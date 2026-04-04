#!/usr/bin/env python3
"""
Optional BLAST validation for sharkmer benchmark amplicons.

Validates amplicon sequences against NCBI nucleotide databases. Uses a local
BLAST database if one is found in /db/ (mounted from ~/db/ on the host),
otherwise falls back to the NCBI remote BLAST API.

Requires network access (for remote fallback) and takes a few minutes per run.

Usage:
    python benchmarks/blast_validate.py benchmarks/results/RESULT.yaml

The script adds a 'blast_hits' field to each product in the result file.
Uses git config user.email for NCBI API identification (remote only).
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib import request, parse, error

import yaml


# NCBI BLAST API endpoint
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
E_VALUE_THRESHOLD = 1e-50
POLL_INTERVAL = 30  # seconds between status checks
MAX_WAIT = 600  # maximum seconds to wait for a single BLAST job

# Local database search path
LOCAL_DB_DIR = Path("/db")


def find_local_blast_db():
    """Discover a local nucleotide BLAST database in /db/.

    Scans /db/ for subdirectories containing .nsq files (nucleotide sequence
    data). Prefers databases with 'core_nt' in the name, then 'nt'.
    Returns the database path (without extension) or None.
    """
    if not LOCAL_DB_DIR.is_dir():
        return None

    # Find all .nsq files (indicates a nucleotide BLAST database)
    nsq_files = list(LOCAL_DB_DIR.rglob("*.nsq"))
    if not nsq_files:
        return None

    # Extract unique database prefixes (strip .nsq or .NN.nsq)
    db_prefixes = set()
    for nsq in nsq_files:
        name = nsq.name
        # Handle multivolume: name.NN.nsq -> name
        # Handle single volume: name.nsq -> name
        stem = name
        if stem.endswith(".nsq"):
            stem = stem[:-4]
        # Strip volume number if present (e.g., core_nt.83 -> keep as-is,
        # it's the db name not a volume number in the old format)
        db_prefixes.add(str(nsq.parent / stem))

    if not db_prefixes:
        return None

    # Prefer core_nt, then nt, then anything else
    for preference in ["core_nt", "nt"]:
        for prefix in sorted(db_prefixes):
            if preference in os.path.basename(prefix):
                return prefix

    # Fall back to first available
    return sorted(db_prefixes)[0]


def check_blastn_available():
    """Check if blastn is installed and accessible."""
    try:
        result = subprocess.run(
            ["blastn", "-version"],
            capture_output=True, text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def verify_local_blast_db(db_path):
    """Verify that a local BLAST database is usable by running a test query.

    Returns True if the database works, False otherwise.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(">test\nATCGATCGATCGATCGATCGATCGATCG\n")
        query_path = f.name

    try:
        result = subprocess.run(
            [
                "blastn",
                "-db", db_path,
                "-query", query_path,
                "-evalue", "10",
                "-max_target_seqs", "1",
                "-outfmt", "6",
            ],
            capture_output=True, text=True, timeout=30
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, Exception):
        return False
    finally:
        os.unlink(query_path)


def blast_sequence_local(sequence, db_path):
    """Run blastn locally against a database. Returns parsed hit or None."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(f">query\n{sequence}\n")
        query_path = f.name

    try:
        result = subprocess.run(
            [
                "blastn",
                "-db", db_path,
                "-query", query_path,
                "-evalue", str(E_VALUE_THRESHOLD),
                "-max_target_seqs", "1",
                "-outfmt", "5",  # XML output
            ],
            capture_output=True, text=True, timeout=120
        )

        if result.returncode != 0:
            print(f"    Local BLAST error: {result.stderr[:200]}")
            return None

        return parse_blast_xml(result.stdout)

    except subprocess.TimeoutExpired:
        print("    Local BLAST timed out")
        return None
    finally:
        os.unlink(query_path)


def get_git_email():
    """Get user email from git config for NCBI API identification."""
    try:
        result = subprocess.run(
            ["git", "config", "user.email"],
            capture_output=True, text=True
        )
        email = result.stdout.strip()
        if email:
            return email
    except Exception:
        pass
    print("WARNING: No git user.email configured. NCBI requests may be throttled.")
    print("  Set with: git config --global user.email you@example.com")
    return "anonymous@example.com"


def submit_blast(sequence, email, program="blastn", database="nt"):
    """Submit a BLAST search and return the request ID (RID)."""
    params = parse.urlencode({
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "EXPECT": str(E_VALUE_THRESHOLD),
        "HITLIST_SIZE": "1",
        "FORMAT_TYPE": "XML",
        "TOOL": "sharkmer_benchmark",
        "EMAIL": email,
    }).encode()

    req = request.Request(BLAST_URL, data=params)
    with request.urlopen(req, timeout=30) as response:
        text = response.read().decode()

    # Parse RID from response
    rid = None
    for line in text.split("\n"):
        if line.strip().startswith("RID = "):
            rid = line.strip().split("= ")[1].strip()
            break

    if not rid:
        raise RuntimeError(f"Failed to get RID from BLAST submission:\n{text[:500]}")

    return rid


def check_blast_status(rid):
    """Check if a BLAST job is complete. Returns 'WAITING', 'READY', or 'UNKNOWN'."""
    params = parse.urlencode({
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo",
    }).encode()

    req = request.Request(BLAST_URL, data=params)
    with request.urlopen(req, timeout=30) as response:
        text = response.read().decode()

    for line in text.split("\n"):
        line = line.strip()
        if line.startswith("Status="):
            return line.split("=")[1].strip()

    return "UNKNOWN"


def get_blast_results(rid):
    """Retrieve BLAST results in XML format."""
    params = parse.urlencode({
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
    }).encode()

    req = request.Request(BLAST_URL, data=params)
    with request.urlopen(req, timeout=60) as response:
        return response.read().decode()


def parse_blast_xml(xml_text):
    """Parse BLAST XML output and extract the top hit."""
    root = ET.fromstring(xml_text)

    iterations = root.findall(".//Iteration")
    if not iterations:
        return None

    iteration = iterations[0]
    hits = iteration.findall(".//Hit")
    if not hits:
        return None

    hit = hits[0]
    hsp = hit.find(".//Hsp")

    if hsp is None:
        return None

    evalue = float(hsp.findtext("Hsp_evalue", "999"))
    if evalue > E_VALUE_THRESHOLD:
        return None

    # Parse hit definition for gene name and taxon
    hit_def = hit.findtext("Hit_def", "")
    accession = hit.findtext("Hit_accession", "")

    identity = int(hsp.findtext("Hsp_identity", "0"))
    align_len = int(hsp.findtext("Hsp_align-len", "1"))
    pct_identity = round(100.0 * identity / align_len, 1)

    return {
        "accession": accession,
        "hit_def": hit_def[:200],  # truncate long descriptions
        "evalue": evalue,
        "pct_identity": pct_identity,
        "align_length": align_len,
    }


def blast_sequence_remote(sequence, email):
    """Submit, wait for, and parse BLAST results via NCBI remote API."""
    rid = submit_blast(sequence, email)
    print(f"    Submitted RID: {rid}")

    elapsed = 0
    while elapsed < MAX_WAIT:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
        status = check_blast_status(rid)
        if status == "READY":
            break
        if status not in ("WAITING",):
            print(f"    Unexpected BLAST status: {status}")
            return None
        print(f"    Waiting... ({elapsed}s)")

    if elapsed >= MAX_WAIT:
        print(f"    BLAST timed out after {MAX_WAIT}s")
        return None

    xml_text = get_blast_results(rid)
    return parse_blast_xml(xml_text)


def submit_blast_batch(fasta_text, email, program="blastn", database="nt"):
    """Submit a multi-FASTA batch and return the RID."""
    params = parse.urlencode({
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": fasta_text,
        "EXPECT": str(E_VALUE_THRESHOLD),
        "HITLIST_SIZE": "1",
        "FORMAT_TYPE": "XML",
        "TOOL": "sharkmer_benchmark",
        "EMAIL": email,
    }).encode()

    req = request.Request(BLAST_URL, data=params)
    with request.urlopen(req, timeout=60) as response:
        text = response.read().decode()

    rid = None
    for line in text.split("\n"):
        if line.strip().startswith("RID = "):
            rid = line.strip().split("= ")[1].strip()
            break

    if not rid:
        raise RuntimeError(f"Failed to get RID from BLAST submission:\n{text[:500]}")

    return rid


def parse_blast_xml_batch(xml_text):
    """Parse BLAST XML output with multiple iterations. Returns dict: query_id -> hit."""
    root = ET.fromstring(xml_text)
    iterations = root.findall(".//Iteration")
    results = {}

    for iteration in iterations:
        # Query ID is in <Iteration_query-def> — we set it via FASTA header
        query_def = iteration.findtext("Iteration_query-def", "").strip()
        # Split to get just the ID (first token)
        query_id = query_def.split()[0] if query_def else ""

        hits = iteration.findall(".//Hit")
        if not hits:
            results[query_id] = None
            continue

        hit = hits[0]
        hsp = hit.find(".//Hsp")
        if hsp is None:
            results[query_id] = None
            continue

        evalue = float(hsp.findtext("Hsp_evalue", "999"))
        if evalue > E_VALUE_THRESHOLD:
            results[query_id] = None
            continue

        hit_def = hit.findtext("Hit_def", "")
        accession = hit.findtext("Hit_accession", "")
        identity = int(hsp.findtext("Hsp_identity", "0"))
        align_len = int(hsp.findtext("Hsp_align-len", "1"))
        pct_identity = round(100.0 * identity / align_len, 1)

        results[query_id] = {
            "accession": accession,
            "hit_def": hit_def[:200],
            "evalue": evalue,
            "pct_identity": pct_identity,
            "align_length": align_len,
        }

    return results


def blast_batch_remote(sequences_with_ids, email):
    """Submit a batch of sequences as multi-FASTA. Returns dict: id -> hit (or None).

    sequences_with_ids: list of (query_id, sequence) tuples.
    """
    fasta_lines = []
    for query_id, seq in sequences_with_ids:
        fasta_lines.append(f">{query_id}")
        fasta_lines.append(seq)
    fasta_text = "\n".join(fasta_lines)

    rid = submit_blast_batch(fasta_text, email)
    print(f"    Submitted batch RID: {rid} ({len(sequences_with_ids)} sequences)")

    elapsed = 0
    while elapsed < MAX_WAIT:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
        status = check_blast_status(rid)
        if status == "READY":
            break
        if status not in ("WAITING",):
            print(f"    Unexpected BLAST status: {status}")
            return {}
        print(f"    Waiting... ({elapsed}s)")

    if elapsed >= MAX_WAIT:
        print(f"    BLAST timed out after {MAX_WAIT}s")
        return {}

    xml_text = get_blast_results(rid)
    return parse_blast_xml_batch(xml_text)


def validate_results(result_path, dry_run=False):
    """Add BLAST validation to amplicons in a benchmark result file.

    Only BLASTs product 0 (the highest-scoring amplicon) per gene to keep
    runtime practical. The primary purpose is confirming gene identity.

    Uses a local BLAST database if available, otherwise falls back to the
    NCBI remote API.
    """
    with open(result_path) as f:
        benchmark = yaml.safe_load(f)

    # Determine BLAST mode: local or remote
    local_db = None
    if check_blastn_available():
        local_db = find_local_blast_db()
        if local_db:
            print(f"Found local BLAST database: {local_db}")
            if verify_local_blast_db(local_db):
                print("  Database verified OK")
            else:
                print("  Database verification failed — falling back to remote API")
                local_db = None

    if local_db:
        mode = "local"
        print(f"Using local BLAST database: {local_db}")
    else:
        mode = "remote"
        if not check_blastn_available():
            print("blastn not found — using NCBI remote API")
        else:
            print("No usable local database — using NCBI remote API")

    email = None
    if mode == "remote":
        email = get_git_email()
        print(f"Using email: {email}")

    results = benchmark.get("results", [])
    total_seqs = sum(
        1 for r in results
        for p in r.get("products", [])
        if p.get("sequences")
    )
    print(f"Total genes to validate (product 0 each): {total_seqs}")

    if dry_run:
        print("Dry run — not submitting to BLAST.")
        return

    if mode == "remote":
        # Batch remote BLAST: collect all sequences with unique IDs, submit in
        # chunks, then map results back to products.
        BATCH_SIZE = 50
        batch = []  # list of (query_id, seq, product_ref, label)
        all_entries = []  # (query_id, product_ref, label)
        for result in results:
            sample = result.get("sample", "?")
            for product in result.get("products", []):
                sequences = product.get("sequences", [])
                if not sequences:
                    continue
                gene = product.get("gene", "?")
                query_id = f"q{len(all_entries)}"
                label = f"{sample} {gene} ({len(sequences[0])} bp)"
                all_entries.append((query_id, sequences[0], product, label))

        print(
            f"Batching {len(all_entries)} sequences into batches of {BATCH_SIZE}"
        )

        for batch_start in range(0, len(all_entries), BATCH_SIZE):
            batch = all_entries[batch_start : batch_start + BATCH_SIZE]
            batch_end = batch_start + len(batch)
            print(
                f"\nBatch {batch_start // BATCH_SIZE + 1}: sequences {batch_start + 1}-{batch_end} of {len(all_entries)}"
            )
            for _, _, _, label in batch:
                print(f"  {label}")

            try:
                sequences_with_ids = [(qid, seq) for qid, seq, _, _ in batch]
                hits_by_id = blast_batch_remote(sequences_with_ids, email)
            except Exception as e:
                print(f"  Batch BLAST error: {e}")
                for _, _, product, _ in batch:
                    product["blast_hit"] = {"error": str(e)}
                time.sleep(1)
                continue

            for query_id, _, product, label in batch:
                hit = hits_by_id.get(query_id)
                if hit:
                    product["blast_hit"] = hit
                    print(
                        f"  {label} -> {hit['accession']} {hit['pct_identity']}% {hit['hit_def'][:60]}"
                    )
                else:
                    product["blast_hit"] = {
                        "accession": None,
                        "hit_def": "no significant hit",
                        "evalue": None,
                        "pct_identity": None,
                    }
                    print(f"  {label} -> no significant hit")

            time.sleep(1)
    else:
        # Local BLAST: one sequence at a time (fast enough)
        seq_count = 0
        for result in results:
            sample = result.get("sample", "?")
            products = result.get("products", [])

            for product in products:
                gene = product.get("gene", "?")
                sequences = product.get("sequences", [])

                if not sequences:
                    continue

                seq = sequences[0]
                seq_count += 1
                print(f"  [{seq_count}/{total_seqs}] {sample} {gene} ({len(seq)} bp)")

                try:
                    hit = blast_sequence_local(seq, local_db)
                    if hit:
                        product["blast_hit"] = hit
                        print(
                            f"    Hit: {hit['accession']} {hit['pct_identity']}% {hit['hit_def'][:80]}"
                        )
                    else:
                        product["blast_hit"] = {
                            "accession": None,
                            "hit_def": "no significant hit",
                            "evalue": None,
                            "pct_identity": None,
                        }
                        print(f"    No significant hit")
                except Exception as e:
                    print(f"    BLAST error: {e}")
                    product["blast_hit"] = {"error": str(e)}

    # Record which BLAST mode was used
    benchmark["blast_mode"] = mode
    if local_db:
        benchmark["blast_db"] = os.path.basename(local_db)

    # Write updated results back
    with open(result_path, "w") as f:
        yaml.dump(benchmark, f, default_flow_style=False, sort_keys=False)

    print(f"\nBLAST validation written to: {result_path}")


def main():
    parser = argparse.ArgumentParser(
        description="BLAST-validate sharkmer benchmark amplicons"
    )
    parser.add_argument(
        "result_file",
        help="Benchmark result YAML file to validate"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Count sequences but don't submit to BLAST"
    )
    args = parser.parse_args()

    result_path = Path(args.result_file)
    if not result_path.exists():
        print(f"File not found: {result_path}")
        sys.exit(1)

    validate_results(result_path, dry_run=args.dry_run)


if __name__ == "__main__":
    main()

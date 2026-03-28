#!/usr/bin/env python3
"""
Optional BLAST validation for sharkmer benchmark amplicons.

Batch-submits all amplicon sequences from a benchmark result to NCBI BLAST,
parses results, and annotates the benchmark YAML with validation data.

Requires network access and takes a few minutes per run.

Usage:
    python benchmarks/blast_validate.py benchmarks/results/RESULT.yaml

The script adds a 'blast_hits' field to each product in the result file.
Uses git config user.email for NCBI API identification.
"""

import argparse
import subprocess
import sys
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


def blast_sequence(sequence, email):
    """Submit, wait for, and parse BLAST results for a single sequence."""
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


def validate_results(result_path, dry_run=False):
    """Add BLAST validation to all amplicons in a benchmark result file."""
    with open(result_path) as f:
        benchmark = yaml.safe_load(f)

    email = get_git_email()
    print(f"Using email: {email}")

    results = benchmark.get("results", [])
    total_seqs = sum(
        len(p.get("sequences", []))
        for r in results
        for p in r.get("products", [])
    )
    print(f"Total sequences to validate: {total_seqs}")

    if dry_run:
        print("Dry run — not submitting to BLAST.")
        return

    seq_count = 0
    for result in results:
        sample = result.get("sample", "?")
        max_reads = result.get("max_reads", "?")
        products = result.get("products", [])

        for product in products:
            gene = product.get("gene", "?")
            sequences = product.get("sequences", [])

            if not sequences:
                continue

            blast_hits = []
            for i, seq in enumerate(sequences):
                seq_count += 1
                print(f"  [{seq_count}/{total_seqs}] {sample} {gene} product {i} ({len(seq)} bp)")

                try:
                    hit = blast_sequence(seq, email)
                    if hit:
                        blast_hits.append({
                            "product_index": i,
                            **hit,
                        })
                        print(f"    Hit: {hit['accession']} {hit['pct_identity']}% {hit['hit_def'][:80]}")
                    else:
                        blast_hits.append({
                            "product_index": i,
                            "accession": None,
                            "hit_def": "no significant hit",
                            "evalue": None,
                            "pct_identity": None,
                        })
                        print(f"    No significant hit")
                except Exception as e:
                    print(f"    BLAST error: {e}")
                    blast_hits.append({
                        "product_index": i,
                        "error": str(e),
                    })

                # NCBI rate limit: max 1 request per second (we already wait 30s in polling)
                time.sleep(1)

            product["blast_hits"] = blast_hits

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

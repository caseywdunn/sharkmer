#!/usr/bin/env python3
"""
Bootstrap reference sequences by running panel validations, collecting
amplicons, and BLASTing them against NCBI nt with taxonomic restriction.

For each panel:
  1. Run sharkmer at the highest declared max_reads for each validation sample
  2. Collect the longest amplicon per gene per sample
  3. BLAST each against NCBI nt (restricted to the sample's taxonomy)
  4. Output a TSV and YAML for evaluation

Usage:
    python scripts/bootstrap_from_runs.py panels/cnidaria.yaml
    python scripts/bootstrap_from_runs.py panels/cnidaria.yaml --reuse-runs
    python scripts/bootstrap_from_runs.py --all
    python scripts/bootstrap_from_runs.py --all --reuse-runs --no-blast
"""

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from urllib import error, parse, request

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from sharkmer_validate import runner

RUNS_DIR = REPO_ROOT / "panels" / "validation_runs"
BOOTSTRAP_DIR = REPO_ROOT / "panels" / "bootstrap"

# ---------------------------------------------------------------------------
# NCBI BLAST with taxonomic restriction
# ---------------------------------------------------------------------------

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
POLL_INTERVAL = 30
MAX_WAIT = 900  # 15 min per job


def get_git_email() -> str:
    try:
        r = runner.subprocess.run(
            ["git", "config", "user.email"], capture_output=True, text=True
        )
        return r.stdout.strip() or "anonymous@example.com"
    except Exception:
        return "anonymous@example.com"


def submit_blast(sequence: str, email: str, taxon: str = "") -> str | None:
    """Submit a BLAST job to NCBI, optionally restricted to a taxon."""
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "nt",
        "QUERY": sequence,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": "5",
        "EXPECT": "1e-10",
        "TOOL": "sharkmer-bootstrap",
        "EMAIL": email,
    }
    if taxon:
        # Use ENTREZ_QUERY for taxonomic restriction.
        params["ENTREZ_QUERY"] = f'"{taxon}"[Organism]'

    data = parse.urlencode(params).encode()
    req = request.Request(BLAST_URL, data=data)
    try:
        with request.urlopen(req, timeout=60) as resp:
            text = resp.read().decode()
    except Exception as e:
        print(f"    Submit failed: {e}")
        return None

    match = re.search(r"RID = (\S+)", text)
    return match.group(1) if match else None


def check_blast_status(rid: str) -> str:
    params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
    url = f"{BLAST_URL}?{parse.urlencode(params)}"
    try:
        with request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode()
    except Exception:
        return "UNKNOWN"

    if "Status=WAITING" in text:
        return "WAITING"
    elif "Status=READY" in text:
        return "READY"
    elif "Status=FAILED" in text:
        return "FAILED"
    return "UNKNOWN"


def get_blast_results(rid: str) -> str:
    params = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid}
    url = f"{BLAST_URL}?{parse.urlencode(params)}"
    with request.urlopen(url, timeout=60) as resp:
        return resp.read().decode()


def parse_blast_xml(xml_text: str) -> dict | None:
    """Extract top hit from BLAST XML."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None

    hits = root.findall(".//Hit")
    if not hits:
        return None

    hit = hits[0]
    hsp = hit.find(".//Hsp")
    if hsp is None:
        return None

    identity = int(hsp.findtext("Hsp_identity", "0"))
    align_len = int(hsp.findtext("Hsp_align-len", "1"))

    return {
        "accession": hit.findtext("Hit_accession", ""),
        "hit_def": hit.findtext("Hit_def", ""),
        "pct_identity": round(100.0 * identity / align_len, 1),
        "align_length": align_len,
    }


def blast_sequence_remote(sequence: str, email: str, taxon: str = "") -> dict | None:
    """Submit, poll, and return the top BLAST hit."""
    # Taxonomic restriction: use genus only for broader matching.
    genus = taxon.split()[0] if taxon else ""

    rid = submit_blast(sequence, email, taxon=genus)
    if not rid:
        return None

    print(f"      RID={rid}", end="", flush=True)
    elapsed = 0
    while elapsed < MAX_WAIT:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL
        status = check_blast_status(rid)
        print(".", end="", flush=True)
        if status == "READY":
            break
        elif status == "FAILED":
            print(" FAILED")
            return None
    else:
        print(" TIMEOUT")
        return None

    print(" done")
    xml = get_blast_results(rid)
    return parse_blast_xml(xml)


# ---------------------------------------------------------------------------
# Collect amplicons from sharkmer runs
# ---------------------------------------------------------------------------


def collect_amplicons_from_runs(
    panel_data: dict, run_dir: Path, panel_name: str
) -> list:
    """Collect best amplicon per (gene, sample) from a run directory.

    Returns list of dicts: {gene, taxon, accession, sequence, length, max_reads}
    """
    validation = panel_data.get("validation") or {}
    samples = validation.get("samples", [])
    gene_names = runner.panel_gene_names(panel_data)

    amplicons = []
    for sample in samples:
        accession = sample["accession"]
        taxon = sample.get("taxon", "")
        max_reads_list = sorted(sample.get("max_reads", []), reverse=True)

        # Try each depth from highest to lowest, take the first successful one.
        for max_reads in max_reads_list:
            k_reads = max_reads // 1000
            prefix = f"{panel_name}_{accession}_{k_reads}k"
            products = runner.parse_fasta_products(prefix, run_dir)
            if not products:
                continue

            for prod in products:
                gene = prod["gene"]
                # Strip panel prefix.
                if gene.startswith(f"{panel_name}_"):
                    gene = gene[len(f"{panel_name}_"):]

                seqs = prod.get("sequences", [])
                if not seqs:
                    continue

                # Only take the first (longest) product per gene per sample.
                # Check if we already have this gene+accession.
                existing = [
                    a for a in amplicons
                    if a["gene"] == gene and a["accession"] == accession
                ]
                if existing:
                    continue

                amplicons.append({
                    "gene": gene,
                    "taxon": taxon,
                    "accession": accession,
                    "sequence": seqs[0],
                    "length": len(seqs[0]),
                    "max_reads": max_reads,
                })

            # If we got products at this depth, don't try lower depths.
            break

    return amplicons


def find_latest_run_dir(panel_name: str) -> Path | None:
    """Find the most recent validation run directory for a panel."""
    candidates = sorted(RUNS_DIR.glob(f"{panel_name}_*"), reverse=True)
    for d in candidates:
        if d.is_dir() and list(d.glob("*.fasta")):
            return d
    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_panel(
    panel_path: Path,
    reuse_runs: bool = False,
    skip_blast: bool = False,
):
    """Run the full bootstrap pipeline for one panel."""
    panel_data = runner.load_panel_yaml(panel_path)
    panel_name = panel_data.get("name", "unknown")
    validation = panel_data.get("validation") or {}
    samples = validation.get("samples", [])

    if not samples:
        print(f"No validation samples in {panel_name}, skipping.")
        return

    print(f"\n{'='*60}")
    print(f"Panel: {panel_name} ({len(samples)} samples)")
    print(f"{'='*60}")

    # Step 1: Get amplicons (run sharkmer or reuse).
    run_dir = None
    if reuse_runs:
        run_dir = find_latest_run_dir(panel_name)
        if run_dir:
            print(f"Reusing existing run: {run_dir}")
        else:
            print(f"No existing runs for {panel_name}, running sharkmer...")

    if run_dir is None:
        runner.build_sharkmer()
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = RUNS_DIR / f"{panel_name}_{stamp}"

        for sample in samples:
            accession = sample["accession"]
            # Run at highest depth only.
            max_reads = max(sample.get("max_reads", [1_000_000]))
            print(f"\n--- {sample.get('taxon', accession)} @ {max_reads // 1000}k ---")
            runner.run_sharkmer(
                panel_path, panel_name, accession, max_reads, run_dir
            )

    # Step 2: Collect amplicons.
    amplicons = collect_amplicons_from_runs(panel_data, run_dir, panel_name)
    print(f"\nCollected {len(amplicons)} amplicons")

    if not amplicons:
        print("No amplicons to process.")
        return

    # Step 3: BLAST with taxonomic restriction.
    if not skip_blast:
        email = get_git_email()
        print(f"\nBLASTing {len(amplicons)} amplicons against NCBI nt...")
        print(f"(taxonomically restricted to genus of each sample)")
        for i, amp in enumerate(amplicons, 1):
            gene = amp["gene"]
            taxon = amp["taxon"]
            accession = amp["accession"]
            print(
                f"\n  [{i}/{len(amplicons)}] {gene} from {taxon} "
                f"({accession}, {amp['length']} bp)"
            )
            hit = blast_sequence_remote(amp["sequence"], email, taxon=taxon)
            if hit:
                amp["blast_hit"] = hit
                print(
                    f"    -> {hit['accession']} {hit['pct_identity']}% "
                    f"{hit['hit_def'][:70]}"
                )
            else:
                amp["blast_hit"] = None
                print("    -> no hit")

    # Step 4: Write outputs.
    BOOTSTRAP_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # TSV for quick evaluation.
    tsv_path = BOOTSTRAP_DIR / f"{panel_name}_{stamp}.tsv"
    with open(tsv_path, "w") as f:
        f.write(
            "panel\tgene\ttaxon\tsample_accession\tlength\t"
            "blast_accession\tblast_identity\tblast_description\n"
        )
        for amp in amplicons:
            hit = amp.get("blast_hit")
            if hit:
                f.write(
                    f"{panel_name}\t{amp['gene']}\t{amp['taxon']}\t"
                    f"{amp['accession']}\t{amp['length']}\t"
                    f"{hit['accession']}\t{hit['pct_identity']}\t"
                    f"{hit['hit_def']}\n"
                )
            else:
                f.write(
                    f"{panel_name}\t{amp['gene']}\t{amp['taxon']}\t"
                    f"{amp['accession']}\t{amp['length']}\t"
                    f"---\t---\t---\n"
                )
    print(f"\nTSV: {tsv_path}")

    # YAML with sequences for reference population.
    yaml_path = BOOTSTRAP_DIR / f"{panel_name}_{stamp}.yaml"
    output = {
        "panel": panel_name,
        "date": datetime.now().strftime("%Y-%m-%d"),
        "amplicons": [],
    }
    for amp in amplicons:
        entry = {
            "gene": amp["gene"],
            "taxon": amp["taxon"],
            "sample_accession": amp["accession"],
            "length": amp["length"],
            "sequence": amp["sequence"],
        }
        hit = amp.get("blast_hit")
        if hit:
            entry["blast"] = hit
        output["amplicons"].append(entry)

    with open(yaml_path, "w") as f:
        yaml.dump(output, f, default_flow_style=False, sort_keys=False, width=4096)
    print(f"YAML: {yaml_path}")

    return amplicons


def main():
    parser = argparse.ArgumentParser(
        description="Bootstrap references by running validations and BLASTing results."
    )
    parser.add_argument(
        "panels", nargs="*",
        help="Panel YAML files to process",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Process all panels that have validation samples",
    )
    parser.add_argument(
        "--reuse-runs", action="store_true",
        help="Reuse existing validation run outputs instead of re-running sharkmer",
    )
    parser.add_argument(
        "--no-blast", action="store_true",
        help="Skip BLAST; just collect amplicons",
    )
    args = parser.parse_args()

    if args.all:
        panel_paths = sorted(runner.PANELS_DIR.glob("*.yaml"))
    elif args.panels:
        panel_paths = [Path(p).resolve() for p in args.panels]
    else:
        parser.print_help()
        sys.exit(1)

    for panel_path in panel_paths:
        if not panel_path.exists():
            print(f"Not found: {panel_path}")
            continue
        run_panel(panel_path, reuse_runs=args.reuse_runs, skip_blast=args.no_blast)

    print("\nDone. Review TSV files in panels/bootstrap/ to evaluate results.")


if __name__ == "__main__":
    main()

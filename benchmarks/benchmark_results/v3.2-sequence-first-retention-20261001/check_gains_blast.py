#!/usr/bin/env python3
"""BLAST the 200 newly retained identities against pinned public references."""

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from scripts.sharkmer_validate import blast_references, reference_provenance
from scripts.sharkmer_validate.reference_targets import logical_gene_name


SAMPLES = {
    "SRR25378184": ("angiospermae", "Liriodendron tulipifera", "005"),
    "SRR31887760": ("insecta", "Drosophila melanogaster", "007"),
    "SRR1057608": ("insecta", "Heliconius pachinus", "008"),
    "SRR27962769": ("insecta", "Gryllus bimaculatus", "009"),
}


def fasta_records(path):
    header, pieces = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(pieces)
            header, pieces = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError(f"Sequence precedes header: {path}")
            pieces.append(line.strip())
    if header is not None:
        yield header, "".join(pieces)


def load_gains():
    path = Path(__file__).with_name("report-v3") / "uncertain_gains.tsv"
    with path.open(newline="") as stream:
        gains = list(csv.DictReader(stream, delimiter="\t"))
    if len(gains) != 200 or len({(r["accession"], r["gene"], r["sha256"]) for r in gains}) != 200:
        raise ValueError("Expected exactly 200 unique gained identities")
    return gains


def load_sequences(gains):
    results = []
    cache = {}
    for gain in gains:
        accession, gene = gain["accession"], gain["gene"]
        panel, taxon, number = SAMPLES[accession]
        key = accession, gene
        if key not in cache:
            glob = (ROOT / "tmp/sequence-retention-20261001/full-v3/runs").glob(
                f"{number}_{panel}_{accession}_1000000_pair1_candidate/outputs/*_{gene}.fasta"
            )
            paths = list(glob)
            if len(paths) != 1:
                raise ValueError(f"Expected one candidate FASTA for {key}: {paths}")
            cache[key] = {
                hashlib.sha256(sequence.encode("ascii")).hexdigest(): (header, sequence)
                for header, sequence in fasta_records(paths[0])
            }
        digest = gain["sha256"]
        if digest not in cache[key]:
            raise ValueError(f"Gained identity absent from candidate FASTA: {key} {digest}")
        header, sequence = cache[key][digest]
        if (len(sequence) != int(gain["length"]) or "assembly_status=uncertain" not in header
                or "marker_detail=summary_only" not in header):
            raise ValueError(f"Gained sequence/header differs: {key} {digest}")
        results.append({"accession": accession, "panel": panel, "taxon": taxon,
                        "gene": gene, "length": len(sequence), "sha256": digest,
                        "header": header, "sequence": sequence})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = load_sequences(load_gains())
    panels = {name: yaml.safe_load((ROOT / "panels" / f"{name}.yaml").read_text())
              for name in {row["panel"] for row in rows}}
    databases = {}
    reference_audits = {}
    for name, panel in panels.items():
        audit = reference_provenance.audit_references(panel)
        if audit["excluded_count"] or not audit["verified_count"]:
            raise ValueError(f"Panel reference audit incomplete: {name}")
        reference_audits[name] = {"verified_count": audit["verified_count"],
                                  "audit_sha256": audit["audit_sha256"],
                                  "catalog": audit["catalog"]}
        path = args.work_dir / name
        path.mkdir(exist_ok=True)
        db = blast_references.build_reference_db(panel, path)
        if db is None:
            raise ValueError(f"BLAST database unavailable: {name}")
        databases[name] = db
    results = []
    for i, row in enumerate(rows, 1):
        panel = panels[row["panel"]]
        target = row["gene"].removeprefix(f"{row['panel']}_")
        expected = logical_gene_name(target, panel)
        available = {logical_gene_name(ref["gene_name"], panel)
                     for ref in reference_provenance.verified_references(panel)}
        if expected in available:
            result = asdict(blast_references.blast_against_references(
                row["sequence"], databases[row["panel"]], target, row["taxon"]))
        else:
            result = asdict(blast_references.RefBlastResult(
                "no_verified_reference", target, row["taxon"],
                expected_logical_gene=expected))
        if result["status"] == "failed_run":
            raise RuntimeError(f"BLAST failed for {row['sha256']}: {result['error']}")
        results.append({key: value for key, value in row.items() if key != "sequence"}
                       | {"reference_match": result})
        if i % 25 == 0:
            print(f"Classified {i}/{len(rows)}", flush=True)
    summary = {
        "schema_version": 1,
        "scope": "200 exact new sequence identities from corrected pair-1 candidate outputs",
        "reference_scope": "verified public-record panel regions in pinned local catalog",
        "criteria": {"evalue": "1e-10", "min_query_coverage_pct": 90,
                     "min_identity_pct": 90, "task": "blastn default"},
        "reference_audits": reference_audits,
        "total": len(results),
        "statuses": dict(Counter(r["reference_match"]["status"] for r in results)),
        "by_target": {
            f"{accession}:{gene}": {
                "count": len(group),
                "statuses": dict(Counter(r["reference_match"]["status"] for r in group)),
            }
            for accession, gene in sorted({(r["accession"], r["gene"]) for r in results})
            for group in [[r for r in results if r["accession"] == accession and r["gene"] == gene]]
        },
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (args.out_dir / "per_sequence.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    with (args.out_dir / "per_sequence.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        columns = ["accession", "gene", "length", "sha256", "status", "matched_gene",
                   "matched_taxon", "matched_accession", "pct_identity", "query_coverage_pct",
                   "reference_coverage_pct", "gap_count", "query_gap_bases", "reference_gap_bases",
                   "sequence_relationship", "alignment_structure_status"]
        writer.writerow(columns)
        for row in results:
            match = row["reference_match"]
            writer.writerow([row.get(key, match.get(key)) for key in columns])
    print(json.dumps(summary["by_target"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

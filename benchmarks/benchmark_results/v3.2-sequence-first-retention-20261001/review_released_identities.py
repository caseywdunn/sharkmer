#!/usr/bin/env python3
"""Check historical release identities directly against corrected raw run receipts."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PRIOR = (ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929"
         / "results/report-v2/released_identity_dispositions.tsv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    directory = args.directory.resolve()
    jobs = {}
    for path in directory.glob("*pair1_candidate.json"):
        job = json.loads(path.read_text())
        if job["status"] != "pass":
            raise ValueError("Failed candidate run")
        jobs[job["invocation"]["input"]] = job
    if len(jobs) != 10:
        raise ValueError("Need ten completed candidate samples")
    with PRIOR.open(newline="") as stream:
        historical = list(csv.DictReader(stream, delimiter="\t"))
    if len(historical) != 86:
        raise ValueError("Historical identity inventory changed")
    records = []
    for row in historical:
        matches = []
        for product in jobs[row["accession"]]["outputs"]["records"]:
            if product["gene"] != row["gene"]:
                continue
            orientation = ("exact" if product["sha256"] == row["sha256"] else "reverse_complement"
                           if product["reverse_complement_sha256"] == row["sha256"] else None)
            if orientation:
                matches.append({"assembly_status": product["assembly_status"],
                                "orientation": orientation, "sha256": product["sha256"]})
        if len(matches) != 1:
            raise ValueError(f"Historical identity has {len(matches)} matches: {row['accession']}/{row['gene']}")
        records.append({"accession": row["accession"], "gene": row["gene"], "scope": row["scope"],
                        "length": int(row["length"]), "historical_sha256": row["sha256"],
                        "prior_disposition": row["disposition"], "current": matches[0]})
    high = [row for row in records if row["scope"] == "high_copy_candidate"]
    deferred = [row for row in records if row["scope"] == "deferred_or_unclassified"]
    if len(high) != 71 or len(deferred) != 15:
        raise ValueError("Historical scope counts changed")
    counts = {scope: dict(Counter(row["current"]["assembly_status"] for row in selected))
              for scope, selected in (("high_copy_candidate", high), ("deferred_or_unclassified", deferred))}
    newly_available = [row for row in records if row["prior_disposition"] == "unobserved"]
    if len(newly_available) != 14 or any(row["current"] != {
            "assembly_status": "uncertain", "orientation": "exact", "sha256": row["historical_sha256"]}
            for row in newly_available):
        raise ValueError("Previously absent identities did not recover exactly as uncertain")
    result = {"status": "pass", "historical_total": len(records), "scope_counts": counts,
              "previously_absent_now_uncertain_exact": newly_available,
              "all_records": records}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "historical_total": len(records),
                      "scope_counts": counts, "formerly_absent": len(newly_available)}))


if __name__ == "__main__":
    main()

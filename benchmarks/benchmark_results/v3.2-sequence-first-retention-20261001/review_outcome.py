#!/usr/bin/env python3
"""Independent raw-output audit of the frozen retention comparisons."""

import argparse
import hashlib
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def reverse_complement(value):
    return value.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def raw_fasta(path):
    rows, header, sequence = [], None, []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if header is not None:
                rows.append((header, "".join(sequence)))
            header, sequence = raw[1:], []
        elif raw.strip():
            need(header is not None, "FASTA sequence before header")
            sequence.append(raw.strip())
    if header is not None:
        rows.append((header, "".join(sequence)))
    return rows


def audit_run(job):
    invocation = job["invocation"]
    folder = Path(invocation["outputs"])
    name = invocation["sample_prefix"]
    stats = yaml.safe_load((folder / (name + ".stats.yaml")).read_text())
    result = {"standard": set(), "uncertain": set(), "uncertain_orientation": [],
              "records": 0, "summary_sources": 0,
              "detail_sources": 0, "target_assessments": {}, "counts": {}}
    for target in stats["pcr_results"]:
        gene = target["gene_name"]
        filename = target.get("output_file")
        rows = raw_fasta(folder / filename) if filename else []
        need(len(rows) == target["n_products"], "FASTA count differs from stats")
        labels = []
        for header, sequence in rows:
            tokens = header.split()
            fields = dict(token.split("=", 1) for token in tokens[1:])
            label = fields.get("assembly_status")
            need(label in ("standard", "uncertain") and sequence and set(sequence) <= set("ACGT"),
                 "Invalid raw FASTA class or sequence")
            if label == "uncertain":
                need(fields.get("repeat_length") == "unresolved"
                     and fields.get("read_support") == "not_evaluated",
                     "Uncertain header omits uncertainty")
                if invocation["version"] == "candidate":
                    need(fields.get("marker_detail") in ("complete", "summary_only"),
                         "New uncertain header omits marker detail")
            result[label].add((gene, digest(sequence)))
            if label == "uncertain":
                result["uncertain_orientation"].append(
                    {"gene": gene, "sha256": digest(sequence),
                     "reverse_complement_sha256": digest(reverse_complement(sequence))})
            labels.append(label)
        need(labels == ["standard"] * target["n_standard_products"]
             + ["uncertain"] * target["n_uncertain_products"], "Class order/count mismatch")
        result["records"] += len(rows)
        assessment = target.get("uncertain_products")
        if assessment:
            result["target_assessments"][gene] = assessment
            for product in assessment["products"]:
                need((gene, product["sequence_sha256"]) in result["uncertain"],
                     "Structured uncertain product absent from FASTA")
                for source in product["sources"]:
                    if source.get("marker_detail") == "summary_only":
                        need(not source["marker_occurrences"] and not source["marker_covered_runs"],
                             "Summary source includes detailed markers")
                        result["summary_sources"] += 1
                    else:
                        result["detail_sources"] += 1
    result["counts"] = {key: stats[key] for key in
                        ("n_reads_read", "n_bases_read", "n_unique_kmers", "n_kmers_above_threshold")
                        if key in stats}
    result["standard"] = sorted(result["standard"])
    result["uncertain"] = sorted(result["uncertain"])
    return result


def audit(directory, mode):
    state = json.loads((directory / "result.json").read_text())
    need(state["status"] == "pass" and state["mode"] == mode, "Comparison did not pass")
    jobs = state["jobs"]
    need(len(jobs) == (60 if mode == "full" else 2), "Wrong invocation count")
    by_cell = {}
    for job in jobs:
        need(job["status"] == "pass", "Failed invocation")
        key = job["invocation"]["input"]
        by_cell.setdefault(key, {}).setdefault(job["invocation"]["version"], []).append(audit_run(job))
    cells = []
    for accession, roles in by_cell.items():
        old, new = (roles[role][0] for role in ("baseline", "candidate"))
        for role, runs in roles.items():
            need(len(runs) == (3 if mode == "full" else 1), "Wrong repetition count")
            need(all({key: row[key] for key in ("standard", "uncertain", "counts")}
                     == {key: runs[0][key] for key in ("standard", "uncertain", "counts")}
                     for row in runs), f"{accession}/{role}: unstable raw output")
        need(old["standard"] == new["standard"] and old["counts"] == new["counts"],
             f"{accession}: standard or input changed")
        need(set(old["uncertain"]) <= set(new["uncertain"]), f"{accession}: lost uncertain output")
        cells.append({"accession": accession, "standard": len(new["standard"]),
                      "old_uncertain": len(old["uncertain"]), "new_uncertain": len(new["uncertain"]),
                      "new_uncertain_identities": sorted(set(new["uncertain"]) - set(old["uncertain"])),
                      "summary_sources": new["summary_sources"], "detail_sources": new["detail_sources"]})
    summary = json.loads((directory / "summary.json").read_text())
    need(summary["status"] == "pass", "Summary did not pass")
    need(summary["totals"] == {"standard": sum(row["standard"] for row in cells),
                               "old_uncertain": sum(row["old_uncertain"] for row in cells),
                               "new_uncertain": sum(row["new_uncertain"] for row in cells),
                               "gained_uncertain": sum(len(row["new_uncertain_identities"]) for row in cells)},
         "Summary totals differ from independently parsed FASTA")
    result = {"status": "pass", "mode": mode, "cells": cells, "totals": summary["totals"]}
    if "SRR31887760" in by_cell:
        dros = next(row for row in cells if row["accession"] == "SRR31887760")
        historical = json.loads((OLD / "review/dros12s_outcome_review.json").read_text())
        expected = {(row["gene"], row["sha256"]) for row in historical["historical_exact_matches"]
                    if row["length"] in (487, 499)}
        need(expected <= set(dros["new_uncertain_identities"]), "Historical Drosophila identities absent")
        oriented = by_cell["SRR31887760"]["candidate"][0]["uncertain_orientation"]
        orientation = [{"gene": gene, "sha256": digest_value,
                        "literal_matches": sum(row["gene"] == gene and row["sha256"] == digest_value
                                               for row in oriented),
                        "reverse_complement_matches": sum(row["gene"] == gene
                                                           and row["reverse_complement_sha256"] == digest_value
                                                           for row in oriented)}
                       for gene, digest_value in sorted(expected)]
        need(all(row["literal_matches"] == 1 for row in orientation),
             "Historical Drosophila sequences lack unique literal matches")
        assessment = by_cell["SRR31887760"]["candidate"][0]["target_assessments"]["insecta_12S"]
        threshold = next(row["retention"] for row in assessment["thresholds"] if row["threshold"] == 2)
        need(threshold["retained_paths"] == 26 and threshold["retained_summary_only_paths"] == 23
             and threshold["observed_not_retained_paths"] == 0 and assessment["marker_detail_truncated"]
             and not assessment["retention_truncated"], "Drosophila retention accounting differs")
        result.update(dros12s_exact_recovered=sorted(expected), dros12s_orientations=orientation,
                      dros12s_threshold2=threshold)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("dros", "plant", "full"), required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.directory.resolve(), args.mode)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "mode": args.mode,
                      "totals": result["totals"]}))


if __name__ == "__main__":
    main()

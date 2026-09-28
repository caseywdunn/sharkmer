#!/usr/bin/env python3
"""Derive compact tables from completed, hash-bound comparison receipts.

This reads saved results only: it neither reads FASTQ inputs nor runs Sharkmer.
The separate results_review auditor checks the raw stats and outcome semantics.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import stat


HERE = Path(__file__).absolute().parent
BUNDLE = HERE.parent
ACCESSIONS = ("SRR27962769", "SRR31887760", "SRR1057608")
ROLES = ("baseline_off", "changed_off", "changed_on")
COUNTS = ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers")
RETENTION_SUMS = (
    "observed_withheld_paths", "retained_paths", "observed_not_retained_paths",
    "retained_sequence_bases", "retained_marker_occurrences", "dropped_oversize_sequence",
    "dropped_oversize_marker_set", "dropped_path_cap", "dropped_sequence_base_cap",
    "dropped_marker_cap", "node_visit_skips",
)
SEARCH_FLAGS = ("node_budget_reached", "dfs_limit_reached", "max_length_reached", "path_quota_reached")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).absolute()
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode), f"Not a regular file: {path}")
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    h = hashlib.sha256()
    with path.open("rb") as handle:
        require(identity(os.fstat(handle.fileno())) == identity(before), f"Replaced before read: {path}")
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
        require(identity(os.fstat(handle.fileno())) == identity(before), f"Changed while reading: {path}")
    require(identity(path.lstat()) == identity(before), f"Replaced after read: {path}")
    return {"path": str(path), "sha256": h.hexdigest(), "size_bytes": before.st_size}


class Evidence:
    def __init__(self):
        self.sources = {}

    def bind(self, path, expected=None):
        actual = fingerprint(path)
        require(expected is None or actual == expected, f"Fingerprint differs: {path}")
        require(actual["path"] not in self.sources or self.sources[actual["path"]] == actual,
                f"Evidence changed between reads: {path}")
        self.sources[actual["path"]] = actual
        return actual

    def json(self, path, expected=None):
        self.bind(path, expected)
        return json.loads(Path(path).read_text())

    def finish(self):
        for row in self.sources.values():
            require(fingerprint(row["path"]) == row, f"Evidence changed before publication: {row['path']}")


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            records.append({"header": line[1:], "sequence": ""})
        elif line:
            require(bool(records), f"Sequence before FASTA header: {path}")
            records[-1]["sequence"] += line
    require(all(row["sequence"] for row in records), f"Empty FASTA sequence: {path}")
    return records


def rc(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_new(path, text):
    with Path(path).open("x") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def tsv(rows, fields):
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: json.dumps(row[key], separators=(",", ":"))
                         if isinstance(row[key], (dict, list, bool)) else row[key] for key in fields})
    return output.getvalue()


def derive(attempt, evidence):
    outer = evidence.json(attempt / "orchestration.json")
    require(outer["status"] == "complete", "Orchestration did not complete")
    execution = evidence.json(attempt / "runs/execution.json", outer["execution"])
    require(execution["status"] == "complete", "Execution did not complete")
    frozen = evidence.json(attempt / "freeze.json", outer["freeze"])
    require(execution["freeze_before"] == execution["freeze_after"] == outer["freeze"], "Freeze bookends differ")
    assessment = evidence.json(attempt / "runs/assessment.json", execution["assessment"])
    require(assessment["status"] == "pass", "Standard parity did not pass")
    losses_path = BUNDLE / "inputs/historical_losses.json"
    losses = evidence.json(losses_path, frozen["sources"][str(losses_path)])["rows"]
    require(len(losses) == 7, "Historical identity set differs")
    require(len(execution["jobs"]) == len(frozen["jobs"]) == 9, "Expected nine runs")
    require([j["id"] for j in execution["jobs"]] == [j["id"] for j in frozen["jobs"]], "Run order differs")
    records = {}
    run_rows = []
    gene_rows = []
    product_rows = []
    for job, item in zip(frozen["jobs"], execution["jobs"]):
        record = evidence.json(attempt / "runs" / (job["id"] + ".receipt.json"), item["receipt"])
        require(item["status"] == record["status"] == "complete", "Incomplete job")
        require(all(record[key] == value for key, value in job.items()), "Receipt job differs from freeze")
        key = (record["accession"], record["role"])
        require(key not in records, "Duplicate sample/role")
        records[key] = record
        sample = frozen["samples"][record["accession"]]
        require(record["prefix_before"] == record["prefix_after"] == sample["expected_prefix"], "Prefix bookends differ")
        output = record["outputs"]
        for row in output["files"].values():
            evidence.bind(row["path"], row)
        standard, uncertain = [], []
        require(len(output["stats"]["pcr_results"]) == 21, "Expected 21 target rows")
        require(len({g["gene_name"] for g in output["stats"]["pcr_results"]}) == 21, "Duplicate target rows")
        for gene in output["stats"]["pcr_results"]:
            name = gene["gene_name"]
            if gene.get("output_file"):
                ordinary = fasta(Path(record["output"]) / gene["output_file"])
            else:
                ordinary = []
            require(len(ordinary) == gene["n_products"] and [len(r["sequence"]) for r in ordinary] == gene.get("product_lengths", []), "Standard FASTA/count/length mismatch")
            standard += [dict(gene=name, **row) for row in ordinary]
            if record["role"] != "changed_on":
                require("uncertain_products" not in gene, "Disabled output has uncertain payload")
                continue
            u = gene["uncertain_products"]
            emitted = fasta(Path(record["output"]) / u["output_file"]) if u.get("output_file") else []
            uncertain += [dict(gene=name, **row) for row in emitted]
            require(len(emitted) == u["n_products"] == u["unique_products"] == len(u["products"]), "Uncertain product count mismatch")
            require([len(r["sequence"]) for r in emitted] == u["product_lengths"], "Uncertain product lengths differ")
            require(u["retained_path_records"] == u["excluded_standard_equivalent_path_records"] + u["exact_duplicate_path_records"] + u["unique_products"], "Retention accounting differs")
            totals = {k: sum(t["retention"][k] for t in u["thresholds"]) for k in RETENTION_SUMS}
            require(totals["retained_paths"] == u["retained_path_records"], "Threshold/gene retained counts differ")
            require(totals["observed_withheld_paths"] == totals["retained_paths"] + totals["observed_not_retained_paths"], "Observed path accounting differs")
            gene_rows.append({
                "accession": record["accession"], "taxon": sample["sample"], "gene": name,
                "standard_count": gene["n_products"], "standard_lengths": gene.get("product_lengths", []),
                **{k: u[k] for k in ("n_products", "product_lengths", "retained_path_records", "excluded_standard_equivalent_path_records", "exact_duplicate_path_records", "retention_truncated", "path_search_evaluated", "search_limits_reached", "limits")},
                "output_file": u.get("output_file"), "retention_totals": totals,
                "search_flags_any_threshold": {k: any(t[k] for t in u["thresholds"]) for k in SEARCH_FLAGS},
                "thresholds": u["thresholds"],
            })
            for row, product in zip(emitted, u["products"]):
                require(row["sequence"] == product["oriented_sequence"], "Uncertain sequence order differs")
                require(hashlib.sha256(row["sequence"].encode()).hexdigest() == product["sequence_sha256"], "Uncertain sequence hash differs")
                product_rows.append({
                    "accession": record["accession"], "gene": name, "header": row["header"],
                    "source_fasta": output["files"][u["output_file"]],
                    **{k: product[k] for k in ("id", "sequence_length", "sequence_sha256")},
                    "sources": [{k: v for k, v in source.items() if k not in ("marker_occurrences", "marker_covered_runs")} for source in product["sources"]],
                })
        require(standard == output["standard"] and uncertain == output["uncertain"], "Captured FASTA differs from raw files")
        run_rows.append({
            "order": len(run_rows) + 1, "accession": record["accession"], "taxon": sample["sample"],
            "role": record["role"], "standard_count": len(standard), "uncertain_count": len(uncertain),
            "wall_seconds": record["wall_seconds"], "maximum_rss_bytes": record["maximum_rss_bytes"],
            "maximum_rss_gib": record["maximum_rss_bytes"] / (1024 ** 3),
            "counts": {k: output["stats"][k] for k in COUNTS},
            "optional_counts": {k: output["stats"][k] for k in ("n_multi_kmers", "n_singleton_kmers") if k in output["stats"]},
            "stage_timings": output["stats"]["stage_timings"], "receipt": item["receipt"],
        })
    require(set(records) == {(a, r) for a in ACCESSIONS for r in ROLES}, "Sample/role set differs")
    samples = []
    for accession in ACCESSIONS:
        rows = [g for g in gene_rows if g["accession"] == accession]
        c = [x for x in assessment["comparisons"] if x["accession"] == accession]
        require(len(c) == 1 and c[0]["passed"], "Missing successful parity result")
        samples.append({**c[0], "taxon": frozen["samples"][accession]["sample"],
                        "expected_prefix": frozen["samples"][accession]["expected_prefix"],
                        "retained_path_records": sum(g["retained_path_records"] for g in rows),
                        "excluded_standard_equivalent_path_records": sum(g["excluded_standard_equivalent_path_records"] for g in rows),
                        "exact_duplicate_path_records": sum(g["exact_duplicate_path_records"] for g in rows),
                        "retention_totals": {k: sum(g["retention_totals"][k] for g in rows) for k in RETENTION_SUMS},
                        "retention_truncated_genes": [g["gene"] for g in rows if g["retention_truncated"]],
                        "search_limits_reached_genes": [g["gene"] for g in rows if g["search_limits_reached"]]})
    historical = []
    matched_fasta = []
    for loss in losses:
        seq = loss["sequence"]
        require(len(seq) == loss["length"] and hashlib.sha256(seq.encode()).hexdigest() == loss["sequence_sha256"], "Historical identity differs")
        matches = {}
        for role in ROLES:
            output = records[(loss["accession"], role)]["outputs"]
            matches[role] = {}
            for kind in ("standard", "uncertain"):
                rows = [r for r in output[kind] if r["gene"] == "insecta_" + loss["gene"]]
                matches[role][kind] = {"exact_oriented": [r["header"] for r in rows if r["sequence"] == seq],
                                       "reverse_complement_only": [r["header"] for r in rows if r["sequence"] != seq and r["sequence"] == rc(seq)]}
        checked = [r for r in assessment["historical_losses"] if (r["accession"], r["gene"], r["sequence_sha256"]) == (loss["accession"], loss["gene"], loss["sequence_sha256"])]
        require(len(checked) == 1 and checked[0]["matches"] == matches, "Historical assessment match differs")
        on = matches["changed_on"]
        categories = [f"{kind}_{orientation}" for kind in ("standard", "uncertain") for orientation in ("exact_oriented", "reverse_complement_only") if on[kind][orientation]]
        row = {k: v for k, v in loss.items() if k not in ("sequence", "historical_header")}
        row.update(classifications=categories or ["bounded_non_observation"], matches=matches,
                   non_observation_scope=checked[0]["non_observation_scope"],
                   gene_summary_key=[loss["accession"], "insecta_" + loss["gene"]])
        historical.append(row)
        for header in on["uncertain"]["exact_oriented"]:
            matched_fasta.append(f">{header}\n{seq}\n")
    return {
        "schema_version": 1, "status": "complete", "standard_parity": "pass",
        "method": "Post-run derivation from hash-bound receipts and their saved FASTA files; no FASTQ scan or new assembly.",
        "interpretation": assessment["interpretation"],
        "scope": {"run_count": 9, "panel_targets": 21, "k": 19, "threads": 2, "chunks": 0,
                  "reads_per_sample": 1000000, "read_threading": "disabled", "run_limits": frozen["limits"],
                  "uncertain_products_are_validated_haplotypes": False,
                  "non_observation": "Bounded search and retention only. Caps are static per-target shares of the full-panel budget; unused shares are not redistributed. Missing exact historical identities cannot be assigned to particular dropped paths.",
                  "performance": "One sequential trio per sample, with prefix verification before every run; descriptive observations only."},
        "samples": samples, "runs": run_rows, "uncertain_genes": gene_rows,
        "uncertain_products": product_rows, "historical_losses": historical,
        "historical_counts": {"total": len(historical), "uncertain_exact_oriented": sum(bool(r["matches"]["changed_on"]["uncertain"]["exact_oriented"]) for r in historical),
                              "uncertain_reverse_complement_only": sum(bool(r["matches"]["changed_on"]["uncertain"]["reverse_complement_only"]) for r in historical),
                              "standard_exact_or_rc": sum(any(r["matches"]["changed_on"]["standard"].values()) for r in historical),
                              "bounded_non_observation": sum(r["classifications"] == ["bounded_non_observation"] for r in historical)},
    }, "".join(matched_fasta)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = ("summary.json", "sample_runs.tsv", "uncertain_genes.tsv", "historical_losses.tsv", "historical_exact_uncertain.fasta", "derivation.json")
    require(not any((args.output_dir / name).exists() for name in names), "Refusing to overwrite derived artifacts")
    evidence = Evidence()
    source = evidence.bind(Path(__file__).absolute())
    summary, matched = derive(args.attempt.absolute(), evidence)
    evidence.finish()
    summary["source_evidence"] = evidence.sources
    payloads = {
        "summary.json": json.dumps(summary, indent=2, sort_keys=True) + "\n",
        "sample_runs.tsv": tsv(summary["runs"], ["order", "accession", "taxon", "role", "standard_count", "uncertain_count", "wall_seconds", "maximum_rss_bytes", "maximum_rss_gib"]),
        "uncertain_genes.tsv": tsv(summary["uncertain_genes"], ["accession", "taxon", "gene", "standard_count", "standard_lengths", "n_products", "product_lengths", "retained_path_records", "excluded_standard_equivalent_path_records", "exact_duplicate_path_records", "retention_truncated", "search_limits_reached", "retention_totals", "limits"]),
        "historical_losses.tsv": tsv(summary["historical_losses"], ["accession", "taxon", "gene", "length", "sequence_sha256", "classifications", "matches", "non_observation_scope"]),
        "historical_exact_uncertain.fasta": matched,
    }
    for name, payload in payloads.items():
        write_new(args.output_dir / name, payload)
    evidence.finish()
    receipt = {"schema_version": 1, "status": "complete", "source": source,
               "source_evidence_bookends_equal": True, "evidence_file_count": len(evidence.sources),
               "outputs": {name: fingerprint(args.output_dir / name) for name in payloads}}
    write_new(args.output_dir / "derivation.json", json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "complete", "historical_counts": summary["historical_counts"], "samples": [{k: row[k] for k in ("accession", "standard_counts", "uncertain_count", "retained_path_records", "excluded_standard_equivalent_path_records", "exact_duplicate_path_records", "retention_truncated_genes")} for row in summary["samples"]]}))


if __name__ == "__main__":
    main()

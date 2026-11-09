#!/usr/bin/env python3
"""Create descriptive tables from a completed, source-bound 60-run comparison.

Uncertain output denotes emitted hypotheses, not independent biological support.
Unique counts use (gene, oriented sequence); strand-invariant counts are explicit.
This posthoc reporter never invokes the assembler or reads benchmark FASTQs.
"""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import statistics


HERE = Path(__file__).resolve().parent
MAX_JSON_BYTES = 128 * 1024 * 1024
LABELS = ("standard", "uncertain")
WALL_SUM_ABS_TOLERANCE_SECONDS = 1e-9


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).absolute()
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode), f"Expected regular file: {path}")
    key = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        require(key(os.fstat(stream.fileno())) == key(before), "File replaced before hashing")
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        require(key(os.fstat(stream.fileno())) == key(before), "File changed while hashing")
    require(key(path.lstat()) == key(before), "File replaced after hashing")
    return {"path": str(path), "size_bytes": before.st_size, "sha256": digest.hexdigest()}


def read_json(path):
    require(Path(path).stat().st_size <= MAX_JSON_BYTES, "JSON artifact exceeds reporter limit")
    return json.loads(Path(path).read_text())


def reverse(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def check_record(row):
    sequence = row["sequence"]
    require(sequence and set(sequence) <= set("ACGT"), "Invalid captured sequence")
    require(row["assembly_status"] in LABELS, "Invalid captured output class")
    require(row["length"] == len(sequence), "Captured length differs")
    require(row["sha256"] == hashlib.sha256(sequence.encode()).hexdigest(), "Captured sequence SHA differs")
    require(row["reverse_complement_sha256"] == hashlib.sha256(reverse(sequence).encode()).hexdigest(),
            "Captured reverse-complement SHA differs")


def counts(records):
    result = {}
    for label in LABELS:
        selected = [row for row in records if row["assembly_status"] == label]
        result[label + "_records"] = len(selected)
        result[label + "_unique_oriented"] = len({(row["gene"], row["sequence"]) for row in selected})
    result["total_records"] = len(records)
    result["total_unique_oriented"] = len({(row["gene"], row["sequence"]) for row in records})
    result["total_unique_strand_invariant"] = len({(row["gene"], min(row["sequence"], reverse(row["sequence"])))
                                                 for row in records})
    return result


def identities(baseline, candidate):
    baseline_unique = {(row["gene"], row["sha256"]): row for row in baseline}
    dispositions = []
    for (gene, digest), original in sorted(baseline_unique.items()):
        evidence = {}
        for label in LABELS:
            for orientation, sequence in (("exact", original["sequence"]), ("rc", reverse(original["sequence"]))):
                evidence[label + "_" + orientation] = [row["sha256"] for row in candidate
                    if row["gene"] == gene and row["assembly_status"] == label and row["sequence"] == sequence]
        disposition = next((name for name in ("standard_exact", "standard_rc", "uncertain_exact", "uncertain_rc")
                            if evidence[name]), "unobserved")
        dispositions.append({"gene": gene, "sha256": digest, "length": original["length"],
                             "scope": original["scope"], "disposition": disposition, "evidence": evidence})
    gains = []
    for key, row in sorted({(row["gene"], row["sha256"], row["assembly_status"]): row for row in candidate}.items()):
        if key[:2] not in baseline_unique:
            gains.append({"gene": row["gene"], "sha256": row["sha256"], "length": row["length"],
                          "assembly_status": row["assembly_status"], "scope": row["scope"],
                          "reverse_complement_in_baseline": (row["gene"], row["reverse_complement_sha256"]) in baseline_unique})
    return dispositions, gains


def resources(before, after, priority_percent):
    require(len(before) == len(after) == 3, "Three paired resource measurements required")
    require(all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value >= 0
                for value in before + after), "Invalid resource observation")
    medians = [statistics.median(values) for values in (before, after)]
    percent = 100 * (medians[1] / medians[0] - 1) if medians[0] else None
    return {"baseline_min": min(before), "baseline_median": medians[0], "baseline_max": max(before),
            "candidate_min": min(after), "candidate_median": medians[1], "candidate_max": max(after),
            "median_absolute_delta": medians[1] - medians[0], "median_ratio": medians[1] / medians[0] if medians[0] else None,
            "median_percent_change": percent, "increase": medians[1] > medians[0],
            "priority_review": percent is not None and percent > priority_percent,
            "paired_ratios": [new / old if old else None for old, new in zip(before, after)]}


def historical_seven(oracle, dispositions):
    outcome = {(row["accession"], row["gene"], row["sha256"]): row for row in dispositions}
    lost = []
    for sample in oracle["sample_expectations"]:
        previous = {(gene, row["sequence_sha256"]) for gene, rows in sample["roles"]["candidate"]["by_gene"].items() for row in rows}
        for gene, rows in sample["roles"]["baseline"]["by_gene"].items():
            for row in rows:
                if row["scope"] != "high_copy_candidate" or (gene, row["sequence_sha256"]) in previous:
                    continue
                key = (sample["accession"], sample["panel"] + "_" + gene, row["sequence_sha256"])
                require(key in outcome, "Historical high-copy loss is absent from released identity table")
                lost.append(dict(outcome[key]))
    require(len(lost) == len({(row["accession"], row["gene"], row["sha256"]) for row in lost}) == 7,
            "Expected exactly seven unique historical high-copy losses")
    return lost


def wall_sum_matches(values, expected):
    """Allow only sub-nanosecond sum roundoff across Python runtime versions."""
    return math.isclose(math.fsum(values), expected, rel_tol=0.0,
                        abs_tol=WALL_SUM_ABS_TOLERANCE_SECONDS)


def derive(execution, summary, oracle):
    require(execution.get("status") == summary.get("status") == oracle.get("status") == "pass", "Incomplete input evidence")
    jobs = execution["jobs"]
    require(len(jobs) == 60 and all(job.get("status") == job["execution"].get("status") == "pass" for job in jobs),
            "Need all sixty successful invocations")
    by_sample = {}
    seen_ids = set()
    for job in jobs:
        invocation = job["invocation"]
        identifier = invocation["invocation_id"]
        require(identifier not in seen_ids, "Duplicate invocation identity")
        seen_ids.add(identifier)
        accession, role, pair = invocation["input"], invocation["version"], invocation["pair_index"]
        require(role in ("baseline", "candidate") and type(pair) is int and pair in (1, 2, 3), "Unexpected role/pair")
        entries = by_sample.setdefault(accession, {})
        require((role, pair) not in entries, "Duplicate sample/role/pair")
        entries[(role, pair)] = job
        for row in job["outputs"]["records"]:
            check_record(row)
    require(len(by_sample) == 10, "Expected ten samples")
    historical = {row["accession"]: row for row in oracle["sample_expectations"]}
    require(len(historical) == len(oracle["sample_expectations"]) == 10 and set(historical) == set(by_sample),
            "Historical oracle sample set differs")
    summary_cells = {cell["accession"]: cell for cell in summary["cells"]}
    require(len(summary_cells) == len(summary["cells"]) == 10 and set(summary_cells) == set(by_sample), "Summary sample set differs")
    tables = {name: [] for name in ("sample_counts", "released_identity_dispositions", "literal_gains", "resources", "paired_resources")}
    total_records = {role: Counter() for role in ("baseline", "candidate")}
    for accession, selected in by_sample.items():
        require(set(selected) == {(role, pair) for role in ("baseline", "candidate") for pair in (1, 2, 3)}, "Missing paired invocation")
        panel = selected[("baseline", 1)]["invocation"]["panel"]
        cell = summary_cells[accession]
        records_by_role = {}
        for role in ("baseline", "candidate"):
            replicates = [selected[(role, pair)] for pair in (1, 2, 3)]
            require(all(job["invocation"]["panel"] == panel for job in replicates), "Sample panel differs")
            for key in ("records", "counts", "pcr_results"):
                require(all(job["outputs"][key] == replicates[0]["outputs"][key] for job in replicates), "Replicate captured outputs differ")
            records = records_by_role[role] = replicates[0]["outputs"]["records"]
            for scope, chosen in (("all", records),
                                  ("high_copy", [row for row in records if row["scope"] == "high_copy_candidate"]),
                                  ("deferred_or_unclassified", [row for row in records if row["scope"] != "high_copy_candidate"])):
                values = counts(chosen)
                tables["sample_counts"].append({"accession": accession, "panel": panel, "version": role, "scope": scope, **values})
                if scope == "all":
                    require(cell["roles"][role]["record_counts"] == {label: values[label + "_records"] for label in LABELS}, "Summary record counts differ")
                    require(cell["roles"][role]["unique_counts"] == {label: values[label + "_unique_oriented"] for label in LABELS}, "Summary unique counts differ")
                    total_records[role].update({label: values[label + "_records"] for label in LABELS})
        require(selected[("baseline", 1)]["outputs"]["counts"] == selected[("candidate", 1)]["outputs"]["counts"], "Read/base/k-mer count parity differs")
        expected_release = {panel + "_" + gene: [row["sequence"] for row in rows]
                            for gene, rows in historical[accession]["roles"]["baseline"]["by_gene"].items()}
        observed_release = {gene: [] for gene in expected_release}
        for row in records_by_role["baseline"]:
            require(row["assembly_status"] == "standard" and row["gene"] in observed_release,
                    "Unexpected released output class or gene")
            observed_release[row["gene"]].append(row["sequence"])
        require(observed_release == expected_release, "Fresh released ordered identities differ from historical oracle")
        released, gains = identities(records_by_role["baseline"], records_by_role["candidate"])
        require(released == cell["released_identity_dispositions"], "Summary released-identity dispositions differ")
        key = lambda row: (row["gene"], row["sha256"], row["assembly_status"])
        require(gains == sorted(cell["literal_gains"], key=key), "Summary literal gains differ")
        tables["released_identity_dispositions"].extend({"accession": accession, "panel": panel, **row} for row in released)
        tables["literal_gains"].extend({"accession": accession, "panel": panel, **row} for row in gains)
        for metric, threshold in (("wall_seconds", 5), ("peak_rss_bytes", 2)):
            values = {role: [selected[(role, pair)]["execution"]["metrics"][metric] for pair in (1, 2, 3)]
                      for role in ("baseline", "candidate")}
            report = resources(values["baseline"], values["candidate"], threshold)
            for role in values:
                expected_metric = {"values": values[role], "median": statistics.median(values[role]),
                                   "min": min(values[role]), "max": max(values[role])}
                require(cell["roles"][role]["metrics"][metric] == expected_metric, "Summary resource observations differ")
            change = cell["resource_changes"][metric]
            require(change == {"absolute_delta": report["median_absolute_delta"], "percent": report["median_percent_change"],
                               "increase": report["increase"], "priority_review": report["priority_review"]}, "Summary resource flags differ")
            tables["resources"].append({"accession": accession, "panel": panel, "metric": metric,
                                        **{key: value for key, value in report.items() if key != "paired_ratios"}})
            for pair in (1, 2, 3):
                old, new = values["baseline"][pair - 1], values["candidate"][pair - 1]
                paired = {"baseline": old, "candidate": new, "absolute_delta": new - old, "ratio": new / old if old else None}
                require(cell["paired_resources"][pair - 1]["pair"] == pair and cell["paired_resources"][pair - 1][metric] == paired,
                        "Summary paired observations differ")
                tables["paired_resources"].append({"accession": accession, "panel": panel, "metric": metric, "pair": pair, **paired})
    require({role: dict(values) for role, values in total_records.items()} == summary["totals"], "Summary total counts differ")
    for role in ("baseline", "candidate"):
        require(wall_sum_matches([row[role + "_median"] for row in tables["resources"] if row["metric"] == "wall_seconds"],
                                 summary["median_wall_sum_seconds"][role]), "Summary sum of sample median walls differs")
    tables["historical_seven"] = historical_seven(oracle, tables["released_identity_dispositions"])
    return tables


def write_tables(tables, output):
    fields = {
        "sample_counts": ["accession", "panel", "version", "scope", "standard_records", "uncertain_records", "total_records",
                          "standard_unique_oriented", "uncertain_unique_oriented", "total_unique_oriented", "total_unique_strand_invariant"],
        "released_identity_dispositions": ["accession", "panel", "gene", "scope", "length", "sha256", "disposition", "evidence"],
        "literal_gains": ["accession", "panel", "gene", "scope", "assembly_status", "length", "sha256", "reverse_complement_in_baseline"],
        "resources": ["accession", "panel", "metric", "baseline_min", "baseline_median", "baseline_max", "candidate_min", "candidate_median",
                      "candidate_max", "median_absolute_delta", "median_ratio", "median_percent_change", "increase", "priority_review"],
        "paired_resources": ["accession", "panel", "metric", "pair", "baseline", "candidate", "absolute_delta", "ratio"],
    }
    fields["historical_seven"] = fields["released_identity_dispositions"]
    for name, rows in tables.items():
        with (output / (name + ".tsv")).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields[name], delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value for key, value in row.items()})


def generate(execution_path, summary_path, oracle_path, output):
    output.mkdir(parents=True, exist_ok=False)
    sources = {key: fingerprint(path) for key, path in (("execution", execution_path), ("summary", summary_path),
                                                        ("historical_identities", oracle_path), ("reporter", Path(__file__)))}
    receipt = {"status": "failed", "sources_before": sources}
    try:
        execution, summary, oracle = (read_json(path) for path in (execution_path, summary_path, oracle_path))
        require(execution["summary"] == sources["summary"], "Execution does not bind supplied summary")
        require(execution["bindings"].get(str(oracle_path.absolute())) == sources["historical_identities"], "Execution does not bind supplied historical oracle")
        for job in execution["jobs"]:
            binding = job["receipt"]
            require(fingerprint(binding["path"]) == binding, "Individual invocation receipt changed")
            require(read_json(binding["path"]) == {key: value for key, value in job.items() if key != "receipt"}, "Invocation receipt differs from execution")
            sources[binding["path"]] = binding
            for collection in (job["outputs"]["files"], job["execution"]["files"]):
                for row in collection.values():
                    require(fingerprint(row["path"]) == row, "Captured output/log changed")
                    sources[row["path"]] = row
        tables = derive(execution, summary, oracle)
        write_tables(tables, output)
        receipt.update(status="pass", numeric_comparison={"wall_sum_method": "math.fsum", "wall_sum_absolute_tolerance_seconds": WALL_SUM_ABS_TOLERANCE_SECONDS, "relative_tolerance": 0.0}, row_counts={key: len(rows) for key, rows in tables.items()},
                       record_totals=summary["totals"],
                       released_identity_dispositions=dict(Counter(row["disposition"] for row in tables["released_identity_dispositions"])),
                       historical_seven_dispositions=dict(Counter(row["disposition"] for row in tables["historical_seven"])),
                       median_wall_sum_seconds=summary["median_wall_sum_seconds"],
                       resource_increases=[row for row in tables["resources"] if row["increase"]],
                       interpretation="Descriptive identity/resource comparison. Uncertain records are emitted hypotheses with biological support unassessed; historical seven is a subset of released identities.",
                       files={path.name: fingerprint(path) for path in sorted(output.glob("*.tsv"))})
    except Exception as error:
        receipt.update(status="failed", error_type=type(error).__name__, error=str(error))
    try:
        receipt["sources_after"] = {key: fingerprint(row["path"]) for key, row in sources.items()}
        require(receipt["sources_after"] == sources, "Report sources changed during extraction")
    except Exception as error:
        receipt.update(status="failed", bookend_error=str(error))
    with (output / "report.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--historical-identities", type=Path, default=HERE.parent / "inputs/historical_identities.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = generate(args.execution, args.summary, args.historical_identities, args.output)
    print(json.dumps({"status": result["status"], "error": result.get("error"), "output": str(args.output)}))
    raise SystemExit(0 if result["status"] == "pass" else 1)

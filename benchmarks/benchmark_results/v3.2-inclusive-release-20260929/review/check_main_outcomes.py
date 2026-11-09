"""Independent actual artifact review of the frozen release comparison.

Reads existing artifacts and fingerprints inputs; never invokes Sharkmer.
Reuses the reviewed format readers for complete transaction validation, then
independently reconstructs identity, scope, historical and resource accounting.
"""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = Path(__file__).resolve().parents[1]
ATTEMPT = ROOT / "tmp/inclusive-release-20260929/comparison-v1"
sys.path.insert(0, str(BUNDLE))
import run_release_comparison as producer


def main():
    execution_path = ATTEMPT / "execution.json"
    execution = json.loads(execution_path.read_text())
    summary = json.loads((ATTEMPT / "summary.json").read_text())
    plan = json.loads((ATTEMPT / "plan.json").read_text())
    oracle = json.loads((BUNDLE / "inputs/historical_identities.json").read_text())
    inventory = json.loads((BUNDLE / "inputs/inventory.json").read_text())
    availability = json.loads((ROOT / "tmp/inclusive-release-20260929/availability.json").read_text())
    assert execution["status"] == summary["status"] == availability["status"] == "pass"
    assert len(execution["jobs"]) == len(plan["schedule"]) == 60
    assert [job["invocation"] for job in execution["jobs"]] == plan["schedule"]
    assert plan["bindings"] == execution["bindings"]
    cpus = execution["cpus"]
    assert len(cpus) == len(set(cpus)) == 2
    topology = {row["cpu"]: (row["socket"], row["core"]) for row in execution["allowed_topology"]}
    assert len({topology[cpu] for cpu in cpus}) == 2
    artifacts = {}

    def bind(row):
        if row["path"] in artifacts:
            assert artifacts[row["path"]] == row
            return
        assert producer.fingerprint(Path(row["path"])) == row, row["path"]
        artifacts[row["path"]] = row

    for row in execution["bindings"].values():
        bind(row)
    bind(execution["plan"])
    bind(execution["summary"])
    assert availability["sources_before"] == availability["sources_after"]
    for row in availability["sources_after"].values():
        bind(row)
    samples = {row["accession"]: row for row in availability["samples"]}
    assert len(samples) == 10
    original_samples = {row["accession"]: row for row in inventory["samples"]}
    for accession, sample in samples.items():
        assert sample["status"] == "verified"
        assert sample["expected_prefix"] == original_samples[accession]["expected_prefix"]
        bind(sample["receipt"])
        receipt = json.loads(Path(sample["receipt"]["path"]).read_text())
        assert receipt["observed_prefix"] == sample["expected_prefix"]
    scope = {panel: {panel + "_" + row["gene"]: row["scope"] for row in values}
             for panel, values in inventory["historical_protocol"]["target_metadata"].items()}
    expectations = {row["accession"]: row for row in oracle["sample_expectations"]}
    grouped = defaultdict(dict)
    for original, job in zip(inventory["historical_schedule"], execution["jobs"]):
        invocation = job["invocation"]
        assert all(invocation[key] == value for key, value in original.items())
        assert job["status"] == job["execution"]["status"] == "pass"
        assert job["execution"]["returncode"] == 0
        assert not job["execution"].get("timed_out")
        bind(job["receipt"])
        assert json.loads(Path(job["receipt"]["path"]).read_text()) == {key: value for key, value in job.items() if key != "receipt"}
        bind(job["input_file"])
        expected_prefix = samples[invocation["input"]]["expected_prefix"]
        assert all(job["input_file"][key] == expected_prefix[key] for key in ("sha256", "size_bytes"))
        for collection in (job["execution"]["files"], job["outputs"]["files"]):
            for row in collection.values():
                bind(row)
        timing = json.loads((Path(invocation["directory"]) / "time.json").read_text())
        assert timing["exit_code"] == 0
        metrics = job["execution"]["metrics"]
        assert all(metrics[key] == value for key, value in timing.items())
        assert metrics["peak_rss_bytes"] == timing["max_rss_kib"] * 1024
        child = json.loads((Path(invocation["directory"]) / "process_constraints.json").read_text())
        assert child["status"] == "pass" and child["affinity"] == cpus
        assert child["rlimit_as_soft_hard"] == [40 * 1024**3] * 2
        assert child["command"][0] == "/usr/bin/time" and child["command"][-len(invocation["command"]):] == invocation["command"]
        recaptured = producer.capture(Path(invocation["outputs"]), invocation, samples[invocation["input"]],
                                      invocation["command"], scope[invocation["panel"]])
        assert recaptured == job["outputs"]
        expected = {invocation["panel"] + "_" + gene: [row["sequence"] for row in rows]
                    for gene, rows in expectations[invocation["input"]]["roles"][invocation["version"]]["by_gene"].items()}
        observed = {gene: [] for gene in expected}
        for row in recaptured["records"]:
            assert hashlib.sha256(row["sequence"].encode()).hexdigest() == row["sha256"]
            assert len(row["sequence"]) == row["length"]
            if row["assembly_status"] == "standard":
                observed[row["gene"]].append(row["sequence"])
        assert observed == expected, (invocation["input"], invocation["version"])
        assert job["historical_identity_check"]["ordered_standard_sequences_identical"]
        grouped[invocation["input"]][(invocation["version"], invocation["pair_index"])] = job
    result_counts = {role: Counter() for role in ("baseline", "candidate")}
    scope_counts = {role: Counter() for role in result_counts}
    dispositions, gains, resources, cells, historical_seven = [], [], [], [], []
    for accession, jobs in grouped.items():
        assert len(jobs) == 6
        rows = {}
        for role in result_counts:
            first = jobs[(role, 1)]["outputs"]
            assert all(jobs[(role, pair)]["outputs"][key] == first[key] for pair in (2, 3) for key in ("records", "counts", "pcr_results"))
            rows[role] = first["records"]
            result_counts[role].update(row["assembly_status"] for row in rows[role])
            scope_counts[role].update((row["scope"], row["assembly_status"]) for row in rows[role])
        assert jobs[("baseline", 1)]["outputs"]["counts"] == jobs[("candidate", 1)]["outputs"]["counts"]
        baseline = {(row["gene"], row["sha256"]): row for row in rows["baseline"]}
        prior_standard = {(gene, product["sequence_sha256"])
                          for gene, products in expectations[accession]["roles"]["candidate"]["by_gene"].items() for product in products}
        panel = jobs[("baseline", 1)]["invocation"]["panel"]
        for key, old in baseline.items():
            candidate = {(row["assembly_status"], row["gene"], row["sha256"]) for row in rows["candidate"]}
            label = "unobserved"
            for status, orientation, digest in (("standard", "exact", old["sha256"]), ("standard", "rc", old["reverse_complement_sha256"]),
                                                 ("uncertain", "exact", old["sha256"]), ("uncertain", "rc", old["reverse_complement_sha256"])):
                if (status, old["gene"], digest) in candidate:
                    label = status + "_" + orientation
                    break
            found = {"accession": accession, "gene": old["gene"], "scope": old["scope"], "sha256": old["sha256"], "length": old["length"], "disposition": label}
            dispositions.append(found)
            if old["scope"] == "high_copy_candidate" and (old["gene"].removeprefix(panel + "_"), old["sha256"]) not in prior_standard:
                historical_seven.append(found)
        unique_candidate = {(row["gene"], row["sha256"], row["assembly_status"]): row for row in rows["candidate"]}
        gains.extend({"accession": accession, **row} for key, row in unique_candidate.items() if key[:2] not in baseline)
        for metric, priority in (("wall_seconds", 5), ("peak_rss_bytes", 2)):
            values = {role: [jobs[(role, pair)]["execution"]["metrics"][metric] for pair in (1, 2, 3)] for role in result_counts}
            medians = {role: statistics.median(values[role]) for role in result_counts}
            old, new = medians["baseline"], medians["candidate"]
            percent = 100 * (new / old - 1) if old else None
            resources.append({"accession": accession, "metric": metric, "values": values, "medians": medians, "absolute_delta": new - old,
                              "percent_change": percent, "increase": new > old, "priority_review": percent is not None and percent > priority})
        cells.append({"accession": accession, "counts": {role: dict(Counter(row["assembly_status"] for row in data)) for role, data in rows.items()}})
    assert result_counts["baseline"] == {"standard": 86}
    assert result_counts["candidate"] == {"standard": 64, "uncertain": 62}
    assert Counter(row["disposition"] for row in dispositions) == {"standard_exact": 64, "uncertain_exact": 8, "unobserved": 14}
    assert len(historical_seven) == 7 and Counter(row["disposition"] for row in historical_seven) == {"uncertain_exact": 5, "unobserved": 2}
    assert len(gains) == 54 and all(row["assembly_status"] == "uncertain" for row in gains)
    medians = {role: sum(row["medians"][role] for row in resources if row["metric"] == "wall_seconds") for role in result_counts}
    assert medians == summary["median_wall_sum_seconds"]
    for path in (Path(__file__), execution_path, BUNDLE / "review/main_execution_review.json", BUNDLE / "review/restoration_v2_review.json"):
        row = producer.fingerprint(path); artifacts[row["path"]] = row
    result = {"schema_version": 1, "status": "pass", "purpose": "independent_release_raw_outcome_review", "date": "2026-09-30",
              "reviewed_artifacts": artifacts, "invocations_recaptured": 60, "historical_baseline_and_candidate_standard_parity": True,
              "counts": {role: dict(count) for role, count in result_counts.items()},
              "scope_counts": {role: {scope + "/" + label: count for (scope, label), count in values.items()} for role, values in scope_counts.items()},
              "released_dispositions": dict(Counter(row["disposition"] for row in dispositions)), "released_identity_rows": dispositions,
              "historical_seven": historical_seven, "new_literal_candidate_identities": len(gains),
              "resources": resources, "median_wall_sum_seconds": medians, "cells": cells,
              "limitations": ["Artifact validity passes; performance regression disposition is separate and remains open.",
                              "All ten sample wall medians increased; four exceed the5% investigation threshold. No RSS increase exceeds2%.",
                              "Release-to-candidate differences cannot be attributed specifically to uncertainty emission from this comparison.",
                              "Uncertain records are bounded alternatives with biological correctness unassessed."]}
    path = BUNDLE / "review/raw_outcome_review.json"
    with path.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True); stream.write("\n")
    print(json.dumps({"status": "pass", "review": producer.fingerprint(path), "counts": result["counts"], "scope_counts": result["scope_counts"],
                      "historical_seven": historical_seven, "priorities": [row for row in resources if row["priority_review"]]}))


if __name__ == "__main__":
    main()

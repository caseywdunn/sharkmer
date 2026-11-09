#!/usr/bin/env python3
"""Independent table reconciliation and final artifact/document binding.

Reuses the completed independent raw-FASTA and scratch-ledger reviews; it does
not rescan biological inputs, rerun Sharkmer, or infer biological correctness.
"""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "tmp/inclusive-release-20260929"
DOCS = ["dev_docs/INCLUSIVE_RELEASE_PROTOCOL.md", "dev_docs/INCLUSIVE_RELEASE_RESULTS.md",
        "dev_docs/DROS12S_RETENTION_RESULTS.md", "dev_docs/HANDOFF.md", "dev_docs/PLAN.md",
        "dev_docs/REVIEW_v3.2.md", "CHANGELOG.md", "AGENTS.md"]


def load(path):
    return json.loads(Path(path).read_text())


def fp(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
            size += len(block)
    return {"path": str(path), "sha256": h.hexdigest(), "size_bytes": size}


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    raw_path = BUNDLE / "review/raw_outcome_review.json"
    dros_path = BUNDLE / "review/dros12s_outcome_review.json"
    resource_path = BUNDLE / "review/resource_assessment.json"
    report_path = BUNDLE / "results/report-v2/report.json"
    execution_path = SCRATCH / "comparison-v1/execution.json"
    raw, dros, resource, report, execution = map(load, [raw_path, dros_path, resource_path, report_path, execution_path])
    assert all(x["status"] == "pass" for x in [raw, dros, resource, report, execution])
    assert raw["invocations_recaptured"] == 60 and len(execution["jobs"]) == 60
    assert resource["priority_flags_resolved"] is False
    assert report["sources_before"] == report["sources_after"]
    inherited = {}
    for review in [raw, dros, resource]:
        for value in review["reviewed_artifacts"].values():
            old = inherited.setdefault(value["path"], value)
            assert old == value, ("conflicting review bindings", value["path"])
    assert inherited[str(execution_path.resolve())] == fp(execution_path)
    for value in report["sources_before"].values():
        if value["path"] in inherited:
            assert inherited[value["path"]] == value
    table = {}
    for name, value in report["files"].items():
        assert fp(value["path"]) == value
        table[name] = rows(Path(value["path"]))
    assert report["record_totals"] == {"baseline": {"standard": 86, "uncertain": 0}, "candidate": {"standard": 64, "uncertain": 62}}
    assert report["released_identity_dispositions"] == raw["released_dispositions"]
    assert report["historical_seven_dispositions"] == dict(Counter(r["disposition"] for r in raw["historical_seven"]))
    assert all(math.isclose(report["median_wall_sum_seconds"][k], raw["median_wall_sum_seconds"][k], rel_tol=0, abs_tol=1e-9) for k in ["baseline", "candidate"])
    for name, expected in [("released_identity_dispositions.tsv", raw["released_identity_rows"]), ("historical_seven.tsv", raw["historical_seven"])]:
        fields = ["accession", "gene", "scope", "length", "sha256", "disposition"]
        normalize = lambda r: tuple(str(r[k]) for k in fields)
        assert Counter(map(normalize, table[name])) == Counter(map(normalize, expected))
    captures = {}
    for job in execution["jobs"]:
        inv = job["invocation"]
        if inv["pair_index"] == 1:
            captures[(inv["input"], inv["version"])] = job["outputs"]["records"]
    assert len(captures) == 20
    assert len(table["sample_counts.tsv"]) == 60
    seen_counts = set()
    for row in table["sample_counts.tsv"]:
        key = (row["accession"], row["version"], row["scope"])
        assert key not in seen_counts
        seen_counts.add(key)
        scope = {"high_copy": "high_copy_candidate"}.get(key[2], key[2])
        records = [r for r in captures[key[:2]] if scope == "all" or r["scope"] == scope]
        for status in ["standard", "uncertain"]:
            selected = [r for r in records if r["assembly_status"] == status]
            assert int(row[status + "_records"]) == len(selected)
            assert int(row[status + "_unique_oriented"]) == len({(r["gene"], r["sha256"]) for r in selected})
        assert int(row["total_records"]) == len(records)
        assert int(row["total_unique_oriented"]) == len({(r["gene"], r["sha256"]) for r in records})
        assert int(row["total_unique_strand_invariant"]) == len({(r["gene"], min(r["sha256"], r["reverse_complement_sha256"])) for r in records})
    expected_gains = set()
    for (accession, version), records in captures.items():
        if version != "candidate":
            continue
        baseline = captures[(accession, "baseline")]
        prior = {(r["gene"], r["sha256"]) for r in baseline}
        for r in records:
            if (r["gene"], r["sha256"]) not in prior:
                assert (r["gene"], r["reverse_complement_sha256"]) not in prior
                expected_gains.add((accession, r["gene"], r["scope"], r["assembly_status"], str(r["length"]), r["sha256"], "False"))
    fields = ["accession", "gene", "scope", "assembly_status", "length", "sha256", "reverse_complement_in_baseline"]
    assert len(expected_gains) == 54
    assert Counter(tuple(r[k] for k in fields) for r in table["literal_gains.tsv"]) == Counter(expected_gains)
    independent_metrics = {(r["accession"], r["metric"]): r for r in raw["resources"]}
    assert len(table["resources.tsv"]) == len(independent_metrics) == 20
    for r in table["resources.tsv"]:
        old = independent_metrics[(r["accession"], r["metric"])]
        for version in ["baseline", "candidate"]:
            assert float(r[version + "_median"]) == old["medians"][version]
            assert float(r[version + "_min"]) == min(old["values"][version])
            assert float(r[version + "_max"]) == max(old["values"][version])
        assert float(r["median_absolute_delta"]) == old["absolute_delta"]
        assert float(r["median_percent_change"]) == old["percent_change"]
        assert r["increase"] == str(old["increase"])
        assert r["priority_review"] == str(old["priority_review"])
    assert len(table["paired_resources.tsv"]) == 60
    seen_pairs = set()
    for r in table["paired_resources.tsv"]:
        key = (r["accession"], r["metric"], int(r["pair"]))
        assert key not in seen_pairs
        seen_pairs.add(key)
        old = independent_metrics[key[:2]]
        baseline, candidate = (old["values"][v][key[2] - 1] for v in ["baseline", "candidate"])
        assert float(r["baseline"]) == baseline and float(r["candidate"]) == candidate
        assert float(r["absolute_delta"]) == candidate - baseline
        assert float(r["ratio"]) == candidate / baseline
    controls = load(BUNDLE / "results/package_controls.json")
    assert controls["status"] == "pass" and controls["sources_before"] == controls["sources_after"]
    for value in controls["sources_after"].values():
        assert fp(value["path"]) == value
    paths = [raw_path, dros_path, resource_path, report_path, execution_path,
             SCRATCH / "comparison-v1/summary.json", SCRATCH / "build-v1/builds.json", SCRATCH / "availability.json",
             SCRATCH / "dros12s/replay-v1/replay.receipt.json", SCRATCH / "dros12s/replay-v1/assessment.json",
             SCRATCH / "dros12s/build-v1/rust-controls.json", BUNDLE / "README.md",
             ROOT / "tmp/inclusive-output-20260929/attempt-v1/checks/checks.json"]
    paths += [ROOT / name for name in DOCS]
    paths += [SCRATCH / "build-v1" / role / name for role in ["baseline", "candidate"] for name in ["build.json", "snapshot.json"]]
    paths += [Path(value["path"]) for value in report["files"].values()]
    paths += list((BUNDLE / "review").glob("*.py")) + list((BUNDLE / "review").glob("*.json"))
    paths += list((BUNDLE / "results").glob("*.py")) + list((BUNDLE / "results").glob("*controls.json"))
    paths += [BUNDLE / "results/report-v1/report.json"]
    output = BUNDLE / "review/outcome_review.json"
    assert not output.exists(), "Final review is immutable; use a versioned review after a change"
    for path in paths:
        value = fp(path)
        inherited[value["path"]] = value
    final = {"schema_version": 1, "status": "pass", "date": "2026-09-30",
             "purpose": "Independent completed outcome, derived tables, and final document review; packaging clearance",
             "review_basis": "Previously completed independent 60-invocation raw FASTA/stats/provenance recapture and 16-ledger scratch review, followed by independent reconciliation of all six report-v2 tables. No biological reruns.",
             "counts": raw["counts"], "scope_counts": raw["scope_counts"], "released_dispositions": raw["released_dispositions"],
             "historical_seven": raw["historical_seven"], "new_literal_candidate_identities": 54,
             "report_row_counts": report["row_counts"], "median_wall_sum_seconds": report["median_wall_sum_seconds"],
             "runtime_priority_flags": [r for r in raw["resources"] if r["priority_review"]],
             "runtime_priority_flags_resolved": False, "release_approval": False,
             "limitations": ["Uncertain products are bounded hypotheses; biological correctness and full-product read support remain unassessed.",
                             "Exact 487/499 bp Drosophila 12S completions are discarded by existing retention; production output and caps remain unchanged.",
                             "All ten median wall times increase, four exceed 5%; no peak-RSS median increase exceeds 2%. Three descriptive replicates cannot isolate inclusive emission overhead.",
                             "Reporter-v1's preserved float-summation failure was resolved in versioned v2 with aggregate-only 1e-9-second absolute tolerance.",
                             "Archive generation and independent member verification follow this receipt; this receipt does not claim those future checks already passed."],
             "reviewed_artifacts": dict(sorted(inherited.items()))}
    with output.open("x") as handle:
        json.dump(final, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": "pass", "review": fp(output), "reviewed_artifact_count": len(inherited)}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Read-only assessment of saved comparison-v1 evidence; never runs Sharkmer."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics

import yaml


ROOT = Path(__file__).resolve().parents[4]
COMPARISON = ROOT / "tmp/inclusive-release-20260929/comparison-v1"
BUILD = ROOT / "tmp/inclusive-release-20260929/build-v1"
FLAGS = {"SRR8617500", "SRR17535371", "SRR25378184", "SRR1057608"}


def assess():
    artifacts = {}

    def read(path, expected=None):
        path = Path(path).resolve()
        data = path.read_bytes()
        row = {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
               "size_bytes": len(data)}
        if expected is not None:
            assert row == expected, str(path)
        key = str(path.relative_to(ROOT))
        assert key not in artifacts or artifacts[key] == row
        artifacts[key] = row
        return data

    def load(path, expected=None):
        return json.loads(read(path, expected))

    read(__file__)
    read(ROOT / "dev_docs/INCLUSIVE_RELEASE_PROTOCOL.md")
    summary = load(COMPARISON / "summary.json")
    execution = load(COMPARISON / "execution.json")
    load(COMPARISON / "plan.json")
    assert summary["status"] == execution["status"] == "pass"
    assert len(execution["jobs"]) == 60
    provenance = {}
    for role in ("baseline", "candidate"):
        snapshot = load(BUILD / role / "snapshot.json")
        build = load(BUILD / role / "build.json")
        assert build["status"] == "complete" and build["returncode"] == 0
        read(build["binary"]["path"], build["binary"])
        for name in ("src/main.rs", "src/stats.rs", "src/pcr/mod.rs", "src/pcr/paths.rs"):
            read(BUILD / role / "source" / name, snapshot["files"][name])
        if role == "candidate":
            read(BUILD / role / "source/src/output.rs", snapshot["files"]["src/output.rs"])
        provenance[role] = {"commit": snapshot["provenance"]["commit"],
                            "binary": build["binary"], "command": build["command"]}

    groups = defaultdict(dict)
    historical_checks = 0
    selected_runs = defaultdict(lambda: defaultdict(list))
    for embedded in execution["jobs"]:
        job = load(embedded["receipt"]["path"], embedded["receipt"])
        assert {k: v for k, v in embedded.items() if k != "receipt"} == job
        inv = job["invocation"]
        accession, role, pair = inv["input"], inv["version"], inv["pair_index"]
        assert job["status"] == job["execution"]["status"] == "pass"
        assert job["execution"]["returncode"] == 0
        hist = job["historical_identity_check"]
        assert hist["ordered_standard_sequences_identical"] is True
        assert not hist["gained_literal_identities"] and not hist["missing_literal_identities"]
        historical_checks += role == "candidate"
        metrics = job["execution"]["metrics"]
        files = job["execution"]["files"]
        raw = load(files["time.json"]["path"], files["time.json"])
        assert all(raw[k] == metrics[k] for k in raw)
        assert metrics["peak_rss_bytes"] == raw["max_rss_kib"] * 1024
        constraints = load(files["process_constraints.json"]["path"], files["process_constraints.json"])
        assert constraints["status"] == "pass" and constraints["affinity"] == [0, 2]
        assert constraints["rlimit_as_soft_hard"] == [42949672960, 42949672960]
        assert constraints["command"][-len(inv["command"]):] == inv["command"]
        assert "--read-threading" not in inv["command"]
        assert accession not in groups or (role, pair) not in groups[accession]
        groups[accession][role, pair] = job
        if accession not in FLAGS:
            continue
        for name in ("stdout.log", "stderr.log"):
            read(files[name]["path"], files[name])
        stats_row = next(row for name, row in job["outputs"]["files"].items()
                         if name.endswith(".stats.yaml"))
        stats = yaml.safe_load(read(stats_row["path"], stats_row))
        assert all(stats[k] == v for k, v in job["outputs"]["counts"].items())
        stages = stats.get("stage_timings")
        if stages is not None:
            assert stages["read_threading_input_s"] is None
        selected_runs[accession][role].append({
            "pair": pair, "receipt": str(Path(embedded["receipt"]["path"]).relative_to(ROOT)),
            "wall_seconds": raw["wall_seconds"], "user_seconds": raw["user_seconds"],
            "system_seconds": raw["system_seconds"],
            "cpu_seconds": raw["user_seconds"] + raw["system_seconds"],
            "peak_rss_bytes": metrics["peak_rss_bytes"],
            "allocator_peak_bytes": job["outputs"]["allocator_peak_bytes"],
            "stage_timings": stages,
            "wall_minus_pipeline_seconds": None if stages is None else
                raw["wall_seconds"] - stages["pipeline_total_s"],
            "stats_size_bytes": stats_row["size_bytes"],
            "fasta_files": sum(name.endswith(".fasta") for name in job["outputs"]["files"]),
        })

    scopes = {role: defaultdict(Counter) for role in ("baseline", "candidate")}
    dispositions, scoped_dispositions, gains = Counter(), defaultdict(Counter), Counter()
    high_copy_losses, cells, gain_cells = [], [], []
    for cell in summary["cells"]:
        acc = cell["accession"]
        jobs = groups[acc]
        assert len(jobs) == 6
        expected_prefix = jobs["baseline", 1]["prefix_before"]
        expected_counts = jobs["baseline", 1]["outputs"]["counts"]
        for job in jobs.values():
            assert job["prefix_before"] == expected_prefix
            assert job["outputs"]["counts"] == expected_counts
        for role, recorded in cell["roles"].items():
            for scope, counts in recorded["scope_unique_counts"].items():
                scopes[role][scope].update(counts)
            for metric, threshold in (("wall_seconds", 5), ("peak_rss_bytes", 2)):
                values = [jobs[role, i]["execution"]["metrics"][metric] for i in (1, 2, 3)]
                assert recorded["metrics"][metric]["values"] == values
                assert recorded["metrics"][metric]["median"] == statistics.median(values)
                base = cell["roles"]["baseline"]["metrics"][metric]["median"]
                new = cell["roles"]["candidate"]["metrics"][metric]["median"]
                change = cell["resource_changes"][metric]
                assert math.isclose(change["percent"], 100 * (new - base) / base)
                assert change["priority_review"] == (change["percent"] > threshold)
        for row in cell["released_identity_dispositions"]:
            dispositions[row["disposition"]] += 1
            scoped_dispositions[row["scope"]][row["disposition"]] += 1
            if row["scope"] == "high_copy_candidate" and row["disposition"] != "standard_exact":
                high_copy_losses.append({"accession": acc, **{k: row[k] for k in
                                        ("gene", "length", "sha256", "disposition")}})
        for row in cell["literal_gains"]:
            gains[row["assembly_status"]] += 1
            assert row["scope"] == "high_copy_candidate"
        gain_cells.append({"accession": acc, "literal_gains": len(cell["literal_gains"]),
                           "classes": dict(Counter(row["assembly_status"] for row in cell["literal_gains"]))})
        cells.append({"accession": acc, "taxon": jobs["candidate", 1]["invocation"]["taxon"],
                      "counts": {role: cell["roles"][role]["unique_counts"] for role in cell["roles"]},
                      "resource_changes": cell["resource_changes"],
                      "wall_seconds": {role: cell["roles"][role]["metrics"]["wall_seconds"]
                                       for role in cell["roles"]}})
    assert {c["accession"] for c in cells if c["resource_changes"]["wall_seconds"]["priority_review"]} == FLAGS
    assert not any(c["resource_changes"]["peak_rss_bytes"]["priority_review"] for c in cells)
    assert sum(dispositions.values()) == 86 and historical_checks == 30
    for acc, roles in selected_runs.items():
        for role, rows in roles.items():
            rows.sort(key=lambda row: row["pair"])

    report = {
        "schema_version": 1, "status": "pass", "review_type": "saved_evidence_resource_assessment",
        "scope": "Read-only numerical/source assessment; no rerun, new assay, production edit, or causal experiment.",
        "priority_flags_resolved": False, "build_provenance": provenance,
        "conditions": {"completed_cli_invocations": 60, "paired_repeats_per_cell": 3,
                       "affinity": [0, 2], "address_space_limit_bytes": 42949672960,
                       "read_threading_enabled": False,
                       "human_input_records": groups["SRR17535371"]["baseline", 1]["prefix_before"]["records"]},
        "sequence_arithmetic": {"totals": summary["totals"], "scope_totals": scopes,
                                "released_dispositions": dispositions, "scope_dispositions": scoped_dispositions,
                                "high_copy_standard_losses": high_copy_losses,
                                "new_literal_identities": gains, "literal_gains_by_cell": gain_cells,
                                "candidate_historical_ordered_standard_checks_passed": historical_checks,
                                "biological_support": "unassessed"},
        "median_wall_sum_seconds": summary["median_wall_sum_seconds"],
        "cells": cells, "flagged_sample_saved_runs": selected_runs,
        "findings": [
            "The four >5% median-wall flags are observed regressions in this batch and remain open. No median-RSS increase exceeds the protocol's 2% priority trigger. These are review triggers, not significance tests or acceptance allowances.",
            "Human median wall increases 1.55 to 4.29 seconds (+2.74 seconds, +176.77%), while summed user+system CPU median increases 0.65 to 0.71 seconds. All five output sequences remain standard; no uncertain sequence is emitted. The complete saved input has 108518 reads, reaching EOF below the one-million-read cap. A small denominator enlarges its percentage but does not erase the absolute slowdown.",
            "Rhopilema also emits zero uncertainty records. Its CPU median increases 64.36 to 64.81 seconds versus wall 43.87 to 46.53 seconds. Therefore uncertain-sequence serialization cannot explain either zero-uncertainty cell on its own.",
            "Liriodendron emits six standard plus 27 uncertain records; CPU median increases 87.43 to 88.79 seconds and wall 52.63 to 56.80 seconds. Heliconius emits six standard plus 13 uncertain records; CPU median increases 71.66 to 72.28 seconds and wall 41.59 to 44.07 seconds. Additional retained payload and search changes are possible costs, not isolated causes.",
            "Candidate stats::run_pcr includes parallel assembly followed by sequential FASTA staging. main.rs times this entire call as pcr_s, so pcr_s is not an assembly-only measurement. pipeline_total_s is captured before stats serialization, transaction commit and reporting. GNU wall minus pipeline_total_s also includes process startup and cleanup and timing quantization; it is not a measured fsync duration.",
            "Candidate output.rs flushes and syncs staged files, synchronizes manifests/directories, hashes output files and transactionally publishes them. Baseline stats.rs writes FASTA files directly. This provides a plausible source-level contributor to wall time outside CPU-heavy assembly, but saved evidence has no syscall/storage/scheduler profiling to establish its contribution.",
            "This compares released v3.1.0 with ed04bde, spanning threshold/admission/repeat-diagnostic and output-publication changes. It is not an inclusive-output on/off ablation, and read threading is disabled. Attribute the measured change to the compared builds, not uniquely to inclusive delivery.",
            "Three paired same-host repetitions show variability: Rhopilema pair 2 is +0.27% and Heliconius pair 3 is +0.047%; the unflagged ERR14009273 candidate has a 47.19-second first repeat versus 25.30/25.28 seconds later. Saved data cannot establish portability or stable causal effect sizes.",
            "CPU seconds sum user and system time across threads and can exceed wall time; subtracting them from wall is not a blocked-time measurement. Peak RSS can be dominated by the count table, so a non-priority RSS result does not establish zero memory overhead.",
            "All54 additional literal identities are high-copy uncertainty records. The seven high-copy standard losses are five exact uncertainty re-emissions and two Drosophila12S identities absent from these normal outputs. Exact identity and greater output count do not validate repeat length, haplotypes or biological correctness. The separate retention observer is outside this resource comparison.",
        ],
        "source_locations": {
            "candidate/src/main.rs": [91, 202, 215, 245, 250, 251],
            "candidate/src/stats.rs": [127, 148, 192],
            "candidate/src/output.rs": [152, 166, 191, 469, 486, 497, 545],
            "baseline/src/stats.rs": [104, 108],
        },
        "reviewed_artifacts": artifacts,
    }
    # Check every consumed artifact again after parsing and calculations.
    for row in list(artifacts.values()):
        read(row["path"], row)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = assess()
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "output": str(args.output),
                      "bound_artifacts": len(result["reviewed_artifacts"]),
                      "open_wall_priority_flags": sorted(FLAGS)}))

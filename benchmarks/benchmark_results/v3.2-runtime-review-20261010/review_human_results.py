#!/usr/bin/env python3
"""Independent receipt, literal-output, timing and saved-trace review; no CLI runs."""
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
import statistics

import yaml

HERE = Path(__file__).resolve().parent
receipts = {}


def bind(path, expected=None):
    path = Path(path).resolve()
    raw = path.read_bytes()
    seen = {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}
    if expected is not None:
        assert seen == {k: expected[k] for k in seen}, str(path)
    receipts[str(path)] = seen
    return raw


def read(path):
    return json.loads(bind(path))


def main():
    protocol = read(HERE / "protocol-v1.json")
    result = read(HERE / "results-v1/result.json")
    analysis = read(HERE / "analysis-v1/analysis.json")
    for source in protocol["sources"]:
        bind(source["path"], source)
    bind(result["protocol"]["path"], result["protocol"])
    bind(HERE / "analyze_human_results.py", analysis["producer"])
    for path, receipt in analysis["source_files"].items():
        bind(path, receipt)
    for name, receipt in analysis["outputs"].items():
        bind(HERE / "analysis-v1" / name, receipt)
    accounting = read(HERE / "slurm-29012047-accounting.json")
    bind(accounting["slurm_log"]["path"], accounting["slurm_log"])
    assert "29012047|COMPLETED|0:0|00:01:17|2|" in accounting["output"]
    assert result["status"] == analysis["status"] == "pass"
    assert result["slurm_job_id"] == analysis["slurm_job_id"] == "29012047"
    assert result["cpus"] == [32, 36]
    assert len({(r["socket"], r["core"]) for r in result["topology"] if r["cpu"] in result["cpus"]}) == 2
    assert result["mounts"]["nfs_output"]["filesystems"][0]["fstype"] == "nfs"
    assert result["mounts"]["local_output"]["filesystems"][0]["fstype"] == "ext4"
    assert [r["schedule"] for r in result["runs"]] == protocol["schedule"]
    observations, stage_medians = [], {}
    rows = list(csv.DictReader((HERE / "analysis-v1/runs.tsv").open(), delimiter="\t"))
    assert len(rows) == len(result["runs"]) == 16
    for run, tab in zip(result["runs"], rows):
        item = run["schedule"]
        saved = HERE / "results-v1" / item["id"]
        assert read(saved / "receipt.json") == run
        assert run["status"] == "pass" and run["returncode"] == 0
        assert run["input_preflight"] == protocol["expected_input"]
        assert run["command"][0] == protocol["binaries"][item["role"]]["path"]
        assert run["command"][-1] == protocol["input"]["path"]
        assert run["command"][run["command"].index("--pcr-panel-file") + 1] == protocol["panel"]["path"]
        time = read(saved / "local_artifacts/time.json")
        assert all(run["metrics"][k] == v for k, v in time.items())
        assert run["metrics"]["cpu_seconds"] == time["user_seconds"] + time["system_seconds"]
        constraints = read(saved / "local_artifacts/process_constraints.json")
        assert constraints["status"] == "pass" and constraints["affinity"] == result["cpus"]
        assert constraints["rlimit_as_soft_hard"] == [12 * 1024**3] * 2
        assert ("/usr/bin/strace" in constraints["command"]) == item["traced"]
        for source in run["preserved_local_files"]:
            bind(source["path"], source)
        output = saved / ("local_artifacts/outputs" if item["location"] == "local" else "outputs")
        for source in run["outputs"]["files"]:
            bind(output / Path(source["path"]).relative_to(run["output_path_during_timing"]), source)
        stats = yaml.safe_load(bind(next(output.glob("*.stats.yaml"))))
        assert stats["command"] == " ".join(run["command"])
        expected = protocol["expected_outputs"][item["role"]]
        assert {k: stats[k] for k in expected["counts"]} == expected["counts"] == run["outputs"]["counts"]
        assert sorted([r["gene_name"], r["status"], r["n_products"]] for r in stats["pcr_results"]) == expected["targets"]
        records = []
        for fasta in sorted(output.glob("*.fasta")):
            gene = fasta.stem.removeprefix(stats["sample"] + "_")
            for record in fasta.read_text().split(">")[1:]:
                header, *parts = record.splitlines()
                sequence = "".join(parts)
                status = re.findall(r"(?:^|\s)assembly_status=(\S+)", header)
                assert status == (["standard"] if item["role"] == "retention" else [])
                records.append({"gene": gene, "assembly_status": "standard", "length": len(sequence), "sha256": hashlib.sha256(sequence.encode()).hexdigest()})
        assert sorted(records, key=lambda r: (r["gene"], r["sha256"])) == expected["records"] == run["outputs"]["records"]
        assert stats.get("stage_timings") == run["outputs"]["stage_timings"]
        assert tab["run_id"] == item["id"] and tab["traced"] == str(item["traced"])
        for key, value in run["metrics"].items():
            assert float(tab[key]) == value
        if item["traced"]:
            for trace in sorted((saved / "local_artifacts").glob("trace.*")):
                for number, line in enumerate(trace.read_text().splitlines(), 1):
                    # This saved batch has no unfinished/resumed calls. Reject them rather than silently omit.
                    assert "unfinished ..." not in line and " resumed>" not in line
                    match = re.fullmatch(r"\d+\.\d+ ([A-Za-z0-9_]+)\(.* <([0-9.]+)>", line)
                    if not match:
                        assert "+++ exited with 0 +++" in line, line
                        continue
                    path = run["output_path_during_timing"]
                    # Independently require a directory boundary rather than a raw substring.
                    assigned = bool(re.search(re.escape(path) + r'(?=[/">])', line))
                    observations.append({"run_id": item["id"], "location": item["location"], "role": item["role"],
                        "thread_trace": trace.name, "line": number, "syscall": match[1],
                        "scope": "product_output" if assigned else "other_or_unresolved", "seconds": float(match[2]),
                        "error": bool(re.search(r"= -1 [A-Z]", line)), "raw": line})
    assert observations == analysis["trace_observations"]
    assert analysis["unparsed_syscall_lines"] == []
    groups = collections.defaultdict(list)
    for row in observations:
        groups[(row["run_id"], row["syscall"], row["scope"])].append(row)
    assert len(groups) == len(analysis["trace_summaries"])
    for row in analysis["trace_summaries"]:
        selected = groups[(row["run_id"], row["syscall"], row["scope"])]
        times = [r["seconds"] for r in selected]
        assert row["calls"] == len(times) and row["sum_call_seconds"] == sum(times)
        assert row["max_call_seconds"] == max(times) and row["median_call_seconds"] == statistics.median(times)
        assert row["errors"] == sum(r["error"] for r in selected)
    normal = [r for r in result["runs"] if not r["schedule"]["traced"]]
    assert len(normal) == analysis["normal_runs"] == 12 and analysis["separate_traced_runs"] == 4
    for row in analysis["normal_medians"]:
        by_role = {}
        for role in ("release", "retention"):
            selected = [r for r in normal if r["schedule"]["role"] == role and r["schedule"]["location"] == row["location"]]
            assert len(selected) == 3
            by_role[role] = statistics.median(r["metrics"]["wall_seconds"] for r in selected)
            assert row[role + "_wall_s"] == by_role[role]
            assert row[role + "_cpu_s"] == statistics.median(r["metrics"]["cpu_seconds"] for r in selected)
            assert row[role + "_peak_rss_kib"] == statistics.median(r["metrics"]["max_rss_kib"] for r in selected)
            if role == "retention":
                stage_medians[row["location"]] = {k: statistics.median(r["outputs"]["stage_timings"][k] for r in selected)
                    for k in ("read_ingest_s", "pcr_s", "pipeline_total_s")}
                stage_medians[row["location"]]["wall_minus_pipeline_s"] = statistics.median(r["metrics"]["wall_seconds"] - r["outputs"]["stage_timings"]["pipeline_total_s"] for r in selected)
        assert row["wall_absolute_delta_s"] == by_role["retention"] - by_role["release"]
        assert row["wall_percent_delta"] == 100 * (by_role["retention"] / by_role["release"] - 1)
    bind(HERE / "HUMAN_FINDINGS.md")
    bind(Path(__file__))
    review = {"status": "pass", "scope": "Independent saved-artifact review; no further assembly or FASTQ assay", "runs": 16,
        "normal_runs": 12, "separate_traced_runs": 4, "literal_products_verified_per_run": 5, "target_outcomes_verified_per_run": 9,
        "raw_trace_calls_independently_parsed": len(observations), "unparsed": 0, "normal_medians": analysis["normal_medians"],
        "retention_stage_medians": stage_medians, "source_receipts": receipts,
        "interpretation": ["Observed output-location sensitivity in this bounded Human batch, with three normal repeats per cell.",
            "Decoded output-path syscall attribution agrees with stricter directory-boundary matching for every saved call.",
            "Current traces support metadata/publication latency; they do not support fsync itself as the dominant delay.",
            "Trace sums overlap across threads and include instrumentation effects; they cannot decompose untraced wall time.",
            "The frozen binary comparison does not identify historical causation or establish final corrected-build performance acceptance.",
            "The other three historical runtime findings and collaborator acceptance remain open; no release authorization."]}
    with (HERE / "independent_human_results_review.json").open("x") as stream:
        json.dump(review, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": "pass", "runs": 16, "trace_calls": len(observations), "source_receipts": len(receipts)}))


if __name__ == "__main__":
    main()

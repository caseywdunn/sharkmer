#!/usr/bin/env python3
"""Summarize normal timings and separate publication syscall diagnostics."""
import argparse
import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import statistics

HERE = Path(__file__).resolve().parent
CALL = re.compile(r"^\d+\.\d+\s+([a-zA-Z0-9_]+)\(")
DURATION = re.compile(r"\s<([0-9.]+)>$")
RESUMED = re.compile(r"^\d+\.\d+\s+<\.\.\. ([a-zA-Z0-9_]+) resumed>")


def fingerprint(path):
    raw = Path(path).read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}


def make_tsv(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def traces_for(run, saved):
    summaries, observations, unparsed = [], [], []
    for path in sorted((saved / "local_artifacts").glob("trace.*")):
        pending = None
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            match = CALL.match(line)
            resumed = RESUMED.match(line)
            if match and line.endswith("<unfinished ...>"):
                if pending is not None:
                    raise ValueError("Nested unfinished syscall in one thread")
                pending = (match[1], line, line_number)
                continue
            if resumed:
                if pending is None or resumed[1] != pending[0]:
                    raise ValueError("Resumed syscall has no matching beginning")
                syscall, first, start_line = pending
                full = first + " " + line
                pending = None
            elif match:
                syscall, full, start_line = match[1], line, line_number
            else:
                continue  # signal and process-exit records are retained in raw logs
            duration = DURATION.search(line)
            if duration is None:
                unparsed.append({"path": str(path), "line": start_line, "text": full})
                continue
            product_path = run["output_path_during_timing"]
            scope = "product_output" if product_path in full else "other_or_unresolved"
            # -yy resolves fd-based synchronization and writes. Relative names
            # without a decoded fd are deliberately not attributed to outputs.
            observations.append({"run_id": run["schedule"]["id"], "location": run["schedule"]["location"],
                "role": run["schedule"]["role"], "thread_trace": path.name, "line": start_line,
                "syscall": syscall, "scope": scope, "seconds": float(duration[1]), "error": bool(re.search(r"= -1 [A-Z]", full)),
                "raw": full})
        if pending is not None:
            unparsed.append({"path": str(path), "line": pending[2], "text": pending[1]})
    groups = collections.defaultdict(list)
    for observation in observations:
        groups[(observation["syscall"], observation["scope"])].append(observation)
    for (syscall, scope), rows in sorted(groups.items()):
        values = [row["seconds"] for row in rows]
        summaries.append({"run_id": run["schedule"]["id"], "location": run["schedule"]["location"], "role": run["schedule"]["role"],
            "syscall": syscall, "scope": scope, "calls": len(values), "sum_call_seconds": sum(values),
            "max_call_seconds": max(values), "median_call_seconds": statistics.median(values), "errors": sum(row["error"] for row in rows)})
    return summaries, observations, unparsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=HERE / "results-v1")
    parser.add_argument("--output", type=Path, default=HERE / "analysis-v1")
    args = parser.parse_args()
    results = args.results.resolve()
    source_files = {}

    def verify(path, expected=None):
        seen = fingerprint(path)
        if expected and seen != {k: expected[k] for k in ("sha256", "size_bytes")}:
            raise ValueError("Preserved source fingerprint differs: " + str(path))
        source_files[str(path)] = seen

    verify(results / "result.json")
    result = json.loads((results / "result.json").read_text())
    if result["status"] != "pass" or len(result["runs"]) != 16:
        raise ValueError("Need all sixteen completed invocations")
    verify(result["protocol"]["path"], result["protocol"])
    rows, trace_summaries, observations, unparsed = [], [], [], []
    for run in result["runs"]:
        item = run["schedule"]
        saved = results / item["id"]
        verify(saved / "receipt.json")
        if json.loads((saved / "receipt.json").read_text()) != run:
            raise ValueError("Per-run receipt differs from batch result")
        for file in run["preserved_local_files"]:
            verify(file["path"], file)
        for file in run["outputs"]["files"]:
            target = Path(file["path"])
            if item["location"] == "local":
                target = saved / "local_artifacts" / "outputs" / target.relative_to(run["output_path_during_timing"])
            verify(target, file)
        stages = run["outputs"]["stage_timings"] or {}
        pipeline = stages.get("pipeline_total_s")
        rows.append({"run_id": item["id"], "block": item["block"], "location": item["location"], "role": item["role"],
            "traced": item["traced"], **run["metrics"], "read_ingest_s": stages.get("read_ingest_s"),
            "count_finalize_s": stages.get("count_finalize_s"), "pcr_including_fasta_staging_s": stages.get("pcr_s"),
            "pipeline_total_s": pipeline, "wall_minus_pipeline_s": run["metrics"]["wall_seconds"] - pipeline if pipeline is not None else None,
            "product_records": len(run["outputs"]["records"]), "input_records": run["outputs"]["counts"]["n_reads_read"],
            "allocator_peak_bytes": run["outputs"]["allocator_peak_bytes"]})
        if item["traced"]:
            summary, detail, remaining = traces_for(run, saved)
            trace_summaries.extend(summary)
            observations.extend(detail)
            unparsed.extend(remaining)
    medians = []
    for location in ("nfs", "local"):
        by_role = {}
        for role in ("release", "retention"):
            selected = [r for r in rows if not r["traced"] and r["location"] == location and r["role"] == role]
            if len(selected) != 3:
                raise ValueError("Incomplete normal cell")
            by_role[role] = {key: statistics.median(r[key] for r in selected) for key in ("wall_seconds", "cpu_seconds", "max_rss_kib")}
        medians.append({"location": location, "release_wall_s": by_role["release"]["wall_seconds"],
            "retention_wall_s": by_role["retention"]["wall_seconds"],
            "wall_absolute_delta_s": by_role["retention"]["wall_seconds"] - by_role["release"]["wall_seconds"],
            "wall_percent_delta": 100 * (by_role["retention"]["wall_seconds"] / by_role["release"]["wall_seconds"] - 1),
            "release_cpu_s": by_role["release"]["cpu_seconds"], "retention_cpu_s": by_role["retention"]["cpu_seconds"],
            "release_peak_rss_kib": by_role["release"]["max_rss_kib"], "retention_peak_rss_kib": by_role["retention"]["max_rss_kib"]})
    args.output.mkdir(parents=True, exist_ok=False)
    payloads = {"runs.tsv": make_tsv(rows), "normal_medians.tsv": make_tsv(medians), "trace_summary.tsv": make_tsv(trace_summaries)}
    report = {"status": "pass", "normal_runs": 12, "separate_traced_runs": 4,
              "normal_medians": medians, "trace_summaries": trace_summaries, "trace_observations": observations,
              "unparsed_syscall_lines": unparsed, "source_files": source_files,
              "producer": fingerprint(Path(__file__)), "slurm_job_id": result["slurm_job_id"],
              "machine": result["machine"], "cpus": result["cpus"], "mounts": result["mounts"],
              "limitations": ["One Human sample, one current node/storage environment, three repeats per normal cell, frozen old binaries.",
                  "No final corrected-build claim or closure of the original four performance findings.",
                  "All timed runs follow complete input/binary/panel reads; this is not cold-cache benchmarking.",
                  "Trace times include tracer overhead and are excluded from all normal medians.",
                  "Syscall durations can overlap across threads. Their sum is not elapsed blocking time or the critical path.",
                  "Only decoded paths containing the exact product-output directory are attributed to product outputs.",
                  "A current output-location effect does not prove what caused the historical September wall delay.",
                  "pcr_s includes FASTA staging; wall minus pipeline includes multiple later steps, not just filesystem synchronization."],
              "outputs": {name: {"sha256": hashlib.sha256(text.encode()).hexdigest(), "size_bytes": len(text.encode())} for name, text in payloads.items()}}
    payloads["analysis.json"] = json.dumps(report, indent=2, sort_keys=True) + "\n"
    for name, text in payloads.items():
        with (args.output / name).open("x") as stream:
            stream.write(text)
    print(json.dumps({"status": "pass", "medians": medians, "unparsed_syscall_lines": len(unparsed)}))


if __name__ == "__main__":
    main()

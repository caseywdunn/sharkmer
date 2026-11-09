#!/usr/bin/env python3
"""Execute the frozen final-build comparison on one dedicated Slurm node."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import signal
import statistics
import subprocess
import sys
import time

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from sharkmer_validate import runner, uncertain

COUNTS = ("n_reads_read", "n_bases_read", "n_subreads_ingested",
          "n_bases_ingested", "n_kmers")
FLAGS = {"wall_seconds": 5.0, "peak_rss_bytes": 2.0}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while data := stream.read(1024 * 1024):
            digest.update(data)
    return {"path": str(path), "sha256": digest.hexdigest(),
            "size_bytes": path.stat().st_size}


def verify(row):
    require(fingerprint(row["path"]) == row, "Frozen file changed: " + row["path"])


def write_new(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def inspect_input(path, expected):
    """Check exact consumed FASTQ prefix, including EOF, outside timing."""
    digest = hashlib.sha256()
    size = bases = 0
    with Path(path).open("rb") as stream:
        for ordinal in range(expected["records"]):
            lines = [stream.readline(4097) for _ in range(4)]
            require(all(line.endswith(b"\n") and len(line) <= 4096 for line in lines),
                    f"Truncated FASTQ at record {ordinal + 1}")
            sequence, quality = lines[1].rstrip(b"\r\n"), lines[3].rstrip(b"\r\n")
            require(lines[0].startswith(b"@") and lines[2].startswith(b"+")
                    and sequence and len(sequence) == len(quality), "Invalid FASTQ record")
            require(all(base in b"ACGTNacgtn" for base in sequence)
                    and all(33 <= value <= 126 for value in quality), "Invalid FASTQ symbols")
            for line in lines:
                digest.update(line)
                size += len(line)
            bases += len(sequence)
        require(not stream.read(1), "Input is longer than frozen consumed prefix")
    observed = {"records": expected["records"], "bases": bases,
                "size_bytes": size, "sha256": digest.hexdigest()}
    require(observed == expected, "Consumed prefix differs: " + str(path))
    return observed


def physical_cpus():
    selected, seen, topology = [], set(), []
    for cpu in sorted(os.sched_getaffinity(0)):
        directory = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        pair = tuple(int((directory / key).read_text())
                     for key in ("physical_package_id", "core_id"))
        topology.append({"cpu": cpu, "socket": pair[0], "core": pair[1]})
        if pair not in seen and len(selected) < 2:
            seen.add(pair)
            selected.append(cpu)
    require(len(selected) == 2, "Need two distinct physical cores")
    os.sched_setaffinity(0, selected)
    return selected, topology


def mount_details(path):
    command = ["/usr/bin/findmnt", "--json", "--target", str(path),
               "--output", "TARGET,SOURCE,FSTYPE,OPTIONS"]
    return json.loads(subprocess.check_output(command, text=True, timeout=10))


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            require(header and header.split(), "Missing FASTA identifier")
            records.append({"id": header.split()[0], "header": header, "sequence": ""})
        elif line:
            require(records, "Sequence precedes FASTA header")
            records[-1]["sequence"] += line
    require(all(row["sequence"] and set(row["sequence"]) <= set("ACGT") for row in records),
            "Invalid emitted bases")
    require(len({row["id"] for row in records}) == len(records), "Duplicate FASTA IDs")
    return records


def header_fields(header):
    fields = {}
    description = header.split(maxsplit=1)
    for token in shlex.split(description[1] if len(description) == 2 else ""):
        if "=" in token:
            key, value = token.split("=", 1)
            require(key not in fields, "Duplicate header field")
            fields[key] = value
    return fields


def capture(directory, sample, role, command, protocol):
    prefix = sample["sample_prefix"]
    stats = runner._parse_stats_yaml(directory / f"{prefix}.stats.yaml")
    expected = sample["expected"][role]
    genes = {entry["gene_name"] for entry in expected["targets"]}
    require(stats["sample"] == prefix and stats["command"] == " ".join(command),
            "Output sample or command differs")
    require(stats["kmer_length"] == 19 and stats["chunks"] == 0, "Output settings differ")
    require(stats["sharkmer_version"] == protocol["versions"][role], "Wrong executable version")
    counts = {key: stats[key] for key in COUNTS}
    require(counts == expected["counts"], "Ingestion/k-mer counts differ from saved oracle")
    require(counts["n_reads_read"] == sample["expected_prefix"]["records"]
            and counts["n_bases_read"] == sample["expected_prefix"]["bases"],
            "Actual consumed records/bases differ")
    entries = stats["pcr_results"]
    require(len(entries) == len(genes) and {row["gene_name"] for row in entries} == genes,
            "Missing, duplicate, or unexpected target result")
    if role == "candidate":
        runner._validate_stats_manifest(stats, prefix, 19, genes, "end-to-end",
                                        expected_command=command)
        manifest = runner._validate_output_manifest(stats, prefix, directory)
        require(manifest["schema_version"] == stats["product_output_schema_version"] == 3
                and stats["product_output_layout"] == "combined_warnings", "Wrong output schema")
        require(stats["input_source"] == {"kind": "local_files", "inputs": [sample["input"]["path"]],
                "paired": False, "max_reads": 1000000}, "Input provenance differs")
        require(stats["stage_timings"]["read_threading_input_s"] is None,
                "Unexpected read threading")
    else:
        require(not any(key in stats for key in ("run_id", "run_status", "product_output_schema_version")),
                "Unexpected released output schema")
    observed, targets, files = [], [], set()
    for entry in entries:
        gene = entry["gene_name"]
        count = entry["n_products"]
        require(type(count) is int and count >= 0, "Invalid count")
        path = directory / f"{prefix}_{gene}.fasta"
        records = fasta(path) if count else []
        if count:
            files.add(path.name)
        require(len(records) == count and [len(row["sequence"]) for row in records]
                == entry.get("product_lengths", []), "Output count or lengths differ")
        target = {"gene_name": gene, "status": entry["status"], "n_products": count,
                  "n_standard_products": entry.get("n_standard_products", count),
                  "n_uncertain_products": entry.get("n_uncertain_products", 0),
                  "product_lengths": entry.get("product_lengths", [])}
        labels = []
        for record in records:
            fields = header_fields(record["header"])
            label = uncertain.assembly_status(record["header"], 3) if role == "candidate" else "standard"
            if role == "baseline":
                require("assembly_status" not in fields, "Unexpected legacy status field")
            labels.append(label)
            observed.append({"gene": gene, "id": record["id"], "length": len(record["sequence"]),
                             "sha256": hashlib.sha256(record["sequence"].encode("ascii")).hexdigest(),
                             "assembly_status": label,
                             "ordinary_fields": fields if label == "standard" else None})
        require(labels == sorted(labels, key=lambda value: value != "standard"),
                "Record class ordering differs")
        require(labels.count("standard") == target["n_standard_products"]
                and labels.count("uncertain") == target["n_uncertain_products"],
                "Class subtotals differ")
        targets.append(target)
    require({path.name for path in directory.glob("*.fasta")} == files,
            "Missing or unreported FASTA")
    require(targets == expected["targets"], "Target order/status/counts differ from saved oracle")
    require(len(observed) == len(expected["records"]), "Sequence inventory length differs")
    for actual, old in zip(observed, expected["records"]):
        require({key: actual[key] for key in old} == old,
                "Sequence/order/ID/class/header differs from saved oracle: " + actual["gene"])
    return {"counts": counts, "records": observed, "targets": targets,
            "stage_timings": stats.get("stage_timings"), "allocator_peak_bytes": stats["peak_memory_bytes"],
            "files": [fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()]}


def run_timed(command, directory, cpus, protocol):
    timed = ["/usr/bin/time", "-o", str(directory / "time.json"), "-f",
             '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}',
             "--", *command]
    limit = protocol["address_space_bytes"]
    launch = ["/usr/bin/prlimit", f"--as={limit}", "--", str(Path(sys.executable).resolve()),
              protocol["launcher"]["path"], "--receipt", str(directory / "constraints.json"),
              "--expected-cpus", ",".join(map(str, cpus)), "--expected-as", str(limit), "--", *timed]
    result = {"command": command, "launch_command": launch,
              "timing_scope": "GNU time directly wraps Sharkmer; input prewarming/validation excluded"}
    start = time.monotonic()
    with (directory / "stdout.log").open("x") as stdout, (directory / "stderr.log").open("x") as stderr:
        process = subprocess.Popen(launch, cwd=ROOT, env=dict(os.environ, LC_ALL="C"),
                                   stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            result["returncode"] = process.wait(timeout=protocol["invocation_timeout_seconds"])
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            result.update(returncode=process.returncode, timed_out=True)
    result["observer_wall_seconds"] = time.monotonic() - start
    result["files"] = [fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()]
    if result["returncode"]:
        return result
    constraints = json.loads((directory / "constraints.json").read_text())
    require(constraints["status"] == "pass" and constraints["affinity"] == cpus
            and constraints["rlimit_as_soft_hard"] == [limit] * 2
            and constraints["command"] == timed, "Actual process constraints differ")
    metrics = json.loads((directory / "time.json").read_text())
    require(metrics["exit_code"] == 0 and all(math.isfinite(value) and value >= 0
            for value in metrics.values()), "Invalid timing result")
    metrics["cpu_seconds"] = metrics["user_seconds"] + metrics["system_seconds"]
    metrics["peak_rss_bytes"] = metrics["max_rss_kib"] * 1024
    result["metrics"] = metrics
    return result


def summarize(runs, protocol):
    require(len(runs) == 32 and all(run["status"] == "pass" for run in runs),
            "Incomplete final comparison")
    performance = [run for run in runs if run["invocation"]["purpose"] == "performance"]
    require(len(performance) == 24, "Wrong performance scope")
    cells = []
    for sample in protocol["performance_samples"]:
        roles = {role: sorted([run for run in performance if run["invocation"]["sample"] == sample
                               and run["invocation"]["role"] == role],
                              key=lambda run: run["invocation"]["pair_index"])
                 for role in ("baseline", "candidate")}
        require(all(len(rows) == 3 for rows in roles.values()), "Incomplete matched triplicate")
        for rows in roles.values():
            require(all(row["outputs"]["records"] == rows[0]["outputs"]["records"]
                        and row["outputs"]["targets"] == rows[0]["outputs"]["targets"]
                        and row["outputs"]["counts"] == rows[0]["outputs"]["counts"] for row in rows),
                    "Repeated output is unstable")
        require(roles["baseline"][0]["outputs"]["counts"] == roles["candidate"][0]["outputs"]["counts"],
                "Paired consumed counts differ")
        keys = lambda rows: {(row["gene"], row["sha256"]) for row in rows}
        old = keys(roles["baseline"][0]["outputs"]["records"])
        new = keys(roles["candidate"][0]["outputs"]["records"])
        require(old <= new, "Historical released sequence was lost")
        cell = {"sample": sample, "released_sequences": len(old), "retained_sequences": len(old & new),
                "resources": {}, "changes": {}, "paired_deltas": []}
        for role, rows in roles.items():
            cell["resources"][role] = {}
            for metric in ("wall_seconds", "cpu_seconds", "user_seconds", "system_seconds", "peak_rss_bytes"):
                values = [row["execution"]["metrics"][metric] for row in rows]
                cell["resources"][role][metric] = {"values": values, "median": statistics.median(values),
                                                   "min": min(values), "max": max(values)}
        for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes"):
            before, after = (cell["resources"][role][metric]["median"] for role in ("baseline", "candidate"))
            percent = 100 * (after / before - 1) if before else None
            cell["changes"][metric] = {"absolute": after - before, "percent": percent,
                                       "investigation_flag": metric in FLAGS and percent is not None and percent > FLAGS[metric]}
        for before, after in zip(roles["baseline"], roles["candidate"]):
            cell["paired_deltas"].append({"pair_index": before["invocation"]["pair_index"],
                **{metric: after["execution"]["metrics"][metric] - before["execution"]["metrics"][metric]
                   for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}})
        cells.append(cell)
    candidates = {run["invocation"]["sample"]: run["outputs"] for run in runs
                  if run["invocation"]["role"] == "candidate"}
    require(len(candidates) == 12, "Need ten historical and two corrected-plant sample inventories")
    historical = [row for key, row in candidates.items() if not key.endswith("_corrected")]
    current = [row for key, row in candidates.items() if key not in protocol["corrected_replaces"]]
    counts = lambda rows: {"records": sum(len(row["records"]) for row in rows),
                          "standard": sum(record["assembly_status"] == "standard" for row in rows for record in row["records"]),
                          "warning": sum(record["assembly_status"] == "uncertain" for row in rows for record in row["records"])}
    require(counts(historical) == {"records": 326, "standard": 64, "warning": 262}, "Historical inventory total differs")
    require(counts(current) == {"records": 287, "standard": 66, "warning": 221}, "Corrected-panel inventory total differs")
    retained, retained_with_warning, total = 0, 0, 0
    for key, outputs in candidates.items():
        if key.endswith("_corrected"):
            continue
        old = {(row["gene"], row["sha256"]) for row in protocol["samples"][key]["expected"]["baseline"]["records"]}
        latest = candidates.get(key + "_corrected", outputs)
        new = {(row["gene"], row["sha256"]) for row in latest["records"]}
        warnings = {(row["gene"], row["sha256"]) for row in latest["records"]
                    if row["assembly_status"] == "uncertain"}
        require(old <= new, "Historical released identity missing from current-panel inventory: " + key)
        total += len(old)
        retained += len(old & new)
        retained_with_warning += len(old & warnings)
    require(total == retained == 86 and retained_with_warning == 22,
            "Released identity retention/classification totals differ")
    return {"status": "pass", "cells": cells, "historical_panel_inventory": counts(historical),
            "corrected_panel_inventory": counts(current), "invocations": len(runs),
            "released_identities": total, "retained_identities": retained,
            "retained_with_warning": retained_with_warning, "lost_identities": total - retained,
            "note": "Exact scoped inventory verification; performance flags require review. Not statistical equivalence or biological validation."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(args.protocol.read_text())
    state = {"status": "failed", "protocol": fingerprint(args.protocol), "runs": [],
             "started": time.time(), "slurm_job_id": os.environ.get("SLURM_JOB_ID")}
    try:
        require(state["slurm_job_id"] and int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) == 2,
                "Requires a dedicated two-CPU Slurm allocation")
        require(protocol["schema_version"] == 1 and len(protocol["schedule"]) == 32,
                "Unsupported protocol or invocation scope")
        state["bindings_before"] = [fingerprint(row["path"]) for row in protocol["bindings"]]
        require(state["bindings_before"] == protocol["bindings"], "Frozen source/tool bindings changed")
        cpus, topology = physical_cpus()
        state.update(cpus=cpus, topology=topology, machine=platform.uname()._asdict(),
                     python=sys.version, python_executable=fingerprint(sys.executable), pyyaml=yaml.__version__,
                     output_mount=mount_details(output))
        require(state["output_mount"]["filesystems"][0]["fstype"].startswith("nfs"),
                "This protocol requires shared NFS output")
        state["input_preflight"] = {}
        for sample in protocol["samples"].values():
            path = sample["input"]["path"]
            if path not in state["input_preflight"]:
                checked = inspect_input(path, sample["expected_prefix"])
                require({"path": path, "sha256": checked["sha256"], "size_bytes": checked["size_bytes"]}
                        == sample["input"], "Input metadata differs from consumed-prefix identity")
                state["input_preflight"][path] = checked
            require(state["input_preflight"][path] == sample["expected_prefix"],
                    "Shared input has conflicting prefix expectations")
        write_new(output / "preflight.json", {key: value for key, value in state.items() if key != "runs"})
        for item in protocol["schedule"]:
            require(time.time() - state["started"] < protocol["total_start_budget_seconds"],
                    "Total start budget exhausted")
            sample, role = protocol["samples"][item["sample"]], item["role"]
            directory = output / item["id"]
            directory.mkdir()
            products = directory / "outputs"
            products.mkdir()
            run = {"status": "failed", "invocation": item, "started": time.time()}
            state["runs"].append(run)
            try:
                run["input_before"] = fingerprint(sample["input"]["path"])
                require(run["input_before"] == sample["input"], "Frozen input changed before invocation")
                for row in (sample["panel"], protocol["binaries"][role]):
                    verify(row)
                command = [protocol["binaries"][role]["path"], "-k", "19", "-t", "2", "--chunks", "0",
                           "--max-reads", "1000000", "-o", str(products) + "/", "-s", sample["sample_prefix"],
                           "--pcr-panel-file", sample["panel"]["path"], sample["input"]["path"]]
                run["execution"] = run_timed(command, directory, cpus, protocol)
                require(run["execution"]["returncode"] == 0, "Timed CLI failed or exceeded its timeout")
                run["outputs"] = capture(products, sample, role, command, protocol)
                run["input_after"] = fingerprint(sample["input"]["path"])
                require(run["input_after"] == run["input_before"], "Input changed during invocation")
                for row in (sample["panel"], protocol["binaries"][role]):
                    verify(row)
                run["status"] = "pass"
            finally:
                run["finished"] = time.time()
                write_new(directory / "receipt.json", run)
            print(json.dumps({"completed": len(state["runs"]), "id": item["id"], "status": run["status"]}), flush=True)
        state["bindings_after"] = [fingerprint(row["path"]) for row in protocol["bindings"]]
        require(state["bindings_after"] == state["bindings_before"], "Frozen source/tool bindings changed during run")
        verify(state["protocol"])
        summary = summarize(state["runs"], protocol)
        write_new(output / "summary.json", summary)
        state["summary"] = fingerprint(output / "summary.json")
        state["status"] = "pass"
    except Exception as error:
        state["error"] = repr(error)
        raise
    finally:
        state["finished"] = time.time()
        write_new(output / "result.json", state)


if __name__ == "__main__":
    main()

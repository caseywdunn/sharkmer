#!/usr/bin/env python3
"""Execute the frozen, bounded Human output-location diagnostic, sequentially."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import statistics
import subprocess
import sys
import tempfile
import time

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ADDRESS_SPACE = 12 * 1024**3
SAMPLE = "human_SRR17535371_1000000"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            h.update(chunk)
    return {"path": str(path.resolve()), "sha256": h.hexdigest(), "size_bytes": path.stat().st_size}


def verify(row):
    require(fingerprint(row["path"]) == row, "Source fingerprint changed: " + row["path"])


def write_new(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def inspect_input(path, expected):
    h = hashlib.sha256()
    bases = size = 0
    with Path(path).open("rb") as stream:
        for ordinal in range(expected["records"]):
            lines = [stream.readline(4097) for _ in range(4)]
            require(all(line.endswith(b"\n") and len(line) <= 4096 for line in lines), "Truncated FASTQ")
            seq, qual = lines[1].rstrip(b"\r\n"), lines[3].rstrip(b"\r\n")
            require(lines[0].startswith(b"@") and lines[2].startswith(b"+") and seq and len(seq) == len(qual), "Invalid FASTQ")
            require(all(base in b"ACGTNacgtn" for base in seq) and all(33 <= q <= 126 for q in qual), "Invalid FASTQ symbols")
            for line in lines:
                h.update(line)
                size += len(line)
            bases += len(seq)
        require(not stream.read(1), "Input contains more than frozen EOF")
    observed = {"records": expected["records"], "bases": bases, "size_bytes": size, "sha256": h.hexdigest()}
    require(observed == expected, "Frozen input prefix differs")
    return observed


def physical_cpus():
    chosen, cores, topology = [], set(), []
    for cpu in sorted(os.sched_getaffinity(0)):
        root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = tuple(int((root / key).read_text()) for key in ("physical_package_id", "core_id"))
        topology.append({"cpu": cpu, "socket": core[0], "core": core[1]})
        if core not in cores and len(chosen) < 2:
            cores.add(core)
            chosen.append(cpu)
    require(len(chosen) == 2, "Need two distinct physical cores")
    os.sched_setaffinity(0, chosen)
    return chosen, topology


def mount_details(path):
    result = subprocess.run(["/usr/bin/findmnt", "--json", "--target", str(path), "--output", "TARGET,SOURCE,FSTYPE,OPTIONS"],
                            capture_output=True, text=True, timeout=10, check=True)
    return json.loads(result.stdout)


def capture_outputs(path, role, protocol, command):
    stats_path = path / (SAMPLE + ".stats.yaml")
    stats = yaml.safe_load(stats_path.read_text())
    require(stats["sample"] == SAMPLE and stats["kmer_length"] == 19 and stats["chunks"] == 0, "Output settings differ")
    require(stats["command"] == " ".join(command), "Recorded CLI differs")
    expected = protocol["expected_outputs"][role]
    counts = {key: stats[key] for key in expected["counts"]}
    require(counts == expected["counts"], "Input/ingestion counts differ")
    records = []
    for fasta in sorted(path.glob("*.fasta")):
        gene = fasta.stem.removeprefix(SAMPLE + "_")
        header = None
        sequence = []
        for line in fasta.read_text().splitlines() + [">"]:
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(sequence)
                    status = re.findall(r"(?:^|\s)assembly_status=(\S+)", header)
                    if role == "retention":
                        require(status == ["standard"], "Unexpected Human product status")
                    else:
                        require(not status, "Unexpected release header schema")
                    require(seq and set(seq) <= set("ACGT"), "Invalid output bases")
                    records.append({"gene": gene, "sha256": hashlib.sha256(seq.encode()).hexdigest(),
                                    "length": len(seq), "assembly_status": "standard"})
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    records.sort(key=lambda row: (row["gene"], row["sha256"]))
    require(records == expected["records"], "Human full sequence identities differ")
    entries = stats["pcr_results"]
    observed_targets = sorted((row["gene_name"], row["status"], row["n_products"]) for row in entries)
    require(observed_targets == [tuple(row) for row in expected["targets"]], "Target set/dispositions differ")
    if role == "retention":
        require(stats["run_status"] == "complete" and stats["product_output_layout"] == "combined_classified", "Unexpected run schema")
        require(stats["stage_timings"]["read_threading_input_s"] is None, "Read threading unexpectedly active")
        require(stats["input_source"] == {"kind": "local_files", "inputs": [protocol["input"]["path"]],
                "paired": False, "max_reads": 1000000}, "Input provenance differs")
    return {"counts": counts, "records": records, "stage_timings": stats.get("stage_timings"),
            "allocator_peak_bytes": stats["peak_memory_bytes"],
            "files": [fingerprint(f) for f in sorted(path.rglob("*")) if f.is_file()]}


def summarize(runs):
    normal = [r for r in runs if not r["schedule"]["traced"]]
    require(len(normal) == 12 and all(r["status"] == "pass" for r in runs), "Incomplete successful diagnostic")
    groups, contrasts = [], []
    for location in ("nfs", "local"):
        for role in ("release", "retention"):
            rows = [r for r in normal if r["schedule"]["location"] == location and r["schedule"]["role"] == role]
            item = {"location": location, "role": role, "n": len(rows)}
            for field in ("wall_seconds", "cpu_seconds", "user_seconds", "system_seconds", "max_rss_kib"):
                values = [r["metrics"][field] for r in rows]
                item[field] = {"values": values, "median": statistics.median(values), "min": min(values), "max": max(values)}
            groups.append(item)
        for block in (1, 2, 3):
            rows = {r["schedule"]["role"]: r for r in normal if r["schedule"]["location"] == location and r["schedule"]["block"] == block}
            before, after = (rows[role]["metrics"] for role in ("release", "retention"))
            contrasts.append({"location": location, "block": block,
                              "wall_delta_seconds": after["wall_seconds"] - before["wall_seconds"],
                              "wall_percent": 100 * (after["wall_seconds"] / before["wall_seconds"] - 1),
                              "cpu_delta_seconds": after["cpu_seconds"] - before["cpu_seconds"]})
    return {"groups": groups, "paired_contrasts": contrasts, "traced_runs_excluded": 4,
            "scope": "Human, warm input, frozen released/retention binaries, one node. Not release acceptance, not an assembly or retention ablation, not historical causation."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=HERE / "protocol-v1.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "failed", "runs": [], "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "protocol": fingerprint(args.protocol), "start_unix_time": time.time()}
    try:
        require(os.environ.get("SLURM_JOB_ID"), "Must run in dedicated Slurm allocation")
        require(int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) == 2, "Need a two-CPU Slurm task")
        protocol = json.loads(args.protocol.read_text())
        for row in protocol["sources"]:
            verify(row)
        cpus, topology = physical_cpus()
        result.update(machine=platform.uname()._asdict(), cpus=cpus, topology=topology,
                      python=sys.version, pyyaml=yaml.__version__)
        local = Path(tempfile.mkdtemp(prefix=f"sharkmer-runtime-{result['slurm_job_id']}-", dir="/tmp"))
        result["local_root"] = str(local)
        result["mounts"] = {"nfs_output": mount_details(output), "local_output": mount_details(local), "input": mount_details(protocol["input"]["path"])}
        require(result["mounts"]["nfs_output"]["filesystems"][0]["fstype"].startswith("nfs"), "Shared output is not NFS")
        require(not result["mounts"]["local_output"]["filesystems"][0]["fstype"].startswith("nfs"), "Temporary output is not node-local")
        write_new(output / "preflight.json", {k: v for k, v in result.items() if k != "runs"})
        environment = dict(os.environ, LC_ALL="C")
        for item in protocol["schedule"]:
            require(time.time() - result["start_unix_time"] < 780, "Internal 13-minute start budget exhausted")
            name = item["id"]
            saved = output / name
            saved.mkdir()
            local_run = local / name
            local_run.mkdir()
            products = (saved if item["location"] == "nfs" else local_run) / "outputs"
            products.mkdir()
            run = {"status": "failed", "schedule": item, "output_path_during_timing": str(products), "start_unix_time": time.time()}
            result["runs"].append(run)
            try:
                # Explicit repeated input/code/panel reads, all outside the timer.
                run["input_preflight"] = inspect_input(protocol["input"]["path"], protocol["expected_input"])
                for row in (protocol["input"], protocol["panel"], protocol["binaries"][item["role"]]):
                    verify(row)
                command = [protocol["binaries"][item["role"]]["path"], "-k", "19", "-t", "2", "--chunks", "0", "--max-reads", "1000000",
                           "-o", str(products) + "/", "-s", SAMPLE, "--pcr-panel-file", protocol["panel"]["path"], protocol["input"]["path"]]
                observed_command = command
                if item["traced"]:
                    observed_command = ["/usr/bin/strace", "-ff", "-ttt", "-T", "-yy", "-s", "0", "-o", str(local_run / "trace"),
                                        "-e", "trace=" + protocol["trace_syscalls"], "--", *command]
                timed = ["/usr/bin/time", "-o", str(local_run / "time.json"), "-f",
                         '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}', "--", *observed_command]
                constraints = local_run / "process_constraints.json"
                launch = ["/usr/bin/prlimit", "--as=" + str(ADDRESS_SPACE), "--", str(Path(sys.executable).resolve()),
                          protocol["launcher"]["path"], "--receipt", str(constraints), "--expected-cpus", ",".join(map(str, cpus)),
                          "--expected-as", str(ADDRESS_SPACE), "--", *timed]
                run.update(command=command, launch_command=launch, timing_scope="traced diagnostic including tracer overhead" if item["traced"] else "GNU time directly wraps Sharkmer")
                start = time.monotonic()
                with (local_run / "stdout.log").open("x") as stdout, (local_run / "stderr.log").open("x") as stderr:
                    process = subprocess.Popen(launch, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr, start_new_session=True)
                    try:
                        run["returncode"] = process.wait(timeout=120 if item["traced"] else 60)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                        run.update(timed_out=True, returncode=process.returncode)
                run["observer_wall_seconds"] = time.monotonic() - start
                require(run["returncode"] == 0, "CLI or tracer failed")
                observed = json.loads(constraints.read_text())
                require(observed["status"] == "pass" and observed["affinity"] == cpus and observed["rlimit_as_soft_hard"] == [ADDRESS_SPACE] * 2 and observed["command"] == timed, "Constraints differ")
                metrics = json.loads((local_run / "time.json").read_text())
                require(metrics["exit_code"] == 0 and all(math.isfinite(v) and v >= 0 for v in metrics.values()), "Invalid GNU time metrics")
                metrics["cpu_seconds"] = metrics["user_seconds"] + metrics["system_seconds"]
                run["metrics"] = metrics
                run["outputs"] = capture_outputs(products, item["role"], protocol, command)
                for row in (protocol["input"], protocol["panel"], protocol["binaries"][item["role"]]):
                    verify(row)
                run["status"] = "pass"
            finally:
                # Preserve all node-local outputs and raw traces outside the timer,
                # including the partial artifacts of any failed invocation.
                shutil.copytree(local_run, saved / "local_artifacts")
                run["preserved_local_files"] = [fingerprint(f) for f in sorted((saved / "local_artifacts").rglob("*")) if f.is_file()]
                run["finished_unix_time"] = time.time()
                write_new(saved / "receipt.json", run)
            print(json.dumps({"completed": len(result["runs"]), "id": name, "status": run["status"]}), flush=True)
        for row in protocol["sources"]:
            verify(row)
        result["summary"] = summarize(result["runs"])
        write_new(output / "summary.json", result["summary"])
        result["status"] = "pass"
    except Exception as error:
        result["error"] = repr(error)
    finally:
        result["finished_unix_time"] = time.time()
        write_new(output / "result.json", result)
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Compare direct counting with the frozen v3.2.1 engineering baseline."""

import argparse
import datetime
import hashlib
import json
import math
import os
import platform
import signal
import statistics
import subprocess
from pathlib import Path
from typing import Any

import yaml

BASELINE_SHA256 = "0d57521d5c6060df3e954d366b72fdccb0b3e40f9532162b400d62b14a59109c"


def digest(path: Path) -> str:
    """Hash an artifact without loading it in memory."""
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def save(path: Path, value: dict[str, Any]) -> None:
    """Replace a receipt only after its complete JSON has been written."""
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def records(directory: Path) -> list[dict[str, Any]]:
    """Read complete record identities and warning reasons from FASTA outputs."""
    observed = []
    for fasta in sorted(directory.glob("*.fasta")):
        for block in fasta.read_text().split(">")[1:]:
            header, *lines = block.splitlines()
            fields = dict(part.split("=", 1) for part in header.split()[1:])
            sequence = "".join(lines)
            observed.append(
                {
                    "id": header.split()[0],
                    "gene": fields.get("target", fields.get("gene")),
                    "length": len(sequence),
                    "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "assembly_status": fields.get("assembly_status", "standard"),
                    "assembly_status_causes": sorted(
                        fields.get("assembly_status_causes", "").split(",")
                    )
                    if fields.get("assembly_status_causes")
                    else [],
                }
            )
            assert int(fields.get("sequence_length", fields.get("length"))) == len(
                sequence
            )
    return sorted(observed, key=lambda row: row["id"])


def invoke(
    command: list[str], work: Path, timeout: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Retain GNU time and output receipts; kill the entire process group on error."""
    work.mkdir(parents=True, exist_ok=False)
    metrics_file = work / "time.json"
    timed = [
        "/usr/bin/time",
        "-f",
        '{"wall_s":%e,"user_s":%U,"system_s":%S,"rss_kib":%M}',
        "-o",
        str(metrics_file),
        *command,
    ]
    (work / "command.json").write_text(json.dumps(timed, indent=2) + "\n")
    with (
        (work / "stdout.log").open("wb") as stdout,
        (work / "stderr.log").open("wb") as stderr,
        subprocess.Popen(
            timed, stdout=stdout, stderr=stderr, start_new_session=True
        ) as child,
    ):
        try:
            returncode = child.wait(timeout=timeout)
        except BaseException:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
            raise
    assert returncode == 0, f"{work}: process exited {returncode}"
    metrics = json.loads(metrics_file.read_text())
    destination = work / "outputs"
    stats_file = next(destination.glob("*.stats.yaml"))
    manifest_file = next(destination.glob("*.manifest.yaml"))
    stats = yaml.safe_load(stats_file.read_text())
    manifest = yaml.safe_load(manifest_file.read_text())
    assert stats["run_status"] == manifest["status"] == "complete"
    assert stats["run_id"] == manifest["run_id"]
    for artifact in manifest["files"]:
        assert digest(destination / artifact["path"]) == artifact["sha256"]
    return metrics, stats


def summarize(runs: list[dict[str, Any]], protocol: dict[str, Any]) -> dict[str, Any]:
    """Apply the previously frozen wall/RSS gates without changing tolerances."""
    gates = protocol["performance_gates"]
    samples = []
    for sample in protocol["samples"]:
        selected = [
            row
            for row in runs
            if row["accession"] == sample["accession"] and row["kind"] == "timed"
        ]
        medians = {}
        for backend in ["baseline", "candidate"]:
            rows = [row for row in selected if row["backend"] == backend]
            assert len(rows) == protocol["timed_replicates"]
            medians[backend] = {
                "wall_s": statistics.median(row["metrics"]["wall_s"] for row in rows),
                "rss_bytes": statistics.median(
                    row["metrics"]["rss_kib"] * 1024 for row in rows
                ),
                "allocator_peak_bytes": statistics.median(
                    row["allocator_peak_bytes"] for row in rows
                ),
                "finalize_s": statistics.median(
                    row["stages"]["count_finalize_s"] for row in rows
                ),
            }
        baseline, candidate = medians["baseline"], medians["candidate"]
        wall_limit = max(
            baseline["wall_s"] * gates["per_sample_median_wall_ratio_max"],
            baseline["wall_s"] + gates["per_sample_wall_absolute_allowance_seconds"],
        )
        rss_limit = max(
            baseline["rss_bytes"] * gates["per_sample_median_rss_ratio_max"],
            baseline["rss_bytes"] + gates["per_sample_rss_absolute_allowance_bytes"],
        )
        samples.append(
            {
                "accession": sample["accession"],
                "medians": medians,
                "wall_ratio": candidate["wall_s"] / baseline["wall_s"],
                "wall_limit_s": wall_limit,
                "rss_limit_bytes": rss_limit,
                "wall_passed": candidate["wall_s"] <= wall_limit,
                "rss_passed": candidate["rss_bytes"] <= rss_limit,
            }
        )
    geomean = math.exp(statistics.mean(math.log(row["wall_ratio"]) for row in samples))
    return {
        "samples": samples,
        "geomean_wall_ratio": geomean,
        "performance_gates_passed": geomean
        <= gates["default_in_memory_geomean_wall_ratio_max"]
        and all(row["wall_passed"] and row["rss_passed"] for row in samples),
        "scope": "paired two-thread performance and all-entry count fingerprints; thread/build agreement is a separate gate",
    }


def verify_inputs(protocol: dict[str, Any]) -> None:
    """Recheck full input prefixes and panels before and after measurement."""
    for sample in protocol["samples"]:
        source = sample["input"]
        assert Path(source["path"]).stat().st_size == source["size_bytes"]
        assert digest(Path(source["path"])) == source["sha256"]
        assert digest(Path(sample["panel_path"])) == sample["panel_sha256"]


def main() -> None:
    """Run each retained case with all attempted invocations preserved."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert os.name == "posix", "Process-group cleanup requires POSIX"
    protocol_file = args.protocol.resolve()
    protocol_sha = digest(protocol_file)
    assert protocol_sha == protocol_file.with_suffix(".sha256").read_text().split()[0]
    protocol = json.loads(protocol_file.read_text())
    assert protocol["timed_replicates"] == 3
    binaries = {
        "baseline": args.baseline.resolve(),
        "candidate": args.candidate.resolve(),
    }
    assert digest(binaries["baseline"]) == BASELINE_SHA256
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {
        "status": "running",
        "protocol_sha256": protocol_sha,
        "harness_sha256": digest(Path(__file__)),
        "binaries": {
            name: {"path": str(path), "sha256": digest(path)}
            for name, path in binaries.items()
        },
        "hostname": platform.node(),
        "platform": platform.platform(),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "started_at_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "runs": [],
    }
    save(output / "comparison.json", receipt)
    try:
        verify_inputs(protocol)
        for sample in protocol["samples"]:
            accession = sample["accession"]
            expected = sorted(sample["expected_records"], key=lambda row: row["id"])
            inventory = None
            fingerprints = []
            schedule = [
                ("verify", "legacy-verify", 0),
                ("verify", "candidate", 0),
                ("warmup", "baseline", 0),
                ("warmup", "candidate", 0),
            ]
            for pair in range(protocol["timed_replicates"]):
                order = (
                    ["baseline", "candidate"]
                    if pair % 2 == 0
                    else ["candidate", "baseline"]
                )
                schedule.extend(("timed", backend, pair) for backend in order)
            for index, (kind, backend, pair) in enumerate(schedule):
                work = output / accession / f"{index:02d}-{kind}-{backend}"
                binary = binaries["baseline" if backend == "baseline" else "candidate"]
                command = [
                    str(binary),
                    "-k",
                    str(protocol["k"]),
                    "-t",
                    str(protocol["threads"]),
                    "--chunks",
                    "0",
                    "--max-reads",
                    str(sample["max_reads"]),
                    "-s",
                    sample["sample"],
                    "-o",
                    str(work / "outputs"),
                    "--pcr-panel-file",
                    sample["panel_path"],
                ]
                if backend != "baseline":
                    command += [
                        "--count-backend",
                        "legacy" if backend == "legacy-verify" else "memory",
                    ]
                if kind == "verify":
                    command += ["--count-fingerprint"]
                command.append(sample["input"]["path"])
                metrics, stats = invoke(command, work, protocol["timeout_seconds"])
                assert {key: stats[key] for key in sample["expected_counts"]} == sample[
                    "expected_counts"
                ]
                observed = records(work / "outputs")
                assert [
                    {key: row[key] for key in reference}
                    for row, reference in zip(observed, expected, strict=True)
                ] == expected
                if inventory is None:
                    inventory = observed
                assert observed == inventory, (
                    f"{accession}: record/warning inventory differs"
                )
                if kind == "verify":
                    fingerprints.append(stats["count_fingerprint"])
                    if len(fingerprints) == 2:
                        assert fingerprints[0] == fingerprints[1], (
                            f"{accession}: count pairs differ"
                        )
                    assert stats["count_fingerprint"]["occurrences"] == stats["n_kmers"]
                row = {
                    "accession": accession,
                    "kind": kind,
                    "backend": backend,
                    "pair": pair,
                    "command": command,
                    "metrics": metrics,
                    "stages": stats["stage_timings"],
                    "allocator_peak_bytes": stats["peak_memory_bytes"],
                    "records": observed,
                    "count_store": stats.get("count_store"),
                    "count_fingerprint": stats.get("count_fingerprint"),
                }
                receipt["runs"].append(row)
                save(output / "comparison.json", receipt)
            assert len(fingerprints) == 2 and fingerprints[0] == fingerprints[1], (
                f"{accession}: count pairs differ"
            )
        verify_inputs(protocol)
        receipt["summary"] = summarize(receipt["runs"], protocol)
        receipt["status"] = (
            "passed"
            if receipt["summary"]["performance_gates_passed"]
            else "performance_gates_failed"
        )
        assert digest(protocol_file) == protocol_sha
        assert digest(Path(__file__)) == receipt["harness_sha256"]
        for backend, binary in binaries.items():
            assert digest(binary) == receipt["binaries"][backend]["sha256"]
    except BaseException as error:
        receipt["status"] = "failed"
        receipt["error"] = repr(error)
        raise
    finally:
        receipt["finished_at_utc"] = datetime.datetime.now(datetime.UTC).isoformat()
        save(output / "comparison.json", receipt)


if __name__ == "__main__":
    main()

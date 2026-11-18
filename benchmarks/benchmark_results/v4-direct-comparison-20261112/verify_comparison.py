#!/usr/bin/env python3
"""Independently reconcile complete paired comparisons with their raw artifacts."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import yaml


def sha(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def inventory(directory: Path) -> list[dict]:
    rows = []
    for path in sorted(directory.glob("*.fasta")):
        header, bases = None, []

        def finish(header: str | None, bases: list[str]) -> None:
            if header is None:
                return
            parts = header.split()
            fields = dict(part.split("=", 1) for part in parts[1:])
            sequence = "".join(bases)
            status = fields.get("assembly_status", "standard")
            status = "warning" if status == "uncertain" else status
            assert int(fields.get("length", fields.get("sequence_length"))) == len(
                sequence
            )
            causes = (
                sorted(fields.get("assembly_status_causes", "").split(","))
                if fields.get("assembly_status_causes")
                else []
            )
            assert status != "warning" or causes
            rows.append(
                {
                    "id": parts[0],
                    "gene": fields.get("target", fields.get("gene")),
                    "length": len(sequence),
                    "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "assembly_status": status,
                    "assembly_status_causes": causes,
                }
            )

        for line in path.read_text().splitlines():
            if line.startswith(">"):
                finish(header, bases)
                header, bases = line[1:], []
            else:
                assert header is not None
                bases.append(line)
        finish(header, bases)
    assert len(rows) == len({row["id"] for row in rows})
    return sorted(rows, key=lambda row: row["id"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    receipt = json.loads((root / "comparison.json").read_text())
    protocol = json.loads(args.protocol.read_text())
    assert sha(args.protocol) == receipt["protocol_sha256"]
    assert len(receipt["runs"]) == 100
    assert protocol["timed_replicates"] == 3 and len(protocol["samples"]) == 10
    for binary in receipt["binaries"].values():
        assert sha(Path(binary["path"])) == binary["sha256"]
    samples = {sample["accession"]: sample for sample in protocol["samples"]}
    seen, reference, fingerprints, timed = set(), {}, {}, {}
    for run in receipt["runs"]:
        accession, backend, kind, pair = (
            run[key] for key in ["accession", "backend", "kind", "pair"]
        )
        token = (accession, backend, kind, pair)
        assert token not in seen
        seen.add(token)
        sample = samples[accession]
        command = run["command"]
        assert command[-1] == sample["input"]["path"]
        assert command[command.index("--max-reads") + 1] == str(sample["max_reads"])
        assert command[command.index("--pcr-panel-file") + 1] == sample["panel_path"]
        assert command[command.index("-k") + 1] == str(protocol["k"])
        assert command[command.index("-t") + 1] == str(protocol["threads"])
        assert ("--count-fingerprint" in command) == (kind == "verify")
        if backend != "baseline":
            assert command[command.index("--count-backend") + 1] == (
                "legacy" if backend == "legacy-verify" else "memory"
            )
        destination = Path(command[command.index("-o") + 1])
        assert destination.is_relative_to(root)
        work = destination.parent
        assert (
            json.loads((work / "command.json").read_text())[-len(command) :] == command
        )
        metrics = json.loads((work / "time.json").read_text())
        assert metrics == run["metrics"]
        stats = yaml.safe_load(
            (destination / f"{sample['sample']}.stats.yaml").read_text()
        )
        manifest = yaml.safe_load(
            (destination / f"{sample['sample']}.manifest.yaml").read_text()
        )
        assert stats["run_status"] == manifest["status"] == "complete"
        assert stats["run_id"] == manifest["run_id"]
        assert stats["stage_timings"] == run["stages"]
        assert stats["peak_memory_bytes"] == run["allocator_peak_bytes"]
        for artifact in manifest["files"]:
            file = destination / artifact["path"]
            assert (
                file.is_file()
                and not file.is_symlink()
                and file.resolve().is_relative_to(destination)
            )
            assert sha(file) == artifact["sha256"]
        for field, expected in sample["expected_counts"].items():
            assert stats[field] == expected
        actual = inventory(destination)
        assert actual == run["records"]
        expected = sorted(sample["expected_records"], key=lambda row: row["id"])
        assert len(actual) == len(expected)
        for row, target in zip(actual, expected, strict=True):
            for field, value in target.items():
                assert row[field] == (
                    "warning"
                    if field == "assembly_status" and value == "uncertain"
                    else value
                )
        assert reference.setdefault(accession, actual) == actual
        if kind == "verify":
            fingerprint = stats["count_fingerprint"]
            assert fingerprint == run["count_fingerprint"]
            assert fingerprint["occurrences"] == stats["n_kmers"] < (1 << 32)
            assert fingerprints.setdefault(accession, fingerprint) == fingerprint
        else:
            assert "count_fingerprint" not in stats
        if backend == "candidate":
            metadata = stats["count_store"]
            assert metadata["count_width_bits"] == 64
            assert metadata["occurrences"] == stats["n_kmers"]
            assert len(bytes.fromhex(metadata["selected_input_sha256"])) == 32
            assert (
                metadata["owned_slot_bytes"]
                <= metadata["peak_slot_bytes"]
                <= metadata["peak_count_input_reserved_bytes"]
            )
        if kind == "timed":
            timed.setdefault((accession, backend), []).append(metrics)
    expected_tokens = set()
    summaries = []
    limits = protocol["performance_gates"]
    for accession in samples:
        expected_tokens.update(
            (accession, backend, kind, 0)
            for kind, backend in [
                ("verify", "legacy-verify"),
                ("verify", "candidate"),
                ("warmup", "baseline"),
                ("warmup", "candidate"),
            ]
        )
        expected_tokens.update(
            (accession, backend, "timed", pair)
            for backend in ["baseline", "candidate"]
            for pair in range(3)
        )
        medians = {}
        for backend in ["baseline", "candidate"]:
            rows = timed[(accession, backend)]
            assert len(rows) == 3
            medians[backend] = {
                "wall": sorted(row["wall_s"] for row in rows)[1],
                "rss": sorted(row["rss_kib"] for row in rows)[1] * 1024,
            }
        baseline, candidate = medians["baseline"], medians["candidate"]
        wall = candidate["wall"] <= max(
            baseline["wall"] * limits["per_sample_median_wall_ratio_max"],
            baseline["wall"] + limits["per_sample_wall_absolute_allowance_seconds"],
        )
        rss = candidate["rss"] <= max(
            baseline["rss"] * limits["per_sample_median_rss_ratio_max"],
            baseline["rss"] + limits["per_sample_rss_absolute_allowance_bytes"],
        )
        summaries.append(
            {
                "accession": accession,
                "wall_ratio": candidate["wall"] / baseline["wall"],
                "wall_passed": wall,
                "rss_passed": rss,
            }
        )
    assert seen == expected_tokens
    geomean = math.prod(row["wall_ratio"] for row in summaries) ** (1 / len(summaries))
    passed = geomean <= limits["default_in_memory_geomean_wall_ratio_max"] and all(
        row["wall_passed"] and row["rss_passed"] for row in summaries
    )
    assert math.isclose(
        geomean, receipt["summary"]["geomean_wall_ratio"], rel_tol=1e-12
    )
    assert passed == receipt["summary"]["performance_gates_passed"]
    assert receipt["status"] == ("passed" if passed else "performance_gates_failed")
    for sample in samples.values():
        assert sha(Path(sample["input"]["path"])) == sample["input"]["sha256"]
        assert sha(Path(sample["panel_path"])) == sample["panel_sha256"]
    output = {
        "verification": "passed",
        "invocations": len(seen),
        "samples": len(samples),
        "records": sum(len(rows) for rows in reference.values()),
        "all_count_fingerprints_agree": True,
        "geomean_wall_ratio": geomean,
        "promotion_gates_passed": passed,
        "sample_performance": summaries,
        "scope": "independent raw-artifact/count/record/gate audit; a rejected candidate remains rejected",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

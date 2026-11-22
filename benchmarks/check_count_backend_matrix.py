#!/usr/bin/env python3
"""Check full count/product equivalence across hashers and thread counts."""

import argparse
import datetime
import json
import os
from pathlib import Path

from compare_count_backends import (
    backend_arguments,
    digest,
    expected_records,
    invoke,
    records,
    save,
    verify_inputs,
)


def verify_backend(metadata: dict, backend: str, threads: int, hasher: str) -> None:
    """Check actual ownership rather than assuming every candidate is serial."""
    assert metadata["hash_algorithm"] == hasher
    if backend == "memory":
        assert metadata["backend"] == "hashbrown_memory"
        assert "parallel" not in metadata
        assert "input_identity_worker" not in metadata
        assert "count_owner_transfer" not in metadata
        return
    assert backend in ("hybrid", "default")
    transfer = metadata.get("count_owner_transfer")
    workers = threads - 2 if threads > 2 else int(transfer is not None)
    assert ("input_identity_worker" in metadata) == (threads > 1)
    if transfer is not None:
        assert threads == 2
        assert transfer["threshold_bytes"] == 16 << 20
        assert transfer["current_table_bytes"] >= transfer["threshold_bytes"]
        assert transfer["selected_records_before_transfer"] > 0
        assert transfer["occurrences_before_transfer"] > 0
        assert (
            metadata["input_identity_worker"]["selected_records"]
            == transfer["selected_records_before_transfer"]
        )
    if workers:
        assert metadata["backend"] == "adaptive_sharded_exact"
        assert metadata["parallel"]["counting_workers"] == workers
        shards = metadata["adaptive_parallel"]["shards"]
        # These fixed calibration inputs fit the unchanged default RAM budget.
        # Native spill correctness is exercised by its separate input matrix.
        assert shards["memory_shards"] == workers
        assert shards["disk_shards"] == 0
    else:
        assert metadata["backend"] == "hashbrown_memory"
        assert "parallel" not in metadata and "adaptive_parallel" not in metadata


def main() -> None:
    """Run the frozen calibration inputs through every supported execution mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--ahash", type=Path, required=True)
    parser.add_argument("--fxhash", type=Path, required=True)
    parser.add_argument(
        "--candidate-backend", choices=("memory", "hybrid", "default"), default="memory"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    reference = json.loads(args.reference.read_text())
    assert reference["protocol_sha256"] == digest(args.protocol)
    # A performance-rejected table can still supply its independently audited
    # exact counts. The reference pair was also checked against legacy counts;
    # total occurrences below u32::MAX exclude legacy saturation here.
    expected = {}
    inventories = {}
    for row in reference["runs"]:
        if row["kind"] == "verify":
            accession = row["accession"]
            assert (
                expected.setdefault(accession, row["count_fingerprint"])
                == row["count_fingerprint"]
            )
            assert inventories.setdefault(accession, row["records"]) == row["records"]
    binaries = {"ahash": args.ahash.resolve(), "fxhash": args.fxhash.resolve()}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {
        "status": "running",
        "protocol_sha256": digest(args.protocol),
        "reference_sha256": digest(args.reference),
        "harness_sha256": digest(Path(__file__)),
        "comparison_helpers_sha256": digest(
            Path(__file__).with_name("compare_count_backends.py")
        ),
        "candidate_backend": args.candidate_backend,
        "binaries": {
            name: {"path": str(path), "sha256": digest(path)}
            for name, path in binaries.items()
        },
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "started_at_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "runs": [],
    }
    save(output / "matrix.json", receipt)
    try:
        verify_inputs(protocol)
        for sample in protocol["samples"]:
            accession = sample["accession"]
            expected_digest = None
            for hasher, binary in binaries.items():
                for threads in [1, 2, 4, 8]:
                    work = output / accession / f"{hasher}-{threads}"
                    command = [
                        str(binary),
                        *backend_arguments(args.candidate_backend),
                        "--count-fingerprint",
                        "-k",
                        str(protocol["k"]),
                        "--threads",
                        str(threads),
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
                        sample["input"]["path"],
                    ]
                    metrics, stats = invoke(command, work, protocol["timeout_seconds"])
                    for field, value in sample["expected_counts"].items():
                        assert stats[field] == value, (
                            accession,
                            hasher,
                            threads,
                            field,
                        )
                    assert stats["count_fingerprint"] == expected[accession]
                    assert stats["count_fingerprint"]["occurrences"] < (1 << 32)
                    observed = records(work / "outputs")
                    assert observed == inventories[accession]
                    targets = expected_records(sample)
                    assert [
                        {key: row[key] for key in target}
                        for row, target in zip(observed, targets, strict=True)
                    ] == targets
                    metadata = stats["count_store"]
                    verify_backend(metadata, args.candidate_backend, threads, hasher)
                    assert (
                        metadata["count_width_bits"] == 64
                        and metadata["occurrences"] == stats["n_kmers"]
                    )
                    assert (
                        metadata["owned_table_bytes"]
                        <= metadata["peak_table_bytes"]
                        <= metadata["peak_count_input_reserved_bytes"]
                    )
                    selected_digest = metadata["selected_input_sha256"]
                    if expected_digest is None:
                        expected_digest = selected_digest
                    assert selected_digest == expected_digest
                    receipt["runs"].append(
                        {
                            "accession": accession,
                            "hasher": hasher,
                            "threads": threads,
                            "command": command,
                            "metrics": metrics,
                            "stages": stats["stage_timings"],
                            "count_store": metadata,
                            "count_fingerprint": stats["count_fingerprint"],
                            "records": observed,
                        }
                    )
                    save(output / "matrix.json", receipt)
        assert len(receipt["runs"]) == 80
        verify_inputs(protocol)
        assert digest(args.protocol) == receipt["protocol_sha256"]
        assert digest(args.reference) == receipt["reference_sha256"]
        assert digest(Path(__file__)) == receipt["harness_sha256"]
        assert (
            digest(Path(__file__).with_name("compare_count_backends.py"))
            == receipt["comparison_helpers_sha256"]
        )
        for hasher, binary in binaries.items():
            assert digest(binary) == receipt["binaries"][hasher]["sha256"]
        receipt["status"] = "passed"
    except BaseException as error:
        receipt["status"] = "failed"
        receipt["error"] = repr(error)
        raise
    finally:
        receipt["finished_at_utc"] = datetime.datetime.now(datetime.UTC).isoformat()
        save(output / "matrix.json", receipt)


if __name__ == "__main__":
    main()

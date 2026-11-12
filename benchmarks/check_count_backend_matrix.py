#!/usr/bin/env python3
"""Check full count/product equivalence across hashers and thread counts."""

import argparse
import datetime
import json
import os
from pathlib import Path

from compare_count_backends import (
    digest,
    expected_records,
    invoke,
    records,
    save,
    verify_inputs,
)


def main() -> None:
    """Run the frozen calibration inputs through every supported execution mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--ahash", type=Path, required=True)
    parser.add_argument("--fxhash", type=Path, required=True)
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
                        "--count-backend",
                        "memory",
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
                    assert (
                        metadata["backend"] == "hashbrown_memory"
                        and metadata["hash_algorithm"] == hasher
                    )
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

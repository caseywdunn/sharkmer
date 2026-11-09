#!/usr/bin/env python3
"""Black-box fail-closed CLI checks using only generated FASTQ fixtures."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "scripts/audit_its_recruitment.py"
MANIFEST = Path(__file__).resolve().parents[1] / "inputs/manifest.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    source_hash = sha(SOURCE.read_bytes())
    manifest = json.loads(MANIFEST.read_bytes())
    sequence = manifest["regions"][0]["classes"][0]["sequence"].encode()
    valid_raw = b"@duplicate\n" + sequence + b"\n+\n" + b"I" * len(sequence) + b"\n"
    cases = []
    with tempfile.TemporaryDirectory(prefix="independent-its-cli-") as directory:
        base = Path(directory)
        for label in ("duplicate_ids", "wrong_hash", "truncated", "malformed_quality", "read_length_bound"):
            raw = valid_raw * 2 if label == "duplicate_ids" else valid_raw
            record_count = 2 if label in ("duplicate_ids", "truncated") else 1
            bases = len(sequence) * record_count
            if label == "malformed_quality":
                raw = b"@bad\n" + sequence + b"\n+\nI\n"
            if label == "read_length_bound":
                raw = b"@long\n" + b"A" * 1001 + b"\n+\n" + b"I" * 1001 + b"\n"
                bases = 1001
            fastq = base / (label + ".fastq")
            fastq.write_bytes(raw)
            expected = {"records": record_count, "bases": bases, "size_bytes": len(raw),
                        "sha256": "0" * 64 if label == "wrong_hash" else sha(raw)}
            identity = base / (label + ".identity.json")
            identity.write_text(json.dumps(expected))
            output = base / (label + ".result.json")
            argv = [sys.executable, str(SOURCE), "--fastq", str(fastq), "--manifest", str(MANIFEST),
                    "--expected-prefix", str(identity), "--output", str(output)]
            process = subprocess.run(argv, capture_output=True, text=True)
            failure = output.with_name(output.name + ".failure.json")
            if label == "duplicate_ids":
                assert process.returncode == 0, process.stderr
                result = json.loads(output.read_text())
                assert result["status"] == "completed" and not failure.exists()
                own = next(region for region in result["regions"] if region["region_id"] == "region_1")
                assert own["unique_records"] == 2 and own["distinct_recruited_read_ids"] == 1
                assert own["classes"]["region_1_class_1"]["distinct_read_ids"] == 1
                assert result["staged_file_before"] == result["staged_file_after"]
                assert result["staged_file_before"]["sha256"] == sha(raw)
                assert result["provenance"]["scanner_sha256"] == source_hash
                assert result["provenance"]["argv"] == argv[1:]
                before = sha(output.read_bytes())
                second = subprocess.run(argv, capture_output=True, text=True)
                assert second.returncode != 0 and sha(output.read_bytes()) == before
                assert not failure.exists()
                cases.append({"case": "existing_output_refused", "status": "passed"})
            else:
                assert process.returncode != 0 and not output.exists() and failure.exists()
                receipt = json.loads(failure.read_text())
                assert receipt["status"] == "failed"
                assert receipt["provenance"]["scanner_sha256"] == source_hash
                assert receipt["provenance"]["staged_file_before"]["sha256"] == sha(raw)
            cases.append({"case": label, "status": "passed", "exit_status": process.returncode})
    assert sha(SOURCE.read_bytes()) == source_hash, "Source changed during review"
    print(json.dumps({"status": "passed", "scanner_sha256": source_hash,
                      "biological_reads_examined": 0, "checks": cases}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

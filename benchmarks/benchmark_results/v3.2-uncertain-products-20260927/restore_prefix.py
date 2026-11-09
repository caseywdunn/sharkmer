#!/usr/bin/env python3
"""Bounded ENA R1 prefix restoration with exact archived-byte verification."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MAX_LINE_BYTES = 4096
MAX_COMPRESSED_BYTES = 512 * 1024 * 1024
TIMEOUT_SECONDS = 1800


def sha(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


class BoundedResponse:
    def __init__(self, response, started):
        self.response, self.started = response, started
        self.bytes = 0

    def read(self, size=-1):
        if size < 0:
            raise ValueError("Unbounded compressed-source read refused")
        if time.monotonic() - self.started > TIMEOUT_SECONDS:
            raise TimeoutError("Prefix restoration exceeded frozen elapsed-time bound")
        if self.bytes + size > MAX_COMPRESSED_BYTES:
            raise ValueError("Compressed-source byte bound exceeded")
        raw = self.response.read(size)
        self.bytes += len(raw)
        return raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", choices=("SRR27962769", "SRR31887760", "SRR1057608"), required=True)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    manifest_path = HERE / "inputs/read_inputs.json"
    manifest_before, source_before = sha(manifest_path), sha(Path(__file__))
    manifest = json.loads(manifest_path.read_text())
    sample = next(item for item in manifest["samples"] if item["accession"] == arguments.accession)
    expected = sample["expected_prefix"]
    path = arguments.output or ROOT / sample["path"]
    partial = path.with_name(path.name + ".partial")
    receipt = path.with_name(path.name + ".restore.json")
    if any(item.exists() or item.is_symlink() for item in (path, partial, receipt)):
        raise ValueError("Restore output, partial or receipt exists; choose a fresh attempt path")
    path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    records = bases = size = 0
    digest = hashlib.sha256()
    transport = {}
    try:
        request = urllib.request.Request(sample["source_url"], headers={"User-Agent": "sharkmer-frozen-prefix-restoration/1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            transport = {"url": response.url, "status": response.status, "headers": dict(response.headers)}
            bounded = BoundedResponse(response, started)
            with gzip.GzipFile(fileobj=bounded) as stream, partial.open("xb") as output:
                for records in range(1, expected["records"] + 1):
                    lines = [stream.readline(MAX_LINE_BYTES + 1) for _ in range(4)]
                    if any(not line.endswith(b"\n") or len(line) > MAX_LINE_BYTES for line in lines):
                        raise ValueError("Malformed, truncated or overlong remote FASTQ record")
                    header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
                    if not header.startswith(b"@") or not separator.startswith(b"+") or not sequence or len(sequence) != len(quality):
                        raise ValueError("Malformed remote FASTQ structure")
                    raw = b"".join(lines)
                    size += len(raw)
                    bases += len(sequence)
                    if size > expected["size_bytes"] or bases > expected["bases"]:
                        raise ValueError("Remote prefix exceeds frozen byte/base counts")
                    output.write(raw)
                    digest.update(raw)
                output.flush()
                os.fsync(output.fileno())
            transport["compressed_bytes_consumed"] = bounded.bytes
        observed = {"records": records, "bases": bases, "size_bytes": size, "sha256": digest.hexdigest()}
        if observed != expected:
            raise ValueError(f"Restored prefix differs from archive: {observed}")
        if sha(manifest_path) != manifest_before or sha(Path(__file__)) != source_before:
            raise ValueError("Restoration source/manifest changed during download")
        os.link(partial, path)
        partial.unlink()
        outcome = {"status": "verified", "accession": arguments.accession, "expected": expected,
                   "observed": observed, "output": str(path.resolve()), "transport": transport,
                   "source_sha256": source_before, "input_manifest_sha256": manifest_before,
                   "elapsed_seconds": time.monotonic() - started,
                   "limits": {"line_bytes": MAX_LINE_BYTES, "compressed_bytes": MAX_COMPRESSED_BYTES, "elapsed_seconds": TIMEOUT_SECONDS},
                   "scope": "Only exact consumed first1M R1; complete compressed upstream checksum and former longer staged-file identity are not verified"}
    except Exception as error:
        outcome = {"status": "failed", "accession": arguments.accession, "error_type": type(error).__name__,
                   "error": str(error), "records": records, "bases": bases, "size_bytes": size,
                   "transport": transport, "source_sha256": source_before,
                   "input_manifest_sha256": manifest_before, "elapsed_seconds": time.monotonic() - started}
        with receipt.open("x") as output:
            json.dump(outcome, output, indent=2, sort_keys=True)
            output.write("\n")
        raise
    with receipt.open("x") as output:
        json.dump(outcome, output, indent=2, sort_keys=True)
        output.write("\n")
    print(json.dumps(outcome, sort_keys=True))


if __name__ == "__main__":
    main()

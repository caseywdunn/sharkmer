#!/usr/bin/env python3
"""Restore the byte-identical archived R1 prefix, without downloading unused R2."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request


URL = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR279/069/SRR27962769/SRR27962769_1.fastq.gz"
EXPECTED = {"records": 1000000, "bases": 150000000, "size_bytes": 365356898,
            "sha256": "6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    path = args.output
    partial = path.with_suffix(path.suffix + ".partial")
    receipt = path.with_suffix(path.suffix + ".restore.json")
    if any(p.exists() or p.is_symlink() for p in (path, partial, receipt)):
        raise ValueError("Restore paths already exist; choose a fresh destination")
    path.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    digest = hashlib.sha256()
    size = bases = records = 0
    try:
        with urllib.request.urlopen(URL, timeout=60) as response:
            with gzip.GzipFile(fileobj=response) as stream, partial.open("xb") as output:
                for records in range(1, EXPECTED["records"] + 1):
                    lines = [stream.readline() for _ in range(4)]
                    if any(not line for line in lines):
                        raise ValueError("Truncated source FASTQ")
                    header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
                    if not header.startswith(b"@") or not separator.startswith(b"+") or len(sequence) != len(quality):
                        raise ValueError("Malformed source FASTQ")
                    content = b"".join(lines)
                    output.write(content)
                    digest.update(content)
                    size += len(content)
                    bases += len(sequence)
                    if time.time() - started > 1800:
                        raise TimeoutError("Restore exceeded 1800 seconds")
        observed = {"records": records, "bases": bases, "size_bytes": size, "sha256": digest.hexdigest()}
        if observed != EXPECTED:
            raise ValueError(f"Restored prefix identity differs: {observed}")
        # A hard link publishes without overwriting any concurrently created file.
        os.link(partial, path)
        partial.unlink()
        outcome = {"status": "verified", "url": URL, "output": str(path.resolve()),
                   "expected": EXPECTED, "observed": observed, "elapsed_seconds": time.time() - started,
                   "scope": "First 1M R1 records only; not original 8M-record staged file",
                   "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        with receipt.open("x") as handle:
            json.dump(outcome, handle, indent=2)
            handle.write("\n")
        print(json.dumps(outcome))
    except Exception as error:
        with receipt.open("x") as handle:
            json.dump({"status": "failed", "url": URL, "error": str(error),
                       "records": records, "size_bytes": size,
                       "elapsed_seconds": time.time() - started}, handle, indent=2)
            handle.write("\n")
        raise


if __name__ == "__main__":
    main()

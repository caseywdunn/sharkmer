#!/usr/bin/env python3
"""Restore exact archived FASTQ prefixes, retaining immutable attempt receipts.

Samples consume the archived ordered sources: complete preceding files, then a
bounded prefix or complete final file. The human sample requires three complete
files; Liriodendron requires complete R1 plus partial R2. Verification runs on a compute node,
outside benchmark timing. No assembler or reference classifier is invoked.

Version 2 correction: the archived compressed_sha256 field fingerprints a local
historical cache gzip, not the upstream ENA representation. Complete upstream
streams must match archived ENA MD5 and size; their newly observed SHA256 is
recorded separately. Exact final raw-prefix SHA256 and counts remain mandatory.
"""

import argparse
from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import stat
import sys
import time
import urllib.parse
import urllib.request


HERE = Path(__file__).resolve().parent
MAX_LINE_BYTES = 4096
MAX_COMPRESSED_BYTES = 512 * 1024 * 1024  # per source, including read-ahead
MAX_PREFIX_BYTES = 512 * 1024 * 1024
MAX_RECORDS = 1_000_000
TIMEOUT_SECONDS = 1800  # each sample, enforced by SIGALRM
NETWORK_TIMEOUT_SECONDS = 60
MAX_MANIFEST_BYTES = 2 * 1024 * 1024


def identity(path):
    value = Path(path).lstat()
    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"Expected regular non-symlink file: {path}")
    return {key: getattr(value, "st_" + key) for key in
            ("dev", "ino", "size", "mtime_ns", "ctime_ns")}


def fingerprint(path):
    path = Path(path).absolute()
    before = identity(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    if identity(path) != before:
        raise ValueError(f"File changed during hashing: {path}")
    return {"path": str(path), "size_bytes": before["size"], "sha256": digest.hexdigest()}


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


@contextmanager
def sample_deadline():
    def expired(_signum, _frame):
        raise TimeoutError("Sample exceeded 1800-second restoration/verification limit")
    prior = signal.signal(signal.SIGALRM, expired)
    signal.alarm(TIMEOUT_SECONDS)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, prior)


class BoundedResponse:
    def __init__(self, response):
        self.response = response
        self.bytes = 0
        self.sha256 = hashlib.sha256()
        self.md5 = hashlib.md5()

    def read(self, size=-1):
        if size < 0:
            raise ValueError("Unbounded compressed read refused")
        raw = self.response.read(min(size, MAX_COMPRESSED_BYTES - self.bytes + 1))
        self.bytes += len(raw)
        self.sha256.update(raw)
        self.md5.update(raw)
        if self.bytes > MAX_COMPRESSED_BYTES:
            raise ValueError("Compressed source byte limit exceeded")
        return raw

    def observation(self):
        return {"size_bytes": self.bytes, "sha256": self.sha256.hexdigest(),
                "md5": self.md5.hexdigest()}


class PrefixCounter:
    def __init__(self, expected):
        self.expected = expected
        self.records = self.bases = self.size_bytes = 0
        self.digest = hashlib.sha256()

    def accept(self, lines):
        if any(not line.endswith(b"\n") or len(line) > MAX_LINE_BYTES for line in lines):
            raise ValueError("Truncated, non-newline-terminated or overlong FASTQ record")
        header, sequence, separator, quality = [line[:-1].removesuffix(b"\r") for line in lines]
        if (not header.startswith(b"@") or not separator.startswith(b"+") or
                not sequence or len(sequence) != len(quality) or
                any(char < 33 or char > 126 for char in quality)):
            raise ValueError("Malformed FASTQ record")
        raw = b"".join(lines)
        self.records += 1
        self.bases += len(sequence)
        self.size_bytes += len(raw)
        if (self.records > self.expected["records"] or self.bases > self.expected["bases"] or
                self.size_bytes > self.expected["size_bytes"]):
            raise ValueError("FASTQ exceeds archived record/base/byte counts")
        self.digest.update(raw)
        return raw

    def observation(self):
        return {"records": self.records, "bases": self.bases,
                "size_bytes": self.size_bytes, "sha256": self.digest.hexdigest()}

    def require_exact(self):
        if self.observation() != self.expected:
            raise ValueError(f"FASTQ differs from archived prefix: {self.observation()}")


def consume(stream, records, counter, output=None, require_eof=False):
    for _ in range(records):
        raw = counter.accept([stream.readline(MAX_LINE_BYTES + 1) for _ in range(4)])
        if output is not None:
            output.write(raw)
    if require_eof and stream.read(1) != b"":
        raise ValueError("Source has additional records/bytes beyond required EOF")


def verify_existing(path, expected):
    before = identity(path)
    if before["size"] != expected["size_bytes"]:
        raise ValueError("Existing FASTQ has wrong byte count")
    counter = PrefixCounter(expected)
    with Path(path).open("rb") as stream:
        consume(stream, expected["records"], counter, require_eof=True)
    if identity(path) != before:
        raise ValueError("Existing FASTQ changed during verification")
    counter.require_exact()
    return counter.observation(), before


def restore(sample, partial, transports, opener=urllib.request.urlopen):
    counter = PrefixCounter(sample["expected_prefix"])
    with partial.open("xb") as output:
        for source in sample["restoration_sources"]:
            transport = {"requested_url": source["url"], "records_required": source["records"],
                         "eof_required": source["eof_required"], "eof_verified": False,
                         "upstream_md5_and_size_verified": False,
                         "historical_local_cache_sha256": source["archived_compressed_sha256"],
                         "historical_cache_sha256_is_upstream_pin": False}
            transports.append(transport)
            request = urllib.request.Request(source["url"], headers={
                "User-Agent": "sharkmer-archived-prefix-restoration/3"})
            with opener(request, timeout=NETWORK_TIMEOUT_SECONDS) as response:
                transport.update({"response_url": response.url, "http_status": response.status,
                                  "headers": {key: response.headers.get(key) for key in
                                              ("Content-Length", "ETag", "Last-Modified")}})
                if response.status != 200:
                    raise ValueError("Expected complete HTTP 200 source response")
                bounded = BoundedResponse(response)
                try:
                    with gzip.GzipFile(fileobj=bounded) as stream:
                        consume(stream, source["records"], counter, output,
                                require_eof=source["eof_required"])
                    transport["eof_verified"] = source["eof_required"]
                    if source["eof_required"]:
                        observed = bounded.observation()
                        required = {"size_bytes": source["archived_full_compressed_size_bytes"],
                                    "md5": source["archived_full_compressed_md5"]}
                        if any(observed[key] != value for key, value in required.items()):
                            raise ValueError("Complete compressed source differs from archived ENA MD5/size")
                        transport["upstream_md5_and_size_verified"] = True
                finally:
                    transport["compressed_bytes_read_including_read_ahead"] = bounded.observation()
        output.flush()
        os.fsync(output.fileno())
    counter.require_exact()
    return counter.observation()


def validate_inventory(document):
    samples = document["samples"]
    if document.get("schema_version") != 1 or not isinstance(samples, list) or not 1 <= len(samples) <= 10:
        raise ValueError("Unsupported inventory schema/sample count")
    seen = set()
    panels = {row["name"] for row in document["panels"]}
    for sample in samples:
        accession = sample["accession"]
        if not re.fullmatch(r"(?:SRR|ERR)[0-9]+", accession) or accession in seen:
            raise ValueError("Invalid or duplicate accession")
        seen.add(accession)
        if sample["panel"] not in panels:
            raise ValueError("Unknown assembly panel")
        expected = sample["expected_prefix"]
        if set(expected) != {"records", "bases", "size_bytes", "sha256"}:
            raise ValueError("Unexpected expected-prefix schema")
        for key, limit in (("records", MAX_RECORDS), ("bases", MAX_PREFIX_BYTES),
                           ("size_bytes", MAX_PREFIX_BYTES)):
            if type(expected[key]) is not int or not 0 < expected[key] <= limit:
                raise ValueError("Invalid or excessive archived count")
        if not re.fullmatch(r"[0-9a-f]{64}", expected["sha256"]):
            raise ValueError("Invalid archived SHA256")
        sources = sample["restoration_sources"]
        if not 1 <= len(sources) <= 3:
            raise ValueError("Invalid number of ordered sources")
        for source_index, source in enumerate(sources):
            parsed = urllib.parse.urlsplit(source["url"])
            if (parsed.scheme != "https" or parsed.netloc != "ftp.sra.ebi.ac.uk" or
                    not parsed.path.startswith("/vol1/fastq/") or parsed.query or parsed.fragment):
                raise ValueError("Unexpected ENA source URL")
            if type(source["records"]) is not int or not 0 < source["records"] <= MAX_RECORDS:
                raise ValueError("Invalid source record count")
            if type(source["eof_required"]) is not bool:
                raise ValueError("Invalid source EOF requirement")
            if source_index < len(sources) - 1 and not source["eof_required"]:
                raise ValueError("Concatenated preceding sources must require complete EOF")
            if source["eof_required"] and (
                    not re.fullmatch(r"[0-9a-f]{64}", source["archived_compressed_sha256"] or "") or
                    not re.fullmatch(r"[0-9a-f]{32}", source["archived_full_compressed_md5"] or "") or
                    type(source["archived_full_compressed_size_bytes"]) is not int or
                    not 0 < source["archived_full_compressed_size_bytes"] <= MAX_COMPRESSED_BYTES):
                raise ValueError("Complete source lacks bounded archived compressed identity")
        if sum(source["records"] for source in sources) != expected["records"]:
            raise ValueError("Ordered source record counts disagree with prefix")


def restore_one(sample, output_dir, frozen_sources, opener=urllib.request.urlopen):
    accession = sample["accession"]
    receipt_path = output_dir / f"{accession}.restore.json"
    final = output_dir / f"{accession}.fastq"
    partial = output_dir / f"{accession}.fastq.partial"
    if any(path.exists() or path.is_symlink() for path in (receipt_path, final, partial)):
        raise ValueError("Attempt output exists; choose a fresh attempt directory")
    started = time.monotonic()
    receipt = {"schema_version": 1, "status": "failed", "accession": accession,
               "panel_id": sample["panel"], "expected_prefix": sample["expected_prefix"],
               "sources_before": frozen_sources, "transport": [], "available_path": None}
    try:
        with sample_deadline():
            current = sample.get("current_file")
            if current:
                available = Path(current["path"]).absolute()
                receipt["mode"] = "verify_existing_exact_prefix"
                observed, file_identity = verify_existing(available, sample["expected_prefix"])
            else:
                receipt["mode"] = "restore_ordered_sources"
                observed = restore(sample, partial, receipt["transport"], opener=opener)
                available = final.absolute()
                file_identity = identity(partial)
            receipt["observed_prefix"] = observed
            receipt["sources_after"] = {key: fingerprint(row["path"]) for key, row in frozen_sources.items()}
            if receipt["sources_after"] != frozen_sources:
                raise ValueError("Restoration source/inventory changed")
            if not current:
                os.link(partial, final)  # no-clobber publication on the same filesystem
                partial.unlink()
                file_identity = identity(final)
            receipt.update(status="verified", available_path=str(available),
                           available_file_identity=file_identity)
    except Exception as error:
        receipt.update(error_type=type(error).__name__, error=str(error))
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        receipt["limits"] = {"elapsed_seconds_per_sample": TIMEOUT_SECONDS,
                             "compressed_bytes_per_source": MAX_COMPRESSED_BYTES,
                             "line_bytes": MAX_LINE_BYTES, "prefix_bytes": MAX_PREFIX_BYTES}
        write_json(receipt_path, receipt)
    return receipt


def run(inventory, output, expected_inventory_sha256, accessions=None):
    output = Path(output).absolute()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    report = {"schema_version": 1, "status": "failed", "samples": [],
              "argv": sys.argv, "python": sys.version, "python_executable": sys.executable,
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
              "scope": "Exact archived consumed FASTQ bytes only; no assembly or biological classification"}
    try:
        frozen = {"inventory": fingerprint(inventory), "helper": fingerprint(__file__)}
        report["sources_before"] = frozen
        if frozen["inventory"]["size_bytes"] > MAX_MANIFEST_BYTES:
            raise ValueError("Inventory exceeds byte limit")
        if frozen["inventory"]["sha256"] != expected_inventory_sha256:
            raise ValueError("Inventory differs from external SHA256 pin")
        document = json.loads(Path(inventory).read_text())
        validate_inventory(document)
        if accessions and (len(set(accessions)) != len(accessions) or
                           set(accessions) - {row["accession"] for row in document["samples"]}):
            raise ValueError("Unknown or repeated requested accession")
        selected = [sample for sample in document["samples"]
                    if not accessions or sample["accession"] in accessions]
        report["complete_inventory_selected"] = len(selected) == len(document["samples"])
        for panel in document["panels"]:
            if fingerprint(panel["assembly_panel"]["path"]) != panel["assembly_panel"]:
                raise ValueError("Assembly panel differs from materialized inventory")
        for sample in selected:
            receipt = restore_one(sample, output, frozen)
            receipt_path = output / f"{sample['accession']}.restore.json"
            report["samples"].append({key: receipt[key] for key in
                                      ("accession", "panel_id", "status", "available_path", "expected_prefix")})
            report["samples"][-1]["receipt"] = fingerprint(receipt_path)
            print(json.dumps(report["samples"][-1]), flush=True)
        report["sources_after"] = {key: fingerprint(row["path"]) for key, row in frozen.items()}
        if report["sources_after"] != frozen:
            raise ValueError("Restoration source/inventory changed")
        if any(sample["status"] != "verified" for sample in report["samples"]):
            raise ValueError("One or more samples failed; preserve this attempt and retry in a fresh directory")
        report["status"] = "pass"
    except Exception as error:
        report.update(error_type=type(error).__name__, error=str(error))
    finally:
        report["elapsed_seconds"] = time.monotonic() - started
        write_json(output / "availability.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "inventory.json")
    parser.add_argument("--expected-inventory-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--accession", action="append")
    args = parser.parse_args()
    result = run(args.inventory, args.output, args.expected_inventory_sha256, args.accession)
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())

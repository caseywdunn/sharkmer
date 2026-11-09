#!/usr/bin/env python3
"""Bounded HTTPS panel/ingestion parity and certificate rejection using Sharkmer.

Run with Python + PyYAML, openssl, outbound GitHub HTTPS, and loopback sockets:
  python https_check.py --binary /absolute/path/sharkmer --output /new/results

This tests the ureq TLS path used for HTTPS panels. ENA metadata also uses
ureq HTTPS, but ordinary ENA fastq_ftp URLs become HTTP in current Sharkmer;
this check makes no claim of HTTPS FASTQ/cache or biological validation.
"""

import argparse
import datetime
import hashlib
import http.server
import json
import os
from pathlib import Path
import re
import shutil
import ssl
import subprocess
import tempfile
import threading
import time
import urllib.request

import yaml


COMMIT = "e914ea728fc5ffc253273560e3de2cfc9395893b"
URL = (
    f"https://raw.githubusercontent.com/caseywdunn/sharkmer/{COMMIT}/"
    "tests/fixtures/test_panel.yaml"
)
PANEL_SHA256 = "6bb5785c055666c0c58cd440315e8d98711310da20f47774ea3aa9f74ac983cb"
READS = 3
K = 19


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-version", default="3.2.1")
    args = parser.parse_args()
    binary, output = args.binary.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {
        "status": "failed",
        "started_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "binary": str(binary),
        "binary_sha256": sha256(binary.read_bytes()),
        "script_sha256": sha256(Path(__file__).read_bytes()),
        "panel_url": URL,
        "panel_sha256": PANEL_SHA256,
        "read_limit": READS,
        "commands": [],
        "scope": "HTTPS panel loading, bounded local ingestion, and untrusted-certificate rejection; no ENA FASTQ/cache or biological validation",
    }

    def execute(label, command, timeout=45, environment=None):
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command, capture_output=True, text=True, timeout=timeout, env=environment
            )
        except subprocess.TimeoutExpired as error:
            for kind in ("stdout", "stderr"):
                contents = getattr(error, kind) or b""
                if isinstance(contents, bytes):
                    contents = contents.decode(errors="replace")
                (output / f"{label}.{kind}.log").write_text(contents)
            receipt["commands"].append({
                "label": label, "command": command, "timed_out": True,
                "timeout_seconds": timeout,
            })
            raise
        for kind in ("stdout", "stderr"):
            (output / f"{label}.{kind}.log").write_text(getattr(completed, kind))
        receipt["commands"].append({
            "label": label,
            "command": command,
            "returncode": completed.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "stdout_sha256": sha256(completed.stdout.encode()),
            "stderr_sha256": sha256(completed.stderr.encode()),
        })
        return completed

    try:
        require(shutil.which("openssl") is not None, "openssl is required")
        version = execute("version", [str(binary), "--version"], timeout=10)
        require(version.returncode == 0, version.stderr)
        require(version.stdout.strip() == f"sharkmer {args.expected_version} (https://github.com/caseywdunn/sharkmer)", "Unexpected binary version")
        receipt["version_output"] = version.stdout.strip()

        # Independently bind the tiny immutable endpoint; Sharkmer then makes its
        # own request below, using its unchanged certificate verification.
        with urllib.request.urlopen(URL, timeout=15, context=ssl.create_default_context()) as response:
            require(response.status == 200 and response.geturl() == URL, "Unexpected panel HTTP response/redirect")
            panel = response.read(4097)
        require(len(panel) <= 4096 and sha256(panel) == PANEL_SHA256, "Pinned HTTPS panel bytes differ")
        local_panel = output / "pinned-panel.yaml"
        local_panel.write_bytes(panel)

        bases = b"".join(hashlib.sha256(f"https-probe-{i}".encode()).digest() for i in range(4))
        sequence = "".join("ACGT"[value % 4] for value in bases[:100])
        reverse = lambda value: value.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        kmers = {min(sequence[i:i + K], reverse(sequence[i:i + K])) for i in range(len(sequence) - K + 1)}
        require(len(kmers) == 82, "Synthetic read unexpectedly repeats a canonical k-mer")
        fastq = output / "four-reads.fastq"
        fastq.write_text("".join(f"@probe_{i}\n{sequence}\n+\n{'I' * 100}\n" for i in range(4)))
        receipt["fastq_sha256"] = sha256(fastq.read_bytes())
        counts = {
            "n_reads_read": 3, "n_bases_read": 300,
            "n_subreads_ingested": 3, "n_bases_ingested": 300, "n_kmers": 246,
        }
        receipt["expected_counts"] = counts

        def command(source, directory):
            return [
                str(binary), "-vv", "--color", "never", "-k", str(K), "-t", "1",
                "--chunks", "0", "--max-reads", str(READS), "-s", "https_probe",
                "--cache-dir", str(output / "read-cache"), "-o", str(directory),
                "--pcr-panel-file", source, str(fastq),
            ]

        results = {}
        for label, source in (("local", str(local_panel)), ("https", URL)):
            directory = output / label
            completed = execute(label, command(source, directory))
            require(completed.returncode == 0, f"{label} run failed: {completed.stderr}")
            stats_path = directory / "https_probe.stats.yaml"
            manifest_path = directory / "https_probe.manifest.yaml"
            stats = yaml.safe_load(stats_path.read_text())
            manifest = yaml.safe_load(manifest_path.read_text())
            require(stats["run_status"] == manifest["status"] == "complete", "Incomplete run")
            require(stats["run_id"] == manifest["run_id"], "Run identity mismatch")
            require({key: stats[key] for key in counts} == counts, "Exact ingestion/count oracle failed")
            require(stats["input_source"]["kind"] == "local_files" and stats["input_source"]["max_reads"] == READS, "Input provenance differs")
            genes = stats["pcr_results"]
            require(len(genes) == 1 and genes[0]["gene_name"] == "test_panel_18S", "Pinned panel target differs")
            require(genes[0]["n_products"] == 0 and not list(directory.glob("*.fasta")), "Unexpected product from 100-base negative input")
            require(len(manifest["files"]) == 1 and manifest["files"][0]["path"] == stats_path.name, "Unexpected output inventory")
            require(manifest["files"][0]["sha256"] == sha256(stats_path.read_bytes()), "Output receipt checksum differs")
            results[label] = genes
            receipt[f"{label}_stats_sha256"] = sha256(stats_path.read_bytes())
            receipt[f"{label}_manifest_sha256"] = sha256(manifest_path.read_bytes())
            if label == "https":
                require(re.search(r"Using ciphersuite TLS13_", completed.stderr) is not None, "Successful Sharkmer request did not record a TLS 1.3 ciphersuite")
                receipt["sharkmer_tls13_observed"] = True
        require(results["local"] == results["https"], "Local/HTTPS panel outcomes differ")
        require(not (output / "read-cache").exists(), "Local-input check unexpectedly created a read cache")

        requests = []

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                requests.append(self.path)
                self.send_response(200)
                self.end_headers()
                self.wfile.write(panel)

            def log_message(self, *_args):
                pass

        with tempfile.TemporaryDirectory(prefix="tls-key-", dir=output) as temporary:
            key, cert = Path(temporary) / "key.pem", output / "untrusted-cert.pem"
            generated = execute("certificate", [
                "openssl", "req", "-x509", "-newkey", "rsa:2048", "-nodes",
                "-keyout", str(key), "-out", str(cert), "-days", "1",
                "-subj", "/CN=localhost", "-addext", "subjectAltName=DNS:localhost,IP:127.0.0.1",
                "-addext", "basicConstraints=critical,CA:FALSE",
                "-addext", "extendedKeyUsage=serverAuth",
            ], timeout=15)
            require(generated.returncode == 0, generated.stderr)
            context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
            context.minimum_version = ssl.TLSVersion.TLSv1_3
            context.load_cert_chain(cert, key)
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
            server.socket = context.wrap_socket(server.socket, server_side=True)
            thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.05}, daemon=True)
            thread.start()
            try:
                url = f"https://127.0.0.1:{server.server_port}/panel.yaml"
                environment = {key: value for key, value in os.environ.items() if key.lower() not in {"http_proxy", "https_proxy", "all_proxy"}}
                environment["NO_PROXY"] = "localhost,127.0.0.1"
                rejected = execute("untrusted", command(url, output / "untrusted"), timeout=15, environment=environment)
                require(rejected.returncode != 0, "Self-signed certificate was accepted")
                require(re.search(r"UnknownIssuer|unknown issuer|self[- ]signed|unknown certificate issuer", rejected.stderr, re.I) is not None, "Failure was not specifically certificate trust rejection")
                require(not requests and not (output / "untrusted").exists(), "Untrusted server reached HTTP or produced output")
                receipt["certificate_rejection"] = {"http_requests": len(requests), "certificate_sha256": sha256(cert.read_bytes())}
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=2)
        require(sha256(binary.read_bytes()) == receipt["binary_sha256"], "Binary changed during check")
        receipt["status"] = "pass"
    except Exception as error:
        receipt["error"] = str(error)
        raise
    finally:
        receipt["finished_at_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "pass", "receipt": str(output / "receipt.json")}))


if __name__ == "__main__":
    main()

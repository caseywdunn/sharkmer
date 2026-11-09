#!/usr/bin/env python3
"""Run one immutable bounded assay and attest input specification/source/output."""

import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from scan_exact import sha


HERE = Path(__file__).resolve().parent


def fingerprint(path):
    return {"sha256": sha(path.read_bytes()), "bytes": path.stat().st_size}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if len(os.sched_getaffinity(0)) != 1 or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("A one-CPU compute allocation is required")
    manifest = json.loads(args.manifest.read_text())
    sources = {Path(p): receipt for p, receipt in manifest["source_files"].items()}
    for path, receipt in sources.items():
        if fingerprint(path) != receipt:
            raise ValueError("Frozen derivation source differs")
    for p in (Path(__file__), HERE / "scan_exact.py", args.manifest):
        sources[p.resolve()] = fingerprint(p)
    args.output.mkdir(parents=True, exist_ok=False)
    command = ["/usr/bin/timeout", "--kill-after=5s", "600s", "/usr/bin/time", "-o", str(args.output / "time.json"),
               "-f", '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}',
               "/usr/bin/prlimit", "--as=2147483648", "--", sys.executable, str(HERE / "scan_exact.py"),
               "--manifest", str(args.manifest.resolve()), "--output", str((args.output / "raw.json").resolve())]
    start = time.monotonic()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (args.output / "stdout.log").open("x") as stdout, (args.output / "stderr.log").open("x") as stderr:
        completed = subprocess.run(command, stdout=stdout, stderr=stderr)
    after = {str(p): fingerprint(p) for p in sources}
    before = {str(p): receipt for p, receipt in sources.items()}
    status = "pass" if completed.returncode == 0 and before == after else "fail"
    receipt = {"status": status, "accession": manifest["accession"], "started_utc": started,
               "returncode": completed.returncode, "command": command, "elapsed_seconds": time.monotonic() - start,
               "source_files_before": before, "source_files_after": after, "sources_unchanged": before == after,
               "slurm_job_id": os.environ["SLURM_JOB_ID"], "cpu_affinity": sorted(os.sched_getaffinity(0)),
               "outputs": {p.name: fingerprint(p) for p in sorted(args.output.iterdir()) if p.is_file()}}
    with (args.output / "execution.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"accession": manifest["accession"], "status": status, "elapsed_seconds": receipt["elapsed_seconds"]}))
    if status != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

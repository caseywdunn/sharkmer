#!/usr/bin/env python3
"""Run one immutable, bounded exact-prefix Drosophila audit."""

import datetime
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = HERE / "run-v1"
    output.mkdir()
    controls = json.loads((HERE / "prefilter_controls.json").read_text())
    if controls["status"] != "pass" or len(os.sched_getaffinity(0)) != 1:
        raise ValueError("Passing prefilter controls and one-CPU affinity are required")
    derivation = json.loads((HERE / "derivation.json").read_text())
    files = {Path(p) for p in derivation["source_files"]}
    files.update(HERE / n for n in ("prepare_pilot.py", "scan_prefiltered.py", "check_prefilter.py",
                 "run_pilot.py", "manifest.json", "expected_prefix.json", "derivation.json", "prefilter_controls.json"))
    before = {str(p): digest(p) for p in sorted(files)}
    if before[str(HERE / "prepare_pilot.py")] != derivation["producer_sha256"]:
        raise ValueError("Preparation source changed")
    for p, value in derivation["source_files"].items():
        if before[p] != value["sha256"]:
            raise ValueError("Frozen source differs")
    for name, expected in controls["source_sha256"].items():
        if before[str(HERE / name)] != expected:
            raise ValueError("Prefilter control source differs")
    command = ["/usr/bin/timeout", "--kill-after=5s", "600s", "/usr/bin/time", "-o", str(output / "time.json"),
               "-f", '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}',
               "/usr/bin/prlimit", "--as=2147483648", "--", sys.executable,
               str(HERE / "scan_prefiltered.py"), "--fastq",
               str(ROOT / "tmp/uncertain-products-20260927/reads/SRR31887760.first1M_R1.retry1.fastq"),
               "--output", str(output / "raw.json")]
    start = time.monotonic()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (output / "stdout.log").open("x") as out, (output / "stderr.log").open("x") as err:
        completed = subprocess.run(command, stdout=out, stderr=err, check=False)
    after = {str(p): digest(p) for p in sorted(files)}
    receipt = {"status": "pass" if completed.returncode == 0 and before == after else "fail",
               "returncode": completed.returncode, "command": command, "started_utc": started,
               "observer_elapsed_seconds": time.monotonic() - start, "source_hashes_before": before,
               "source_hashes_after": after, "sources_unchanged": before == after,
               "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "cpu_affinity": sorted(os.sched_getaffinity(0)),
               "address_space_limit_bytes": 2147483648, "time_limit_seconds": 600,
               "outputs": {p.name: {"sha256": digest(p), "bytes": p.stat().st_size}
                           for p in sorted(output.iterdir()) if p.is_file()}}
    with (output / "execution.json").open("x") as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: receipt[k] for k in ("status", "returncode", "observer_elapsed_seconds", "sources_unchanged")}))
    if receipt["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

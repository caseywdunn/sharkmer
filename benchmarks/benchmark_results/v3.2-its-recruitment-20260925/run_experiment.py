#!/usr/bin/env python3
"""Freeze and execute the bounded same-input diagnostic with immutable receipts."""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = Path(__file__).resolve().parent
EXPECTED_SHA = "6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794"
SCANNER_SHA = "5c6268e840a626eb980713da8b0a1992030c50bad0e258241dcce703b3f0a75d"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def fingerprint(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular nonsymlink file: {path}")
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, field) != getattr(after, field) for field in fields):
        raise ValueError(f"File changed during fingerprinting: {path}")
    return {"size_bytes": after.st_size, "sha256": digest.hexdigest()}


def write_json(path, payload):
    path = Path(path)
    raw = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf8")
    descriptor, temporary = tempfile.mkstemp(prefix=".its-run-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(raw)
            output.flush()
            os.fsync(output.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def source_paths():
    paths = [ROOT / relative for relative in (
        "dev_docs/ITS_READ_RECRUITMENT_PROTOCOL.md", "scripts/audit_its_recruitment.py",
        "scripts/audit_read_bridges.py", "tests/test_its_recruitment.py", "tests/test_read_bridges.py")]
    paths.extend(EVIDENCE.glob("*.py"))
    for directory in ("inputs", "review", "restoration"):
        paths.extend(path for path in (EVIDENCE / directory).rglob("*")
                     if path.is_file() and "__pycache__" not in path.parts)
    return sorted(set(paths))


def freeze(fastq, run_dir, freeze_path):
    if freeze_path.exists() or run_dir.exists():
        raise ValueError("Freeze/output path exists; choose fresh paths")
    raw = fingerprint(fastq)
    if raw != {"size_bytes": 365356898, "sha256": EXPECTED_SHA}:
        raise ValueError("Input is not the exact archived R1 prefix")
    if fingerprint(ROOT / "scripts/audit_its_recruitment.py")["sha256"] != SCANNER_SHA:
        raise ValueError("Scanner differs from independently reviewed source")
    manifest = EVIDENCE / "inputs/manifest.json"
    expected = EVIDENCE / "inputs/originals/inputs/expected_prefix.json"
    control = EVIDENCE / "inputs/originals/inputs/co1_control_manifest.json"
    run_dir.mkdir(parents=True)
    commands = [{"id": "recruitment", "argv": [sys.executable,
        str(ROOT / "scripts/audit_its_recruitment.py"), "--fastq", str(fastq),
        "--manifest", str(manifest), "--expected-prefix", str(expected),
        "--output", str(run_dir / "recruitment.json")]}]
    for substitutions in (0, 1):
        identifier = f"co1_substitutions_{substitutions}"
        commands.append({"id": identifier, "argv": [sys.executable,
            str(ROOT / "scripts/audit_read_bridges.py"), "--fastq", str(fastq),
            "--events", str(control), "--expected-prefix", str(expected),
            "--max-substitutions", str(substitutions), "--output", str(run_dir / f"{identifier}.json")]})
    sources = {str(path.relative_to(ROOT)): fingerprint(path) for path in source_paths()}
    for relative in sources:
        snapshot = run_dir / "source_snapshot" / relative
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        with snapshot.open("xb") as output:
            output.write((ROOT / relative).read_bytes())
        if fingerprint(snapshot) != sources[relative]:
            raise ValueError("Source changed while snapshotting")
    plan = {"schema_version": 1, "frozen_utc": now(), "baseline": "956e771437a2410a7bcbb1504b033c61dd1e11b6",
            "branch": "experiment/its-read-threading", "python": sys.version,
            "root": str(ROOT), "run_dir": str(run_dir), "fastq": str(fastq),
            "input_identity": raw, "sources": sources, "commands": commands,
            "timeout_seconds_per_command": 1800, "address_space_bytes": 40 * 1024 ** 3,
            "statement": "Protocol, sources, inputs, synthetic review, and argv frozen before new biological scans"}
    write_json(freeze_path, plan)
    print(json.dumps({"freeze": str(freeze_path), **fingerprint(freeze_path), "sources": len(sources)}))


def verify(plan):
    if fingerprint(plan["fastq"]) != plan["input_identity"]:
        raise ValueError("Input changed from freeze")
    for relative, identity in plan["sources"].items():
        if fingerprint(ROOT / relative) != identity:
            raise ValueError(f"Source changed from freeze: {relative}")


def execute(freeze_path):
    freeze_identity = fingerprint(freeze_path)
    plan = json.loads(freeze_path.read_text())
    run_dir = Path(plan["run_dir"])
    write_json(run_dir / "started.json", {"started_utc": now(), "freeze": freeze_identity,
               "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "python": sys.version})
    failed = False
    outcomes = []
    for command in plan["commands"]:
        identifier = command["id"]
        started = time.monotonic()
        timing = run_dir / f"{identifier}.time.txt"
        argv = ["/usr/bin/time", "-v", "-o", str(timing), *command["argv"]]
        receipt = {"id": identifier, "argv": argv, "started_utc": now(), "freeze": freeze_identity,
                   "exit_code": None, "launched": False}
        def limits():
            resource.setrlimit(resource.RLIMIT_AS, (plan["address_space_bytes"], plan["address_space_bytes"]))
        process = None
        try:
            verify(plan)
            if fingerprint(freeze_path) != freeze_identity:
                raise ValueError("Freeze changed before command execution")
            if plan.get("python", sys.version) != sys.version:
                raise ValueError("Python version differs from freeze")
            if list(run_dir.glob(identifier + ".*")):
                raise ValueError("Refusing to overwrite existing command output, timing, or receipt")
            receipt["input_before"] = fingerprint(plan["fastq"])
            with (run_dir / f"{identifier}.stdout.txt").open("x") as stdout, (run_dir / f"{identifier}.stderr.txt").open("x") as stderr:
                process = subprocess.Popen(argv, cwd=ROOT, stdout=stdout, stderr=stderr,
                                           preexec_fn=limits, start_new_session=True)
                receipt["launched"] = True
                try:
                    receipt["exit_code"] = process.wait(timeout=plan["timeout_seconds_per_command"])
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
                    receipt.update(exit_code=124, timeout=True)
        except (Exception, KeyboardInterrupt) as error:
            if process is not None and process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            receipt.update(error_type=type(error).__name__, error=str(error))
        receipt["elapsed_seconds"] = time.monotonic() - started
        receipt["finished_utc"] = now()
        try:
            verify(plan)
            if fingerprint(freeze_path) != freeze_identity:
                raise ValueError("Freeze changed during execution")
            receipt["source_input_bookends_match"] = True
        except Exception as error:
            receipt.update(source_input_bookends_match=False, verification_error=str(error))
        try:
            receipt["input_after"] = fingerprint(plan["fastq"])
        except Exception as error:
            receipt["input_after_error"] = str(error)
        receipt["outputs"] = {}
        for path in sorted(run_dir.glob(identifier + ".*")):
            try:
                receipt["outputs"][path.name] = fingerprint(path)
            except Exception as error:
                receipt.setdefault("output_errors", {})[path.name] = str(error)
        receipt["status"] = "complete" if (receipt["exit_code"] == 0 and
            receipt["source_input_bookends_match"] and not any(key in receipt for key in
            ("error", "input_after_error", "output_errors"))) else "failed"
        try:
            write_json(run_dir / f"{identifier}.receipt.json", receipt)
        except Exception as error:
            receipt.update(status="failed", receipt_publication_error=str(error))
        outcomes.append(receipt)
        print(json.dumps({"id": identifier, "status": receipt["status"], "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)
        failed |= receipt["status"] != "complete"
        if not receipt["source_input_bookends_match"] or not receipt["launched"] or receipt.get("error_type") == "KeyboardInterrupt":
            break
    write_json(run_dir / "execution.json", {"status": "failed" if failed else "complete", "freeze": freeze_identity,
               "finished_utc": now(), "commands_planned": [item["id"] for item in plan["commands"]],
               "commands": outcomes})
    return 1 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--fastq", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute:
        return execute(args.freeze.resolve())
    if not args.fastq or not args.run_dir:
        parser.error("--fastq and --run-dir are required when creating a freeze")
    freeze(args.fastq.resolve(), args.run_dir.resolve(), args.freeze.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Package completed inclusive-output evidence, then reopen and hash every member.

Only repository source/evidence payloads are included. External benchmark FASTQ
inputs, compiled binaries, Cargo targets and the toolchain remain external; their
identities remain in the original unmodified verification/check receipts.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import tarfile


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SCRATCH = ROOT / "tmp/inclusive-output-20260929"
SKIP_DIRS = {"target", "toolchain", "cargo", "rustup", "__pycache__", ".git"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).absolute()
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode), f"Not a regular file: {path}")
    identity = lambda row: (row.st_dev, row.st_ino, row.st_size, row.st_mtime_ns, row.st_ctime_ns)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        require(identity(os.fstat(stream.fileno())) == identity(before), f"Replaced before hashing: {path}")
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        require(identity(os.fstat(stream.fileno())) == identity(before), f"Changed while hashing: {path}")
    require(identity(path.lstat()) == identity(before), f"Replaced after hashing: {path}")
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": before.st_size}


class Selection:
    def __init__(self):
        self.files = {}

    def add(self, path, reason, expected=None):
        path = Path(path).absolute()
        relative = path.relative_to(ROOT)
        row = fingerprint(path)
        require(expected is None or row == expected, f"Bound identity differs: {path}")
        name = "payload/" + relative.as_posix()
        if name in self.files:
            require(self.files[name]["file"] == row, f"Identity changed during selection: {path}")
            self.files[name]["reasons"] = sorted(set(self.files[name]["reasons"]) | {reason})
        else:
            self.files[name] = {"file": row, "reasons": [reason]}
        return row

    def bound(self, row, reason):
        return self.add(row["path"], reason, row)

    def document(self, path, reason, expected=None):
        self.add(path, reason, expected)
        return json.loads(Path(path).read_text())

    def bookend(self):
        for value in self.files.values():
            require(fingerprint(value["file"]["path"]) == value["file"], "Selected evidence changed")


def collect_tree(selection, directory, excluded_paths):
    excluded = []
    for parent, dirs, files in os.walk(directory, followlinks=False):
        parent = Path(parent)
        for name in list(dirs):
            path = parent / name
            if name in SKIP_DIRS or path.is_symlink():
                dirs.remove(name)
                excluded.append({"path": str(path), "reason": "generated directory or symlink"})
        for name in sorted(files):
            path = (parent / name).absolute()
            if path in excluded_paths:
                continue
            if (name == "sharkmer" or name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".pyc", ".tar", ".tar.gz", ".tgz", ".zip"))):
                excluded.append({"path": str(path), "reason": "raw input, compiled binary, cache or archive"})
                continue
            require(not path.is_symlink(), f"Unexpected evidence symlink: {path}")
            with path.open("rb") as stream:
                if stream.read(4) == b"\x7fELF":
                    excluded.append({"path": str(path), "reason": "compiled ELF binary"})
                    continue
            selection.add(path, "new bundle or local attempt artifact")
    return excluded


def collect(verification_path, review_path, docs, excluded_paths):
    selection = Selection()
    verification = selection.document(verification_path, "complete six-run verification")
    require(verification.get("status") == "pass", "Verification did not pass")
    require(verification.get("source_before") and verification["source_before"] == verification.get("source_after"), "Source bookends are absent or differ")
    checks = selection.document(verification["checks"]["path"], "complete compiled checks", verification["checks"])
    require(checks.get("status") == "pass" and checks.get("sources_before") == checks.get("sources_after") == verification["source_before"], "Checked source bindings differ")
    expected_steps = {"format", "clippy-default", "release-tests-default", "clippy-fxhashmap", "release-tests-fxhashmap", "python-tests", "reference-audit"}
    require(len(checks.get("steps", [])) == 7 and {step["name"] for step in checks["steps"]} == expected_steps
            and all(step["status"] == "pass" for step in checks["steps"]), "All seven checks must pass")
    for row in verification["source_before"].values():
        selection.bound(row, "exact checked source payload")
    require(verification.get("extra_sources"), "Extra source bindings are absent")
    for row in verification["extra_sources"].values():
        selection.bound(row, "extra source/input/prior-outcome binding")
    for step in checks["steps"]:
        selection.bound(step["log"], "compiled-check log")
    for name in ("comparison_controls", "checks_step", "build_step"):
        require(verification[name]["status"] == "pass", f"{name} did not pass")
        selection.bound(verification[name]["log"], name)
    expected_runs = {(a, mode) for a in ("SRR27962769", "SRR31887760", "SRR1057608") for mode in ("default", "standard_only")}
    jobs = verification.get("jobs", [])
    require(len(jobs) == 6 and {(row["accession"], row["mode"]) for row in jobs} == expected_runs, "Expected six distinct runs")
    for job in jobs:
        require(job.get("status") == "pass" and job.get("prefix_before") == job.get("prefix_after"), "Run failed or prefix bookends differ")
        saved = selection.document(job["receipt"]["path"], "individual run receipt", job["receipt"])
        require(saved == {key: value for key, value in job.items() if key != "receipt"}, "Run receipt differs from verification")
        for row in job["outputs"]["files"].values():
            selection.bound(row, "six-run output")
        previous = selection.document(job["prior_receipt"]["path"], "prior exact comparison receipt", job["prior_receipt"])
        for row in previous["outputs"]["files"].values():
            selection.bound(row, "prior exact comparison output")
    review = selection.document(review_path, "final independent outcome review")
    require(review.get("status") == "pass", "Independent outcome review did not pass")
    require(review.get("verification") == fingerprint(verification_path), "Outcome review does not bind this verification receipt")
    reviewed = review.get("reviewed_artifacts", {})
    rows = list(reviewed.values()) if isinstance(reviewed, dict) else reviewed
    require(isinstance(rows, list) and all(isinstance(row, dict) and "path" in row for row in rows),
            "Outcome review lacks documented artifact fingerprints")
    reviewed_by_path = {row["path"]: row for row in rows}
    require(len(reviewed_by_path) == len(rows), "Duplicate reviewed artifact paths")
    for path in [HERE / "README.md", *docs]:
        path = Path(path).absolute()
        require(str(path) in reviewed_by_path, f"Final documentation was not reviewed: {path}")
        selection.add(path, "reviewed final documentation", reviewed_by_path[str(path)])
    excluded = collect_tree(selection, HERE, excluded_paths)
    excluded += collect_tree(selection, SCRATCH, excluded_paths)
    selection.bookend()
    return selection, excluded


def verify_archive(path, expected_manifest):
    """Reopen the finished archive and independently hash actual member streams."""
    expected = expected_manifest["members"]
    seen = set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            require(member.isfile() and member.name not in seen, "Non-file or duplicate archive member")
            seen.add(member.name)
            stream = archive.extractfile(member)
            if member.name == "MANIFEST.json":
                require(json.load(stream) == expected_manifest, "Embedded inventory differs")
                continue
            require(member.name in expected, f"Unexpected archive member: {member.name}")
            digest, size = hashlib.sha256(), 0
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
                size += len(block)
            row = expected[member.name]["file"]
            require(size == member.size == row["size_bytes"] and digest.hexdigest() == row["sha256"], f"Archive hash differs: {member.name}")
    require(seen == set(expected) | {"MANIFEST.json"}, "Archive member set differs")
    return fingerprint(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--outcome-review", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--doc", type=Path, action="append", default=[])
    args = parser.parse_args()
    archive_path, receipt_path = args.archive.absolute(), args.receipt.absolute()
    require(archive_path != receipt_path and not archive_path.exists() and not receipt_path.exists(), "Archive and receipt must be fresh paths")
    selection, excluded = collect(args.verification.absolute(), args.outcome_review.absolute(), args.doc,
                                  {archive_path, receipt_path})
    manifest = {"schema_version": 1, "purpose": "inclusive-output completed source and result evidence",
                "members": dict(sorted(selection.files.items())), "excluded": excluded,
                "verification": fingerprint(args.verification), "outcome_review": fingerprint(args.outcome_review),
                "scope": "Exact checked source payloads, extra bindings, six current output sets, directly compared prior outputs, all local attempts and final new-bundle documentation/reviews. External read inputs, compiler/runtime and built executables remain represented by their recorded identities only."}
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with archive_path.open("xb") as destination:
        with tarfile.open(fileobj=destination, mode="w:gz", format=tarfile.PAX_FORMAT) as archive:
            info = tarfile.TarInfo("MANIFEST.json")
            info.size = len(manifest_bytes)
            archive.addfile(info, io.BytesIO(manifest_bytes))
            for name, row in manifest["members"].items():
                info = tarfile.TarInfo(name)
                info.size = row["file"]["size_bytes"]
                info.mode = 0o644
                with Path(row["file"]["path"]).open("rb") as stream:
                    archive.addfile(info, stream)
        destination.flush()
        os.fsync(destination.fileno())
    archive_identity = verify_archive(archive_path, manifest)
    selection.bookend()
    receipt = {"schema_version": 1, "status": "pass", "archive": archive_identity,
               "independent_reopen_hash_verification": "pass", "source_bookends_equal": True,
               "member_count": len(selection.files) + 1, "payload_bytes": sum(row["file"]["size_bytes"] for row in selection.files.values()),
               "embedded_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
               "verification": manifest["verification"], "outcome_review": manifest["outcome_review"]}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    with receipt_path.open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()

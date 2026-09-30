#!/usr/bin/env python3
"""Package completed release evidence; reopen and independently hash every member.

Expanded exact source snapshots are included. Their redundant source archives,
external benchmark FASTQs/partials, built executables, targets and toolchains are
excluded and retain their original receipt identities. No packaging occurs on import.
"""

import argparse
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tarfile


ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "tmp/inclusive-release-20260929"
PRIOR_HELPER = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-output-20260929/package_evidence.py"
SPEC = importlib.util.spec_from_file_location("reviewed_packaging_primitives", PRIOR_HELPER)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
fingerprint, require = BASE.fingerprint, BASE.require
SKIP_DIRS = {"target", "toolchain", "cargo", "rustup", "__pycache__", ".git"}
MAX_FILE_BYTES = 256 * 1024 * 1024
MAX_PAYLOAD_BYTES = 512 * 1024 * 1024


class Selection(BASE.Selection):
    def __init__(self):
        super().__init__()
        self.by_path = {}
        self.excluded = {}
        self.source_paths = set()

    def add(self, path, reason, expected=None):
        path = Path(path).absolute()
        require(path.stat().st_size <= MAX_FILE_BYTES, f"Unexpectedly large evidence file: {path}")
        if str(path) in self.by_path:
            row = self.by_path[str(path)]
            require(expected is None or row == expected, f"Conflicting bound identity: {path}")
            name = "payload/" + path.relative_to(ROOT).as_posix()
            self.files[name]["reasons"] = sorted(set(self.files[name]["reasons"]) | {reason})
            return row
        row = super().add(path, reason, expected)
        self.by_path[str(path)] = row
        require(sum(value["file"]["size_bytes"] for value in self.files.values()) <= MAX_PAYLOAD_BYTES,
                "Evidence payload exceeds packaging bound")
        return row

    def omit(self, path, reason, expected=None):
        if expected is not None or str(path) not in self.excluded:
            self.excluded[str(path)] = {"path": str(path), "reason": reason,
                                        **({"recorded_identity": expected} if expected else {})}


def exclusion(path, source_fixture=False):
    path = Path(path).absolute()
    if not path.is_relative_to(ROOT):
        return "external tool/input: identity retained in original receipt"
    if any(part in SKIP_DIRS for part in path.relative_to(ROOT).parts):
        return "build target, toolchain or cache"
    name = path.name.lower()
    if name.endswith((".tar", ".tar.gz", ".tgz", ".zip")):
        return "nested archive omitted; source snapshots are included as expanded payloads"
    if not source_fixture and (".fastq" in name or name.endswith(".fq") or ".fq." in name):
        return "external benchmark FASTQ or failed transfer partial"
    if name == "sharkmer" or name.endswith(".pyc"):
        return "compiled executable or Python bytecode"
    require(not path.is_symlink(), f"Evidence symlink refused: {path}")
    with path.open("rb") as stream:
        if stream.read(4) == b"\x7fELF":
            return "compiled ELF executable"
    return None


def bound(selection, row, reason):
    path = Path(row["path"]).absolute()
    omitted = exclusion(path, path in selection.source_paths)
    if omitted:
        selection.omit(path, omitted, row)
    else:
        selection.bound(row, reason)


def references(selection, value, reason):
    """Collect fingerprint rows within an already selected receipt, without following JSON recursively."""
    if isinstance(value, dict):
        if {"path", "sha256", "size_bytes"} <= set(value):
            row = {key: value[key] for key in ("path", "sha256", "size_bytes")}
            bound(selection, row, reason)
        else:
            for child in value.values():
                references(selection, child, reason)
    elif isinstance(value, list):
        for child in value:
            references(selection, child, reason)


def source_tree(selection, root, files, reason):
    root = Path(root).absolute()
    require(files, "Empty checked source snapshot")
    for relative, original in files.items():
        path = (root / relative).absolute()
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts, "Unsafe source member path")
        expected = {"path": str(path), "size_bytes": original["size_bytes"], "sha256": original["sha256"]}
        require("path" not in original or original["path"] == str(path), "Source map path disagrees with root")
        selection.source_paths.add(path)
        require(exclusion(path, True) is None, f"Unexpected executable/archive in checked source map: {path}")
        selection.bound(expected, reason)


def collect_tree(selection, directory, excluded_paths):
    require(directory.is_dir() and not directory.is_symlink(), f"Missing/invalid evidence directory: {directory}")
    for parent, dirs, files in os.walk(directory, followlinks=False):
        parent = Path(parent)
        for name in list(dirs):
            path = parent / name
            if name in SKIP_DIRS or path.is_symlink():
                dirs.remove(name)
                selection.omit(path, "generated directory or symlink")
        for name in sorted(files):
            path = (parent / name).absolute()
            if path in excluded_paths:
                continue
            omitted = exclusion(path, path in selection.source_paths)
            if omitted:
                selection.omit(path, omitted)
            else:
                selection.add(path, "new bundle or local attempt artifact, including failures")


def reviewed_paths(review):
    require(review.get("status") == "pass", "Final independent outcome review did not pass")
    values = review.get("reviewed_artifacts", {})
    rows = list(values.values()) if isinstance(values, dict) else values
    require(isinstance(rows, list) and rows, "Missing reviewed artifact map")
    result = {row["path"]: row for row in rows}
    require(len(result) == len(rows), "Duplicate reviewed artifact paths")
    return result


def collect(args, excluded_paths):
    selection = Selection()
    documents = {}
    for name in ("execution", "report", "builds", "availability", "dros_receipt", "dros_controls", "outcome_review"):
        path = getattr(args, name).absolute()
        documents[name] = selection.document(path, name.replace("_", " "))
    execution, report, builds = (documents[name] for name in ("execution", "report", "builds"))
    require(execution.get("status") == report.get("status") == builds.get("status") == "pass", "Main execution/report/build gate failed")
    for path in (args.builds, args.availability):
        require(execution["bindings"].get(str(path.absolute())) == fingerprint(path),
                "Supplied builds/availability are not bound by main execution")
    jobs = execution.get("jobs", [])
    require(len(jobs) == 60 and all(job.get("status") == "pass" for job in jobs), "Need sixty successful main invocations")
    require(len({job["invocation"]["invocation_id"] for job in jobs}) == 60, "Duplicate main invocation")
    summary = selection.document(execution["summary"]["path"], "complete comparison summary", execution["summary"])
    require(summary.get("status") == "pass", "Summary failed")
    require(report["sources_before"] == report["sources_after"] and
            report["sources_before"]["execution"] == fingerprint(args.execution) and
            report["sources_before"]["summary"] == execution["summary"], "Posthoc report does not bind execution/summary")
    mandatory_review = [args.execution, Path(execution["summary"]["path"]), args.report, args.builds,
                        args.availability, args.dros_receipt, args.dros_controls, HERE / "README.md", *args.doc]

    require(set(builds["versions"]) == {"baseline", "candidate"}, "Wrong production build roles")
    for role, version in builds["versions"].items():
        snapshot = selection.document(version["snapshot"]["path"], f"{role} source snapshot", version["snapshot"])
        build = selection.document(version["build"]["path"], f"{role} build receipt", version["build"])
        require(build.get("status") == "complete" and build["returncode"] == 0, "Production build failed")
        require(snapshot["files"] == build["source_before"] == build["source_after"], "Production source bookends differ")
        source_tree(selection, snapshot["source"], snapshot["files"], f"exact {role} checked source payload")
        references(selection, build, f"{role} build lineage")
        references(selection, snapshot, f"{role} source export lineage")
        mandatory_review.extend(Path(version[key]["path"]) for key in ("snapshot", "build"))

    dros, controls = documents["dros_receipt"], documents["dros_controls"]
    require(dros.get("status") in ("pass", "complete") and controls.get("status") == "pass", "Focused replay/checks failed")
    require(controls.get("steps") and all(step.get("status") == "pass" for step in controls["steps"]),
            "Focused compiled check step failed or is missing")
    require(dros["before"] == dros["after"] == dros["after_audit"], "Focused replay source bookends differ")
    require(controls["source_before"] == controls["source_after"], "Focused checked sources differ")
    frozen = dros["source_freeze"]
    source_tree(selection, frozen["source"], frozen["files"], "exact focused diagnostic source payload")
    require({key: {field: row[field] for field in ("size_bytes", "sha256")} for key, row in controls["source_before"].items()}
            == frozen["files"], "Focused replay differs from checked source map")
    references(selection, dros, "focused replay source/control/input/output lineage")
    references(selection, controls, "focused compiled-check lineage")
    assessment = selection.document(dros["assessment"]["path"], "focused audit assessment", dros["assessment"])
    mandatory_review.append(Path(dros["assessment"]["path"]))

    available = documents["availability"]
    require(available.get("status") == "pass" and available.get("complete_inventory_selected") is True,
            "Incomplete merged input availability")
    require(available.get("sources_before") and available["sources_before"] == available.get("sources_after"),
            "Restoration merge source bookends are absent or differ")
    require(len(available["samples"]) == len({row["accession"] for row in available["samples"]}) == 10,
            "Wrong available sample set")
    references(selection, available, "restoration merge and attempt lineage")
    for attempt in available["source_attempts"]:
        saved = selection.document(attempt["receipt"]["path"], "source restoration attempt", attempt["receipt"])
        references(selection, saved, "successful and failed restoration receipts")

    for job in jobs:
        saved = selection.document(job["receipt"]["path"], "individual main invocation", job["receipt"])
        require(saved == {key: value for key, value in job.items() if key != "receipt"}, "Embedded main invocation differs from saved receipt")
        references(selection, job["outputs"]["files"], "main run output")
        references(selection, job["execution"]["files"], "main timing/constraint/log evidence")
    references(selection, execution["bindings"], "frozen main source/input/build bindings")
    references(selection, report, "posthoc tables and checked report sources")
    for path in args.evidence:
        extra = selection.document(path, "additional checked evidence")
        references(selection, extra, "additional checked evidence lineage")

    review_map = reviewed_paths(documents["outcome_review"])
    for path in mandatory_review:
        path = Path(path).absolute()
        require(str(path) in review_map, f"Outcome review does not bind required artifact: {path}")
        selection.bound(review_map[str(path)], "final reviewed outcome/documentation")
    for row in review_map.values():
        bound(selection, row, "independently reviewed artifact")
    selection.add(PRIOR_HELPER, "reused reviewed packaging primitives")
    selection.add(__file__, "packaging source")
    collect_tree(selection, HERE, excluded_paths)
    collect_tree(selection, SCRATCH, excluded_paths)
    selection.bookend()
    return selection


def sha_listing(manifest):
    return "".join(f"{row['file']['sha256']}  {name}\n" for name, row in sorted(manifest["members"].items())).encode()


def verify_archive(path, manifest):
    expected = manifest["members"]
    seen = set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            require(member.isfile() and member.name not in seen, "Duplicate or non-file archive member")
            seen.add(member.name)
            stream = archive.extractfile(member)
            if member.name == "MANIFEST.json":
                require(json.load(stream) == manifest, "Embedded manifest differs")
            elif member.name == "SHA256SUMS":
                require(stream.read() == sha_listing(manifest), "Embedded SHA256 listing differs")
            else:
                require(member.name in expected, "Unexpected archive member")
                digest, size = hashlib.sha256(), 0
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
                    size += len(block)
                row = expected[member.name]["file"]
                require(size == member.size == row["size_bytes"] and digest.hexdigest() == row["sha256"],
                        f"Archived payload hash/size differs: {member.name}")
    require(seen == set(expected) | {"MANIFEST.json", "SHA256SUMS"}, "Archive member set differs")
    return fingerprint(path)


def publish(selection, archive_path, receipt_path, anchors):
    require(archive_path != receipt_path and all(not path.exists() and not path.is_symlink() for path in (archive_path, receipt_path)),
            "Archive and receipt must be fresh distinct paths")
    manifest = {"schema_version": 1, "purpose": "inclusive release source and outcome evidence",
                "members": dict(sorted(selection.files.items())), "excluded": list(selection.excluded.values()), "anchors": anchors,
                "scope": "Exact expanded source snapshots, run outputs/logs, restoration lineage, focused replay, reviews and documentation; external inputs and tools remain represented by receipts."}
    manifest_bytes = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    selection.bookend()
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with archive_path.open("xb") as destination:
        with tarfile.open(fileobj=destination, mode="w:gz", format=tarfile.PAX_FORMAT) as archive:
            for name, raw in (("MANIFEST.json", manifest_bytes), ("SHA256SUMS", sha_listing(manifest))):
                member = tarfile.TarInfo(name)
                member.size, member.mode = len(raw), 0o644
                archive.addfile(member, io.BytesIO(raw))
            for name, row in manifest["members"].items():
                member = tarfile.TarInfo(name)
                member.size, member.mode = row["file"]["size_bytes"], 0o644
                with Path(row["file"]["path"]).open("rb") as stream:
                    archive.addfile(member, stream)
        destination.flush()
        os.fsync(destination.fileno())
    verified = verify_archive(archive_path, manifest)
    selection.bookend()
    receipt = {"schema_version": 1, "status": "pass", "archive": verified, "anchors": anchors,
               "independent_reopen_hash_verification": "pass", "source_bookends_equal": True,
               "member_count": len(selection.files) + 2, "payload_bytes": sum(row["file"]["size_bytes"] for row in selection.files.values()),
               "embedded_manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest()}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    with receipt_path.open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("execution", "report", "builds", "availability", "dros-receipt", "dros-controls", "outcome-review", "archive", "receipt"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--doc", type=Path, action="append", default=[])
    parser.add_argument("--evidence", type=Path, action="append", default=[])
    args = parser.parse_args()
    args.archive, args.receipt = args.archive.absolute(), args.receipt.absolute()
    require(args.archive != args.receipt and all(not path.exists() and not path.is_symlink() for path in (args.archive, args.receipt)),
            "Archive and receipt must be fresh distinct paths")
    selection = collect(args, {args.archive, args.receipt})
    anchors = {name: fingerprint(getattr(args, name)) for name in ("execution", "report", "builds", "availability", "dros_receipt", "dros_controls", "outcome_review")}
    print(json.dumps(publish(selection, args.archive, args.receipt, anchors), sort_keys=True))


if __name__ == "__main__":
    main()

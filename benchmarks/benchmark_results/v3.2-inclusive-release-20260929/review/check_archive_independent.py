#!/usr/bin/env python3
"""Independent archive verification; no imports from either packaging helper."""
import hashlib
import json
import tarfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = Path(__file__).resolve().parents[1]
SCRATCH = ROOT / "tmp/inclusive-release-20260929"


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": size}


def load(path):
    return json.loads(Path(path).read_text())


def main():
    receipt_path = BUNDLE / "archive_receipt.json"
    receipt = load(receipt_path)
    archive_path = BUNDLE / "evidence.tar.gz"
    before = fingerprint(archive_path)
    assert receipt["status"] == "pass" and receipt["archive"] == before
    with tarfile.open(archive_path, "r:gz") as archive:
        manifest_bytes = archive.extractfile("MANIFEST.json").read()
    assert hashlib.sha256(manifest_bytes).hexdigest() == receipt["embedded_manifest_sha256"]
    manifest = json.loads(manifest_bytes)
    assert manifest["anchors"] == receipt["anchors"]
    expected = manifest["members"]
    by_original = {}
    for name, row in expected.items():
        original = row["file"]["path"]
        assert original not in by_original, ("duplicate original path", original)
        by_original[original] = row["file"]
        assert name == "payload/" + str(Path(original).relative_to(ROOT))

    def required(row):
        assert by_original.get(row["path"]) == row, ("missing or changed required payload", row)

    def required_path(path):
        row = fingerprint(path)
        required(row)
        return row

    source_paths = set()
    source_counts = {}
    builds = load(SCRATCH / "build-v1/builds.json")
    for role, version in builds["versions"].items():
        required(version["snapshot"])
        required(version["build"])
        snapshot = load(version["snapshot"]["path"])
        build = load(version["build"]["path"])
        assert snapshot["files"] == build["source_before"] == build["source_after"]
        for relative, row in snapshot["files"].items():
            path = str((Path(snapshot["source"]) / relative).resolve())
            required({"path": path, "sha256": row["sha256"], "size_bytes": row["size_bytes"]})
            source_paths.add(path)
        source_counts[role] = len(snapshot["files"])
    dros = load(SCRATCH / "dros12s/replay-v1/replay.receipt.json")
    frozen = dros["source_freeze"]
    for relative, row in frozen["files"].items():
        path = str((Path(frozen["source"]) / relative).resolve())
        required({"path": path, "sha256": row["sha256"], "size_bytes": row["size_bytes"]})
        source_paths.add(path)
    source_counts["dros_scratch"] = len(frozen["files"])
    checks_path = ROOT / "tmp/inclusive-output-20260929/attempt-v1/checks/checks.json"
    checks = load(checks_path)
    required_path(checks_path)
    assert checks["sources_before"] == checks["sources_after"]
    prior_source_aliases = {}
    for row in checks["sources_before"].values():
        if row["path"] not in by_original:
            relative = Path(row["path"]).relative_to(ROOT)
            alias = str((SCRATCH / "build-v1/candidate/source" / relative).resolve())
            required({"path": alias, "sha256": row["sha256"], "size_bytes": row["size_bytes"]})
            prior_source_aliases[row["path"]] = by_original[alias]
        else:
            required(row)
        source_paths.add(row["path"])
    source_counts["prior_full_checks"] = len(checks["sources_before"])
    fixture_paths = {path for path in source_paths if path in by_original and any(path.endswith(ext) for ext in (".fastq", ".fastq.gz", ".fq", ".fq.gz"))}

    checksums = "".join(f"{expected[name]['file']['sha256']}  {name}\n" for name in sorted(expected)).encode()
    seen = set()
    total = 0
    first_bytes = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            parts = PurePosixPath(member.name).parts
            assert member.isfile() and member.name not in seen
            assert not PurePosixPath(member.name).is_absolute() and ".." not in parts
            seen.add(member.name)
            stream = archive.extractfile(member)
            if member.name == "MANIFEST.json":
                assert stream.read() == manifest_bytes
                continue
            if member.name == "SHA256SUMS":
                assert stream.read() == checksums
                continue
            assert member.name in expected
            row = expected[member.name]["file"]
            digest = hashlib.sha256()
            size = 0
            first = b""
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                if not first:
                    first = block[:4]
                digest.update(block)
                size += len(block)
            assert size == member.size == row["size_bytes"]
            assert digest.hexdigest() == row["sha256"]
            assert first != b"\x7fELF", member.name
            original = Path(row["path"])
            assert not {"target", "toolchain", "cargo-home", "rustup-home", "__pycache__", ".git"}.intersection(original.parts)
            assert not str(original).endswith((".tar", ".tar.gz", ".tgz", ".zip"))
            if str(original).endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
                assert str(original) in fixture_paths, ("external raw input payload", original)
            assert original.name != "sharkmer"
            total += size
    assert seen == set(expected) | {"MANIFEST.json", "SHA256SUMS"}
    assert len(seen) == receipt["member_count"] and total == receipt["payload_bytes"]
    for row in receipt["anchors"].values():
        required(row)
        assert fingerprint(row["path"]) == row
    outcome = load(BUNDLE / "review/outcome_review.json")
    assert outcome["status"] == "pass"
    for path, row in outcome["reviewed_artifacts"].items():
        if path in by_original:
            assert by_original[path] == row
    for doc in ["dev_docs/INCLUSIVE_RELEASE_PROTOCOL.md", "dev_docs/INCLUSIVE_RELEASE_RESULTS.md",
                "dev_docs/DROS12S_RETENTION_RESULTS.md", "dev_docs/HANDOFF.md", "dev_docs/PLAN.md",
                "dev_docs/REVIEW_v3.2.md", "CHANGELOG.md", "AGENTS.md"]:
        required_path(ROOT / doc)
    required_path(BUNDLE / "README.md")
    execution = load(SCRATCH / "comparison-v1/execution.json")
    assert len(execution["jobs"]) == 60
    main_outputs = set()
    for job in execution["jobs"]:
        assert job["status"] == "pass"
        required(job["receipt"])
        for collection in [job["outputs"]["files"], job["execution"]["files"]]:
            for row in collection.values():
                required(row)
                main_outputs.add(row["path"])
    report = load(BUNDLE / "results/report-v2/report.json")
    for row in report["files"].values():
        required(row)
    required(dros["assessment"])
    assert load(dros["assessment"]["path"])["status"] == "pass"
    availability = load(SCRATCH / "availability.json")
    for attempt in availability["source_attempts"]:
        required(attempt["receipt"])
    for path in [BUNDLE / "results/report-v1/report.json", BUNDLE / "results/report_comparison.py",
                 BUNDLE / "inputs/restore_inputs.py", BUNDLE / "inputs/restore_inputs_v2.py"]:
        required_path(path)
    sidecars = [p for p in by_original if p.startswith(str(SCRATCH / "dros12s/replay-v1")) and Path(p).name.startswith("completion_ledger_")]
    assert len(sidecars) == 16
    assert fingerprint(archive_path) == before
    result = {"schema_version": 1, "status": "pass", "date": "2026-09-30",
              "purpose": "Independent archive member/hash, required output/source/document coverage and exclusion review",
              "archive": before, "packaging_receipt": fingerprint(receipt_path),
              "checker": fingerprint(__file__), "outcome_review": fingerprint(BUNDLE / "review/outcome_review.json"),
              "member_count": len(seen), "payload_count": len(expected), "payload_bytes": total,
              "expanded_source_counts": source_counts, "retained_source_fastq_fixtures": sorted(fixture_paths),
              "prior_checked_source_equivalent_payloads": prior_source_aliases,
              "main_invocations": 60, "main_output_and_execution_files": len(main_outputs), "complete_dros_ledgers": len(sidecars),
              "checks": ["Every payload and exact member set match the embedded manifest and SHA256SUMS; no links, duplicate or unsafe member names.",
                         "Both exact production source trees, complete scratch diagnostic source and all prior full-check sources are included.",
                         "All sixty invocation receipts, output/timing/log/constraint files, six report tables, reviewed documents and anchor fingerprints are included.",
                         "Failed restoration/report attempts and versioned corrected sources are preserved.",
                         "External raw FASTQs, executable payloads, target/toolchain/cache trees and nested archives are excluded; exact checked-source FASTQ fixtures are retained.",
                         "Archive digest is unchanged across the independent reopen and coverage checks."],
              "limitations": ["Evidence integrity pass does not resolve the four runtime flags, validate uncertain biological hypotheses, or authorize release."]}
    output = BUNDLE / "review/archive_review.json"
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": "pass", "review": fingerprint(output), "archive": before}))


if __name__ == "__main__":
    main()

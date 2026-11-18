#!/usr/bin/env python3
"""Retain every attempted paired comparison and its exact source bindings."""

import hashlib
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

ROOT = Path("tmp/v4-direct-ingestion-20261110")
DESTINATION = Path("benchmarks/benchmark_results/v4-direct-comparison-20261112")


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def main() -> None:
    audit = json.loads((ROOT / "verification-v2.json").read_text())
    assert audit["verification"] == "passed" and not audit["promotion_gates_passed"]
    DESTINATION.mkdir(parents=True, exist_ok=False)
    archives = {}
    for version, commit in [("v1", "27ea351"), ("v2", "724db89")]:
        source = ROOT / f"paired-{version}"
        bindings = ROOT / f"candidate-{version}-source.json"
        for path, expected in json.loads(bindings.read_text()).items():
            blob = subprocess.check_output(["git", "show", f"{commit}:{path}"])
            assert hashlib.sha256(blob).hexdigest() == expected, (version, path)
        shutil.copy2(bindings, DESTINATION / bindings.name)
        members = {
            str(path.relative_to(source)): {
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            }
            for path in sorted(source.rglob("*"))
            if path.is_file() and not path.is_symlink()
        }
        archive = DESTINATION / f"raw-results-{version}.tar.gz"
        with tarfile.open(archive, "w:gz", compresslevel=6) as tar:
            for relative in members:
                tar.add(source / relative, arcname=relative, recursive=False)
        # Re-open the retained artifact rather than trusting archive creation.
        with tarfile.open(archive, "r:gz") as tar:
            assert set(tar.getnames()) == set(members)
            for member in tar:
                stream = tar.extractfile(member)
                assert stream is not None
                data = stream.read()
                assert len(data) == members[member.name]["bytes"]
                assert (
                    hashlib.sha256(data).hexdigest() == members[member.name]["sha256"]
                )
        (DESTINATION / f"artifact-manifest-{version}.json").write_text(
            json.dumps(members, indent=2) + "\n"
        )
        (DESTINATION / f"harness-{version}.py").write_bytes(
            subprocess.check_output(
                ["git", "show", f"{commit}:benchmarks/compare_count_backends.py"]
            )
        )
        archives[version] = {
            "source_commit": subprocess.check_output(
                ["git", "rev-parse", commit], text=True
            ).strip(),
            "members": len(members),
            "raw_bytes": sum(row["bytes"] for row in members.values()),
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": digest(archive),
        }
    for name in [
        "candidate-v1-rejection.json",
        "verification-v2.json",
        "verify_comparison.py",
        "verify.sh",
        "batch_compare.sh",
        "batch_compare_v2.sh",
    ]:
        shutil.copy2(ROOT / name, DESTINATION / name)
    for name in ["paired-32947879.out", "paired-33016469.out", "verify-33120716.out"]:
        shutil.copy2(ROOT / "logs" / name, DESTINATION / name)
    (DESTINATION / "archives.json").write_text(json.dumps(archives, indent=2) + "\n")
    print(json.dumps(archives, indent=2))


if __name__ == "__main__":
    main()

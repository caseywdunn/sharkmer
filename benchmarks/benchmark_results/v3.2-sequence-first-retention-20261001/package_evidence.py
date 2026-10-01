#!/usr/bin/env python3
"""Archive the completed assessment without copying large FASTQ inputs."""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def add_tree(files, label, base):
    for path in base.rglob("*"):
        if path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc":
            name = f"{label}/{path.relative_to(base).as_posix()}"
            if name in files:
                raise ValueError(f"Duplicate archive name: {name}")
            files[name] = path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    scratch = args.scratch.resolve()
    full = json.loads((scratch / "full-v3/result.json").read_text())
    if full["status"] != "pass" or len(full["jobs"]) != 60:
        raise ValueError("Full comparison is incomplete")
    for name in ("dros-v3", "plant-v3", "optout-v1"):
        if json.loads((scratch / name / "result.json").read_text())["status"] != "pass":
            raise ValueError(f"{name} is incomplete")
    checks = scratch / "checks-v3"
    if (checks / "sources.before.sha256").read_bytes() != (checks / "sources.after.sha256").read_bytes():
        raise ValueError("Compiled-check sources changed")

    files = {}
    for name in ("checks-v1", "checks-v2", "checks-v3", "dros-v2", "dros-v3", "plant-v3",
                 "optout-v1", "full-v1", "full-v3", "slurm"):
        add_tree(files, f"scratch/{name}", scratch / name)
    for path in HERE.rglob("*"):
        if (path.is_file() and "__pycache__" not in path.parts and path.suffix != ".pyc"
                and path.name not in ("evidence.tar.gz", "archive_receipt.json", "SHA256SUMS")):
            files[f"bundle/{path.relative_to(HERE).as_posix()}"] = path
    for name in ("README.md", "AGENTS.md", "dev_docs/SEQUENCE_RETENTION_PROTOCOL.md",
                 "dev_docs/SEQUENCE_RETENTION_RESULTS.md", "dev_docs/HANDOFF.md",
                 "dev_docs/PLAN.md", "src/pcr/mod.rs", "src/pcr/paths.rs",
                 "src/pcr/uncertain.rs", "scripts/sharkmer_validate/uncertain.py",
                 "tests/test_uncertain_output.py"):
        files[f"source/{name}"] = ROOT / name
    files["inputs/availability.json"] = ROOT / "tmp/inclusive-release-20260929/availability.json"
    files["inputs/inventory.json"] = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/inventory.json"
    files["inputs/released_identity_dispositions.tsv"] = (ROOT
        / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929/results/report-v2/released_identity_dispositions.tsv")

    manifest = {name: {"sha256": digest(path), "size_bytes": path.stat().st_size}
                for name, path in sorted(files.items())}
    manifest_bytes = (json.dumps({"schema_version": 1, "files": manifest}, indent=2, sort_keys=True)
                      + "\n").encode()
    archive = HERE / "evidence.tar.gz"
    with archive.open("xb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as compressed:
            with tarfile.open(fileobj=compressed, mode="w|") as tar:
                for name, path in sorted(files.items()):
                    info = tarfile.TarInfo(name)
                    info.size = path.stat().st_size
                    info.mtime = 0
                    info.mode = 0o644
                    with path.open("rb") as stream:
                        tar.addfile(info, stream)
                info = tarfile.TarInfo("manifest.json")
                info.size = len(manifest_bytes)
                info.mtime = 0
                info.mode = 0o644
                tar.addfile(info, io.BytesIO(manifest_bytes))

    observed = set()
    with tarfile.open(archive, mode="r:gz") as tar:
        for member in tar:
            if member.name == "manifest.json":
                if tar.extractfile(member).read() != manifest_bytes:
                    raise ValueError("Archive manifest changed")
                continue
            if member.name not in manifest or member.name in observed:
                raise ValueError("Unexpected/duplicate archive member")
            observed.add(member.name)
            hasher = hashlib.sha256()
            total = 0
            stream = tar.extractfile(member)
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                hasher.update(chunk)
                total += len(chunk)
            if hasher.hexdigest() != manifest[member.name]["sha256"] or total != manifest[member.name]["size_bytes"]:
                raise ValueError(f"Corrupt archive member: {member.name}")
    if observed != set(manifest):
        raise ValueError("Missing archive member")
    receipt = {"status": "pass", "schema_version": 1, "archive": archive.name,
               "sha256": digest(archive), "size_bytes": archive.stat().st_size,
               "member_count": len(manifest), "member_bytes": sum(row["size_bytes"] for row in manifest.values())}
    (HERE / "archive_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    (HERE / "SHA256SUMS").write_text(f"{receipt['sha256']}  {archive.name}\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()

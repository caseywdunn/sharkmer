#!/usr/bin/env python3
"""Package the completed experiment without raw FASTQ or overwritten evidence."""

import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile


BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
RUN = ROOT / "tmp/read-histories-20260926/run-v1"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    archive = BASE / "evidence.tar.gz"
    index = BASE / "ARCHIVE_CONTENTS.json"
    sums = BASE / "SHA256SUMS"
    if any(path.exists() for path in (archive, index, sums)):
        raise ValueError("Archive outputs exist; do not overwrite")
    execution = json.loads((RUN / "execution.json").read_text())
    if execution["status"] != "complete":
        raise ValueError("Execution is not complete")
    review = json.loads((BASE / "results_review/provenance_ledger_review.retry1.json").read_text())
    if review["status"] != "pass":
        raise ValueError("Independent provenance/ledger review did not pass")
    final_review = json.loads((BASE / "results_review/final_results_review.json").read_text())
    if final_review["status"] != "pass":
        raise ValueError("Final interpretation review did not pass")
    for relative, expected in final_review["reviewed_artifacts"].items():
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or (ROOT / path).is_symlink():
            raise ValueError("Invalid reviewed artifact path")
        raw = (ROOT / path).read_bytes()
        if {"size_bytes": len(raw), "sha256": sha(raw)} != expected:
            raise ValueError(f"Artifact changed after final review: {relative}")
    files = {}
    for prefix, directory in (("experiment", BASE), ("run-v1", RUN)):
        for path in sorted(directory.rglob("*")):
            if "__pycache__" in path.parts or not path.is_file():
                continue
            if path.is_symlink() or path.suffix in (".fastq", ".fq", ".pyc"):
                raise ValueError(f"Unexpected evidence member: {path}")
            name = f"{prefix}/{path.relative_to(directory).as_posix()}"
            files[name] = path.read_bytes()
    for relative in ("dev_docs/READ_HISTORY_EXPERIMENT.md", "dev_docs/READ_HISTORY_RESULTS.md",
                     "dev_docs/HANDOFF.md", "dev_docs/PLAN.md", "dev_docs/REVIEW_v3.2.md",
                     "scripts/audit_read_histories.py", "tests/test_read_histories.py"):
        files["repository_snapshots/" + relative] = (ROOT / relative).read_bytes()
    members = [{"path": name, "size_bytes": len(data), "sha256": sha(data)} for name, data in sorted(files.items())]
    with archive.open("xb") as output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as bundle:
                for item in members:
                    info = tarfile.TarInfo(item["path"])
                    info.size = item["size_bytes"]
                    info.mode = 0o644
                    bundle.addfile(info, io.BytesIO(files[item["path"]]))
    with tarfile.open(archive, "r:gz") as bundle:
        actual = []
        for item in bundle.getmembers():
            data = bundle.extractfile(item).read()
            actual.append({"path": item.name, "size_bytes": len(data), "sha256": sha(data)})
        if actual != members:
            raise ValueError("Archive read-back verification failed")
    with index.open("x") as output:
        json.dump({"archive_format": "tar.gz", "members": members}, output, indent=2)
        output.write("\n")
    archive_sha = sha(archive.read_bytes())
    with sums.open("x") as output:
        output.write(f"{archive_sha}  evidence.tar.gz\n{sha(index.read_bytes())}  ARCHIVE_CONTENTS.json\n")
    print(json.dumps({"members": len(members), "archive_bytes": archive.stat().st_size,
                      "archive_sha256": archive_sha, "read_back_verified": True}))


if __name__ == "__main__":
    main()

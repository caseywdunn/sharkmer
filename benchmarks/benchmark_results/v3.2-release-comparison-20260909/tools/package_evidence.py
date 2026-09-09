import gzip
import hashlib
import json
import shutil
import tarfile
from pathlib import Path


workspace = Path("/tmp/sharkmer-release-comparison")
destination = Path("/home/claude/repos/sharkmer/benchmarks/benchmark_results/v3.2-release-comparison-20260909")
destination.mkdir(parents=True, exist_ok=True)
archive_sources = [
    "execution",
    "normalized",
    "ena",
    "supplemental-synthetic-20260909T162310Z",
]
inventory = []
with (destination / "evidence.tar.gz").open("wb") as archive_output:
    with gzip.GzipFile(filename="", mode="wb", fileobj=archive_output, mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w|") as archive:
            for source_name in archive_sources:
                for path in sorted((workspace / source_name).rglob("*")):
                    if path.is_symlink():
                        raise ValueError(f"Unexpected evidence symlink: {path}")
                    if not path.is_file():
                        continue
                    relative = path.relative_to(workspace).as_posix()
                    content = path.read_bytes()
                    inventory.append({"path": relative, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})
                    information = archive.gettarinfo(str(path), arcname=relative)
                    information.uid = information.gid = 0
                    information.uname = information.gname = ""
                    information.mtime = 0
                    information.mode = 0o644
                    with path.open("rb") as input_stream:
                        archive.addfile(information, input_stream)
(destination / "ARCHIVE_CONTENTS.json").write_text(json.dumps(inventory, indent=2) + "\n")
copies = {
    "builds.json": "builds.json",
    "protocol.json": "protocol.json",
    "inputs.json": "inputs.json",
    "normalized/comparison.json": "comparison.json",
    "analysis/analysis.json": "analysis.json",
    "analysis/analysis.md": "tables.md",
    "supplemental-synthetic-20260909T162310Z/summary.json": "synthetic.json",
    "staging.log": "staging.log",
    "driver.log": "driver.log",
    "postprocess.log": "postprocess.log",
}
tool_names = [
    "build_versions.py",
    "stage_inputs.py",
    "driver.py",
    "test_driver.py",
    "postprocess_legacy.py",
    "test_postprocess_legacy.py",
    "analyze_comparison.py",
    "supplemental_synthetic_probes.py",
    "package_evidence.py",
]
copies.update({name: f"tools/{name}" for name in tool_names})
for source_name, target_name in copies.items():
    target = destination / target_name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(workspace / source_name, target)
print(f"Packaged {len(inventory)} unchanged evidence files in {destination}")

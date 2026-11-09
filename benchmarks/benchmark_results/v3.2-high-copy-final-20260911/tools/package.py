import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, action="append", default=[])
    parser.add_argument("--extra-dir", type=Path, action="append", default=[])
    arguments = parser.parse_args()
    destination = arguments.destination
    if destination.exists():
        raise ValueError("Do not overwrite archived evidence")
    destination.mkdir(parents=True)
    inventory = []
    with (destination / "evidence.tar.gz").open("wb") as output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|") as archive:
                for source in (arguments.execution, arguments.analysis, *arguments.extra_dir):
                    for path in sorted(source.rglob("*")):
                        if path.is_symlink():
                            raise ValueError(f"Unexpected evidence symlink: {path}")
                        if not path.is_file():
                            continue
                        relative = f"{source.name}/{path.relative_to(source).as_posix()}"
                        inventory.append({"path": relative, "bytes": path.stat().st_size, "sha256": digest(path)})
                        information = archive.gettarinfo(str(path), arcname=relative)
                        information.uid = information.gid = 0
                        information.uname = information.gname = ""
                        information.mtime = 0
                        information.mode = 0o644
                        with path.open("rb") as input_stream:
                            archive.addfile(information, input_stream)
    (destination / "ARCHIVE_CONTENTS.json").write_text(json.dumps(inventory, indent=2) + "\n")
    for source in (arguments.execution / "comparison.json", arguments.execution / "provenance.json"):
        shutil.copyfile(source, destination / source.name)
    for path in sorted(arguments.analysis.iterdir()):
        if path.is_file():
            shutil.copyfile(path, destination / path.name)
    for build_directory in arguments.build_dir:
        target = destination / "builds" / build_directory.name
        target.mkdir(parents=True)
        for source in sorted(build_directory.iterdir()):
            if source.is_file() and source.suffix in (".json", ".jsonl", ".stderr"):
                shutil.copyfile(source, target / source.name)
    tool_directory = destination / "tools"
    tool_directory.mkdir()
    for name in ("package.py", "analyze.py", "test_benchmark.py", "test_analyze.py"):
        source = Path(__file__).with_name(name)
        shutil.copyfile(source, tool_directory / name)
    if arguments.extra_dir:
        shutil.copyfile(Path(__file__).with_name("controls.py"), tool_directory / "controls.py")
        shutil.copyfile(Path(__file__).with_name("diagnostics.py"), tool_directory / "diagnostics.py")
        shutil.copyfile(Path(__file__).with_name("reclassify_audit.py"), tool_directory / "reclassify_audit.py")
    files = sorted(path for path in destination.rglob("*") if path.is_file())
    (destination / "SHA256SUMS").write_text("".join(f"{digest(path)}  {path.relative_to(destination)}\n" for path in files))
    print(f"Archived {len(inventory)} evidence files to {destination}")


if __name__ == "__main__":
    main()

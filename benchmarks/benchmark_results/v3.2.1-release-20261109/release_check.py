#!/usr/bin/env python3
"""Run finite release gates; this does not run biological panel campaigns.

Requires Cargo with rustfmt/Clippy, Python with PyYAML/jsonschema, and BLAST+.
Use --only to run the same gates separately in CI, or --offline with cached Cargo
dependencies. The default checks both supported hash-map backends.
"""

import argparse
import datetime
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional


ROOT = Path(__file__).resolve().parents[1]
CHECKS = ("metadata", "fmt", "clippy", "rust", "python", "references", "panels", "fixtures")


def require(condition: bool, message: str) -> None:
    """Stop the check with an actionable failure."""
    if not condition:
        raise ValueError(message)


def string_field(section: str, name: str) -> str:
    """Read one plain quoted Cargo field without requiring Python 3.11 tomllib."""
    values = re.findall(rf'^\s*{name}\s*=\s*"([^"\n]+)"\s*(?:#.*)?$', section, re.M)
    require(len(values) == 1, f"Expected one plain quoted Cargo {name} field")
    return values[0]


def check_metadata(tag: Optional[str] = None) -> None:
    """Require Cargo/lock agreement and changelog entries for the current version."""
    manifest = (ROOT / "Cargo.toml").read_text()
    package = re.search(r"(?ms)^\[package\]\s*\n(.*?)(?=^\[|\Z)", manifest)
    require(package is not None, "Cargo.toml has no [package] section")
    name = string_field(package.group(1), "name")
    version = string_field(package.group(1), "version")
    require(
        re.fullmatch(r"\d+\.\d+\.\d+(?:-dev)?", version) is not None,
        f"Unsupported release/development version: {version}",
    )
    lock_packages = re.split(
        r"(?m)^\[\[package\]\]\s*$", (ROOT / "Cargo.lock").read_text()
    )[1:]
    locked = [
        section for section in lock_packages if string_field(section, "name") == name
    ]
    require(
        len(locked) == 1 and string_field(locked[0], "version") == version,
        "Cargo.toml and the root package in Cargo.lock disagree",
    )

    changelog = (ROOT / "CHANGELOG.md").read_text()
    if version.endswith("-dev"):
        require(tag is None, "A release tag cannot use a development version")
        require(
            re.search(r"(?m)^## \[Unreleased\]\s*$", changelog) is not None,
            "Development version requires an Unreleased changelog section",
        )
    else:
        recipe_versions = re.findall(
            r'^\{% set version = "([^"\n]+)" %\}$',
            (ROOT / "meta.yaml").read_text(),
            re.M,
        )
        require(
            recipe_versions == [version],
            "Cargo.toml and the version in meta.yaml disagree",
        )
        releases = re.findall(
            r"(?m)^## \[(\d+\.\d+\.\d+)\] - (\d{4}-\d{2}-\d{2})\s*$",
            changelog,
        )
        require(
            bool(releases) and releases[0][0] == version,
            f"First dated CHANGELOG release must be [{version}]",
        )
        datetime.date.fromisoformat(releases[0][1])
        section = re.search(
            rf"(?ms)^## \[{re.escape(version)}\] - [^\n]+\n(.*?)(?=^## |\Z)",
            changelog,
        )
        require(
            section is not None
            and re.search(r"(?m)^- \S", section.group(1)) is not None,
            f"CHANGELOG [{version}] must contain release entries",
        )
        if tag is not None:
            require(tag == f"v{version}", f"Tag {tag!r} does not match v{version}")
    print(
        f"Metadata agrees: {name} {version}" + (f", tag {tag}" if tag else ""),
        flush=True,
    )


def run(
    command: list[str], capture: bool = False, timeout: Optional[int] = None
) -> subprocess.CompletedProcess:
    """Run a gate in the checkout and expose diagnostics on failure."""
    print("+ " + shlex.join(command), flush=True)
    result = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=capture, timeout=timeout
    )
    if result.returncode and capture:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
    result.check_returncode()
    return result


def cargo_flags(offline: bool) -> list[str]:
    return ["--locked"] + (["--offline"] if offline else [])


def build_binary(offline: bool) -> Path:
    """Obtain Cargo's actual default-backend executable, respecting target dirs."""
    result = run(
        [
            "cargo", "build", "--release", "--bin", "sharkmer",
            "--message-format=json", *cargo_flags(offline),
        ],
        capture=True,
    )
    print(result.stderr, end="", file=sys.stderr)
    artifacts = [json.loads(line) for line in result.stdout.splitlines()]
    binaries = [
        row["executable"]
        for row in artifacts
        if row.get("reason") == "compiler-artifact"
        and row.get("executable")
        and row["target"]["name"] == "sharkmer"
    ]
    require(len(binaries) == 1, "Cargo did not report exactly one sharkmer binary")
    return Path(binaries[0])


def normalize_yaml(value: Any) -> Any:
    """Match JSON Schema strings for YAML's implicit date values."""
    if isinstance(value, dict):
        return {str(key): normalize_yaml(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_yaml(item) for item in value]
    if isinstance(value, (datetime.date, datetime.datetime)):
        return value.isoformat()
    return value


def check_panels(binary: Path) -> None:
    """Check schemas, embedded exports, and loaders without running assemblies."""
    import jsonschema
    import yaml

    schema = json.loads((ROOT / "schemas/panel/v2.json").read_text())
    validator_type = jsonschema.validators.validator_for(schema)
    validator_type.check_schema(schema)
    validator = validator_type(schema)
    panels = sorted((ROOT / "panels").glob("*.yaml"))
    require(len(panels) == 9, "Expected all nine built-in panels")
    fixture = ROOT / "tests/fixtures/ERR571460_100k_R1.fastq.gz"
    with tempfile.TemporaryDirectory(prefix="sharkmer-panels-") as directory:
        for panel in panels + [ROOT / "panels/examples/reference.yaml"]:
            data = normalize_yaml(yaml.safe_load(panel.read_text()))
            validator.validate(data)
            modes = [["--pcr-panel-file", str(panel)]]
            if panel in panels:
                exported = run(
                    [str(binary), "--export-panel", panel.stem],
                    capture=True,
                    timeout=30,
                )
                require(
                    normalize_yaml(yaml.safe_load(exported.stdout)) == data,
                    f"Embedded/exported panel differs: {panel.stem}",
                )
                modes.append(["--pcr-panel", panel.stem])
            for index, mode in enumerate(modes):
                output = Path(directory) / f"{panel.stem}-{index}"
                run(
                    [
                        str(binary), "-k", "19", "-t", "2", "--chunks", "0",
                        "-s", "schema_probe", "-o", str(output), "--dry-run",
                        *mode, str(fixture),
                    ],
                    capture=True,
                    timeout=30,
                )
                require(
                    not any(output.rglob("*")),
                    f"Dry run wrote outputs: {panel.stem}",
                )
    print(
        "Panels passed: 10 schemas, 9 embedded exports, 19 loader dry runs",
        flush=True,
    )


def check_fixtures(binary: Path) -> None:
    """Check the pinned full-sequence CLI oracle; warnings are not biological truth."""
    import yaml

    truth = yaml.safe_load((ROOT / "benchmarks/known_truth.yaml").read_text())
    require(len(truth["samples"]) == 1, "Expected one pinned known-answer fixture")
    oracle = truth["samples"][0]
    require(
        oracle["accession"] == "ERR571460" and oracle["panel"] == "cnidaria",
        "Unexpected known-answer fixture",
    )
    fixture = ROOT / "tests/fixtures/ERR571460_100k_R1.fastq.gz"
    require(
        hashlib.sha256(fixture.read_bytes()).hexdigest() == oracle["input_sha256"],
        "Known-answer fixture checksum differs",
    )
    expected = {}
    for gene, records in oracle["expected"]["end-to-end"]["genes"].items():
        require(len(records) == 1, f"Expected one pinned product for {gene}")
        expected[f"cnidaria_{gene}"] = (records[0]["length"], records[0]["sha256"])
    with tempfile.TemporaryDirectory(prefix="sharkmer-fixture-") as directory:
        output = Path(directory)
        run(
            [
                str(binary), "-k", "31", "-t", "2", "--chunks", "0",
                "-s", "fixture", "-o", str(output), "--pcr-panel", "cnidaria",
                str(fixture),
            ],
            capture=True,
            timeout=300,
        )
        stats = yaml.safe_load((output / "fixture.stats.yaml").read_text())
        require(
            stats["n_reads_read"] == 100000 and stats["run_status"] == "complete",
            "Known-answer fixture read count or completion differs",
        )
        observed = {}
        for fasta in output.glob("*.fasta"):
            gene = fasta.stem.removeprefix("fixture_")
            for block in fasta.read_text().split(">"):
                if not block:
                    continue
                header, *lines = block.splitlines()
                if "assembly_status=warning" in header.split():
                    continue
                require(
                    gene not in observed,
                    f"Duplicate ordinary fixture product: {gene}",
                )
                sequence = "".join(lines)
                observed[gene] = (len(sequence), hashlib.sha256(sequence.encode()).hexdigest())
        require(observed == expected, "Known-answer full-sequence products differ")
    print(
        "Fixture passed: complete ordinary 18S and 28S_2 sequence hashes; "
        "Rust tests cover the counting oracle",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only", nargs="+", choices=CHECKS,
        help="Run selected gates, in the supplied order",
    )
    parser.add_argument(
        "--backend", choices=("both", "default", "fxhashmap"), default="both",
        help="Backends for Rust tests and Clippy (default: both)",
    )
    parser.add_argument(
        "--offline", action="store_true", help="Require cached Cargo dependencies"
    )
    parser.add_argument(
        "--tag", help="Require this exact vX.Y.Z tag in the metadata gate"
    )
    args = parser.parse_args()
    selected = args.only or CHECKS
    if args.tag and "metadata" not in selected:
        parser.error("--tag requires the metadata gate")
    backend_flags = {
        "default": [],
        "fxhashmap": ["--no-default-features", "--features", "fxhashmap"],
    }
    backends = backend_flags if args.backend == "both" else [args.backend]
    binary = None
    for check in selected:
        print(f"\n== {check} ==", flush=True)
        if check == "metadata":
            check_metadata(args.tag)
        elif check == "fmt":
            run(["cargo", "fmt", "--all", "--check"])
        elif check in ("clippy", "rust"):
            for backend in backends:
                command = (
                    ["cargo", "clippy", "--all-targets"]
                    if check == "clippy"
                    else ["cargo", "test", "--release"]
                )
                command += cargo_flags(args.offline) + backend_flags[backend]
                if check == "clippy":
                    command += ["--", "-D", "warnings"]
                run(command)
        elif check == "python":
            require(
                all(shutil.which(tool) for tool in ("blastn", "makeblastdb")),
                "BLAST+ (blastn and makeblastdb) is required so integration tests do not skip",
            )
            run([
                sys.executable, "-m", "unittest", "discover", "-s", "tests",
                "-p", "test_*.py", "-v",
            ])
        elif check == "references":
            run([sys.executable, "scripts/audit_panel_references.py", "--verify-existing"])
        else:
            if binary is None:
                binary = build_binary(args.offline)
            if check == "panels":
                check_panels(binary)
            else:
                check_fixtures(binary)
    print("\nRelease checks passed: " + ", ".join(selected), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print(f"Release check failed: {error}", file=sys.stderr)
        sys.exit(1)

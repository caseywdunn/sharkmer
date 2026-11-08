#!/usr/bin/env python3
"""Compare annotation-only changes against the frozen October 11 executable.

The eight invocations use the tiny repeat/clean fixture from the Rust CLI
regressions. This is output compatibility evidence, not biological validation
or a performance benchmark. Run inside the project's compute allocation.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASELINE_SHA256 = "3a6505da2750982c333b077a065ad9e0cd0e6607c6963fa7140a0b44da71a899"
SAMPLE = "annotation_parity"
CAUSES = {
    "omitted_self_loop_marker_occurrences": "homopolymer_length_unresolved",
    "cyclic_scc_marker_occurrences": "repeat_structure_unresolved",
    "retained_collision_marker_occurrences": "high_coverage_repeat_suspected",
}
MODES = {
    "inclusive": [],
    "standard_only": ["--no-emit-uncertain"],
    "inclusive_diagnostics": ["--diagnose-withheld-paths"],
    "standard_only_diagnostics": [
        "--no-emit-uncertain", "--diagnose-withheld-paths"
    ],
}


def require(condition: bool, message: str) -> None:
    """Fail without allowing Python optimization to disable a check."""
    if not condition:
        raise ValueError(message)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fasta(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            tokens = line[1:].split()
            fields: dict[str, str] = {}
            for token in tokens[1:]:
                require("=" in token, f"Malformed header token in {path}: {token}")
                key, value = token.split("=", 1)
                require(key not in fields, f"Duplicate header field in {path}: {key}")
                fields[key] = value
            records.append({
                "id": tokens[0], "tokens": tokens[1:], "fields": fields,
                "sequence": "",
            })
        elif line:
            require(bool(records), f"Sequence before header in {path}")
            records[-1]["sequence"] += line
    return records


def verify_outputs(directory: Path, expected_schema: int) -> dict[str, Any]:
    stats_path = directory / f"{SAMPLE}.stats.yaml"
    stats = yaml.safe_load(stats_path.read_text())
    manifest = yaml.safe_load((directory / f"{SAMPLE}.manifest.yaml").read_text())
    require(stats["sample"] == manifest["sample"] == SAMPLE, "Wrong sample")
    require(stats["run_id"] == manifest["run_id"], "Run identity mismatch")
    require(stats["run_status"] == manifest["status"] == "complete", "Incomplete run")
    require(manifest["schema_version"] == expected_schema, "Wrong manifest schema")
    require(stats["product_output_schema_version"] == expected_schema, "Wrong stats schema")
    require(stats["product_output_layout"] == {
        2: "combined_classified", 3: "combined_warnings"
    }[expected_schema], "Wrong output layout")
    receipts = manifest["files"]
    names = [receipt["path"] for receipt in receipts]
    expected_files = {stats_path.name}
    for target in stats["pcr_results"]:
        if target.get("output_file"):
            expected_files.add(target["output_file"])
    require(len(names) == len(set(names)) and set(names) == expected_files,
            "Manifest file set mismatch")
    for receipt in receipts:
        name = receipt["path"]
        require(Path(name).name == name, "Nonlocal receipt filename")
        require(digest(directory / name) == receipt["sha256"], f"Bad receipt: {name}")

    records_by_target: dict[str, list[dict[str, Any]]] = {}
    for target in stats["pcr_results"]:
        name = target["gene_name"]
        records = fasta(directory / target["output_file"]) if target.get("output_file") else []
        require(len(records) == target["n_products"], f"Count mismatch: {name}")
        require([len(record["sequence"]) for record in records] == target.get("product_lengths", []),
                f"Length mismatch: {name}")
        require(len({record["id"] for record in records}) == len(records),
                f"Duplicate sequence identifiers: {name}")
        standard_count = target["n_standard_products"]
        warning_count = target["n_uncertain_products"]
        require(standard_count + warning_count == len(records), f"Class counts mismatch: {name}")
        assessment = target.get("uncertain_products")
        products = assessment["products"] if assessment else []
        require(len(products) == warning_count, f"Evidence count mismatch: {name}")
        if assessment:
            require(assessment["schema_version"] == expected_schema, "Wrong evidence schema")
        for index, record in enumerate(records):
            fields = record["fields"]
            is_warning = index >= standard_count
            expected_status = ("uncertain" if is_warning else "standard") if expected_schema == 2 else (
                "warning" if is_warning else None
            )
            require(fields.get("assembly_status") == expected_status,
                    f"Wrong class/order: {record['id']}")
            record["warning"] = is_warning
            if not is_warning:
                require("assembly_status_causes" not in fields, "Orphan warning reasons")
                continue
            ordinal = index - standard_count
            product = products[ordinal]
            require(record["id"] == product["id"] and record["sequence"] == product["oriented_sequence"],
                    "FASTA/evidence identity mismatch")
            sequence_hash = hashlib.sha256(record["sequence"].encode("ascii")).hexdigest()
            require(product["sequence_sha256"] == sequence_hash, "Bad sequence evidence digest")
            causes = sorted({
                cause for source in product["sources"] for field, cause in CAUSES.items()
                if source[field] > 0
            })
            require(bool(causes), "Warning has no underlying cause")
            if expected_schema == 3:
                expected_fields = {
                    "sample": SAMPLE, "gene": name, "product": str(ordinal),
                    "length": str(len(record["sequence"])), "assembly_status": "warning",
                    "assembly_status_causes": ",".join(causes),
                }
                require(fields == expected_fields, f"Unexpected warning header: {record['id']}")
                require(product["assembly_status_causes"] == causes,
                        "Structured causes do not match source evidence")
        records_by_target[name] = records
    return {"stats": stats, "records": records_by_target}


def normalized_stats(stats: dict[str, Any]) -> dict[str, Any]:
    normalized = copy.deepcopy(stats)
    # Runtime provenance is captured per invocation, not compared as evidence.
    for field in ("command", "run_id", "peak_memory_bytes", "stage_timings",
                  "product_output_schema_version", "product_output_layout"):
        normalized.pop(field, None)
    for target in normalized["pcr_results"]:
        assessment = target.get("uncertain_products")
        if assessment:
            assessment.pop("schema_version")
            for product in assessment["products"]:
                product.pop("description")
                product.pop("assembly_status_causes", None)
    # All source markers, limits, counters, and hidden diagnostics stay intact.
    return normalized


def compare(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    require(normalized_stats(before["stats"]) == normalized_stats(after["stats"]),
            "Structured stats/evidence changed beyond allowed annotation fields")
    require(before["records"].keys() == after["records"].keys(), "Target set changed")
    totals = {"records": 0, "standard": 0, "warning": 0}
    for target, old_records in before["records"].items():
        new_records = after["records"][target]
        require(len(old_records) == len(new_records), f"Record count changed: {target}")
        for old, new in zip(old_records, new_records):
            for field in ("id", "sequence", "warning"):
                require(old[field] == new[field], f"{field} or record order changed: {target}")
            if not old["warning"]:
                original = [token for token in old["tokens"] if token != "assembly_status=standard"]
                require(original == new["tokens"], "Released ordinary header fields changed")
            totals["records"] += 1
            totals["warning" if old["warning"] else "standard"] += 1
    return {"passed": True, **totals, "full_structured_evidence_equal": True}


def run_one(binary: Path, directory: Path, flags: list[str], fixture: dict[str, Any],
            timeout: int) -> dict[str, Any]:
    directory.mkdir(parents=True)
    command = [str(binary), "-k", "19", "--threads", "1", "--chunks", "0",
               "-s", SAMPLE, "-o", str(directory)]
    for primer in fixture["primers"]:
        command.extend(["--pcr-primers", primer])
    command.extend([*flags, str(HERE / "inputs" / "reads.fastq")])
    started = datetime.now(timezone.utc).isoformat()
    tick = time.monotonic()
    returncode = None
    timed_out = False
    with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
        try:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
            returncode = completed.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    receipt = {
        "command": command, "started_utc": started,
        "wall_seconds": time.monotonic() - tick, "returncode": returncode,
        "timeout_seconds": timeout, "timed_out": timed_out,
    }
    write_json(directory / "invocation.json", receipt)
    require(not timed_out and returncode == 0, f"Invocation failed or timed out: {directory}")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=REPO / "benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/results-v1/bin/sharkmer")
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    result: dict[str, Any] = {"purpose": "annotation compatibility, not biological/performance validation", "passed": False,
                              "modes": {}, "invocations": []}
    try:
        binaries = {"baseline": args.baseline.resolve(), "candidate": args.candidate.resolve()}
        binding_paths = {**binaries, "harness": Path(__file__).resolve(),
                         "input": HERE / "inputs/reads.fastq", "fixture": HERE / "inputs/fixture.json"}
        bindings = {name: {"path": str(path), "sha256": digest(path)} for name, path in binding_paths.items()}
        result["bindings"] = bindings
        require(bindings["baseline"]["sha256"] == BASELINE_SHA256, "Baseline is not the frozen October 11 binary")
        fixture = json.loads((HERE / "inputs/fixture.json").read_text())
        require(bindings["input"]["sha256"] == fixture["input_sha256"], "Stored fixture input changed")
        for mode, flags in MODES.items():
            observed = {}
            for label, schema in (("baseline", 2), ("candidate", 3)):
                directory = output / mode / label
                receipt = run_one(binaries[label], directory, flags, fixture, args.timeout)
                result["invocations"].append({"mode": mode, "label": label, **receipt})
                observed[label] = verify_outputs(directory, schema)
                require(observed[label]["stats"]["n_reads_read"] == fixture["reads"],
                        "Fixture read count changed")
                clean = observed[label]["records"]["clean"]
                require(len(clean) == 1 and clean[0]["sequence"] == fixture["clean_sequence"],
                        "Clean positive-control sequence was not reconstructed exactly")
            result["modes"][mode] = compare(observed["baseline"], observed["candidate"])
            require(result["modes"][mode]["standard"] == 1, "Clean positive control missing")
            require((result["modes"][mode]["warning"] > 0) == ("standard_only" not in mode),
                    "Repeat warning control or opt-out failed")
            write_json(output / "result.json", result)
        require(all(digest(path) == bindings[name]["sha256"] for name, path in binding_paths.items()),
                "Executable, input, or harness changed during comparison")
        result["bindings_unchanged"] = True
        result["passed"] = True
    except (OSError, ValueError, KeyError, TypeError, IndexError, yaml.YAMLError) as error:
        result["failure"] = str(error)
    write_json(output / "result.json", result)
    print(json.dumps({"passed": result["passed"], "result": str(output / "result.json"),
                      "failure": result.get("failure")}, sort_keys=True))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

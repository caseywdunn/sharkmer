#!/usr/bin/env python3
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import yaml


REPOSITORY = Path("/home/claude/repos/sharkmer")
FROZEN_SUPPLEMENTAL = Path("/tmp/sharkmer-release-comparison/supplemental_synthetic_probes.py")
FROZEN_RUNNER_ROOT = Path("/tmp/sharkmer-release-comparison/sources/candidate/scripts")
SUPPLEMENTAL_SHA256 = "40ba6e4d69c5d2d72c41f5c3510c9af528208fced33f38cf0957c6afd2bbda68"
FIXTURE = REPOSITORY / "tests/fixtures/ERR571460_100k_R1.fastq.gz"
KNOWN_TRUTH = REPOSITORY / "benchmarks/known_truth.yaml"
PANEL = REPOSITORY / "panels/cnidaria.yaml"
TIMEOUT_SECONDS = 1200


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def receipt(path):
    path = Path(path).resolve()
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Receipt source is not a regular file: {path}")
    return {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def load_builds(path):
    document = json.loads(Path(path).read_text())
    versions = document.get("versions") if isinstance(document, dict) else None
    if set(versions or ()) != {"baseline", "candidate"}:
        raise ValueError("Build manifest must contain exactly baseline and candidate")
    for name, build in versions.items():
        binary = Path(build.get("binary_path", ""))
        expected = build.get("binary_sha256")
        if binary.is_symlink() or not binary.is_file() or sha256(binary) != expected:
            raise ValueError(f"Binary receipt invalid for {name}")
    return document


def load_frozen_sources():
    if sha256(FROZEN_SUPPLEMENTAL) != SUPPLEMENTAL_SHA256:
        raise ValueError("Frozen supplemental probe source changed")
    runner_path = FROZEN_RUNNER_ROOT / "sharkmer_validate/runner.py"
    if not runner_path.is_file():
        raise ValueError("Frozen validation runner is unavailable")
    sys.path.insert(0, str(FROZEN_RUNNER_ROOT))
    try:
        from sharkmer_validate import runner
    finally:
        sys.path.pop(0)
    return load_module("frozen_supplemental", FROZEN_SUPPLEMENTAL), runner, runner_path


def as_text(value):
    return value.decode(errors="replace") if isinstance(value, bytes) else value or ""


def write_logs(directory, completed, timeout):
    (directory / "stdout.log").write_text(as_text(completed.stdout))
    (directory / "stderr.log").write_text(as_text(completed.stderr))
    write_json(directory / "execution.json", {
        "returncode": completed.returncode,
        "timed_out": timeout,
        "timeout_seconds": TIMEOUT_SECONDS,
    })


def parse_products(runner, stats, sample, directory):
    transaction = runner._validate_output_manifest(stats, sample, directory)
    files = [entry["output_file"] for entry in stats.get("pcr_results", []) if entry.get("status") == "success"]
    groups = runner.parse_fasta_products(sample, directory, files)
    if runner._validate_output_manifest(stats, sample, directory) != transaction:
        raise ValueError("Output transaction changed during product parsing")
    return groups


def parse_legacy_fasta(path, sample):
    records, header, sequence = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "sequence": "".join(sequence), "length": len("".join(sequence))})
            header, sequence = line[1:], []
        else:
            sequence.append(line)
    if header is not None:
        records.append({"header": header, "sequence": "".join(sequence), "length": len("".join(sequence))})
    gene = Path(path).stem.removeprefix(f"{sample}_")
    return {"gene": gene, "n_products": len(records), "lengths": [record["length"] for record in records], "products": records, "output_file": Path(path).name}


def validate_legacy(runner, stats, sample, directory, command, expected_genes):
    if any(field in stats for field in ("run_id", "run_status", "output_manifest", "input_source", "stage_timings")):
        raise ValueError("Legacy output contains current-only transaction fields")
    if stats.get("sample") != sample or stats.get("kmer_length") != int(command[command.index("-k") + 1]) or stats.get("command") != " ".join(command):
        raise ValueError("Legacy stats do not match command")
    results = stats.get("pcr_results")
    if not isinstance(results, list) or {entry.get("gene_name") for entry in results} != expected_genes:
        raise ValueError("Legacy stats gene set differs from expected genes")
    output_files = []
    for entry in results:
        if entry.get("status") not in {"success", "fail"} or type(entry.get("n_products")) is not int:
            raise ValueError("Legacy PCR status or product count is invalid")
        if entry["status"] == "success":
            file_name = f"{sample}_{entry['gene_name']}.fasta"
            if entry["n_products"] <= 0 or not isinstance(entry.get("product_lengths"), list):
                raise ValueError("Legacy successful PCR result is invalid")
            output_files.append(file_name)
        elif entry["n_products"] != 0:
            raise ValueError("Legacy failed PCR result reports products")
    expected_files = {f"{sample}.stats.yaml", *output_files}
    observed_files = {path.name for path in directory.iterdir() if path.is_file() and not path.is_symlink()}
    if observed_files != expected_files:
        raise ValueError("Legacy output file set differs from stats")
    groups = [parse_legacy_fasta(directory / file_name, sample) for file_name in output_files]
    by_file = {group["output_file"]: group for group in groups}
    for entry in results:
        if entry["status"] == "success":
            group = by_file.get(f"{sample}_{entry['gene_name']}.fasta")
            if group is None or group["n_products"] != entry["n_products"] or group["lengths"] != entry["product_lengths"]:
                raise ValueError("Legacy FASTA differs from stats")
    return groups


def invoke(runner, binary, version, label, command, expected_genes, scope, directory):
    directory.mkdir(parents=True)
    artifacts = directory.parent / f"{directory.name}.artifacts"
    artifacts.mkdir()
    write_json(artifacts / "command.json", command)
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
        write_logs(artifacts, completed, False)
    except subprocess.TimeoutExpired as error:
        completed = subprocess.CompletedProcess(command, None, as_text(error.stdout), as_text(error.stderr))
        write_logs(artifacts, completed, True)
        return {"version": version, "label": label, "command": command, "success": False, "failure": "timeout", "raw_output_directory": str(directory), "artifacts": str(artifacts)}
    result = {"version": version, "label": label, "command": command, "returncode": completed.returncode, "raw_output_directory": str(directory), "artifacts": str(artifacts)}
    if completed.returncode != 0:
        result.update({"success": False, "failure": "nonzero_exit"})
        return result
    try:
        sample = command[command.index("-s") + 1]
        stats = runner._parse_stats_yaml(directory / f"{sample}.stats.yaml")
        if version == "baseline":
            groups = validate_legacy(runner, stats, sample, directory, command, expected_genes)
        else:
            runner._validate_stats_manifest(stats, sample, int(command[command.index("-k") + 1]), expected_genes, scope, expected_command=command)
            groups = parse_products(runner, stats, sample, directory)
        result.update({"success": True, "stats": stats, "products": groups})
    except (KeyError, OSError, TypeError, ValueError, yaml.YAMLError) as error:
        result.update({"success": False, "failure": "invalid_current_manifest", "message": str(error)})
    return result


def synthetic_expectation(case, result):
    if not result.get("success"):
        return False
    observed = [product["sequence"] for group in result["products"] for product in group["products"]]
    if case["case"] in {"threshold", "repeat-a18"}:
        return observed == case["expected_sequences"]
    diagnostics = json.dumps(result["stats"].get("pcr_results", [])).lower()
    return not observed and "repeat" in diagnostics


def known_truth_expected():
    document = yaml.safe_load(KNOWN_TRUTH.read_text())
    sample = document["samples"][0]
    return sample["expected"]


def exact_products(result):
    observed = {}
    for group in result.get("products", []):
        gene = group["gene"].removeprefix("cnidaria_")
        observed[gene] = sorted((product["length"], hashlib.sha256(product["sequence"].encode()).hexdigest()) for product in group["products"])
    return observed


def main():
    parser = argparse.ArgumentParser(description="Run bounded high-copy and pinned exact controls")
    parser.add_argument("--builds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--synthetic-only", action="store_true")
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise SystemExit(f"Output already exists: {arguments.output}")
    builds = load_builds(arguments.builds)
    supplemental, runner, runner_path = load_frozen_sources()
    expected = known_truth_expected()
    if sha256(FIXTURE) != yaml.safe_load(KNOWN_TRUTH.read_text())["samples"][0]["input_sha256"]:
        raise ValueError("Pinned ERR571460 fixture checksum differs from known truth")
    arguments.output.mkdir(parents=True)
    receipts = {name: receipt(path) for name, path in {
        "builds": arguments.builds, "controls": __file__, "supplemental": FROZEN_SUPPLEMENTAL,
        "frozen_runner": runner_path, "known_truth": KNOWN_TRUTH, "fixture": FIXTURE, "panel": PANEL,
    }.items()}
    write_json(arguments.output / "receipts.json", receipts)
    input_directory = arguments.output / "synthetic_inputs"
    input_directory.mkdir()
    cases = [supplemental.make_threshold_input(input_directory), *supplemental.make_repeat_inputs(input_directory)]
    if len(cases[0]["expected_sequences"]) != 1 or len(cases[0]["expected_sequences"][0]) != 180:
        raise ValueError("Frozen threshold control no longer specifies the 180 bp product")
    results = []
    for version, build in builds["versions"].items():
        binary = Path(build["binary_path"]).resolve()
        for case in cases:
            directory = arguments.output / "synthetic" / version / case["case"]
            command = [str(binary), "-s", case["sample"], "-o", str(directory), "-k", "19", "-t", "2", "--pcr-primers", case["primer"], str(case["input_path"])]
            result = invoke(runner, binary, version, case["case"], command, {"target"}, "end-to-end", directory)
            result["input"] = receipt(case["input_path"])
            result["candidate_expectation_pass"] = synthetic_expectation(case, result) if version == "candidate" else None
            results.append(result)
    candidate = Path(builds["versions"]["candidate"]["binary_path"]).resolve()
    for label, kmer_length, threaded in (() if arguments.synthetic_only else (("count-k19", 19, False), ("pcr-k31", 31, False), ("pcr-k31-read-threading", 31, True))):
        directory = arguments.output / "known_truth" / label
        sample = f"known_{label}"
        command = [str(candidate), "-k", str(kmer_length), "-t", "2", "--max-reads", "100000", "-o", str(directory), "-s", sample]
        scope = "counting-only" if label == "count-k19" else "end-to-end"
        panel_data = runner.load_panel_yaml(PANEL)
        prefix = panel_data.get("gene_prefix") or panel_data["name"]
        expected_genes = set() if scope == "counting-only" else {f"{prefix}_{gene}" for gene in runner.panel_gene_names(panel_data)}
        if scope == "end-to-end":
            command.extend(["--pcr-panel-file", str(PANEL)])
        if threaded:
            command.append("--read-threading")
        command.append(str(FIXTURE))
        result = invoke(runner, candidate, "candidate", label, command, expected_genes, scope, directory)
        if result.get("success"):
            if scope == "counting-only":
                result["exact_oracle_pass"] = result["stats"].get("n_kmers") == expected["counting-only"]["n_kmers"]
            else:
                desired = {gene: sorted((entry["length"], entry["sha256"]) for entry in products) for gene, products in expected["end-to-end"]["genes"].items()}
                result["exact_oracle_pass"] = exact_products(result) == desired
        else:
            result["exact_oracle_pass"] = False
        result["input"] = receipt(FIXTURE)
        results.append(result)
    failed = [result for result in results if not result.get("success") or result.get("candidate_expectation_pass") is False or result.get("exact_oracle_pass") is False]
    write_json(arguments.output / "summary.json", {"schema_version": 1, "label": "bounded supplemental controls; not timing or biological calibration", "builds": builds, "receipts": receipts, "results": results, "failed_labels": [result["label"] for result in failed]})
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()

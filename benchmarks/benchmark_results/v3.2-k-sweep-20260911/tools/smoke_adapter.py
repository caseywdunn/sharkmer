#!/usr/bin/env python3
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
DRIVER_PATH = ROOT / "driver.py"
BUILDS_PATH = Path("/tmp/sharkmer-high-copy-20260911/final-builds.json")
INPUT_PATH = Path("/tmp/sharkmer-output-smoke.fastq")
PANEL_PATH = Path("/tmp/sharkmer-release-comparison/sources/candidate/panels/insecta.yaml")


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def run_smoke(output_root):
    driver = load_module("k_sweep_smoke_driver", DRIVER_PATH)
    driver.verify_dependencies()
    if output_root.exists():
        raise ValueError("Smoke output exists; preserve the prior run")
    output_root.mkdir(parents=True)
    builds_document = driver.load_json(BUILDS_PATH)
    protocol = driver.load_json(ROOT / "discovery-protocol.json")
    input_lines = INPUT_PATH.read_text().splitlines()
    if len(input_lines) != 4 or len(input_lines[1]) != len(input_lines[3]):
        raise ValueError("Smoke input must contain exactly one valid FASTQ record")
    results = []
    for selected_k in protocol["allowed_k"]:
        helper = driver.load_module(
            f"smoke_helper_{selected_k}", driver.EXPECTED_DEPENDENCIES["helper"][0]
        )
        legacy = driver.load_module(
            f"smoke_legacy_{selected_k}",
            driver.EXPECTED_DEPENDENCIES["legacy_adapter"][0],
        )
        driver.install_dynamic_adapters(helper, legacy, selected_k)
        builds = driver.verify_builds(helper, protocol, builds_document)
        validator_build = driver.load_json(
            driver.EXPECTED_DEPENDENCIES["original_builds"][0]
        )["versions"]["candidate"]
        runner, unused_blast, unused_validator = helper.load_validator(
            validator_build["source_export"]
        )
        panel_data = runner.load_panel_yaml(PANEL_PATH)
        prefix = panel_data.get("gene_prefix") or panel_data["name"]
        panel = {
            "path": str(PANEL_PATH),
            "sha256": driver.sha256_file(PANEL_PATH),
            "data": panel_data,
            "output_prefix": prefix,
            "stats_genes": {
                f"{prefix}_{gene}" for gene in runner.panel_gene_names(panel_data)
            },
        }
        for version in ("baseline", "candidate"):
            cell_root = output_root / f"k{selected_k}" / version
            output_dir = cell_root / "output"
            output_dir.mkdir(parents=True)
            sample = f"insecta_smoke_k{selected_k}_{version}"
            command = [
                builds[version]["binary_path"],
                "-k",
                str(selected_k),
                "-t",
                "2",
                "--chunks",
                "0",
                "--max-reads",
                "1",
                "-o",
                str(output_dir) + "/",
                "-s",
                sample,
                "--pcr-panel-file",
                str(PANEL_PATH),
                str(INPUT_PATH),
            ]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=120)
            (cell_root / "stdout.log").write_text(completed.stdout)
            (cell_root / "stderr.log").write_text(completed.stderr)
            if completed.returncode != 0:
                raise ValueError(f"Smoke command failed for k{selected_k}/{version}")
            stats_path = output_dir / f"{sample}.stats.yaml"
            stats = helper.parse_stats(stats_path)
            invocation = {
                "version": version,
                "depth": 1,
                "sample_prefix": sample,
            }
            input_record = {
                "path": str(INPUT_PATH),
                "prefixes": {1: {"records": 1, "bases": len(input_lines[1])}},
            }
            expected_version = builds[version]["binary_version"].split()[1]
            pcr_results = helper.validate_shared_stats(
                stats,
                invocation,
                panel,
                input_record,
                command,
                expected_version,
            )
            if version == "baseline":
                genes, completion = helper.validate_legacy_outputs(
                    output_dir, invocation, panel, stats, pcr_results
                )
            else:
                genes, completion = helper.validate_current_outputs(
                    output_dir, invocation, panel, stats, pcr_results, command, runner
                )
            wrong_k = 23 if selected_k == 19 else 19
            wrong_helper = driver.load_module(
                f"wrong_smoke_helper_{selected_k}_{version}",
                driver.EXPECTED_DEPENDENCIES["helper"][0],
            )
            wrong_legacy = driver.load_module(
                f"wrong_smoke_legacy_{selected_k}_{version}",
                driver.EXPECTED_DEPENDENCIES["legacy_adapter"][0],
            )
            driver.install_dynamic_adapters(wrong_helper, wrong_legacy, wrong_k)
            wrong_k_rejected = False
            try:
                wrong_helper.validate_shared_stats(
                    stats,
                    invocation,
                    panel,
                    input_record,
                    command,
                    expected_version,
                )
            except ValueError as error:
                wrong_k_rejected = "selected sweep k" in str(error)
            if not wrong_k_rejected:
                raise ValueError("Deliberate wrong-k stats validation was not rejected")
            results.append(
                {
                    "k": selected_k,
                    "version": version,
                    "command": command,
                    "returncode": completed.returncode,
                    "stats": helper.file_receipt(stats_path),
                    "outputs": driver.inventory_receipts(helper, output_dir),
                    "n_genes": len(genes),
                    "completion": completion,
                    "wrong_k_rejected": True,
                }
            )
    summary = {
        "schema_version": 1,
        "purpose": "Nonbenchmark dynamic-k adapter smoke; timings are not evidence",
        "driver": {
            "path": str(DRIVER_PATH),
            "sha256": driver.sha256_file(DRIVER_PATH),
        },
        "builds": {
            "path": str(BUILDS_PATH),
            "sha256": driver.sha256_file(BUILDS_PATH),
        },
        "input": {
            "path": str(INPUT_PATH),
            "sha256": driver.sha256_file(INPUT_PATH),
        },
        "panel": {
            "path": str(PANEL_PATH),
            "sha256": driver.sha256_file(PANEL_PATH),
        },
        "results": results,
    }
    (output_root / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run_smoke(arguments.output)


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    main()

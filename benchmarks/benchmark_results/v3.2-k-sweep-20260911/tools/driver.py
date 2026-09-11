#!/usr/bin/env python3
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


HIGH_COPY = Path("/tmp/sharkmer-high-copy-20260911")
RELEASE_COMPARISON = Path("/tmp/sharkmer-release-comparison")
EXPECTED_DEPENDENCIES = {
    "helper": (
        RELEASE_COMPARISON / "execution/receipts/driver.py",
        "f0d7eee4631819921eadb6b4dddf05b422c5e6eaf4c554db82259cd1c41bec97",
    ),
    "legacy_adapter": (
        RELEASE_COMPARISON / "postprocess_legacy.py",
        "f9e58a665c530b71c96ad7fb79508ba1598df7119feb7c52ee403aeb3e790aef",
    ),
    "high_copy_driver": (
        HIGH_COPY / "benchmark.py",
        "10d199dcaa7610d29dfde61d92b4646d8abc974065b1cee6f5ef02d800ae7b1d",
    ),
    "high_copy_protocol": (
        HIGH_COPY / "final-protocol.json",
        "8b9e1ad915f0e1ebde5fe904d2e938d408c218bb6ff17c7da73a7986677f7f5b",
    ),
    "high_copy_builds": (
        HIGH_COPY / "final-builds.json",
        "c705b80b311f4874cb3258ce8a203dfa3f7478ac771ad6e250fd7630e70e9c79",
    ),
    "original_protocol": (
        RELEASE_COMPARISON / "protocol.json",
        "097485889647071ce0ec7de56ceb9030f90a6ad27ddacbf5cb37cc2a37cbb825",
    ),
    "inputs": (
        RELEASE_COMPARISON / "inputs.json",
        "4763732cebacd43d3e295afe17ff903bdd1a6fc07449de24a03cd263b692cacf",
    ),
    "original_builds": (
        RELEASE_COMPARISON / "builds.json",
        "475518df001573fd5f4a05316a16d6b29e4bac4d402bd49530a4207d9f32cec7",
    ),
}
EXPECTED_ORDERING = {
    "outer_order": ["pair_index", "protocol_sample_order"],
    "k_rotation": "left_by_(pair_index_minus_1_plus_zero_based_sample_index)_modulo_k_count",
    "version_order": "baseline_first_when_(pair_index_plus_zero_based_sample_index_plus_allowed_k_index)_modulo_2_equals_1",
    "same_k_version_pair_contiguous": True,
}
ALLOWED_K_VALUES = {19, 23, 27, 31}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as input_stream:
        for block in iter(lambda: input_stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    with Path(path).open() as input_stream:
        value = json.load(input_stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def require_string(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Expected nonblank string for {label}")
    return value


def verify_dependencies():
    for label, (path, expected_sha256) in EXPECTED_DEPENDENCIES.items():
        if not path.is_file() or path.is_symlink() or sha256_file(path) != expected_sha256:
            raise ValueError(f"Frozen dependency changed: {label}")


def install_dynamic_adapters(helper, legacy, selected_k):
    original_shared = helper.validate_shared_stats
    original_legacy = helper.validate_legacy_outputs
    original_current = helper.validate_current_outputs

    def validate_shared(stats, invocation, panel, input_record, command, expected_version):
        if stats.get("kmer_length") != selected_k:
            raise ValueError("Stats k differs from selected sweep k")
        adapted_stats = dict(stats)
        adapted_stats["kmer_length"] = 19
        results = original_shared(
            adapted_stats, invocation, panel, input_record, command, expected_version
        )
        if invocation["version"] != "baseline":
            expected_source = {
                "kind": "local_files",
                "inputs": [input_record["path"]],
                "paired": False,
                "max_reads": invocation["depth"],
            }
            if stats.get("input_source") != expected_source:
                raise ValueError("Current stats input source differs from frozen invocation")
        return results

    def validate_legacy(output_dir, invocation, panel, stats, pcr_results):
        sample = invocation["sample_prefix"]
        expected_genes = {
            f"{sample}_{entry['gene_name']}.fasta": entry["gene_name"]
            for entry in pcr_results
            if entry["status"] == "success"
        }
        original_parser = helper.parse_fasta_file
        helper.parse_fasta_file = lambda path: legacy.parse_legacy_fasta(
            path, sample, expected_genes[Path(path).name]
        )
        try:
            return original_legacy(output_dir, invocation, panel, stats, pcr_results)
        finally:
            helper.parse_fasta_file = original_parser

    def validate_current(output_dir, invocation, panel, stats, pcr_results, command, runner):
        runner._validate_stats_manifest(
            stats,
            invocation["sample_prefix"],
            selected_k,
            panel["stats_genes"],
            "end-to-end",
            expected_command=command,
        )
        adapted_stats = dict(stats)
        adapted_stats["kmer_length"] = 19
        return original_current(
            output_dir, invocation, panel, adapted_stats, pcr_results, command, runner
        )

    helper.validate_shared_stats = validate_shared
    helper.validate_legacy_outputs = validate_legacy
    helper.validate_current_outputs = validate_current


def validate_protocol(protocol):
    if protocol.get("schema_version") != 1:
        raise ValueError("Unsupported protocol schema")
    require_string(protocol.get("purpose"), "purpose")
    allowed_k = protocol.get("allowed_k")
    if (
        not isinstance(allowed_k, list)
        or not allowed_k
        or any(type(value) is not int or value not in ALLOWED_K_VALUES for value in allowed_k)
        or len(allowed_k) != len(set(allowed_k))
    ):
        raise ValueError("allowed_k must be a unique nonempty list drawn from 19, 23, 27, 31")
    if protocol.get("ordering") != EXPECTED_ORDERING:
        raise ValueError("Protocol ordering differs from the preregistered interleave")
    original = load_json(EXPECTED_DEPENDENCIES["high_copy_protocol"][0])
    expected_settings = {key: value for key, value in original["settings"].items() if key != "k"}
    if protocol.get("settings") != expected_settings:
        raise ValueError("Non-k settings differ from the frozen high-copy comparison")
    if protocol.get("baseline_role") != original["baseline_role"]:
        raise ValueError("Baseline role differs from the frozen comparison")
    if protocol.get("expected_commits") != original["expected_commits"]:
        raise ValueError("Expected commits differ from the frozen comparison")
    samples = protocol.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ValueError("Protocol samples must be a nonempty list")
    original_samples = {
        (sample["panel"], sample["input"], sample["depth"]): sample
        for sample in original["samples"]
    }
    seen_cells = set()
    for sample in samples:
        if not isinstance(sample, dict):
            raise ValueError("Protocol sample is not an object")
        cell = (sample.get("panel"), sample.get("input"), sample.get("depth"))
        if cell in seen_cells or sample != original_samples.get(cell):
            raise ValueError("Sample differs from or duplicates the frozen high-copy comparison")
        seen_cells.add(cell)
    selected_panels = {sample["panel"] for sample in samples}
    expected_metadata = {
        panel: original["target_metadata"][panel] for panel in sorted(selected_panels)
    }
    if protocol.get("target_metadata") != expected_metadata:
        raise ValueError("Target metadata differs from the frozen high-copy comparison")
    fixed_reference = protocol.get("fixed_reference")
    if not isinstance(fixed_reference, dict):
        raise ValueError("Fixed k19 reference metadata is missing")
    if (
        fixed_reference.get("protocol_sha256")
        != EXPECTED_DEPENDENCIES["high_copy_protocol"][1]
        or fixed_reference.get("baseline_k") != 19
    ):
        raise ValueError("Fixed k19 reference identity differs from the frozen comparison")
    diagnostic_path = Path(fixed_reference.get("analysis", "")) / "diagnostic-audit.json"
    if (
        not diagnostic_path.is_file()
        or sha256_file(diagnostic_path) != fixed_reference.get("diagnostic_audit_sha256")
    ):
        raise ValueError("Fixed k19 diagnostic audit changed")
    limitations = protocol.get("limitations")
    if not isinstance(limitations, list) or not limitations or not all(
        isinstance(value, str) and value.strip() for value in limitations
    ):
        raise ValueError("Protocol limitations must be explicit")
    return allowed_k


def schedule_for_protocol(protocol):
    allowed_k = protocol["allowed_k"]
    baseline_role = protocol["baseline_role"]
    schedule = []
    maximum_pairs = max(sample["pairs"] for sample in protocol["samples"])
    for pair_index in range(1, maximum_pairs + 1):
        for sample_index, sample in enumerate(protocol["samples"]):
            if pair_index > sample["pairs"]:
                continue
            rotation = (pair_index - 1 + sample_index) % len(allowed_k)
            rotated_k = allowed_k[rotation:] + allowed_k[:rotation]
            for selected_k in rotated_k:
                allowed_index = allowed_k.index(selected_k)
                baseline_first = (
                    pair_index + sample_index + allowed_index
                ) % 2 == 1
                order = (
                    (baseline_role, "candidate")
                    if baseline_first
                    else ("candidate", baseline_role)
                )
                cell = f"{sample['panel']}/{sample['input']}/{sample['depth']}"
                for order_position, version in enumerate(order, start=1):
                    schedule.append(
                        {
                            "invocation_id": (
                                f"{sample_index:03d}_{cell.replace('/', '_')}_k{selected_k}_"
                                f"pair{pair_index}_{version}"
                            ),
                            "cell": cell,
                            "cell_index": sample_index,
                            "panel": sample["panel"],
                            "input": sample["input"],
                            "taxon": sample["taxon"],
                            "depth": sample["depth"],
                            "k": selected_k,
                            "pair_index": pair_index,
                            "order_position": order_position,
                            "version": version,
                        }
                    )
    return schedule


def verify_builds(helper, protocol, builds_document):
    builds = builds_document.get("versions")
    if not isinstance(builds, dict) or set(builds) != {"baseline", "candidate"}:
        raise ValueError("Expected baseline and candidate builds")
    if {
        version: build.get("commit") for version, build in builds.items()
    } != protocol["expected_commits"]:
        raise ValueError("Build commits differ from protocol")
    helper.EXPECTED_COMMITS = protocol["expected_commits"]
    verified = helper.verify_builds(builds_document)
    compared = list(verified.values())
    for field in (
        "rustc",
        "cargo",
        "cargo_artifact_features",
        "cargo_artifact_profile",
    ):
        if compared[0][field] != compared[1][field]:
            raise ValueError(f"Build configurations differ: {field}")
    return verified


def prepare_context(protocol_path, builds_path, selected_k):
    verify_dependencies()
    protocol = load_json(protocol_path)
    allowed_k = validate_protocol(protocol)
    if selected_k not in allowed_k:
        raise ValueError("Selected k is not preregistered")
    if sha256_file(builds_path) != EXPECTED_DEPENDENCIES["high_copy_builds"][1]:
        raise ValueError("Build receipt differs from the frozen clean builds")
    helper = load_module(
        f"k_sweep_helper_{selected_k}", EXPECTED_DEPENDENCIES["helper"][0]
    )
    legacy = load_module(
        f"k_sweep_legacy_{selected_k}", EXPECTED_DEPENDENCIES["legacy_adapter"][0]
    )
    install_dynamic_adapters(helper, legacy, selected_k)
    builds = verify_builds(helper, protocol, load_json(builds_path))
    validator_build = load_json(EXPECTED_DEPENDENCIES["original_builds"][0])["versions"]["candidate"]
    if helper.source_tree_sha256(Path(validator_build["source_export"])) != validator_build[
        "source_tree_sha256_after_build"
    ]:
        raise ValueError("Frozen validator source changed")
    runner, blast, validator = helper.load_validator(validator_build["source_export"])
    input_document = load_json(EXPECTED_DEPENDENCIES["inputs"][0])
    selected_inputs = {sample["input"] for sample in protocol["samples"]}
    input_document["inputs"] = [
        record for record in input_document["inputs"] if record["id"] in selected_inputs
    ]
    inputs = helper.verify_inputs(input_document)
    if set(inputs) != selected_inputs:
        raise ValueError("Selected input set is incomplete")
    original_protocol = load_json(EXPECTED_DEPENDENCIES["original_protocol"][0])
    panel_records = {
        record["name"]: record for record in original_protocol["panels"]
    }
    panels = {}
    for panel_name in {sample["panel"] for sample in protocol["samples"]}:
        panel_record = panel_records[panel_name]
        if helper.sha256_file(panel_record["path"]) != panel_record["sha256"]:
            raise ValueError("Frozen panel changed")
        panel_data = runner.load_panel_yaml(Path(panel_record["path"]))
        metadata_names = [entry["gene"] for entry in protocol["target_metadata"][panel_name]]
        if len(metadata_names) != len(set(metadata_names)) or set(metadata_names) != set(
            runner.panel_gene_names(panel_data)
        ):
            raise ValueError("Target metadata differs from panel genes")
        prefix = panel_data.get("gene_prefix") or panel_data["name"]
        panels[panel_name] = {
            **panel_record,
            "data": panel_data,
            "output_prefix": prefix,
            "stats_genes": {
                f"{prefix}_{gene}" for gene in runner.panel_gene_names(panel_data)
            },
        }
    settings = {**protocol["settings"], "k": selected_k}
    helper.validate_cpu_list(settings["cpu_list"], settings["threads"])
    return {
        "protocol": protocol,
        "helper": helper,
        "runner": runner,
        "blast": blast,
        "validator": validator,
        "builds": builds,
        "inputs": inputs,
        "input_document": input_document,
        "panels": panels,
        "settings": settings,
    }


def fixed_k_schedule(master_schedule, selected_k):
    return [copy.deepcopy(item) for item in master_schedule if item["k"] == selected_k]


def freeze_execution_receipts(context, protocol_path, builds_path, output_root):
    helper = context["helper"]
    receipts_dir = output_root / "receipts"
    receipts = {
        "driver": helper.freeze_receipt(__file__, receipts_dir / "driver.py"),
        "protocol": helper.freeze_receipt(protocol_path, receipts_dir / "protocol.json"),
        "builds": helper.freeze_receipt(builds_path, receipts_dir / "builds.json"),
    }
    for label, (source, unused_sha256) in EXPECTED_DEPENDENCIES.items():
        suffix = source.suffix or ".receipt"
        receipts[label] = helper.freeze_receipt(
            source, receipts_dir / f"dependency-{label}{suffix}"
        )
    fixed_reference = context["protocol"]["fixed_reference"]
    reference_paths = {
        "k19_execution_provenance": Path(fixed_reference["execution"]) / "provenance.json",
        "k19_execution_comparison": Path(fixed_reference["execution"]) / "comparison.json",
        "k19_analysis": Path(fixed_reference["analysis"]) / "analysis.json",
        "k19_diagnostic_audit": Path(fixed_reference["analysis"]) / "diagnostic-audit.json",
    }
    for label, source in reference_paths.items():
        receipts[label] = helper.freeze_receipt(
            source, receipts_dir / f"reference-{label}.json"
        )
    return receipts


def ensure_result_k(result, invocation, selected_k):
    signature = result.get("signature")
    if not isinstance(signature, dict) or signature.get("settings", {}).get("k") != selected_k:
        raise ValueError(f"Result signature has wrong k: {invocation['invocation_id']}")
    signed_invocation = signature.get("invocation")
    if not isinstance(signed_invocation, dict) or signed_invocation.get("k") != selected_k:
        raise ValueError(f"Signed invocation has wrong k: {invocation['invocation_id']}")
    command = result.get("binary_command")
    if not isinstance(command, list) or command.count("-k") != 1:
        raise ValueError(f"Result command lacks an unambiguous k: {invocation['invocation_id']}")
    command_index = command.index("-k")
    if command_index + 1 >= len(command) or command[command_index + 1] != str(selected_k):
        raise ValueError(f"Result command has wrong k: {invocation['invocation_id']}")
    if result.get("timing_status") == "complete":
        stats_path = Path(require_string(result.get("stats_path"), "result stats_path"))
        stats = context_free_stats(stats_path)
        if stats.get("kmer_length") != selected_k:
            raise ValueError(f"Preserved stats has wrong k: {invocation['invocation_id']}")


def context_free_stats(path):
    import yaml

    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Stats path is invalid: {path}")
    value = yaml.safe_load(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Stats is not an object: {path}")
    return value


def inventory_receipts(helper, directory):
    receipts = []
    for path in sorted(Path(directory).rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink in preserved output: {path}")
        if path.is_file():
            receipt = helper.file_receipt(path)
            receipt["relative_path"] = str(path.relative_to(directory))
            receipts.append(receipt)
    return receipts


def write_measurement_index(context, output_root, schedule, results):
    helper = context["helper"]
    frozen_results = output_root / "measurement-results"
    result_receipts = {}
    for invocation in schedule:
        identity = invocation["invocation_id"]
        source = output_root / "results" / f"{identity}.json"
        receipt = helper.freeze_receipt(source, frozen_results / f"{identity}.json")
        result_receipts[identity] = receipt
    index = {
        "schema_version": 1,
        "status": "complete",
        "completed_at": helper.utc_now(),
        "k": context["settings"]["k"],
        "schedule": schedule,
        "result_receipts": result_receipts,
        "raw_artifact_receipts": inventory_receipts(helper, output_root / "attempts"),
        "statuses": {identity: result["status"] for identity, result in results.items()},
        "timing_boundary": "no_blast_or_reference_database_work_before_this_index",
    }
    helper.atomic_write_json(output_root / "measurement-index.json", index)
    return index


def measure_execution(context, protocol_path, builds_path, output_root, schedule):
    if output_root.exists():
        raise ValueError(f"Output exists: {output_root}")
    output_root.mkdir(parents=True)
    receipts = freeze_execution_receipts(context, protocol_path, builds_path, output_root)
    provenance = {
        "schema_version": 1,
        "sweep_id": f"protocol-sha256:{sha256_file(protocol_path)}",
        "k": context["settings"]["k"],
        "protocol": context["protocol"],
        "frozen_receipts": receipts,
        "builds": {"schema_version": 1, "versions": context["builds"]},
        "validator": context["validator"],
        "inputs": context["input_document"],
        "schedule": schedule,
        "started_at": context["helper"].utc_now(),
        "timing_boundary": "all sweep measurements before any BLAST",
        "machine": context["runner"].get_machine_info(),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
    }
    context["helper"].atomic_write_json(output_root / "provenance.json", provenance)
    context["helper"].preflight_panels(
        context["builds"], context["panels"], output_root
    )
    results = {}
    for invocation in schedule:
        result = context["helper"].run_invocation(
            invocation,
            context["builds"][invocation["version"]],
            context["panels"][invocation["panel"]],
            context["inputs"][invocation["input"]],
            context["settings"],
            output_root,
            context["runner"],
        )
        ensure_result_k(result, invocation, context["settings"]["k"])
        results[invocation["invocation_id"]] = result
        wall_time = result.get("execution", {}).get("wall_time_s")
        print(invocation["invocation_id"], result["status"], wall_time, flush=True)
    return write_measurement_index(context, output_root, schedule, results)


def verify_receipt(receipt):
    if not isinstance(receipt, dict):
        return False
    path = Path(receipt.get("path", ""))
    return (
        path.is_file()
        and not path.is_symlink()
        and receipt.get("size_bytes") == path.stat().st_size
        and receipt.get("sha256") == sha256_file(path)
    )


def validate_master_measurement(protocol, master_path):
    master_path = Path(master_path).resolve()
    master = load_json(master_path)
    if master.get("status") != "complete":
        raise ValueError("Master measurement schedule is incomplete")
    if master.get("schedule") != schedule_for_protocol(protocol):
        raise ValueError("Master measurement schedule differs from protocol")
    receipts = master.get("receipts")
    if not isinstance(receipts, dict) or set(receipts) != {"driver", "protocol", "builds"}:
        raise ValueError("Master measurement receipts are incomplete")
    if not all(verify_receipt(receipt) for receipt in receipts.values()):
        raise ValueError("Master measurement receipt changed")
    if sha256_file(__file__) != receipts["driver"]["sha256"]:
        raise ValueError("Active classification driver differs from master measurement driver")
    expected_children = {f"k{selected_k}" for selected_k in protocol["allowed_k"]}
    children = master.get("child_measurement_indexes")
    if not isinstance(children, dict) or set(children) != expected_children:
        raise ValueError("Master child measurement set differs from protocol")
    for child_name, receipt in children.items():
        expected_path = master_path.parent / child_name / "measurement-index.json"
        if Path(receipt.get("path", "")).resolve() != expected_path.resolve() or not verify_receipt(receipt):
            raise ValueError(f"Child measurement index changed: {child_name}")
    return master


def tool_receipt(name):
    resolved = shutil.which(name)
    if resolved is None:
        raise ValueError(f"Required classification tool unavailable: {name}")
    path = Path(resolved).resolve()
    completed = subprocess.run(
        [str(path), "-version"], capture_output=True, text=True, timeout=30
    )
    if completed.returncode != 0:
        raise ValueError(f"Unable to identify classification tool: {name}")
    return {
        "name": name,
        "path": str(path),
        "sha256": sha256_file(path),
        "version_stdout": completed.stdout,
        "version_stderr": completed.stderr,
        "returncode": completed.returncode,
    }


def validate_measurement(context, output_root):
    provenance = load_json(output_root / "provenance.json")
    if provenance.get("k") != context["settings"]["k"]:
        raise ValueError("Execution provenance has wrong k")
    frozen_receipts = provenance.get("frozen_receipts", {})
    for receipt in frozen_receipts.values():
        if not verify_receipt(receipt):
            raise ValueError("Frozen execution receipt changed")
    if sha256_file(__file__) != frozen_receipts.get("driver", {}).get("sha256"):
        raise ValueError("Active classification driver differs from measured driver")
    index_path = output_root / "measurement-index.json"
    index = load_json(index_path)
    if index.get("status") != "complete" or index.get("k") != context["settings"]["k"]:
        raise ValueError("Measurement index is incomplete or has wrong k")
    expected_schedule = fixed_k_schedule(
        schedule_for_protocol(context["protocol"]), context["settings"]["k"]
    )
    master_path = provenance.get("master_schedule_path")
    if master_path is not None:
        master = validate_master_measurement(context["protocol"], master_path)
        child_name = f"k{context['settings']['k']}"
        child_receipt = master["child_measurement_indexes"][child_name]
        if (
            Path(child_receipt["path"]).resolve() != index_path.resolve()
            or child_receipt["sha256"] != sha256_file(index_path)
        ):
            raise ValueError("Linked child measurement index differs from master")
    if index.get("schedule") != expected_schedule or provenance.get("schedule") != expected_schedule:
        raise ValueError("Measurement schedule differs from protocol")
    results = {}
    for invocation in expected_schedule:
        identity = invocation["invocation_id"]
        receipt = index.get("result_receipts", {}).get(identity)
        if not verify_receipt(receipt):
            raise ValueError(f"Frozen measurement result changed: {identity}")
        current_path = output_root / "results" / f"{identity}.json"
        if sha256_file(current_path) != receipt["sha256"]:
            raise ValueError(f"Source result changed before classification: {identity}")
        result = load_json(current_path)
        ensure_result_k(result, invocation, context["settings"]["k"])
        context["helper"].validate_preserved_result(
            result,
            result["signature"],
            output_root,
            identity,
        )
        results[identity] = result
    if inventory_receipts(context["helper"], output_root / "attempts") != index.get(
        "raw_artifact_receipts"
    ):
        raise ValueError("Raw measurement artifacts changed before classification")
    return index_path, index, expected_schedule, results


def classification_identity(context, output_root, measurement_index_path, source_results):
    tools = {name: tool_receipt(name) for name in ("blastn", "makeblastdb")}
    return {
        "measurement_index": {
            "path": str(measurement_index_path),
            "sha256": sha256_file(measurement_index_path),
        },
        "source_results": {
            identity: {
                "path": str(output_root / "measurement-results" / f"{identity}.json"),
                "sha256": sha256_file(
                    output_root / "measurement-results" / f"{identity}.json"
                ),
            }
            for identity in sorted(source_results)
        },
        "build_binaries": {
            version: {
                "path": build["binary_path"],
                "sha256": sha256_file(build["binary_path"]),
                "commit": build["commit"],
                "source_archive_sha256": build["source_archive_sha256"],
            }
            for version, build in context["builds"].items()
        },
        "validator_source_tree_sha256": context["helper"].source_tree_sha256(
            Path(RELEASE_COMPARISON / "sources/candidate")
        ),
        "tools": tools,
    }


def database_receipts(helper, output_root):
    return inventory_receipts(helper, output_root / "reference_databases")


def classify_execution(context, output_root):
    index_path, unused_index, schedule, results = validate_measurement(context, output_root)
    comparison_path = output_root / "comparison.json"
    classification_path = output_root / "classification-provenance.json"
    if comparison_path.exists() or classification_path.exists():
        raise ValueError("Classification output already exists; preserve it")
    before = classification_identity(context, output_root, index_path, results)
    databases, references = context["helper"].prepare_reference_databases(
        context["panels"], context["blast"], output_root
    )
    before["reference_databases"] = database_receipts(context["helper"], output_root)
    context["helper"].atomic_write_json(
        classification_path,
        {
            "schema_version": 1,
            "status": "in_progress",
            "k": context["settings"]["k"],
            "started_at": context["helper"].utc_now(),
            "before": before,
            "references": references,
        },
    )
    for invocation in schedule:
        identity = invocation["invocation_id"]
        results[identity] = context["helper"].finalize_classification(
            invocation,
            results[identity],
            context["panels"][invocation["panel"]],
            output_root,
            context["blast"],
            databases[invocation["panel"]],
        )
        if any(
            product.get("reference_match", {}).get("status") == "failed_run"
            for gene in results[identity].get("genes", [])
            for product in gene.get("products", [])
        ):
            results[identity].update(
                {
                    "status": "failed",
                    "classification_status": "failed",
                    "failure": "per_product_blast_failed_run",
                }
            )
            context["helper"].atomic_write_json(
                output_root / "results" / f"{identity}.json", results[identity]
            )
    after = classification_identity(context, output_root, index_path, results)
    after["reference_databases"] = database_receipts(context["helper"], output_root)
    for key in ("measurement_index", "source_results", "build_binaries", "validator_source_tree_sha256", "tools"):
        if before[key] != after[key]:
            raise ValueError(f"Classification provenance changed: {key}")
    if before["reference_databases"] != after["reference_databases"]:
        raise ValueError("Reference database changed during classification")
    summary = context["helper"].summarize_results(schedule, results)
    summary["k"] = context["settings"]["k"]
    summary["protocol"] = context["protocol"]
    summary["references"] = references
    summary["limitations"] = context["protocol"]["limitations"]
    summary["release_gate"] = "Not decided; this sweep is descriptive and does not authorize release."
    context["helper"].atomic_write_json(comparison_path, summary)
    completed = {
        "schema_version": 1,
        "status": "complete",
        "k": context["settings"]["k"],
        "started_at": load_json(classification_path)["started_at"],
        "completed_at": context["helper"].utc_now(),
        "before": before,
        "after": after,
        "references": references,
        "classified_results": {
            identity: context["helper"].file_receipt(
                output_root / "results" / f"{identity}.json"
            )
            for identity in sorted(results)
        },
        "comparison": context["helper"].file_receipt(comparison_path),
    }
    context["helper"].atomic_write_json(classification_path, completed)
    return bool(summary["failures"])


def measure_sweep(arguments):
    protocol = load_json(arguments.protocol)
    allowed_k = validate_protocol(protocol)
    if arguments.output.exists():
        raise ValueError(f"Output exists: {arguments.output}")
    arguments.output.mkdir(parents=True)
    master_schedule = schedule_for_protocol(protocol)
    helper_context = prepare_context(arguments.protocol, arguments.builds, allowed_k[0])
    master_receipts = {
        "driver": helper_context["helper"].freeze_receipt(
            __file__, arguments.output / "receipts/driver.py"
        ),
        "protocol": helper_context["helper"].freeze_receipt(
            arguments.protocol, arguments.output / "receipts/protocol.json"
        ),
        "builds": helper_context["helper"].freeze_receipt(
            arguments.builds, arguments.output / "receipts/builds.json"
        ),
    }
    helper_context["helper"].atomic_write_json(
        arguments.output / "master-schedule.json",
        {
            "schema_version": 1,
            "sweep_id": f"protocol-sha256:{sha256_file(arguments.protocol)}",
            "status": "in_progress",
            "started_at": helper_context["helper"].utc_now(),
            "receipts": master_receipts,
            "schedule": master_schedule,
        },
    )
    contexts = {}
    for selected_k in allowed_k:
        context = prepare_context(arguments.protocol, arguments.builds, selected_k)
        child = arguments.output / f"k{selected_k}"
        child.mkdir()
        receipts = freeze_execution_receipts(
            context, arguments.protocol, arguments.builds, child
        )
        child_schedule = fixed_k_schedule(master_schedule, selected_k)
        context["helper"].atomic_write_json(
            child / "provenance.json",
            {
                "schema_version": 1,
                "sweep_id": f"protocol-sha256:{sha256_file(arguments.protocol)}",
                "k": selected_k,
                "protocol": protocol,
                "frozen_receipts": receipts,
                "builds": {"schema_version": 1, "versions": context["builds"]},
                "validator": context["validator"],
                "inputs": context["input_document"],
                "schedule": child_schedule,
                "master_schedule_path": str(arguments.output / "master-schedule.json"),
                "started_at": context["helper"].utc_now(),
                "timing_boundary": "all sweep measurements before any BLAST",
                "machine": context["runner"].get_machine_info(),
                "cpu_affinity": sorted(os.sched_getaffinity(0)),
            },
        )
        context["helper"].preflight_panels(
            context["builds"], context["panels"], child
        )
        contexts[selected_k] = context
    results_by_k = {selected_k: {} for selected_k in allowed_k}
    for invocation in master_schedule:
        selected_k = invocation["k"]
        context = contexts[selected_k]
        child = arguments.output / f"k{selected_k}"
        result = context["helper"].run_invocation(
            invocation,
            context["builds"][invocation["version"]],
            context["panels"][invocation["panel"]],
            context["inputs"][invocation["input"]],
            context["settings"],
            child,
            context["runner"],
        )
        ensure_result_k(result, invocation, selected_k)
        results_by_k[selected_k][invocation["invocation_id"]] = result
        print(
            invocation["invocation_id"],
            result["status"],
            result.get("execution", {}).get("wall_time_s"),
            flush=True,
        )
    child_indexes = {}
    for selected_k in allowed_k:
        child = arguments.output / f"k{selected_k}"
        index = write_measurement_index(
            contexts[selected_k],
            child,
            fixed_k_schedule(master_schedule, selected_k),
            results_by_k[selected_k],
        )
        child_indexes[f"k{selected_k}"] = contexts[selected_k]["helper"].file_receipt(
            child / "measurement-index.json"
        )
    master = load_json(arguments.output / "master-schedule.json")
    master.update(
        {
            "status": "complete",
            "completed_at": helper_context["helper"].utc_now(),
            "child_measurement_indexes": child_indexes,
        }
    )
    helper_context["helper"].atomic_write_json(
        arguments.output / "master-schedule.json", master
    )
    return any(
        result["status"] != "timed_complete"
        for results in results_by_k.values()
        for result in results.values()
    )


def measure_one(arguments):
    context = prepare_context(arguments.protocol, arguments.builds, arguments.k)
    schedule = fixed_k_schedule(
        schedule_for_protocol(context["protocol"]), arguments.k
    )
    index = measure_execution(
        context, arguments.protocol, arguments.builds, arguments.output, schedule
    )
    return any(status != "timed_complete" for status in index["statuses"].values())


def classify_one(arguments):
    provenance = load_json(arguments.execution / "provenance.json")
    selected_k = provenance.get("k")
    if type(selected_k) is not int:
        raise ValueError("Execution lacks fixed k")
    protocol_receipt = provenance.get("frozen_receipts", {}).get("protocol", {})
    builds_receipt = provenance.get("frozen_receipts", {}).get("builds", {})
    context = prepare_context(
        Path(protocol_receipt.get("path", "")),
        Path(builds_receipt.get("path", "")),
        selected_k,
    )
    return classify_execution(context, arguments.execution)


def classify_sweep(arguments):
    master_path = arguments.sweep_root / "master-schedule.json"
    receipt_protocol = load_json(master_path).get("receipts", {}).get("protocol", {})
    protocol = load_json(Path(receipt_protocol.get("path", "")))
    validate_protocol(protocol)
    master = validate_master_measurement(protocol, master_path)
    failures = False
    for child_name, receipt in master.get("child_measurement_indexes", {}).items():
        if not verify_receipt(receipt):
            raise ValueError(f"Child measurement index changed: {child_name}")
        child_arguments = argparse.Namespace(execution=arguments.sweep_root / child_name)
        failures = classify_one(child_arguments) or failures
    return failures


def main():
    parser = argparse.ArgumentParser(
        description="Frozen released-v3.1.0/current k-sweep measurement harness"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    measure = commands.add_parser("measure")
    measure.add_argument("--protocol", type=Path, required=True)
    measure.add_argument("--builds", type=Path, required=True)
    measure.add_argument("--k", type=int, required=True)
    measure.add_argument("--output", type=Path, required=True)
    sweep = commands.add_parser("measure-sweep")
    sweep.add_argument("--protocol", type=Path, required=True)
    sweep.add_argument("--builds", type=Path, required=True)
    sweep.add_argument("--output", type=Path, required=True)
    classify = commands.add_parser("classify")
    classify.add_argument("--execution", type=Path, required=True)
    classify_all = commands.add_parser("classify-sweep")
    classify_all.add_argument("--sweep-root", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.command == "measure":
        return measure_one(arguments)
    if arguments.command == "measure-sweep":
        return measure_sweep(arguments)
    if arguments.command == "classify":
        return classify_one(arguments)
    return classify_sweep(arguments)


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    raise SystemExit(main())

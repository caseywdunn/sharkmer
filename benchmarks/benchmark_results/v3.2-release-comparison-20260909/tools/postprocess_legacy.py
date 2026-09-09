#!/usr/bin/env python3

import argparse
import copy
import hashlib
import importlib.util
import json
import re
import shlex
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

LEGACY_MEDIAN = re.compile(r"(?:0|[1-9][0-9]*)(?:\.5)?")
NONNEGATIVE_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?")
INTEGER = re.compile(r"0|[1-9][0-9]*")
EXPECTED_LEGACY_ERROR = "ValueError: FASTA output lacks current kmer median header: "
REQUIRED_HEADER_FIELDS = {
    "sample",
    "gene",
    "product",
    "length",
    "kmer_count_mean",
    "kmer_count_median",
    "kmer_count_min",
    "kmer_count_max",
    "score",
}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as input_stream:
        while True:
            block = input_stream.read(1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)


def load_json(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"JSON input must be a regular non-symlink file: {path}")
    with path.open() as input_stream:
        value = json.load(input_stream)
    if not isinstance(value, dict):
        raise ValueError(f"JSON input must contain an object: {path}")
    return value


def load_frozen_driver(execution_root, suite_provenance):
    receipt = suite_provenance.get("frozen_receipts", {}).get("driver")
    if not isinstance(receipt, dict):
        raise ValueError("Source suite lacks a frozen driver receipt")
    driver_path = Path(receipt.get("path", ""))
    expected_path = execution_root / "receipts" / "driver.py"
    if driver_path.resolve() != expected_path.resolve():
        raise ValueError("Frozen driver receipt path differs from source suite")
    if driver_path.is_symlink() or not driver_path.is_file():
        raise ValueError("Frozen driver receipt is not a regular file")
    if sha256_file(driver_path) != receipt.get("sha256"):
        raise ValueError("Frozen driver receipt checksum mismatch")
    sys.dont_write_bytecode = True
    specification = importlib.util.spec_from_file_location("frozen_release_comparison_driver", driver_path)
    if specification is None or specification.loader is None:
        raise ValueError("Cannot load frozen release-comparison driver")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def parse_header(header):
    tokens = shlex.split(header)
    if not tokens:
        raise ValueError("FASTA header is empty")
    identifier = tokens[0]
    fields = {}
    for token in tokens[1:]:
        if token.count("=") != 1:
            raise ValueError(f"Malformed FASTA header token: {token}")
        key, value = token.split("=", 1)
        if not key or not value or key in fields:
            raise ValueError(f"Invalid or duplicate FASTA header field: {key}")
        fields[key] = value
    if set(fields) != REQUIRED_HEADER_FIELDS:
        raise ValueError(f"FASTA header fields differ from legacy format: {sorted(fields)}")
    return identifier, fields


def parse_nonnegative_decimal(raw_value, field):
    if NONNEGATIVE_DECIMAL.fullmatch(raw_value) is None:
        raise ValueError(f"FASTA {field} is not a nonnegative decimal: {raw_value}")
    try:
        return Decimal(raw_value)
    except InvalidOperation as error:
        raise ValueError(f"FASTA {field} is invalid: {raw_value}") from error


def parse_legacy_median(raw_value):
    if LEGACY_MEDIAN.fullmatch(raw_value) is None:
        raise ValueError(f"Legacy kmer_count_median is not an integer or .5 value: {raw_value}")
    numeric = Decimal(raw_value)
    if numeric == numeric.to_integral_value():
        return int(numeric), raw_value
    return float(numeric), raw_value


def parse_legacy_fasta(path, sample, gene):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Legacy FASTA must be a regular non-symlink file: {path}")
    records = []
    header = None
    sequence_parts = []
    with path.open() as input_stream:
        for raw_line in input_stream:
            line = raw_line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_parts)))
                header = line[1:]
                sequence_parts = []
            elif line:
                if header is None:
                    raise ValueError(f"FASTA sequence precedes header: {path}")
                sequence_parts.append(line)
    if header is not None:
        records.append((header, "".join(sequence_parts)))
    if not records:
        raise ValueError(f"Legacy FASTA is empty: {path}")
    products = []
    for record_order, (record_header, sequence) in enumerate(records):
        if not sequence or re.fullmatch(r"[ACGT]+", sequence) is None:
            raise ValueError(f"Legacy FASTA contains invalid sequence: {path}")
        identifier, fields = parse_header(record_header)
        for field in ("product", "length", "kmer_count_min", "kmer_count_max"):
            if INTEGER.fullmatch(fields[field]) is None:
                raise ValueError(f"FASTA {field} is not a nonnegative integer: {fields[field]}")
        product_index = int(fields["product"])
        length = int(fields["length"])
        minimum = int(fields["kmer_count_min"])
        maximum = int(fields["kmer_count_max"])
        median, median_raw = parse_legacy_median(fields["kmer_count_median"])
        median_decimal = Decimal(median_raw)
        mean_decimal = parse_nonnegative_decimal(fields["kmer_count_mean"], "kmer_count_mean")
        score_decimal = parse_nonnegative_decimal(fields["score"], "score")
        score = int(score_decimal) if score_decimal == score_decimal.to_integral_value() else float(score_decimal)
        if product_index != record_order:
            raise ValueError(f"Legacy product indices are not contiguous from zero: {path}")
        if fields["sample"] != sample or fields["gene"] != gene:
            raise ValueError(f"Legacy FASTA sample or gene header differs from stats: {path}")
        if identifier != f"{sample}_{gene}_{product_index}":
            raise ValueError(f"Legacy FASTA record identifier differs from metadata: {path}")
        if length != len(sequence):
            raise ValueError(f"Legacy FASTA header length differs from sequence: {path}")
        if not Decimal(minimum) <= median_decimal <= Decimal(maximum):
            raise ValueError(f"Legacy median lies outside kmer count range: {path}")
        if not Decimal(minimum) <= mean_decimal <= Decimal(maximum):
            raise ValueError(f"Legacy mean lies outside kmer count range: {path}")
        products.append(
            {
                "header": record_header,
                "product_index": product_index,
                "kmer_count_median": median,
                "kmer_count_median_raw": median_raw,
                "score": score,
                "score_raw": fields["score"],
                "length": length,
                "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "sequence": sequence,
                "record_order": record_order,
            }
        )
    return products


def should_correct_legacy_result(result):
    invocation = result.get("signature", {}).get("invocation", {})
    execution = result.get("execution", {})
    return (
        invocation.get("version") == "baseline"
        and result.get("status") == "failed"
        and result.get("failure") == "invalid_output"
        and isinstance(result.get("error"), str)
        and result["error"].startswith(EXPECTED_LEGACY_ERROR)
        and execution.get("returncode") == 0
        and execution.get("timed_out") is False
        and execution.get("orphaned_process_group_cleaned") is False
        and execution.get("measurement_complete") is True
        and execution.get("measurement_matches_exit") is True
        and execution.get("input_changed_during_invocation") is False
        and execution.get("binary_changed_during_invocation") is False
        and execution.get("panel_changed_during_invocation") is False
    )


def validate_correctable_legacy_result(result, invocation, panel, input_record, build, frozen_driver):
    attempt_dir = Path(result["attempt_dir"])
    output_dir = attempt_dir / "output"
    sample = invocation["sample_prefix"]
    stats_path = output_dir / f"{sample}.stats.yaml"
    stats = frozen_driver.parse_stats(stats_path)
    binary_command = result.get("binary_command")
    if not isinstance(binary_command, list) or not binary_command:
        raise ValueError("Correctable legacy result lacks the exact binary command")
    pcr_results = frozen_driver.validate_shared_stats(
        stats,
        invocation,
        panel,
        input_record,
        binary_command,
        "3.1.0",
    )
    if result["signature"].get("binary_sha256") != build["binary_sha256"]:
        raise ValueError("Correctable legacy result uses a different baseline binary")
    expected_files = {stats_path.name}
    genes = []
    for entry in pcr_results:
        gene = frozen_driver.normalized_gene_name(entry["gene_name"], panel)
        if entry["status"] == "success":
            output_name = f"{sample}_{entry['gene_name']}.fasta"
            expected_files.add(output_name)
            products = parse_legacy_fasta(output_dir / output_name, sample, entry["gene_name"])
            lengths = [product["length"] for product in products]
            if len(products) != entry["n_products"] or lengths != entry["product_lengths"]:
                raise ValueError(f"Legacy FASTA product lengths disagree with stats for {gene}")
            genes.append({"gene": gene, "recovered": True, "n_products": len(products), "products": products})
        else:
            genes.append(
                {
                    "gene": gene,
                    "recovered": False,
                    "n_products": 0,
                    "products": [],
                    "failure_reason": entry.get("failure_reason"),
                }
            )
    observed_files = set()
    for path in output_dir.iterdir():
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Legacy output contains a non-regular entry: {path.name}")
        observed_files.add(path.name)
    if observed_files != expected_files:
        raise ValueError(
            f"Legacy output file set mismatch: expected {sorted(expected_files)}, got {sorted(observed_files)}"
        )
    wall_time_s = result["execution"]["wall_time_s"]
    if not isinstance(wall_time_s, (int, float)) or wall_time_s <= 0:
        raise ValueError("Correctable legacy result lacks a positive GNU-time elapsed measurement")
    return {
        "stats": stats,
        "stats_path": stats_path,
        "genes": genes,
        "raw_output_files": frozen_driver.raw_output_files(output_dir),
        "metrics": frozen_driver.normalized_metrics(stats, wall_time_s),
        "final_output_bytes": sum(path.stat().st_size for path in output_dir.iterdir()),
    }


def original_lineage(result_path, original_result, correction):
    return {
        "original_result_path": str(result_path),
        "original_result_sha256": sha256_file(result_path),
        "original_status": original_result.get("status"),
        "original_timing_status": original_result.get("timing_status"),
        "original_failure": original_result.get("failure"),
        "original_error": original_result.get("error"),
        "normalization": correction,
    }


def gene_output_signature(genes):
    signature = []
    for gene in genes:
        products = []
        for product in gene.get("products", []):
            products.append(
                {
                    "header": product.get("header"),
                    "product_index": product.get("product_index"),
                    "kmer_count_median": product.get("kmer_count_median"),
                    "length": product.get("length"),
                    "sha256": product.get("sha256"),
                    "record_order": product.get("record_order"),
                }
            )
        signature.append(
            {
                "gene": gene.get("gene"),
                "recovered": gene.get("recovered"),
                "n_products": gene.get("n_products"),
                "failure_reason": gene.get("failure_reason"),
                "products": products,
            }
        )
    return signature


def revalidate_completed_result(result, invocation, panel, input_record, build, frozen_driver, runner):
    if result.get("status") != "complete" or result.get("timing_status") != "complete":
        raise ValueError(f"Result is not complete: {invocation['invocation_id']}")
    output_dir = Path(result["attempt_dir"]) / "output"
    stats_path = output_dir / f"{invocation['sample_prefix']}.stats.yaml"
    stats = frozen_driver.parse_stats(stats_path)
    binary_command = result.get("binary_command")
    if not isinstance(binary_command, list) or not binary_command:
        raise ValueError(f"Complete result lacks binary command: {invocation['invocation_id']}")
    expected_version = build["binary_version"].split()[1]
    pcr_results = frozen_driver.validate_shared_stats(
        stats,
        invocation,
        panel,
        input_record,
        binary_command,
        expected_version,
    )
    if invocation["version"] == "baseline":
        derived_genes, derived_completion = frozen_driver.validate_legacy_outputs(
            output_dir,
            invocation,
            panel,
            stats,
            pcr_results,
        )
    else:
        derived_genes, derived_completion = frozen_driver.validate_current_outputs(
            output_dir,
            invocation,
            panel,
            stats,
            pcr_results,
            binary_command,
            runner,
        )
    if gene_output_signature(derived_genes) != gene_output_signature(result.get("genes", [])):
        raise ValueError(f"Complete result genes differ from preserved raw outputs: {invocation['invocation_id']}")
    if result.get("completion") != derived_completion:
        raise ValueError(f"Complete result completion evidence differs: {invocation['invocation_id']}")
    if result.get("stats_path") != str(stats_path) or result.get("stats_sha256") != sha256_file(stats_path):
        raise ValueError(f"Complete result stats receipt differs: {invocation['invocation_id']}")
    expected_metrics = frozen_driver.normalized_metrics(stats, result["execution"]["wall_time_s"])
    if result.get("metrics") != expected_metrics:
        raise ValueError(f"Complete result metrics differ from stats and GNU time: {invocation['invocation_id']}")
    for gene in result.get("genes", []):
        for product in gene.get("products", []):
            if "sequence" in product or not isinstance(product.get("reference_match"), dict):
                raise ValueError(f"Complete result lacks finalized classification: {invocation['invocation_id']}")


def normalize_legacy_adapter_failure(
    result_path,
    original_result,
    invocation,
    panel,
    input_record,
    build,
    frozen_driver,
):
    derived = validate_correctable_legacy_result(
        original_result,
        invocation,
        panel,
        input_record,
        build,
        frozen_driver,
    )
    normalized = copy.deepcopy(original_result)
    normalized.pop("failure", None)
    normalized.pop("error", None)
    normalized.update(
        {
            "status": "normalized_legacy_pending_classification",
            "timing_status": "complete",
            "completion": {
                "completion_evidence": "exit_zero_revalidated_legacy_stats_fasta_exact_fresh_directory",
                "manifest_available": False,
            },
            "stats_path": str(derived["stats_path"]),
            "stats_sha256": sha256_file(derived["stats_path"]),
            "metrics": derived["metrics"],
            "genes": derived["genes"],
            "raw_output_files": derived["raw_output_files"],
            "final_output_bytes": derived["final_output_bytes"],
            "lineage": original_lineage(
                result_path,
                original_result,
                "legacy_v3.1_integer_or_half_median_adapter_correction",
            ),
            "correction_note": (
                "The original Sharkmer invocation exited zero with valid timing and outputs; "
                "the frozen comparison adapter rejected a valid v3.1.0 .5 median header."
            ),
        }
    )
    return normalized


def verify_source_suite(execution_root):
    suite_path = execution_root / "suite-provenance.json"
    schedule_path = execution_root / "schedule.json"
    comparison_path = execution_root / "comparison.json"
    suite = load_json(suite_path)
    schedule_document = load_json(schedule_path)
    comparison = load_json(comparison_path)
    if comparison.get("invocations") != 114:
        raise ValueError("Source comparison is incomplete")
    frozen_driver = load_frozen_driver(execution_root, suite)
    receipt_paths = {}
    for name in ("protocol", "builds", "inputs"):
        receipt = suite.get("frozen_receipts", {}).get(name)
        if not isinstance(receipt, dict):
            raise ValueError(f"Source suite lacks frozen {name} receipt")
        receipt_path = Path(receipt.get("path", ""))
        expected_path = execution_root / "receipts" / f"{name}.json"
        if (
            receipt_path.is_symlink()
            or not receipt_path.is_file()
            or receipt_path.resolve() != expected_path.resolve()
            or sha256_file(receipt_path) != receipt.get("sha256")
        ):
            raise ValueError(f"Source suite frozen {name} receipt mismatch")
        receipt_paths[name] = receipt_path
    builds = frozen_driver.verify_builds(frozen_driver.load_json(receipt_paths["builds"]))
    protocol_document = frozen_driver.load_json(receipt_paths["protocol"])
    validator_root = protocol_document["validator_root"]
    runner, blast, validator_provenance = frozen_driver.load_validator(validator_root)
    inputs = frozen_driver.verify_inputs(frozen_driver.load_json(receipt_paths["inputs"]))
    protocol = frozen_driver.verify_protocol(protocol_document, inputs, runner)
    schedule = frozen_driver.build_schedule(protocol)
    if schedule != schedule_document.get("schedule") or len(schedule) != 114:
        raise ValueError("Source schedule differs from frozen protocol")
    expected_result_names = {f"{invocation['invocation_id']}.json" for invocation in schedule}
    results_dir = execution_root / "results"
    observed_result_names = {path.name for path in results_dir.iterdir() if path.is_file() and not path.is_symlink()}
    if any(path.is_symlink() or not path.is_file() for path in results_dir.iterdir()):
        raise ValueError("Source results directory contains a non-regular entry")
    if observed_result_names != expected_result_names:
        raise ValueError("Source result set is incomplete or contains unexpected files")
    observed_reference_provenance = {
        panel_name: blast.reference_checksums(panel["data"])
        for panel_name, panel in protocol["panels"].items()
    }
    if observed_reference_provenance != suite.get("reference_checksums"):
        raise ValueError("Current validator reference checksums differ from source suite")
    return {
        "suite": suite,
        "suite_path": suite_path,
        "schedule_document": schedule_document,
        "schedule_path": schedule_path,
        "comparison": comparison,
        "comparison_path": comparison_path,
        "driver": frozen_driver,
        "runner": runner,
        "blast": blast,
        "validator_provenance": validator_provenance,
        "builds": builds,
        "inputs": inputs,
        "protocol": protocol,
        "schedule": schedule,
        "results_dir": results_dir,
    }


def classify_corrected_results(normalized_results, schedule, protocol, blast, frozen_driver, output_root):
    corrected_panels = {
        invocation["panel"]
        for invocation in schedule
        if normalized_results[invocation["invocation_id"]].get("status")
        == "normalized_legacy_pending_classification"
    }
    if not corrected_panels:
        return
    selected_panels = {name: protocol["panels"][name] for name in corrected_panels}
    databases, reference_provenance = frozen_driver.prepare_reference_databases(
        selected_panels,
        blast,
        output_root,
    )
    source_reference_provenance = {
        name: blast.reference_checksums(protocol["panels"][name]["data"])
        for name in corrected_panels
    }
    if reference_provenance != source_reference_provenance:
        raise ValueError("Reference checksums changed before corrected legacy classification")
    for invocation in schedule:
        invocation_id = invocation["invocation_id"]
        result = normalized_results[invocation_id]
        if result.get("status") != "normalized_legacy_pending_classification":
            continue
        genes = copy.deepcopy(result["genes"])
        try:
            genes = frozen_driver.evaluate_products(
                genes,
                protocol["panels"][invocation["panel"]],
                invocation["taxon"],
                databases[invocation["panel"]],
                blast,
            )
            result.update(
                {
                    "status": "complete",
                    "classification_status": "complete",
                    "classification_provenance": "current_validator_after_timing_legacy_adapter_correction",
                    "genes": genes,
                }
            )
        except Exception as error:
            result.update(
                {
                    "status": "failed",
                    "classification_status": "failed",
                    "failure": "classification_failed_after_legacy_adapter_correction",
                    "classification_error": f"{type(error).__name__}: {error}",
                }
            )


def run_postprocess(execution_root, output_root):
    execution_root = Path(execution_root).resolve()
    output_root = Path(output_root).resolve()
    if output_root == execution_root or execution_root in output_root.parents:
        raise ValueError("Postprocessed output must be outside the source execution tree")
    if output_root.exists():
        raise ValueError(f"Postprocessed output already exists: {output_root}")
    source = verify_source_suite(execution_root)
    output_root.mkdir(parents=True)
    results_output = output_root / "results"
    receipts_output = output_root / "receipts"
    results_output.mkdir()
    receipts_output.mkdir()
    frozen_driver = source["driver"]
    receipts = {
        "postprocessor": frozen_driver.freeze_receipt(__file__, receipts_output / "postprocess_legacy.py"),
        "source_suite": frozen_driver.freeze_receipt(source["suite_path"], receipts_output / "source-suite-provenance.json"),
        "source_schedule": frozen_driver.freeze_receipt(source["schedule_path"], receipts_output / "source-schedule.json"),
        "source_comparison": frozen_driver.freeze_receipt(source["comparison_path"], receipts_output / "source-comparison.json"),
    }
    normalized_results = {}
    corrections = []
    for invocation in source["schedule"]:
        invocation_id = invocation["invocation_id"]
        result_path = source["results_dir"] / f"{invocation_id}.json"
        original_result = load_json(result_path)
        if original_result.get("status") == "timed_complete":
            raise ValueError(f"Source timing/classification matrix is unfinished: {invocation_id}")
        invocation_with_sample = {
            **invocation,
            "sample_prefix": f"{invocation['panel']}_{invocation['input']}_{invocation['depth']}",
        }
        panel = source["protocol"]["panels"][invocation["panel"]]
        input_record = source["inputs"][invocation["input"]]
        build = source["builds"][invocation["version"]]
        expected_signature = frozen_driver.invocation_signature(
            invocation_with_sample,
            build,
            panel,
            input_record,
            source["protocol"]["settings"],
        )
        frozen_driver.validate_preserved_result(
            original_result,
            expected_signature,
            execution_root,
            invocation_id,
        )
        if should_correct_legacy_result(original_result):
            normalized = normalize_legacy_adapter_failure(
                result_path,
                original_result,
                invocation_with_sample,
                panel,
                input_record,
                build,
                frozen_driver,
            )
            corrections.append(invocation_id)
        else:
            if original_result.get("status") == "complete":
                revalidate_completed_result(
                    original_result,
                    invocation_with_sample,
                    panel,
                    input_record,
                    build,
                    frozen_driver,
                    source["runner"],
                )
            normalized = copy.deepcopy(original_result)
            normalized["lineage"] = original_lineage(
                result_path,
                original_result,
                "unchanged",
            )
        normalized_results[invocation_id] = normalized
    classify_corrected_results(
        normalized_results,
        source["schedule"],
        source["protocol"],
        source["blast"],
        frozen_driver,
        output_root,
    )
    for invocation_id, result in normalized_results.items():
        frozen_driver.atomic_write_json(results_output / f"{invocation_id}.json", result)
    comparison = frozen_driver.summarize_results(source["schedule"], normalized_results)
    comparison["source_comparison_sha256"] = sha256_file(source["comparison_path"])
    comparison["legacy_adapter_corrections"] = corrections
    comparison["normalization_note"] = (
        "Only exact v3.1.0 adapter failures caused by valid nonnegative integer-or-half median headers were corrected; "
        "all other original failures remain failures."
    )
    frozen_driver.atomic_write_json(output_root / "comparison.json", comparison)
    provenance = {
        "schema_version": 1,
        "source_execution_root": str(execution_root),
        "output_root": str(output_root),
        "receipts": receipts,
        "source_suite_sha256": sha256_file(source["suite_path"]),
        "source_schedule_sha256": sha256_file(source["schedule_path"]),
        "source_comparison_sha256": sha256_file(source["comparison_path"]),
        "validator": source["validator_provenance"],
        "legacy_adapter_corrections": corrections,
        "original_results_preserved": True,
        "sharkmer_invocations_rerun": False,
        "classification_scope": "parser-corrected legacy products only; all classification remains outside timing",
    }
    frozen_driver.atomic_write_json(output_root / "provenance.json", provenance)
    return 1 if comparison["failures"] else 0


def main():
    parser = argparse.ArgumentParser(description="Normalize valid v3.1.0 half-median outputs after release timing")
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        return run_postprocess(arguments.execution_root, arguments.output_root)
    except Exception as error:
        print(f"postprocess error: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

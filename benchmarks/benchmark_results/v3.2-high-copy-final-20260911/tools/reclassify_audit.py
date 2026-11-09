#!/usr/bin/env python3
import argparse
import copy
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


BLAST_DIRECTORY = Path("/tmp/sharkmer-benchmark-env/bin")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def regular_receipt(path):
    path = Path(path).resolve()
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected regular file: {path}")
    return {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}


def tree_receipts(directory):
    directory = Path(directory).resolve()
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError(f"Expected regular directory: {directory}")
    receipts = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink() or not path.is_file():
            if path.is_dir() and not path.is_symlink():
                continue
            raise ValueError(f"Database entry is not a regular file: {path}")
        receipts.append({"path": str(path.relative_to(directory)), "sha256": sha256(path), "size_bytes": path.stat().st_size})
    if not receipts:
        raise ValueError(f"Reference database is empty: {directory}")
    return receipts


def version_receipt(path):
    path = Path(path).resolve()
    completed = subprocess.run([str(path), "-version"], capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise ValueError(f"Version command failed for {path}: {completed.stderr}")
    return {**regular_receipt(path), "command": [str(path), "-version"], "stdout": completed.stdout, "stderr": completed.stderr}


def parse_fasta(path):
    records, header, sequence = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "sequence": "".join(sequence)})
            header, sequence = line[1:], []
        else:
            sequence.append(line)
    if header is not None:
        records.append({"header": header, "sequence": "".join(sequence)})
    return records


def product_sequences(result):
    invocation = result["signature"]["invocation"]
    output = Path(result["attempt_dir"]) / "output"
    values = copy.deepcopy(result["genes"])
    for gene_result in values:
        if not gene_result["products"]:
            continue
        gene = gene_result["gene"]
        path = output / f"{invocation['sample_prefix']}_{invocation['panel']}_{gene}.fasta"
        original = {product["product_index"]: product for product in gene_result["products"]}
        parsed = {}
        for index, record in enumerate(parse_fasta(path)):
            fields = dict(token.split("=", 1) for token in record["header"].split()[1:] if "=" in token)
            product_index = int(fields.get("product", index))
            parsed[product_index] = record
        if set(parsed) != set(original):
            raise ValueError(f"Saved FASTA product indices differ for {invocation['invocation_id']}/{gene}")
        for product_index, product in original.items():
            record = parsed[product_index]
            if product["header"] != record["header"] or product["length"] != len(record["sequence"]) or product["sha256"] != hashlib.sha256(record["sequence"].encode()).hexdigest():
                raise ValueError(f"Saved FASTA product differs from result receipt for {invocation['invocation_id']}/{gene}/{product_index}")
            product["sequence"] = record["sequence"]
    return values


def panels(helper, runner, provenance):
    protocol_path = Path(provenance["frozen_receipts"]["original_protocol"]["path"])
    protocol = load_json(protocol_path)
    result = {}
    for record in protocol["panels"]:
        if record["name"] not in {item["panel"] for item in provenance["schedule"]}:
            continue
        if sha256(record["path"]) != record["sha256"]:
            raise ValueError(f"Frozen panel changed: {record['name']}")
        result[record["name"]] = {**record, "data": runner.load_panel_yaml(Path(record["path"]))}
    return result


def record_difference(original, audited):
    differences = []
    for original_gene, audited_gene in zip(original, audited, strict=True):
        if original_gene["gene"] != audited_gene["gene"] or len(original_gene["products"]) != len(audited_gene["products"]):
            raise ValueError("Gene structure changed during classification audit")
        for original_product, audited_product in zip(original_gene["products"], audited_gene["products"], strict=True):
            if original_product["sha256"] != audited_product["sha256"]:
                raise ValueError("Product hash changed during classification audit")
            if original_product.get("reference_match") != audited_product.get("reference_match"):
                differences.append({"gene": original_gene["gene"], "product_index": original_product["product_index"], "sha256": original_product["sha256"], "original": original_product.get("reference_match"), "audited": audited_product.get("reference_match")})
    return differences


def main():
    parser = argparse.ArgumentParser(description="Post-timing BLAST reclassification audit")
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    execution = arguments.execution.resolve()
    if arguments.output.exists():
        raise SystemExit(f"Output already exists: {arguments.output}")
    provenance = load_json(execution / "provenance.json")
    helper_path = Path(provenance["frozen_receipts"]["helper"]["path"])
    if sha256(helper_path) != provenance["frozen_receipts"]["helper"]["sha256"]:
        raise ValueError("Frozen measurement helper receipt changed")
    validator_root = Path(provenance["validator"]["root"])
    runner_path = validator_root / "scripts/sharkmer_validate/runner.py"
    blast_path = validator_root / "scripts/sharkmer_validate/blast_references.py"
    if sha256(runner_path) != provenance["validator"]["runner_sha256"] or sha256(blast_path) != provenance["validator"]["blast_references_sha256"]:
        raise ValueError("Frozen validator receipt changed")
    blastn = (BLAST_DIRECTORY / "blastn").resolve()
    makeblastdb = (BLAST_DIRECTORY / "makeblastdb").resolve()
    before_tools = {"blastn": version_receipt(blastn), "makeblastdb": version_receipt(makeblastdb)}
    before_databases = {directory.name: tree_receipts(directory) for directory in sorted((execution / "reference_databases").iterdir())}
    arguments.output.mkdir(parents=True)
    shutil.copyfile(__file__, arguments.output / "reclassify_audit.py")
    receipt = {"script": regular_receipt(arguments.output / "reclassify_audit.py"), "helper": regular_receipt(helper_path), "runner": regular_receipt(runner_path), "blast_references": regular_receipt(blast_path), "blast_before": before_tools, "databases_before": before_databases}
    write_json(arguments.output / "receipts.json", receipt)
    previous_path = os.environ.get("PATH", "")
    os.environ["PATH"] = f"{BLAST_DIRECTORY}:{previous_path}"
    helper = load_module("frozen_measurement_helper", helper_path)
    runner, blast, _ = helper.load_validator(validator_root)
    panel_map = panels(helper, runner, provenance)
    source_hashes = {path.stem: sha256(path) for path in sorted((execution / "results").glob("*.json"))}
    derived, failures, differences = [], [], []
    derived_directory = arguments.output / "derived"
    derived_directory.mkdir()
    for invocation_id, source_hash in source_hashes.items():
        path = execution / "results" / f"{invocation_id}.json"
        try:
            result = load_json(path)
            if result.get("status") != "complete" or result.get("classification_status") != "complete":
                raise ValueError("Source result is not completely classified")
            invocation = result["signature"]["invocation"]
            classified = helper.evaluate_products(product_sequences(result), panel_map[invocation["panel"]], invocation["taxon"], execution / "reference_databases" / invocation["panel"] / "ref_db", blast)
            changed = record_difference(result["genes"], classified)
            item = {"invocation_id": invocation_id, "source_result_sha256": source_hash, "panel": invocation["panel"], "taxon": invocation["taxon"], "products": [{"gene": gene["gene"], "products": [{key: product[key] for key in ("product_index", "length", "sha256", "reference_match")} for product in gene["products"]]} for gene in classified], "differences": changed}
            write_json(derived_directory / f"{invocation_id}.json", item)
            derived.append(item)
            differences.extend({"invocation_id": invocation_id, **difference} for difference in changed)
        except Exception as error:
            failures.append({"invocation_id": invocation_id, "source_result_sha256": source_hash, "error": f"{type(error).__name__}: {error}"})
    after_tools = {"blastn": version_receipt(blastn), "makeblastdb": version_receipt(makeblastdb)}
    after_databases = {directory.name: tree_receipts(directory) for directory in sorted((execution / "reference_databases").iterdir())}
    source_after = {path.stem: sha256(path) for path in sorted((execution / "results").glob("*.json"))}
    unchanged = source_hashes == source_after
    write_json(arguments.output / "summary.json", {"schema_version": 1, "label": "post-timing reclassification audit; it does not retroactively pin the original classification executables", "execution": str(execution), "source_result_hashes_before": source_hashes, "source_result_hashes_after": source_after, "source_results_unchanged": unchanged, "blast_before": before_tools, "blast_after": after_tools, "databases_before": before_databases, "databases_after": after_databases, "derived_results": len(derived), "failures": failures, "differences": differences, "passed": unchanged and before_tools == after_tools and before_databases == after_databases and not failures and not differences})
    raise SystemExit(0 if unchanged and before_tools == after_tools and before_databases == after_databases and not failures and not differences else 1)


if __name__ == "__main__":
    main()

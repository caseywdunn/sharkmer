#!/usr/bin/env python3

import argparse
import copy
import hashlib
import importlib
import json
import os
import re
import resource
import shlex
import signal
import statistics
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import yaml


HEX_SHA256 = re.compile(r"[0-9a-f]{64}")
HEX_COMMIT = re.compile(r"[0-9a-f]{40}")
EXPECTED_COMMITS = {
    "baseline": "5a664680c91ad59f32b8c2a847b8fef37f34a0ae",
    "candidate": "ba64f573048b6c19a528278028810af5bd475b81",
}


def sha256_file(path, byte_limit=None):
    digest = hashlib.sha256()
    remaining = byte_limit
    with Path(path).open("rb") as input_stream:
        while remaining is None or remaining > 0:
            block_size = 1024 * 1024 if remaining is None else min(1024 * 1024, remaining)
            block = input_stream.read(block_size)
            if not block:
                break
            digest.update(block)
            if remaining is not None:
                remaining -= len(block)
    if remaining not in (None, 0):
        raise ValueError(f"File ended before requested byte limit: {path}")
    return digest.hexdigest()


def source_tree_sha256(directory):
    directory = Path(directory)
    digest = hashlib.sha256()
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            value = f"symlink:{os.readlink(path)}"
        elif path.is_file():
            value = sha256_file(path)
        else:
            continue
        digest.update(str(path.relative_to(directory)).encode())
        digest.update(b"\0")
        digest.update(value.encode())
        digest.update(b"\n")
    return digest.hexdigest()


def load_json(path):
    with Path(path).open() as input_stream:
        value = json.load(input_stream)
    if not isinstance(value, dict):
        raise ValueError(f"JSON document must be an object: {path}")
    return value


def atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("x") as output_stream:
        json.dump(value, output_stream, indent=2, sort_keys=True)
        output_stream.write("\n")
        output_stream.flush()
        os.fsync(output_stream.fileno())
    os.replace(temporary, path)


def freeze_receipt(source_path, destination_path):
    source_path = verify_regular_file(source_path, "receipt source")
    destination_path = Path(destination_path)
    content = source_path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    if destination_path.exists() or destination_path.is_symlink():
        destination_path = verify_regular_file(destination_path, "frozen receipt")
        if destination_path.read_bytes() != content:
            raise ValueError(f"Frozen receipt differs from current input: {destination_path.name}")
    else:
        with destination_path.open("xb") as output_stream:
            output_stream.write(content)
            output_stream.flush()
            os.fsync(output_stream.fileno())
        destination_path.chmod(0o444)
    return {"path": str(destination_path.resolve()), "sha256": digest, "size_bytes": len(content)}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def require_mapping(value, context):
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def require_list(value, context):
    if not isinstance(value, list):
        raise ValueError(f"{context} must be a list")
    return value


def require_string(value, context):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a nonblank string")
    return value


def verify_sha256(value, context):
    if not isinstance(value, str) or HEX_SHA256.fullmatch(value) is None:
        raise ValueError(f"{context} must be a lowercase SHA-256")


def verify_regular_file(path, context):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{context} must be a regular non-symlink file: {path}")
    return path


def verify_builds(builds_document):
    if builds_document.get("schema_version") != 1:
        raise ValueError("builds.json schema_version must be 1")
    versions = require_mapping(builds_document.get("versions"), "builds.json versions")
    if set(versions) != {"baseline", "candidate"}:
        raise ValueError("builds.json versions must contain exactly baseline and candidate")
    verified = {}
    for version_name, raw_build in versions.items():
        build = require_mapping(raw_build, f"build {version_name}")
        commit = require_string(build.get("commit"), f"build {version_name} commit")
        if HEX_COMMIT.fullmatch(commit) is None:
            raise ValueError(f"build {version_name} commit must be a full Git commit")
        if commit != EXPECTED_COMMITS[version_name]:
            raise ValueError(f"build {version_name} commit differs from the frozen comparison")
        binary_path = verify_regular_file(build.get("binary_path"), f"build {version_name} binary")
        verify_sha256(build.get("binary_sha256"), f"build {version_name} binary_sha256")
        if sha256_file(binary_path) != build["binary_sha256"]:
            raise ValueError(f"build {version_name} binary checksum mismatch")
        source_archive = verify_regular_file(build.get("source_archive"), f"build {version_name} source archive")
        verify_sha256(build.get("source_archive_sha256"), f"build {version_name} source archive checksum")
        if sha256_file(source_archive) != build["source_archive_sha256"]:
            raise ValueError(f"build {version_name} source archive checksum mismatch")
        source_export = Path(require_string(build.get("source_export"), f"build {version_name} source export"))
        if source_export.is_symlink() or not source_export.is_dir():
            raise ValueError(f"build {version_name} source export must be a directory")
        if build.get("source_export_pristine_after_build") is not True:
            raise ValueError(f"build {version_name} source export was not attested pristine")
        before_hash = build.get("source_tree_sha256_before_build")
        after_hash = build.get("source_tree_sha256_after_build")
        verify_sha256(before_hash, f"build {version_name} pre-build source checksum")
        verify_sha256(after_hash, f"build {version_name} post-build source checksum")
        if before_hash != after_hash or source_tree_sha256(source_export) != after_hash:
            raise ValueError(f"build {version_name} source export is not pristine")
        cargo_lock = verify_regular_file(source_export / "Cargo.lock", f"build {version_name} Cargo.lock")
        verify_sha256(build.get("cargo_lock_sha256"), f"build {version_name} Cargo.lock checksum")
        if sha256_file(cargo_lock) != build["cargo_lock_sha256"]:
            raise ValueError(f"build {version_name} Cargo.lock checksum mismatch")
        cargo_toml = verify_regular_file(source_export / "Cargo.toml", f"build {version_name} Cargo.toml")
        build_command = require_list(build.get("build_command"), f"build {version_name} command")
        if not all(flag in build_command for flag in ("--locked", "--offline", "--release")):
            raise ValueError(f"build {version_name} was not attested as locked offline release build")
        features = set(require_list(build.get("cargo_artifact_features"), f"build {version_name} features"))
        if features != {"ahashmap", "default"}:
            raise ValueError(f"build {version_name} does not use default ahashmap features")
        if build.get("source_binding") != "built_from_verified_clean_git_archive_with_locked_dependencies":
            raise ValueError(f"build {version_name} source binding is not explicit")
        if Path(build.get("build_cwd", "")).resolve() != source_export.resolve():
            raise ValueError(f"build {version_name} working directory differs from source export")
        build_environment = require_mapping(build.get("build_environment"), f"build {version_name} environment")
        for variable in ("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC_WRAPPER", "RUSTC_WORKSPACE_WRAPPER"):
            if build_environment.get(variable) is not None:
                raise ValueError(f"build {version_name} uses unsupported {variable}")
        observed_version = subprocess.run(
            [str(binary_path), "--version"], capture_output=True, text=True, check=True
        ).stdout.strip()
        if observed_version != build.get("binary_version"):
            raise ValueError(f"build {version_name} binary version changed")
        if version_name == "baseline":
            if build.get("release_tag") != "v3.1.0" or build.get("release_tag_commit") != commit:
                raise ValueError("baseline build is not attested to the v3.1.0 tag commit")
        verified[version_name] = {
            **build,
            "binary_path": str(binary_path.resolve()),
            "source_archive": str(source_archive.resolve()),
            "source_export": str(source_export.resolve()),
            "cargo_toml_sha256": sha256_file(cargo_toml),
        }
    return verified


def scan_fastq(path, prefixes):
    requested_depths = sorted(prefixes)
    prefix_results = {}
    digest = hashlib.sha256()
    records = 0
    bases = 0
    bytes_read = 0
    with Path(path).open("rb") as input_stream:
        magic = input_stream.read(2)
        input_stream.seek(0)
        if magic == b"\x1f\x8b":
            raise ValueError(f"Timed input must be uncompressed FASTQ: {path}")
        while True:
            header = input_stream.readline()
            if not header:
                break
            sequence = input_stream.readline()
            plus = input_stream.readline()
            quality = input_stream.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Incomplete FASTQ record {records + 1}: {path}")
            sequence_value = sequence.rstrip(b"\r\n")
            quality_value = quality.rstrip(b"\r\n")
            if not header.startswith(b"@") or not plus.startswith(b"+"):
                raise ValueError(f"Invalid FASTQ structure at record {records + 1}: {path}")
            if len(sequence_value) != len(quality_value):
                raise ValueError(f"FASTQ sequence/quality length mismatch at record {records + 1}: {path}")
            for line in (header, sequence, plus, quality):
                digest.update(line)
                bytes_read += len(line)
            records += 1
            bases += len(sequence_value)
            if records in requested_depths:
                prefix_results[records] = {
                    "records": records,
                    "bases": bases,
                    "size_bytes": bytes_read,
                    "sha256": digest.copy().hexdigest(),
                }
    for requested_depth in requested_depths:
        if requested_depth > records:
            prefix_results[requested_depth] = {
                "records": records,
                "bases": bases,
                "size_bytes": bytes_read,
                "sha256": digest.hexdigest(),
            }
    return {
        "total_records": records,
        "total_bases": bases,
        "size_bytes": bytes_read,
        "sha256": digest.hexdigest(),
        "prefixes": prefix_results,
    }


def verify_inputs(inputs_document):
    if inputs_document.get("schema_version") != 1:
        raise ValueError("inputs.json schema_version must be 1")
    preparation = require_mapping(inputs_document.get("preparation"), "inputs preparation")
    if preparation.get("record_order") != "ENA URL order, sequential, unpaired":
        raise ValueError("inputs must attest ENA URL order with sequential unpaired records")
    verified = {}
    for raw_input in require_list(inputs_document.get("inputs"), "inputs"):
        input_record = require_mapping(raw_input, "input record")
        input_id = require_string(input_record.get("id"), "input id")
        if input_id in verified:
            raise ValueError(f"Duplicate input id: {input_id}")
        accession = require_string(input_record.get("accession"), f"input {input_id} accession")
        if accession != input_id:
            raise ValueError(f"input {input_id} id and accession differ")
        input_path = verify_regular_file(input_record.get("path"), f"input {input_id} path")
        verify_sha256(input_record.get("sha256"), f"input {input_id} checksum")
        prefixes_raw = require_mapping(input_record.get("prefixes"), f"input {input_id} prefixes")
        prefixes = {}
        for depth_text, prefix in prefixes_raw.items():
            try:
                requested_depth = int(depth_text)
            except (TypeError, ValueError) as error:
                raise ValueError(f"input {input_id} prefix key must be an integer string") from error
            if requested_depth <= 0:
                raise ValueError(f"input {input_id} prefix depth must be positive")
            prefix = require_mapping(prefix, f"input {input_id} prefix {requested_depth}")
            verify_sha256(prefix.get("sha256"), f"input {input_id} prefix {requested_depth} checksum")
            prefixes[requested_depth] = prefix
        observed = scan_fastq(input_path, prefixes)
        for field in ("sha256", "size_bytes", "total_records", "total_bases"):
            if input_record.get(field) != observed[field]:
                raise ValueError(f"input {input_id} {field} mismatch")
        for requested_depth, expected in prefixes.items():
            if expected != observed["prefixes"][requested_depth]:
                raise ValueError(f"input {input_id} prefix {requested_depth} mismatch")
        source_urls = require_list(input_record.get("source_urls_ordered"), f"input {input_id} source URLs")
        if not source_urls or any(not isinstance(url, str) or not url for url in source_urls):
            raise ValueError(f"input {input_id} source URLs are missing")
        verified[input_id] = {
            **input_record,
            "path": str(input_path.resolve()),
            "prefixes": prefixes,
        }
    return verified


def load_validator(validator_root):
    validator_root = Path(validator_root).resolve()
    scripts_path = validator_root / "scripts"
    verify_regular_file(scripts_path / "sharkmer_validate" / "runner.py", "validator runner")
    verify_regular_file(scripts_path / "sharkmer_validate" / "blast_references.py", "validator BLAST module")
    scripts_text = str(scripts_path)
    if scripts_text not in sys.path:
        sys.path.insert(0, scripts_text)
    sys.dont_write_bytecode = True
    runner = importlib.import_module("sharkmer_validate.runner")
    blast = importlib.import_module("sharkmer_validate.blast_references")
    return runner, blast, {
        "root": str(validator_root),
        "runner_sha256": sha256_file(scripts_path / "sharkmer_validate" / "runner.py"),
        "blast_references_sha256": sha256_file(scripts_path / "sharkmer_validate" / "blast_references.py"),
    }


def validate_cpu_list(cpu_list, threads):
    if not isinstance(cpu_list, list) or len(cpu_list) != threads:
        raise ValueError("cpu_list must contain exactly one CPU per Sharkmer thread")
    if any(type(cpu_number) is not int or cpu_number < 0 for cpu_number in cpu_list):
        raise ValueError("cpu_list values must be nonnegative integers")
    if len(set(cpu_list)) != len(cpu_list):
        raise ValueError("cpu_list contains duplicates")
    if hasattr(os, "sched_getaffinity"):
        allowed = os.sched_getaffinity(0)
        if not set(cpu_list).issubset(allowed):
            raise ValueError(f"cpu_list is outside current process affinity: {sorted(allowed)}")
    physical_cores = []
    for cpu_number in cpu_list:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu_number}/topology")
        package_path = topology / "physical_package_id"
        core_path = topology / "core_id"
        if not package_path.is_file() or not core_path.is_file():
            raise ValueError(f"Cannot verify physical-core identity for CPU {cpu_number}")
        physical_cores.append((int(package_path.read_text()), int(core_path.read_text())))
    if len(set(physical_cores)) != len(physical_cores):
        raise ValueError("cpu_list contains sibling threads from the same physical core")
    return physical_cores


def verify_protocol(protocol_document, verified_inputs, runner):
    if protocol_document.get("schema_version") != 1:
        raise ValueError("protocol.json schema_version must be 1")
    comparison_id = require_string(protocol_document.get("comparison_id"), "comparison_id")
    if re.fullmatch(r"[A-Za-z0-9_.-]+", comparison_id) is None:
        raise ValueError("comparison_id contains unsafe characters")
    settings = require_mapping(protocol_document.get("settings"), "protocol settings")
    exact_settings = {
        "k": 19,
        "threads": 2,
        "chunks": 0,
        "read_threading": False,
        "paired": False,
        "timeout_seconds": 1800,
        "address_space_limit_bytes": 40 * 1024 * 1024 * 1024,
        "prewarm_each_invocation": True,
        "input_format": "uncompressed_four_line_fastq",
    }
    for key, expected in exact_settings.items():
        if settings.get(key) != expected:
            raise ValueError(f"protocol setting {key} must equal {expected!r}")
    physical_cores = validate_cpu_list(settings.get("cpu_list"), settings["threads"])
    panels = {}
    for raw_panel in require_list(protocol_document.get("panels"), "protocol panels"):
        panel = require_mapping(raw_panel, "panel")
        panel_name = require_string(panel.get("name"), "panel name")
        if re.fullmatch(r"[A-Za-z0-9_.-]+", panel_name) is None:
            raise ValueError(f"Panel name contains unsafe characters: {panel_name}")
        if panel_name in panels:
            raise ValueError(f"Duplicate panel: {panel_name}")
        panel_path = verify_regular_file(panel.get("path"), f"panel {panel_name} path")
        verify_sha256(panel.get("sha256"), f"panel {panel_name} checksum")
        if sha256_file(panel_path) != panel["sha256"]:
            raise ValueError(f"panel {panel_name} checksum mismatch")
        panel_data = runner.load_panel_yaml(panel_path)
        if panel_data.get("name") != panel_name:
            raise ValueError(f"panel {panel_name} internal name differs")
        native_genes = runner.panel_gene_names(panel_data)
        if not native_genes or len(native_genes) != len(set(native_genes)):
            raise ValueError(f"panel {panel_name} has missing or duplicate gene names")
        output_prefix = panel_data.get("gene_prefix") or panel_data["name"]
        panels[panel_name] = {
            **panel,
            "path": str(panel_path.resolve()),
            "data": panel_data,
            "native_genes": native_genes,
            "output_prefix": output_prefix,
            "stats_genes": {f"{output_prefix}_{gene}" for gene in native_genes},
        }
    samples = []
    seen_samples = set()
    for raw_sample in require_list(protocol_document.get("samples"), "protocol samples"):
        sample = require_mapping(raw_sample, "sample")
        panel_name = require_string(sample.get("panel"), "sample panel")
        input_id = require_string(sample.get("input"), "sample input")
        taxon = require_string(sample.get("taxon"), "sample taxon")
        if panel_name not in panels or input_id not in verified_inputs:
            raise ValueError(f"sample {panel_name}/{input_id} references unknown panel or input")
        if re.fullmatch(r"[A-Za-z0-9_.-]+", input_id) is None:
            raise ValueError(f"Input id contains unsafe characters: {input_id}")
        identity = (panel_name, input_id)
        if identity in seen_samples:
            raise ValueError(f"Duplicate sample: {panel_name}/{input_id}")
        seen_samples.add(identity)
        depths = require_list(sample.get("depths"), f"sample {panel_name}/{input_id} depths")
        if not depths or any(type(depth) is not int or depth <= 0 for depth in depths):
            raise ValueError(f"sample {panel_name}/{input_id} depths must be positive integers")
        if len(depths) != len(set(depths)):
            raise ValueError(f"sample {panel_name}/{input_id} has duplicate depths")
        for depth in depths:
            if depth not in verified_inputs[input_id]["prefixes"]:
                raise ValueError(f"input {input_id} has no verified prefix for depth {depth}")
        samples.append({**sample, "depths": sorted(depths)})
    if len(samples) != 13:
        raise ValueError("protocol must contain exactly 13 panel/sample combinations")
    primary_depth = 1_000_000
    expected_deep_depths = [primary_depth, 2_000_000, 4_000_000, 8_000_000]
    if any(primary_depth not in sample["depths"] for sample in samples):
        raise ValueError("Every protocol sample must include the primary depth")
    deep_samples = [sample for sample in samples if sample["depths"] != [primary_depth]]
    if len(deep_samples) != 6 or any(sample["depths"] != expected_deep_depths for sample in deep_samples):
        raise ValueError("Protocol must define exactly six samples at 1M, 2M, 4M, and 8M")
    schedule_settings = require_mapping(protocol_document.get("schedule"), "protocol schedule")
    if schedule_settings.get("primary_depth") != 1_000_000:
        raise ValueError("primary_depth must be 1000000")
    if schedule_settings.get("primary_pairs") != 3 or schedule_settings.get("deeper_pairs") != 1:
        raise ValueError("protocol requires three primary pairs and one deeper pair")
    if schedule_settings.get("depth_order") != "ascending; all primary cells before deeper cells":
        raise ValueError("protocol depth_order must put all primary cells before deeper cells")
    order_seed_raw = schedule_settings.get("order_seed")
    if not isinstance(order_seed_raw, (str, int)) or isinstance(order_seed_raw, bool):
        raise ValueError("schedule order_seed must be a string or integer")
    order_seed = str(order_seed_raw)
    return {
        "comparison_id": comparison_id,
        "settings": settings,
        "physical_cores": physical_cores,
        "panels": panels,
        "samples": samples,
        "schedule": schedule_settings,
        "order_seed": order_seed,
    }


def build_schedule(protocol):
    schedule = []
    cell_index = 0
    primary_depth = protocol["schedule"]["primary_depth"]
    cells = []
    for sample in protocol["samples"]:
        if primary_depth in sample["depths"]:
            cells.append((sample, primary_depth))
    deeper_depths = sorted(
        {depth for sample in protocol["samples"] for depth in sample["depths"] if depth != primary_depth}
    )
    for depth in deeper_depths:
        for sample in protocol["samples"]:
            if depth in sample["depths"]:
                cells.append((sample, depth))
    for sample, depth in cells:
        pair_count = protocol["schedule"]["primary_pairs"] if depth == primary_depth else protocol["schedule"]["deeper_pairs"]
        identity = f"{sample['panel']}/{sample['input']}/{depth}"
        seed_digest = hashlib.sha256(f"{protocol['order_seed']}:{identity}".encode()).digest()
        initial_order = ["baseline", "candidate"] if seed_digest[0] % 2 == 0 else ["candidate", "baseline"]
        for pair_index in range(pair_count):
            version_order = initial_order if pair_index % 2 == 0 else list(reversed(initial_order))
            for order_position, version_name in enumerate(version_order):
                safe_identity = identity.replace("/", "_")
                invocation_id = f"{cell_index:03d}_{safe_identity}_pair{pair_index + 1}_{version_name}"
                schedule.append(
                    {
                        "invocation_id": invocation_id,
                        "cell_index": cell_index,
                        "cell": identity,
                        "panel": sample["panel"],
                        "input": sample["input"],
                        "taxon": sample["taxon"],
                        "depth": depth,
                        "pair_index": pair_index + 1,
                        "order_position": order_position + 1,
                        "version": version_name,
                    }
                )
        cell_index += 1
    return schedule


def parse_header_fields(header):
    fields = {}
    for token in shlex.split(header):
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def parse_fasta_file(path):
    path = verify_regular_file(path, "FASTA output")
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
        raise ValueError(f"FASTA output is empty: {path}")
    parsed = []
    seen_indices = set()
    for record_number, (record_header, sequence) in enumerate(records):
        if not sequence or re.fullmatch(r"[ACGT]+", sequence) is None:
            raise ValueError(f"FASTA output contains invalid sequence: {path}")
        fields = parse_header_fields(record_header)
        product_text = fields.get("product")
        median_text = fields.get("kmer_count_median")
        if product_text is None or not product_text.isdigit():
            raise ValueError(f"FASTA output lacks current product index header: {path}")
        if median_text is None or not median_text.isdigit():
            raise ValueError(f"FASTA output lacks current kmer median header: {path}")
        product_index = int(product_text)
        if product_index in seen_indices:
            raise ValueError(f"FASTA output has duplicate product index: {path}")
        seen_indices.add(product_index)
        parsed.append(
            {
                "header": record_header,
                "product_index": product_index,
                "kmer_count_median": int(median_text),
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "sequence": sequence,
                "record_order": record_number,
            }
        )
    return parsed


def parse_stats(path):
    path = verify_regular_file(path, "stats output")
    try:
        value = yaml.safe_load(path.read_text())
    except yaml.YAMLError as error:
        raise ValueError(f"Stats YAML is malformed: {error}") from error
    return require_mapping(value, "stats YAML")


def validate_shared_stats(stats, invocation, panel, input_record, binary_command, expected_version):
    sample_prefix = invocation["sample_prefix"]
    expected_prefix = input_record["prefixes"][invocation["depth"]]
    if stats.get("sample") != sample_prefix:
        raise ValueError("Stats sample differs from invocation")
    if stats.get("kmer_length") != 19 or stats.get("chunks") != 0:
        raise ValueError("Stats k/chunks differ from protocol")
    if stats.get("command") != " ".join(binary_command):
        raise ValueError("Stats command differs from exact invocation")
    if stats.get("sharkmer_version") != expected_version:
        raise ValueError("Stats version differs from selected binary")
    if stats.get("n_reads_read") != expected_prefix["records"]:
        raise ValueError("Stats record count differs from frozen input prefix")
    if stats.get("n_bases_read") != expected_prefix["bases"]:
        raise ValueError("Stats base count differs from frozen input prefix")
    for field in ("n_subreads_ingested", "n_bases_ingested", "n_kmers", "peak_memory_bytes"):
        if type(stats.get(field)) is not int or stats[field] < 0:
            raise ValueError(f"Stats field {field} is missing or invalid")
    pcr_results = require_list(stats.get("pcr_results"), "stats pcr_results")
    observed_genes = {entry.get("gene_name") for entry in pcr_results if isinstance(entry, dict)}
    if len(pcr_results) != len(panel["stats_genes"]) or observed_genes != panel["stats_genes"]:
        raise ValueError("Stats PCR gene set differs from frozen panel")
    for entry in pcr_results:
        require_mapping(entry, "PCR result")
        status = entry.get("status")
        if status not in {"success", "fail"}:
            raise ValueError(f"Invalid PCR status for {entry.get('gene_name')}")
        if type(entry.get("n_products")) is not int or entry["n_products"] < 0:
            raise ValueError(f"Invalid product count for {entry.get('gene_name')}")
        if status == "success":
            lengths = require_list(entry.get("product_lengths"), "PCR product lengths")
            if (
                entry["n_products"] == 0
                or len(lengths) != entry["n_products"]
                or any(type(length) is not int or length <= 0 for length in lengths)
            ):
                raise ValueError(f"Successful PCR result has invalid lengths for {entry.get('gene_name')}")
        elif entry["n_products"] != 0:
            raise ValueError(f"Failed PCR result reports products for {entry.get('gene_name')}")
    return pcr_results


def normalized_gene_name(stats_gene, panel):
    prefix = f"{panel['output_prefix']}_"
    if not stats_gene.startswith(prefix):
        raise ValueError(f"Stats gene does not use frozen panel prefix: {stats_gene}")
    gene = stats_gene[len(prefix):]
    if not gene:
        raise ValueError("Stats gene name is empty after removing panel prefix")
    return gene


def validate_legacy_outputs(output_dir, invocation, panel, stats, pcr_results):
    if any(field in stats for field in ("run_id", "run_status", "output_manifest", "input_source", "stage_timings")):
        raise ValueError("Legacy output unexpectedly contains current-only provenance fields")
    if not str(stats.get("sharkmer_version", "")).startswith("3.1.0"):
        raise ValueError("Legacy adapter only accepts Sharkmer 3.1.0 stats")
    sample_prefix = invocation["sample_prefix"]
    expected_files = {f"{sample_prefix}.stats.yaml"}
    gene_results = []
    for entry in pcr_results:
        gene = normalized_gene_name(entry["gene_name"], panel)
        if entry["status"] == "success":
            output_name = f"{sample_prefix}_{entry['gene_name']}.fasta"
            expected_files.add(output_name)
            products = parse_fasta_file(output_dir / output_name)
            if len(products) != entry["n_products"] or [product["length"] for product in products] != entry["product_lengths"]:
                raise ValueError(f"Legacy FASTA disagrees with stats for {gene}")
            gene_results.append({"gene": gene, "recovered": True, "n_products": len(products), "products": products})
        else:
            gene_results.append(
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
            raise ValueError(f"Legacy output contains non-regular entry: {path.name}")
        observed_files.add(path.name)
    if observed_files != expected_files:
        raise ValueError(f"Legacy output file set mismatch: expected {sorted(expected_files)}, got {sorted(observed_files)}")
    return gene_results, {
        "completion_evidence": "exit_zero_validated_stats_fasta_exact_fresh_directory",
        "manifest_available": False,
    }


def validate_current_outputs(output_dir, invocation, panel, stats, pcr_results, binary_command, runner):
    sample_prefix = invocation["sample_prefix"]
    runner._validate_stats_manifest(
        stats,
        sample_prefix,
        19,
        panel["stats_genes"],
        "end-to-end",
        expected_command=binary_command,
    )
    transaction = runner._validate_output_manifest(stats, sample_prefix, output_dir)
    output_files = [entry["output_file"] for entry in pcr_results if entry["status"] == "success"]
    parsed_groups = runner.parse_fasta_products(sample_prefix, output_dir, output_files)
    if runner._validate_output_manifest(stats, sample_prefix, output_dir) != transaction:
        raise ValueError("Current output transaction changed during parsing")
    parsed_by_file = {group["output_file"]: group for group in parsed_groups}
    gene_results = []
    for entry in pcr_results:
        gene = normalized_gene_name(entry["gene_name"], panel)
        if entry["status"] == "success":
            group = parsed_by_file.get(entry["output_file"])
            if group is None:
                raise ValueError(f"Current manifest lacks successful FASTA for {gene}")
            products = []
            for product in group["products"]:
                products.append(
                    {
                        **product,
                        "sha256": hashlib.sha256(product["sequence"].encode()).hexdigest(),
                        "record_order": len(products),
                    }
                )
            if len(products) != entry["n_products"] or [product["length"] for product in products] != entry["product_lengths"]:
                raise ValueError(f"Current FASTA disagrees with stats for {gene}")
            gene_results.append({"gene": gene, "recovered": True, "n_products": len(products), "products": products})
        else:
            gene_results.append(
                {
                    "gene": gene,
                    "recovered": False,
                    "n_products": 0,
                    "products": [],
                    "failure_reason": entry.get("failure_reason"),
                }
            )
    return gene_results, {
        "completion_evidence": "exit_zero_validated_current_transaction",
        "manifest_available": True,
        "run_id": stats.get("run_id"),
        "manifest": stats.get("output_manifest"),
    }


def evaluate_products(gene_results, panel, taxon, db_path, blast):
    run_results = [{"success": True, "genes": gene_results}]
    reference_genes = {reference["gene_name"] for reference in blast.extract_references(panel["data"])}
    blast.blast_all_products(
        run_results,
        db_path,
        taxon,
        skip_blast=False,
        reference_genes=reference_genes,
    )
    for gene_result in gene_results:
        for product in gene_result.get("products", []):
            product.pop("sequence", None)
    return gene_results


def parse_gnu_time(path):
    result = {"raw_path": str(path), "peak_rss_bytes": None, "exit_status": None, "wall_time_s": None}
    if not Path(path).is_file():
        return result
    for line in Path(path).read_text(errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("Maximum resident set size (kbytes):"):
            result["peak_rss_bytes"] = int(stripped.rsplit(":", 1)[1].strip()) * 1024
        elif stripped.startswith("Exit status:"):
            result["exit_status"] = int(stripped.rsplit(":", 1)[1].strip())
        elif stripped.startswith("Elapsed (wall clock) time"):
            elapsed_text = stripped.split("):", 1)[1].strip()
            components = elapsed_text.split(":")
            if len(components) == 2:
                result["wall_time_s"] = int(components[0]) * 60 + float(components[1])
            elif len(components) == 3:
                result["wall_time_s"] = int(components[0]) * 3600 + int(components[1]) * 60 + float(components[2])
    return result


def file_receipt(path):
    path = verify_regular_file(path, "artifact")
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def inventory_untrusted_directory(directory):
    inventory = []
    for path in sorted(Path(directory).iterdir()):
        if path.is_symlink():
            inventory.append({"path": path.name, "kind": "symlink", "target": os.readlink(path)})
        elif path.is_file():
            inventory.append({"path": path.name, "kind": "file", "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
        elif path.is_dir():
            inventory.append({"path": path.name, "kind": "directory"})
        else:
            inventory.append({"path": path.name, "kind": "other"})
    return inventory


def process_limit(address_space_limit_bytes):
    resource.setrlimit(resource.RLIMIT_AS, (address_space_limit_bytes, address_space_limit_bytes))


def process_group_exists(process_group_id):
    try:
        os.killpg(process_group_id, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def stop_process_group(process_group_id):
    if not process_group_exists(process_group_id):
        return False
    try:
        os.killpg(process_group_id, signal.SIGTERM)
    except ProcessLookupError:
        return True
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        if not process_group_exists(process_group_id):
            return True
        time.sleep(0.05)
    if process_group_exists(process_group_id):
        try:
            os.killpg(process_group_id, signal.SIGKILL)
        except ProcessLookupError:
            pass
    return True


def run_timed(command, stdout_path, stderr_path, time_path, timeout_seconds, address_space_limit_bytes, cpu_list, environment):
    timed_command = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_path),
        "/usr/bin/taskset",
        "--cpu-list",
        ",".join(str(cpu_number) for cpu_number in cpu_list),
        *command,
    ]
    started = time.monotonic()
    timed_out = False
    orphaned_group = False
    with Path(stdout_path).open("wb") as stdout_stream, Path(stderr_path).open("wb") as stderr_stream:
        process = subprocess.Popen(
            timed_command,
            stdout=stdout_stream,
            stderr=stderr_stream,
            env=environment,
            start_new_session=True,
            preexec_fn=lambda: process_limit(address_space_limit_bytes),
        )
        try:
            returncode = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            stop_process_group(process.pid)
            process.wait()
            returncode = process.returncode
    if not timed_out and process_group_exists(process.pid):
        orphaned_group = stop_process_group(process.pid)
    observer_wall_time_s = time.monotonic() - started
    gnu_time = parse_gnu_time(time_path)
    measurement_complete = (
        gnu_time["wall_time_s"] is not None
        and gnu_time["wall_time_s"] > 0
        and gnu_time["peak_rss_bytes"] is not None
        and gnu_time["peak_rss_bytes"] >= 0
        and gnu_time["exit_status"] is not None
    )
    return {
        "timed_command": timed_command,
        "wall_time_s": gnu_time["wall_time_s"],
        "wall_time_provenance": "gnu_time_elapsed",
        "observer_wall_time_s": observer_wall_time_s,
        "returncode": returncode,
        "timed_out": timed_out,
        "orphaned_process_group_cleaned": orphaned_group,
        "measurement_complete": measurement_complete,
        "gnu_time": gnu_time,
    }


def prewarm_prefix(input_record, depth):
    prefix = input_record["prefixes"][depth]
    started = time.monotonic()
    observed_sha256 = sha256_file(input_record["path"], prefix["size_bytes"])
    if observed_sha256 != prefix["sha256"]:
        raise ValueError(f"Input prefix changed before invocation: {input_record['id']} @ {depth}")
    return {
        "method": "sequential_read_and_sha256_outside_timing",
        "size_bytes": prefix["size_bytes"],
        "sha256": observed_sha256,
        "wall_time_s": time.monotonic() - started,
    }


def raw_output_files(output_dir):
    files = []
    for path in sorted(Path(output_dir).iterdir()):
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Output contains non-regular entry: {path}")
        files.append({"path": path.name, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return files


def directory_size_bytes(directory):
    total = 0
    for path in Path(directory).rglob("*"):
        if path.is_file() and not path.is_symlink():
            total += path.stat().st_size
    return total


def receipt_matches(receipt):
    if not isinstance(receipt, dict):
        return False
    try:
        path = verify_regular_file(receipt.get("path"), "preserved artifact")
    except (TypeError, ValueError):
        return False
    return receipt.get("size_bytes") == path.stat().st_size and receipt.get("sha256") == sha256_file(path)


def validate_preserved_result(result, signature, output_root, invocation_id):
    if result.get("status") not in {"failed", "timed_complete", "complete"}:
        raise ValueError(f"Preserved result has invalid status for {invocation_id}")
    if result.get("signature") != signature:
        raise ValueError(f"Preserved result signature differs for {invocation_id}")
    attempt_id = require_string(result.get("attempt_id"), "preserved result attempt_id")
    expected_attempt_dir = (Path(output_root) / "attempts" / invocation_id / attempt_id).resolve()
    attempt_dir = Path(require_string(result.get("attempt_dir"), "preserved result attempt_dir")).resolve()
    if attempt_dir != expected_attempt_dir or not attempt_dir.is_dir() or attempt_dir.is_symlink():
        raise ValueError(f"Preserved result attempt directory is invalid for {invocation_id}")
    logs = require_mapping(result.get("logs"), "preserved result logs")
    if not all(receipt_matches(logs.get(name)) for name in ("stdout", "stderr")):
        raise ValueError(f"Preserved result logs changed for {invocation_id}")
    if logs.get("gnu_time") is not None and not receipt_matches(logs["gnu_time"]):
        raise ValueError(f"Preserved result GNU time log changed for {invocation_id}")
    output_dir = attempt_dir / "output"
    if inventory_untrusted_directory(output_dir) != result.get("unvalidated_output_inventory"):
        raise ValueError(f"Preserved result output inventory changed for {invocation_id}")
    if result.get("timing_status") != "complete":
        return result
    stats_path = verify_regular_file(result.get("stats_path"), "preserved result stats")
    if stats_path.resolve().parent != (attempt_dir / "output").resolve():
        raise ValueError(f"Preserved result stats path is outside its attempt for {invocation_id}")
    if result.get("stats_sha256") != sha256_file(stats_path):
        raise ValueError(f"Preserved result stats changed for {invocation_id}")
    expected_outputs = result.get("raw_output_files")
    if not isinstance(expected_outputs, list) or raw_output_files(output_dir) != expected_outputs:
        raise ValueError(f"Preserved result outputs changed for {invocation_id}")
    return result


def normalized_metrics(stats, wall_time_s):
    metrics = {
        "n_reads_read": stats.get("n_reads_read"),
        "n_bases_read": stats.get("n_bases_read"),
        "n_subreads_ingested": stats.get("n_subreads_ingested"),
        "n_bases_ingested": stats.get("n_bases_ingested"),
        "n_kmers": stats.get("n_kmers"),
        "count_table_capacity": stats.get("count_table_capacity"),
        "allocator_peak_bytes": stats.get("peak_memory_bytes"),
        "stage_timings": stats.get("stage_timings"),
        "availability": {
            "count_table_capacity": "sharkmer_stats" if "count_table_capacity" in stats else None,
            "allocator_peak_bytes": "sharkmer_peak_alloc" if "peak_memory_bytes" in stats else None,
            "stage_timings": "sharkmer_stats" if "stage_timings" in stats else None,
        },
    }
    metrics["end_to_end_input_mbp_s"] = stats["n_bases_read"] / 1_000_000 / wall_time_s
    metrics["end_to_end_ingested_kmers_s"] = stats["n_kmers"] / wall_time_s
    metrics["availability"]["end_to_end_input_mbp_s"] = "gnu_time_elapsed_and_sharkmer_stats"
    metrics["availability"]["end_to_end_ingested_kmers_s"] = "gnu_time_elapsed_and_sharkmer_stats"
    return metrics


def invocation_signature(invocation, build, panel, input_record, settings):
    prefix = input_record["prefixes"][invocation["depth"]]
    return {
        "invocation": invocation,
        "binary_sha256": build["binary_sha256"],
        "source_commit": build["commit"],
        "panel_sha256": panel["sha256"],
        "input_sha256": input_record["sha256"],
        "input_subset": {
            "method": "first_records",
            "requested_records": invocation["depth"],
            "actual_records": prefix["records"],
            "actual_bases": prefix["bases"],
            "size_bytes": prefix["size_bytes"],
            "sha256": prefix["sha256"],
            "source_exhausted": prefix["records"] < invocation["depth"],
        },
        "settings": settings,
    }


def run_invocation(invocation, build, panel, input_record, settings, output_root, runner):
    invocation_id = invocation["invocation_id"]
    result_path = output_root / "results" / f"{invocation_id}.json"
    sample_prefix = f"{invocation['panel']}_{invocation['input']}_{invocation['depth']}"
    invocation = {**invocation, "sample_prefix": sample_prefix}
    signature = invocation_signature(invocation, build, panel, input_record, settings)
    if result_path.is_file():
        existing = load_json(result_path)
        return validate_preserved_result(existing, signature, output_root, invocation_id)
    attempt_id = uuid.uuid4().hex
    attempt_dir = output_root / "attempts" / invocation_id / attempt_id
    output_dir = attempt_dir / "output"
    temporary_dir = attempt_dir / "tmp"
    output_dir.mkdir(parents=True)
    temporary_dir.mkdir()
    full_input_started = time.monotonic()
    full_input_sha256 = sha256_file(input_record["path"])
    full_input_verification = {
        "method": "full_file_sha256_outside_timing_before_prefix_prewarm",
        "sha256": full_input_sha256,
        "size_bytes": Path(input_record["path"]).stat().st_size,
        "wall_time_s": time.monotonic() - full_input_started,
    }
    if full_input_sha256 != input_record["sha256"]:
        raise ValueError(f"Full input changed before invocation: {input_record['id']}")
    prewarm = prewarm_prefix(input_record, invocation["depth"])
    binary_command = [
        build["binary_path"],
        "-k",
        str(settings["k"]),
        "-t",
        str(settings["threads"]),
        "--chunks",
        str(settings["chunks"]),
        "--max-reads",
        str(invocation["depth"]),
        "-o",
        str(output_dir) + "/",
        "-s",
        sample_prefix,
        "--pcr-panel-file",
        panel["path"],
        input_record["path"],
    ]
    stdout_path = attempt_dir / "stdout.log"
    stderr_path = attempt_dir / "stderr.log"
    time_path = attempt_dir / "gnu-time.txt"
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    environment["TMPDIR"] = str(temporary_dir)
    started_at = utc_now()
    execution = run_timed(
        binary_command,
        stdout_path,
        stderr_path,
        time_path,
        settings["timeout_seconds"],
        settings["address_space_limit_bytes"],
        settings["cpu_list"],
        environment,
    )
    post_run_input_sha256 = sha256_file(input_record["path"], input_record["prefixes"][invocation["depth"]]["size_bytes"])
    if post_run_input_sha256 != prewarm["sha256"]:
        execution["input_changed_during_invocation"] = True
    else:
        execution["input_changed_during_invocation"] = False
    artifact_integrity = {
        "binary_sha256": sha256_file(build["binary_path"]),
        "panel_sha256": sha256_file(panel["path"]),
    }
    execution["binary_changed_during_invocation"] = artifact_integrity["binary_sha256"] != build["binary_sha256"]
    execution["panel_changed_during_invocation"] = artifact_integrity["panel_sha256"] != panel["sha256"]
    execution["measurement_matches_exit"] = execution["gnu_time"]["exit_status"] == execution["returncode"]
    log_receipts = {
        "stdout": file_receipt(stdout_path),
        "stderr": file_receipt(stderr_path),
        "gnu_time": file_receipt(time_path) if time_path.is_file() else None,
    }
    common = {
        "schema_version": 1,
        "signature": signature,
        "attempt_id": attempt_id,
        "attempt_dir": str(attempt_dir),
        "started_at": started_at,
        "finished_at": utc_now(),
        "full_input_verification": full_input_verification,
        "prewarm": prewarm,
        "binary_command": binary_command,
        "environment": {"LC_ALL": "C", "TMPDIR": str(temporary_dir)},
        "resource_limits": {
            "address_space_limit_bytes": settings["address_space_limit_bytes"],
            "address_space_limit_kind": "RLIMIT_AS",
            "timeout_seconds": settings["timeout_seconds"],
            "cpu_list": settings["cpu_list"],
        },
        "execution": execution,
        "logs": log_receipts,
        "post_run_input_prefix_sha256": post_run_input_sha256,
        "post_run_artifact_integrity": artifact_integrity,
        "unvalidated_output_inventory": inventory_untrusted_directory(output_dir),
        "temporary_directory_final_bytes": directory_size_bytes(temporary_dir),
        "temporary_directory_measurement": "final_size_not_peak",
    }
    integrity_failure = any(
        execution[field]
        for field in ("input_changed_during_invocation", "binary_changed_during_invocation", "panel_changed_during_invocation")
    )
    if (
        execution["timed_out"]
        or execution["returncode"] != 0
        or integrity_failure
        or execution["orphaned_process_group_cleaned"]
        or not execution["measurement_complete"]
        or not execution["measurement_matches_exit"]
    ):
        if integrity_failure:
            failure = "artifact_changed"
        elif execution["orphaned_process_group_cleaned"]:
            failure = "orphaned_process_group"
        elif execution["timed_out"]:
            failure = "timeout"
        elif execution["returncode"] != 0:
            failure = "nonzero_exit"
        elif not execution["measurement_complete"] or not execution["measurement_matches_exit"]:
            failure = "invalid_gnu_time_measurement"
        result = {**common, "status": "failed", "timing_status": "failed", "failure": failure}
        atomic_write_json(result_path, result)
        return result
    try:
        stats_path = output_dir / f"{sample_prefix}.stats.yaml"
        stats = parse_stats(stats_path)
        expected_version = build["binary_version"].split()[1]
        pcr_results = validate_shared_stats(
            stats, invocation, panel, input_record, binary_command, expected_version
        )
        if invocation["version"] == "baseline":
            gene_results, completion = validate_legacy_outputs(output_dir, invocation, panel, stats, pcr_results)
        else:
            gene_results, completion = validate_current_outputs(
                output_dir, invocation, panel, stats, pcr_results, binary_command, runner
            )
        result = {
            **common,
            "status": "timed_complete",
            "timing_status": "complete",
            "completion": completion,
            "stats_path": str(stats_path),
            "stats_sha256": sha256_file(stats_path),
            "metrics": normalized_metrics(stats, execution["wall_time_s"]),
            "genes": gene_results,
            "raw_output_files": raw_output_files(output_dir),
            "final_output_bytes": sum(path.stat().st_size for path in output_dir.iterdir() if path.is_file()),
        }
    except Exception as error:
        result = {
            **common,
            "status": "failed",
            "timing_status": "failed",
            "failure": "invalid_output",
            "error": f"{type(error).__name__}: {error}",
        }
    atomic_write_json(result_path, result)
    return result


def finalize_classification(invocation, result, panel, output_root, blast, db_path):
    if result.get("status") in {"complete", "failed"}:
        return result
    if result.get("status") != "timed_complete" or result.get("timing_status") != "complete":
        raise ValueError(f"Invocation is not ready for classification: {invocation['invocation_id']}")
    result_path = Path(output_root) / "results" / f"{invocation['invocation_id']}.json"
    genes = copy.deepcopy(result["genes"])
    try:
        classified_genes = evaluate_products(genes, panel, invocation["taxon"], db_path, blast)
        finalized = {
            **result,
            "status": "complete",
            "classification_status": "complete",
            "classification_finished_at": utc_now(),
            "classification_provenance": "current_validator_after_all_timed_invocations",
            "genes": classified_genes,
        }
    except Exception as error:
        finalized = {
            **result,
            "status": "failed",
            "classification_status": "failed",
            "classification_finished_at": utc_now(),
            "failure": "classification_failed",
            "classification_error": f"{type(error).__name__}: {error}",
        }
    atomic_write_json(result_path, finalized)
    return finalized


def product_signature(result):
    products = []
    for gene_result in result.get("genes", []):
        for product in gene_result.get("products", []):
            products.append((gene_result["gene"], product["product_index"], product["length"], product["sha256"]))
    return sorted(products)


def classification_signature(result):
    classifications = []
    for gene_result in result.get("genes", []):
        for product in gene_result.get("products", []):
            classifications.append(
                (
                    gene_result["gene"],
                    product["product_index"],
                    product["sha256"],
                    product.get("reference_match", {}).get("status"),
                )
            )
    return sorted(classifications)


def summarize_results(schedule, results):
    by_pair = {}
    for invocation in schedule:
        result = results[invocation["invocation_id"]]
        pair_key = (invocation["cell"], invocation["pair_index"])
        by_pair.setdefault(pair_key, {})[invocation["version"]] = result
    comparisons = []
    for (cell, pair_index), versions in sorted(by_pair.items()):
        baseline = versions.get("baseline")
        candidate = versions.get("candidate")
        complete = baseline is not None and candidate is not None and baseline.get("status") == candidate.get("status") == "complete"
        comparison = {"cell": cell, "pair_index": pair_index, "both_complete": complete}
        if complete:
            baseline_metrics = baseline["metrics"]
            candidate_metrics = candidate["metrics"]
            count_fields = ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers")
            comparison.update(
                {
                    "counts_identical": all(baseline_metrics[field] == candidate_metrics[field] for field in count_fields),
                    "products_identical": product_signature(baseline) == product_signature(candidate),
                    "classifications_identical": classification_signature(baseline) == classification_signature(candidate),
                    "baseline_wall_time_s": baseline["execution"]["wall_time_s"],
                    "candidate_wall_time_s": candidate["execution"]["wall_time_s"],
                    "candidate_speedup_ratio": baseline["execution"]["wall_time_s"]
                    / candidate["execution"]["wall_time_s"],
                    "baseline_peak_rss_bytes": baseline["execution"]["gnu_time"]["peak_rss_bytes"],
                    "candidate_peak_rss_bytes": candidate["execution"]["gnu_time"]["peak_rss_bytes"],
                    "baseline_input_mbp_s": baseline_metrics["end_to_end_input_mbp_s"],
                    "candidate_input_mbp_s": candidate_metrics["end_to_end_input_mbp_s"],
                    "baseline_ingested_kmers_s": baseline_metrics["end_to_end_ingested_kmers_s"],
                    "candidate_ingested_kmers_s": candidate_metrics["end_to_end_ingested_kmers_s"],
                }
            )
        comparisons.append(comparison)
    timing_groups = {}
    for invocation in schedule:
        result = results[invocation["invocation_id"]]
        if result.get("status") != "complete":
            continue
        key = (invocation["cell"], invocation["version"])
        timing_groups.setdefault(key, []).append(result["execution"]["wall_time_s"])
    timings = []
    for (cell, version_name), values in sorted(timing_groups.items()):
        timings.append(
            {
                "cell": cell,
                "version": version_name,
                "n": len(values),
                "values_s": values,
                "median_s": statistics.median(values),
                "minimum_s": min(values),
                "maximum_s": max(values),
            }
        )
    failures = [invocation_id for invocation_id, result in results.items() if result.get("status") != "complete"]
    return {
        "schema_version": 1,
        "generated_at": utc_now(),
        "invocations": len(schedule),
        "failures": failures,
        "pair_comparisons": comparisons,
        "timings": timings,
        "limitations": [
            "One-million-record cells have three paired repetitions; deeper cells have one paired run.",
            "Input files are prewarmed outside every timed invocation; operating-system cache residency is not guaranteed.",
            "Peak RSS is GNU time process RSS; the 40 GiB RLIMIT_AS ceiling is an address-space limit, not an RSS cap.",
            "BLAST classification uses the current validator after timing and is not a biological truth claim.",
            "Legacy 3.1.0 lacks current source-plan, stage-timing, capacity, and output-transaction fields; unavailable metrics remain null.",
        ],
    }


def preflight_panels(builds, panels, output_root):
    preflight_dir = output_root / "preflight"
    preflight_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for version_name, build in builds.items():
        for panel_name, panel in panels.items():
            command = [build["binary_path"], "--validate-panels", "--pcr-panel-file", panel["path"]]
            completed = subprocess.run(command, capture_output=True, text=True, timeout=120)
            result = {
                "version": version_name,
                "panel": panel_name,
                "command": command,
                "returncode": completed.returncode,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
            }
            results.append(result)
            atomic_write_json(preflight_dir / f"{version_name}_{panel_name}.json", result)
            if completed.returncode != 0:
                raise ValueError(f"Panel preflight failed for {version_name}/{panel_name}")
    return results


def prepare_reference_databases(panels, blast, output_root):
    has_references = any(blast.extract_references(panel["data"]) for panel in panels.values())
    if has_references and (not blast.check_blastn_available() or not blast.check_makeblastdb_available()):
        raise ValueError("Current per-product validation requires blastn and makeblastdb")
    databases = {}
    provenance = {}
    for panel_name, panel in panels.items():
        db_dir = output_root / "reference_databases" / panel_name
        db_dir.mkdir(parents=True, exist_ok=True)
        db_path = blast.build_reference_db(panel["data"], db_dir)
        if blast.extract_references(panel["data"]) and db_path is None:
            raise ValueError(f"Reference database build failed for panel {panel_name}")
        databases[panel_name] = db_path
        provenance[panel_name] = blast.reference_checksums(panel["data"])
    return databases, provenance


def run_driver(arguments):
    initial_protocol = load_json(arguments.protocol)
    output_root_text = str(arguments.output_root) if arguments.output_root else require_string(initial_protocol.get("output_root"), "output_root")
    output_root = Path(output_root_text).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    frozen_receipts = {
        "driver": freeze_receipt(__file__, output_root / "receipts" / "driver.py"),
        "protocol": freeze_receipt(arguments.protocol, output_root / "receipts" / "protocol.json"),
        "builds": freeze_receipt(arguments.builds, output_root / "receipts" / "builds.json"),
        "inputs": freeze_receipt(arguments.inputs, output_root / "receipts" / "inputs.json"),
    }
    protocol_document = load_json(frozen_receipts["protocol"]["path"])
    builds_document = load_json(frozen_receipts["builds"]["path"])
    inputs_document = load_json(frozen_receipts["inputs"]["path"])
    builds = verify_builds(builds_document)
    validator_root = require_string(protocol_document.get("validator_root"), "validator_root")
    if Path(validator_root).resolve() != Path(builds["candidate"]["source_export"]).resolve():
        raise ValueError("Validator root must be the frozen candidate source export")
    runner, blast, validator_provenance = load_validator(validator_root)
    inputs = verify_inputs(inputs_document)
    protocol = verify_protocol(protocol_document, inputs, runner)
    candidate_panels = Path(builds["candidate"]["source_export"]).resolve() / "panels"
    for panel in protocol["panels"].values():
        if Path(panel["path"]).resolve().parent != candidate_panels:
            raise ValueError("Frozen panels must come directly from the candidate source export")
    schedule = build_schedule(protocol)
    if len(schedule) != 114:
        raise ValueError(f"Frozen protocol produced {len(schedule)} invocations instead of 114")
    schedule_document = {
        "schema_version": 1,
        "comparison_id": protocol["comparison_id"],
        "created_at": utc_now(),
        "protocol_sha256": frozen_receipts["protocol"]["sha256"],
        "builds_sha256": frozen_receipts["builds"]["sha256"],
        "inputs_sha256": frozen_receipts["inputs"]["sha256"],
        "schedule": schedule,
    }
    schedule_path = output_root / "schedule.json"
    if schedule_path.exists():
        existing = load_json(schedule_path)
        comparable_existing = {key: value for key, value in existing.items() if key != "created_at"}
        comparable_new = {key: value for key, value in schedule_document.items() if key != "created_at"}
        if comparable_existing != comparable_new:
            raise ValueError("Existing schedule differs from requested immutable schedule")
    else:
        atomic_write_json(schedule_path, schedule_document)
    preflight = preflight_panels(builds, protocol["panels"], output_root)
    reference_provenance = {
        panel_name: blast.reference_checksums(panel["data"])
        for panel_name, panel in protocol["panels"].items()
    }
    if any(blast.extract_references(panel["data"]) for panel in protocol["panels"].values()):
        if not blast.check_blastn_available() or not blast.check_makeblastdb_available():
            raise ValueError("Current per-product validation requires blastn and makeblastdb")
    suite_provenance = {
        "schema_version": 1,
        "comparison_id": protocol["comparison_id"],
        "created_at": utc_now(),
        "effective_output_root": str(output_root),
        "driver_path": frozen_receipts["driver"]["path"],
        "driver_sha256": frozen_receipts["driver"]["sha256"],
        "protocol": protocol_document,
        "builds": builds,
        "inputs": inputs_document,
        "frozen_receipts": frozen_receipts,
        "validator": validator_provenance,
        "physical_cores": protocol["physical_cores"],
        "reference_checksums": reference_provenance,
        "preflight": preflight,
        "timing_boundary": "all_sharkmer_invocations_complete_before_reference_database_build_or_blast_classification",
    }
    suite_provenance_path = output_root / "suite-provenance.json"
    if suite_provenance_path.exists():
        existing_provenance = load_json(suite_provenance_path)
        comparable_existing = {key: value for key, value in existing_provenance.items() if key != "created_at"}
        comparable_new = {key: value for key, value in suite_provenance.items() if key != "created_at"}
        if comparable_existing != comparable_new:
            raise ValueError("Existing suite provenance differs from frozen comparison")
    else:
        atomic_write_json(suite_provenance_path, suite_provenance)
    results = {}
    for invocation in schedule:
        panel = protocol["panels"][invocation["panel"]]
        input_record = inputs[invocation["input"]]
        build = builds[invocation["version"]]
        result = run_invocation(
            invocation,
            build,
            panel,
            input_record,
            protocol["settings"],
            output_root,
            runner,
        )
        results[invocation["invocation_id"]] = result
        print(f"{invocation['invocation_id']}: timing {result['status']}", flush=True)
    databases, observed_reference_provenance = prepare_reference_databases(protocol["panels"], blast, output_root)
    if observed_reference_provenance != reference_provenance:
        raise ValueError("Reference provenance changed after timed invocations")
    for invocation in schedule:
        invocation_id = invocation["invocation_id"]
        panel = protocol["panels"][invocation["panel"]]
        results[invocation_id] = finalize_classification(
            invocation,
            results[invocation_id],
            panel,
            output_root,
            blast,
            databases[invocation["panel"]],
        )
        print(f"{invocation_id}: classification {results[invocation_id]['status']}", flush=True)
    summary = summarize_results(schedule, results)
    atomic_write_json(output_root / "comparison.json", summary)
    return 1 if summary["failures"] else 0


def main():
    parser = argparse.ArgumentParser(description="Compare released Sharkmer 3.1.0 with a current candidate")
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--builds", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    arguments = parser.parse_args()
    try:
        return run_driver(arguments)
    except Exception as error:
        print(f"driver error: {type(error).__name__}: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

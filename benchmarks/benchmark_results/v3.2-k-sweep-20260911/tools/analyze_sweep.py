#!/usr/bin/env python3
import argparse
import collections
import hashlib
import json
import statistics
from pathlib import Path

import yaml


ROLES = {"baseline", "candidate"}
COUNT_NAMES = (
    "n_reads_read",
    "n_bases_read",
    "n_subreads_ingested",
    "n_bases_ingested",
    "n_kmers",
)


def load_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def json_digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def target_scopes(protocol):
    registered = {}
    for panel, entries in protocol.get("target_metadata", {}).items():
        if not isinstance(panel, str) or not isinstance(entries, list):
            raise ValueError("Invalid target metadata")
        panel_targets = registered.setdefault(panel, {})
        for entry in entries:
            if not isinstance(entry, dict):
                raise ValueError("Invalid target metadata entry")
            gene, scope = entry.get("gene"), entry.get("scope")
            if not isinstance(gene, str) or scope not in {"high_copy_candidate", "deferred_or_unclassified"}:
                raise ValueError("Target metadata needs a registered gene and scope")
            if gene in panel_targets:
                raise ValueError(f"Duplicate target metadata gene: {panel}/{gene}")
            panel_targets[gene] = entry
    return registered


def product_key(product):
    gene, digest = product.get("gene"), product.get("sha256")
    if not isinstance(gene, str) or not isinstance(digest, str) or len(digest) != 64:
        raise ValueError("Product requires normalized gene and SHA-256")
    return gene, digest


def product_map(product_values):
    products, duplicates = {}, []
    for product in product_values:
        key = product_key(product)
        if key in products:
            duplicates.append(product)
        else:
            products[key] = product
    return products, duplicates


def reference_signature(product):
    return json.dumps(product.get("reference_match"), sort_keys=True, separators=(",", ":"))


def comparison(anchor_values, observed_values):
    anchors, anchor_duplicates = product_map(anchor_values)
    observed, observed_duplicates = product_map(observed_values)
    retained_keys = sorted(set(anchors) & set(observed))
    lost = [anchors[key] for key in sorted(set(anchors) - set(observed))]
    gained = [observed[key] for key in sorted(set(observed) - set(anchors))]
    same_gene_changes = []
    for product in lost:
        alternatives = [candidate for candidate in observed.values() if candidate["gene"] == product["gene"]]
        if alternatives:
            same_gene_changes.append({"anchor": product, "observed_same_gene": alternatives})
    scope_status_counts = {}
    for scope in ("high_copy_candidate", "deferred_or_unclassified"):
        scope_status_counts[scope] = {
            "exact_retained": sum(product.get("scope") == scope for product in (anchors[key] for key in retained_keys)),
            "lost": sum(product.get("scope") == scope for product in lost),
            "gained": sum(product.get("scope") == scope for product in gained),
        }
    return {
        "exact_retained": [anchors[key] for key in retained_keys],
        "lost": lost,
        "gained": gained,
        "changed_endpoint_or_sequence_unresolved": same_gene_changes,
        "anchor_duplicates": anchor_duplicates,
        "observed_duplicates": observed_duplicates,
        "scope_status_counts": scope_status_counts,
    }


def normalized_products(result, registered_targets, panel):
    if panel not in registered_targets or not isinstance(result.get("genes"), list):
        raise ValueError("Result has an unknown panel or no normalized genes")
    values = []
    for gene_result in result["genes"]:
        if not isinstance(gene_result, dict):
            raise ValueError("Invalid gene result")
        gene, product_values = gene_result.get("gene"), gene_result.get("products")
        if gene not in registered_targets[panel] or not isinstance(product_values, list):
            raise ValueError("Result product is absent from registered target metadata")
        for product in product_values:
            if not isinstance(product, dict):
                raise ValueError("Invalid product")
            value = dict(product)
            value["gene"] = gene
            value["scope"] = registered_targets[panel][gene]["scope"]
            product_key(value)
            values.append(value)
    return values


def result_metrics(result):
    execution, metrics = result.get("execution"), result.get("metrics")
    if not isinstance(execution, dict) or not isinstance(metrics, dict):
        raise ValueError("Result lacks execution or count metrics")
    if execution.get("returncode") != 0 or execution.get("timed_out") is not False:
        raise ValueError("Result execution is incomplete")
    wall_time = execution.get("wall_time_s")
    peak_rss = execution.get("gnu_time", {}).get("peak_rss_bytes")
    if not isinstance(wall_time, (int, float)) or wall_time <= 0 or type(peak_rss) is not int or peak_rss < 0:
        raise ValueError("Result wall time or RSS is invalid")
    counts = {name: metrics.get(name) for name in COUNT_NAMES}
    if any(type(value) is not int or value < 0 for value in counts.values()):
        raise ValueError("Result count metrics are invalid")
    return {"wall_time_s": wall_time, "peak_rss_bytes": peak_rss, "counts": counts}


def raw_stats_diagnostics(result, role):
    stats_path = Path(result.get("stats_path", ""))
    expected_digest = result.get("stats_sha256")
    if not isinstance(expected_digest, str) or not stats_path.is_file() or stats_path.is_symlink():
        raise ValueError("Result lacks a regular checksummed stats file")
    contents = stats_path.read_bytes()
    if hashlib.sha256(contents).hexdigest() != expected_digest:
        raise ValueError("Stats receipt differs from preserved raw stats")
    receipts = result.get("raw_output_files")
    if not isinstance(receipts, list) or not any(
        receipt.get("path") == stats_path.name and receipt.get("sha256") == expected_digest
        for receipt in receipts
        if isinstance(receipt, dict)
    ):
        raise ValueError("Raw output inventory lacks matching stats receipt")
    stats = yaml.safe_load(contents)
    pcr_results = stats.get("pcr_results") if isinstance(stats, dict) else None
    if not isinstance(pcr_results, list):
        raise ValueError("Raw stats lacks PCR result diagnostics")
    diagnostics = []
    for gene_result in pcr_results:
        if not isinstance(gene_result, dict) or not isinstance(gene_result.get("gene_name"), str):
            raise ValueError("Raw stats has malformed PCR result")
        values = gene_result.get("threshold_diagnostics")
        if values is not None and not isinstance(values, list):
            raise ValueError("Raw stats has malformed threshold diagnostics")
        diagnostics.append({"gene_name": gene_result["gene_name"], "threshold_diagnostics": values})
    return {"stats_path": str(stats_path), "stats_sha256": expected_digest, "genes": diagnostics}


def execution_directories(root):
    children = sorted(path for path in root.glob("k*/provenance.json") if path.parent.is_dir())
    if children:
        return [path.parent for path in children]
    if (root / "provenance.json").is_file():
        return [root]
    raise ValueError(f"No execution provenance beneath: {root}")


def invocation_kmer_length(invocation, default_kmer_length=None):
    value = invocation.get("kmer_length", invocation.get("k", default_kmer_length))
    if type(value) is not int or value < 2:
        raise ValueError("Signed invocation requires kmer_length")
    return value


def read_execution(root, protocol, require_diagnostics, default_kmer_length=None):
    registered_targets = target_scopes(protocol)
    all_records = {}
    for execution in execution_directories(root):
        provenance = load_json(execution / "provenance.json")
        if provenance.get("protocol") != protocol:
            raise ValueError(f"Execution protocol differs from frozen protocol: {execution}")
        schedule = provenance.get("schedule")
        if not isinstance(schedule, list):
            raise ValueError(f"Execution lacks a schedule: {execution}")
        expected = {entry.get("invocation_id"): entry for entry in schedule if isinstance(entry, dict)}
        result_paths = {path.stem: path for path in (execution / "results").glob("*.json")}
        if len(expected) != len(schedule) or None in expected or set(expected) != set(result_paths):
            raise ValueError(f"Results do not exactly match frozen schedule: {execution}")
        for identity, invocation in expected.items():
            result = load_json(result_paths[identity])
            signed = result.get("signature", {}).get("invocation")
            if not isinstance(signed, dict) or any(signed.get(name) != value for name, value in invocation.items()):
                raise ValueError(f"Signed invocation differs: {identity}")
            if result.get("status") != "complete" or result.get("classification_status") != "complete":
                raise ValueError(f"Unfinished classification: {identity}")
            role, panel, cell, pair_index = invocation.get("version"), invocation.get("panel"), invocation.get("cell"), invocation.get("pair_index")
            if role not in ROLES or panel not in registered_targets or not isinstance(cell, str) or type(pair_index) is not int:
                raise ValueError(f"Invalid signed invocation identity: {identity}")
            diagnostics = raw_stats_diagnostics(result, role) if require_diagnostics else None
            kmer_length = invocation_kmer_length(invocation, default_kmer_length)
            key = cell, kmer_length, pair_index, role
            if key in all_records:
                raise ValueError(f"Duplicate execution record: {key}")
            all_records[key] = {
                "identity": identity,
                "invocation": invocation,
                "panel": panel,
                "signature": result.get("signature"),
                "metrics": result_metrics(result),
                "products": normalized_products(result, registered_targets, panel),
                "threshold_diagnostics": diagnostics,
            }
    return all_records


def high_copy_keys(product_values):
    return {key for key, product in product_map(product_values)[0].items() if product["scope"] == "high_copy_candidate"}


def confirmed_keys(product_values):
    return {
        key
        for key, product in product_map(product_values)[0].items()
        if isinstance(product.get("reference_match"), dict) and product["reference_match"].get("status") == "confirmed_product"
    }


def is_confirmed(product):
    return isinstance(product.get("reference_match"), dict) and product["reference_match"].get("status") == "confirmed_product"


def classification_stable(product_collections):
    return len({tuple(sorted((product_key(product), reference_signature(product)) for product in products)) for products in product_collections}) == 1


def source_identity(record):
    signature = record["signature"]
    return {
        name: signature.get(name)
        for name in ("source_commit", "binary_sha256", "input_sha256", "input_subset", "panel_sha256", "settings")
    }


def validate_discovery_matrix(records, protocol):
    samples = protocol.get("samples")
    expected_commits = protocol.get("expected_commits")
    allowed_k = protocol.get("allowed_k")
    if not isinstance(samples, list) or not isinstance(expected_commits, dict):
        raise ValueError("Discovery protocol lacks samples or expected commits")
    expected = set()
    for sample in samples:
        if not isinstance(sample, dict):
            raise ValueError("Discovery protocol has invalid sample")
        panel, input_name, depth, pairs = sample.get("panel"), sample.get("input"), sample.get("depth"), sample.get("pairs")
        if not isinstance(panel, str) or not isinstance(input_name, str) or type(depth) is not int or type(pairs) is not int or pairs < 1:
            raise ValueError("Discovery protocol has invalid sample identity")
        cell = f"{panel}/{input_name}/{depth}"
        for kmer_length in allowed_k:
            for pair_index in range(1, pairs + 1):
                for role in ROLES:
                    expected.add((cell, kmer_length, pair_index, role))
    if set(records) != expected:
        raise ValueError("Discovery records differ from the preregistered sample/k/pair matrix")
    binary_hashes = collections.defaultdict(set)
    for (_, kmer_length, _, role), record in records.items():
        signature = record["signature"]
        if not isinstance(signature, dict):
            raise ValueError("Result lacks source signature")
        if signature.get("source_commit") != expected_commits.get(role):
            raise ValueError("Result source commit differs from version role")
        if signature.get("settings", {}).get("k") != kmer_length:
            raise ValueError("Result signature k differs from scheduled k")
        if not isinstance(signature.get("binary_sha256"), str):
            raise ValueError("Result lacks binary receipt")
        if not isinstance(signature.get("input_sha256"), str) or not isinstance(signature.get("input_subset"), dict):
            raise ValueError("Result lacks fixed input receipts")
        if not isinstance(signature.get("panel_sha256"), str):
            raise ValueError("Result lacks panel receipt")
        binary_hashes[role].add(signature["binary_sha256"])
    if any(len(binary_hashes[role]) != 1 for role in ROLES):
        raise ValueError("Version role used multiple binary receipts")


def validate_fixed_inputs(discovery_records, fixed_records):
    for (cell, _, pair_index, role), record in discovery_records.items():
        fixed = fixed_records.get((cell, 19, pair_index, role))
        if fixed is None:
            raise ValueError(f"Fixed input anchor missing: {cell}/{pair_index}/{role}")
        current_signature, fixed_signature = record["signature"], fixed["signature"]
        for name in ("source_commit", "binary_sha256", "input_sha256", "input_subset", "panel_sha256"):
            if current_signature.get(name) != fixed_signature.get(name):
                raise ValueError(f"Input or panel receipt differs from fixed k19: {cell}/{pair_index}/{role}")
        for name in ("n_reads_read", "n_bases_read"):
            if record["metrics"]["counts"][name] != fixed["metrics"]["counts"][name]:
                raise ValueError(f"Input count differs from fixed k19: {cell}/{pair_index}/{role}")


def fixed_reference(protocol):
    reference = protocol.get("fixed_reference")
    if not isinstance(reference, dict):
        raise ValueError("Cross-k protocol lacks fixed reference")
    root = Path(reference.get("execution", ""))
    provenance = load_json(root / "provenance.json")
    source_protocol = provenance.get("protocol")
    if not isinstance(source_protocol, dict):
        raise ValueError("Fixed reference lacks its frozen protocol")
    expected_digest = reference.get("protocol_sha256")
    receipt_path = root / "receipts" / "protocol.json"
    if not isinstance(expected_digest, str) or not receipt_path.is_file():
        raise ValueError("Fixed reference lacks its protocol receipt")
    if hashlib.sha256(receipt_path.read_bytes()).hexdigest() != expected_digest:
        raise ValueError("Fixed reference protocol receipt digest differs")
    if load_json(receipt_path) != source_protocol:
        raise ValueError("Fixed reference protocol receipt differs from provenance")
    diagnostic_path = Path(reference.get("analysis", "")) / "diagnostic-audit.json"
    diagnostic_digest = reference.get("diagnostic_audit_sha256")
    if not isinstance(diagnostic_digest, str) or not diagnostic_path.is_file():
        raise ValueError("Fixed reference lacks diagnostic-audit receipt")
    if hashlib.sha256(diagnostic_path.read_bytes()).hexdigest() != diagnostic_digest:
        raise ValueError("Fixed reference diagnostic-audit receipt differs")
    baseline_k = reference.get("baseline_k")
    if type(baseline_k) is not int:
        raise ValueError("Fixed reference lacks baseline k")
    records = read_execution(root, source_protocol, require_diagnostics=False, default_kmer_length=baseline_k)
    candidate_high = set()
    losses = {}
    for (cell, kmer_length, pair_index, role), record in records.items():
        if role == "candidate":
            candidate_high.update(high_copy_keys(record["products"]))
        if kmer_length != baseline_k or role != "baseline":
            continue
        candidate = records.get((cell, kmer_length, pair_index, "candidate"))
        if candidate is None:
            raise ValueError(f"Fixed reference lacks candidate pair: {cell}/{pair_index}")
        baseline_products = product_map(record["products"])[0]
        candidate_products = product_map(candidate["products"])[0]
        for key in high_copy_keys(record["products"]) - high_copy_keys(candidate["products"]):
            losses[cell, key[0], key[1]] = baseline_products[key]
    if len(losses) != reference.get("expected_lost_high_copy_products"):
        raise ValueError(f"Fixed reference lost-product count differs: {len(losses)}")
    if len(candidate_high) != reference.get("expected_candidate_high_copy_products_all_ten_samples"):
        raise ValueError(f"Fixed reference candidate high-copy count differs: {len(candidate_high)}")
    its_digest = reference.get("exact_its2_sha256")
    its_keys = [key for _, gene, digest in losses if gene == "ITS_2" and digest == its_digest]
    if len(its_keys) != 1:
        raise ValueError("Fixed reference lacks the exact preregistered ITS_2 loss")
    return records, losses, ("ITS_2", its_digest)


def aggregate_cell(cell, kmer_length, pairs, fixed_records, known_losses):
    rows = []
    values_by_role = collections.defaultdict(list)
    for pair_index, role_records in sorted(pairs.items()):
        if set(role_records) != ROLES:
            raise ValueError(f"Incomplete version pair: {cell}/k{kmer_length}/{pair_index}")
        fixed_baseline = fixed_records.get((cell, 19, pair_index, "baseline"))
        fixed_candidate = fixed_records.get((cell, 19, pair_index, "candidate"))
        if fixed_baseline is None or fixed_candidate is None:
            raise ValueError(f"Fixed k19 pair missing: {cell}/{pair_index}")
        baseline, candidate = role_records["baseline"], role_records["candidate"]
        same_k = comparison(baseline["products"], candidate["products"])
        candidate_vs_released = comparison(fixed_baseline["products"], candidate["products"])
        candidate_vs_current = comparison(fixed_candidate["products"], candidate["products"])
        baseline_vs_released = comparison(fixed_baseline["products"], baseline["products"])
        baseline_vs_current = comparison(fixed_candidate["products"], baseline["products"])
        candidate_keys = set(product_map(candidate["products"])[0])
        restored = sorted(
            (gene, digest)
            for loss_cell, gene, digest in known_losses
            if loss_cell == cell and (gene, digest) in candidate_keys
        )
        rows.append({
            "pair_index": pair_index,
            "baseline": baseline["metrics"],
            "candidate": candidate["metrics"],
            "source_identity": {"baseline": source_identity(baseline), "candidate": source_identity(candidate)},
            "same_k_count_parity": baseline["metrics"]["counts"] == candidate["metrics"]["counts"],
            "same_k_membership": same_k,
            "baseline_vs_k19_released": baseline_vs_released,
            "baseline_vs_k19_candidate": baseline_vs_current,
            "candidate_vs_k19_released": candidate_vs_released,
            "candidate_vs_k19_candidate": candidate_vs_current,
            "candidate_exact_restored_known_losses": restored,
            "threshold_diagnostics": {
                "baseline": baseline["threshold_diagnostics"],
                "candidate": candidate["threshold_diagnostics"],
            },
        })
        values_by_role["baseline"].append(baseline["products"])
        values_by_role["candidate"].append(candidate["products"])
    stability = {role: classification_stable(values) for role, values in values_by_role.items()}
    return {
        "cell": cell,
        "kmer_length": kmer_length,
        "repetitions": len(rows),
        "pairs": rows,
        "median_wall_time_s": {role: statistics.median(row[role]["wall_time_s"] for row in rows) for role in ROLES},
        "median_peak_rss_bytes": {role: statistics.median(row[role]["peak_rss_bytes"] for row in rows) for role in ROLES},
        "repeat_sequence_and_classification_stable": stability,
    }


def option_summary(kmer_length, cells, its_key, global_stability, global_count_parity):
    rows = [row for cell in cells for row in cell["pairs"]]
    restored_every_replicate = set()
    retained_every_replicate = set()
    for cell in cells:
        cell_restorations = [set(row["candidate_exact_restored_known_losses"]) for row in cell["pairs"]]
        stable_restorations = set.intersection(*cell_restorations) if cell_restorations else set()
        restored_every_replicate.update((cell["cell"], gene, digest) for gene, digest in stable_restorations)
        cell_retention = [high_copy_keys(row["candidate_vs_k19_released"]["exact_retained"]) for row in cell["pairs"]]
        stable_retention = set.intersection(*cell_retention) if cell_retention else set()
        retained_every_replicate.update((cell["cell"], gene, digest) for gene, digest in stable_retention)
    confirmed_losses = [
        product
        for row in rows
        for product in row["candidate_vs_k19_candidate"]["lost"]
        if is_confirmed(product)
    ]
    stable = all(all(cell["repeat_sequence_and_classification_stable"].values()) for cell in cells)
    count_parity = all(row["same_k_count_parity"] for row in rows)
    eligible = kmer_length > 19 and global_stability and global_count_parity and bool(restored_every_replicate) and not confirmed_losses
    return {
        "kmer_length": kmer_length,
        "all_discovery_invocations_valid": True,
        "same_version_sequence_and_classification_stable": stable,
        "all_discovery_sequence_and_classification_stable": global_stability,
        "all_discovery_same_k_count_parity": global_count_parity,
        "same_k_count_parity_every_replicate": count_parity,
        "exact_restored_known_losses_every_replicate": sorted(restored_every_replicate),
        "exact_its2_ak281180_restored_every_replicate": any((gene, digest) == its_key for _, gene, digest in restored_every_replicate),
        "confirmed_k19_candidate_sequences_lost": confirmed_losses,
        "k19_released_high_copy_exact_retained_every_replicate": sorted(retained_every_replicate),
        "eligible": eligible,
    }


def analyze(execution):
    root = Path(execution).resolve()
    directories = execution_directories(root)
    provenance = load_json(directories[0] / "provenance.json")
    protocol = provenance.get("protocol")
    if not isinstance(protocol, dict) or protocol.get("baseline_role") != "baseline":
        raise ValueError("Discovery lacks the preregistered baseline protocol")
    allowed_k = protocol.get("allowed_k")
    if not isinstance(allowed_k, list) or 19 not in allowed_k or any(type(value) is not int for value in allowed_k):
        raise ValueError("Discovery protocol has invalid allowed k values")
    discovery_records = read_execution(root, protocol, require_diagnostics=True)
    observed_k = {key[1] for key in discovery_records}
    if set(allowed_k) != observed_k:
        raise ValueError(f"Discovery k values differ from protocol: {sorted(observed_k)}")
    validate_discovery_matrix(discovery_records, protocol)
    fixed_records, known_losses, its_key = fixed_reference(protocol)
    validate_fixed_inputs(discovery_records, fixed_records)
    grouped = collections.defaultdict(lambda: collections.defaultdict(dict))
    for (cell, kmer_length, pair_index, role), record in discovery_records.items():
        grouped[cell, kmer_length][pair_index][role] = record
    cell_results = [
        aggregate_cell(cell, kmer_length, pairs, fixed_records, known_losses)
        for (cell, kmer_length), pairs in sorted(grouped.items())
    ]
    global_stability = all(
        all(cell["repeat_sequence_and_classification_stable"].values())
        for cell in cell_results
    )
    global_count_parity = all(
        row["same_k_count_parity"]
        for cell in cell_results
        for row in cell["pairs"]
    )
    options = [
        option_summary(
            kmer_length,
            [cell for cell in cell_results if cell["kmer_length"] == kmer_length],
            its_key,
            global_stability,
            global_count_parity,
        )
        for kmer_length in sorted(observed_k)
        if kmer_length > 19
    ]
    options.sort(key=lambda option: (
        not option["eligible"],
        not option["exact_its2_ak281180_restored_every_replicate"],
        -len(option["exact_restored_known_losses_every_replicate"]),
        -len(option["k19_released_high_copy_exact_retained_every_replicate"]),
        option["kmer_length"],
    ))
    return {
        "schema_version": 1,
        "execution": str(root),
        "protocol_sha256": json_digest(protocol),
        "fixed_reference_execution": protocol["fixed_reference"]["execution"],
        "known_lost_high_copy_products": [known_losses[key] for key in sorted(known_losses)],
        "cells": cell_results,
        "confirmation_options_ranked": options,
        "interpretation": "Exact gene-plus-sequence membership is primary. Changed endpoint or sequence hashes remain unresolved evidence and never count as recovery. Eligibility is preregistered, requires every replicate, and does not choose by wall time or RSS.",
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze preregistered cross-k recovery sweep")
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise SystemExit(f"Output already exists: {arguments.output}")
    analysis = analyze(arguments.execution)
    arguments.output.mkdir()
    write_json(arguments.output / "analysis.json", analysis)
    print(arguments.output / "analysis.json")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import collections
import hashlib
import json
import re
import statistics
from decimal import Decimal, InvalidOperation
from pathlib import Path
COUNT_FIELDS = (
    "n_reads_read",
    "n_bases_read",
    "n_subreads_ingested",
    "n_bases_ingested",
    "n_kmers",
)
def load_json(path):
    with Path(path).open() as input_stream:
        value = json.load(input_stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value
def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as input_stream:
        for block in iter(lambda: input_stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
def require_mapping(value, label):
    if not isinstance(value, dict):
        raise ValueError(f"Expected mapping for {label}")
    return value
def header_fields(header, label):
    if not isinstance(header, str):
        raise ValueError(f"Product header missing for {label}")
    fields = {}
    for token in header.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    required = {"gene", "length", "kmer_count_median", "kmer_count_min"}
    if not required <= fields.keys():
        raise ValueError(f"Product header lacks required support fields for {label}")
    return fields
def product_evidence(panel, gene, product):
    product = require_mapping(product, f"product for {gene}")
    sha256 = product.get("sha256")
    length = product.get("length")
    if not isinstance(sha256, str) or re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
        raise ValueError(f"Product SHA-256 invalid for {gene}")
    if type(length) is not int or length <= 0:
        raise ValueError(f"Product length invalid for {gene}")
    fields = header_fields(product.get("header"), gene)
    try:
        header_length = int(fields["length"])
        median = Decimal(fields["kmer_count_median"])
        minimum = int(fields["kmer_count_min"])
    except (InvalidOperation, ValueError) as error:
        raise ValueError(f"Product support header invalid for {gene}") from error
    if fields["gene"] != f"{panel}_{gene}" or header_length != length or not median.is_finite() or median < 0 or minimum < 0:
        raise ValueError(f"Product support header disagrees with product for {gene}")
    reference_match = require_mapping(product.get("reference_match"), f"reference match for {gene}")
    if not isinstance(reference_match.get("status"), str) or reference_match["status"] == "failed_run":
        raise ValueError(f"Product classification missing for {gene}")
    return {
        "gene": gene,
        "sha256": sha256,
        "length": length,
        "kmer_count_median": fields["kmer_count_median"],
        "kmer_count_min": minimum,
        "reference_match": reference_match,
    }
def result_products(result, panel, scopes):
    values = []
    gene_results = result.get("genes")
    if not isinstance(gene_results, list):
        raise ValueError("Gene results missing")
    for gene_result in gene_results:
        gene_result = require_mapping(gene_result, "gene result")
        gene = gene_result.get("gene")
        if gene not in scopes:
            raise ValueError(f"Gene is absent from preregistered target metadata: {gene}")
        products = gene_result.get("products")
        if not isinstance(products, list):
            raise ValueError(f"Products missing for {gene}")
        if gene_result.get("n_products") != len(products):
            raise ValueError(f"Product count disagrees for {gene}")
        for product in products:
            evidence = product_evidence(panel, gene, product)
            evidence["scope"] = scopes[gene]["scope"]
            values.append(evidence)
    return values
def membership(values):
    unique, duplicates = {}, []
    for value in values:
        key = (value["gene"], value["sha256"])
        if key in unique:
            duplicates.append(value)
        else:
            unique[key] = value
    return unique, duplicates
def metrics(result):
    execution = require_mapping(result.get("execution"), "execution")
    gnu_time = require_mapping(execution.get("gnu_time"), "GNU time")
    if (
        execution.get("returncode") != 0
        or execution.get("timed_out") is not False
        or execution.get("orphaned_process_group_cleaned") is not False
        or execution.get("measurement_complete") is not True
        or execution.get("measurement_matches_exit") is not True
        or any(execution.get(field) is not False for field in ("input_changed_during_invocation", "binary_changed_during_invocation", "panel_changed_during_invocation"))
        or gnu_time.get("exit_status", 0) != 0
    ):
        raise ValueError("Execution integrity check failed")
    wall_time = execution.get("wall_time_s")
    peak_rss = gnu_time.get("peak_rss_bytes")
    if not isinstance(wall_time, (int, float)) or wall_time <= 0 or type(peak_rss) is not int or peak_rss < 0:
        raise ValueError("Timing or RSS is invalid")
    values = require_mapping(result.get("metrics"), "metrics")
    counts = {field: values.get(field) for field in COUNT_FIELDS}
    if any(type(value) is not int or value < 0 for value in counts.values()):
        raise ValueError("Aggregate count metric is invalid")
    return {"wall_time_s": wall_time, "peak_rss_bytes": peak_rss, "counts": counts}
def summary(values):
    return {"n": len(values), "median": statistics.median(values), "minimum": min(values), "maximum": max(values)}
def validate_result(result, invocation):
    signed = require_mapping(require_mapping(result.get("signature"), "signature").get("invocation"), "signed invocation")
    if any(signed.get(key) != value for key, value in invocation.items()):
        raise ValueError(f"Signed invocation differs for {invocation['invocation_id']}")
    expected_sample = f"{invocation['panel']}_{invocation['input']}_{invocation['depth']}"
    if signed.get("sample_prefix") != expected_sample:
        raise ValueError(f"Signed sample differs for {invocation['invocation_id']}")
    if result.get("status") != "complete" or result.get("timing_status") != "complete" or result.get("classification_status") != "complete":
        raise ValueError(f"Result is unfinished or unclassified: {invocation['invocation_id']}")
    if result.get("classification_error") or result.get("error"):
        raise ValueError(f"Result contains an evaluation error: {invocation['invocation_id']}")
def scope_map(protocol):
    scopes = {}
    for panel, entries in require_mapping(protocol.get("target_metadata"), "target_metadata").items():
        if not isinstance(panel, str) or not isinstance(entries, list) or not entries:
            raise ValueError(f"Invalid target metadata for {panel}")
        for entry in entries:
            entry = require_mapping(entry, f"target metadata for {panel}")
            gene = entry.get("gene")
            scope = entry.get("scope")
            if not isinstance(gene, str) or scope not in {"high_copy_candidate", "deferred_or_unclassified"}:
                raise ValueError(f"Invalid target metadata for {panel}")
            if gene in scopes.setdefault(panel, {}):
                raise ValueError(f"Conflicting target metadata for {panel}/{gene}")
            scopes[panel][gene] = entry
    return scopes
def roles(protocol, schedule):
    baseline = protocol.get("baseline_role")
    versions = {item.get("version") for item in schedule if isinstance(item, dict)}
    candidate = versions - {baseline}
    if not isinstance(baseline, str) or len(candidate) != 1 or not all(isinstance(version, str) for version in versions):
        raise ValueError("Protocol must specify exactly one baseline role and one candidate role")
    return baseline, candidate.pop()
def analyze(execution):
    execution = Path(execution).resolve()
    provenance = load_json(execution / "provenance.json")
    protocol = require_mapping(provenance.get("protocol"), "frozen protocol")
    receipt = require_mapping(require_mapping(provenance.get("frozen_receipts"), "frozen receipts").get("protocol"), "protocol receipt")
    receipt_path = Path(receipt.get("path", "")).resolve()
    if receipt_path != (execution / "receipts" / "protocol.json").resolve() or sha256_file(receipt_path) != receipt.get("sha256") or load_json(receipt_path) != protocol:
        raise ValueError("Frozen protocol receipt is invalid")
    scopes = scope_map(protocol)
    schedule = provenance.get("schedule")
    if not isinstance(schedule, list) or not schedule:
        raise ValueError("Execution schedule is missing")
    baseline_role, candidate_role = roles(protocol, schedule)
    expected = {item.get("invocation_id"): item for item in schedule if isinstance(item, dict)}
    if len(expected) != len(schedule) or None in expected:
        raise ValueError("Execution schedule has duplicate or invalid identities")
    result_paths = {path.stem: path for path in (execution / "results").glob("*.json")}
    if set(result_paths) != set(expected):
        raise ValueError("Results do not exactly match the frozen schedule")
    cells = collections.defaultdict(lambda: collections.defaultdict(dict))
    for invocation_id, invocation in expected.items():
        if invocation.get("version") not in {baseline_role, candidate_role}:
            raise ValueError(f"Unexpected version for {invocation_id}")
        result = load_json(result_paths[invocation_id])
        validate_result(result, invocation)
        cell = invocation.get("cell")
        if not isinstance(cell, str) or invocation["version"] in cells[cell][invocation.get("pair_index")]:
            raise ValueError(f"Invalid pair identity for {invocation_id}")
        cells[cell][invocation["pair_index"]][invocation["version"]] = result
    analyzed_cells = []
    for cell, pairs in sorted(cells.items()):
        pair_rows = []
        panel = next(iter(pairs.values()))[baseline_role]["signature"]["invocation"]["panel"]
        if panel not in scopes:
            raise ValueError(f"Panel metadata missing for {panel}")
        for pair_index, versions in sorted(pairs.items()):
            if set(versions) != {baseline_role, candidate_role}:
                raise ValueError(f"Incomplete pair for {cell}/{pair_index}")
            if any(versions[version]["signature"]["invocation"].get("panel") != panel for version in (baseline_role, candidate_role)):
                raise ValueError(f"Mixed panels in {cell}/{pair_index}")
            baseline_metrics, candidate_metrics = metrics(versions[baseline_role]), metrics(versions[candidate_role])
            baseline_products = result_products(versions[baseline_role], panel, scopes[panel])
            candidate_products = result_products(versions[candidate_role], panel, scopes[panel])
            baseline_membership, baseline_duplicates = membership(baseline_products)
            candidate_membership, candidate_duplicates = membership(candidate_products)
            lost = [baseline_membership[key] for key in sorted(set(baseline_membership) - set(candidate_membership))]
            gained = [candidate_membership[key] for key in sorted(set(candidate_membership) - set(baseline_membership))]
            pair_rows.append({
                "pair_index": pair_index,
                "baseline_role": baseline_role,
                "candidate_role": candidate_role,
                "baseline": baseline_metrics,
                "candidate": candidate_metrics,
                "counts_identical": baseline_metrics["counts"] == candidate_metrics["counts"],
                "wall_ratio_candidate_over_baseline": candidate_metrics["wall_time_s"] / baseline_metrics["wall_time_s"],
                "rss_ratio_candidate_over_baseline": candidate_metrics["peak_rss_bytes"] / baseline_metrics["peak_rss_bytes"] if baseline_metrics["peak_rss_bytes"] else None,
                "lost_products": lost,
                "gained_products": gained,
                "baseline_duplicate_products": baseline_duplicates,
                "candidate_duplicate_products": candidate_duplicates,
            })
        version_products = {version: [result_products(pair[version], panel, scopes[panel]) for pair in pairs.values()] for version in (baseline_role, candidate_role)}
        stability = {version: {"sequence_membership_stable": len({tuple(sorted(membership(products)[0])) for products in outputs}) == 1, "classification_stable": len({tuple(sorted((value["gene"], value["sha256"], value["reference_match"].get("status")) for value in products)) for products in outputs}) == 1} for version, outputs in version_products.items()}
        losses = [product for pair in pair_rows for product in pair["lost_products"]]
        count_sums = {version: {field: sum(pair[version]["counts"][field] for pair in pair_rows) for field in COUNT_FIELDS} for version in ("baseline", "candidate")}
        analyzed_cells.append({"cell": cell, "pairs": pair_rows, "repetitions": len(pair_rows), "baseline_wall_time_s": summary([pair["baseline"]["wall_time_s"] for pair in pair_rows]), "candidate_wall_time_s": summary([pair["candidate"]["wall_time_s"] for pair in pair_rows]), "baseline_peak_rss_bytes": summary([pair["baseline"]["peak_rss_bytes"] for pair in pair_rows]), "candidate_peak_rss_bytes": summary([pair["candidate"]["peak_rss_bytes"] for pair in pair_rows]), "aggregate_count_parity": {"all_pairs_identical": all(pair["counts_identical"] for pair in pair_rows), "baseline": count_sums["baseline"], "candidate": count_sums["candidate"], "identical": count_sums["baseline"] == count_sums["candidate"]}, "same_version_stability": {"baseline": stability[baseline_role], "candidate": stability[candidate_role]}, "losses_by_scope": {scope: [product for product in losses if product["scope"] == scope] for scope in ("high_copy_candidate", "deferred_or_unclassified")}, "all_unfiltered_losses": losses})
    total_baseline = sum(cell["baseline_wall_time_s"]["median"] for cell in analyzed_cells)
    total_candidate = sum(cell["candidate_wall_time_s"]["median"] for cell in analyzed_cells)
    return {"schema_version": 1, "execution": str(execution), "baseline_role": baseline_role, "candidate_role": candidate_role, "protocol_priority": protocol.get("priority"), "review_criteria": protocol.get("review_criteria"), "cells": analyzed_cells, "sum_per_cell_median_wall_seconds": {"baseline": total_baseline, "candidate": total_candidate, "candidate_over_baseline_ratio": total_candidate / total_baseline}, "interpretation": "Descriptive comparison only: membership changes are neither automatic false positives nor automatic true gains, and this analysis does not authorize release or claim no regression."}
def write_report(path, analysis):
    lines = ["# High-copy comparison analysis", "", analysis["interpretation"], "", "| Cell | Repetitions | Baseline median wall s | Candidate median wall s | Wall ratio | Baseline median RSS | Candidate median RSS | RSS ratio | Count parity |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"]
    for cell in analysis["cells"]:
        pairs = cell["pairs"]
        lines.append(f"| {cell['cell']} | {cell['repetitions']} | {cell['baseline_wall_time_s']['median']:.3f} | {cell['candidate_wall_time_s']['median']:.3f} | {statistics.median(pair['wall_ratio_candidate_over_baseline'] for pair in pairs):.3f} | {cell['baseline_peak_rss_bytes']['median']} | {cell['candidate_peak_rss_bytes']['median']} | {statistics.median(pair['rss_ratio_candidate_over_baseline'] for pair in pairs if pair['rss_ratio_candidate_over_baseline'] is not None):.3f} | {all(pair['counts_identical'] for pair in pairs)} |")
    totals = analysis["sum_per_cell_median_wall_seconds"]
    lines.extend(["", f"Sum of per-cell medians: baseline {totals['baseline']:.3f}s; candidate {totals['candidate']:.3f}s; ratio {totals['candidate_over_baseline_ratio']:.3f}.", "", "## Losses", ""])
    for cell in analysis["cells"]:
        for scope, products in cell["losses_by_scope"].items():
            if products:
                lines.append(f"- `{cell['cell']}` {scope}: " + ", ".join(f"{product['gene']} {product['length']}bp sha256={product['sha256']} median={product['kmer_count_median']} min={product['kmer_count_min']} status={product['reference_match']['status']}" for product in products))
    lines.extend(["", "## Gains", ""])
    for cell in analysis["cells"]:
        for pair in cell["pairs"]:
            if pair["gained_products"]:
                lines.append(f"- `{cell['cell']}` pair {pair['pair_index']}: " + ", ".join(f"{product['gene']} {product['length']}bp sha256={product['sha256']} median={product['kmer_count_median']} min={product['kmer_count_min']} status={product['reference_match']['status']}" for product in pair["gained_products"]))
    Path(path).write_text("\n".join(lines) + "\n")
def main():
    parser = argparse.ArgumentParser(description="Analyze completed high-copy release comparison results")
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise SystemExit(f"Output already exists: {arguments.output}")
    analysis = analyze(arguments.execution)
    arguments.output.mkdir(parents=True)
    with (arguments.output / "analysis.json").open("x") as output_stream:
        json.dump(analysis, output_stream, indent=2, sort_keys=True)
        output_stream.write("\n")
    write_report(arguments.output / "report.md", analysis)
if __name__ == "__main__":
    main()

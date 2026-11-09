#!/usr/bin/env python3

import argparse
import collections
import hashlib
import json
import statistics
from pathlib import Path


PRIMARY_DEPTH = 1_000_000
VERSIONS = ("baseline", "candidate")
COUNT_FIELDS = (
    "n_reads_read",
    "n_bases_read",
    "n_subreads_ingested",
    "n_bases_ingested",
    "n_kmers",
)


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


def write_json(path, value):
    with Path(path).open("x") as output_stream:
        json.dump(value, output_stream, indent=2, sort_keys=True)
        output_stream.write("\n")


def product_multiset(result):
    products = collections.Counter()
    for gene in result.get("genes", []):
        for product in gene.get("products", []):
            products[(gene.get("gene"), product.get("length"), product.get("sha256"))] += 1
    return products


def indexed_product_signature(result):
    products = []
    for gene in result.get("genes", []):
        for product in gene.get("products", []):
            products.append(
                (
                    gene.get("gene"),
                    product.get("product_index"),
                    product.get("length"),
                    product.get("sha256"),
                )
            )
    return sorted(products)


def classification_multiset(result):
    classifications = collections.Counter()
    for gene in result.get("genes", []):
        for product in gene.get("products", []):
            status = product.get("reference_match", {}).get("status")
            classifications[(gene.get("gene"), product.get("length"), product.get("sha256"), status)] += 1
    return classifications


def expand_counter(counter):
    values = []
    for value, count in sorted(counter.items()):
        values.extend([list(value)] * count)
    return values


def subset_classifications(result, selected_products):
    remaining = selected_products.copy()
    statuses = collections.Counter()
    for gene in result.get("genes", []):
        for product in gene.get("products", []):
            key = (gene["gene"], product["length"], product["sha256"])
            if remaining[key]:
                statuses[product["reference_match"]["status"]] += 1
                remaining[key] -= 1
    if any(remaining.values()):
        raise ValueError("Classification subset is not present in its source result")
    return dict(sorted(statuses.items()))


def status_summary(result):
    if result.get("status") != "complete":
        return {
            "observed_result_status": result.get("status"),
            "failure": result.get("failure"),
            "error": result.get("error") or result.get("classification_error"),
            "genes": None,
            "gene_status_counts": None,
            "classification_status_counts": None,
        }
    genes = []
    gene_status_counts = collections.Counter()
    classification_status_counts = collections.Counter()
    for gene in result.get("genes", []):
        recovered = bool(gene.get("recovered"))
        gene_status_counts["recovered" if recovered else "failed"] += 1
        genes.append(
            {
                "gene": gene.get("gene"),
                "recovered": recovered,
                "n_products": gene.get("n_products"),
                "failure_reason": gene.get("failure_reason"),
            }
        )
        for product in gene.get("products", []):
            classification_status_counts[product.get("reference_match", {}).get("status")] += 1
    return {
        "observed_result_status": "complete",
        "genes": genes,
        "gene_status_counts": dict(sorted(gene_status_counts.items())),
        "classification_status_counts": dict(sorted(classification_status_counts.items(), key=lambda item: str(item[0]))),
    }


def numeric_summary(values):
    usable = [value for value in values if isinstance(value, (int, float))]
    if not usable:
        return None
    return {
        "n": len(usable),
        "median": statistics.median(usable),
        "minimum": min(usable),
        "maximum": max(usable),
    }


def run_metrics(result):
    if result.get("status") != "complete":
        return None
    metrics = result.get("metrics")
    execution = result.get("execution")
    if not isinstance(metrics, dict) or not isinstance(execution, dict):
        return None
    gnu_time = execution.get("gnu_time") if isinstance(execution.get("gnu_time"), dict) else {}
    return {
        "wall_time_s": execution.get("wall_time_s"),
        "peak_rss_bytes": gnu_time.get("peak_rss_bytes"),
        "counts": {field: metrics.get(field) for field in COUNT_FIELDS},
    }


def compare_pair(baseline, candidate):
    comparison = {
        "baseline_status": baseline.get("status"),
        "candidate_status": candidate.get("status"),
        "both_complete": baseline.get("status") == candidate.get("status") == "complete",
    }
    baseline_metrics = run_metrics(baseline)
    candidate_metrics = run_metrics(candidate)
    if not comparison["both_complete"]:
        comparison.update(
            {
                "baseline_observation": status_summary(baseline),
                "candidate_observation": status_summary(candidate),
                "baseline_wall_time_s": baseline_metrics["wall_time_s"] if baseline_metrics else None,
                "candidate_wall_time_s": candidate_metrics["wall_time_s"] if candidate_metrics else None,
                "baseline_peak_rss_bytes": baseline_metrics["peak_rss_bytes"] if baseline_metrics else None,
                "candidate_peak_rss_bytes": candidate_metrics["peak_rss_bytes"] if candidate_metrics else None,
            }
        )
        return comparison
    baseline_products = product_multiset(baseline)
    candidate_products = product_multiset(candidate)
    baseline_classifications = classification_multiset(baseline)
    candidate_classifications = classification_multiset(candidate)
    baseline_wall = baseline_metrics["wall_time_s"]
    candidate_wall = candidate_metrics["wall_time_s"]
    comparison.update(
        {
            "counts_identical": all(
                baseline_metrics["counts"][field] == candidate_metrics["counts"][field]
                for field in COUNT_FIELDS
            ),
            "count_values": {
                "baseline": baseline_metrics["counts"],
                "candidate": candidate_metrics["counts"],
                "n_kmers_label": "accepted k-mer occurrences",
            },
            "products_identical_ignoring_index": baseline_products == candidate_products,
            "product_index_ordering_identical": indexed_product_signature(baseline) == indexed_product_signature(candidate),
            "classifications_identical_ignoring_index": baseline_classifications == candidate_classifications,
            "retained_products": expand_counter(baseline_products & candidate_products),
            "lost_products": expand_counter(baseline_products - candidate_products),
            "gained_products": expand_counter(candidate_products - baseline_products),
            "lost_classifications": subset_classifications(baseline, baseline_products - candidate_products),
            "gained_classifications": subset_classifications(candidate, candidate_products - baseline_products),
            "retained_baseline_classifications": subset_classifications(baseline, baseline_products & candidate_products),
            "retained_candidate_classifications": subset_classifications(candidate, baseline_products & candidate_products),
            "baseline_status_summary": status_summary(baseline),
            "candidate_status_summary": status_summary(candidate),
            "baseline_wall_time_s": baseline_wall,
            "candidate_wall_time_s": candidate_wall,
            "candidate_speedup_ratio": baseline_wall / candidate_wall if baseline_wall and candidate_wall else None,
            "baseline_peak_rss_bytes": baseline_metrics["peak_rss_bytes"],
            "candidate_peak_rss_bytes": candidate_metrics["peak_rss_bytes"],
        }
    )
    return comparison


def repeat_stability(runs):
    completed = [result for result in runs if result.get("status") == "complete"]
    statuses = [result.get("status") for result in runs]
    if len(completed) != len(runs):
        return {
            "assessable": False,
            "completed_pairs": len(completed),
            "observed_statuses": statuses,
        }
    product_signatures = [product_multiset(result) for result in completed]
    classification_signatures = [classification_multiset(result) for result in completed]
    return {
        "assessable": True,
        "completed_pairs": len(completed),
        "sequence_multiset_stable": all(signature == product_signatures[0] for signature in product_signatures[1:]),
        "classification_multiset_stable": all(signature == classification_signatures[0] for signature in classification_signatures[1:]),
    }


def cell_analysis(cell, pairs, primary):
    pair_indexes = sorted(pairs)
    pair_comparisons = [
        {"pair_index": pair_index, **compare_pair(pairs[pair_index]["baseline"], pairs[pair_index]["candidate"])}
        for pair_index in pair_indexes
    ]
    representative_pair_index = pair_indexes[0]
    representative = pair_comparisons[0]
    actual_records = {
        version: pairs[representative_pair_index][version]["signature"]["input_subset"]["actual_records"]
        for version in VERSIONS
    }
    result = {
        "cell": cell,
        "depth": pairs[representative_pair_index]["baseline"]["signature"]["input_subset"]["requested_records"],
        "actual_records": actual_records,
        "pairs": pair_comparisons,
        "representative_pair_index": representative_pair_index,
        "representative_pair": representative,
    }
    if primary:
        for version in VERSIONS:
            version_runs = [pairs[pair_index][version] for pair_index in pair_indexes]
            metrics = [run_metrics(run) for run in version_runs]
            result[f"{version}_wall_time_s"] = numeric_summary([metric["wall_time_s"] for metric in metrics if metric])
            result[f"{version}_peak_rss_bytes"] = numeric_summary([metric["peak_rss_bytes"] for metric in metrics if metric])
            result[f"{version}_repeat_stability"] = repeat_stability(version_runs)
        paired_ratios = [pair["candidate_speedup_ratio"] for pair in pair_comparisons if pair.get("both_complete")]
        result["paired_speedup_ratio"] = numeric_summary(paired_ratios)
    return result


def markdown_value(value):
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def markdown_range(summary):
    if not summary:
        return "—"
    return f"{summary['median']:.3f} [{summary['minimum']:.3f}, {summary['maximum']:.3f}]"


def markdown_product_counts(pair):
    if not pair.get("both_complete"):
        return "—", "—", "—"
    return (
        str(len(pair["retained_products"])),
        str(len(pair["lost_products"])),
        str(len(pair["gained_products"])),
    )


def write_markdown(path, analysis):
    lines = ["# Released 3.1.0 vs Candidate Analysis", "", "This is calibration/regression analysis, not held-out biological validation.", ""]
    lines.extend(["## Primary 1M Timing", "", "Median [min, max] wall seconds and RSS are across three paired runs. `n_kmers` means accepted k-mer occurrences.", "", "| Cell | Actual records | Baseline wall s | Candidate wall s | Baseline RSS bytes | Candidate RSS bytes | Paired speedup |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
    for cell in analysis["primary_cells"]:
        actual_records = cell["actual_records"]["baseline"]
        lines.append(
            f"| {cell['cell']} | {actual_records} | {markdown_range(cell['baseline_wall_time_s'])} | {markdown_range(cell['candidate_wall_time_s'])} | {markdown_range(cell['baseline_peak_rss_bytes'])} | {markdown_range(cell['candidate_peak_rss_bytes'])} | {markdown_range(cell['paired_speedup_ratio'])} |"
        )
    lines.extend(["", "| Version | Sum of available per-cell median wall seconds | Cells with medians | Missing cells |", "| --- | ---: | ---: | ---: |"])
    for version, totals in analysis["primary_median_sums"].items():
        lines.append(f"| {version} | {totals['sum_per_cell_medians_seconds']:.3f} | {totals['cells_with_medians']} | {totals['cells_missing_medians']} |")
    lines.extend(["", "## Deeper Descriptive Runs", "", "Each 2M/4M/8M cell has one paired run; values are descriptive, not replicated timing estimates.", "", "| Cell | Actual records | Baseline status/time/RSS | Candidate status/time/RSS | Retained | Lost | Gained |", "| --- | ---: | --- | --- | ---: | ---: | ---: |"])
    for cell in analysis["deeper_cells"]:
        pair = cell["representative_pair"]
        baseline = f"{pair['baseline_status']}/{markdown_value(pair.get('baseline_wall_time_s'))}/{markdown_value(pair.get('baseline_peak_rss_bytes'))}"
        candidate = f"{pair['candidate_status']}/{markdown_value(pair.get('candidate_wall_time_s'))}/{markdown_value(pair.get('candidate_peak_rss_bytes'))}"
        retained, lost, gained = markdown_product_counts(pair)
        lines.append(f"| {cell['cell']} | {cell['actual_records']['baseline']} | {baseline} | {candidate} | {retained} | {lost} | {gained} |")
    lines.extend(["", "## Per-Cell Product and Status Changes", "", "Product comparisons use the `(gene, length, SHA-256)` multiset, so biological gains/losses do not depend on output index. The JSON includes per-gene recovery/failure reasons and classification status counts.", "", "| Cell | Baseline status | Candidate status | Retained | Lost | Gained | Product multiset identical | Index ordering identical |", "| --- | --- | --- | ---: | ---: | ---: | --- | --- |"])
    for cell in [*analysis["primary_cells"], *analysis["deeper_cells"]]:
        pair = cell["representative_pair"]
        retained, lost, gained = markdown_product_counts(pair)
        lines.append(f"| {cell['cell']} | {pair['baseline_status']} | {pair['candidate_status']} | {retained} | {lost} | {gained} | {markdown_value(pair.get('products_identical_ignoring_index'))} | {markdown_value(pair.get('product_index_ordering_identical'))} |")
    lines.extend(["", "## Primary Repeat Stability", "", "Sequence and classification stability compare product multisets within each version; failed observations are not interpreted as zero-product runs.", "", "| Cell | Version | Assessable | Sequence stable | Classification stable | Observed statuses |", "| --- | --- | --- | --- | --- | --- |"])
    for cell in analysis["primary_cells"]:
        for version in VERSIONS:
            stability = cell[f"{version}_repeat_stability"]
            lines.append(f"| {cell['cell']} | {version} | {stability['assessable']} | {markdown_value(stability.get('sequence_multiset_stable'))} | {markdown_value(stability.get('classification_multiset_stable'))} | {', '.join(stability.get('observed_statuses', []))} |")
    lines.extend(["", "## Observed Failures", "", "Failures remain failures; no missing result is converted to a zero-product observation.", "", "| Invocation | Status | Failure | Error |", "| --- | --- | --- | --- |"])
    for failure in analysis["failures"]:
        lines.append(f"| {failure['invocation_id']} | {failure['status']} | {failure.get('failure') or ''} | {failure.get('error') or ''} |")
    if not analysis["failures"]:
        lines.append("| None | — | — | — |")
    lines.extend(["", "## Representative Totals", "", "Totals use primary pair 1 once per primary cell plus each deeper cell once; three primary repetitions are not tripled.", "", "| Version | Complete cells | Products | Classification counts |", "| --- | ---: | ---: | --- |"])
    for version, totals in analysis["representative_totals"].items():
        lines.append(f"| {version} | {totals['complete_cells']} | {totals['products']} | {json.dumps(totals['classification_status_counts'], sort_keys=True)} |")
    lines.extend(["", "## Primary-Only Product Totals", "", "One representative run per primary cell; no deeper cells or repeated runs included.", "", "| Version | Complete cells | Products | Classification counts |", "| --- | ---: | ---: | --- |"])
    for version, totals in analysis["primary_representative_totals"].items():
        lines.append(f"| {version} | {totals['complete_cells']} | {totals['products']} | {json.dumps(totals['classification_status_counts'], sort_keys=True)} |")
    Path(path).write_text("\n".join(lines) + "\n")


def analyze(normalized_root, source_execution):
    normalized_root = Path(normalized_root)
    source_execution = Path(source_execution)
    schedule_document = load_json(source_execution / "schedule.json")
    source_comparison_path = source_execution / "comparison.json"
    source_comparison = load_json(source_comparison_path)
    schedule = schedule_document.get("schedule")
    if not isinstance(schedule, list):
        raise ValueError("Source schedule is missing")
    if source_comparison.get("invocations") != len(schedule):
        raise ValueError("Source comparison is unfinished")
    normalized_comparison_path = normalized_root / "comparison.json"
    normalized_comparison = load_json(normalized_comparison_path)
    if normalized_comparison.get("source_comparison_sha256") != sha256_file(source_comparison_path):
        raise ValueError("Normalized comparison does not attest the source comparison")
    results_directory = normalized_root / "results"
    results = {}
    for invocation in schedule:
        invocation_id = invocation["invocation_id"]
        result_path = results_directory / f"{invocation_id}.json"
        results[invocation_id] = load_json(result_path)
        if results[invocation_id].get("status") not in {"complete", "failed"}:
            raise ValueError(f"Normalized result is not final: {invocation_id}")
    unexpected = {path.name for path in results_directory.glob("*.json")} - {f"{entry['invocation_id']}.json" for entry in schedule}
    if unexpected:
        raise ValueError(f"Normalized results contain unexpected files: {sorted(unexpected)}")
    cells = collections.defaultdict(lambda: collections.defaultdict(dict))
    failures = []
    for invocation in schedule:
        result = results[invocation["invocation_id"]]
        cells[invocation["cell"]][invocation["pair_index"]][invocation["version"]] = result
        if result.get("status") != "complete":
            failures.append({"invocation_id": invocation["invocation_id"], "status": result.get("status"), "failure": result.get("failure"), "error": result.get("error") or result.get("classification_error")})
    for cell, pairs in cells.items():
        for pair_index, versions in pairs.items():
            if set(versions) != set(VERSIONS):
                raise ValueError(f"Cell {cell} pair {pair_index} lacks a version")
    primary_cells = []
    deeper_cells = []
    for cell, pairs in sorted(cells.items()):
        depth = next(iter(next(iter(pairs.values())).values()))["signature"]["input_subset"]["requested_records"]
        analysis = cell_analysis(cell, pairs, depth == PRIMARY_DEPTH)
        (primary_cells if depth == PRIMARY_DEPTH else deeper_cells).append(analysis)
    primary_median_sums = {}
    for version in VERSIONS:
        summaries = [cell[f"{version}_wall_time_s"] for cell in primary_cells]
        primary_median_sums[version] = {
            "sum_per_cell_medians_seconds": sum(summary["median"] for summary in summaries if summary),
            "cells_with_medians": sum(summary is not None for summary in summaries),
            "cells_missing_medians": sum(summary is None for summary in summaries),
        }
    representative_totals = {version: {"complete_cells": 0, "products": 0, "classification_status_counts": {}} for version in VERSIONS}
    primary_representative_totals = {version: {"complete_cells": 0, "products": 0, "classification_status_counts": {}} for version in VERSIONS}
    for cell in [*primary_cells, *deeper_cells]:
        pair_index = cell["representative_pair_index"]
        for version in VERSIONS:
            result = cells[cell["cell"]][pair_index][version]
            if result.get("status") != "complete":
                continue
            representative_totals[version]["complete_cells"] += 1
            representative_totals[version]["products"] += sum(len(gene.get("products", [])) for gene in result.get("genes", []))
            counts = collections.Counter(representative_totals[version]["classification_status_counts"])
            for classification in expand_counter(classification_multiset(result)):
                counts[classification[3]] += 1
            representative_totals[version]["classification_status_counts"] = dict(sorted(counts.items(), key=lambda item: str(item[0])))
            if cell["depth"] == PRIMARY_DEPTH:
                primary_representative_totals[version]["complete_cells"] += 1
                primary_representative_totals[version]["products"] += sum(len(gene.get("products", [])) for gene in result.get("genes", []))
                primary_counts = collections.Counter(primary_representative_totals[version]["classification_status_counts"])
                for classification in expand_counter(classification_multiset(result)):
                    primary_counts[classification[3]] += 1
                primary_representative_totals[version]["classification_status_counts"] = dict(sorted(primary_counts.items(), key=lambda item: str(item[0])))
    return {
        "schema_version": 1,
        "source_execution": str(source_execution.resolve()),
        "source_schedule_sha256": sha256_file(source_execution / "schedule.json"),
        "normalized_root": str(normalized_root.resolve()),
        "normalized_comparison_sha256": sha256_file(normalized_comparison_path),
        "primary_cells": primary_cells,
        "deeper_cells": deeper_cells,
        "primary_median_sums": primary_median_sums,
        "representative_totals": representative_totals,
        "primary_representative_totals": primary_representative_totals,
        "failures": failures,
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze finalized normalized release-comparison results")
    parser.add_argument("--normalized-root", type=Path, required=True)
    parser.add_argument("--source-execution", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output_dir.exists():
        raise SystemExit(f"Analysis output already exists: {arguments.output_dir}")
    analysis = analyze(arguments.normalized_root, arguments.source_execution)
    arguments.output_dir.mkdir(parents=True)
    write_json(arguments.output_dir / "analysis.json", analysis)
    write_markdown(arguments.output_dir / "analysis.md", analysis)


if __name__ == "__main__":
    main()

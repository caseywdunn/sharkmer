#!/usr/bin/env python3
import argparse
import json
import statistics
from pathlib import Path


HIGH_COPY = "high_copy_candidate"


def load_json(path):
    value = json.loads(Path(path).read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def product_key(product):
    gene, digest = product.get("gene"), product.get("sha256")
    if not isinstance(gene, str) or not isinstance(digest, str):
        raise ValueError("Membership product lacks gene or SHA-256")
    return gene, digest


def observed_products(comparison):
    values = comparison.get("exact_retained", []) + comparison.get("gained", [])
    if not all(isinstance(value, dict) for value in values):
        raise ValueError("Membership comparison has malformed products")
    products = {}
    for product in values:
        key = product_key(product)
        if key in products:
            raise ValueError("Membership comparison duplicates an exact product")
        products[key] = product
    return products


def high_copy_products(comparison):
    return {
        key: product
        for key, product in observed_products(comparison).items()
        if product.get("scope") == HIGH_COPY
    }


def median(values):
    if not values:
        raise ValueError("Cannot calculate an empty median")
    return statistics.median(values)


def percent_change(value, reference):
    if reference == 0:
        return None
    return (value - reference) * 100.0 / reference


def representative_pair(cell):
    pairs = cell.get("pairs")
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("Cell lacks paired measurements")
    return pairs[0]


def same_k_differences(cell):
    differences = []
    for pair in cell["pairs"]:
        comparison = pair["same_k_membership"]
        lost = comparison.get("lost", [])
        gained = comparison.get("gained", [])
        if lost or gained:
            differences.append({"pair_index": pair["pair_index"], "lost": lost, "gained": gained})
    return differences


def its2_diagnostics(cell):
    evidence = []
    for pair in cell["pairs"]:
        diagnostics = pair["threshold_diagnostics"]["candidate"]
        genes = diagnostics.get("genes", []) if isinstance(diagnostics, dict) else []
        matching = [
            gene
            for gene in genes
            if isinstance(gene, dict)
            and isinstance(gene.get("gene_name"), str)
            and gene["gene_name"].endswith("_ITS_2")
        ]
        evidence.append({"pair_index": pair["pair_index"], "matches": matching})
    return evidence


def summarize(analysis):
    cells = analysis.get("cells")
    options = analysis.get("confirmation_options_ranked")
    if not isinstance(cells, list) or not isinstance(options, list):
        raise ValueError("Analysis lacks cells or confirmation options")
    by_k = {}
    for cell in cells:
        kmer_length = cell.get("kmer_length")
        if type(kmer_length) is not int or not isinstance(cell.get("cell"), str):
            raise ValueError("Analysis cell has invalid identity")
        by_k.setdefault(kmer_length, []).append(cell)
    if 19 not in by_k:
        raise ValueError("Analysis lacks k19 anchor")
    k19_cells = {cell["cell"]: cell for cell in by_k[19]}
    per_k = {}
    for kmer_length, k_cells in sorted(by_k.items()):
        if set(cell["cell"] for cell in k_cells) != set(k19_cells):
            raise ValueError(f"k{kmer_length} cells differ from k19")
        totals = {"baseline": 0, "candidate": 0}
        resource_rows = []
        parity_differences = []
        cross_k_current_differences = []
        stable = True
        count_parity = True
        diagnostics = []
        for cell in sorted(k_cells, key=lambda value: value["cell"]):
            pair = representative_pair(cell)
            totals["baseline"] += len(high_copy_products(pair["baseline_vs_k19_released"]))
            totals["candidate"] += len(high_copy_products(pair["candidate_vs_k19_released"]))
            differences = same_k_differences(cell)
            parity_differences.extend({"cell": cell["cell"], **difference} for difference in differences)
            current_comparison = pair["candidate_vs_k19_candidate"]
            if current_comparison.get("lost") or current_comparison.get("gained"):
                cross_k_current_differences.append({
                    "cell": cell["cell"],
                    "lost": current_comparison["lost"],
                    "gained": current_comparison["gained"],
                    "changed_endpoint_or_sequence_unresolved": current_comparison["changed_endpoint_or_sequence_unresolved"],
                })
            stable = stable and all(cell.get("repeat_sequence_and_classification_stable", {}).values())
            count_parity = count_parity and all(pair_value.get("same_k_count_parity") is True for pair_value in cell["pairs"])
            k19_pair = representative_pair(k19_cells[cell["cell"]])
            resource_rows.append({
                "cell": cell["cell"],
                "candidate_wall_time_s": cell["median_wall_time_s"]["candidate"],
                "candidate_wall_vs_own_k19_percent": percent_change(
                    cell["median_wall_time_s"]["candidate"],
                    k19_cells[cell["cell"]]["median_wall_time_s"]["candidate"],
                ),
                "candidate_peak_rss_bytes": cell["median_peak_rss_bytes"]["candidate"],
                "candidate_rss_vs_own_k19_percent": percent_change(
                    cell["median_peak_rss_bytes"]["candidate"],
                    k19_cells[cell["cell"]]["median_peak_rss_bytes"]["candidate"],
                ),
                "candidate_wall_vs_baseline_same_k_percent": percent_change(
                    cell["median_wall_time_s"]["candidate"],
                    cell["median_wall_time_s"]["baseline"],
                ),
                "candidate_rss_vs_baseline_same_k_percent": percent_change(
                    cell["median_peak_rss_bytes"]["candidate"],
                    cell["median_peak_rss_bytes"]["baseline"],
                ),
                "representative_same_k_counts": {
                    "baseline": pair["baseline"]["counts"],
                    "candidate": pair["candidate"]["counts"],
                },
                "k19_pair_index_used_for_anchor": k19_pair["pair_index"],
            })
            if cell["cell"].startswith("insecta/SRR27962769/"):
                diagnostics.append({"cell": cell["cell"], "pairs": its2_diagnostics(cell)})
        per_k[str(kmer_length)] = {
            "high_copy_exact_products_total": totals,
            "same_k_sequence_parity": not parity_differences,
            "same_k_sequence_differences": parity_differences,
            "candidate_cross_k19_current_sequence_differences": cross_k_current_differences,
            "same_k_count_parity_every_pair": count_parity,
            "repeat_sequence_and_classification_stable": stable,
            "candidate_resources_by_cell": resource_rows,
            "candidate_wall_sum_cell_medians_s": sum(row["candidate_wall_time_s"] for row in resource_rows),
            "baseline_wall_sum_cell_medians_s": sum(cell["median_wall_time_s"]["baseline"] for cell in k_cells),
            "candidate_rss_median_cell_medians_bytes": median(row["candidate_peak_rss_bytes"] for row in resource_rows),
            "baseline_rss_median_cell_medians_bytes": median(cell["median_peak_rss_bytes"]["baseline"] for cell in k_cells),
            "candidate_its2_raw_threshold_diagnostics": diagnostics,
        }
    options_by_k = {option.get("kmer_length"): option for option in options if isinstance(option, dict)}
    selection = {}
    for kmer_length in sorted(by_k):
        if kmer_length == 19:
            continue
        option = options_by_k.get(kmer_length)
        if option is None:
            raise ValueError(f"Analysis lacks confirmation option for k{kmer_length}")
        selection[str(kmer_length)] = {
            "eligible": option.get("eligible"),
            "exact_restored_known_losses_every_replicate": option.get("exact_restored_known_losses_every_replicate"),
            "exact_its2_ak281180_restored_every_replicate": option.get("exact_its2_ak281180_restored_every_replicate"),
            "confirmed_k19_candidate_sequences_lost": option.get("confirmed_k19_candidate_sequences_lost"),
            "k19_released_high_copy_exact_retained_every_replicate": option.get("k19_released_high_copy_exact_retained_every_replicate"),
        }
    expected_totals = {"k19_released_high_copy_products": 33, "k19_current_high_copy_products": 26}
    observed_k19 = per_k["19"]["high_copy_exact_products_total"]
    return {
        "schema_version": 1,
        "analysis_source": analysis.get("execution"),
        "reference_totals_expected": expected_totals,
        "reference_totals_observed": {
            "k19_released_high_copy_products": observed_k19["baseline"],
            "k19_current_high_copy_products": observed_k19["candidate"],
        },
        "reference_totals_match_expected": {
            "k19_released_high_copy_products": observed_k19["baseline"] == expected_totals["k19_released_high_copy_products"],
            "k19_current_high_copy_products": observed_k19["candidate"] == expected_totals["k19_current_high_copy_products"],
        },
        "per_k": per_k,
        "selection_evidence": selection,
        "caveats": [
            "Exact gene-plus-sequence membership is the recovery criterion; a same-gene sequence change remains unresolved rather than a rescue.",
            "Threshold, node, SCC, and seed-adjacent diagnostics are observations, not proof of a single cause for a retained or missing product.",
            "The frozen AK281180 probe found no same-strand repeated 18-, 22-, 26-, or 30-mers in that reference sequence alone; it does not exclude read-graph SCCs, other reads, errors, or orientation effects.",
            "No product absence establishes biological falsehood, and no gain establishes biological truth; this is a fixed-input calibration summary, not release approval.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description="Summarize approved cross-k sweep analysis")
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise SystemExit(f"Output already exists: {arguments.output}")
    summary = summarize(load_json(arguments.analysis))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(arguments.output, summary)
    print(arguments.output)


if __name__ == "__main__":
    main()

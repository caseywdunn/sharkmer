#!/usr/bin/env python3
"""
Compare two sharkmer benchmark result files.

Reports performance changes, product stability, and flags regressions.

Usage:
    python benchmarks/compare.py benchmarks/results/baseline.yaml benchmarks/results/current.yaml
"""

import argparse
import sys
from pathlib import Path

import yaml


def load_results(path):
    with open(path) as f:
        return yaml.safe_load(f)


def results_by_sample(benchmark):
    """Index results by sample name for easy lookup."""
    return {r["sample"]: r for r in benchmark.get("results", [])}


def compare_sequences(products_a, products_b):
    """Compare product sequences between two runs."""
    genes_a = {p["gene"]: p for p in products_a}
    genes_b = {p["gene"]: p for p in products_b}

    all_genes = sorted(set(list(genes_a.keys()) + list(genes_b.keys())))
    changes = []

    for gene in all_genes:
        pa = genes_a.get(gene)
        pb = genes_b.get(gene)

        if pa is None and pb is not None:
            changes.append(f"  {gene}: NEW (+{pb['n_products']} products)")
        elif pa is not None and pb is None:
            changes.append(f"  {gene}: LOST (was {pa['n_products']} products)")
        elif pa is not None and pb is not None:
            if pa["n_products"] != pb["n_products"]:
                changes.append(
                    f"  {gene}: product count {pa['n_products']} → {pb['n_products']}"
                )
            elif pa.get("sequences") and pb.get("sequences"):
                if pa["sequences"] != pb["sequences"]:
                    # Check lengths first
                    if pa["lengths"] != pb["lengths"]:
                        changes.append(
                            f"  {gene}: lengths changed {pa['lengths']} → {pb['lengths']}"
                        )
                    else:
                        changes.append(
                            f"  {gene}: sequences differ (same lengths)"
                        )

    return changes


def compare(path_a, path_b):
    """Compare two benchmark files and print a report."""
    a = load_results(path_a)
    b = load_results(path_b)

    label_a = a.get("sharkmer_version", "A")
    label_b = b.get("sharkmer_version", "B")
    commit_a = a.get("git_commit", "?")
    commit_b = b.get("git_commit", "?")

    print(f"Comparing: {label_a} ({commit_a}) vs {label_b} ({commit_b})")
    print(f"  Baseline: {path_a}")
    print(f"  Current:  {path_b}")
    print()

    samples_a = results_by_sample(a)
    samples_b = results_by_sample(b)

    all_samples = sorted(set(list(samples_a.keys()) + list(samples_b.keys())))

    regressions = []
    improvements = []

    # Header
    print(f"{'Sample':<30} {'Time A':>8} {'Time B':>8} {'Delta':>8}  {'Products':>10}  Sequence changes")
    print("-" * 100)

    for sample in all_samples:
        ra = samples_a.get(sample)
        rb = samples_b.get(sample)

        if ra is None:
            print(f"{sample:<30} {'---':>8} {rb.get('wall_time_s', '?'):>8.1f}s {'NEW':>8}")
            continue
        if rb is None:
            print(f"{sample:<30} {ra.get('wall_time_s', '?'):>8.1f}s {'---':>8} {'MISSING':>8}")
            regressions.append(f"{sample}: missing from current run")
            continue

        # Handle failed runs
        if ra.get("status") == "failed" or rb.get("status") == "failed":
            status_a = ra.get("status", "ok")
            status_b = rb.get("status", "ok")
            print(f"{sample:<30} {status_a:>8} {status_b:>8}")
            if status_a != "failed" and status_b == "failed":
                regressions.append(f"{sample}: now failing")
            continue

        time_a = ra.get("wall_time_s", 0)
        time_b = rb.get("wall_time_s", 0)

        if time_a > 0:
            delta_pct = ((time_b - time_a) / time_a) * 100
            delta_str = f"{delta_pct:+.1f}%"
        else:
            delta_str = "---"

        products_a = ra.get("products", [])
        products_b = rb.get("products", [])
        n_genes_a = len(products_a)
        n_genes_b = len(products_b)
        product_str = f"{n_genes_a} → {n_genes_b}" if n_genes_a != n_genes_b else f"{n_genes_a}"

        # Check sequences
        seq_changes = compare_sequences(products_a, products_b)
        seq_str = f"{len(seq_changes)} changes" if seq_changes else "stable"

        # Flag regressions
        flag = ""
        if time_a > 0 and ((time_b - time_a) / time_a) > 0.10:
            flag = " ⚠"
            regressions.append(f"{sample}: wall time +{delta_pct:.1f}%")
        elif time_a > 0 and ((time_b - time_a) / time_a) < -0.10:
            improvements.append(f"{sample}: wall time {delta_pct:.1f}%")

        if seq_changes:
            flag += " ⚠"
            for change in seq_changes:
                regressions.append(f"{sample}: {change.strip()}")

        print(f"{sample:<30} {time_a:>7.1f}s {time_b:>7.1f}s {delta_str:>8}  {product_str:>10}  {seq_str}{flag}")

    print()

    if regressions:
        print("REGRESSIONS:")
        for r in regressions:
            print(f"  ⚠ {r}")
        print()

    if improvements:
        print("IMPROVEMENTS:")
        for i in improvements:
            print(f"  ✓ {i}")
        print()

    if not regressions and not improvements:
        print("No significant changes detected.")

    return len(regressions)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two sharkmer benchmark results"
    )
    parser.add_argument("baseline", help="Baseline benchmark YAML file")
    parser.add_argument("current", help="Current benchmark YAML file")
    args = parser.parse_args()

    n_regressions = compare(args.baseline, args.current)
    sys.exit(1 if n_regressions > 0 else 0)


if __name__ == "__main__":
    main()

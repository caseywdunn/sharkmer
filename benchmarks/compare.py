#!/usr/bin/env python3
"""
Compare two sharkmer benchmark result files.

Reports performance changes, product stability, and flags regressions.
Handles per-sample multi-read-count sweeps.

Usage:
    python benchmarks/compare.py benchmarks/results/baseline.yaml benchmarks/results/current.yaml
"""

import argparse
import hashlib
import sys
from pathlib import Path

import yaml


def load_results(path):
    with open(path) as f:
        return yaml.safe_load(f)


def results_by_key(benchmark):
    """Index results by (sample, max_reads) for easy lookup."""
    index = {}
    for r in benchmark.get("results", []):
        key = (r["sample"], r.get("max_reads", 0))
        index[key] = r
    return index


def seq_fingerprint(seq):
    """Short hash for display when sequences differ."""
    return hashlib.md5(seq.encode()).hexdigest()[:8]


def compare_sequences(products_a, products_b):
    """Compare product sequences between two runs.

    Returns (changes, details) where changes is a list of one-line summaries
    and details is a list of verbose diff information.
    """
    genes_a = {p["gene"]: p for p in products_a}
    genes_b = {p["gene"]: p for p in products_b}

    all_genes = sorted(set(list(genes_a.keys()) + list(genes_b.keys())))
    changes = []
    details = []

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
                    f"  {gene}: product count {pa['n_products']} -> {pb['n_products']}"
                )

            seqs_a = pa.get("sequences", [])
            seqs_b = pb.get("sequences", [])

            if seqs_a and seqs_b and seqs_a != seqs_b:
                if pa["lengths"] != pb["lengths"]:
                    changes.append(
                        f"  {gene}: lengths changed {pa['lengths']} -> {pb['lengths']}"
                    )
                else:
                    # Same lengths but different content — show fingerprints
                    n_diff = sum(1 for a, b in zip(seqs_a, seqs_b) if a != b)
                    changes.append(
                        f"  {gene}: {n_diff}/{len(seqs_a)} sequences differ (same lengths)"
                    )

                # Detailed per-product diff
                for i, (sa, sb) in enumerate(zip(seqs_a, seqs_b)):
                    if sa != sb:
                        details.append(
                            f"  {gene} product {i}: "
                            f"{len(sa)}bp {seq_fingerprint(sa)} -> "
                            f"{len(sb)}bp {seq_fingerprint(sb)}"
                        )

    return changes, details


def compare(path_a, path_b, verbose=False):
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

    index_a = results_by_key(a)
    index_b = results_by_key(b)

    all_keys = sorted(set(list(index_a.keys()) + list(index_b.keys())))

    regressions = []
    improvements = []
    all_details = []

    # Header
    print(
        f"{'Sample':<30} {'Reads':>8} {'Time A':>8} {'Time B':>8} "
        f"{'Delta':>8}  {'Products':>10}  Sequence changes"
    )
    print("-" * 110)

    for sample, max_reads in all_keys:
        ra = index_a.get((sample, max_reads))
        rb = index_b.get((sample, max_reads))
        k_reads = f"{max_reads // 1000}k" if max_reads else "?"

        if ra is None:
            time_b = rb.get("wall_time_s", "?")
            time_str = f"{time_b:.1f}s" if isinstance(time_b, (int, float)) else time_b
            print(f"{sample:<30} {k_reads:>8} {'---':>8} {time_str:>8} {'NEW':>8}")
            continue
        if rb is None:
            time_a = ra.get("wall_time_s", "?")
            time_str = f"{time_a:.1f}s" if isinstance(time_a, (int, float)) else time_a
            print(f"{sample:<30} {k_reads:>8} {time_str:>8} {'---':>8} {'MISSING':>8}")
            regressions.append(f"{sample} ({k_reads}): missing from current run")
            continue

        # Handle failed runs
        if ra.get("status") == "failed" or rb.get("status") == "failed":
            status_a = ra.get("status", "ok")
            status_b = rb.get("status", "ok")
            print(f"{sample:<30} {k_reads:>8} {status_a:>8} {status_b:>8}")
            if status_a != "failed" and status_b == "failed":
                regressions.append(f"{sample} ({k_reads}): now failing")
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
        product_str = f"{n_genes_a} -> {n_genes_b}" if n_genes_a != n_genes_b else f"{n_genes_a}"

        # Check sequences
        seq_changes, seq_details = compare_sequences(products_a, products_b)
        seq_str = f"{len(seq_changes)} changes" if seq_changes else "stable"

        # Flag regressions
        flag = ""
        if time_a > 0 and ((time_b - time_a) / time_a) > 0.10:
            flag = " !"
            regressions.append(f"{sample} ({k_reads}): wall time +{delta_pct:.1f}%")
        elif time_a > 0 and ((time_b - time_a) / time_a) < -0.10:
            improvements.append(f"{sample} ({k_reads}): wall time {delta_pct:.1f}%")

        if n_genes_a > n_genes_b:
            regressions.append(
                f"{sample} ({k_reads}): genes lost ({n_genes_a} -> {n_genes_b})"
            )
            flag += " !"

        if seq_changes:
            flag += " !"
            for change in seq_changes:
                regressions.append(f"{sample} ({k_reads}): {change.strip()}")
            all_details.extend(seq_details)

        print(
            f"{sample:<30} {k_reads:>8} {time_a:>7.1f}s {time_b:>7.1f}s "
            f"{delta_str:>8}  {product_str:>10}  {seq_str}{flag}"
        )

    print()

    if all_details and verbose:
        print("SEQUENCE DETAILS:")
        for d in all_details:
            print(d)
        print()

    if regressions:
        print("REGRESSIONS:")
        for r in regressions:
            print(f"  ! {r}")
        print()

    if improvements:
        print("IMPROVEMENTS:")
        for i in improvements:
            print(f"  + {i}")
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
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Show detailed per-product sequence diffs"
    )
    args = parser.parse_args()

    n_regressions = compare(args.baseline, args.current, verbose=args.verbose)
    sys.exit(1 if n_regressions > 0 else 0)


if __name__ == "__main__":
    main()

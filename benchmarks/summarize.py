#!/usr/bin/env python3
"""
Generate a human-readable summary of benchmark results.

Produces a Markdown file with one table per panel. Each row is a
(sample, max_reads) combination, each column is a gene. Cells show:

    ✓  product recovered, BLAST hit confirms identity
    ?  product recovered, no significant BLAST hit
    ·  no product recovered

Usage:
    python benchmarks/summarize.py benchmarks/results/RESULT.yaml
    python benchmarks/summarize.py benchmarks/results/RESULT.yaml -o summary.md
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import yaml


# Status icons (plain Unicode, renders everywhere)
HIT = "\u2713"    # ✓  product + BLAST hit
NO_BLAST = "?"    # ?  product, no BLAST hit
EMPTY = "\u00b7"  # ·  no product


def load_results(path):
    with open(path) as f:
        return yaml.safe_load(f)


def strip_panel_prefix(gene, panel):
    """Remove panel prefix from gene name for compact column headers.

    e.g. 'cnidaria_18S' with panel 'cnidaria' -> '18S'
         'bacteria_16S-515F-806R' with panel 'bacteria' -> '16S-515F-806R'
    """
    for p in panel.split(", "):
        prefix = p + "_"
        if gene.startswith(prefix):
            return gene[len(prefix):]
    return gene


def classify_product(product):
    """Classify a product as HIT, NO_BLAST, or EMPTY."""
    if not product:
        return EMPTY

    blast = product.get("blast_hit")
    if blast is None:
        # No BLAST annotation at all (BLAST not run)
        return NO_BLAST

    if blast.get("accession") is not None:
        return HIT
    elif blast.get("error"):
        return NO_BLAST
    else:
        return NO_BLAST


def summarize(result_path, output_path=None):
    """Generate a summary from a benchmark result file."""
    benchmark = load_results(result_path)
    results = benchmark.get("results", [])

    version = benchmark.get("sharkmer_version", "?")
    commit = benchmark.get("git_commit", "?")
    date = benchmark.get("date", "?")

    # Group results by panel
    # panel -> [(sample, max_reads, {gene: product})]
    panel_data = defaultdict(list)

    for r in results:
        if r.get("status") == "failed":
            continue

        sample = r["sample"]
        max_reads = r.get("max_reads", 0)
        panel = r.get("panel", "unknown")
        products = {p["gene"]: p for p in r.get("products", [])}

        panel_data[panel].append((sample, max_reads, products))

    # Build output
    lines = []
    lines.append(f"# Benchmark Summary")
    lines.append(f"")
    lines.append(f"- **Version**: {version}")
    lines.append(f"- **Commit**: {commit}")
    lines.append(f"- **Date**: {date}")
    lines.append(f"- **Source**: `{Path(result_path).name}`")
    lines.append(f"")
    lines.append(f"Legend: `{HIT}` product + BLAST hit &ensp; "
                 f"`{NO_BLAST}` product, no BLAST hit &ensp; "
                 f"`{EMPTY}` no product")
    lines.append(f"")

    for panel in sorted(panel_data.keys()):
        rows = panel_data[panel]

        # Collect all genes seen in this panel, preserve order from first appearance
        seen_genes = []
        for _, _, products in rows:
            for gene in products:
                if gene not in seen_genes:
                    seen_genes.append(gene)

        # Short gene names (strip panel prefix)
        short_names = [strip_panel_prefix(g, panel) for g in seen_genes]

        lines.append(f"## {panel}")
        lines.append(f"")

        # Header
        header = "| Sample | Reads |"
        for sn in short_names:
            header += f" {sn} |"
        lines.append(header)

        # Separator
        sep = "| --- | --- |"
        for _ in short_names:
            sep += " --- |"
        lines.append(sep)

        # Data rows
        for sample, max_reads, products in rows:
            k_reads = f"{max_reads // 1000}k"
            row = f"| {sample} | {k_reads} |"
            for gene in seen_genes:
                product = products.get(gene)
                icon = classify_product(product)
                row += f" {icon} |"
            lines.append(row)

        lines.append(f"")

    text = "\n".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(text)
        print(f"Summary written to: {output_path}")
    else:
        print(text)

    return output_path or "(stdout)"


def main():
    parser = argparse.ArgumentParser(
        description="Generate human-readable benchmark summary"
    )
    parser.add_argument(
        "result_file",
        help="Benchmark result YAML file"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file (default: stdout)"
    )
    args = parser.parse_args()

    result_path = Path(args.result_file)
    if not result_path.exists():
        print(f"File not found: {result_path}")
        sys.exit(1)

    summarize(result_path, output_path=args.output)


if __name__ == "__main__":
    main()

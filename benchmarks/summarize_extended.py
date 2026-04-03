#!/usr/bin/env python3
"""
Generate an extended benchmark summary with primer binding analysis.

Cross-references benchmark results with primer_binding.yaml to annotate
each gene × sample cell with primer match information. Highlights failures
where primers match perfectly — these are the actionable graph traversal
issues rather than primer design problems.

Usage:
    python benchmarks/summarize_extended.py benchmarks/results/RESULT.yaml
    python benchmarks/summarize_extended.py benchmarks/results/RESULT.yaml -o extended.md
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import yaml


# Status icons
HIT = "\u2713"    # ✓  product + BLAST hit
NO_BLAST = "?"    # ?  product, no BLAST hit
EMPTY = "\u00b7"  # ·  no product

PANELS_DIR = Path(__file__).parent.parent / "panels"
PRIMER_BINDING_PATH = Path(__file__).parent.parent / "dev_docs" / "primer_binding.yaml"
CONFIG_PATH = Path(__file__).parent / "config.yaml"

# Gene type classification for categorized tables
GENE_TYPES = {
    # Cnidaria panel
    "cnidaria_16S": "mt", "cnidaria_CO1": "mt",
    "cnidaria_18S": "rRNA", "cnidaria_28S": "rRNA",
    "cnidaria_ITS": "rRNA", "cnidaria_ITS-v2": "rRNA",
    "cnidaria_28S-v2": "rRNA",
    "cnidaria_EF1A": "nuclear",
    # Insecta panel
    "insecta_12S": "mt", "insecta_16S": "mt",
    "insecta_CO1": "mt", "insecta_CO2": "mt",
    "insecta_CytB": "mt", "insecta_ND1": "mt",
    "insecta_ND4": "mt", "insecta_ND5": "mt",
    "insecta_16S-v2": "mt", "insecta_CO1-v2": "mt",
    "insecta_CO2-v2": "mt", "insecta_NADH": "mt",
    "insecta_18S": "rRNA", "insecta_28S": "rRNA",
    "insecta_ITS": "rRNA", "insecta_18S-v2": "rRNA",
    "insecta_28S-v2": "rRNA", "insecta_ITS-v2": "rRNA",
    "insecta_ITS-v3": "rRNA",
    "insecta_EF1g": "nuclear", "insecta_Fz4": "nuclear",
    "insecta_Gpdh": "nuclear", "insecta_Pgi": "nuclear",
    "insecta_Yp2": "nuclear",
    # Teleostei panel
    "teleostei_18S": "rRNA", "teleostei_16S": "mt",
    "teleostei_CO1": "mt", "teleostei_CytB": "mt",
    "teleostei_12S": "mt",
    # Human panel — all mitochondrial
    "human_mt1404-3947": "mt", "human_mt3734-6739": "mt",
    "human_mt6511-9220": "mt", "human_mt8910-10648": "mt",
    "human_mt10360-12226": "mt", "human_mt11977-13830": "mt",
    "human_mt13477-15349": "mt", "human_mt14898-151": "mt",
    "human_mt16488-1677": "mt",
    # Bacteria panel — all 16S rRNA
    "bacteria_16S-27F-338R": "rRNA", "bacteria_16S-V2f-V3r": "rRNA",
    "bacteria_16S-341F-785R": "rRNA",
    "bacteria_16S-PRK341F-PRK806R": "rRNA",
    "bacteria_16S-515F-806R": "rRNA",
    "bacteria_16S-515F-806RB": "rRNA",
    "bacteria_16S-515F-Y-926R": "rRNA",
    "bacteria_16S-B969F-BA1406R": "rRNA",
    "bacteria_16S-799F-1391R": "rRNA",
    "bacteria_16S-967F-1391R": "rRNA",
    "bacteria_16S-799F-1193R": "rRNA",
    "bacteria_16S-68F-783Rabc": "rRNA",
    "bacteria_16S-68F-518R": "rRNA",
    "bacteria_16S-341F-783Rabc": "rRNA",
    # Angiospermae panel — all chloroplast/rRNA
    "angiospermae_psbA-trnH": "plastid",
    "angiospermae_rpl36-infA-rps8": "plastid",
    "angiospermae_trnK-rps16": "plastid",
    "angiospermae_trnV-atpE": "plastid",
    "angiospermae_trnC-ycf6": "plastid",
    "angiospermae_ycf6-psbM": "plastid",
    "angiospermae_psbM-trnD": "plastid",
    "angiospermae_atpB-rbcL": "plastid",
    "angiospermae_trnL-F": "plastid",
    "angiospermae_ITS": "rRNA",
}

GENE_TYPE_LABELS = {
    "mt": "mitochondrial",
    "rRNA": "rRNA",
    "nuclear": "nuclear",
    "plastid": "plastid",
}


def load_results(path):
    with open(path) as f:
        return yaml.safe_load(f)


def load_panel_genes(panel_str):
    genes = []
    for panel_name in panel_str.split(", "):
        panel_file = PANELS_DIR / f"{panel_name}.yaml"
        if not panel_file.exists():
            continue
        with open(panel_file) as f:
            panel = yaml.safe_load(f)
        for primer in panel.get("primers", []):
            prefixed = f"{panel_name}_{primer['gene_name']}"
            if prefixed not in genes:
                genes.append(prefixed)
    return genes


def load_primer_binding():
    """Load primer binding analysis. Returns dict keyed by
    '{gene}__{sample}' -> {forward_mismatches, reverse_mismatches, notes}
    """
    if not PRIMER_BINDING_PATH.exists():
        return {}
    with open(PRIMER_BINDING_PATH) as f:
        data = yaml.safe_load(f)
    return data.get("primer_binding", {})


def load_config():
    """Load benchmark config for genome size info."""
    if not CONFIG_PATH.exists():
        return {}
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def strip_panel_prefix(gene, panel):
    for p in panel.split(", "):
        prefix = p + "_"
        if gene.startswith(prefix):
            return gene[len(prefix):]
    return gene


def classify_product(product):
    if not product:
        return EMPTY
    blast = product.get("blast_hit")
    if blast is None:
        return NO_BLAST
    if blast.get("accession") is not None:
        return HIT
    return NO_BLAST


def lookup_primer_binding(primer_binding, gene, sample):
    """Look up primer binding info for a gene/sample pair.

    The primer_binding.yaml keys are like 'cnidaria_16S__Porites_lutea'.
    Returns (fwd_mm, rev_mm, notes) or None if no data.
    """
    key = f"{gene}__{sample}"
    entry = primer_binding.get(key)
    if entry is None:
        return None
    fwd = entry.get("forward_mismatches")
    rev = entry.get("reverse_mismatches")
    return (fwd, rev, entry.get("notes", ""))


def format_mm(fwd, rev):
    """Format mismatch counts as compact string like '0+0' or '2+1'."""
    f = str(fwd) if fwd is not None else "?"
    r = str(rev) if rev is not None else "?"
    return f"{f}+{r}"


def is_perfect_match(fwd, rev):
    """True if both primer directions have 0 mismatches."""
    return fwd == 0 and rev == 0


def summarize_extended(result_path, output_path=None):
    benchmark = load_results(result_path)
    results = benchmark.get("results", [])
    primer_binding = load_primer_binding()
    config = load_config()

    version = benchmark.get("sharkmer_version", "?")
    commit = benchmark.get("git_commit", "?")
    date = benchmark.get("date", "?")

    # Group results by panel
    panel_data = defaultdict(list)
    for r in results:
        if r.get("status") == "failed":
            continue
        sample = r["sample"]
        max_reads = r.get("max_reads", 0)
        panel = r.get("panel", "unknown")
        products = {p["gene"]: p for p in r.get("products", [])}
        panel_data[panel].append((sample, max_reads, products))

    # Collect genome sizes from config
    genome_sizes = {}
    if config and "sample" in config:
        for sname, sdata in config["sample"].items():
            if isinstance(sdata, dict) and "genome" in sdata:
                genome_sizes[sname] = sdata["genome"].get("size_mb")

    lines = []
    lines.append("# Extended Benchmark Summary")
    lines.append("")
    lines.append(f"- **Version**: {version}")
    lines.append(f"- **Commit**: {commit}")
    lines.append(f"- **Date**: {date}")
    lines.append(f"- **Source**: `{Path(result_path).name}`")
    lines.append("")
    lines.append("## Legend")
    lines.append("")
    lines.append("Recovery status:")
    lines.append(f"- `{HIT}` product recovered, BLAST confirms identity")
    lines.append(f"- `{NO_BLAST}` product recovered, no BLAST validation")
    lines.append(f"- `{EMPTY}` no product recovered")
    lines.append("")
    lines.append("Primer match (fwd+rev mismatches):")
    lines.append("- **0+0** perfect match in both directions")
    lines.append("- **2+1** 2 forward mismatches, 1 reverse mismatch")
    lines.append("- **?+?** no public reference sequence available")
    lines.append("- *blank* gene recovered (primer binding not assessed)")
    lines.append("")
    lines.append("Gene type: mt = mitochondrial, rRNA = ribosomal RNA, "
                 "nuc = nuclear, cp = plastid")
    lines.append("")

    # Collect all failures with perfect primer matches for the analysis section
    perfect_match_failures = []

    for panel in sorted(panel_data.keys()):
        rows = panel_data[panel]

        seen_genes = load_panel_genes(panel)
        if not seen_genes:
            for _, _, products in rows:
                for gene in products:
                    if gene not in seen_genes:
                        seen_genes.append(gene)

        short_names = [strip_panel_prefix(g, panel) for g in seen_genes]
        gene_types = [GENE_TYPES.get(g, "?") for g in seen_genes]

        lines.append(f"## {panel}")
        lines.append("")

        # Gene type row + header
        header = "| Sample | Reads |"
        type_row = "| | |"
        for sn, gt in zip(short_names, gene_types):
            header += f" {sn} |"
            label = {"mt": "mt", "rRNA": "rRNA", "nuclear": "nuc",
                     "plastid": "cp", "?": "?"}.get(gt, "?")
            type_row += f" *{label}* |"

        sep = "| --- | --- |"
        for _ in short_names:
            sep += " :---: |"

        lines.append(header)
        lines.append(sep)
        lines.append(type_row)

        for sample, max_reads, products in rows:
            k_reads = f"{max_reads // 1000}k"
            row = f"| {sample} | {k_reads} |"

            for gene, gt in zip(seen_genes, gene_types):
                product = products.get(gene)
                icon = classify_product(product)
                binding = lookup_primer_binding(primer_binding, gene, sample)

                if icon == EMPTY and binding is not None:
                    fwd_mm, rev_mm, notes = binding
                    mm_str = format_mm(fwd_mm, rev_mm)
                    row += f" {icon} {mm_str} |"

                    if fwd_mm is not None and rev_mm is not None:
                        if isinstance(fwd_mm, int) and isinstance(rev_mm, int):
                            if fwd_mm + rev_mm <= 1:
                                perfect_match_failures.append({
                                    "panel": panel,
                                    "gene": gene,
                                    "short_name": strip_panel_prefix(gene, panel),
                                    "sample": sample,
                                    "reads": k_reads,
                                    "fwd_mm": fwd_mm,
                                    "rev_mm": rev_mm,
                                    "gene_type": gt,
                                    "notes": notes,
                                    "genome_mb": genome_sizes.get(sample),
                                })
                elif icon == EMPTY:
                    row += f" {icon} |"
                else:
                    row += f" {icon} |"

            lines.append(row)

        lines.append("")

    # ── Analysis section: perfect-match failures ──
    lines.append("---")
    lines.append("")
    lines.append("## Analysis: failures with perfect or near-perfect primer matches")
    lines.append("")
    lines.append("These are genes where primers bind with 0-1 total mismatches "
                 "but no product is recovered at 1M reads. These represent graph "
                 "traversal failures, not primer design problems.")
    lines.append("")

    if perfect_match_failures:
        lines.append("| Panel | Gene | Type | Sample | Genome | Fwd+Rev mm | Notes |")
        lines.append("| --- | --- | --- | --- | --- | :---: | --- |")

        for f in perfect_match_failures:
            genome = f"{f['genome_mb']:.0f} Mb" if f["genome_mb"] else "?"
            mm = format_mm(f["fwd_mm"], f["rev_mm"])
            # Truncate notes for table
            notes = f["notes"].strip().replace("\n", " ")
            if len(notes) > 120:
                notes = notes[:117] + "..."
            gene_label = {"mt": "mt", "rRNA": "rRNA", "nuclear": "nuc",
                          "plastid": "cp"}.get(f["gene_type"], "?")
            lines.append(
                f"| {f['panel']} | {f['short_name']} | {gene_label} | "
                f"{f['sample']} | {genome} | {mm} | {notes} |"
            )
        lines.append("")
    else:
        lines.append("*No perfect-match failures found.*")
        lines.append("")

    # Summary statistics
    lines.append("## Recovery statistics by gene type")
    lines.append("")

    # Collect stats across all panels
    type_stats = defaultdict(lambda: {"recovered": 0, "total": 0})
    for panel, rows in panel_data.items():
        seen_genes = load_panel_genes(panel)
        if not seen_genes:
            for _, _, products in rows:
                for gene in products:
                    if gene not in seen_genes:
                        seen_genes.append(gene)

        for sample, max_reads, products in rows:
            for gene in seen_genes:
                gt = GENE_TYPES.get(gene, "unknown")
                product = products.get(gene)
                icon = classify_product(product)
                type_stats[gt]["total"] += 1
                if icon != EMPTY:
                    type_stats[gt]["recovered"] += 1

    lines.append("| Gene type | Recovered | Total | Rate |")
    lines.append("| --- | ---: | ---: | ---: |")
    for gt in ["mt", "rRNA", "nuclear", "plastid", "unknown"]:
        if gt not in type_stats:
            continue
        s = type_stats[gt]
        rate = f"{100 * s['recovered'] / s['total']:.0f}%" if s["total"] else "—"
        label = GENE_TYPE_LABELS.get(gt, gt)
        lines.append(f"| {label} | {s['recovered']} | {s['total']} | {rate} |")

    total_rec = sum(s["recovered"] for s in type_stats.values())
    total_all = sum(s["total"] for s in type_stats.values())
    total_rate = f"{100 * total_rec / total_all:.0f}%" if total_all else "—"
    lines.append(f"| **total** | **{total_rec}** | **{total_all}** | **{total_rate}** |")
    lines.append("")

    text = "\n".join(lines)

    if output_path:
        with open(output_path, "w") as f:
            f.write(text)
        print(f"Extended summary written to: {output_path}")
    else:
        print(text)


def main():
    parser = argparse.ArgumentParser(
        description="Generate extended benchmark summary with primer binding analysis"
    )
    parser.add_argument("result_file", help="Benchmark result YAML file")
    parser.add_argument("-o", "--output", help="Output file (default: stdout)")
    args = parser.parse_args()

    result_path = Path(args.result_file)
    if not result_path.exists():
        print(f"File not found: {result_path}")
        sys.exit(1)

    summarize_extended(result_path, output_path=args.output)


if __name__ == "__main__":
    main()

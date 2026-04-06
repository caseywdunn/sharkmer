"""Depth-focused markdown report generation.

Both validation and benchmarks produce reports with the same sections. This
module generates them from result dicts (as built by results.build_result()).
"""

from datetime import datetime
from pathlib import Path

from . import primer_analysis, runner


# ---------------------------------------------------------------------------
# Per-panel report (used by both validation and benchmark)
# ---------------------------------------------------------------------------


def write_panel_report(
    result: dict,
    panel_data: dict,
    sample_results: list,
    report_path: Path,
    gene_filter: list | None = None,
):
    """Write a full markdown report for one panel.

    result: the result dict from results.build_result()
    panel_data: the panel YAML data
    sample_results: list of (sample_block, runs) for primer analysis
    """
    lines = []
    panel_name = result.get("panel", "unknown")
    panel_version = result.get("panel_version", "unversioned")
    sharkmer_version = result.get("sharkmer_version", "?")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Section 1: Header
    lines.append(f"# Panel validation: {panel_name}")
    lines.append("")
    lines.append(f"- **Panel version**: `{panel_version}`")
    lines.append(f"- **sharkmer version**: `{sharkmer_version}`")
    lines.append(f"- **Date**: {now}")
    machine = result.get("machine", {})
    if machine:
        lines.append(
            f"- **Machine**: {machine.get('os', '?')}, "
            f"{machine.get('cpu_cores', '?')} cores, "
            f"{machine.get('total_ram_gb', '?')} GB RAM"
        )
    if gene_filter:
        lines.append(f"- **Gene filter**: {', '.join(gene_filter)}")
    lines.append("")

    declared_genes = runner.panel_gene_names(panel_data)
    if gene_filter:
        considered_genes = [g for g in declared_genes if g in gene_filter]
    else:
        considered_genes = declared_genes

    # Section 2: Depth-recovery matrix (one per sample)
    for sample_entry in result.get("samples", []):
        accession = sample_entry["accession"]
        taxon = sample_entry.get("taxon", "")
        depths = sample_entry.get("depths", [])

        heading = f"## {taxon} ({accession})" if taxon else f"## {accession}"
        lines.append(heading)
        lines.append("")

        # Check for failed runs.
        failed = [d for d in depths if not d.get("success", True)]
        if failed:
            failed_reads = ", ".join(
                f"{d['max_reads'] // 1000}k" for d in failed
            )
            lines.append(f"Failed runs: {failed_reads}")
            lines.append("")

        successful = [d for d in depths if d.get("success", True)]
        if not successful:
            lines.append("_No successful runs._")
            lines.append("")
            continue

        # Sort depths ascending for left-to-right reading.
        successful.sort(key=lambda d: d["max_reads"])

        # Build gene -> {max_reads: gene_result} map.
        depth_data: dict[str, dict[int, dict]] = {}
        for depth in successful:
            for gene_result in depth.get("genes", []):
                gene = gene_result["gene"]
                depth_data.setdefault(gene, {})[depth["max_reads"]] = gene_result

        # Find best reference match per gene (from highest depth).
        best_ref: dict[str, dict] = {}
        for gene in considered_genes:
            gene_depths = depth_data.get(gene, {})
            for depth in reversed(successful):
                gr = gene_depths.get(depth["max_reads"])
                if gr and gr.get("recovered") and gr.get("reference_match"):
                    best_ref[gene] = gr["reference_match"]
                    break

        # Render table.
        depth_headers = [f"{d['max_reads'] // 1000}k" for d in successful]
        header = "| Gene | " + " | ".join(depth_headers) + " | Best ref match |"
        sep = "|------|" + "|".join(["---:"] * len(successful)) + "|----------------|"
        lines.append(header)
        lines.append(sep)

        for gene in considered_genes:
            gene_depths = depth_data.get(gene, {})
            cells = []
            for depth in successful:
                gr = gene_depths.get(depth["max_reads"])
                if gr and gr.get("recovered"):
                    length = gr.get("length")
                    ref = gr.get("reference_match")
                    if ref and ref.get("pct_identity") is not None:
                        cells.append(f"{length}bp ({ref['pct_identity']}%)")
                    elif length:
                        cells.append(f"{length}bp")
                    else:
                        cells.append("?")
                else:
                    cells.append("---")

            # Best ref match column.
            ref = best_ref.get(gene)
            if ref and ref.get("matched_taxon"):
                pct = ref.get("pct_identity")
                pct_str = f" {pct}%" if pct is not None else ""
                target = " ON" if ref.get("on_target") else " OFF-TARGET"
                ref_str = f"{ref['matched_taxon']}{pct_str}{target}"
            elif ref and ref.get("error"):
                ref_str = f"error: {ref['error'][:30]}"
            else:
                ref_str = "(no reference)"

            lines.append(
                f"| {gene} | " + " | ".join(cells) + f" | {ref_str} |"
            )

        lines.append("")

        # Wall time summary for this sample.
        times = [
            f"{d['max_reads'] // 1000}k: {d.get('wall_time_s', '?')}s"
            for d in successful
        ]
        lines.append(f"Wall times: {', '.join(times)}")
        lines.append("")

    # Section 3: Cross-sample summary at highest depth.
    lines.extend(_cross_sample_summary(result, considered_genes))

    # Section 4: Primer binding analysis.
    if sample_results:
        analyses = primer_analysis.analyze_primer_bindings(
            panel_data, sample_results, considered_genes
        )
        if analyses:
            lines.extend(_format_binding_section(analyses))

    # Section 5: Reference match details.
    lines.extend(_reference_details(result, considered_genes))

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Report written to: {report_path}")


# ---------------------------------------------------------------------------
# Cross-sample summary
# ---------------------------------------------------------------------------


def _cross_sample_summary(result: dict, considered_genes: list) -> list:
    """Genes x samples pass/fail grid at highest depth."""
    samples = result.get("samples", [])
    if len(samples) < 2:
        return []

    lines = []
    lines.append("## Cross-sample summary (highest depth)")
    lines.append("")

    # Column headers: sample labels.
    labels = []
    for s in samples:
        taxon = s.get("taxon", "")
        label = taxon if taxon else s["accession"]
        # Truncate long taxon names.
        if len(label) > 20:
            label = label[:17] + "..."
        labels.append(label)

    header = "| Gene | " + " | ".join(labels) + " |"
    sep = "|------|" + "|".join([":---:"] * len(samples)) + "|"
    lines.append(header)
    lines.append(sep)

    for gene in considered_genes:
        cells = []
        for s in samples:
            depths = s.get("depths", [])
            successful = [d for d in depths if d.get("success", True)]
            if not successful:
                cells.append("---")
                continue
            # Take highest depth.
            best = max(successful, key=lambda d: d["max_reads"])
            gene_results = {g["gene"]: g for g in best.get("genes", [])}
            gr = gene_results.get(gene)
            if gr and gr.get("recovered"):
                ref = gr.get("reference_match")
                if ref and ref.get("on_target"):
                    cells.append("PASS")
                elif ref and ref.get("pct_identity") is not None:
                    if not ref.get("on_target"):
                        cells.append("OFF-TARGET")
                    else:
                        cells.append("PASS")
                else:
                    cells.append("recovered")
            else:
                cells.append("---")
        lines.append(f"| {gene} | " + " | ".join(cells) + " |")

    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Primer binding section formatting
# ---------------------------------------------------------------------------


def _format_binding_section(analyses: list) -> list:
    """Render primer binding analysis for the markdown report."""
    lines = []
    lines.append("## Primer binding analysis")
    lines.append("")
    lines.append(
        "For each gene, the first and last `trim` bases of each recovered "
        "amplicon are compared against the user-specified primer sequence "
        "(3'-trimmed to the match window). Reverse primer bindings are shown "
        "reverse-complemented so they appear in the same orientation as the "
        "primer was written in the panel. "
        "Use this section to decide whether a primer's degeneracy should be "
        "reduced (only a subset of coded bases is actually seen) or widened "
        "(an off-code base was absorbed by sharkmer's `--mismatches` "
        "tolerance)."
    )
    lines.append("")

    for a in analyses:
        gene = a["gene_name"]
        trim = a["trim"]
        lines.append(f"### {gene}")
        lines.append("")
        if a["missing_samples"]:
            lines.append(
                f"**Not recovered** in: {', '.join(a['missing_samples'])}. "
                "This may indicate insufficient primer degeneracy, "
                "insufficient read coverage, or the target taxon lacks the "
                "locus."
            )
            lines.append("")
        if not a["forward"]["per_sample"] and not a["reverse"]["per_sample"]:
            lines.append(
                "_(no samples recovered this gene; skipping alignment)_"
            )
            lines.append("")
            continue

        for which, title in (
            ("forward", "Forward primer"),
            ("reverse", "Reverse primer"),
        ):
            info = a[which]
            full = a[f"{which}_full"]
            lines.append(
                f"**{title}** — spec `{full}` (trim={trim}, match window "
                f"`{info['spec']}`)"
            )
            lines.append("")
            lines.append("```")
            lines.append(f"  {'spec':<14}{info['spec']}")
            marker = []
            for pos in info["position_analysis"]:
                if pos["status"] == "fixed":
                    marker.append("|")
                elif pos["status"] in ("fully_utilised", "could_reduce"):
                    marker.append(".")
                else:
                    marker.append("x")
            lines.append(f"  {'':<14}{''.join(marker)}")
            for s in info["per_sample"]:
                lines.append(f"  {s['accession']:<14}{s['observed']}")
            lines.append("```")
            lines.append("")
            for v in info["verdict_lines"]:
                lines.append(v)
            lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Reference match details
# ---------------------------------------------------------------------------


def _reference_details(result: dict, considered_genes: list) -> list:
    """Detailed reference match table."""
    rows = []
    for s in result.get("samples", []):
        accession = s["accession"]
        depths = s.get("depths", [])
        successful = [d for d in depths if d.get("success", True)]
        if not successful:
            continue
        best = max(successful, key=lambda d: d["max_reads"])
        for gr in best.get("genes", []):
            gene = gr.get("gene")
            if gene not in considered_genes:
                continue
            if not gr.get("recovered"):
                continue
            ref = gr.get("reference_match")
            if ref is None:
                continue
            rows.append(
                {
                    "sample": accession,
                    "gene": gene,
                    "matched_taxon": ref.get("matched_taxon", "---"),
                    "matched_accession": ref.get("matched_accession", "---"),
                    "pct_identity": ref.get("pct_identity"),
                    "align_length": ref.get("align_length"),
                    "on_target": ref.get("on_target", False),
                }
            )

    if not rows:
        return []

    lines = []
    lines.append("## Reference match details")
    lines.append("")
    lines.append(
        "| Sample | Gene | Ref taxon | Ref accession | Identity | "
        "Align len | On-target |"
    )
    lines.append(
        "|--------|------|-----------|---------------|----------|"
        "-----------|-----------|"
    )
    for r in rows:
        pct = f"{r['pct_identity']:.1f}%" if r["pct_identity"] is not None else "---"
        alen = str(r["align_length"]) if r["align_length"] is not None else "---"
        target = "yes" if r["on_target"] else "NO"
        lines.append(
            f"| {r['sample']} | {r['gene']} | {r['matched_taxon']} | "
            f"{r['matched_accession']} | {pct} | {alen} | {target} |"
        )
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Benchmark combined summary (cross-panel)
# ---------------------------------------------------------------------------


def write_benchmark_summary(
    panel_results: list,
    summary_path: Path,
):
    """Write a combined benchmark summary spanning multiple panels.

    panel_results is a list of result dicts (one per panel).
    """
    lines = []
    lines.append("# Benchmark Summary")
    lines.append("")

    if panel_results:
        r0 = panel_results[0]
        lines.append(f"- **sharkmer version**: `{r0.get('sharkmer_version', '?')}`")
        lines.append(f"- **Commit**: `{r0.get('git_commit', '?')}`")
        lines.append(f"- **Date**: {r0.get('date', '?')}")
        machine = r0.get("machine", {})
        if machine:
            lines.append(
                f"- **Machine**: {machine.get('os', '?')}, "
                f"{machine.get('cpu_cores', '?')} cores, "
                f"{machine.get('total_ram_gb', '?')} GB RAM"
            )
    lines.append("")

    # One section per panel with a compact recovery matrix.
    for result in panel_results:
        panel_name = result.get("panel", "unknown")
        panel_version = result.get("panel_version", "?")
        lines.append(f"## {panel_name} v{panel_version}")
        lines.append("")

        samples = result.get("samples", [])
        if not samples:
            lines.append("_No samples._")
            lines.append("")
            continue

        # Collect all genes across samples.
        all_genes = []
        for s in samples:
            for d in s.get("depths", []):
                for g in d.get("genes", []):
                    if g["gene"] not in all_genes:
                        all_genes.append(g["gene"])

        # Table: rows = (sample, depth), columns = genes.
        header = "| Sample | Reads | " + " | ".join(all_genes) + " |"
        sep = "| --- | ---: | " + " | ".join(["---"] * len(all_genes)) + " |"
        lines.append(header)
        lines.append(sep)

        for s in samples:
            label = s.get("taxon") or s["accession"]
            if len(label) > 20:
                label = label[:17] + "..."
            for d in sorted(s.get("depths", []), key=lambda x: x["max_reads"]):
                if not d.get("success", True):
                    continue
                k_reads = f"{d['max_reads'] // 1000}k"
                gene_map = {g["gene"]: g for g in d.get("genes", [])}
                cells = []
                for gene in all_genes:
                    gr = gene_map.get(gene)
                    if gr and gr.get("recovered"):
                        length = gr.get("length", "?")
                        ref = gr.get("reference_match")
                        if ref and ref.get("on_target"):
                            cells.append(f"{length}bp")
                        elif ref and ref.get("pct_identity") is not None:
                            cells.append(f"{length}bp*")
                        else:
                            cells.append(f"{length}bp?")
                    else:
                        cells.append("---")
                lines.append(
                    f"| {label} | {k_reads} | " + " | ".join(cells) + " |"
                )
        lines.append("")
        lines.append("Legend: `Nbp` = on-target, `Nbp*` = off-target, `Nbp?` = no reference, `---` = not recovered")
        lines.append("")

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Benchmark summary written to: {summary_path}")

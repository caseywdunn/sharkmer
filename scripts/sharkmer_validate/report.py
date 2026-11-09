"""Depth-focused markdown report generation.

Both validation and benchmarks produce reports with the same sections. This
module generates them from result dicts (as built by results.build_result()).
"""

from datetime import datetime
from pathlib import Path

from . import blast_references, primer_analysis, runner
from .reference_targets import logical_gene_name, target_logical_genes


# ---------------------------------------------------------------------------
# Three-position scoring code
# ---------------------------------------------------------------------------
#
# Each gene × sample is scored with a 3-character code:
#
#   Position 1 — Recovery:   `-` not recovered, `+` recovered
#   Position 2 — Reference:  `-` no verified ref for this gene,
#                             `+` verified ref for other taxa only,
#                             `*` verified ref for this taxon
#   Position 3 — Alignment:  `-` no sufficient same-gene support,
#                             `+` gene support, other taxon,
#                             `*` gene and expected-taxon support
#
# Possible codes:
#   ---  not recovered, no references for this gene
#   -+-  not recovered, refs exist for other species
#   -*-  not recovered, ref exists for this species
#   +--  recovered, no references for this gene
#   ++-  recovered, refs for other species, no hit
#   +++  recovered, refs for other species, hit same gene different species
#   +*-  recovered, ref for this species, no hit (suspicious)
#   +*+  recovered, ref for this species, hit different species (unexpected)
#   +**  recovered, reference supports same gene and expected taxon


SCORE_LEGEND = (
    "**Scoring** — three positions: recovery / verified-reference availability / alignment support.\n"
    "\n"
    "| Code | Meaning |\n"
    "|------|---------|\n"
    "| `+**` | Recovered; alignment supports the gene and expected taxon |\n"
    "| `+*+` | Recovered; expected-taxon reference exists but strongest gene support is another taxon |\n"
    "| `+*-` | Recovered; no sufficient unambiguous same-gene alignment |\n"
    "| `+++` | Recovered; alignment supports the gene using another taxon |\n"
    "| `++-` | Recovered; no sufficient unambiguous same-gene alignment |\n"
    "| `+--` | Recovered; no verified reference for this gene |\n"
    "| `-*-` | Not recovered; verified reference exists for this taxon |\n"
    "| `-+-` | Not recovered; verified references exist for other taxa |\n"
    "| `---` | Not recovered; no verified reference for this gene |\n"
    "\n"
    "Position 1: `-` no product, `+` product recovered. "
    "Position 2: `-` no verified reference for this gene, "
    "`+` verified reference for other taxa, `*` verified reference for this taxon. "
    "Position 3: `-` insufficient, ambiguous, or conflicting evidence; "
    "`+` same-gene support from another taxon; `*` same-gene and expected-taxon support.\n"
    "Identity and query-coverage gates are recorded analysis criteria, not universal truth. "
    "Reference provenance validates the public source region, not its panel gene annotation. "
    "Reference alignment does not establish sample haplotype truth or read support.\n"
)


def _build_ref_availability(
    panel_data: dict | None = None, reference_summary: dict | None = None
) -> dict:
    """Build a map of reference availability per gene.

    Returns: {gene_name: {taxon1, taxon2, ...}} — set of taxa that have
    a reference for each gene. Genes with no references are absent.
    """
    ref_map = {}
    if isinstance(reference_summary, dict):
        references = reference_summary.get("verified", [])
        target_mapping = reference_summary.get("target_logical_genes", {})
        for reference in references:
            ref_map.setdefault(reference["gene"], set()).add(reference["taxon"])
            reference_logical_gene = reference.get("logical_gene") or logical_gene_name(
                reference["gene"]
            )
            for target, logical_gene in target_mapping.items():
                if logical_gene == reference_logical_gene:
                    ref_map.setdefault(target, set()).add(reference["taxon"])
        return ref_map
    target_mapping = target_logical_genes(panel_data or {})
    for reference in blast_references.extract_references(panel_data or {}):
        ref_map.setdefault(reference["gene_name"], set()).add(reference["taxon"])
        reference_logical_gene = logical_gene_name(
            reference["gene_name"], panel_data=panel_data
        )
        for target, logical_gene in target_mapping.items():
            if logical_gene == reference_logical_gene:
                ref_map.setdefault(target, set()).add(reference["taxon"])
    return ref_map


def _score_gene(
    recovered: bool,
    gene: str,
    sample_taxon: str,
    ref_match: dict | None,
    ref_availability: dict,
) -> str:
    """Compute the 3-position score for a gene × sample."""
    # Position 1: recovery
    if not recovered:
        p1 = "-"
    else:
        p1 = "+"

    # Position 2: reference availability for this gene
    expected_logical_gene = (
        ref_match.get("expected_logical_gene")
        if isinstance(ref_match, dict)
        else logical_gene_name(gene)
    ) or logical_gene_name(gene)
    gene_refs = set().union(
        *(
            taxa
            for target, taxa in ref_availability.items()
            if target == gene
            or logical_gene_name(target) == expected_logical_gene
        )
    )
    if not gene_refs:
        p2 = "-"
    elif sample_taxon in gene_refs:
        p2 = "*"
    else:
        p2 = "+"

    # Position 3: BLAST result (only meaningful if recovered)
    if not recovered or ref_match is None:
        p3 = "-"
    elif (
        ref_match.get("status") == "gene_supported_expected_taxon"
        and (
            ref_match.get("matched_logical_gene")
            or logical_gene_name(ref_match.get("matched_gene"))
        )
        == expected_logical_gene
        and ref_match.get("all_products_expected_taxon_supported", True)
    ):
        p3 = "*"
    elif (
        ref_match.get("status") == "gene_supported_other_taxon"
        and (
            ref_match.get("matched_logical_gene")
            or logical_gene_name(ref_match.get("matched_gene"))
        )
        == expected_logical_gene
    ):
        p3 = "+"
    else:
        p3 = "-"

    return f"{p1}{p2}{p3}"


def _reference_provenance_summary(result: dict) -> list:
    references = result.get("provenance", {}).get("references")
    if not isinstance(references, dict):
        return []
    catalog = references.get("catalog") or {}
    lines = ["## Reference evidence provenance", ""]
    lines.append(
        f"- **Verified reference entries**: {references.get('verified_count', 0)}"
    )
    lines.append(
        f"- **Excluded/unverified entries**: {references.get('excluded_count', 0)}"
    )
    lines.append(
        f"- **Reference catalog**: status `{catalog.get('status', 'unavailable')}`, "
        f"SHA-256 `{catalog.get('sha256') or 'unavailable'}`, path `{catalog.get('path', 'unavailable')}`"
    )
    biological_truth = references.get("biological_truth")
    if biological_truth:
        lines.append(f"- **Evidence scope**: {biological_truth}")
    target_mapping = references.get("target_logical_genes") or {}
    grouped_targets = {}
    for target, logical_gene in target_mapping.items():
        grouped_targets.setdefault(logical_gene, []).append(target)
    reviewed_groups = {
        logical_gene: sorted(targets)
        for logical_gene, targets in grouped_targets.items()
        if len(targets) > 1
    }
    if reviewed_groups:
        rendered_groups = "; ".join(
            f"`{logical_gene}`: {', '.join(f'`{target}`' for target in targets)}"
            for logical_gene, targets in sorted(reviewed_groups.items())
        )
        lines.append(f"- **Reviewed logical target groups**: {rendered_groups}")
        lines.append(
            "- **Primer-region support**: not established by logical target grouping"
        )
    excluded = references.get("excluded") or []
    if excluded:
        reason_counts = {}
        for entry in excluded:
            reason = entry.get("reason", "unspecified")
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
        lines.extend(["", "| Exclusion reason | Entries |", "|------------------|--------:|"])
        for reason, count in sorted(reason_counts.items()):
            lines.append(f"| {reason} | {count} |")
    lines.extend(
        [
            "",
            "Public-region verification establishes source provenance only. Gene labels remain panel annotations; alignment does not establish sample haplotype truth or read support.",
            "",
        ]
    )
    return lines


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
    executable = result.get("provenance", {}).get("executable") or {}
    if executable:
        source_observation = executable.get("workspace_source_observation", {})
        lines.append(
            f"- **Executable**: `{executable.get('path', '?')}` "
            f"(sha256 `{executable.get('sha256', '?')}`)"
        )
        lines.append(
            f"- **Source observation**: `{source_observation.get('revision', '?')}`; "
            f"dirty={source_observation.get('dirty', '?')}; "
            f"tree sha256 `{source_observation.get('tree_sha256', '?')}`"
        )
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

    stored_reference_summary = result.get("provenance", {}).get("references")
    ref_availability = _build_ref_availability(
        panel_data,
        stored_reference_summary if isinstance(stored_reference_summary, dict) else None,
    )

    lines.extend(_reference_provenance_summary(result))

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

        # Render table.
        depth_headers = [f"{d['max_reads'] // 1000}k" for d in successful]
        header = "| Gene | " + " | ".join(depth_headers) + " | Score |"
        sep = "|------|" + "|".join(["---:"] * len(successful)) + "|:-----:|"
        lines.append(header)
        lines.append(sep)

        for gene in considered_genes:
            gene_depths = depth_data.get(gene, {})
            cells = []
            for depth in successful:
                gr = gene_depths.get(depth["max_reads"])
                if gr and gr.get("evaluation_status") == "not_evaluated":
                    cells.append("N/E")
                elif gr and gr.get("recovered"):
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

            # Score from highest successful depth.
            best_depth = successful[-1]
            gr = gene_depths.get(best_depth["max_reads"])
            if gr and gr.get("evaluation_status") == "not_evaluated":
                score = "N/E"
            elif gr and gr.get("recovered"):
                score = _score_gene(
                    True, gene, taxon,
                    gr.get("reference_match"), ref_availability,
                )
            else:
                score = _score_gene(
                    False, gene, taxon, None, ref_availability,
                )

            lines.append(
                f"| {gene} | " + " | ".join(cells) + f" | `{score}` |"
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
    lines.extend(_cross_sample_summary(result, considered_genes, ref_availability))

    # Section 4: Primer binding analysis.
    if sample_results:
        analyses = primer_analysis.analyze_primer_bindings(
            panel_data, sample_results, considered_genes
        )
        if analyses:
            lines.extend(_format_binding_section(analyses))

    # Section 5: Reference match details.
    lines.extend(_reference_details(result, considered_genes))

    # Section 6: Performance summary.
    lines.extend(_performance_summary(result))

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Report written to: {report_path}")


# ---------------------------------------------------------------------------
# Cross-sample summary
# ---------------------------------------------------------------------------


def _cross_sample_summary(
    result: dict, considered_genes: list, ref_availability: dict,
) -> list:
    """Genes x samples scoring grid at highest depth."""
    samples = result.get("samples", [])
    if not samples:
        return []

    lines = []
    lines.append("## Cross-sample summary (highest depth)")
    lines.append("")

    # Column headers: sample labels.
    labels = []
    for s in samples:
        taxon = s.get("taxon", "")
        label = taxon if taxon else s["accession"]
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
            taxon = s.get("taxon", "")
            depths = s.get("depths", [])
            successful = [d for d in depths if d.get("success", True)]
            if not successful:
                score = _score_gene(False, gene, taxon, None, ref_availability)
                cells.append(f"`{score}`")
                continue
            best = max(successful, key=lambda d: d["max_reads"])
            gene_results = {g["gene"]: g for g in best.get("genes", [])}
            gr = gene_results.get(gene)
            if gr and gr.get("evaluation_status") == "not_evaluated":
                cells.append("`N/E`")
                continue
            if gr and gr.get("recovered"):
                score = _score_gene(
                    True, gene, taxon,
                    gr.get("reference_match"), ref_availability,
                )
            else:
                score = _score_gene(
                    False, gene, taxon, None, ref_availability,
                )
            cells.append(f"`{score}`")
        lines.append(f"| {gene} | " + " | ".join(cells) + " |")

    lines.append("")
    lines.append(SCORE_LEGEND)
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
        taxon = s.get("taxon", "")
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
            for product in gr.get("products", []):
                ref = product.get("reference_match")
                if ref is None:
                    continue
                rows.append(
                    {
                        "sample": accession,
                        "sample_taxon": taxon,
                        "gene": gene,
                        "logical_gene": ref.get("expected_logical_gene")
                        or gr.get("logical_gene")
                        or logical_gene_name(gene),
                        "product": product.get("product_index"),
                        "status": ref.get("status", "unknown"),
                        "target_support": ref.get("target_support", "unavailable"),
                        "sequence_relationship": ref.get("sequence_relationship", "unavailable"),
                        "matched_gene": ref.get("matched_gene", "---"),
                        "matched_logical_gene": ref.get(
                            "matched_logical_gene", "---"
                        ),
                        "matched_taxon": ref.get("matched_taxon", "---"),
                        "matched_accession": ref.get("matched_accession", "---"),
                        "pct_identity": ref.get("pct_identity"),
                        "query_coverage_pct": ref.get("query_coverage_pct"),
                        "query_unaligned_bases": ref.get("query_unaligned_bases"),
                        "unmatched_query_regions": ref.get("unmatched_query_regions"),
                        "reference_coverage_pct": ref.get("reference_coverage_pct"),
                        "unmatched_reference_regions": ref.get("unmatched_reference_regions"),
                        "gap_count": ref.get("gap_count"),
                        "haplotype_truth": ref.get("haplotype_truth", "not_established"),
                        "read_support": ref.get("read_support", "not_evaluated"),
                        "primer_region_support": ref.get(
                            "primer_region_support", "not_established"
                        ),
                    }
                )

    if not rows:
        return []

    lines = []
    lines.append("## Reference match details")
    lines.append("")
    lines.append(
        "| Sample | Target | Logical target | Product | Status | Target support | Sequence relationship | "
        "Matched reference target | Matched logical target | Sample taxon | Ref taxon | Ref accession | Identity | Query coverage | "
        "Unmatched query regions | Reference coverage | Unmatched reference regions | "
        "Gaps | Primer-region support | Haplotype truth | Read support |"
    )
    lines.append(
        "|--------|--------|----------------|--------:|--------|----------------|-----------------------|"
        "--------------------------|------------------------|-------------|-----------|---------------|----------:|---------------:|"
        "-------------------------|-------------------:|-----------------------------|"
        "-----:|----------------------|-----------------|--------------|"
    )
    for r in rows:
        pct = f"{r['pct_identity']:.1f}%" if r["pct_identity"] is not None else "---"
        coverage = (
            f"{r['query_coverage_pct']:.1f}%"
            if r["query_coverage_pct"] is not None
            else "---"
        )
        reference_coverage = (
            f"{r['reference_coverage_pct']:.1f}%"
            if r["reference_coverage_pct"] is not None
            else "---"
        )
        unmatched_query = ", ".join(
            f"{start}-{end}" for start, end in r["unmatched_query_regions"] or []
        ) or "none"
        unmatched_reference = ", ".join(
            f"{start}-{end}" for start, end in r["unmatched_reference_regions"] or []
        ) or "none"
        lines.append(
            f"| {r['sample']} | {r['gene']} | {r['logical_gene']} | {r['product']} | {r['status']} | "
            f"{r['target_support']} | {r['sequence_relationship']} | "
            f"{r['matched_gene']} | {r['matched_logical_gene']} | {r['sample_taxon']} | "
            f"{r['matched_taxon']} | {r['matched_accession']} | {pct} | {coverage} | "
            f"{unmatched_query} | {reference_coverage} | {unmatched_reference} | "
            f"{r['gap_count'] if r['gap_count'] is not None else '---'} | "
            f"{r['primer_region_support']} | {r['haplotype_truth']} | {r['read_support']} |"
        )
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Performance summary
# ---------------------------------------------------------------------------


def _format_bytes(n: int | None) -> str:
    """Format bytes as a human-readable string."""
    if n is None:
        return "---"
    if n < 1024:
        return f"{n} B"
    elif n < 1024 ** 2:
        return f"{n / 1024:.0f} KB"
    elif n < 1024 ** 3:
        return f"{n / 1024 ** 2:.0f} MB"
    else:
        return f"{n / 1024 ** 3:.1f} GB"


def _format_count(n: int | None) -> str:
    """Format a large number with commas."""
    if n is None:
        return "---"
    return f"{n:,}"


def _performance_summary(result: dict) -> list:
    """Performance table: timing, allocator, input, and counting metrics."""
    rows = []
    for s in result.get("samples", []):
        accession = s["accession"]
        taxon = s.get("taxon", "")
        label = taxon if taxon else accession
        if len(label) > 25:
            label = label[:22] + "..."
        for d in s.get("depths", []):
            if not d.get("success", True):
                continue
            stats = d.get("run_stats", {})
            rows.append({
                "sample": label,
                "reads": d["max_reads"],
                "wall_time_s": d.get("wall_time_s"),
                "peak_mem": stats.get("peak_memory_bytes"),
                "n_reads": stats.get("n_reads_read"),
                "n_subreads": stats.get("n_subreads_ingested"),
                "n_bases": stats.get("n_bases_read"),
                "n_kmers": stats.get("n_kmers"),
                "stage_times": stats.get("stage_timings", {}),
                "peak_rss": d.get("peak_rss_bytes"),
                "table_capacity": stats.get("count_table_capacity"),
                "temp_disk": d.get("temp_disk_final_bytes"),
            })

    if not rows:
        return []

    lines = []
    lines.append("## Performance")
    lines.append("")
    lines.append(
        "`Records ingested` is the legacy `n_subreads_ingested` field: it currently "
        "duplicates successfully ingested FASTQ records, while `N` only breaks kmer windows."
    )
    lines.append("")
    lines.append(
        "| Sample | Max reads | Wall time | Count time | PCR time | Allocator peak | "
        "Reads read | Records ingested | Bases read | Kmer occurrences | Mbp/s | k-mers/s | "
        "Peak RSS | Table capacity | Final run disk |"
    )
    lines.append(
        "|--------|----------:|----------:|-----------:|---------:|---------------:|"
        "------------:|------------------:|-----------:|-----------------:|------:|----------:|"
        "---------:|---------------:|---------------:|"
    )
    for r in rows:
        k_reads = f"{r['reads'] // 1000}k"
        wall = f"{r['wall_time_s']}s" if r["wall_time_s"] is not None else "---"
        counting_seconds = sum(
            r["stage_times"].get(key) or 0
            for key in ("read_ingest_s", "count_finalize_s")
        )
        mbps = (
            f"{r['n_bases'] / 1_000_000 / counting_seconds:.2f}"
            if r["n_bases"] is not None and counting_seconds > 0
            else "---"
        )
        kmers_per_second = (
            f"{r['n_kmers'] / counting_seconds:.0f}"
            if r["n_kmers"] is not None and counting_seconds > 0
            else "---"
        )
        count_time = f"{counting_seconds:.6f}s" if counting_seconds > 0 else "---"
        pcr_seconds = r["stage_times"].get("pcr_s")
        pcr_time = f"{pcr_seconds:.6f}s" if pcr_seconds is not None else "---"
        lines.append(
            f"| {r['sample']} | {k_reads} | {wall} | {count_time} | {pcr_time} | "
            f"{_format_bytes(r['peak_mem'])} | "
            f"{_format_count(r['n_reads'])} | "
            f"{_format_count(r['n_subreads'])} | "
            f"{_format_count(r['n_bases'])} | "
            f"{_format_count(r['n_kmers'])} | {mbps} | {kmers_per_second} | "
            f"{_format_bytes(r['peak_rss'])} | {_format_count(r['table_capacity'])} | "
            f"{_format_bytes(r['temp_disk'])} |"
        )
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Benchmark combined summary (cross-panel)
# ---------------------------------------------------------------------------


def write_benchmark_summary(
    panel_results: list,
    summary_path: Path,
    panel_data_map: dict | None = None,
):
    """Write a combined benchmark summary spanning multiple panels.

    panel_results is a list of result dicts (one per panel).
    panel_data_map: optional {panel_name: panel_data} for ref availability.
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

    # One section per panel with a compact scoring matrix.
    for result in panel_results:
        panel_name = result.get("panel", "unknown")
        panel_version = result.get("panel_version", "?")
        lines.append(f"## {panel_name} v{panel_version}")
        lines.append("")
        references = result.get("provenance", {}).get("references")
        if isinstance(references, dict):
            catalog = references.get("catalog") or {}
            lines.append(
                f"Verified references: {references.get('verified_count', 0)}; "
                f"excluded/unverified: {references.get('excluded_count', 0)}; "
                f"catalog status `{catalog.get('status', 'unavailable')}`; "
                f"catalog SHA-256 `{catalog.get('sha256') or 'unavailable'}`."
            )
            lines.append("")

        samples = result.get("samples", [])
        if not samples:
            lines.append("_No samples._")
            lines.append("")
            continue

        # Get ref availability if we have panel data.
        ref_availability = {}
        reference_summary = result.get("provenance", {}).get("references")
        if isinstance(reference_summary, dict):
            ref_availability = _build_ref_availability(
                reference_summary=reference_summary
            )
        elif panel_data_map and panel_name in panel_data_map:
            ref_availability = _build_ref_availability(panel_data_map[panel_name])

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
            taxon = s.get("taxon", "")
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
                    if gr and gr.get("evaluation_status") == "not_evaluated":
                        cells.append("`N/E`")
                    elif gr and gr.get("recovered"):
                        length = gr.get("length", "?")
                        score = _score_gene(
                            True, gene, taxon,
                            gr.get("reference_match"), ref_availability,
                        )
                        cells.append(f"{length}bp `{score}`")
                    else:
                        score = _score_gene(
                            False, gene, taxon, None, ref_availability,
                        )
                        cells.append(f"`{score}`")
                lines.append(
                    f"| {label} | {k_reads} | " + " | ".join(cells) + " |"
                )
        lines.append("")

    lines.append(SCORE_LEGEND)
    lines.append("")

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Benchmark summary written to: {summary_path}")

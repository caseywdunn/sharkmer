#!/usr/bin/env python3
"""
Summarize a sharkmer parameter sweep.

Walks panels/validation_results/ for sweep result YAMLs (files whose body
contains a `sweep_label` field, written by validate_panel.py --label),
groups them by knob name and value, and produces a single Markdown report
with one section per knob:

  1. A "winners across panels" header table — one row per panel, one
     column per swept value, cells = total points across all samples in
     that panel × value cell. Points are computed from the existing
     3-position score code in sharkmer_validate.report:

       +**  recovered, ref for this species, BLAST hit same gene same species   3 pts
       +++  recovered, hit same gene in a different species                     2 pts
       +--  recovered, no references for this gene (no validation possible)     1 pt
       everything else (suspicious or not-recovered)                            0 pts

  2. Per-panel detail sub-tables — one per panel, rows = (sample, gene),
     columns = swept values, cells = the raw 3-position score string.

Panels with no `references` block in their YAML (currently bacteria and
metazoa) cannot be BLAST-validated, so their cells in the header table use
"genes recovered" as a fallback metric. This is noted in the table header.

Usage:
    python scripts/sweep_summary.py
    python scripts/sweep_summary.py --results-dir panels/validation_results
    python scripts/sweep_summary.py --output panels/validation_reports/sweep.md
"""

import argparse
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# Make the sharkmer_validate package importable.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sharkmer_validate import results, runner  # noqa: E402
from sharkmer_validate.report import _build_ref_availability, _score_gene  # noqa: E402

DEFAULT_RESULTS_DIR = runner.REPO_ROOT / "panels" / "validation_results"
DEFAULT_REPORTS_DIR = runner.REPO_ROOT / "panels" / "validation_reports"
PANELS_DIR = runner.REPO_ROOT / "panels"


# ---------------------------------------------------------------------------
# Knob registry
# ---------------------------------------------------------------------------
#
# Maps the underscore form (as it appears in sweep labels) to the canonical
# hyphenated form (as it appears in CLI flags and report headers). This
# table is the only place that needs to grow when a new knob is added to
# validate_all_panels_slurm.sh's SWEEPS array. Order does not matter as
# long as no key is a prefix of another, which is true for the current set.

KNOB_LABEL_TO_DISPLAY: dict[str, str] = {
    "k": "k",
    "max_primer_kmers": "max-primer-kmers",
    "high_coverage_ratio": "high-coverage-ratio",
    "tip_coverage_fraction": "tip-coverage-fraction",
    "node_budget_global": "node-budget-global",
    "max_node_visits": "max-node-visits",
    "max_dfs_states": "max-dfs-states",
}

# Sort key for sweep values within a knob: numeric where possible, then
# lexicographic. "0.05" < "0.1" < "0.2" works as a float; "k=17" < "k=19"
# works as an int.
def _value_sort_key(v: str):
    try:
        return (0, float(v))
    except ValueError:
        return (1, v)


# ---------------------------------------------------------------------------
# Sweep label parsing
# ---------------------------------------------------------------------------


def parse_sweep_label(label: str) -> tuple[str, str] | None:
    """Parse `sweep_<knob_underscore>_<value>` into (display_knob, value).

    Returns None for labels that don't match a known knob — those are
    skipped with a warning so a stale or hand-edited label doesn't crash
    the summary.
    """
    if not label.startswith("sweep_"):
        return None
    rest = label[len("sweep_") :]
    for underscore_knob, display_knob in KNOB_LABEL_TO_DISPLAY.items():
        prefix = underscore_knob + "_"
        if rest.startswith(prefix):
            value = rest[len(prefix) :]
            return display_knob, value
    return None


# ---------------------------------------------------------------------------
# Result discovery and grouping
# ---------------------------------------------------------------------------


def discover_sweep_results(results_dir: Path) -> list[tuple[Path, dict]]:
    """Walk results_dir for sweep result YAMLs.

    Returns [(path, loaded_result_dict)]. Files without a `sweep_label`
    field are silently ignored — those are non-sweep validation runs.
    """
    out = []
    for path in sorted(results_dir.glob("sweep_*.yaml")):
        try:
            result = results.load_result(path)
        except Exception as e:
            print(f"WARNING: failed to load {path.name}: {e}")
            continue
        if "sweep_label" not in result:
            continue
        out.append((path, result))
    return out


def group_by_knob_panel_value(
    sweep_results: list[tuple[Path, dict]],
) -> dict[str, dict[str, dict[str, dict]]]:
    """Group results by (knob -> panel -> value -> result).

    A given (knob, panel, value) triple should have exactly one result. If
    multiple files match (e.g. the sweep was re-run), the most recent one
    by file timestamp wins.
    """
    grouped: dict[str, dict[str, dict[str, tuple[Path, dict]]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    for path, result in sweep_results:
        label = result.get("sweep_label", "")
        parsed = parse_sweep_label(label)
        if parsed is None:
            print(
                f"WARNING: skipping {path.name}: could not parse sweep_label '{label}'"
            )
            continue
        knob, value = parsed
        panel = result.get("panel", "unknown")
        existing = grouped[knob][panel].get(value)
        if existing is None or path.stat().st_mtime > existing[0].stat().st_mtime:
            grouped[knob][panel][value] = (path, result)
    # Strip the path tuple, keep just the dicts.
    return {
        knob: {
            panel: {value: pair[1] for value, pair in by_value.items()}
            for panel, by_value in by_panel.items()
        }
        for knob, by_panel in grouped.items()
    }


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


SCORE_POINTS = {
    "+**": 3,  # confirmed: same gene, same species
    "+++": 2,  # right gene, different species
    "+--": 1,  # recovered, no possible validation
    # Everything else (suspicious or not recovered) is 0.
}


def points_for_score(score: str) -> int:
    return SCORE_POINTS.get(score, 0)


def load_panel_data(panel_name: str) -> dict | None:
    """Locate and load the panel YAML for `panel_name`. Returns None if not found."""
    candidate = PANELS_DIR / f"{panel_name}.yaml"
    if not candidate.exists():
        return None
    try:
        return runner.load_panel_yaml(candidate)
    except Exception:
        return None


def score_result(result: dict, panel_data: dict | None) -> list[dict]:
    """Score every (sample, gene, depth) row in a result dict.

    Returns a flat list of {sample, taxon, gene, depth, score, recovered, points}.
    """
    if panel_data is None:
        ref_availability = {}
    else:
        ref_availability = _build_ref_availability(panel_data)

    rows = []
    for sample in result.get("samples", []):
        accession = sample.get("accession", "?")
        taxon = sample.get("taxon", "")
        for depth in sample.get("depths", []):
            max_reads = depth.get("max_reads")
            for gene_entry in depth.get("genes", []):
                gene = gene_entry.get("gene", "?")
                recovered = gene_entry.get("recovered", False)
                ref_match = gene_entry.get("reference_match")
                score = _score_gene(
                    recovered=recovered,
                    gene=gene,
                    sample_taxon=taxon,
                    ref_match=ref_match,
                    ref_availability=ref_availability,
                )
                rows.append(
                    {
                        "accession": accession,
                        "taxon": taxon,
                        "gene": gene,
                        "max_reads": max_reads,
                        "score": score,
                        "recovered": recovered,
                        "points": points_for_score(score),
                    }
                )
    return rows


def has_references(panel_data: dict | None) -> bool:
    """Whether the panel YAML has any reference sequences declared."""
    if panel_data is None:
        return False
    return any(
        block.get("sequences") for block in (panel_data.get("references") or [])
    )


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------


def render_winners_table(
    knob: str,
    by_panel_value: dict[str, dict[str, dict]],
    sorted_values: list[str],
    panel_has_refs: dict[str, bool],
) -> list[str]:
    """Render the per-knob "winners across panels" table.

    Cells are total points for ref-having panels, total recovered-gene
    count for ref-less panels (with a † suffix to flag the different metric).
    """
    lines: list[str] = []
    header_cells = ["Panel"] + [f"`{knob}={v}`" for v in sorted_values]
    lines.append("| " + " | ".join(header_cells) + " |")
    lines.append("|" + "|".join(["---"] * len(header_cells)) + "|")

    panel_totals_by_value = defaultdict(int)
    n_panels_for_value = defaultdict(int)

    for panel in sorted(by_panel_value.keys()):
        by_value = by_panel_value[panel]
        ref_panel = panel_has_refs.get(panel, False)
        cells: list[str] = []
        for value in sorted_values:
            result = by_value.get(value)
            if result is None:
                cells.append("—")
                continue
            panel_data = load_panel_data(panel)
            rows = score_result(result, panel_data)
            if ref_panel:
                total = sum(r["points"] for r in rows)
                cells.append(str(total))
            else:
                total = sum(1 for r in rows if r["recovered"])
                cells.append(f"{total}†")
            panel_totals_by_value[value] += total
            n_panels_for_value[value] += 1
        lines.append(f"| {panel} | " + " | ".join(cells) + " |")

    # Footer: per-value totals across all panels.
    totals_row = ["**TOTAL**"] + [
        str(panel_totals_by_value.get(v, 0)) for v in sorted_values
    ]
    lines.append("| " + " | ".join(totals_row) + " |")

    # Winner annotation.
    if panel_totals_by_value:
        best_value = max(sorted_values, key=lambda v: panel_totals_by_value.get(v, -1))
        lines.append("")
        lines.append(
            f"Winning value (highest TOTAL): **`{knob}={best_value}`**"
            f" ({panel_totals_by_value[best_value]} points across "
            f"{n_panels_for_value[best_value]} panels)."
        )
    if any(not v for v in panel_has_refs.values()):
        lines.append("")
        lines.append(
            "† panels without a `references` block (no BLAST validation possible). "
            "Cell shows raw recovered-gene count, not on-target points; cannot be "
            "directly compared with the other rows."
        )
    return lines


def render_panel_detail_table(
    knob: str,
    panel: str,
    by_value: dict[str, dict],
    sorted_values: list[str],
    panel_data: dict | None,
) -> list[str]:
    """Render the per-(sample, gene) score grid for one (knob, panel)."""
    lines: list[str] = []

    # Build the grid: row key = (accession, taxon, gene), col key = value.
    grid: dict[tuple[str, str, str], dict[str, str]] = defaultdict(dict)
    for value in sorted_values:
        result = by_value.get(value)
        if result is None:
            continue
        rows = score_result(result, panel_data)
        for r in rows:
            key = (r["accession"], r["taxon"], r["gene"])
            grid[key][value] = r["score"]

    if not grid:
        lines.append("_No data._")
        return lines

    # Stable row ordering: sort by sample (accession), then gene.
    sorted_keys = sorted(grid.keys(), key=lambda k: (k[1] or k[0], k[2]))

    header = ["Sample", "Gene"] + [f"`{v}`" for v in sorted_values]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "|".join(["---"] * len(header)) + "|")
    for accession, taxon, gene in sorted_keys:
        cells: list[str] = []
        for value in sorted_values:
            score = grid[(accession, taxon, gene)].get(value, "—")
            cells.append(f"`{score}`")
        sample_label = taxon or accession
        lines.append(f"| {sample_label} | {gene} | " + " | ".join(cells) + " |")
    return lines


def render_summary(
    grouped: dict[str, dict[str, dict[str, dict]]],
    n_files: int,
) -> str:
    """Render the full sweep summary Markdown."""
    out: list[str] = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out.append("# Sharkmer parameter sweep summary")
    out.append("")
    out.append(f"- **Generated**: {now}")
    out.append(f"- **Result files scanned**: {n_files}")
    out.append(f"- **Knobs found**: {len(grouped)}")
    out.append("")
    out.append("## Scoring")
    out.append("")
    out.append(
        "Each (sample, gene) cell carries the same 3-position code as the "
        "per-panel validation report. The cross-panel header table converts "
        "those codes to points and sums them per (panel, knob value):"
    )
    out.append("")
    out.append("| Code | Meaning | Points |")
    out.append("|------|---------|-------:|")
    out.append("| `+**` | Recovered, confirmed: same gene, same species | 3 |")
    out.append("| `+++` | Recovered, hit same gene in a different species | 2 |")
    out.append("| `+--` | Recovered, no references for this gene (no validation possible) | 1 |")
    out.append("| anything else | Suspicious or not recovered | 0 |")
    out.append("")

    for knob in sorted(grouped.keys()):
        by_panel_value = grouped[knob]

        # Collect the union of swept values for this knob across all panels.
        all_values: set[str] = set()
        for by_value in by_panel_value.values():
            all_values.update(by_value.keys())
        sorted_values = sorted(all_values, key=_value_sort_key)

        # Look up which panels have reference DBs (for the metric note).
        panel_has_refs = {
            panel: has_references(load_panel_data(panel))
            for panel in by_panel_value.keys()
        }

        out.append(f"## Knob: `{knob}`")
        out.append("")
        out.append("### Winners across panels")
        out.append("")
        out.extend(
            render_winners_table(knob, by_panel_value, sorted_values, panel_has_refs)
        )
        out.append("")
        out.append("### Per-panel detail")
        out.append("")
        for panel in sorted(by_panel_value.keys()):
            panel_data = load_panel_data(panel)
            out.append(f"#### {panel}")
            out.append("")
            out.extend(
                render_panel_detail_table(
                    knob, panel, by_panel_value[panel], sorted_values, panel_data
                )
            )
            out.append("")
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Summarize a sharkmer parameter sweep.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"Directory containing sweep result YAMLs (default: {DEFAULT_RESULTS_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output Markdown path. Default: "
            "panels/validation_reports/sweep_summary_TIMESTAMP.md"
        ),
    )
    args = parser.parse_args()

    if not args.results_dir.exists():
        print(f"Results directory not found: {args.results_dir}")
        sys.exit(1)

    sweep_results = discover_sweep_results(args.results_dir)
    if not sweep_results:
        print(f"No sweep result YAMLs found in {args.results_dir}.")
        print("Run validate_all_panels_slurm.sh --sweep first.")
        sys.exit(1)

    print(f"Loaded {len(sweep_results)} sweep result files.")
    grouped = group_by_knob_panel_value(sweep_results)
    print(f"Knobs: {sorted(grouped.keys())}")
    md = render_summary(grouped, n_files=len(sweep_results))

    if args.output is None:
        DEFAULT_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = DEFAULT_REPORTS_DIR / f"sweep_summary_{stamp}.md"
    else:
        out_path = args.output
        out_path.parent.mkdir(parents=True, exist_ok=True)

    out_path.write_text(md)
    print(f"Sweep summary written to: {out_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Validate a sharkmer panel against its declared ENA/SRA samples.

Reads the panel's `validation.samples` block, runs sharkmer against each
sample at each declared read depth, BLASTs recovered amplicons against the
panel's gold-standard reference sequences, and emits a markdown report plus
a detailed YAML result file to panels/validation_results/.

With --genes, only runs the specified gene(s). Useful when iterating on
a single primer pair without re-running the whole panel.

Sweep mode: pass `-k`, `--extra-args`, and `--label` to forward extra
sharkmer CLI arguments and tag the resulting output files. Combined with
`--max-reads-tier highest` this is the per-cell entry point used by
scripts/validate_all_panels_slurm.sh --sweep.

Requires the sharkmer-bench conda environment (for makeblastdb/blastn).

Usage:
    python scripts/validate_panel.py panels/cnidaria.yaml
    python scripts/validate_panel.py panels/cnidaria.yaml --genes 16S CO1
    python scripts/validate_panel.py panels/cnidaria.yaml --no-blast
    python scripts/validate_panel.py panels/cnidaria.yaml \\
        -k 21 --label sweep_k_21 --max-reads-tier highest
    python scripts/validate_panel.py panels/cnidaria.yaml \\
        --extra-args "--max-primer-kmers 40" \\
        --label sweep_max_primer_kmers_40 --max-reads-tier highest
"""

import argparse
import datetime
import json
import shlex
import sys
import tempfile
from pathlib import Path

from ruamel.yaml import YAML

# Ensure sharkmer_validate package is importable.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sharkmer_validate import runner, blast_references, results, report  # noqa: E402

REPORTS_DIR = runner.REPO_ROOT / "panels" / "validation_reports"
RUNS_DIR = runner.REPO_ROOT / "panels" / "validation_runs"


# ---------------------------------------------------------------------------
# Panel loading with ruamel.yaml (for --genes filtering temp panel)
# ---------------------------------------------------------------------------


def load_panel(panel_path: Path):
    """Load a panel with ruamel.yaml for round-tripping."""
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 4096
    with open(panel_path) as f:
        data = yaml.load(f)
    return yaml, data


def build_filtered_panel(panel_path: Path, gene_filter):
    """If gene_filter is non-empty, write a temp panel YAML containing only
    those genes. Returns (path_to_use, temp_path_or_none).
    """
    if not gene_filter:
        return panel_path, None

    yaml, data = load_panel(panel_path)
    all_genes = runner.panel_gene_names(data)
    missing = [g for g in gene_filter if g not in all_genes]
    if missing:
        raise SystemExit(
            f"--genes filter references genes not in panel: {missing}\n"
            f"Available: {all_genes}"
        )
    data["primers"] = [
        p for p in data["primers"] if runner.derive_gene_name(p) in gene_filter
    ]

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    import tempfile as _tf

    tmp = _tf.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False, dir=str(RUNS_DIR)
    )
    yaml.dump(data, tmp)
    tmp.close()
    return Path(tmp.name), Path(tmp.name)


# ---------------------------------------------------------------------------
# JSON Schema validation
# ---------------------------------------------------------------------------

_SCHEMA_PATH = runner.REPO_ROOT / "schemas" / "panel" / "v2.json"


def _to_json_compat(obj):
    """Recursively convert YAML-deserialized values to JSON-compatible types.

    PyYAML parses bare dates (e.g. 2026-04-05) as datetime.date objects.
    The JSON Schema expects strings, so we convert them to ISO format strings
    before validation.
    """
    if isinstance(obj, dict):
        return {str(k): _to_json_compat(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_json_compat(v) for v in obj]
    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()
    return obj


def check_json_schema(panel_data: dict, panel_path: Path) -> None:
    """Validate panel_data against the v2 JSON Schema.

    Only runs for schema_version "2" panels. Skips silently if the jsonschema
    package is not installed or if the schema file is missing.
    Exits with a non-zero status if validation fails.
    """
    if str(panel_data.get("schema_version")) != "2":
        return

    try:
        import jsonschema
    except ImportError:
        print(
            "  [schema] jsonschema package not installed — skipping JSON Schema check.\n"
            "  Install with: conda install jsonschema  (or pip install jsonschema)"
        )
        return

    if not _SCHEMA_PATH.exists():
        print(f"  [schema] Schema file not found at {_SCHEMA_PATH} — skipping check.")
        return

    with open(_SCHEMA_PATH) as f:
        schema = json.load(f)

    # Normalize date objects to ISO strings before validation (PyYAML parses
    # bare YAML dates such as 2026-04-05 as datetime.date, not str).
    panel_normalized = _to_json_compat(panel_data)

    try:
        jsonschema.validate(instance=panel_normalized, schema=schema)
        print("  [schema] JSON Schema validation passed.")
    except jsonschema.ValidationError as e:
        location = " → ".join(str(p) for p in e.absolute_path) or "(root)"
        raise SystemExit(
            f"Panel {panel_path.name} fails JSON Schema validation:\n"
            f"  Location: {location}\n"
            f"  {e.message}"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Validate a sharkmer panel against its declared samples."
    )
    parser.add_argument("panel", help="Path to the panel YAML file")
    parser.add_argument(
        "--genes",
        nargs="+",
        help="Only validate these gene names (panel-native, without prefix).",
    )
    parser.add_argument(
        "--no-blast",
        action="store_true",
        help="Skip BLAST validation against reference sequences.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPORTS_DIR,
        help=f"Directory for markdown reports (default: {REPORTS_DIR})",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=runner.K,
        help=(
            f"Kmer length passed to sharkmer (default: {runner.K}). "
            "Override for parameter sweeps."
        ),
    )
    parser.add_argument(
        "--extra-args",
        type=str,
        default="",
        help=(
            "Additional CLI arguments forwarded to sharkmer, as a single "
            "shell-quoted string. Parsed with shlex. Example: "
            '--extra-args "--max-primer-kmers 40 --high-coverage-ratio 5.0".'
        ),
    )
    parser.add_argument(
        "--label",
        type=str,
        default=None,
        help=(
            "Sweep tag prepended to result/report filenames and stored in "
            "the result YAML so sweep_summary.py can group runs. Example: "
            "--label sweep_max_primer_kmers_40."
        ),
    )
    parser.add_argument(
        "--max-reads-tier",
        choices=["all", "highest", "lowest"],
        default="all",
        help=(
            "Which max_reads tier(s) to run per sample. Default 'all' runs "
            "every depth declared in the panel; 'highest' and 'lowest' run "
            "only the corresponding tier (used by sweep mode to reduce job "
            "count)."
        ),
    )
    args = parser.parse_args()

    panel_path = Path(args.panel).resolve()
    if not panel_path.exists():
        print(f"Panel file not found: {panel_path}")
        sys.exit(1)

    # Parse extra_args once at the top so a malformed string fails fast.
    extra_args = shlex.split(args.extra_args) if args.extra_args else []

    # Override the module-level default k so all consumers
    # (run_sharkmer cmd, results.parameters.k, primer_analysis trim cap)
    # see the swept value. validate_panel.py is the only entry point that
    # uses this code per process invocation, so the global mutation is
    # safe and confined.
    if args.k != runner.K:
        runner.K = args.k

    runner.build_sharkmer()

    # Load panel data (read-only copy for metadata).
    panel_data = runner.load_panel_yaml(panel_path)
    panel_name = panel_data.get("name")
    if not panel_name:
        print(f"Panel file missing required 'name' field: {panel_path}")
        sys.exit(1)

    check_json_schema(panel_data, panel_path)

    validation = panel_data.get("validation") or {}
    samples = validation.get("samples") or []
    if not samples:
        print(
            f"Panel '{panel_name}' has no validation.samples declared.\n"
            "Add samples to the panel file before running the validator."
        )
        sys.exit(1)

    gene_filter = list(args.genes) if args.genes else None

    # If gene filter is set, write a filtered temp panel.
    panel_path_for_run, temp_panel_path = build_filtered_panel(
        panel_path, gene_filter
    )

    try:
        sharkmer_version = runner.clean_sharkmer_version(
            runner.get_sharkmer_version()
        )
        print(f"sharkmer version: {sharkmer_version}")
        print(
            f"panel: {panel_name} v{runner.get_panel_version(panel_data)}"
        )
        print(f"samples: {len(samples)}")
        if gene_filter:
            print(f"gene filter: {', '.join(gene_filter)}")
        print()

        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = RUNS_DIR / f"{panel_name}_{stamp}"

        # Build reference BLAST DB (if references exist in panel).
        tmpdir = Path(tempfile.mkdtemp(prefix="sharkmer_refs_"))
        ref_db = None
        if not args.no_blast:
            ref_db = blast_references.build_reference_db(panel_data, tmpdir)
            if ref_db:
                print(f"Reference BLAST DB built: {ref_db}")
            else:
                refs = blast_references.extract_references(panel_data)
                if not refs:
                    print("No reference sequences in panel; skipping BLAST.")
                # else: tool not available, warning already printed
        blast_mode = "references" if ref_db else "none"

        # Run sharkmer for each sample x max_reads.
        sample_results = []
        for sample_block in samples:
            accession = sample_block["accession"]
            declared_reads = sorted(
                sample_block.get("max_reads", [1_000_000]), reverse=True
            )
            # Filter by tier. `highest` and `lowest` keep just one tier each;
            # `all` keeps every declared depth (the existing behavior).
            if args.max_reads_tier == "highest":
                max_reads_list = [declared_reads[0]]
            elif args.max_reads_tier == "lowest":
                max_reads_list = [declared_reads[-1]]
            else:
                max_reads_list = declared_reads
            print(f"=== {accession} ({len(max_reads_list)} depths) ===")
            runs = []
            for max_reads in max_reads_list:
                run = runner.run_sharkmer(
                    panel_path_for_run,
                    panel_name,
                    accession,
                    max_reads,
                    run_dir,
                    extra_args=extra_args or None,
                    k=args.k,
                )
                runs.append(run)

            # BLAST this sample's amplicons against references.
            taxon = sample_block.get("taxon", "")
            blast_references.blast_all_products(
                runs, ref_db, sample_taxon=taxon, skip_blast=args.no_blast
            )

            sample_results.append((sample_block, runs))
            print()

        # Build and write result YAML.
        result = results.build_result(
            panel_path,
            panel_data,
            sample_results,
            sharkmer_version,
            blast_mode=blast_mode,
            extra_args=extra_args or None,
            sweep_label=args.label,
        )
        result_name = results.result_filename(
            panel_data, sharkmer_version, stamp, label=args.label
        )
        result_path = results.RESULTS_DIR / result_name
        results.write_result(result, result_path)

        # Write markdown report.
        report_name = result_name.replace(".yaml", ".md")
        report_path = args.output_dir / report_name
        report.write_panel_report(
            result,
            panel_data,
            sample_results,
            report_path,
            gene_filter=gene_filter,
        )

    finally:
        if temp_panel_path is not None and temp_panel_path.exists():
            temp_panel_path.unlink()
        # Clean up temp BLAST DB directory.
        import shutil

        if tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()

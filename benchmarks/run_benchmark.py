#!/usr/bin/env python3
"""
Regression benchmark for sharkmer.

Reads benchmarks/benchmark.yaml to determine which panel + accession + depth
combinations to run. Sample metadata (taxon, taxonomy) is resolved from the
panel YAML files. Produces per-panel result YAMLs and markdown reports, plus
a combined cross-panel summary.

Usage:
    python benchmarks/run_benchmark.py
    python benchmarks/run_benchmark.py --samples Xenia_sp Agalma_elegans
    python benchmarks/run_benchmark.py --panels cnidaria insecta
    python benchmarks/run_benchmark.py --no-blast

Run from the repo root directory with the sharkmer-bench conda environment.
"""

import argparse
import hashlib
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

# Ensure sharkmer_validate package is importable.
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from sharkmer_validate import runner, blast_references, results, report  # noqa: E402

BENCHMARK_CONFIG = REPO_ROOT / "benchmarks" / "benchmark.yaml"
RUNS_DIR = REPO_ROOT / "benchmarks" / "benchmark_runs"
BENCHMARK_RESULTS_DIR = REPO_ROOT / "benchmarks" / "benchmark_results"


def load_benchmark_config(config_path: Path = BENCHMARK_CONFIG) -> list:
    """Load benchmark.yaml and return the samples list."""
    with open(config_path) as f:
        data = yaml.safe_load(f)
    return data.get("samples", [])


def resolve_sample_metadata(panel_data: dict, accession: str) -> dict | None:
    """Look up a validation sample in panel_data by accession.

    Returns the sample block (taxon, taxonomy, etc.) or None if not found.
    """
    validation = panel_data.get("validation") or {}
    for sample in validation.get("samples", []):
        if sample.get("accession") == accession:
            return sample
    return None


def evaluate_known_truth(
    run: dict,
    expected: dict | None,
    benchmark_scope: str,
    expected_input_sha256: str | None = None,
) -> dict | None:
    if not expected:
        return None
    checks = []
    if expected_input_sha256 is not None:
        observed_sha256 = run.get("input", {}).get("sha256")
        checks.append(
            {
                "check": "input_sha256",
                "expected": expected_input_sha256,
                "observed": observed_sha256,
                "passed": observed_sha256 == expected_input_sha256,
            }
        )
    if not run.get("success"):
        checks.append({"check": "run_success", "passed": False})
    elif benchmark_scope == "counting-only":
        expected_kmers = expected.get("counting-only", {}).get("n_kmers")
        if expected_kmers is not None:
            observed_kmers = run.get("run_stats", {}).get("n_kmers")
            checks.append(
                {
                    "check": "n_kmers",
                    "expected": expected_kmers,
                    "observed": observed_kmers,
                    "passed": observed_kmers == expected_kmers,
                }
            )
    else:
        observed = {
            gene["gene"]: [
                {
                    "length": product.get("length"),
                    "sha256": hashlib.sha256(product.get("sequence", "").encode()).hexdigest(),
                }
                for product in gene.get("products", [])
            ]
            for gene in run.get("genes", [])
            if gene.get("recovered")
        }
        expected_genes = expected.get("end-to-end", {}).get("genes", {})
        checks.append(
            {
                "check": "recovered_gene_set",
                "expected": sorted(expected_genes),
                "observed": sorted(observed),
                "passed": set(observed) == set(expected_genes),
            }
        )
        for gene, expected_products in expected_genes.items():
            checks.append(
                {
                    "check": f"products:{gene}",
                    "expected": expected_products,
                    "observed": observed.get(gene),
                    "passed": observed.get(gene) == expected_products,
                }
            )
    return {"passed": bool(checks) and all(check["passed"] for check in checks), "checks": checks}


def run_benchmark(
    panel_filter: list | None = None,
    sample_filter: list | None = None,
    threads: int = runner.THREADS,
    max_reads_override: list | None = None,
    run_blast: bool = True,
    executable: Path | None = None,
    benchmark_scope: str = "end-to-end",
    cache_mode: str = "warm",
    diagnostic_graphs: bool = False,
    config_path: Path = BENCHMARK_CONFIG,
    k: int = runner.K,
):
    """Run the benchmark suite."""
    executable_provenance = runner.build_sharkmer(executable=executable)
    selected_executable = Path(executable_provenance["path"])

    sharkmer_version = runner.clean_sharkmer_version(
        runner.get_sharkmer_version(selected_executable)
    )
    git_commit = runner.get_git_commit()
    machine_info = runner.get_machine_info()

    print(f"sharkmer version: {sharkmer_version}")
    print(f"git commit: {git_commit}")
    print(f"machine: {machine_info}")

    # Load benchmark config and resolve panel data.
    bench_samples = load_benchmark_config(config_path)

    # Load all referenced panels once.
    panel_cache: dict[str, tuple[Path, dict]] = {}
    for entry in bench_samples:
        panel_name = entry["panel"]
        if panel_name not in panel_cache:
            panel_path = runner.PANELS_DIR / f"{panel_name}.yaml"
            if not panel_path.exists():
                print(f"WARNING: panel file not found: {panel_path}")
                continue
            panel_cache[panel_name] = (panel_path, runner.load_panel_yaml(panel_path))

    # Filter and group by panel.
    by_panel: dict[str, list] = {}
    for entry in bench_samples:
        panel_name = entry["panel"]
        accession = entry["accession"]

        if panel_filter and panel_name not in panel_filter:
            continue

        if panel_name not in panel_cache:
            continue

        panel_path, panel_data = panel_cache[panel_name]
        sample_meta = resolve_sample_metadata(panel_data, accession)

        # Allow filtering by taxon name (underscored) or accession.
        if sample_filter:
            taxon = (sample_meta or {}).get("taxon", accession)
            sample_id = taxon.replace(" ", "_")
            if accession not in sample_filter and sample_id not in sample_filter:
                continue

        max_reads = entry.get("max_reads", [1_000_000])

        by_panel.setdefault(panel_name, []).append(
            {
                "panel_path": panel_path,
                "panel_data": panel_data,
                "accession": accession,
                "max_reads": max_reads,
                "sample_meta": sample_meta
                or {"accession": accession, "taxon": entry.get("taxon", "")},
                "notes": entry.get("notes"),
                "input": (config_path.parent / entry["input"]).resolve()
                if entry.get("input")
                else None,
                "expected": entry.get("expected"),
                "input_sha256": entry.get("input_sha256"),
            }
        )

    if not by_panel:
        print("No benchmark samples matched the filters.")
        sys.exit(1)

    total = sum(len(items) for items in by_panel.values())
    print(f"Panels: {len(by_panel)}, benchmark entries: {total}")
    print()

    stamp = runner.unique_run_id()
    all_panel_results = []
    known_truth_failures = 0
    execution_failures = 0

    for panel_name, items in sorted(by_panel.items()):
        panel_path = items[0]["panel_path"]
        panel_data = items[0]["panel_data"]
        panel_version = runner.get_panel_version(panel_data)

        print(f"=== Panel: {panel_name} v{panel_version} ===")

        # Build reference BLAST DB.
        tmpdir = Path(tempfile.mkdtemp(prefix=f"sharkmer_refs_{panel_name}_"))
        ref_db = None
        if run_blast and benchmark_scope == "end-to-end":
            ref_db = blast_references.build_reference_db(panel_data, tmpdir)
            if ref_db:
                print(f"  Reference DB built: {ref_db}")
        blast_mode = "references" if ref_db else "none"

        run_dir = RUNS_DIR / f"benchmark_{panel_name}_{stamp}"

        sample_results = []
        for item in items:
            accession = item["accession"]
            sample_meta = item["sample_meta"]
            taxon = sample_meta.get("taxon", "")
            max_reads_list = max_reads_override or item["max_reads"]
            max_reads_sorted = sorted(max_reads_list, reverse=True)

            print(
                f"\n--- {taxon or accession} ({len(max_reads_sorted)} depths) ---"
            )
            runs = []
            for max_reads in max_reads_sorted:
                run = runner.run_sharkmer(
                    panel_path,
                    panel_name,
                    accession,
                    max_reads,
                    run_dir,
                    threads=threads,
                    dump_graph=diagnostic_graphs,
                    executable=selected_executable,
                    benchmark_scope=benchmark_scope,
                    cache_mode=cache_mode,
                    input_path=item["input"],
                    k=k,
                )
                known_truth = evaluate_known_truth(
                    run, item["expected"], benchmark_scope, item["input_sha256"]
                )
                if not run.get("success", False):
                    execution_failures += 1
                if known_truth is not None:
                    run["known_truth"] = known_truth
                    if not known_truth["passed"]:
                        known_truth_failures += 1
                runs.append(run)

            # BLAST against references.
            blast_references.blast_all_products(
                runs,
                ref_db,
                sample_taxon=taxon,
                skip_blast=not run_blast,
                reference_genes={
                    reference["gene_name"]
                    for reference in blast_references.extract_references(panel_data)
                },
            )

            sample_results.append((sample_meta, runs))

        # Build and write per-panel result YAML.
        result = results.build_result(
            panel_path,
            panel_data,
            sample_results,
            sharkmer_version,
            blast_mode=blast_mode,
            machine_info=machine_info,
            executable_provenance=executable_provenance,
            evaluated_genes=runner.panel_gene_names(panel_data)
            if benchmark_scope == "end-to-end"
            else [],
            run_id=stamp,
        )
        result_name = results.result_filename(panel_data, sharkmer_version, stamp)
        result_path = BENCHMARK_RESULTS_DIR / result_name
        results.write_result(result, result_path)

        # Write per-panel markdown report.
        report_name = result_name.replace(".yaml", ".md")
        report_path = BENCHMARK_RESULTS_DIR / report_name
        report.write_panel_report(
            result, panel_data, sample_results, report_path
        )

        all_panel_results.append(result)

        shutil.rmtree(tmpdir, ignore_errors=True)
        print()

    # Write combined benchmark summary.
    if all_panel_results:
        summary_name = (
            f"benchmark_{sharkmer_version}_{git_commit}_{stamp}.summary.md"
        )
        summary_path = BENCHMARK_RESULTS_DIR / summary_name
        report.write_benchmark_summary(all_panel_results, summary_path)

    print("Benchmark complete.")
    if execution_failures:
        raise SystemExit(f"{execution_failures} benchmark execution(s) failed; artifacts preserved")
    if known_truth_failures:
        raise SystemExit(f"{known_truth_failures} known-truth benchmark run(s) failed")


def main():
    parser = argparse.ArgumentParser(
        description="Run sharkmer regression benchmark"
    )
    parser.add_argument(
        "--panels",
        nargs="+",
        help="Only run these panels (default: all in benchmark.yaml)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=BENCHMARK_CONFIG,
        help=f"Benchmark dataset definition (default: {BENCHMARK_CONFIG})",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=runner.K,
        help=f"Kmer length (default: {runner.K})",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        help="Only run these samples (by taxon name or accession)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=runner.THREADS,
        help=f"Number of threads (default: {runner.THREADS})",
    )
    parser.add_argument(
        "--max-reads",
        type=int,
        nargs="+",
        help="Override max_reads for all samples",
    )
    parser.add_argument(
        "--no-blast",
        action="store_true",
        help="Skip BLAST validation of amplicons",
    )
    parser.add_argument(
        "--executable",
        type=Path,
        help="Use and fingerprint this executable instead of building the workspace binary",
    )
    parser.add_argument(
        "--scope",
        choices=["end-to-end", "counting-only"],
        default="end-to-end",
        help="Benchmark the full pipeline or counting only (default: end-to-end)",
    )
    parser.add_argument(
        "--cache-mode",
        choices=["warm", "cold"],
        default="warm",
        help="Reuse the sharkmer read cache or use an empty per-run cache",
    )
    parser.add_argument(
        "--diagnostic-graphs",
        action="store_true",
        help="Include graph dump overhead and record it in run parameters",
    )
    args = parser.parse_args()

    run_benchmark(
        panel_filter=args.panels,
        sample_filter=args.samples,
        threads=args.threads,
        max_reads_override=args.max_reads,
        run_blast=not args.no_blast,
        executable=args.executable,
        benchmark_scope=args.scope,
        cache_mode=args.cache_mode,
        diagnostic_graphs=args.diagnostic_graphs,
        config_path=args.config.resolve(),
        k=args.k,
    )


if __name__ == "__main__":
    main()

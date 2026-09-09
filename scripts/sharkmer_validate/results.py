"""Read and write validation/benchmark result YAML files.

Both validation and benchmarks produce the same result format. Structure is
depth-first: panel -> samples -> depths -> genes.
"""

import hashlib
from datetime import datetime
from pathlib import Path

import yaml

from . import blast_references, runner

RESULTS_DIR = runner.REPO_ROOT / "panels" / "validation_results"


def build_result(
    panel_path: Path,
    panel_data: dict,
    sample_results: list,
    sharkmer_version: str,
    blast_mode: str = "none",
    machine_info: dict | None = None,
    extra_args: list | None = None,
    sweep_label: str | None = None,
    executable_provenance: dict | None = None,
    evaluated_genes: list[str] | None = None,
    run_id: str | None = None,
) -> dict:
    """Build the result dict from sample_results.

    sample_results is a list of (sample_block, runs) tuples, where each
    run is a dict returned by runner.run_sharkmer().

    `extra_args` and `sweep_label` are recorded so that the sweep summary
    script can group result files by knob/value without having to parse
    filenames. `extra_args` is the shlex-parsed list of CLI args that
    were forwarded to sharkmer; `sweep_label` is a free-form tag like
    `sweep_max_primer_kmers_40` that identifies which sweep cell this
    run belongs to.
    """
    panel_name = panel_data.get("name", "unknown")
    panel_version = runner.get_panel_version(panel_data)

    runs_flat = [run for _, runs in sample_results for run in runs]
    actual_parameter_sets = [run.get("actual_parameters", {}) for run in runs_flat]
    parameters: dict = {"runs": actual_parameter_sets}
    reference_genes = {
        reference["gene_name"] for reference in blast_references.extract_references(panel_data)
    }
    if extra_args:
        parameters["extra_args"] = list(extra_args)

    result = {
        "panel": panel_name,
        "panel_version": panel_version,
        "sharkmer_version": sharkmer_version,
        "git_commit": runner.get_git_commit(),
        "date": datetime.now().strftime("%Y-%m-%d"),
        "run_id": run_id,
        "parameters": parameters,
        "blast_mode": blast_mode,
        "rustc_version": runner.get_rustc_version(),
        "provenance": {
            "executable": executable_provenance,
            "panel": {
                "path": str(panel_path.resolve()),
                "sha256": runner._sha256_file(panel_path),
            },
            "references": blast_references.reference_checksums(panel_data),
        },
    }

    if sweep_label:
        result["sweep_label"] = sweep_label

    if machine_info is None:
        machine_info = runner.get_machine_info()
    result["machine"] = machine_info
    result["metric_availability"] = {
        "stage_times": "sharkmer_stats",
        "allocator_peak_bytes": "sharkmer_peak_alloc",
        "peak_rss_bytes": "gnu_time_linux_or_null",
        "temp_disk_peak_bytes": None,
        "temp_disk_final_bytes": "invocation_directory_scan",
        "os_page_cache_state": None,
    }

    samples = []
    for sample_block, runs in sample_results:
        accession = sample_block["accession"]
        taxon = sample_block.get("taxon", "")

        depths = []
        for run in sorted(runs, key=lambda r: r["max_reads"]):
            gene_results = []
            for prod in run.get("genes", []):
                gene_entry = {
                    "gene": prod["gene"],
                    "recovered": prod.get("recovered", bool(prod.get("products"))),
                    "length": prod["lengths"][0] if prod.get("lengths") else None,
                    "n_products": prod.get("n_products", 0),
                    "kmer_count_median": prod.get("kmer_count_median"),
                    "evaluation_status": "recovered"
                    if prod.get("recovered", bool(prod.get("products")))
                    else "no_product",
                    "reference_status": "available"
                    if prod["gene"] in reference_genes
                    else "no_reference",
                }

                if prod.get("failure_reason"):
                    gene_entry["failure_reason"] = prod["failure_reason"]

                product_entries = []
                for product in prod.get("products", []):
                    product_entry = {
                        "product_index": product.get("product_index"),
                        "length": product.get("length"),
                        "sha256": hashlib.sha256(product.get("sequence", "").encode()).hexdigest(),
                        "kmer_count_median": product.get("kmer_count_median"),
                        "reference_match": product.get("reference_match"),
                    }
                    product_entries.append(product_entry)
                if product_entries:
                    gene_entry["products"] = product_entries

                ref_match = prod.get("reference_match")
                if ref_match:
                    gene_entry["reference_match"] = ref_match

                gene_results.append(gene_entry)

            recovered_genes = {g["gene"] for g in gene_results}
            all_genes = runner.panel_gene_names(panel_data)
            for gene in all_genes:
                if gene not in recovered_genes:
                    if evaluated_genes is not None and gene not in evaluated_genes:
                        evaluation_status = "not_evaluated"
                    elif not run.get("success", False):
                        evaluation_status = "failed_run"
                    else:
                        evaluation_status = "no_product"
                    gene_results.append(
                        {
                            "gene": gene,
                            "recovered": False,
                            "evaluation_status": evaluation_status,
                            "reference_status": "available"
                            if gene in reference_genes
                            else "no_reference",
                        }
                    )

            depth_entry = {
                "max_reads": run["max_reads"],
                "wall_time_s": run.get("wall_time_s"),
                "success": run.get("success", False),
                "genes": gene_results,
                "command": run.get("command"),
                "actual_parameters": run.get("actual_parameters"),
                "input": run.get("input"),
                "logs": run.get("logs"),
                "failure": run.get("failure"),
                "peak_rss_bytes": run.get("peak_rss_bytes"),
                "peak_rss_provenance": run.get("peak_rss_provenance"),
                "temp_disk_final_bytes": run.get("temp_disk_bytes"),
                "known_truth": run.get("known_truth"),
            }
            run_stats = run.get("run_stats")
            if run_stats:
                depth_entry["run_stats"] = run_stats
            depths.append(depth_entry)

        samples.append(
            {
                "accession": accession,
                "taxon": taxon,
                "depths": depths,
            }
        )

    result["samples"] = samples
    return result


def write_result(result: dict, output_path: Path) -> Path:
    """Write a result dict to YAML."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        yaml.dump(result, f, default_flow_style=False, sort_keys=False)
    print(f"Results written to: {output_path}")
    return output_path


def load_result(path: Path) -> dict:
    """Load a result YAML file."""
    with open(path) as f:
        return yaml.safe_load(f)


def result_filename(
    panel_data: dict,
    sharkmer_version: str,
    timestamp: str,
    label: str | None = None,
) -> str:
    """Generate a result filename from panel metadata.

    If `label` is provided (e.g. `sweep_max_primer_kmers_40`), it is
    prepended to the filename so concurrent sweep runs that share a
    second-resolution timestamp do not clobber each other.
    """
    panel_name = panel_data.get("name", "unknown")
    panel_version = runner.get_panel_version(panel_data)
    safe_version = sharkmer_version.replace(" ", "_")
    base = f"{panel_name}_{panel_version}_{safe_version}_{timestamp}.yaml"
    if label:
        return f"{label}_{base}"
    return base

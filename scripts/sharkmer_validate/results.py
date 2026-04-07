"""Read and write validation/benchmark result YAML files.

Both validation and benchmarks produce the same result format. Structure is
depth-first: panel -> samples -> depths -> genes.
"""

from datetime import datetime
from pathlib import Path

import yaml

from . import runner

RESULTS_DIR = runner.REPO_ROOT / "panels" / "validation_results"


def build_result(
    panel_path: Path,
    panel_data: dict,
    sample_results: list,
    sharkmer_version: str,
    blast_mode: str = "none",
    machine_info: dict | None = None,
) -> dict:
    """Build the result dict from sample_results.

    sample_results is a list of (sample_block, runs) tuples, where each
    run is a dict returned by runner.run_sharkmer().
    """
    panel_name = panel_data.get("name", "unknown")
    panel_version = runner.get_panel_version(panel_data)

    result = {
        "panel": panel_name,
        "panel_version": panel_version,
        "sharkmer_version": sharkmer_version,
        "git_commit": runner.get_git_commit(),
        "date": datetime.now().strftime("%Y-%m-%d"),
        "parameters": {
            "k": runner.K,
            "threads": runner.THREADS,
        },
        "blast_mode": blast_mode,
        "rustc_version": runner.get_rustc_version(),
        "hash_backend": "ahashmap",  # default feature flag
        "build_profile": "release",
    }

    if machine_info is None:
        machine_info = runner.get_machine_info()
    result["machine"] = machine_info

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
                    "recovered": True,
                    "length": prod["lengths"][0] if prod.get("lengths") else None,
                    "n_products": prod.get("n_products", 0),
                    "kmer_count_median": prod.get("kmer_count_median"),
                }

                ref_match = prod.get("reference_match")
                if ref_match:
                    gene_entry["reference_match"] = ref_match

                gene_results.append(gene_entry)

            # Add entries for genes not recovered at this depth.
            recovered_genes = {g["gene"] for g in gene_results}
            all_genes = runner.panel_gene_names(panel_data)
            for gene in all_genes:
                if gene not in recovered_genes:
                    gene_results.append(
                        {
                            "gene": gene,
                            "recovered": False,
                        }
                    )

            depth_entry = {
                "max_reads": run["max_reads"],
                "wall_time_s": run.get("wall_time_s"),
                "success": run.get("success", False),
                "genes": gene_results,
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


def result_filename(panel_data: dict, sharkmer_version: str, timestamp: str) -> str:
    """Generate a result filename from panel metadata."""
    panel_name = panel_data.get("name", "unknown")
    panel_version = runner.get_panel_version(panel_data)
    safe_version = sharkmer_version.replace(" ", "_")
    return f"{panel_name}_{panel_version}_{safe_version}_{timestamp}.yaml"

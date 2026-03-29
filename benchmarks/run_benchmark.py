#!/usr/bin/env python3
"""
Regression benchmark for sharkmer.

Downloads SRA datasets (if not cached), runs sharkmer at multiple coverage
levels, and collects results into a YAML file for comparison across versions.

Usage:
    python benchmarks/run_benchmark.py [--samples SAMPLE1 SAMPLE2 ...]
    python benchmarks/run_benchmark.py --all

Run from the repo root directory.
"""

import argparse
import concurrent.futures
import glob
import os
import platform
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blast_validate import validate_results
from summarize import summarize


# --- Configuration ---

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "benchmarks" / "config.yaml"
SHARKMER_BIN = REPO_ROOT / "target" / "release" / "sharkmer"
DATA_DIR = REPO_ROOT / "benchmarks" / "data"
OUTPUT_DIR = REPO_ROOT / "benchmarks" / "output"
RESULTS_DIR = REPO_ROOT / "benchmarks" / "results"

K = 31
THREADS = 8
DEFAULT_MAX_READS = [1_000_000]


def load_config():
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    return config["sample"]


def get_max_reads_for_sample(sample_config):
    """Return the max_reads list for a sample, sorted descending."""
    reads = sample_config.get("max_reads", DEFAULT_MAX_READS)
    return sorted(reads, reverse=True)


def get_sharkmer_version():
    result = subprocess.run(
        [str(SHARKMER_BIN), "--version"],
        capture_output=True, text=True
    )
    return result.stdout.strip()


def get_git_commit():
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True
    )
    return result.stdout.strip()


def get_rustc_version():
    result = subprocess.run(
        ["rustc", "--version"],
        capture_output=True, text=True
    )
    return result.stdout.strip()


def get_machine_info():
    info = {
        "os": f"{platform.system()} {platform.release()}",
        "cpu_model": platform.processor() or platform.machine(),
        "cpu_cores": os.cpu_count(),
    }
    # Try to get total RAM
    try:
        if platform.system() == "Darwin":
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True
            )
            info["total_ram_gb"] = round(int(result.stdout.strip()) / (1024**3), 1)
        elif platform.system() == "Linux":
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal"):
                        kb = int(line.split()[1])
                        info["total_ram_gb"] = round(kb / (1024**2), 1)
                        break
    except Exception:
        info["total_ram_gb"] = None
    return info


def find_sample_data(sample_name, sample_config):
    """Find local data files for a sample.

    Data files are named {accession}.fastq in benchmarks/data/.
    Returns a list of FASTQ file paths, or None if any are missing.
    """
    accessions = sample_config.get("reads", [])
    if not accessions:
        print(f"  WARNING: No reads defined for {sample_name}, skipping")
        return None

    all_fastq_paths = []
    for accession in accessions:
        fastq_path = DATA_DIR / f"{accession}.fastq"
        if fastq_path.exists():
            all_fastq_paths.append(fastq_path)
        else:
            print(f"  ERROR: data file not found: {fastq_path}")
            print(f"    Download with: benchmarks/sra_download.sh {accession} <nreads>")
            print(f"    Then rename to {fastq_path}")
            return None

    return all_fastq_paths


def run_sharkmer(sample_name, fastq_paths, arguments, max_reads, threads=THREADS):
    """Run sharkmer and return wall time in seconds."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    k_reads = max_reads // 1000
    sample_prefix = f"{sample_name}_{k_reads}k"

    cmd = [
        str(SHARKMER_BIN),
        "-k", str(K),
        "-t", str(threads),
        "--max-reads", str(max_reads),
        "--dump-graph",
        "-o", str(OUTPUT_DIR) + "/",
        "-s", sample_prefix,
    ]

    # Parse the arguments string (e.g., "--pcr-panel cnidaria --pcr-panel bacteria")
    if arguments:
        cmd.extend(arguments.split())

    # Add input files
    for fq in fastq_paths:
        cmd.append(str(fq))

    print(f"  Running: {' '.join(cmd)}")

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall_time = time.time() - start_time

    if result.returncode != 0:
        print(f"  ERROR running sharkmer: {result.stderr[-500:]}")
        return None, wall_time

    # Save log
    log_path = OUTPUT_DIR / f"{sample_prefix}.log"
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    return sample_prefix, wall_time


def parse_stats_file(stats_path):
    """Parse the current .stats format (key\tvalue per line)."""
    stats = {}
    if not stats_path.exists():
        return stats
    with open(stats_path) as f:
        for line in f:
            parts = line.strip().split("\t", 1)
            if len(parts) == 2:
                key, value = parts
                try:
                    stats[key] = int(value)
                except ValueError:
                    stats[key] = value
    return stats


def parse_fasta_products(sample_prefix):
    """Find and parse all FASTA output files for a sample."""
    products = []
    pattern = str(OUTPUT_DIR / f"{sample_prefix}_*.fasta")
    for fasta_path in sorted(glob.glob(pattern)):
        gene_name = Path(fasta_path).stem.replace(f"{sample_prefix}_", "")
        sequences = []
        current_header = None
        current_seq = []

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append("".join(current_seq))
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header is not None:
                sequences.append("".join(current_seq))

        # Parse kmer stats from header if available
        kmer_count_median = None
        if current_header:
            median_match = re.search(r"median (\d+)", current_header)
            if median_match:
                kmer_count_median = int(median_match.group(1))

        products.append({
            "gene": gene_name,
            "n_products": len(sequences),
            "lengths": [len(s) for s in sequences],
            "kmer_count_median": kmer_count_median,
            "sequences": sequences,
        })

    return products


def collect_sample_result(sample_name, sample_config, sample_prefix, wall_time):
    """Collect all results for a single sample."""
    # Try YAML stats first (v2.0+), fall back to old TSV format
    stats_yaml_path = OUTPUT_DIR / f"{sample_prefix}.stats.yaml"
    stats_path = OUTPUT_DIR / f"{sample_prefix}.stats"
    if stats_yaml_path.exists():
        with open(stats_yaml_path) as f:
            stats = yaml.safe_load(f) or {}
    else:
        stats = parse_stats_file(stats_path)
    products = parse_fasta_products(sample_prefix)

    # Determine panel from arguments
    arguments = sample_config.get("arguments", "")
    panels = re.findall(r"--pcr-panel\s+(\S+)", arguments)
    panel = ", ".join(panels) if panels else "none"

    genes_amplified = len(products)

    result = {
        "sample": sample_name,
        "panel": panel,
        "wall_time_s": round(wall_time, 1),
        "n_reads_read": stats.get("n_reads_read"),
        "n_bases_read": stats.get("n_bases_read"),
        "n_kmers": stats.get("n_kmers"),
        "n_unique_kmers": None,  # Not in current stats format
        "genes_amplified": genes_amplified,
        "products": products,
    }

    return result


def run_benchmark(samples_to_run=None, threads=THREADS, max_reads_override=None,
                  run_blast=True):
    """Run the full benchmark suite."""
    config = load_config()

    if samples_to_run is None:
        samples_to_run = list(config.keys())

    # Build sharkmer if needed
    print("Building sharkmer (release)...")
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=str(REPO_ROOT),
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Build failed: {result.stderr}")
        sys.exit(1)

    sharkmer_version = get_sharkmer_version()
    git_commit = get_git_commit()
    machine_info = get_machine_info()

    print(f"sharkmer version: {sharkmer_version}")
    print(f"git commit: {git_commit}")
    print(f"machine: {machine_info}")
    print(f"samples: {len(samples_to_run)}")
    print()

    results = []
    for sample_name in samples_to_run:
        if sample_name not in config:
            print(f"WARNING: {sample_name} not found in config, skipping")
            continue

        sample_config = config[sample_name]

        # Find data files (one per accession, reused across sweep levels)
        fastq_paths = find_sample_data(sample_name, sample_config)
        if fastq_paths is None:
            continue

        # Determine sweep levels
        if max_reads_override is not None:
            sample_reads = sorted(max_reads_override, reverse=True)
        else:
            sample_reads = get_max_reads_for_sample(sample_config)

        for max_reads in sample_reads:
            k_reads = max_reads // 1000
            print(f"=== {sample_name} ({k_reads}k reads) ===")

            # Run sharkmer
            arguments = sample_config.get("arguments", "")
            sample_prefix, wall_time = run_sharkmer(
                sample_name, fastq_paths, arguments, max_reads, threads=threads
            )
            if sample_prefix is None:
                print(f"  FAILED, skipping result collection")
                results.append({
                    "sample": sample_name,
                    "max_reads": max_reads,
                    "status": "failed",
                    "wall_time_s": round(wall_time, 1),
                })
                continue

            # Collect results
            sample_result = collect_sample_result(
                sample_name, sample_config, sample_prefix, wall_time
            )
            sample_result["max_reads"] = max_reads
            results.append(sample_result)

            print(f"  Completed in {wall_time:.1f}s, {sample_result['genes_amplified']} genes amplified")
            print()

    # Assemble the benchmark result
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    date_str = datetime.now().strftime("%Y-%m-%d")
    safe_version = sharkmer_version.split("(")[0].strip().replace(" ", "_")
    filename = f"{date_str}_{safe_version}_{git_commit}.yaml"
    result_path = RESULTS_DIR / filename

    benchmark_result = {
        "sharkmer_version": sharkmer_version,
        "git_commit": git_commit,
        "date": date_str,
        "machine": machine_info,
        "rustc_version": get_rustc_version(),
        "hash_backend": "ahashmap",  # default feature flag
        "build_profile": "release",
        "parameters": {
            "k": K,
            "threads": threads,
        },
        "results": results,
    }

    with open(result_path, "w") as f:
        yaml.dump(benchmark_result, f, default_flow_style=False, sort_keys=False)

    print(f"Benchmark results written to: {result_path}")

    # BLAST validation
    if run_blast:
        print()
        print("Running BLAST validation...")
        try:
            validate_results(result_path)
        except Exception as e:
            print(f"BLAST validation failed: {e}")
            print("Results saved without BLAST annotations. Run manually with:")
            print(f"  python benchmarks/blast_validate.py {result_path}")
    else:
        print("Skipping BLAST validation (--no-blast). Run manually with:")
        print(f"  python benchmarks/blast_validate.py {result_path}")

    # Generate human-readable summary
    summary_path = result_path.with_suffix(".summary.md")
    try:
        summarize(result_path, output_path=summary_path)
    except Exception as e:
        print(f"Summary generation failed: {e}")

    return result_path


def main():
    parser = argparse.ArgumentParser(
        description="Run sharkmer regression benchmark"
    )
    parser.add_argument(
        "--samples", nargs="+",
        help="Specific samples to run (default: all)"
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Run all samples from config"
    )
    parser.add_argument(
        "--threads", type=int, default=THREADS,
        help=f"Number of threads (default: {THREADS})"
    )
    parser.add_argument(
        "--max-reads", type=int, nargs="+",
        help="Override max_reads for all samples (ignores per-sample config)"
    )
    parser.add_argument(
        "--no-blast", action="store_true",
        help="Skip BLAST validation of amplicons"
    )
    args = parser.parse_args()

    if args.samples:
        run_benchmark(args.samples, threads=args.threads,
                      max_reads_override=args.max_reads,
                      run_blast=not args.no_blast)
    else:
        run_benchmark(threads=args.threads,
                      max_reads_override=args.max_reads,
                      run_blast=not args.no_blast)


if __name__ == "__main__":
    main()

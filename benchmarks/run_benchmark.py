#!/usr/bin/env python3
"""
Regression benchmark for sharkmer.

Downloads SRA datasets (if not cached), runs sharkmer at fixed coverage,
and collects results into a YAML file for comparison across versions.

Usage:
    python benchmarks/run_benchmark.py [--samples SAMPLE1 SAMPLE2 ...]
    python benchmarks/run_benchmark.py --all

Run from the repo root directory.
"""

import argparse
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


# --- Configuration ---

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "benchmarks" / "config.yaml"
SHARKMER_BIN = REPO_ROOT / "sharkmer" / "target" / "release" / "sharkmer"
DATA_DIR = REPO_ROOT / "benchmarks" / "data"
OUTPUT_DIR = REPO_ROOT / "benchmarks" / "output"
RESULTS_DIR = REPO_ROOT / "benchmarks" / "results"

MAX_READS = 1_000_000
K = 31
THREADS = 8


def load_config():
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    return config["sample"]


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


def download_sample(sample_name, sample_config):
    """Download SRA data by streaming from ENA, keeping only MAX_READS reads.

    Uses benchmarks/sra_download.sh which streams ENA-hosted FASTQ files
    and stops as soon as N reads have been written, avoiding downloading
    the full run. Returns a list of FASTQ file paths (one per mate for
    paired-end, one for single-end).
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    accessions = sample_config.get("reads", [])
    if not accessions:
        print(f"  WARNING: No reads defined for {sample_name}, skipping")
        return None

    download_script = REPO_ROOT / "benchmarks" / "sra_download.sh"

    all_fastq_paths = []
    for accession in accessions:
        # Check if already downloaded
        existing = sorted(DATA_DIR.glob(f"{accession}_*_{MAX_READS}.fastq"))
        if existing:
            print(f"  Using cached data for {accession}: {[p.name for p in existing]}")
            all_fastq_paths.extend(existing)
            continue

        print(f"  Downloading {accession} ({MAX_READS} reads per mate from ENA)...")
        result = subprocess.run(
            [str(download_script), accession, str(MAX_READS)],
            capture_output=True, text=True,
            cwd=str(DATA_DIR),
        )
        if result.returncode != 0:
            print(f"  ERROR downloading {accession}:")
            print(f"    {result.stderr.strip()}")
            return None

        # Print download log
        for line in result.stderr.strip().split("\n"):
            print(f"    {line}")

        # Find the downloaded files
        downloaded = sorted(DATA_DIR.glob(f"{accession}_*_{MAX_READS}.fastq"))
        if not downloaded:
            print(f"  ERROR: no FASTQ files found after download for {accession}")
            return None

        all_fastq_paths.extend(downloaded)

    return all_fastq_paths


def run_sharkmer(sample_name, fastq_paths, arguments, threads=THREADS):
    """Run sharkmer and return wall time in seconds.

    fastq_paths is a list of FASTQ files (e.g. [R1.fastq, R2.fastq] for
    paired-end, or [reads.fastq] for single-end). All are passed as
    positional arguments to sharkmer.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_prefix = f"{sample_name}_1000k"

    cmd = [
        str(SHARKMER_BIN),
        "-k", str(K),
        "-t", str(threads),
        "--max-reads", str(MAX_READS),
        "-o", str(OUTPUT_DIR) + "/",
        "-s", sample_prefix,
    ]

    # Parse the arguments string (e.g., "--pcr cnidaria --pcr bacteria")
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
    stats_path = OUTPUT_DIR / f"{sample_prefix}.stats"
    stats = parse_stats_file(stats_path)
    products = parse_fasta_products(sample_prefix)

    # Determine panel from arguments
    arguments = sample_config.get("arguments", "")
    panels = re.findall(r"--pcr\s+(\S+)", arguments)
    panel = ", ".join(panels) if panels else "none"

    genes_amplified = len(products)
    genes_failed = 0  # Can't easily determine from current output format

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


def run_benchmark(samples_to_run=None, threads=THREADS):
    """Run the full benchmark suite."""
    config = load_config()

    if samples_to_run is None:
        samples_to_run = list(config.keys())

    # Build sharkmer if needed
    print("Building sharkmer (release)...")
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=str(REPO_ROOT / "sharkmer"),
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
        print(f"=== {sample_name} ===")

        # Download
        fastq_paths = download_sample(sample_name, sample_config)
        if fastq_paths is None:
            continue

        # Run sharkmer
        arguments = sample_config.get("arguments", "")
        sample_prefix, wall_time = run_sharkmer(
            sample_name, fastq_paths, arguments, threads=threads
        )
        if sample_prefix is None:
            print(f"  FAILED, skipping result collection")
            results.append({
                "sample": sample_name,
                "status": "failed",
                "wall_time_s": round(wall_time, 1),
            })
            continue

        # Collect results
        sample_result = collect_sample_result(
            sample_name, sample_config, sample_prefix, wall_time
        )
        results.append(sample_result)

        print(f"  Completed in {wall_time:.1f}s, {sample_result['genes_amplified']} genes amplified")
        print()

    # Assemble the benchmark result
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    date_str = datetime.now().strftime("%Y-%m-%d")
    filename = f"{date_str}_{sharkmer_version.replace(' ', '_')}_{git_commit}.yaml"
    result_path = RESULTS_DIR / filename

    benchmark_result = {
        "sharkmer_version": sharkmer_version,
        "git_commit": git_commit,
        "date": date_str,
        "machine": machine_info,
        "rustc_version": get_rustc_version(),
        "hash_backend": "fxhashmap",  # default feature flag
        "build_profile": "release",
        "parameters": {
            "max_reads": MAX_READS,
            "k": K,
            "threads": threads,
        },
        "results": results,
    }

    with open(result_path, "w") as f:
        yaml.dump(benchmark_result, f, default_flow_style=False, sort_keys=False)

    print(f"Benchmark results written to: {result_path}")
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
    args = parser.parse_args()

    if args.samples:
        run_benchmark(args.samples, threads=args.threads)
    else:
        run_benchmark(threads=args.threads)


if __name__ == "__main__":
    main()

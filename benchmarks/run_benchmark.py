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
CONFIG_PATH = REPO_ROOT / "tests" / "config.yaml"
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
    """Download SRA data using fasterq-dump if not already cached."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    fastq_path = DATA_DIR / f"{sample_name}.fastq"
    if fastq_path.exists():
        print(f"  Using cached data: {fastq_path}")
        return fastq_path

    accessions = sample_config.get("reads", [])
    if not accessions:
        print(f"  WARNING: No reads defined for {sample_name}, skipping")
        return None

    temp_dir = DATA_DIR / f"{sample_name}_temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    for accession in accessions:
        print(f"  Downloading {accession}...")
        result = subprocess.run(
            ["fasterq-dump", accession, "-O", str(temp_dir), "-e", "4"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  ERROR downloading {accession}: {result.stderr}")
            return None

    # Concatenate all fastq files
    print(f"  Concatenating FASTQ files...")
    with open(fastq_path, "w") as outfile:
        for fq in sorted(temp_dir.glob("*.fastq")):
            with open(fq) as infile:
                for line in infile:
                    outfile.write(line)

    # Clean up temp files
    for fq in temp_dir.glob("*.fastq"):
        fq.unlink()
    temp_dir.rmdir()

    return fastq_path


def run_sharkmer(sample_name, fastq_path, arguments, threads=THREADS):
    """Run sharkmer and return wall time in seconds."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_prefix = f"{sample_name}_1000k"
    output_prefix = OUTPUT_DIR / sample_prefix

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

    # Add input file
    cmd.append(str(fastq_path))

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
        fastq_path = download_sample(sample_name, sample_config)
        if fastq_path is None:
            continue

        # Run sharkmer
        arguments = sample_config.get("arguments", "")
        sample_prefix, wall_time = run_sharkmer(
            sample_name, fastq_path, arguments, threads=threads
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

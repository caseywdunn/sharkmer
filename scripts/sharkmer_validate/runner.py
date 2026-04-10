"""Run sharkmer and parse its output.

Consolidates execution logic formerly split across benchmarks/run_benchmark.py
and scripts/validate_panel.py into a single module.
"""

import glob
import os
import platform
import re
import subprocess
import time
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SHARKMER_BIN = REPO_ROOT / "target" / "release" / "sharkmer"
PANELS_DIR = REPO_ROOT / "panels"
DATA_DIR = REPO_ROOT / "benchmarks" / "data"
CACHE_DIR = REPO_ROOT / "benchmarks" / "data" / "cache"

K = 19  # Match sharkmer default
THREADS = 8
DEFAULT_MAX_READS = [1_000_000]


# ---------------------------------------------------------------------------
# Sharkmer version / git helpers
# ---------------------------------------------------------------------------


def get_sharkmer_version() -> str:
    result = subprocess.run(
        [str(SHARKMER_BIN), "--version"], capture_output=True, text=True
    )
    return result.stdout.strip()


def clean_sharkmer_version(raw: str) -> str:
    """Extract just the version number from sharkmer --version output.

    Example: "sharkmer 3.0.0-dev (https://...)" -> "3.0.0-dev"
    """
    s = raw.split("(")[0].strip()
    parts = s.split()
    if len(parts) >= 2 and parts[0].lower() == "sharkmer":
        return parts[1]
    return s or raw


def get_git_commit() -> str:
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def get_rustc_version() -> str:
    result = subprocess.run(["rustc", "--version"], capture_output=True, text=True)
    return result.stdout.strip()


def get_machine_info() -> dict:
    info = {
        "os": f"{platform.system()} {platform.release()}",
        "cpu_model": platform.processor() or platform.machine(),
        "cpu_cores": os.cpu_count(),
    }
    try:
        if platform.system() == "Darwin":
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True
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


def build_sharkmer():
    """Build sharkmer in release mode. Exits on failure."""
    if SHARKMER_BIN.exists():
        return
    print("Building sharkmer (release)...")
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"Build failed: {result.stderr}")
        raise SystemExit(1)


# ---------------------------------------------------------------------------
# Panel loading helpers
# ---------------------------------------------------------------------------


def load_panel_yaml(panel_path: Path) -> dict:
    """Load panel YAML with PyYAML (read-only, no round-tripping)."""
    with open(panel_path) as f:
        return yaml.safe_load(f)


def panel_gene_names(panel_data: dict) -> list:
    return [p["gene_name"] for p in panel_data.get("primers", [])]


def discover_panels(panels_dir: Path = None) -> list:
    """Return list of (path, data) for all panel YAML files in panels_dir."""
    if panels_dir is None:
        panels_dir = PANELS_DIR
    panels = []
    for path in sorted(panels_dir.glob("*.yaml")):
        try:
            data = load_panel_yaml(path)
            if data and data.get("primers"):
                panels.append((path, data))
        except Exception as e:
            print(f"WARNING: skipping {path}: {e}")
    return panels


def get_panel_version(panel_data: dict) -> str:
    return str(panel_data.get("version", "unversioned"))


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def parse_fasta_products(sample_prefix: str, output_dir: Path) -> list:
    """Find and parse all FASTA output files for a sample run."""
    products = []
    pattern = str(output_dir / f"{sample_prefix}_*.fasta")
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

        kmer_count_median = None
        if current_header:
            median_match = re.search(r"median (\d+)", current_header)
            if median_match:
                kmer_count_median = int(median_match.group(1))

        products.append(
            {
                "gene": gene_name,
                "n_products": len(sequences),
                "lengths": [len(s) for s in sequences],
                "kmer_count_median": kmer_count_median,
                "sequences": sequences,
            }
        )

    return products


def parse_stats_file(stats_path: Path) -> dict:
    """Parse legacy .stats format (key\\tvalue per line)."""
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


def _parse_stats_yaml(stats_path: Path) -> dict:
    """Parse the sharkmer .stats.yaml file for performance metrics."""
    if not stats_path.exists():
        return {}
    try:
        with open(stats_path) as f:
            data = yaml.safe_load(f) or {}
        return {
            "n_reads_read": data.get("n_reads_read"),
            "n_bases_read": data.get("n_bases_read"),
            "n_kmers": data.get("n_kmers"),
            "peak_memory_bytes": data.get("peak_memory_bytes"),
        }
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Sharkmer execution
# ---------------------------------------------------------------------------


def run_sharkmer(
    panel_path: Path,
    panel_name: str,
    accession: str,
    max_reads: int,
    output_dir: Path,
    threads: int = THREADS,
    dump_graph: bool = False,
    extra_args: list = None,
    k: int | None = None,
) -> dict:
    """Run sharkmer once for a (panel, accession, max_reads) combination.

    Returns a dict with:
      sample_prefix, accession, max_reads, wall_time_s, success, genes
    Gene names are returned without the panel prefix.

    `k` defaults to the module-level `K` (so existing callers are unchanged).
    Sweep callers pass an explicit value to override.
    """
    if k is None:
        k = K
    k_reads = max_reads // 1000
    sample_prefix = f"{panel_name}_{accession}_{k_reads}k"
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(SHARKMER_BIN),
        "-k",
        str(k),
        "-t",
        str(threads),
        "--max-reads",
        str(max_reads),
        "-o",
        str(output_dir) + "/",
        "-s",
        sample_prefix,
        "--pcr-panel-file",
        str(panel_path),
    ]

    if dump_graph:
        cmd.append("--dump-graph")

    if extra_args:
        cmd.extend(extra_args)

    # Prefer local FASTQ; otherwise stream via --ena.
    local_fq = DATA_DIR / f"{accession}.fastq"
    if local_fq.exists():
        cmd.append(str(local_fq))
        source = "local"
    else:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--ena", accession, "--cache-dir", str(CACHE_DIR)])
        source = "ena"

    print(f"  [{source}] running sharkmer for {accession} @ {k_reads}k reads...")
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - start

    if result.returncode != 0:
        print(f"  ERROR: sharkmer failed ({wall:.1f}s)")
        print(f"  stderr tail: {result.stderr[-500:]}")
        return {
            "sample_prefix": sample_prefix,
            "accession": accession,
            "max_reads": max_reads,
            "wall_time_s": round(wall, 1),
            "success": False,
            "genes": [],
        }

    # Persist log for provenance.
    log_path = output_dir / f"{sample_prefix}.log"
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    products = parse_fasta_products(sample_prefix, output_dir)

    # Strip panel prefix from gene names.
    prefix_to_strip = f"{panel_name}_"
    stripped = []
    for p in products:
        gene = p["gene"]
        if gene.startswith(prefix_to_strip):
            gene = gene[len(prefix_to_strip) :]
        stripped.append({**p, "gene": gene})

    # Parse stats YAML for performance data.
    stats_path = output_dir / f"{sample_prefix}.stats.yaml"
    run_stats = _parse_stats_yaml(stats_path)

    print(f"  completed in {wall:.1f}s, {len(stripped)} genes amplified")
    return {
        "sample_prefix": sample_prefix,
        "accession": accession,
        "max_reads": max_reads,
        "wall_time_s": round(wall, 1),
        "success": True,
        "genes": stripped,
        "run_stats": run_stats,
    }

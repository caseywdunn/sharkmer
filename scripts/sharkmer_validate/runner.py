"""Run sharkmer and parse its output.

Consolidates execution logic formerly split across benchmarks/run_benchmark.py
and scripts/validate_panel.py into a single module.
"""

import hashlib
import json
import os
import platform
import re
import shlex
import subprocess
import time
import uuid
from datetime import datetime, timezone
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


def unique_run_id() -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    return f"{timestamp}_{uuid.uuid4().hex[:8]}"


# ---------------------------------------------------------------------------
# Sharkmer version / git helpers
# ---------------------------------------------------------------------------


def get_sharkmer_version(executable: Path = SHARKMER_BIN) -> str:
    result = subprocess.run(
        [str(executable), "--version"], capture_output=True, text=True
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


def get_git_commit(full: bool = False) -> str:
    command = ["git", "-C", str(REPO_ROOT), "rev-parse"]
    if not full:
        command.append("--short")
    command.append("HEAD")
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as input_file:
        for block in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_tree_sha256() -> str:
    result = subprocess.run(
        [
            "git", "-C", str(REPO_ROOT), "ls-files", "--cached", "--others",
            "--exclude-standard", "-z",
        ],
        capture_output=True,
        check=True,
    )
    digest = hashlib.sha256()
    for relative_bytes in result.stdout.split(b"\0"):
        if not relative_bytes:
            continue
        relative_path = relative_bytes.decode()
        path = REPO_ROOT / relative_path
        digest.update(relative_bytes)
        digest.update(b"\0")
        if path.exists():
            with open(path, "rb") as source_file:
                for block in iter(lambda: source_file.read(1024 * 1024), b""):
                    digest.update(block)
        else:
            digest.update(b"<deleted>")
        digest.update(b"\0")
    return digest.hexdigest()


def _git_dirty() -> bool:
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--untracked-files=all"],
        capture_output=True,
        text=True,
        check=True,
    )
    return bool(result.stdout.strip())


def get_rustc_version() -> str:
    try:
        result = subprocess.run(
            ["rustc", "--version"], capture_output=True, text=True
        )
    except FileNotFoundError:
        return "unavailable"
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


def build_sharkmer(
    executable: Path | None = None,
    build_flags: list[str] | None = None,
) -> dict:
    """Build the workspace binary or fingerprint an explicitly selected binary."""
    explicit = executable is not None
    selected = Path(executable).resolve() if explicit else SHARKMER_BIN.resolve()
    cargo_flags = list(build_flags or [])
    build_command = None
    artifact_features = None
    if not explicit:
        build_command = [
            "cargo", "build", "--release", "--message-format=json-render-diagnostics",
            *cargo_flags,
        ]
        print("Building sharkmer (release)...")
        result = subprocess.run(
            build_command,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"Build failed: {result.stderr}")
            raise SystemExit(1)
        artifacts = []
        for line in result.stdout.splitlines():
            try:
                message = json.loads(line)
            except json.JSONDecodeError:
                continue
            if (
                message.get("reason") == "compiler-artifact"
                and message.get("target", {}).get("name") == "sharkmer"
                and message.get("executable")
            ):
                artifacts.append(message)
        if not artifacts:
            raise SystemExit("Cargo did not report a sharkmer executable artifact")
        artifact = artifacts[-1]
        selected = Path(artifact["executable"]).resolve()
        artifact_features = artifact.get("features", [])
    if not selected.is_file():
        raise SystemExit(f"sharkmer executable not found: {selected}")
    return {
        "path": str(selected),
        "sha256": _sha256_file(selected),
        "size_bytes": selected.stat().st_size,
        "selected_explicitly": explicit,
        "built_for_run": not explicit,
        "build_command": build_command,
        "build_flags": cargo_flags,
        "build_environment": {
            key: value
            for key, value in os.environ.items()
            if key == "RUSTFLAGS"
            or key == "CARGO_ENCODED_RUSTFLAGS"
            or key.startswith("CARGO_PROFILE_RELEASE_")
        },
        "build_profile": "release" if not explicit else None,
        "cargo_artifact_features": artifact_features,
        "hash_backend": "fxhashmap"
        if artifact_features and "fxhashmap" in artifact_features
        else ("ahashmap" if artifact_features and "ahashmap" in artifact_features else "unknown"),
        "workspace_source_observation": {
            "revision": get_git_commit(full=True),
            "dirty": _git_dirty(),
            "tree_sha256": _source_tree_sha256(),
            "relationship_to_explicit_binary": "unknown" if explicit else "built_immediately_before_fingerprint",
        },
        "fingerprinted_at": datetime.now(timezone.utc).isoformat(),
    }


# ---------------------------------------------------------------------------
# Panel loading helpers
# ---------------------------------------------------------------------------


def load_panel_yaml(panel_path: Path) -> dict:
    """Load panel YAML with PyYAML (read-only, no round-tripping)."""
    with open(panel_path) as f:
        return yaml.safe_load(f)


def derive_gene_name(primer: dict) -> str:
    """Derive the output gene name from structured primer fields.

    Mirrors the Rust ``derive_gene_name`` logic in preconfigured.rs:
    - gene only          → gene
    - gene + region      → gene-region
    - gene + index       → gene_index
    - gene + region + index → gene-region_index
    """
    gene = primer["gene"]
    region = primer.get("region")
    index = primer.get("index")
    name = f"{gene}-{region}" if region is not None else gene
    if index is not None:
        name = f"{name}_{index}"
    return name


def panel_gene_names(panel_data: dict) -> list:
    return [derive_gene_name(p) for p in panel_data.get("primers", [])]


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
    return str(panel_data.get("panel_version", "unversioned"))


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def _parse_header_fields(header: str) -> dict:
    fields = {}
    for token in shlex.split(header):
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def parse_fasta_products(
    sample_prefix: str,
    output_dir: Path,
    output_files: list[str] | None = None,
) -> list:
    """Parse FASTA products, optionally restricted to a stats manifest."""
    products = []
    if output_files is None:
        fasta_paths = sorted(output_dir.glob(f"{sample_prefix}_*.fasta"))
    else:
        fasta_paths = [output_dir / file_name for file_name in output_files]
    for fasta_path in fasta_paths:
        if not fasta_path.is_file() or fasta_path.parent.resolve() != output_dir.resolve():
            raise ValueError(f"Manifest FASTA is missing or outside run directory: {fasta_path}")
        gene_name = fasta_path.stem.replace(f"{sample_prefix}_", "", 1)
        sequence_records = []
        current_header = None
        current_seq = []

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequence = "".join(current_seq)
                        fields = _parse_header_fields(current_header)
                        sequence_records.append(
                            {
                                "header": current_header,
                                "sequence": sequence,
                                "length": len(sequence),
                                "product_index": int(fields["product"])
                                if fields.get("product", "").isdigit()
                                else len(sequence_records),
                                "kmer_count_median": int(fields["kmer_count_median"])
                                if fields.get("kmer_count_median", "").isdigit()
                                else None,
                            }
                        )
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header is not None:
                sequence = "".join(current_seq)
                fields = _parse_header_fields(current_header)
                sequence_records.append(
                    {
                        "header": current_header,
                        "sequence": sequence,
                        "length": len(sequence),
                        "product_index": int(fields["product"])
                        if fields.get("product", "").isdigit()
                        else len(sequence_records),
                        "kmer_count_median": int(fields["kmer_count_median"])
                        if fields.get("kmer_count_median", "").isdigit()
                        else None,
                    }
                )

        products.append(
            {
                "gene": gene_name,
                "n_products": len(sequence_records),
                "lengths": [record["length"] for record in sequence_records],
                "kmer_count_median": sequence_records[0]["kmer_count_median"]
                if sequence_records
                else None,
                "sequences": [record["sequence"] for record in sequence_records],
                "products": sequence_records,
                "output_file": fasta_path.name,
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
    """Parse the complete sharkmer stats manifest."""
    if not stats_path.exists():
        raise ValueError(f"Current-run stats manifest is missing: {stats_path}")
    try:
        with open(stats_path) as f:
            data = yaml.safe_load(f) or {}
    except Exception as error:
        raise ValueError(f"Current-run stats manifest is malformed: {error}") from error
    if not isinstance(data.get("pcr_results", []), list):
        raise ValueError("Current-run stats manifest has invalid pcr_results")
    return data


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _cache_manifest(cache_dir: Path) -> list[dict]:
    entries = []
    if not cache_dir.exists():
        return entries
    for meta_path in sorted(cache_dir.glob("*.meta.yaml")):
        try:
            metadata = yaml.safe_load(meta_path.read_text()) or {}
        except Exception as error:
            entries.append({"metadata_file": meta_path.name, "error": str(error)})
            continue
        entries.append(
            {
                "metadata_file": meta_path.name,
                "data_file": meta_path.name.replace(".meta.yaml", ".fastq.gz"),
                "url": metadata.get("url"),
                "sha256": metadata.get("sha256"),
                "n_reads": metadata.get("n_reads"),
                "complete": metadata.get("complete"),
            }
        )
    return entries


def _run_with_rss(command: list[str], rss_path: Path) -> tuple[subprocess.CompletedProcess, int | None, str]:
    if platform.system() == "Linux" and Path("/usr/bin/time").is_file():
        wrapped = ["/usr/bin/time", "-v", "-o", str(rss_path), *command]
        completed = subprocess.run(
            wrapped,
            capture_output=True,
            text=True,
            env={**os.environ, "LC_ALL": "C"},
        )
        rss_kib = None
        if rss_path.exists():
            match = re.search(
                r"Maximum resident set size \(kbytes\):\s*(\d+)",
                rss_path.read_text(),
            )
            if match:
                rss_kib = int(match.group(1))
        return completed, rss_kib * 1024 if rss_kib is not None else None, "proc_time_linux"
    return subprocess.run(command, capture_output=True, text=True), None, "unavailable"


def _option_value(command: list[str], option: str) -> str | None:
    positions = [index for index, value in enumerate(command) if value == option]
    if not positions:
        return None
    position = positions[-1]
    return command[position + 1] if position + 1 < len(command) else None


def has_execution_failures(sample_results: list) -> bool:
    return any(not run.get("success", False) for _, runs in sample_results for run in runs)


def _validate_extra_args(extra_args: list[str]):
    reserved = {
        "-k", "-t", "--threads", "--max-reads", "-o", "--outdir", "-s", "--sample",
        "--pcr-panel", "--pcr-panel-file", "--pcr-primers", "--ena", "--cache-dir",
        "--no-cache", "--paired",
    }
    for argument in extra_args:
        lower = argument.lower()
        if argument in reserved:
            raise ValueError(f"--extra-args cannot override provenance option {argument}")
        if (
            argument == "-"
            or lower.startswith(("http://", "https://", "ftp://"))
            or lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
            or Path(argument).exists()
        ):
            raise ValueError(f"--extra-args cannot add an untracked input source: {argument}")


def _validate_stats_manifest(
    stats: dict,
    sample_prefix: str,
    expected_k: int,
    expected_genes: set[str],
    benchmark_scope: str,
    expected_command: list[str] | None = None,
):
    if stats.get("sample") != sample_prefix:
        raise ValueError("Stats sample does not match this invocation")
    if stats.get("kmer_length") != expected_k:
        raise ValueError("Stats kmer_length does not match the effective command")
    if expected_command is not None and stats.get("command") != " ".join(expected_command):
        raise ValueError("Stats command does not match the current invocation")
    pcr_results = stats.get("pcr_results", [])
    if not isinstance(pcr_results, list):
        raise ValueError("Stats pcr_results is missing or invalid")
    if benchmark_scope == "end-to-end" and "pcr_results" not in stats:
        raise ValueError("End-to-end stats manifest is missing pcr_results")
    observed_genes = {entry.get("gene_name") for entry in pcr_results}
    if benchmark_scope == "end-to-end" and observed_genes != expected_genes:
        raise ValueError(
            f"Stats gene set mismatch: expected {sorted(expected_genes)}, got {sorted(observed_genes)}"
        )
    if benchmark_scope == "counting-only" and pcr_results:
        raise ValueError("Counting-only run unexpectedly reported PCR results")
    output_files = []
    for entry in pcr_results:
        status = entry.get("status")
        if status not in {"success", "fail"}:
            raise ValueError(f"Invalid PCR status for {entry.get('gene_name')}: {status}")
        if status == "success":
            if not entry.get("output_file") or not entry.get("n_products"):
                raise ValueError(f"Successful gene lacks output manifest data: {entry.get('gene_name')}")
            output_files.append(entry["output_file"])
        elif entry.get("n_products") != 0:
            raise ValueError(f"Failed gene reports products: {entry.get('gene_name')}")
    if len(output_files) != len(set(output_files)):
        raise ValueError("Stats manifest contains duplicate FASTA paths")


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
    executable: Path = SHARKMER_BIN,
    benchmark_scope: str = "end-to-end",
    cache_mode: str = "warm",
    input_path: Path | None = None,
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
    invocation_dir = output_dir / f"{sample_prefix}_{uuid.uuid4().hex[:12]}"
    invocation_dir.mkdir()

    cmd = [
        str(Path(executable).resolve()),
        "-k",
        str(k),
        "-t",
        str(threads),
        "--max-reads",
        str(max_reads),
        "-o",
        str(invocation_dir) + "/",
        "-s",
        sample_prefix,
    ]

    if benchmark_scope == "end-to-end":
        cmd.extend(["--pcr-panel-file", str(panel_path)])
    elif benchmark_scope != "counting-only":
        raise ValueError(f"Unknown benchmark scope: {benchmark_scope}")

    if dump_graph and benchmark_scope == "end-to-end":
        cmd.append("--dump-graph")

    if extra_args:
        _validate_extra_args(extra_args)
        cmd.extend(extra_args)

    # Prefer local FASTQ; otherwise stream via --ena.
    local_candidates = [
        Path(input_path) if input_path else DATA_DIR / f"{accession}.fastq",
        DATA_DIR / f"{accession}.fastq.gz",
    ]
    local_fq = next((candidate.resolve() for candidate in local_candidates if candidate.exists()), None)
    if local_fq is not None:
        cmd.append(str(local_fq))
        source = "local"
        input_provenance = {
            "source": source,
            "path": str(local_fq.resolve()),
            "sha256": _sha256_file(local_fq),
            "size_bytes": local_fq.stat().st_size,
            "subset": {"method": "first_reads", "max_reads": max_reads},
        }
        cache_dir = CACHE_DIR
        cache_before = []
    else:
        cache_dir = (
            invocation_dir / "cold_cache" if cache_mode == "cold" else CACHE_DIR
        )
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_before = _cache_manifest(cache_dir)
        cmd.extend(["--ena", accession, "--cache-dir", str(cache_dir)])
        source = "ena"
        input_provenance = {
            "source": source,
            "accession": accession,
            "subset": {"method": "first_reads", "max_reads": max_reads},
            "cache_mode": cache_mode,
            "cache_state_before": cache_before,
        }

    print(f"  [{source}] running sharkmer for {accession} @ {k_reads}k reads...")
    start = time.monotonic()
    rss_path = invocation_dir / "resource_usage.txt"
    result, peak_rss_bytes, rss_method = _run_with_rss(cmd, rss_path)
    wall = time.monotonic() - start

    stdout_path = invocation_dir / f"{sample_prefix}.stdout.log"
    stderr_path = invocation_dir / f"{sample_prefix}.stderr.log"
    stdout_path.write_text(result.stdout)
    stderr_path.write_text(result.stderr)
    if source == "ena":
        input_provenance["cache_state_after"] = _cache_manifest(cache_dir)

    common = {
        "sample_prefix": sample_prefix,
        "accession": accession,
        "max_reads": max_reads,
        "wall_time_s": wall,
        "invocation_dir": str(invocation_dir),
        "command": cmd,
        "actual_parameters": {
            "k": int(_option_value(cmd, "-k") or k),
            "threads": int(_option_value(cmd, "-t") or threads),
            "max_reads": int(_option_value(cmd, "--max-reads") or max_reads),
            "options": cmd[1:],
            "benchmark_scope": benchmark_scope,
            "graph_dump": dump_graph and benchmark_scope == "end-to-end",
            "cache_mode": cache_mode,
            "effective_panel": {
                "path": str(panel_path.resolve()),
                "sha256": _sha256_file(panel_path),
            }
            if benchmark_scope == "end-to-end"
            else None,
        },
        "input": input_provenance,
        "logs": {
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
            "resource_usage": str(rss_path) if rss_path.exists() else None,
        },
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_provenance": rss_method,
    }

    if result.returncode != 0:
        print(f"  ERROR: sharkmer failed ({wall:.1f}s)")
        print(f"  stderr tail: {result.stderr[-500:]}")
        return {
            **common,
            "success": False,
            "genes": [],
            "failure": {
                "kind": "failed_run",
                "returncode": result.returncode,
                "message": result.stderr[-500:],
            },
            "temp_disk_bytes": _directory_size(invocation_dir),
        }

    stats_path = invocation_dir / f"{sample_prefix}.stats.yaml"
    try:
        run_stats = _parse_stats_yaml(stats_path)
        effective_k = int(_option_value(cmd, "-k") or k)
        expected_genes = {
            f"{panel_name}_{gene}" for gene in panel_gene_names(load_panel_yaml(panel_path))
        }
        _validate_stats_manifest(
            run_stats,
            sample_prefix,
            effective_k,
            expected_genes,
            benchmark_scope,
            expected_command=cmd,
        )
        input_provenance["observed_n_reads"] = run_stats.get("n_reads_read")
        input_provenance["observed_n_bases"] = run_stats.get("n_bases_read")
        source_plan = run_stats.get("input_source")
        if not isinstance(source_plan, dict):
            raise ValueError("Stats manifest lacks actual input-source provenance")
        input_provenance["observed_source_plan"] = source_plan
        if source == "local":
            observed_inputs = [
                str(Path(path).resolve())
                for path in source_plan.get("inputs", [])
            ]
            if observed_inputs != [str(local_fq)]:
                raise ValueError("Stats input source does not match the fingerprinted local input")
        else:
            observed_paths = [Path(path).resolve() for path in source_plan.get("inputs", [])]
            if source_plan.get("kind") != "cached_remote" or (
                run_stats.get("n_reads_read", 0) > 0 and not observed_paths
            ):
                raise ValueError("Stats remote input-source provenance is missing or unexpected")
            cache_entries = {
                entry.get("data_file"): entry
                for entry in input_provenance.get("cache_state_after", [])
                if entry.get("data_file")
            }
            selected_entries = []
            for observed_path in observed_paths:
                entry = cache_entries.get(observed_path.name)
                if entry is None or not entry.get("sha256") or not observed_path.is_file():
                    raise ValueError(
                        f"No checksummed cache provenance for consumed input {observed_path}"
                    )
                observed_sha256 = _sha256_file(observed_path)
                if observed_sha256 != entry["sha256"]:
                    raise ValueError(f"Consumed cache input checksum mismatch: {observed_path}")
                selected_entries.append({**entry, "path": str(observed_path)})
            input_provenance["selected_cache_entries"] = selected_entries
        output_files = [
            gene_result["output_file"]
            for gene_result in run_stats.get("pcr_results", [])
            if gene_result.get("status") == "success"
        ]
        products = parse_fasta_products(sample_prefix, invocation_dir, output_files)
    except (KeyError, TypeError, ValueError) as error:
        print(f"  ERROR: sharkmer output manifest invalid: {error}")
        return {
            **common,
            "success": False,
            "genes": [],
            "failure": {"kind": "invalid_manifest", "message": str(error)},
            "temp_disk_bytes": _directory_size(invocation_dir),
        }

    prefix_to_strip = f"{panel_name}_"
    stripped = []
    products_by_file = {product["output_file"]: product for product in products}
    for gene_result in run_stats.get("pcr_results", []):
        gene = gene_result.get("gene_name", "")
        if gene.startswith(prefix_to_strip):
            gene = gene[len(prefix_to_strip) :]
        if gene_result.get("status") == "success":
            output_file = gene_result.get("output_file")
            if not output_file or output_file not in products_by_file:
                error = f"Successful gene {gene} has no current-run FASTA in manifest"
                return {
                    **common,
                    "success": False,
                    "genes": [],
                    "failure": {"kind": "invalid_manifest", "message": error},
                    "temp_disk_bytes": _directory_size(invocation_dir),
                }
            stripped.append({**products_by_file[output_file], "gene": gene, "recovered": True})
            parsed = stripped[-1]
            expected_count = gene_result.get("n_products")
            expected_lengths = gene_result.get("product_lengths", [])
            if parsed["n_products"] != expected_count or parsed["lengths"] != expected_lengths:
                error = f"FASTA content disagrees with stats manifest for gene {gene}"
                return {
                    **common,
                    "success": False,
                    "genes": [],
                    "failure": {"kind": "invalid_manifest", "message": error},
                    "temp_disk_bytes": _directory_size(invocation_dir),
                }
        else:
            stripped.append(
                {
                    "gene": gene,
                    "recovered": False,
                    "n_products": 0,
                    "products": [],
                    "sequences": [],
                    "lengths": [],
                    "failure_reason": gene_result.get("failure_reason"),
                }
            )

    amplified = sum(1 for gene in stripped if gene.get("recovered"))
    print(f"  completed in {wall:.1f}s, {amplified} genes amplified")
    return {
        **common,
        "success": True,
        "genes": stripped,
        "run_stats": run_stats,
        "temp_disk_bytes": _directory_size(invocation_dir),
    }

#!/usr/bin/env python3
"""Frozen finite prerelease checks; no release, upload, or acceptance action."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
import tomllib

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PLANT = ROOT / "benchmarks/benchmark_results/v3.2-plant-assay-correction-20261010"
LAUNCHER = ROOT / "benchmarks/benchmark_results/v3.2-uncertain-products-20260927/launch_bounded.py"
BENCH = Path("/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin")
LIMIT = 12 * 1024**3


def require(value, message):
    if not value:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as stream:
        while data := stream.read(1024 * 1024):
            h.update(data)
    return {"path": str(path.resolve()), "sha256": h.hexdigest(), "size_bytes": path.stat().st_size}


def dump(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def verify(row):
    require(fingerprint(row["path"]) == row, "Changed frozen file: " + row["path"])


def source_inventory():
    # Include tracked source and new nonignored source, without traversing
    # ignored historical panel-validation outputs or Python caches.
    names = subprocess.check_output(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z", "--",
        "src", "scripts", "panels", "schemas", "tests", "Cargo.toml", "Cargo.lock", "build.rs", "benchmarks/known_truth.yaml"], cwd=ROOT).decode().split("\0")
    selected = {ROOT / name for name in names if name and "__pycache__" not in Path(name).parts and not name.endswith(".pyc")}
    return {str(path.relative_to(ROOT)): fingerprint(path) for path in sorted(selected)}


def normalize(value):
    if isinstance(value, dict):
        return {str(k): normalize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [normalize(v) for v in value]
    if isinstance(value, (datetime.date, datetime.datetime)):
        return value.isoformat()
    return value


def helper_gate(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    if args.gate == "panels":
        import importlib.metadata
        import jsonschema
        schema = json.loads((ROOT / "schemas/panel/v2.json").read_text())
        validator = jsonschema.validators.validator_for(schema)
        validator.check_schema(schema)
        panels = sorted((ROOT / "panels").glob("*.yaml"))
        require(len(panels) == 9, "Need all nine built-in panels")
        rows = []
        fixture = str(ROOT / "tests/fixtures/ERR571460_100k_R1.fastq.gz")
        for panel in panels + [ROOT / "panels/examples/reference.yaml"]:
            data = normalize(yaml.safe_load(panel.read_text()))
            validator(schema).validate(data)
            row = {"panel": fingerprint(panel), "schema": "pass", "commands": []}
            if panel in panels:
                command = [str(args.binary), "--export-panel", panel.stem]
                run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
                (output / (panel.stem + ".export.yaml")).write_text(run.stdout)
                (output / (panel.stem + ".export.stderr")).write_text(run.stderr)
                require(run.returncode == 0 and normalize(yaml.safe_load(run.stdout)) == data, "Embedded/exported panel differs: " + panel.stem)
                row["commands"].append({"command": command, "returncode": run.returncode})
                modes = [["--pcr-panel", panel.stem], ["--pcr-panel-file", str(panel)]]
            else:
                modes = [["--pcr-panel-file", str(panel)]]
            for index, mode in enumerate(modes):
                product_path = output / f"{panel.stem}-{index}-dry-output"
                command = [str(args.binary), "-k", "19", "-t", "2", "--chunks", "0", "-s", "schema_probe", "-o", str(product_path), "--dry-run", *mode, fixture]
                run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=30)
                (output / f"{panel.stem}-{index}.stdout").write_text(run.stdout)
                (output / f"{panel.stem}-{index}.stderr").write_text(run.stderr)
                require(run.returncode == 0 and not any(product_path.rglob("*")), "Panel loader/dry run failed or wrote outputs: " + str(panel))
                row["commands"].append({"command": command, "returncode": run.returncode})
            rows.append(row)
        dump(output / "verification.json", {"status": "pass", "jsonschema_version": importlib.metadata.version("jsonschema"), "panels": rows,
             "builtin_exports": 9, "schema_checks": 10, "loader_dry_runs": 19, "binary": fingerprint(args.binary)})
    elif args.gate == "references":
        command = [sys.executable, str(ROOT / "scripts/audit_panel_references.py"), "--verify-existing"]
        run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=300)
        (output / "stdout.json").write_text(run.stdout)
        (output / "stderr.log").write_text(run.stderr)
        require(run.returncode == 0, "Reference audit failed")
        result = json.loads(run.stdout)
        require(result["panel_count"] == 10 and result["total_references"] == result["verified_references"] == 158 and result["failure_count"] == result["excluded_references"] == 0, "Expected all 158 active references to verify")
        dump(output / "verification.json", {"status": "pass", "command": command, "audit": result})
    elif args.gate == "fixtures":
        truth = yaml.safe_load((ROOT / "benchmarks/known_truth.yaml").read_text())
        require(len(truth["samples"]) == 1, "Expected one pinned fixture")
        oracle = truth["samples"][0]
        require(oracle["panel"] == "cnidaria" and oracle["accession"] == "ERR571460", "Wrong fixture oracle")
        fixture = ROOT / "tests/fixtures/ERR571460_100k_R1.fastq.gz"
        require(fingerprint(fixture)["sha256"] == oracle["input_sha256"], "Fixture differs")
        genes = oracle["expected"]["end-to-end"]["genes"]
        require(set(genes) == {"18S", "28S_2"} and all(len(rows) == 1 for rows in genes.values()), "Expected two single-record fixture products")
        expected = {"cnidaria_" + gene: (rows[0]["length"], rows[0]["sha256"]) for gene, rows in genes.items()}
        rows = []
        for k in (19, 31):
            sample = "fixture_k" + str(k)
            directory = output / sample
            directory.mkdir()
            command = [str(args.binary), "-k", str(k), "-t", "2", "--chunks", "0", "-s", sample, "-o", str(directory), str(fixture)]
            if k == 31:
                command += ["--pcr-panel", "cnidaria"]
            run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=300)
            (output / (sample + ".stdout")).write_text(run.stdout)
            (output / (sample + ".stderr")).write_text(run.stderr)
            require(run.returncode == 0, "Fixture CLI failed")
            stats = yaml.safe_load((directory / (sample + ".stats.yaml")).read_text())
            require(stats["n_reads_read"] == 100000 and stats["run_status"] == "complete", "Fixture ingestion or completion differs")
            if k == 19:
                require(stats["n_kmers"] == oracle["expected"]["counting-only"]["n_kmers"], "k19 count oracle failed")
            standards, records = {}, []
            for path in sorted(directory.glob("*.fasta")):
                gene = path.stem.removeprefix(sample + "_")
                for block in path.read_text().split(">"):
                    if not block:
                        continue
                    header, *sequence = block.splitlines()
                    sequence = "".join(sequence)
                    row = {"gene": gene, "length": len(sequence), "sha256": hashlib.sha256(sequence.encode()).hexdigest(), "header": header}
                    records.append(row)
                    if "assembly_status=standard" in header:
                        require(gene not in standards, "Unexpected duplicate standard fixture product")
                        standards[gene] = (row["length"], row["sha256"])
            if k == 31:
                require(standards == expected, "k31 full standard product identities differ")
            rows.append({"command": command, "counts": {key: stats[key] for key in ("n_reads_read", "n_bases_read", "n_subreads_ingested", "n_bases_ingested", "n_kmers")}, "all_records": records})
        dump(output / "verification.json", {"status": "pass", "fixture": fingerprint(fixture), "oracle": fingerprint(ROOT / "benchmarks/known_truth.yaml"), "runs": rows,
             "scope": "Regression count/full-standard-sequence oracle. Any uncertain records are listed separately, not used as biological truth."})
    print(json.dumps({"status": "pass", "gate": args.gate, "output": str(output)}))


def freeze(path):
    bindings = {Path(__file__).resolve(), HERE / "run_final_checks.sbatch", LAUNCHER,
                ROOT / "tmp/uncertain-products-20260927/toolchain/activate.sh", PLANT / "run_comparison.py", PLANT / "protocol-v3.json"}
    bindings.update(p for p in (PLANT / "proposal-v2").rglob("*") if p.is_file())
    for folder in ("v3.2-inclusive-release-20260929", "v3.2-uncertain-products-20260927"):
        bindings.update((ROOT / "benchmarks/benchmark_results" / folder).glob("*.py"))
    bindings.add(ROOT / "benchmarks/benchmark_results/v3.2-release-comparison-20260909/tools/postprocess_legacy.py")
    bindings.update([BENCH / "python", BENCH / "blastn", BENCH / "makeblastdb", Path("/usr/bin/time"), Path("/usr/bin/prlimit")])
    compiler_bin = ROOT / "tmp/uncertain-products-20260927/toolchain/rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin"
    bindings.update(compiler_bin / name for name in ("cargo", "rustc", "rustfmt", "cargo-fmt", "clippy-driver", "cargo-clippy"))
    for config_root in (ROOT / ".cargo", ROOT / "tmp/uncertain-products-20260927/toolchain/cargo"):
        bindings.update(path for path in (config_root / "config", config_root / "config.toml") if path.is_file())
    bindings.update(path for path in (ROOT / "rust-toolchain", ROOT / "rust-toolchain.toml") if path.is_file())
    deps = ROOT / "tmp/prerelease-checks-20261010/python-deps"
    bindings.update(p for p in deps.rglob("*") if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    lock = tomllib.loads((ROOT / "Cargo.lock").read_text())
    packages = {(row["name"], row["version"]) for row in lock["package"]}
    require(("rand", "0.8.6") in packages and ("rustls-webpki", "0.103.13") in packages, "Patched dependency versions absent")
    require(sys.version_info[:2] == (3, 14), "Use sharkmer-bench Python 3.14 for rpds ABI")
    import jsonschema
    import importlib.metadata
    require(importlib.metadata.version("jsonschema") == "4.26.0", "Expected pinned JSON Schema installation")
    dump(path, {"status": "frozen_for_review", "source_files": source_inventory(), "bindings": [fingerprint(p) for p in sorted(bindings)],
         "dependency_versions": {"rand": "0.8.6", "rustls-webpki": "0.103.13", "jsonschema": "4.26.0"},
         "resources": {"cpus": 2, "physical_cores": 2, "slurm_memory_gib": 16, "slurm_seconds": 7200, "child_address_space_bytes": LIMIT},
         "plant_output": "results-v1/plant", "gate_order": ["build_default", "plant_eight", "test_default", "test_fxhash", "clippy_default", "clippy_fxhash", "fmt", "python_all", "references", "panels", "fixtures"],
         "fixture_scope": "k19 count and k31 counts/full standard hashes from benchmarks/known_truth.yaml; Rust integration also tests independent Jellyfish histogram oracle on both backends.",
         "scope": "Finite local quality and assay checks; does not authorize release or establish held-out truth/performance acceptance."})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--plan", type=Path, default=HERE / "plan-v2.json")
    parser.add_argument("--output", type=Path, default=HERE / "results-v1")
    parser.add_argument("--gate", choices=("references", "panels", "fixtures"))
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    if args.freeze:
        freeze(args.freeze.resolve())
        print(json.dumps({"status": "frozen_for_review", "plan": fingerprint(args.freeze)}))
        return 0
    if args.gate:
        helper_gate(args)
        return 0
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "failed", "gates": [], "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "started": time.time(), "plan": fingerprint(args.plan)}
    try:
        require(os.environ.get("SLURM_JOB_ID") and os.environ.get("SLURM_CPUS_PER_TASK") == "2", "Need dedicated two-CPU Slurm allocation")
        plan = json.loads(args.plan.read_text())
        require(source_inventory() == plan["source_files"], "Production source changed since review")
        for row in plan["bindings"]:
            verify(row)
        selected, cores, topology = [], set(), []
        for cpu in sorted(os.sched_getaffinity(0)):
            top = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
            pair = tuple(int((top / field).read_text()) for field in ("physical_package_id", "core_id"))
            topology.append({"cpu": cpu, "socket": pair[0], "core": pair[1]})
            if pair not in cores and len(selected) < 2:
                selected.append(cpu)
                cores.add(pair)
        require(len(selected) == 2, "Need two distinct physical cores")
        os.sched_setaffinity(0, selected)
        require(Path(sys.executable).resolve() == (BENCH / "python").resolve(), "Wrong Python/ABI")
        require(all(Path(shutil.which(name)).resolve() == (BENCH / name).resolve() for name in ("blastn", "makeblastdb")), "BLAST tools missing from required environment")
        result.update(machine=platform.uname()._asdict(), cpus=selected, topology=topology, python=sys.version,
                      build_environment={key: os.environ.get(key) for key in ("CARGO_HOME", "RUSTUP_HOME", "CARGO_TARGET_DIR", "CARGO_BUILD_JOBS", "RUST_TEST_THREADS", "RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "RUSTC_WRAPPER", "PYTHONPATH")},
                      tool_versions={name: subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip() for name, command in
                        {"cargo": ["cargo", "--version"], "rustc": ["rustc", "--version"], "blastn": ["blastn", "-version"], "makeblastdb": ["makeblastdb", "-version"]}.items()})
        dump(output / "sources.before.json", plan["source_files"])
        dump(output / "preflight.json", {k: v for k, v in result.items() if k != "gates"})
        binary = output / "bin/sharkmer"
        commands = [
            ("build_default", ["cargo", "build", "--release", "--locked", "--offline", "-j", "2"], 1200),
            ("plant_eight", [sys.executable, str(PLANT / "run_comparison.py"), "--binary", str(binary), "--scope", "all", "--output", str(output / "plant")], 3600),
            ("test_default", ["cargo", "test", "--release", "--locked", "--offline", "-j", "2"], 1200),
            ("test_fxhash", ["cargo", "test", "--release", "--locked", "--offline", "-j", "2", "--no-default-features", "--features", "fxhashmap"], 1200),
            ("clippy_default", ["cargo", "clippy", "--all-targets", "--locked", "--offline", "-j", "2", "--", "-D", "warnings"], 1200),
            ("clippy_fxhash", ["cargo", "clippy", "--all-targets", "--locked", "--offline", "-j", "2", "--no-default-features", "--features", "fxhashmap", "--", "-D", "warnings"], 1200),
            ("fmt", ["cargo", "fmt", "--all", "--check"], 120),
            ("python_all", [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py", "-v"], 600),
            *[(name, [sys.executable, str(Path(__file__).resolve()), "--gate", name, "--binary", str(binary), "--output", str(output / name)], 900 if name == "fixtures" else 600) for name in ("references", "panels", "fixtures")]]
        require([row[0] for row in commands] == plan["gate_order"], "Gate order differs")
        for name, command, cap in commands:
            remaining = int(7000 - (time.time() - result["started"]))
            require(remaining > 0, "Finite suite wall budget exhausted")
            gate = output / "gates" / name
            gate.mkdir(parents=True)
            row = {"name": name, "command": command, "status": "failed", "timeout_seconds": min(cap, remaining)}
            result["gates"].append(row)
            timed = ["/usr/bin/time", "-o", str(gate / "time.json"), "-f", '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}', "--", *command]
            launch = ["/usr/bin/prlimit", "--as=" + str(LIMIT), "--", sys.executable, str(LAUNCHER), "--receipt", str(gate / "constraints.json"),
                      "--expected-cpus", ",".join(map(str, selected)), "--expected-as", str(LIMIT), "--", *timed]
            row["launch_command"] = launch
            start = time.monotonic()
            with (gate / "stdout.log").open("x") as stdout, (gate / "stderr.log").open("x") as stderr:
                child = subprocess.Popen(launch, cwd=ROOT, stdout=stdout, stderr=stderr, start_new_session=True)
                try:
                    row["returncode"] = child.wait(timeout=row["timeout_seconds"])
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                    row.update(timed_out=True, returncode=child.returncode)
            row["observer_wall_seconds"] = time.monotonic() - start
            try:
                if (gate / "time.json").exists():
                    row["metrics"] = json.loads((gate / "time.json").read_text())
                if row["returncode"] == 0:
                    constraints = json.loads((gate / "constraints.json").read_text())
                    require(constraints["status"] == "pass" and constraints["affinity"] == selected and constraints["rlimit_as_soft_hard"] == [LIMIT, LIMIT], "Actual gate constraints differ")
                    if name == "build_default":
                        binary.parent.mkdir()
                        shutil.copy2(Path(os.environ["CARGO_TARGET_DIR"]) / "release/sharkmer", binary)
                        result["frozen_binary"] = fingerprint(binary)
                        dump(output / "frozen_binary.json", {"binary": result["frozen_binary"], "build_command": command, "default_features": True, "source_plan": result["plan"]})
                    if name == "python_all":
                        log = (gate / "stderr.log").read_text()
                        require(re.search(r"test_real_database_build_and_query_binds_logical_gene_metadata.*\.\.\. ok", log), "Real BLAST test did not run successfully")
                        require("skipped=" not in log and re.search(r"Ran [1-9][0-9]* tests", log), "Python suite incomplete or skipped tests")
                    row["status"] = "pass"
            except Exception as error:
                row["validation_error"] = repr(error)
            row["files"] = [fingerprint(p) for p in sorted(gate.iterdir()) if p.is_file()]
            dump(gate / "receipt.json", row)
            print(json.dumps({"gate": name, "status": row["status"], "wall_seconds": row.get("metrics", {}).get("wall_seconds")}), flush=True)
            require(row["status"] == "pass", "Gate failed: " + name)
            if binary.exists():
                verify(result["frozen_binary"])
        require(len(result["gates"]) == 11, "Incomplete finite gate suite")
        result["status"] = "pass"
    except Exception as error:
        result["error"] = repr(error)
    finally:
        after = source_inventory()
        dump(output / "sources.after.json", after)
        result["source_files_unchanged"] = "plan" in locals() and after == plan["source_files"]
        if not result["source_files_unchanged"]:
            result["status"] = "failed"
        if "plan" in locals():
            try:
                for row in plan["bindings"]:
                    verify(row)
                result["bindings_unchanged"] = True
            except Exception as error:
                result.update(status="failed", bindings_unchanged=False, binding_error=repr(error))
        result["finished"] = time.time()
        if "plan" in locals():
            result["not_run_gates"] = [name for name in plan["gate_order"] if name not in {row["name"] for row in result["gates"]}]
        dump(output / "result.json", result)
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())

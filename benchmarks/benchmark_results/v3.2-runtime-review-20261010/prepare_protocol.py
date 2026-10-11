#!/usr/bin/env python3
"""Freeze exact sources and deterministic schedule; does not execute Sharkmer."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location("diagnostic", HERE / "run_human_diagnostic.py")
diagnostic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostic)


def main():
    source_paths = {HERE / "run_human_diagnostic.py", HERE / "run_human_diagnostic.sbatch", Path(__file__).resolve(),
                    ROOT / "tmp/inclusive-release-20260929/build-v1/baseline/build.json",
                    ROOT / "tmp/inclusive-release-20260929/build-v1/baseline/snapshot.json",
                    ROOT / "tmp/sequence-retention-20261001/checks-v3/binary.sha256",
                    ROOT / "tmp/sequence-retention-20261001/checks-v3/sources.before.sha256",
                    ROOT / "tmp/sequence-retention-20261001/checks-v3/sources.after.sha256",
                    ROOT / "tmp/sequence-retention-20261001/checks-v3/build.log",
                    Path("/usr/bin/time"), Path("/usr/bin/prlimit"), Path("/usr/bin/strace"), Path("/usr/bin/findmnt")}
    expected, binaries = {}, {}
    receipts = {
        "release": next((ROOT / "tmp/inclusive-release-20260929/comparison-v1").glob("*SRR17535371_*pair1_baseline.json")),
        "retention": next((ROOT / "tmp/sequence-retention-20261001/full-v3").glob("*SRR17535371_*pair1_candidate.json"))}
    known_binary_hashes = {"release": "e751e5eeca9150495e8d74d95bf27e4109466cad4ccf862e470fe271351d5899",
                           "retention": "5e18783ce66a4aa9a971fde78c7c620e9c204ea50c5cda656e808bff47e85c81"}
    for role, path in receipts.items():
        source_paths.add(path)
        data = json.loads(path.read_text())
        diagnostic.require(data["status"] == "pass", "Historical receipt failed")
        command = data["invocation"]["command"]
        binaries[role] = diagnostic.fingerprint(command[0])
        diagnostic.require(binaries[role]["sha256"] == known_binary_hashes[role], "Frozen binary differs")
        expected[role] = {"counts": data["outputs"]["counts"],
                          "records": sorted([{k: r[k] for k in ("gene", "sha256", "length", "assembly_status")} for r in data["outputs"]["records"]], key=lambda r: (r["gene"], r["sha256"])),
                          "targets": sorted((r["gene_name"], r["status"], r["n_products"]) for r in data["outputs"]["pcr_results"])}
        source_paths.add(Path(command[0]))
        current_input = diagnostic.fingerprint(command[-1])
        current_panel = diagnostic.fingerprint(command[command.index("--pcr-panel-file") + 1])
        current_prefix = data["prefix_before"]
        if role == "release":
            frozen_input, frozen_panel, frozen_prefix = current_input, current_panel, current_prefix
        else:
            diagnostic.require((current_input, current_panel, current_prefix) == (frozen_input, frozen_panel, frozen_prefix), "Historical input/panel differs")
        # Exercise output parsing against the actual archived results, with no
        # assembly. This catches legacy/current schema mistakes before Slurm.
        fixture_protocol = {"expected_outputs": expected, "input": frozen_input}
        diagnostic.capture_outputs(Path(data["invocation"]["outputs"]), role, fixture_protocol, command)
    diagnostic.require(expected["release"]["records"] == expected["retention"]["records"], "Human records differ across roles")
    diagnostic.inspect_input(frozen_input["path"], frozen_prefix)
    launcher = diagnostic.fingerprint(ROOT / "benchmarks/benchmark_results/v3.2-uncertain-products-20260927/launch_bounded.py")
    source_paths.update(Path(row["path"]) for row in (frozen_input, frozen_panel, launcher))
    blocks = [
        [("nfs", "release"), ("nfs", "retention"), ("local", "retention"), ("local", "release")],
        [("local", "release"), ("local", "retention"), ("nfs", "retention"), ("nfs", "release")],
        [("nfs", "retention"), ("nfs", "release"), ("local", "release"), ("local", "retention")]]
    schedule = []
    for block, rows in enumerate(blocks, 1):
        for location, role in rows:
            schedule.append({"id": f"{len(schedule) + 1:02d}_block{block}_{location}_{role}", "block": block,
                             "location": location, "role": role, "traced": False})
    for location, role in [("nfs", "release"), ("local", "retention"), ("nfs", "retention"), ("local", "release")]:
        schedule.append({"id": f"{len(schedule) + 1:02d}_traced_{location}_{role}", "block": None,
                         "location": location, "role": role, "traced": True})
    protocol = {"status": "frozen_for_review", "design": "Human 2 binaries x 2 output locations; three blocks; then four separate traced diagnostics",
        "scope": "Investigate whether the current environment reproduces output-location dependence. Does not establish the cause of historical waits or close the other three flags.",
        "lineage": {"release": "v3.1.0 commit 5a664680c91ad59f32b8c2a847b8fef37f34a0ae, frozen release build",
                    "retention": "Frozen checks-v3 executable used in full-v3 retention comparison; precedes the 2026-10-10 panel/dependency corrections",
                    "interpretation": "Whole release-to-retention comparison, not retention-only ablation. Original Human flag was release-to-earlier-inclusive."},
        "input": frozen_input, "expected_input": frozen_prefix, "panel": frozen_panel, "binaries": binaries,
        "launcher": launcher, "expected_outputs": expected, "schedule": schedule,
        "resources": {"slurm_cpus": 2, "distinct_physical_cores": True, "slurm_memory_gib": 16,
                      "process_address_space_bytes": diagnostic.ADDRESS_SPACE, "slurm_timeout_seconds": 900,
                      "normal_run_timeout_seconds": 60, "traced_run_timeout_seconds": 120, "internal_start_budget_seconds": 780},
        "cache_policy": "No assembly warm-ups and no cache dropping. Validate and read complete input, executable and panel before every invocation outside timer. Warm-input intent; no claim of guaranteed resident pages or cold filesystem state.",
        "timing_boundary": "Normal GNU time wraps Sharkmer only. Constraint check, hashing, output parsing and all artifact copying are outside. Log/time files always node-local. Only product output location varies between location cells.",
        "trace_syscalls": "fsync,fdatasync,rename,renameat,renameat2,openat,close,write,pwrite64,writev,newfstatat,unlink,unlinkat,mkdir,mkdirat",
        "trace_policy": "Four extra runs, after normal measurements. Trace publication/metadata syscalls with per-call duration and decoded paths. No trace result enters timing medians; concurrent syscall sums are not critical-path time.",
        "preservation": "Copy all local outputs/logs/traces to each persistent run directory after timing, including failed runs. Preserve NFS outputs in place. No deletion of evidence.",
        "analysis": ["Report each untraced wall/user/system/RSS and available stage measurement; absolute and percent paired deltas plus medians/ranges.",
                     "Require all five full product identities, nine target outcomes and ingestion counts to match frozen per-role receipts.",
                     "Report mount types/options, node/CPU topology, process constraints and exact producer/source hashes.",
                     "Describe fsync/publication traces as diagnostics for this run, not proof of the 2026-09-29 cause.",
                     "Human is short; timing quantization and fixed publication costs can dominate percentages. Three repeats do not establish equivalence or acceptance."],
        "sources": [diagnostic.fingerprint(path) for path in sorted(source_paths)]}
    diagnostic.write_new(HERE / "protocol-v1.json", protocol)
    diagnostic.write_new(HERE / "protocol_preflight.json", {"status": "pass", "human_input": frozen_prefix,
        "normal_runs": 12, "traced_runs": 4, "fixture_output_parsers_pass": ["release", "retention"],
        "no_sharkmer_executions": True, "protocol": diagnostic.fingerprint(HERE / "protocol-v1.json")})
    print(json.dumps({"status": "pass", "protocol": diagnostic.fingerprint(HERE / "protocol-v1.json")}))


if __name__ == "__main__":
    main()

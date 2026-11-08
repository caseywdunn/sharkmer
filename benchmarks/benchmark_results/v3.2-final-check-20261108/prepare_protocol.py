#!/usr/bin/env python3
"""Freeze existing inputs/oracles and current source without running Sharkmer."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import yaml

from run_check import fingerprint, require, verify, write_new, ROOT, HERE

PERFORMANCE = ["SRR17535371", "SRR8617500", "SRR25378184", "SRR1057608"]
CORRECTED = ["SRR25378184", "ERR14009273"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "inputs/input_inventory.json")
    parser.add_argument("--output", type=Path, default=HERE / "protocol-v1.json")
    args = parser.parse_args()
    inventory = json.loads(args.inventory.read_text())
    require(inventory["schema_version"] == 1 and len(inventory["samples"]) == 12,
            "Incomplete input/oracle inventory")
    candidate = inventory["candidate"]
    checked = json.loads(Path(candidate["validation_receipt"]["path"]).read_text())
    require(checked["status"] == "passed" and checked["source_unchanged"]
            and checked["source_before"] == checked["source_after"] == candidate["source_bindings"]
            and checked["candidate_sha256"] == candidate["binary"]["sha256"],
            "Candidate does not match the frozen passing source validation")
    require(all(row["exit_code"] == 0 for row in checked["gates"])
            and len(checked["gates"]) == 8, "Candidate validation gates incomplete")
    require(candidate["binary"]["sha256"] == "7ef94fedafec49300146c0acf52086915fa78a82a8786dc11413e399a339dd9a",
            "Unexpected final candidate executable")
    baseline = inventory["baseline"]
    build = json.loads(Path(baseline["build"]["path"]).read_text())
    snapshot = json.loads(Path(baseline["snapshot"]["path"]).read_text())
    require(build["status"] == "complete" and build["returncode"] == 0
            and build["binary"] == baseline["binary"]
            and snapshot["provenance"]["commit"] == baseline["commit"]
            and build["source_before"] == build["source_after"] == snapshot["files"]
            and set(build["artifact"]["features"]) == {"default", "ahashmap"},
            "Released executable does not match the saved clean build")
    bindings = dict(inventory["source_bindings"])
    sources = subprocess.check_output([
        "git", "ls-files", "-z", "--", "src", "scripts", "panels", "schemas",
        "Cargo.toml", "Cargo.lock", "build.rs"], cwd=ROOT).decode().split("\0")
    for name in sources:
        if name:
            bindings[name] = fingerprint(ROOT / name)
    for name, digest in candidate["source_bindings"].items():
        row = fingerprint(ROOT / name)
        require(row["sha256"] == digest, "Current source differs from validated candidate: " + name)
        bindings[name] = row
    launcher = ROOT / "benchmarks/benchmark_results/v3.2-uncertain-products-20260927/launch_bounded.py"
    required = [Path(__file__), HERE / "run_check.py", HERE / "check_controls.py",
                HERE / "prepare_inputs.py", HERE / "review_results.py",
                HERE / "controls.json", args.inventory, launcher,
                Path(sys.executable), Path("/usr/bin/time"), Path("/usr/bin/prlimit"),
                Path("/usr/bin/findmnt")]
    required.extend(path for path in Path(yaml.__file__).parent.rglob("*")
                    if path.is_file() and path.suffix in (".py", ".so"))
    for filename in ("README.md", "batch.sh"):
        path = HERE / filename
        require(path.exists(), "Protocol prose and batch wrapper must exist before freezing")
        required.append(path)
    for path in required:
        bindings[str(path)] = fingerprint(path)
    for role in (baseline, candidate):
        for key in ("binary", "build", "snapshot", "validation_receipt"):
            if key in role:
                bindings[role[key]["path"]] = role[key]
    for row in bindings.values():
        verify(row)
    for sample in inventory["samples"].values():
        verify(sample["panel"])
        bindings[sample["panel"]["path"]] = sample["panel"]
        # The dedicated compute job checks full FASTQ hashes and content.
        require(Path(sample["input"]["path"]).stat().st_size == sample["input"]["size_bytes"],
                "Frozen input size changed")
    schedule = []
    for block in (1, 2, 3):
        for sample in PERFORMANCE:
            roles = ("baseline", "candidate") if block != 2 else ("candidate", "baseline")
            for role in roles:
                schedule.append({"id": f"{len(schedule) + 1:02d}_{sample}_pair{block}_{role}",
                                 "sample": sample, "role": role, "pair_index": block,
                                 "purpose": "performance"})
    for key in inventory["samples"]:
        if key not in PERFORMANCE and not key.endswith("_corrected"):
            schedule.append({"id": f"{len(schedule) + 1:02d}_{key}_inventory", "sample": key,
                             "role": "candidate", "pair_index": None, "purpose": "inventory"})
    for key in CORRECTED:
        schedule.append({"id": f"{len(schedule) + 1:02d}_{key}_corrected", "sample": key + "_corrected",
                         "role": "candidate", "pair_index": None, "purpose": "corrected_plant"})
    require(len(schedule) == 32 and len({row["id"] for row in schedule}) == 32, "Invalid schedule")
    unique = {row["path"]: row for row in bindings.values()}
    protocol = {
        "schema_version": 1,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_scope": "Current production, panels, schemas, validators, and saved validation source inventory; docs-only edits do not alter this binary.",
        "candidate_build_provenance": candidate, "baseline_build_provenance": baseline,
        "input_inventory": fingerprint(args.inventory), "samples": inventory["samples"],
        "binaries": {"baseline": baseline["binary"], "candidate": candidate["binary"]},
        "versions": {"baseline": "3.1.0", "candidate": "3.2.0-dev"},
        "launcher": fingerprint(launcher), "bindings": [unique[key] for key in sorted(unique)],
        "address_space_bytes": 40 * 1024**3, "invocation_timeout_seconds": 900,
        "total_start_budget_seconds": 6900, "performance_samples": PERFORMANCE,
        "corrected_replaces": CORRECTED, "schedule": schedule,
        "cache_policy": "No assembly warmups/cache dropping. Entire input, panel and executable are verified outside each timing; warm-input intent, no claim of resident pages or controlled OS cache.",
        "output_location": "Fresh shared-NFS directories under one result root; both roles use the same mount and two physical cores.",
        "review_flags": {"median_wall_percent_above": 5, "median_peak_rss_percent_above": 2,
                         "meaning": "Investigation priorities, not statistical-significance tests or automatic acceptance thresholds."},
        "failure_policy": "Any identity/class/order/count/provenance mismatch, timeout or nonzero exit fails the check; preserve partial receipts and do not change the frozen oracle to match outputs.",
    }
    write_new(args.output, protocol)
    print(json.dumps(fingerprint(args.output)), flush=True)


if __name__ == "__main__":
    main()

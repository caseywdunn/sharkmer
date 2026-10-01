#!/usr/bin/env python3
"""Replay fixed inputs against old inclusive and sequence-first binaries."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OLD = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929"
sys.path.insert(0, str(OLD))
spec = importlib.util.spec_from_file_location("inclusive_comparison", OLD / "run_release_comparison.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def dump(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def record_keys(rows, label):
    return {(row["gene"], row["sha256"]) for row in rows if row["assembly_status"] == label}


def assess(jobs, full):
    by_sample = {}
    for job in jobs:
        by_sample.setdefault(job["invocation"]["input"], {}).setdefault(job["invocation"]["version"], []).append(job)
    cells = []
    for accession, roles in by_sample.items():
        old, new = (roles[role][0]["outputs"] for role in ("baseline", "candidate"))
        old_standard, new_standard = (record_keys(rows["records"], "standard") for rows in (old, new))
        old_uncertain, new_uncertain = (record_keys(rows["records"], "uncertain") for rows in (old, new))
        require(old_standard == new_standard, f"{accession}: changed standard identities")
        require(old_uncertain <= new_uncertain, f"{accession}: lost previously included uncertain identity")
        require(old["counts"] == new["counts"], f"{accession}: changed input or k-mer counts")
        if full:
            require(len(roles["baseline"]) == len(roles["candidate"]) == 3, "Incomplete repetitions")
        gains = sorted(new_uncertain - old_uncertain)
        cell = {"accession": accession, "standard": len(new_standard), "old_uncertain": len(old_uncertain),
                "new_uncertain": len(new_uncertain), "uncertain_gains": [
                    {"gene": gene, "sha256": digest,
                     "length": next(row["length"] for row in new["records"] if row["gene"] == gene and row["sha256"] == digest)}
                    for gene, digest in gains], "resources": {}}
        for role in ("baseline", "candidate"):
            role_jobs = sorted(roles[role], key=lambda job: job["invocation"]["pair_index"])
            for job in role_jobs[1:]:
                require(job["outputs"]["records"] == role_jobs[0]["outputs"]["records"]
                        and job["outputs"]["counts"] == role_jobs[0]["outputs"]["counts"],
                        f"{accession}/{role}: unstable output")
            cell["resources"][role] = {
                metric: [job["execution"]["metrics"][metric] for job in role_jobs]
                for metric in ("wall_seconds", "peak_rss_bytes")}
        cells.append(cell)
    result = {"status": "pass", "cells": cells,
              "totals": {"standard": sum(cell["standard"] for cell in cells),
                         "old_uncertain": sum(cell["old_uncertain"] for cell in cells),
                         "new_uncertain": sum(cell["new_uncertain"] for cell in cells),
                         "gained_uncertain": sum(len(cell["uncertain_gains"]) for cell in cells)}}
    if "SRR31887760" in by_sample:
        cell = next(cell for cell in cells if cell["accession"] == "SRR31887760")
        historical = json.loads((OLD / "review/dros12s_outcome_review.json").read_text())
        recovered = []
        for row in historical["historical_exact_matches"]:
            if row["length"] not in (487, 499):
                continue
            require({"gene": row["gene"], "sha256": row["sha256"], "length": row["length"]}
                    in cell["uncertain_gains"], f"Drosophila {row['length']} bp sequence not recovered exactly")
            recovered.append(row)
        result["dros12s_recovered"] = recovered
        assessment = next(entry["uncertain_products"]
                          for entry in by_sample["SRR31887760"]["candidate"][0]["outputs"]["pcr_results"]
                          if entry["gene_name"] == "insecta_12S")
        result["dros12s_retention"] = {"retained_path_records": assessment["retained_path_records"],
                                      "unique_products": assessment["unique_products"],
                                      "marker_detail_truncated": assessment["marker_detail_truncated"],
                                      "retention_truncated": assessment["retention_truncated"],
                                      "thresholds": [row["retention"] for row in assessment["thresholds"]
                                                     if row.get("retention") is not None]}
    if full:
        require(len(cells) == 10 and len(jobs) == 60, "Incomplete ten-sample comparison")
        result["paired"] = comparison.summarize(jobs)
        result["resource_flags"] = [
            {"accession": cell["accession"],
             "wall_percent": cell["resource_changes"]["wall_seconds"]["percent"],
             "rss_percent": cell["resource_changes"]["peak_rss_bytes"]["percent"]}
            for cell in result["paired"]["cells"]]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("dros", "plant", "full"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--old-binary", type=Path, required=True)
    parser.add_argument("--new-binary", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    state = {"status": "failed", "mode": args.mode, "jobs": [], "slurm_job_id": os.getenv("SLURM_JOB_ID")}
    try:
        cpus, topology = comparison.physical_cpus()
        inventory = json.loads((OLD / "inputs/inventory.json").read_text())
        available = json.loads((ROOT / "tmp/inclusive-release-20260929/availability.json").read_text())
        require(available["status"] == "pass" and available["complete_inventory_selected"], "Input inventory incomplete")
        samples = {row["accession"]: row for row in available["samples"]}
        require(len(samples) == 10, "Expected ten verified inputs")
        for row in inventory["samples"]:
            sample = samples[row["accession"]]
            require(sample["status"] == "verified" and sample["expected_prefix"] == row["expected_prefix"]
                    and sample["panel_id"] == row["panel"], "Input receipt mismatch")
        historical = json.loads((ROOT / "benchmarks/benchmark_results/v3.2-high-copy-final-20260911/provenance.json").read_text())
        require(inventory["historical_schedule"] == historical["schedule"], "Comparison schedule changed")
        panels = {row["name"]: row["assembly_panel"] for row in inventory["panels"]}
        scopes = {panel: {panel + "_" + row["gene"]: row["scope"] for row in rows}
                  for panel, rows in historical["protocol"]["target_metadata"].items()}
        binaries = {"baseline": comparison.fingerprint(args.old_binary.resolve()),
                    "candidate": comparison.fingerprint(args.new_binary.resolve())}
        prior_build = json.loads((ROOT / "tmp/inclusive-release-20260929/build-v1/candidate/build.json").read_text())
        require(binaries["baseline"] == prior_build["binary"], "Old inclusive binary differs from verified build")
        expected_new = (ROOT / "tmp/sequence-retention-20261001/checks-v3/binary.sha256").read_text().split()[0]
        require(binaries["candidate"]["sha256"] == expected_new, "New binary differs from checked build")
        bindings = {str(path.resolve()): comparison.fingerprint(path.resolve()) for path in
                    (Path(__file__), OLD / "run_release_comparison.py", OLD / "inputs/inventory.json",
                     ROOT / "tmp/inclusive-release-20260929/availability.json",
                     ROOT / "tmp/sequence-retention-20261001/checks-v3/binary.sha256",
                     ROOT / "dev_docs/SEQUENCE_RETENTION_PROTOCOL.md")}
        for row in list(panels.values()) + list(binaries.values()):
            comparison.verify(row)
            bindings[row["path"]] = row
        state.update(bindings=bindings, cpus=cpus, topology=topology,
                     machine=platform.uname()._asdict())
        focused = {"dros": "SRR31887760", "plant": "SRR25378184"}
        schedule = [row for row in inventory["historical_schedule"]
                    if args.mode == "full" or (row["input"] == focused[args.mode]
                                               and row["pair_index"] == 1)]
        require(len(schedule) == (60 if args.mode == "full" else 2), "Wrong schedule size")
        plan = []
        for original in schedule:
            invocation = dict(original)
            invocation["sample_prefix"] = f"{original['panel']}_{original['input']}_{original['depth']}"
            directory = output / "runs" / original["invocation_id"]
            invocation.update(directory=str(directory), outputs=str(directory / "outputs"))
            invocation["command"] = [binaries[original["version"]]["path"], "-k", "19", "-t", "2", "--chunks", "0",
                                     "--max-reads", "1000000", "-o", invocation["outputs"] + "/",
                                     "-s", invocation["sample_prefix"], "--pcr-panel-file",
                                     panels[original["panel"]]["path"], samples[original["input"]]["available_path"]]
            plan.append(invocation)
        dump(output / "plan.json", {"schedule": plan, "bindings": bindings, "cpus": cpus,
                                    "topology": topology, "mode": args.mode})
        state["plan"] = comparison.fingerprint(output / "plan.json")
        environment = dict(os.environ, LC_ALL="C")
        for invocation in plan:
            sample = samples[invocation["input"]]
            job = {"status": "failed", "invocation": invocation}
            state["jobs"].append(job)
            try:
                input_path = Path(sample["available_path"])
                job["prefix_before"] = comparison.inspect_prefix(input_path, sample["expected_prefix"])
                job["input_file"] = comparison.fingerprint(input_path)
                comparison.verify(binaries[invocation["version"]])
                comparison.verify(panels[invocation["panel"]])
                job["execution"] = comparison.run_timed(invocation, cpus, environment)
                require(job["execution"]["status"] == "pass", "Timed invocation failed")
                modern = dict(invocation, version="candidate")
                job["outputs"] = comparison.capture(Path(invocation["outputs"]), modern, sample,
                                                      invocation["command"], scopes[invocation["panel"]])
                comparison.verify(job["input_file"])
                comparison.verify(binaries[invocation["version"]])
                comparison.verify(panels[invocation["panel"]])
                job["status"] = "pass"
            finally:
                receipt = output / (invocation["invocation_id"] + ".json")
                dump(receipt, job)
                job["receipt"] = comparison.fingerprint(receipt)
            print(json.dumps({"completed": len(state["jobs"]), "invocation": invocation["invocation_id"],
                              "status": job["status"]}), flush=True)
        result = assess(state["jobs"], args.mode == "full")
        dump(output / "summary.json", result)
        state["summary"] = comparison.fingerprint(output / "summary.json")
        state["status"] = "pass"
    except Exception as error:
        state["error"] = repr(error)
        raise
    finally:
        dump(output / "result.json", state)


if __name__ == "__main__":
    main()

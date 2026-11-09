#!/usr/bin/env python3
"""Compare old and corrected diagnostic-only output on fixed Drosophila input."""

import importlib.util
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT / "benchmarks/benchmark_results/v3.2-inclusive-release-20260929"
sys.path.insert(0, str(OLD))
spec = importlib.util.spec_from_file_location("inclusive_comparison", OLD / "run_release_comparison.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def main():
    output = ROOT / "tmp/sequence-retention-20261001/optout-v1"
    output.mkdir(parents=True, exist_ok=False)
    available = json.loads((ROOT / "tmp/inclusive-release-20260929/availability.json").read_text())
    sample = next(row for row in available["samples"] if row["accession"] == "SRR31887760")
    inventory = json.loads((OLD / "inputs/inventory.json").read_text())
    panel = next(row["assembly_panel"] for row in inventory["panels"] if row["name"] == "insecta")
    historical = json.loads((ROOT / "benchmarks/benchmark_results/v3.2-high-copy-final-20260911/provenance.json").read_text())
    scopes = {"insecta_" + row["gene"]: row["scope"] for row in historical["protocol"]["target_metadata"]["insecta"]}
    binaries = {
        "old": ROOT / "tmp/inclusive-release-20260929/build-v1/candidate/sharkmer",
        "new": ROOT / "tmp/sequence-retention-20261001/checks-v3/sharkmer",
    }
    cpus, topology = comparison.physical_cpus()
    result = {"status": "failed", "slurm_job_id": os.getenv("SLURM_JOB_ID"),
              "cpus": cpus, "topology": topology, "jobs": {},
              "input": comparison.fingerprint(Path(sample["available_path"])),
              "panel": comparison.fingerprint(Path(panel["path"])),
              "binaries": {role: comparison.fingerprint(path) for role, path in binaries.items()}}
    try:
        for role, binary in binaries.items():
            directory = output / role
            name = "insecta_SRR31887760_1000000_optout"
            invocation = {"version": "candidate", "sample_prefix": name,
                          "directory": str(directory), "outputs": str(directory / "outputs")}
            invocation["command"] = [str(binary), "-k", "19", "-t", "2", "--chunks", "0",
                                     "--max-reads", "1000000", "--no-emit-uncertain",
                                     "--diagnose-withheld-paths", "-o", invocation["outputs"] + "/",
                                     "-s", name, "--pcr-panel-file", panel["path"], sample["available_path"]]
            comparison.inspect_prefix(Path(sample["available_path"]), sample["expected_prefix"])
            execution = comparison.run_timed(invocation, cpus, dict(os.environ, LC_ALL="C"))
            if execution["status"] != "pass":
                raise ValueError(f"{role} opt-out invocation failed")
            captured = comparison.capture(Path(invocation["outputs"]), invocation, sample,
                                          invocation["command"], scopes)
            if any(row["assembly_status"] != "standard" for row in captured["records"]):
                raise ValueError("Opt-out output includes uncertain record")
            result["jobs"][role] = {"invocation": invocation, "execution": execution,
                                    "output": captured}
        old, new = (result["jobs"][role]["output"] for role in ("old", "new"))
        for key in ("records", "counts", "pcr_results"):
            if old[key] != new[key]:
                raise ValueError(f"Diagnostic-only {key} changed")
        hidden = [row["withheld_path_diagnostics"] for target in new["pcr_results"]
                  for row in target["threshold_diagnostics"] if row.get("withheld_path_diagnostics")]
        if not hidden or any(row["schema_version"] != 1 or "marker_detail_truncated" in row
                             for row in hidden):
            raise ValueError("Diagnostic-only payload format changed")
        for role, receipt in result["binaries"].items():
            comparison.verify(receipt)
        comparison.verify(result["input"])
        comparison.verify(result["panel"])
        result.update(status="pass", standard_records=len(new["records"]),
                      hidden_threshold_payloads=len(hidden))
    finally:
        (output / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "standard_records": result["standard_records"],
                      "hidden_threshold_payloads": result["hidden_threshold_payloads"]}))


if __name__ == "__main__":
    main()

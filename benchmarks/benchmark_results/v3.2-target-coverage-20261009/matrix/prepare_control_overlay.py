#!/usr/bin/env python3
"""Preserve newly measured read control evidence in the full target table."""

import json
from pathlib import Path

from build_matrix import HERE, identity, sha


def main():
    source = HERE.parent / "reads/results-v1/summary.json"
    summary = json.loads(source.read_text())
    inventory_path = HERE / "sequence_inventory.json"
    inventory = {identity(r): r for r in json.loads(inventory_path.read_text())}
    controls = [c for sample in summary["samples"] for c in sample["positive_controls"]]
    if len(controls) != 8 or len({identity(c) for c in controls}) != 8:
        raise ValueError("Expected eight distinct controls")
    records = []
    for control in controls:
        candidate = inventory[identity(control)]
        if not control["positive_control_pass"] or candidate["assembly_status"] != "standard" or candidate["length"] != control["length"]:
            raise ValueError("Control identity or occurrence check failed")
        records.append({"accession": control["accession"], "gene": control["gene"], "sha256": control["sha256"],
                        "read_evidence_status": "positive_control_local_exact_q20_occurrences_observed",
                        "read_evidence_links": [str(source), str(HERE.parent / "reads/run-v1" / control["accession"] / "raw.json")],
                        "read_evidence_details": control})
    result = {"sequence_evidence": records,
              "source_hashes": {str(p): sha(p.read_bytes()) for p in (source, inventory_path, Path(__file__))},
              "scope": "Exact Q20 local tile occurrences in eight controls; no full-product truth inferred"}
    with (HERE / "control_evidence_overlay.json").open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "controls": len(records)}))


if __name__ == "__main__":
    main()

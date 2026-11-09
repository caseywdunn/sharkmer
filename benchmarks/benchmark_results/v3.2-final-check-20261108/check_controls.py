#!/usr/bin/env python3
"""Check saved legacy readers and finite summary failure cases without assembly."""

import argparse
import copy
import json
from pathlib import Path

from prepare_protocol import PERFORMANCE, CORRECTED
from run_check import capture, fingerprint, summarize, write_new, HERE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inventory = json.loads((HERE / "inputs/input_inventory.json").read_text())
    legacy = []
    for key, sample in inventory["samples"].items():
        if "baseline" not in sample["expected"]:
            continue
        receipt = json.loads(Path(sample["source_receipts"]["baseline"]["path"]).read_text())
        invocation = receipt["invocation"]
        capture(Path(invocation["outputs"]), sample, "baseline", invocation["command"],
                {"versions": {"baseline": "3.1.0"}})
        legacy.append(key)
    schedule = []
    for pair in (1, 2, 3):
        for key in PERFORMANCE:
            for role in ("baseline", "candidate"):
                schedule.append({"sample": key, "role": role, "pair_index": pair,
                                 "purpose": "performance"})
    for key in inventory["samples"]:
        if key not in PERFORMANCE:
            schedule.append({"sample": key, "role": "candidate", "pair_index": None,
                             "purpose": "inventory"})
    runs = []
    for item in schedule:
        metrics = {"wall_seconds": 1.1 if item["role"] == "candidate" else 1.0,
                   "cpu_seconds": 1.0, "user_seconds": 0.9, "system_seconds": 0.1,
                   "peak_rss_bytes": 10000}
        runs.append({"status": "pass", "invocation": item,
                     "outputs": copy.deepcopy(inventory["samples"][item["sample"]]["expected"][item["role"]]),
                     "execution": {"metrics": metrics}})
    protocol = {"samples": inventory["samples"], "performance_samples": PERFORMANCE,
                "corrected_replaces": CORRECTED}
    result = summarize(runs, protocol)
    assert result["retained_identities"] == 86 and result["retained_with_warning"] == 22
    assert all(cell["changes"]["wall_seconds"]["investigation_flag"]
               and not cell["changes"]["peak_rss_bytes"]["investigation_flag"] for cell in result["cells"])
    failures = []
    controls = (
        ("lost sequence", lambda rows: rows[0]["outputs"]["records"][0].update(sha256="changed")),
        ("unstable classes", lambda rows: rows[1]["outputs"]["records"][0].update(assembly_status="uncertain")),
        ("incomplete run", lambda rows: rows.pop()),
    )
    for name, mutate in controls:
        changed = copy.deepcopy(runs)
        mutate(changed)
        try:
            summarize(changed, protocol)
        except ValueError as error:
            failures.append({"case": name, "rejected": True, "error": str(error)})
        else:
            raise AssertionError("Accepted invalid control: " + name)
    write_new(args.output, {"status": "pass", "legacy_capture_replays": legacy,
                           "synthetic_summary_totals": {key: result[key] for key in
                               ("released_identities", "retained_identities", "retained_with_warning",
                                "historical_panel_inventory", "corrected_panel_inventory")},
                           "resource_flag_control": "10% wall flagged, equal RSS unflagged",
                           "negative_controls": failures,
                           "sources": [fingerprint(HERE / name) for name in
                               ("run_check.py", "prepare_protocol.py", "check_controls.py", "inputs/input_inventory.json")]})


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Combine immutable restoration receipts without copying or rescanning FASTQs.

All ten selected prefixes must have exactly one successful receipt. Failed
attempts remain bound as provenance. This checks current regular-file identity;
the benchmark's pre-execution hash checks still verify the complete FASTQ bytes.
"""

import argparse
import json
from pathlib import Path

from restore_inputs import fingerprint, identity, write_json


HERE = Path(__file__).resolve().parent


def merge(inventory_path, attempt_paths, output):
    if output.exists() or output.is_symlink():
        raise ValueError("Merged availability exists; choose a fresh output")
    inventory_binding = fingerprint(inventory_path)
    inventory = json.loads(inventory_path.read_text())
    expected = {row["accession"]: row for row in inventory["samples"]}
    if len(expected) != 10:
        raise ValueError("Expected ten distinct frozen samples")
    bindings = {"inventory": inventory_binding, "merger": fingerprint(__file__),
                "shared_helper": fingerprint(HERE / "restore_inputs.py")}
    chosen, attempts = {}, []

    def bind(label, row):
        if fingerprint(row["path"]) != row:
            raise ValueError(f"Changed source/receipt: {label}")
        bindings[label] = row

    for attempt_index, path in enumerate(attempt_paths):
        binding = fingerprint(path)
        bind(f"attempt_{attempt_index}", binding)
        attempt = json.loads(path.read_text())
        if attempt.get("schema_version") != 1 or attempt.get("status") not in ("pass", "failed"):
            raise ValueError("Unsupported restoration attempt")
        source_before = attempt["sources_before"]
        if source_before != attempt["sources_after"] or source_before["inventory"] != inventory_binding:
            raise ValueError("Attempt source bookends or inventory disagree")
        if Path(source_before["helper"]["path"]) not in (
                HERE / "restore_inputs.py", HERE / "restore_inputs_v2.py"):
            raise ValueError("Unexpected restoration helper")
        bind(f"helper_{attempt_index}", source_before["helper"])
        seen = set()
        for row in attempt["samples"]:
            accession = row["accession"]
            if accession not in expected or accession in seen:
                raise ValueError("Unexpected or repeated sample in attempt")
            seen.add(accession)
            sample = expected[accession]
            if row["expected_prefix"] != sample["expected_prefix"] or row["panel_id"] != sample["panel"]:
                raise ValueError("Attempt altered historical sample identity")
            bind(f"receipt_{attempt_index}_{accession}", row["receipt"])
            receipt = json.loads(Path(row["receipt"]["path"]).read_text())
            for key in ("accession", "panel_id", "status", "available_path", "expected_prefix"):
                if receipt[key] != row[key]:
                    raise ValueError("Receipt and availability summary disagree")
            if receipt["sources_before"] != source_before:
                raise ValueError("Sample receipt has different restoration sources")
            if row["status"] == "failed":
                continue
            if row["status"] != "verified" or accession in chosen:
                raise ValueError("Ambiguous or invalid successful sample receipt")
            if (receipt["sources_after"] != source_before or
                    receipt["observed_prefix"] != sample["expected_prefix"]):
                raise ValueError("Successful sample did not verify exact raw prefix")
            observed_identity = identity(row["available_path"])
            if (observed_identity != receipt["available_file_identity"] or
                    observed_identity["size"] != sample["expected_prefix"]["size_bytes"]):
                raise ValueError("Verified FASTQ file changed after restoration")
            chosen[accession] = row
        attempts.append({"receipt": binding, "status": attempt["status"],
                         "samples": [{"accession": row["accession"], "status": row["status"],
                                      "receipt": row["receipt"]} for row in attempt["samples"]]})
    if set(chosen) != set(expected):
        raise ValueError(f"Missing verified prefixes: {sorted(set(expected) - set(chosen))}")
    for panel in inventory["panels"]:
        bind(f"panel_{panel['name']}", panel["assembly_panel"])
    after = {key: fingerprint(row["path"]) for key, row in bindings.items()}
    if after != bindings:
        raise ValueError("Sources or receipts changed during merge")
    result = {"schema_version": 1, "status": "pass", "complete_inventory_selected": True,
              "sources_before": bindings, "sources_after": after, "source_attempts": attempts,
              "samples": [chosen[row["accession"]] for row in inventory["samples"]],
              "scope": "Receipt merge and current file-identity checks; no new FASTQ scan or biological inference"}
    write_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "inventory.json")
    parser.add_argument("--attempt", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = merge(args.inventory, args.attempt, args.output)
    print(json.dumps({"status": result["status"], "samples": len(result["samples"]),
                      "output": fingerprint(args.output)}, sort_keys=True))

#!/usr/bin/env python3
"""Recheck final joined artifacts and persist source hashes, including links."""

import argparse
import json
from pathlib import Path

from build_matrix import HERE, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    validation_path = args.directory / "validation.json"
    validation = json.loads(validation_path.read_text())
    if validation["status"] != "pass":
        raise ValueError("Joined output was not validated")
    bound = validation["source_files"]
    verified = {}
    for name, expected in bound.items():
        path = Path(name)
        data = path.read_bytes()
        if sha(data) != expected["sha256"] or len(data) != expected["size_bytes"]:
            raise ValueError("Bound source changed after joining: " + name)
        verified[name] = {**expected, "binding": "join_input_receipt"}
    for name, expected in validation["outputs"].items():
        data = (args.directory / name).read_bytes()
        if sha(data) != expected["sha256"] or len(data) != expected["size_bytes"]:
            raise ValueError("Final joined output changed: " + name)
    producer = HERE / "join_evidence_v2.py"
    if sha(producer.read_bytes()) != validation["producer_sha256"]:
        raise ValueError("Final join producer changed")
    base_validation = json.loads((HERE / "validation.json").read_text())
    base_producer = HERE / "build_matrix.py"
    if sha(base_producer.read_bytes()) != base_validation["producer_sha256"]:
        raise ValueError("Frozen base matrix producer changed")
    # Supplemental control overlay explicitly binds its underlying summary.
    control = json.loads((HERE / "control_evidence_overlay.json").read_text())
    for name, expected in control["source_hashes"].items():
        data = Path(name).read_bytes()
        if sha(data) != expected:
            raise ValueError("Control overlay source changed")
        verified[name] = {"sha256": expected, "size_bytes": len(data), "binding": "control_overlay_receipt"}
    inventory = json.loads((args.directory / "sequence_inventory.json").read_text())
    linked = set()
    for row in inventory:
        linked.update(row["read_evidence_links"])
        linked.update(row["reference_evidence_links"])
    for name in sorted(linked):
        path = Path(name)
        if not path.is_file():
            raise ValueError("Missing evidence link: " + name)
        if name not in verified:
            data = path.read_bytes()
            verified[name] = {"sha256": sha(data), "size_bytes": len(data), "binding": "linked_artifact_snapshot"}
    for path in (validation_path, producer, base_producer, Path(__file__)):
        data = path.read_bytes()
        verified[str(path)] = {"sha256": sha(data), "size_bytes": len(data), "binding": "final_join_verification"}
    result = {"status": "pass", "sources": verified, "source_count": len(verified),
              "source_binding_scope": "Receipt-bound inputs rechecked; additional evidence links are explicitly labeled artifact snapshots.",
              "joined_output_count": len(validation["outputs"]), "emitted_identities": len(inventory)}
    with (args.directory / "joined_source_hashes.json").open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "source_count": len(verified), "joined_outputs": len(validation["outputs"])}))


if __name__ == "__main__":
    main()

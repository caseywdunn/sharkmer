#!/usr/bin/env python3
"""Summarize the frozen local occurrence pilot without assigning haplotypes."""

import hashlib
import json
from pathlib import Path


BASE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    run = BASE / "run-v1"
    receipt = json.loads((run / "execution.json").read_text())
    if receipt["status"] != "pass" or not receipt["sources_unchanged"]:
        raise ValueError("Pilot execution or source attestation failed")
    for name, expected in receipt["outputs"].items():
        path = run / name
        if digest(path) != expected["sha256"] or path.stat().st_size != expected["bytes"]:
            raise ValueError(f"Changed execution output: {name}")
    for name, expected in receipt["source_hashes_after"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Changed frozen source: {name}")
    raw = json.loads((run / "raw.json").read_text())
    derivation = json.loads((BASE / "derivation.json").read_text())
    expected_prefix = json.loads((BASE / "expected_prefix.json").read_text())
    if raw["input_subset"] != expected_prefix or raw["minimum_quality"] != 20 or raw["max_substitutions"] != 0:
        raise ValueError("Input or assay differs from the pilot specification")
    complement = str.maketrans("ACGT", "TGCA")
    rows, its_ids = [], set()
    for event in raw["events"]:
        reads = [read for read in raw["read_evidence"]
                 if any(obs["event_id"] == event["id"] for obs in read["observations"])]
        if len(reads) != event["records"]:
            raise ValueError("Event record count differs from retained ledger")
        ids = {read["header"].split()[0] for read in reads}
        is_its = event["id"].startswith("its_local_")
        if is_its:
            its_ids.update(ids)
        sequences = {min(read["sequence"], read["sequence"].translate(complement)[::-1])
                     for read in reads}
        rows.append({"event_id": event["id"], "role": "ITS_local_occurrence" if is_its else "CO1_control",
                     "window_sequence": event["window_sequence"],
                     "window_length": len(event["window_sequence"]),
                     "exact_q20_records": len(reads), "distinct_read_ids": len(ids),
                     "distinct_full_read_sequences_rc_collapsed": len(sequences),
                     "compatible_frozen_candidates": len({placement["candidate_id"]
                                                          for placement in event["template_placements"]})})
    controls = [row["exact_q20_records"] for row in rows if row["role"] == "CO1_control"]
    if controls != [27, 0]:
        raise ValueError(f"Archived CO1 controls did not reproduce: {controls}")
    its = [row for row in rows if row["role"] == "ITS_local_occurrence"]
    result = {
        "status": "completed_local_occurrence_pilot",
        "verified_input_subset": raw["input_subset"],
        "metrics": json.loads((run / "time.json").read_text()),
        "its_windows_assayed": len(its),
        "its_windows_observed": sum(row["exact_q20_records"] > 0 for row in its),
        "distinct_read_ids_across_its_windows": len(its_ids),
        "control_counts_reproduced": controls,
        "events": rows,
        "limitations": derivation["limitations"] + [
            "High-quality exact occurrence does not exclude sequencing errors or establish biological polymorphism.",
            "Distinct read IDs or sequences do not certify independent molecules. Counts across overlapping events cannot be summed as independent observations.",
        ],
        "next_question": "Which combinations of observed local strings are linked within individual reads, and which full hypotheses remain indistinguishable or unlinked?",
        "source_receipts": {name: digest(BASE / name) for name in
                            ("derivation.json", "run-v1/execution.json", "run-v1/raw.json", "summarize_pilot.py")},
    }
    with (BASE / "summary.json").open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("its_windows_assayed", "its_windows_observed", "distinct_read_ids_across_its_windows", "control_counts_reproduced")}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify saved evidence and summarize every assayed Drosophila hypothesis."""

import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

from scan_prefiltered import ExactPrefilter, HERE, load_scanner


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def span_union(intervals):
    total, end = 0, -1
    for left, right in sorted(intervals):
        total += max(0, right - max(left, end))
        end = max(end, right)
    return total


def main():
    paths = [HERE / name for name in ("derivation.json", "manifest.json", "expected_prefix.json", "prefilter_controls.json",
             "run-v1/raw.json", "run-v1/execution.json", "run-v1/time.json")]
    derivation, manifest, prefix, controls, raw, execution, timing = [json.loads(p.read_text()) for p in paths]
    if execution["status"] != "pass" or not execution["sources_unchanged"] or raw["input_subset"] != prefix:
        raise ValueError("Execution or consumed-prefix validation failed")
    for path, checksum in execution["source_hashes_before"].items():
        if digest(Path(path)) != checksum or execution["source_hashes_after"][path] != checksum:
            raise ValueError("Source changed after execution")
    for name, receipt in execution["outputs"].items():
        path = HERE / "run-v1" / name
        if digest(path) != receipt["sha256"] or path.stat().st_size != receipt["bytes"]:
            raise ValueError("Execution output differs")
    scanner = load_scanner()
    prepared = scanner.prepare_events(manifest, 0)
    events = {e["id"]: e for e in raw["events"]}
    expected = {e["id"]: e for e in prepared}
    if set(events) != set(expected):
        raise ValueError("Event inventory differs")
    for identifier, event in events.items():
        if any(event[k] != expected[identifier][k] for k in ("window_sequence", "window_start", "window_end", "template_placements")):
            raise ValueError("Prepared event differs from saved evidence")
    # Verify every stored qualifying read against the brute-force scanner.
    recounted = Counter()
    for read in raw["read_evidence"]:
        observed = scanner.observations_for_read(read["sequence"].encode(), read["quality"].encode(), prepared, 0)
        if observed != read["observations"]:
            raise ValueError("Saved read observations differ from direct scanner")
        recounted.update(o["event_id"] for o in observed)
    if any(recounted[e["id"]] != e["records"] for e in raw["events"]):
        raise ValueError("Event count differs from saved read ledger")
    # Verify the supplementary quality-unfiltered exact occurrence counter
    # independently of its quality-qualified delegate on synthetic reads.
    adapter, reference_raw = ExactPrefilter(scanner), Counter()
    synthetic = [b"N" * 101]
    for event in prepared:
        for _, pattern in event["patterns"]:
            synthetic.extend([pattern, b"NN" + pattern + b"NN", pattern[:19] + b"N" + pattern[20:]])
    for sequence in synthetic:
        for event in prepared:
            reference_raw[event["id"]] += any(sequence.find(pattern) >= 0 for _, pattern in event["patterns"])
        if adapter(sequence, b"!" * len(sequence), prepared, 0):
            raise ValueError("Q0 synthetic sequence produced primary Q20 evidence")
    if any(reference_raw[e["id"]] != adapter.raw_counts[e["id"]] for e in prepared):
        raise ValueError("Supplementary raw count differs from independent direct count")
    for identifier, values in derivation["expected_co1_controls"].items():
        if any(events[identifier][key] != value for key, value in values.items()):
            raise ValueError("Unchanged CO1 control count differs")
    current_ids = {c["candidate_id"] for c in derivation["candidates"]}
    candidate_sequences = {c["id"]: c["sequence"] for c in manifest["candidates"]}
    rows = []
    historical_ids = set(derivation["historical_12s_aliases"].values())
    for candidate in sorted(derivation["candidates"], key=lambda r: r["candidate_id"]):
        identifier = candidate["candidate_id"]
        if digest_sequence(candidate_sequences[identifier]) != candidate["sha256"]:
            raise ValueError("Candidate identity differs")
        memberships = [(window, member) for window in derivation["windows"] for member in window["memberships"]
                       if member["candidate_id"] == identifier]
        tiles = [(w, m) for w, m in memberships if m["role"] == "backbone_tile"]
        starts = sorted(m["window_start"] for _, m in tiles)
        expected_starts = sorted(set(range(0, candidate["length"] - 60, 21)) | {candidate["length"] - 61})
        if starts != expected_starts:
            raise ValueError("All-candidate backbone schedule incomplete")
        for window, membership in memberships:
            if candidate_sequences[identifier][membership["window_start"]:membership["window_end"]] != window["sequence"]:
                raise ValueError("Collapsed membership has wrong sequence or coordinates")
        nearest = {w["event_id"] for w, m in memberships if m["role"] == "nearest_local_difference"}
        observed_tiles = [(w, m) for w, m in tiles if events[w["event_id"]]["records"]]
        row = {"candidate_id": identifier, "length": candidate["length"], "sha256": candidate["sha256"],
               "gained": candidate["gained"], "historical_loss": identifier in historical_ids,
               "backbone_tiles": len(tiles), "q20_observed_backbone_tiles": len(observed_tiles),
               "unfiltered_observed_backbone_tiles": sum(events[w["event_id"]]["unfiltered_exact_records"] > 0 for w, _ in tiles),
               "q20_observed_span_union": span_union((m["window_start"], m["window_end"]) for _, m in observed_tiles),
               "distinct_nearest_window_strings": len(nearest),
               "q20_observed_nearest_window_strings": sum(events[e]["records"] > 0 for e in nearest),
               "unfiltered_observed_nearest_window_strings": sum(events[e]["unfiltered_exact_records"] > 0 for e in nearest),
               "whole_product_status": "unresolved_no_sample_local_corroboration_in_this_assay"}
        rows.append(row)
    if len(rows) != 26 or sum(r["gained"] for r in rows) != 23 or len(historical_ids) != 3:
        raise ValueError("Scope mismatch")
    primary_structural = [e for e in raw["events"] if e["role"] == "primary_frozen_repeat_interval"]
    event_rows = [{"event_id": e["id"], "candidate_representative": e["candidate_id"],
                   "window_length": len(e["window_sequence"]), "callability": e["assay_callability"],
                   "q20_exact_records": e["records"], "unfiltered_exact_records": e["unfiltered_exact_records"],
                   "finite_context_placements": len(e["template_placements"]),
                   "role": e["role"]} for e in raw["events"]]
    result = {"status": "completed_bounded_local_occurrence_check", "current_12s_candidates": len(rows),
              "gained_12s_candidates": 23, "historical_recovered_candidates": [r for r in rows if r["historical_loss"]],
              "per_candidate": rows, "unique_current_12s_local_strings": len(derivation["windows"]),
              "all_candidate_tile_memberships": sum(r["backbone_tiles"] for r in rows),
              "nearest_endpoint_memberships": sum(m["role"] == "nearest_local_difference" for w in derivation["windows"] for m in w["memberships"]),
              "q20_current_12s_strings_observed": sum(events[w["event_id"]]["records"] > 0 for w in derivation["windows"]),
              "unfiltered_current_12s_strings_observed": sum(events[w["event_id"]]["unfiltered_exact_records"] > 0 for w in derivation["windows"]),
              "historical_structural_events": [{"id": e["id"], "window_length": len(e["window_sequence"]),
                                                "callability": e["assay_callability"], "records": e["records"]} for e in primary_structural],
              "positive_control_events_reproduced": len(derivation["expected_co1_controls"]),
              "positive_control_q20_count_range": [min(events[e]["records"] for e in derivation["expected_co1_controls"]),
                                                    max(events[e]["records"] for e in derivation["expected_co1_controls"])],
              "verified_input_subset": prefix, "metrics": timing, "prefilter": raw["prefilter"],
              "verification": {"q20_ledger_reads_bruteforce_rechecked": len(raw["read_evidence"]),
                               "prefilter_observation_equivalence_cases": controls["exact_observation_equality_cases"],
                               "supplementary_raw_counter_equivalence_cases": len(synthetic),
                               "all_candidate_tile_schedules_checked": 26,
                               "all_control_counts_and_footprints_reproduced": True},
              "interpretation": "None of the 26 current 12S hypotheses receives local exact-string corroboration from this bounded first-million R1 assay. Removing the quality threshold does not explain the 12S zeros. All 16 retained CO1 controls reproduce. Whole-product structure, source locus and biological validity remain unresolved; do not infer absence or reject hypotheses solely from this assay.",
              "limitations": derivation["limitations"] + [
                  "The supplementary unfiltered CO1 counts include some low-quality matches; those supplemental-only reads were not retained. The primary Q20 matching read ledger is complete.",
                  "Seed-only matching records are a performance counter across 12S and CO1 contexts, not a validated 12S k-mer-support measure.",
                  "This assay excludes mismatches and indels and does not reconstruct full repeats or establish why these graph paths were produced."],
              "source_receipts": {str(p.relative_to(HERE)): digest(p) for p in paths + [Path(__file__)]}}
    outputs = {"summary.json": json.dumps(result, indent=2, sort_keys=True) + "\n"}
    for name, records in (("per_candidate.tsv", rows), ("per_event.tsv", event_rows)):
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
        outputs[name] = buffer.getvalue()
    for name, contents in outputs.items():
        with (HERE / name).open("x") as stream:
            stream.write(contents)
    print(json.dumps({"candidates": len(rows), "local_strings": len(derivation["windows"]),
                      "q20_positive_strings": result["q20_current_12s_strings_observed"],
                      "unfiltered_positive_strings": result["unfiltered_current_12s_strings_observed"],
                      "positive_controls_reproduced": result["positive_control_events_reproduced"],
                      "verification": result["verification"]}))


def digest_sequence(value):
    return hashlib.sha256(value.encode()).hexdigest()


if __name__ == "__main__":
    main()

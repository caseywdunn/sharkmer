#!/usr/bin/env python3
"""Join complete local counts to candidate schedules without asserting phase."""

import argparse
from collections import defaultdict
import csv
import io
import json
from pathlib import Path

from scan_exact_v2 import sha, rc


HERE = Path(__file__).resolve().parent


def fingerprint(path):
    return {"sha256": sha(path.read_bytes()), "bytes": path.stat().st_size}


def tsv(rows):
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v for k, v in row.items()} for row in rows)
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifests = sorted((HERE / "manifests-v2").glob("*.json"))
    if len(manifests) != 4:
        raise ValueError("Expected four manifests")
    per_candidate, per_event, per_contrast, sample_reports, evidence_rows = [], [], [], [], []
    sources = {str(Path(__file__)): fingerprint(Path(__file__))}
    for manifest_path in manifests:
        manifest = json.loads(manifest_path.read_text())
        accession = manifest["accession"]
        run_dir = args.runs / accession
        receipt_path = run_dir / "execution.json"
        raw_path = run_dir / "raw.json"
        time_path = run_dir / "time.json"
        receipt = json.loads(receipt_path.read_text())
        raw = json.loads(raw_path.read_text())
        if receipt["status"] != "pass" or raw["verified_prefix"] != manifest["expected_prefix"]:
            raise ValueError("Run or prefix verification failed")
        for p, expected in receipt["source_files_before"].items():
            if fingerprint(Path(p)) != expected:
                raise ValueError("Execution source changed")
        for name, expected in receipt["outputs"].items():
            if fingerprint(run_dir / name) != expected:
                raise ValueError("Run output changed")
        for p in (manifest_path, receipt_path, raw_path, time_path):
            sources[str(p.resolve())] = fingerprint(p)
        counts = {e["id"]: e for e in raw["events"]}
        candidates = {c["id"]: c for c in manifest["candidates"]}
        by_candidate = defaultdict(list)
        by_pattern = {e["sequence"]: e["id"] for e in manifest["events"]}
        for event in manifest["events"]:
            count = counts[event["id"]]
            for membership in event["memberships"]:
                by_candidate[membership["candidate_id"]].append((event, membership))
            per_event.append({"accession": accession, "event_id": event["id"], "length": len(event["sequence"]),
                              "sequence_sha256": event["sha256"], "raw_records": count["raw_records"],
                              "q20_records": count["q20_records"], "raw_occurrences": count["raw_occurrences"],
                              "q20_occurrences": count["q20_occurrences"], "callability": count["callability"],
                              "eligible_length_records": count["eligible_length_records"],
                              "context_candidate_count": len({x["candidate_id"] for x in count["candidate_placements"]}),
                              "retained_q20_records": count["retained_q20_records"],
                              "retained_raw_only_records": count["retained_raw_only_records"],
                              "example_ledger_complete_for_event": count["example_ledger_complete_for_event"]})
        control_reports = []
        for identifier in manifest["focused_candidate_ids"] + manifest["positive_control_ids"]:
            candidate = candidates[identifier]
            memberships = by_candidate[identifier]
            tile_members = [(event, m) for event, m in memberships if m["role"].endswith("_tile")]
            tile_ids = {event["id"] for event, _ in tile_members}
            contrast_ids = {event["id"] for event, m in memberships if m["role"] == "nearest_difference_endpoint"}
            event_ids = tile_ids | contrast_ids
            covered_q20, covered_raw = set(), set()
            for event, membership in tile_members:
                if counts[event["id"]]["q20_records"]:
                    covered_q20.update(range(membership["start"], membership["end"]))
                if counts[event["id"]]["raw_records"]:
                    covered_raw.update(range(membership["start"], membership["end"]))
            observed = {identifier: counts[identifier] for identifier in event_ids if counts[identifier]["q20_records"]}
            row = {"accession": accession, "gene": candidate["gene"], "sha256": candidate["sha256"],
                   "length": len(candidate["sequence"]), "role": candidate["role"],
                   "scheduled_tiles": len(tile_members), "distinct_tile_strings": len(tile_ids),
                   "q20_observed_tile_strings": sum(counts[i]["q20_records"] > 0 for i in tile_ids),
                   "raw_observed_tile_strings": sum(counts[i]["raw_records"] > 0 for i in tile_ids),
                   "q20_observed_tile_span_union_bases": len(covered_q20),
                   "raw_observed_tile_span_union_bases": len(covered_raw),
                   "distinct_nearest_endpoint_strings": len(contrast_ids),
                   "q20_observed_nearest_endpoint_strings": sum(counts[i]["q20_records"] > 0 for i in contrast_ids),
                   "raw_observed_nearest_endpoint_strings": sum(counts[i]["raw_records"] > 0 for i in contrast_ids),
                   "overlength_nearest_endpoint_strings": sum(counts[i]["eligible_length_records"] == 0 for i in contrast_ids),
                   "incomplete_flank_endpoint_memberships": sum(x["candidate_id"] == identifier for x in manifest["uncallable_endpoints"]),
                   "q20_observed_strings_unique_to_candidate_in_finite_context": sum(
                       {x["candidate_id"] for x in count["candidate_placements"]} == {identifier} for count in observed.values()),
                   "local_read_status": "local_exact_q20_occurrences_observed" if observed else
                       "raw_only_local_exact_occurrences" if any(counts[i]["raw_records"] for i in event_ids) else "no_tested_local_exact_occurrence",
                   "whole_product_validation": "not_established"}
            if candidate["role"] == "positive_control":
                row["positive_control_pass"] = row["q20_observed_tile_strings"] > 0
                control_reports.append(row)
            else:
                per_candidate.append(row)
                evidence_rows.append({"accession": accession, "gene": candidate["gene"], "sha256": candidate["sha256"],
                                      "read_evidence_status": row["local_read_status"],
                                      "read_evidence_links": [str(raw_path.resolve()), str((args.output / "per_candidate.tsv").resolve())],
                                      "read_evidence_details": row})
        for contrast in manifest["nearest_contrasts"]:
            row = {"accession": accession, **contrast}
            for side, identifier, end_key in (("candidate", contrast["candidate_id"], "candidate_end"),
                                              ("partner", contrast["partner_id"], "partner_end")):
                sequence = candidates[identifier]["sequence"]
                start, end = contrast["start"] - 21, contrast[end_key] + 21
                if start < 0 or end > len(sequence):
                    row[side + "_callability"] = "missing_complete_flanks"
                    row[side + "_q20_records"] = None
                    row[side + "_raw_records"] = None
                else:
                    piece = sequence[start:end].encode()
                    event_id = by_pattern[min(piece, rc(piece)).decode()]
                    event = counts[event_id]
                    row[side + "_callability"] = event["callability"]
                    row[side + "_q20_records"] = event["q20_records"] if event["eligible_length_records"] else None
                    row[side + "_raw_records"] = event["raw_records"] if event["eligible_length_records"] else None
            per_contrast.append(row)
        if not all(x["positive_control_pass"] for x in control_reports):
            raise ValueError("A predefined positive control lacks Q20 observations")
        sample_reports.append({"accession": accession, "events": len(manifest["events"]),
                               "verified_prefix": raw["verified_prefix"], "read_length_histogram": raw["read_length_histogram"],
                               "positive_controls": control_reports, "metrics": json.loads(time_path.read_text()),
                               "ledger": raw["ledger"], "matching_records": raw["records_with_any_exact_pattern"]})
    targets = []
    for accession, gene in sorted({(r["accession"], r["gene"]) for r in per_candidate}):
        rows = [r for r in per_candidate if (r["accession"], r["gene"]) == (accession, gene)]
        targets.append({"accession": accession, "gene": gene, "candidates": len(rows),
                        "candidates_with_any_q20_local_string": sum(r["local_read_status"] == "local_exact_q20_occurrences_observed" for r in rows),
                        "candidates_with_any_q20_tile": sum(r["q20_observed_tile_strings"] > 0 for r in rows),
                        "candidates_with_all_tiles_q20_observed": sum(r["q20_observed_tile_strings"] == r["distinct_tile_strings"] for r in rows),
                        "candidates_with_q20_nearest_endpoint": sum(r["q20_observed_nearest_endpoint_strings"] > 0 for r in rows),
                        "candidates_without_any_raw_local_string": sum(r["local_read_status"] == "no_tested_local_exact_occurrence" for r in rows),
                        "whole_product_validation": "not_established"})
    if len(per_candidate) != 154 or len(targets) != 6:
        raise ValueError("Focused scope incomplete")
    result = {"status": "complete", "scope": "All 154 candidates in the six previously unrecovered targets",
              "targets": targets, "samples": sample_reports, "source_files": sources,
              "limitations": ["Exact local evidence does not establish locus specificity, biological truth, full-product linkage or independent molecules.",
                              "Tiles overlap and may be shared between candidate contexts. Covered span unions are marginal coverage only.",
                              "Raw/Q20 counts cover every expected input record; example ledgers are explicitly capped and selected.",
                              "Incomplete-flank and overlength windows remain uncallable; exact zero counts do not establish sample absence."]}
    payloads = {"summary.json": json.dumps(result, indent=2, sort_keys=True) + "\n",
                "per_candidate.tsv": tsv(per_candidate), "per_candidate.json": json.dumps(per_candidate, indent=2, sort_keys=True) + "\n",
                "per_event.tsv": tsv(per_event), "contrasts.tsv": tsv(per_contrast), "target_summary.tsv": tsv(targets),
                "evidence_overlay.json": json.dumps({"sequence_evidence": evidence_rows}, indent=2, sort_keys=True) + "\n"}
    args.output.mkdir(parents=True, exist_ok=False)
    for name, content in payloads.items():
        (args.output / name).write_text(content)
    with (args.output / "checksums.json").open("x") as stream:
        json.dump({name: fingerprint(args.output / name) for name in payloads}, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(targets))


if __name__ == "__main__":
    main()

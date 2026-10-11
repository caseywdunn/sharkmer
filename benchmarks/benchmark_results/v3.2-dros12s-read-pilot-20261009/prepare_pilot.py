#!/usr/bin/env python3
"""Freeze a bounded all-candidate Drosophila 12S occurrence assay."""

import csv
import hashlib
import importlib.util
import io
import json
import tarfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RETENTION = HERE.parent / "v3.2-sequence-first-retention-20261001"
HISTORICAL = HERE.parent / "v3.2-read-backed-audit-20260913"
TRIAGE = HERE.parent / "v3.2-gain-triage-20261009"
MEMBER = "scratch/full-v3/runs/007_insecta_SRR31887760_1000000_pair1_candidate/outputs/insecta_SRR31887760_1000000_insecta_12S.fasta"


def digest(value):
    return hashlib.sha256(value).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def write_new(name, value):
    with (HERE / name).open("x") as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True) + "\n")


def main():
    sources, members = {}, {}

    def read(path):
        value = path.read_bytes()
        sources[str(path)] = {"sha256": digest(value), "bytes": len(value)}
        return value

    def archive(directory, selected):
        value = read(directory / "evidence.tar.gz")
        checksum = read(directory / "SHA256SUMS").decode().split()[0]
        if digest(value) != checksum:
            raise ValueError("Archive checksum differs")
        result = {}
        with tarfile.open(fileobj=io.BytesIO(value), mode="r:gz") as contents:
            for name in selected:
                result[name] = contents.extractfile(name).read()
                members[f"{directory.name}/{name}"] = {
                    "sha256": digest(result[name]), "bytes": len(result[name])}
        return result

    archived = archive(RETENTION, [MEMBER])
    old_names = ["scan_specs/SRR31887760.scanner_manifest.json",
                 "scan_specs/SRR31887760.expected_prefix.json",
                 "execution_replication_outputs/SRR31887760.exact.json"]
    old = archive(HISTORICAL, old_names)
    old_manifest, prefix, old_result = [json.loads(old[n]) for n in old_names]
    triage_path = TRIAGE / "triage_saved_gains.py"
    read(triage_path)
    triage = module("triage", triage_path)
    fasta = triage.parse_fasta(archived[MEMBER].decode())
    mappings = json.loads(read(RETENTION / "blast-gains-v3/per_sequence.json"))
    gained = {r["sha256"]: r["query_id"] for r in mappings
              if r["accession"] == "SRR31887760" and r["gene"] == "insecta_12S"}
    restoration = json.loads(read(HERE.parent / "v3.2-uncertain-products-20260927/inputs/restoration/SRR31887760.first1M_R1.retry1.fastq.restore.json"))
    if prefix != restoration["expected"] or prefix != restoration["observed"]:
        raise ValueError("Historical and restored consumed-prefix identities differ")
    candidates, identities = [], []
    for name, sequence in fasta.items():
        sha = digest(sequence.encode())
        identifier = gained.get(sha, "current_" + sha[:16])
        if f"sequence_sha256={sha}" not in archived[MEMBER].decode():
            raise ValueError("FASTA header hash not found")
        candidates.append({"id": identifier, "sequence": sequence,
                           "sequence_sha256": sha})
        identities.append({"candidate_id": identifier, "source_fasta_id": name,
                           "length": len(sequence), "sha256": sha,
                           "gained": sha in gained})
    candidates.sort(key=lambda c: c["id"])
    if len(candidates) != 26 or sum(r["gained"] for r in identities) != 23:
        raise ValueError("Expected all 26 current candidates including 23 gains")
    by_sha = {c["sequence_sha256"]: c for c in candidates}
    aliases = {}
    for candidate in old_manifest["candidates"]:
        sha = candidate["sequence_sha256"]
        if sha in by_sha:
            aliases[candidate["id"]] = by_sha[sha]["id"]
        else:
            candidates.append(candidate)
    old_12s = [c for c in old_manifest["candidates"] if len(c["sequence"]) in (475, 487, 499)]
    if len(old_12s) != 3 or not all(c["sequence_sha256"] in by_sha for c in old_12s):
        raise ValueError("Three historical 12S identities not recovered")
    group = defaultdict(list)
    comparisons = []

    def add(sequence, candidate, start, end, role, **extra):
        group[sequence].append({"candidate_id": candidate["id"], "window_start": start,
                                "window_end": end, "role": role, **extra})

    for candidate in by_sha.values():
        sequence = candidate["sequence"]
        starts = sorted(set(range(0, len(sequence) - 60, 21)) | {len(sequence) - 61})
        for start in starts:
            add(sequence[start:start + 61], candidate, start, start + 61, "backbone_tile")
        choices = []
        for other in by_sha.values():
            if other == candidate:
                continue
            window = triage.differing_window(sequence, other["sequence"])
            choices.append(((window["max_allele_length"], window["total_allele_length"],
                             other["sequence_sha256"]), other, window))
        _, other, window = min(choices, key=lambda x: x[0])
        row = {"query_id": candidate["id"], "neighbor_id": other["id"], **window}
        comparisons.append(row)
        for side, member in (("query", candidate), ("neighbor", other)):
            start, end = window[f"{side}_start"] - 21, window[f"{side}_end"] + 21
            if start < 0 or end > len(member["sequence"]):
                raise ValueError("Nearest difference lacks full 21-base flanks")
            add(member["sequence"][start:end], member, start, end,
                "nearest_local_difference", query_id=candidate["id"],
                neighbor_id=other["id"], endpoint_role=side)
    events, windows = [], []
    for ordinal, (sequence, memberships) in enumerate(sorted(group.items()), 1):
        representative = min(memberships, key=lambda r: (r["candidate_id"], r["window_start"], r["role"]))
        identifier = f"dros_local_{ordinal:03d}"
        events.append({"id": identifier, "candidate_id": representative["candidate_id"],
                       "start": representative["window_start"] + 21,
                       "end": representative["window_end"] - 21,
                       "comparison_group": "current_12S_local_occurrence",
                       "role": "collapsed_local_occurrence"})
        windows.append({"event_id": identifier, "sequence": sequence,
                        "length": len(sequence), "memberships": memberships})
    co1 = next(c for c in old_manifest["candidates"] if len(c["sequence"]) == 352)
    controls = [e for e in old_manifest["events"] if e["candidate_id"] == co1["id"]]
    structural = [e for e in old_manifest["events"]
                  if e["candidate_id"] in aliases and e["role"] == "primary_frozen_repeat_interval"]
    for event in structural:
        events.append({**event, "candidate_id": aliases[event["candidate_id"]]})
    events.extend(controls)
    if len(controls) != 16 or len(structural) != 6 or len(events) > 400:
        raise ValueError("Unexpected control, structural or event count")
    manifest = {"schema_version": 1, "candidates": candidates, "events": events}
    scanner_path = ROOT / "scripts/audit_read_bridges.py"
    read(scanner_path)
    scanner = module("bridge_scanner", scanner_path)
    prepared = scanner.prepare_events(manifest, 0)
    old_prepared = {e["id"]: e for e in scanner.prepare_events(old_manifest, 0)}
    for e in prepared[-16:]:
        old_e = old_prepared[e["id"]]
        if any(e[key] != old_e[key] for key in ("window_sequence", "window_start", "window_end", "template_placements")):
            raise ValueError("CO1 control placement changed")
    for window, event in zip(windows, prepared):
        assert event["window_sequence"] == window["sequence"]
        window["finite_context_template_placements"] = event["template_placements"]
    expected_controls = {e["id"]: {k: e[k] for k in ("records", "unique_records", "ambiguous_records", "distinct_ids", "distinct_footprints")}
                         for e in old_result["events"] if e["candidate_id"] == co1["id"]}
    derivation = {"schema_version": 1, "producer_sha256": digest(Path(__file__).read_bytes()),
                  "source_files": sources, "archive_members": members,
                  "candidates": identities, "historical_12s_aliases": aliases,
                  "nearest_comparisons": comparisons, "windows": windows,
                  "events": len(events), "unique_local_windows": len(windows),
                  "current_12s_candidates": len(by_sha), "gained_12s_candidates": len(gained),
                  "control_events": len(controls), "old_structural_events": len(structural),
                  "window_length_counts": {str(n): sum(e["window_end"] - e["window_start"] == n for e in prepared)
                                           for n in sorted({e["window_end"] - e["window_start"] for e in prepared})},
                  "over_read_length_events": [e["id"] for e in prepared if e["window_end"] - e["window_start"] > 101],
                  "backbone_schedule": "61 bases; starts 0,21,42,... while complete; add final end-aligned window; collapse identical strings",
                  "expected_co1_controls": expected_controls,
                  "limitations": [
                      "Local exact strings are occurrence assays; nearest differences are one literal representation and do not localize repeat indels uniquely.",
                      "Overlapping tiles and observations are not independent; a union of observed spans is not whole-product phase or truth.",
                      "Placements and uniqueness inspect only these 28 frozen candidate contexts in both orientations; no genome-wide target specificity is established.",
                      "A callable window with no Q20 exact match in first-million R1 is uncorroborated, not biologically absent or false.",
                      "Windows longer than 101 bases are uncallable in this fixed input, regardless of zero count.",
                      "Supplementary unfiltered exact counts diagnose sensitivity to base quality and do not replace the Q20 primary assay.",
                      "All current 12S records are uncertain hypotheses; no standard 12S positive control exists in this comparison.",
                      "61-base backbone tiles share strings across full hypotheses; even complete tile occurrence would not establish their combination."]}
    write_new("manifest.json", manifest)
    write_new("expected_prefix.json", prefix)
    write_new("derivation.json", derivation)
    write_new("old_control_manifest.json", {"schema_version": 1, "candidates": old_manifest["candidates"], "events": controls})
    print(json.dumps({k: derivation[k] for k in ("events", "unique_local_windows", "current_12s_candidates", "window_length_counts", "over_read_length_events")}))


if __name__ == "__main__":
    main()

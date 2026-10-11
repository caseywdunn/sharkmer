#!/usr/bin/env python3
"""Freeze exact local windows for the saved Gryllus ITS gain pilot."""

import csv
import hashlib
import importlib.util
import io
import json
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TRIAGE = HERE.parent / "v3.2-gain-triage-20261009"
PRIOR = HERE.parent / "v3.2-its-recruitment-20260925/inputs/originals/inputs"


def digest(value):
    return hashlib.sha256(value).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    sources = {}

    def read(path, expected=None):
        path = Path(path)
        value = path.read_bytes()
        observed = {"path": str(path.resolve()), "bytes": len(value), "sha256": digest(value)}
        if expected and any(observed[key] != expected[key] for key in ("bytes", "sha256")):
            raise ValueError(f"Source identity differs: {path}")
        sources[str(path.resolve())] = observed
        return value

    report = json.loads(read(TRIAGE / "report.json"))
    input_bytes = {key: read(value["path"], value) for key, value in report["source_files"].items()}
    window_spec = report["result_files"]["nearest_windows.tsv"]
    rows = list(csv.DictReader(io.StringIO(read(window_spec["path"], window_spec).decode()), delimiter="\t"))
    triage_path = TRIAGE / "triage_saved_gains.py"
    read(triage_path)
    if sources[str(triage_path.resolve())]["sha256"] != report["producer_sha256"]:
        raise ValueError("Triage producer identity differs")
    triage = load_module("saved_gain_triage", triage_path)
    records = triage.bind_inputs(
        list(csv.DictReader(io.StringIO(input_bytes["gains"].decode()), delimiter="\t")),
        json.loads(input_bytes["mapping"]), triage.parse_fasta(input_bytes["fasta"].decode()))
    selected = {r["query_id"]: r for r in records
                if r["accession"] == "SRR27962769" and r["gene"] == "insecta_ITS_2"}
    rows = [r for r in rows if r["accession"] == "SRR27962769" and r["gene"] == "insecta_ITS_2"]
    if len(selected) != 32 or len(rows) != 32 or {r["query_id"] for r in rows} != set(selected):
        raise ValueError("Expected exactly 32 selected gain queries and nearest windows")
    control_bytes = read(PRIOR / "co1_control_manifest.json")
    controls = json.loads(control_bytes)
    if len(controls["candidates"]) != 5 or len(controls["events"]) != 2:
        raise ValueError("Unexpected frozen control manifest")
    prefix_bytes = read(PRIOR / "expected_prefix.json")
    expected_prefix = json.loads(prefix_bytes)
    if expected_prefix != {"records": 1000000, "bases": 150000000, "size_bytes": 365356898,
                           "sha256": "6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794"}:
        raise ValueError("Unexpected consumed-prefix identity")

    candidates = [{"id": q, "sequence": r["sequence"], "sequence_sha256": r["sha256"],
                   "role": "saved_gained_uncertain_hypothesis"} for q, r in sorted(selected.items())]
    candidates.extend(controls["candidates"])
    grouped = defaultdict(list)
    for row in rows:
        first, second = selected[row["query_id"]], selected[row["neighbor_query_id"]]
        for key, rec in (("query", first), ("neighbor", second)):
            if row[f"{key}_sha256"] != rec["sha256"] or int(row[f"{key}_length"]) != rec["length"]:
                raise ValueError("Window identity differs from bound query")
        reconstructed = triage.differing_window(first["sequence"], second["sequence"])
        if any(str(value) != row[key] for key, value in reconstructed.items()):
            raise ValueError("Nearest-window coordinates or alleles differ")
        for role, rec in (("query", first), ("neighbor", second)):
            start, end = int(row[f"{role}_start"]), int(row[f"{role}_end"])
            if start < 21 or end + 21 > rec["length"]:
                raise ValueError("Missing complete 21-base flanks")
            sequence = rec["sequence"][start - 21:end + 21]
            if len(sequence) > 150:
                raise ValueError("Local window exceeds fixed read length")
            grouped[sequence].append({
                "candidate_id": rec["query_id"], "candidate_sha256": rec["sha256"],
                "start": start, "end": end, "allele": row[f"{role}_allele"],
                "window_start": start - 21, "window_end": end + 21,
                "nearest_query_id": row["query_id"], "nearest_neighbor_id": row["neighbor_query_id"],
                "endpoint_role": role,
            })
    events, window_derivation = [], []
    for ordinal, sequence in enumerate(sorted(grouped), 1):
        members = sorted(grouped[sequence], key=lambda r: (
            r["candidate_id"], r["start"], r["end"], r["nearest_query_id"], r["endpoint_role"]))
        representative = members[0]
        identifier = f"its_local_{ordinal:02d}"
        events.append({"id": identifier, "candidate_id": representative["candidate_id"],
                       "start": representative["start"], "end": representative["end"],
                       "comparison_group": "ITS_gain_local_windows",
                       "role": "local_window_occurrence_only"})
        window_derivation.append({"event_id": identifier, "sequence": sequence,
                                  "sequence_sha256": digest(sequence.encode()),
                                  "length": len(sequence), "memberships": members})
    events.extend(controls["events"])
    if len(events) > 64:
        raise ValueError("Event budget exceeded")
    manifest = {
        "schema_version": 1, "candidates": candidates, "events": events,
        "purpose": "Exact Q20 local window occurrence pilot; no full-hypothesis support claim",
    }
    scanner_path = ROOT / "scripts/audit_read_bridges.py"
    read(scanner_path)
    scanner = load_module("frozen_bridge_scanner", scanner_path)
    prepared = scanner.prepare_events(manifest, max_substitutions=0)
    old_controls = scanner.prepare_events(controls, max_substitutions=0)
    if any(event["callability"] != "assayed" for event in prepared):
        raise ValueError("All pilot events must be callable")
    if manifest["events"][-2:] != controls["events"]:
        raise ValueError("CO1 event definitions changed")
    for old, new in zip(old_controls, prepared[-2:]):
        for key in ("window_sequence", "window_start", "window_end", "template_placements"):
            if old[key] != new[key]:
                raise ValueError("CO1 windows or finite-context ambiguity changed")
    for event, window in zip(prepared, window_derivation):
        if event["window_sequence"] != window["sequence"]:
            raise ValueError("Scanner window differs from derivation")
        window["finite_context_template_placements"] = event["template_placements"]
        window["unique_template_placement"] = event["unique_template_placement"]
    sequence_list = [r["sequence"] for r in selected.values()]
    common_prefix = 0
    while len({s[common_prefix] for s in sequence_list}) == 1:
        common_prefix += 1
    common_suffix = 0
    while len({s[-common_suffix - 1] for s in sequence_list}) == 1:
        common_suffix += 1
    derivation = {
        "schema_version": 1, "producer_sha256": digest(Path(__file__).read_bytes()),
        "source_files": sources, "selected_queries": len(selected), "nearest_window_rows": len(rows),
        "endpoint_memberships": sum(len(x) for x in grouped.values()),
        "unique_local_windows": len(grouped), "control_events": len(controls["events"]),
        "total_events": len(events), "candidate_context_count": len(candidates),
        "coordinate_system": "zero_based_half_open_in_emitted_oriented_sequences",
        "representative_rule": "Identical forward window strings collapse; smallest candidate ID then coordinates represents each window",
        "windows": window_derivation,
        "global_gain_shared_prefix": common_prefix, "global_gain_shared_suffix": common_suffix,
        "global_variable_extent_with_21bp_flanks_lengths": sorted({len(s) - common_prefix - common_suffix + 42 for s in sequence_list}),
        "limitations": [
            "ITS read matches establish only occurrence of the local string; multiple full hypotheses share each string.",
            "Scanner unique-record and corroborated-local-bridge fields must not be interpreted as full-hypothesis support.",
            "Candidate IDs on collapsed events identify coordinate representatives; memberships preserve every original nearest-window endpoint.",
            "Finite-context placements inspect only the 32 gains and five old control candidates; target specificity is not established.",
            "The exact assay requires Q20 or better throughout both 21-base flanks and the intervening allele, in either read orientation.",
            "No matching read in this fixed first-million R1 prefix does not prove absence; mismatches, gaps, sampling and low quality limit sensitivity.",
            "Overlapping windows may share reads and cannot be counted as independent evidence.",
            "Nearest-window indel coordinates are one representation; the literal spanning string is assayed without asserting a unique mutation position.",
            "The entire variable extent of the 32 gains is longer than a 150-base read; local matches cannot establish whole-product phase or biological correctness.",
        ],
    }
    outputs = {"manifest.json": (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode(),
               "expected_prefix.json": prefix_bytes}
    derivation["result_files"] = {name: {"sha256": digest(value), "bytes": len(value)} for name, value in outputs.items()}
    outputs["derivation.json"] = (json.dumps(derivation, indent=2, sort_keys=True) + "\n").encode()
    for name in outputs:
        if (HERE / name).exists() or (HERE / name).is_symlink():
            raise ValueError(f"Refusing to overwrite {name}")
    for name, value in outputs.items():
        with (HERE / name).open("xb") as handle:
            handle.write(value)
    print(json.dumps({"windows": len(grouped), "events": len(events), "candidates": len(candidates),
                      "endpoint_memberships": derivation["endpoint_memberships"],
                      "output": str(HERE)}))


if __name__ == "__main__":
    main()

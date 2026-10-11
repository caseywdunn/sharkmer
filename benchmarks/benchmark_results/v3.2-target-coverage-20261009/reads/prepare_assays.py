#!/usr/bin/env python3
"""Freeze deterministic six-target windows and same-sample standard controls."""

import argparse
from collections import defaultdict
import json
from pathlib import Path

from scan_exact import rc, sha


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FOCUS = {
    "SRR25378184": {"angiospermae_psbA-trnH": 49, "angiospermae_trnV-atpE": 41, "angiospermae_atpB-rbcL": 1},
    "ERR14009273": {"angiospermae_trnL-F": 1},
    "SRR1057608": {"insecta_ND4": 44},
    "SRR27962769": {"insecta_18S_1": 18},
}
CONTROL_GENES = {
    "SRR25378184": ["angiospermae_psbM-trnD", "angiospermae_rpl36-infA-rps8"],
    "ERR14009273": ["angiospermae_ITS", "angiospermae_atpB-rbcL"],
    "SRR1057608": ["insecta_CO1_1", "insecta_CO2_1"],
    "SRR27962769": ["insecta_CO1_1", "insecta_CO2_1"],
}


def candidate_id(row):
    return row["accession"] + "|" + row["gene"] + "|" + row["sha256"]


def tile_starts(length):
    if length < 61:
        return []
    return sorted(set(range(0, length - 60, 21)) | {length - 61})


def difference(a, b):
    left = 0
    while left < min(len(a), len(b)) and a[left] == b[left]:
        left += 1
    suffix = 0
    while suffix < min(len(a), len(b)) - left and a[-1 - suffix] == b[-1 - suffix]:
        suffix += 1
    return left, len(a) - suffix, len(b) - suffix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    availability_path = ROOT / "tmp/inclusive-release-20260929/availability.json"
    source_paths = [args.inventory, availability_path, Path(__file__), HERE / "scan_exact.py"]
    sources = {str(p.resolve()): {"sha256": sha(p.read_bytes()), "bytes": p.stat().st_size} for p in source_paths}
    inventory = json.loads(args.inventory.read_text())
    availability = {x["accession"]: x for x in json.loads(availability_path.read_text())["samples"]}
    manifests = {}
    for accession, genes in FOCUS.items():
        context = [r for r in inventory if r["accession"] == accession]
        focus = [r for r in context if r["gene"] in genes]
        for gene, expected_count in genes.items():
            selected = [r for r in focus if r["gene"] == gene]
            if len(selected) != expected_count or any(not r["new_vs_release"] for r in selected):
                raise ValueError("Focused candidate inventory differs")
        controls = []
        for gene in CONTROL_GENES[accession]:
            selected = sorted((r for r in context if r["gene"] == gene and r["assembly_status"] == "standard"), key=lambda r: r["sha256"])
            if not selected or not selected[0]["historical_retained"]:
                raise ValueError("Missing unchanged standard control")
            controls.append(selected[0])
        patterns = defaultdict(list)
        uncallable = []
        contrasts = []

        def add(row, start, end, role, details=None):
            membership = {"candidate_id": candidate_id(row), "gene": row["gene"],
                          "start": start, "end": end, "role": role, **(details or {})}
            if start < 0 or end > len(row["sequence"]):
                uncallable.append({**membership, "reason": "missing_complete_21base_flanks"})
                return
            pattern = row["sequence"][start:end].encode("ascii")
            canonical = min(pattern, rc(pattern)).decode("ascii")
            membership["pattern_strand_in_candidate"] = "+" if canonical == pattern.decode("ascii") else "-"
            patterns[canonical].append(membership)

        for row in focus + controls:
            role = "candidate_backbone_tile" if row in focus else "positive_control_tile"
            starts = tile_starts(len(row["sequence"]))
            if not starts:
                raise ValueError("A candidate is too short for the preregistered tiles")
            for start in starts:
                add(row, start, start + 61, role)
        for row in focus:
            others = [x for x in focus if x["gene"] == row["gene"] and x["sha256"] != row["sha256"]]
            if not others:
                continue
            def rank(other):
                start, end_a, end_b = difference(row["sequence"], other["sequence"])
                return max(end_a - start, end_b - start), end_a + end_b - 2 * start, other["sha256"]
            other = min(others, key=rank)
            start, end_a, end_b = difference(row["sequence"], other["sequence"])
            if row["sequence"][:start] + other["sequence"][start:end_b] + row["sequence"][end_a:] != other["sequence"]:
                raise ValueError("Nearest difference does not reconstruct its partner")
            contrast_id = candidate_id(row) + "->" + candidate_id(other)
            contrasts.append({"id": contrast_id, "candidate_id": candidate_id(row), "partner_id": candidate_id(other),
                              "start": start, "candidate_end": end_a, "partner_end": end_b,
                              "definition": "maximal exact prefix and nonoverlapping suffix; not an alignment or unique indel placement"})
            add(row, start - 21, end_a + 21, "nearest_difference_endpoint", {"contrast_id": contrast_id, "endpoint": "candidate"})
            add(other, start - 21, end_b + 21, "nearest_difference_endpoint", {"contrast_id": contrast_id, "endpoint": "partner"})
        events = [{"id": "window_" + sha(seq.encode())[:16], "sha256": sha(seq.encode()),
                   "sequence": seq, "memberships": members} for seq, members in sorted(patterns.items())]
        if len({e["id"] for e in events}) != len(events):
            raise ValueError("Event ID collision")
        sample = availability[accession]
        if sample["status"] != "verified" or not Path(sample["available_path"]).is_file():
            raise ValueError("Original input unavailable")
        control_ids = {candidate_id(r) for r in controls}
        focused_ids = {candidate_id(r) for r in focus}
        manifest = {
            "schema_version": 1, "accession": accession, "expected_prefix": sample["expected_prefix"],
            "fastq_path": sample["available_path"], "source_files": sources,
            "focus_targets": genes, "control_policy": "Two predefined genes per sample; lexicographically first full hash among unchanged standard records",
            "positive_control_ids": sorted(control_ids),
            "positive_control_acceptance": "At least one Q20 exact 61-base tile observed for each control; exact counts not preregistered for these new tile schedules",
            "focused_candidate_ids": sorted(focused_ids),
            "candidates": [{"id": candidate_id(r), "gene": r["gene"], "sha256": r["sha256"], "sequence": r["sequence"],
                            "assembly_status": r["assembly_status"], "role": "focused" if candidate_id(r) in focused_ids else
                            "positive_control" if candidate_id(r) in control_ids else "placement_context"} for r in context],
            "events": events, "uncallable_endpoints": uncallable, "nearest_contrasts": contrasts,
            "protocol": {"tile_length": 61, "tile_step": 21, "last_tile_end_aligned": True, "contrast_flank": 21,
                         "minimum_phred": 20, "supplementary_quality_unfiltered_counts": True, "exact_only": True,
                         "both_orientations": True, "all_overlapping_occurrences": True,
                         "maximum_seconds": 600, "maximum_address_space_bytes": 2147483648, "cpus": 1},
            "limits": ["Local strings can be shared among hypotheses and other loci; whole-product linkage and biological validity are not established.",
                       "Candidate context covers only the current outputs for this sample, not an entire genome or reference universe.",
                       "A zero count is bounded non-observation; windows exceeding all read lengths and incomplete flanks are uncallable.",
                       "Shared read examples are capped separately from complete counts; no molecule or footprint independence inferred."]}
        manifests[accession] = manifest
    args.output.mkdir(parents=True, exist_ok=True)
    for accession, manifest in manifests.items():
        with (args.output / (accession + ".json")).open("x") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
        print(json.dumps({"accession": accession, "focus_candidates": len(manifest["focused_candidate_ids"]),
                          "events": len(manifest["events"]), "uncallable_endpoints": len(manifest["uncallable_endpoints"])}))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify archived evidence, restore fresh, and derive unchanged local targets."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile


EXPECTED_ARCHIVE_SHA256 = "14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3"
ORIGINALS = (
    "inputs/candidate.json",
    "inputs/expected_prefix.json",
    "inputs/co1_control_manifest.json",
    "inputs/prior_read_manifest.json",
    "inputs/source_manifest.json",
    "local_assays/region_1.scanner_manifest.json",
    "local_assays/region_2.scanner_manifest.json",
)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def publish_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write("\n")


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def verify_restore(archive_dir, restore):
    """Check both outer checksums and every regular member before any extraction."""
    if restore.exists():
        raise ValueError(f"Restore directory must be fresh: {restore}")
    for line in (archive_dir / "SHA256SUMS").read_text().splitlines():
        expected, filename = line.split(maxsplit=1)
        filename = filename.lstrip("*")
        if filename not in {"ARCHIVE_CONTENTS.json", "evidence.tar.gz"}:
            raise ValueError(f"Unexpected outer checksum filename: {filename}")
        if sha((archive_dir / filename).read_bytes()) != expected:
            raise ValueError(f"Outer checksum mismatch: {filename}")
    archive_path = archive_dir / "evidence.tar.gz"
    if sha(archive_path.read_bytes()) != EXPECTED_ARCHIVE_SHA256:
        raise ValueError("Archive differs from handoff-pinned identity")
    manifest = json.loads((archive_dir / "ARCHIVE_CONTENTS.json").read_text())
    expected = {row["path"]: row for row in manifest["members"]}
    if len(expected) != len(manifest["members"]):
        raise ValueError("Duplicate archive member in manifest")
    verified = {}
    with tarfile.open(archive_path) as archive:
        for member in archive.getmembers():
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts or not path.parts:
                raise ValueError(f"Unsafe archive path: {member.name}")
            if member.isdir():
                continue
            if not member.isfile() or member.name in verified:
                raise ValueError(f"Unexpected nonregular or duplicate member: {member.name}")
            data = archive.extractfile(member).read()
            expected_row = expected.get(member.name)
            if expected_row is None or len(data) != expected_row["size_bytes"] or sha(data) != expected_row["sha256"]:
                raise ValueError(f"Member checksum/size mismatch: {member.name}")
            verified[member.name] = data
    if set(verified) != set(expected):
        raise ValueError("Archive member set differs from manifest")
    restore.mkdir(parents=True)
    for name, data in verified.items():
        target = restore / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(data)
    return verified


def derive(verified, output, archive_dir, restore):
    regions = []
    contexts = []
    co1 = json.loads(verified["inputs/co1_control_manifest.json"])
    full_candidates = {candidate["id"]: candidate for candidate in co1["candidates"]}
    its = json.loads(verified["inputs/candidate.json"])["record"]
    its_id = next(key for key, row in full_candidates.items() if row["sequence_sha256"] == its["sequence_sha256"])
    if full_candidates[its_id]["sequence"] != its["oriented_sequence"]:
        raise ValueError("Historical ITS differs between archived manifests")
    for region_id in ("region_1", "region_2"):
        source_name = f"local_assays/{region_id}.scanner_manifest.json"
        original = json.loads(verified[source_name])
        classes = []
        for candidate in original["candidates"]:
            sequence = candidate["sequence"]
            if sha(sequence.encode()) != candidate["sequence_sha256"]:
                raise ValueError("Archived class sequence hash mismatch")
            event = next(event for event in original["events"] if event["candidate_id"] == candidate["id"])
            if event["start"] != 21 or event["end"] != len(sequence) - 21:
                raise ValueError("Unexpected archived ITS flank span")
            self_context = []
            # Equality is checked in both orientations without assuming genomic truth.
            for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                start = its["oriented_sequence"].find(oriented)
                while start >= 0:
                    self_context.append({"decoy_id": its_id, "strand": strand, "start": start, "end": start + len(sequence)})
                    start = its["oriented_sequence"].find(oriented, start + 1)
            row = {"id": candidate["id"], "sequence": sequence, "sequence_sha256": candidate["sequence_sha256"], "window_start": 0, "candidate_id": candidate["id"], "kind": "historical_local_window" if self_context else "extra_cycle_local_window", "self_context": self_context}
            classes.append(row)
            contexts.append({"class_id": row["id"], "source_member": source_name, "source_candidate_id": candidate["id"], "source_event_id": event["id"], "source_interval": [0, len(sequence)], "self_context": self_context})
        regions.append({"id": region_id, "classes": classes})
    co1_classes = []
    for event in co1["events"]:
        candidate = full_candidates[event["candidate_id"]]
        start, end = event["start"] - 21, event["end"] + 21
        sequence = candidate["sequence"][start:end]
        if not 0 <= start < end <= len(candidate["sequence"]):
            raise ValueError("Invalid CO1 window")
        self_context = [{"decoy_id": candidate["id"], "strand": "+", "start": start, "end": end}]
        co1_classes.append({"id": event["id"], "sequence": sequence, "sequence_sha256": sha(sequence.encode()), "window_start": start, "candidate_id": candidate["id"], "kind": event["role"], "self_context": self_context})
        contexts.append({"class_id": event["id"], "source_member": "inputs/co1_control_manifest.json", "source_candidate_id": candidate["id"], "source_event_id": event["id"], "source_interval": [start, end], "self_context": self_context})
    regions.append({"id": "co1_control", "classes": co1_classes})
    decoys = []
    for candidate in co1["candidates"]:
        if sha(candidate["sequence"].encode()) != candidate["sequence_sha256"]:
            raise ValueError("Archived full-product sequence hash mismatch")
        decoys.append({key: candidate[key] for key in ("id", "sequence", "sequence_sha256")})
    manifest = {"schema_version": 1, "regions": regions, "decoys": decoys}
    publish_json(output / "manifest.json", manifest)
    originals = []
    for name in ORIGINALS:
        dest = output / "originals" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("xb") as stream:
            stream.write(verified[name])
        originals.append({"path": str(dest.relative_to(output)), "archive_member": name, "size_bytes": len(verified[name]), "sha256": sha(verified[name])})
    provenance = {"schema_version": 1, "scope": "Unchanged archived sequence hypotheses and controls; no new read outcomes inspected.", "archive_path": str(archive_dir / "evidence.tar.gz"), "archive_sha256": EXPECTED_ARCHIVE_SHA256, "archive_contents_sha256": sha((archive_dir / "ARCHIVE_CONTENTS.json").read_bytes()), "verified_members": len(verified), "restore_path": str(restore), "script_sha256": sha(Path(__file__).read_bytes()), "derived_manifest_sha256": sha((output / "manifest.json").read_bytes()), "originals": originals, "classes": contexts, "decoy_scope": "All five original CO1-control competitors: historical ITS, 16S, 12S, retained CO1, lost CO1. These are hypotheses, not independent references.", "coordinate_semantics": "Zero-based half-open; self_context strand expresses class orientation relative to decoy.", "unassayed_limitations": "Extra-cycle sequences have no asserted full-product self-context; local alignments do not establish repeat copy number or global phase."}
    publish_json(output / "provenance.json", provenance)
    print(json.dumps({"archive_members_verified": len(verified), "manifest_sha256": provenance["derived_manifest_sha256"], "classes": [{"id": row["id"], "length": len(row["sequence"]), "self_context": row["self_context"]} for region in regions for row in region["classes"]], "decoys": len(decoys)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--restore", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verified = verify_restore(args.archive_dir, args.restore)
    derive(verified, args.output, args.archive_dir, args.restore)


if __name__ == "__main__":
    main()

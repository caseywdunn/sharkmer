#!/usr/bin/env python3
"""Independently verify the archived, restored, and derived diagnostic inputs.

Does not import the derivation helper or inspect new reads/results.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile


ARCHIVE_SHA = "14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3"
CONTENTS_SHA = "dcf9821b93b1b993115890953b6c45a3ab82e7dda090e988100f2523a9deda10"
ITS_ID = "SRR27962769_23163ca76f97f5d1"
EXPECTED_DECOYS = {
    ITS_ID: ("ITS_2", 978),
    "SRR27962769_66686929b8d17dfa": ("16S_2", 534),
    "SRR27962769_8a88bc972180963c": ("12S", 423),
    "SRR27962769_e2d376bc5b0f09a9": ("CO1_1", 352),
    "SRR27962769_e86a7ccf9e78adaf": ("CO1_1", 373),
}
EXPECTED_CONTEXTS = {
    "region_1_class_1": [(ITS_ID, "-", 72, 133)],
    "region_1_class_2": [],
    "region_2_class_1": [(ITS_ID, "+", 255, 317)],
    "region_2_class_2": [],
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def rc(sequence):
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[base]
                   for base in reversed(sequence))


def contexts(sequence, candidate_id, reference):
    return [(candidate_id, strand, start, start + len(sequence))
            for strand, query in (("+", sequence), ("-", rc(sequence)))
            for start in range(len(reference) - len(sequence) + 1)
            if reference[start:start + len(sequence)] == query]


def verify(args):
    archive_path = args.archive_dir / "evidence.tar.gz"
    contents_path = args.archive_dir / "ARCHIVE_CONTENTS.json"
    require(digest(archive_path.read_bytes()) == ARCHIVE_SHA, "Pinned archive hash")
    require(digest(contents_path.read_bytes()) == CONTENTS_SHA, "Pinned member manifest hash")
    sums = {}
    for line in (args.archive_dir / "SHA256SUMS").read_text().splitlines():
        expected_hash, name = line.split(maxsplit=1)
        require(name not in sums, "Duplicate SHA256SUMS entry")
        sums[name] = expected_hash
    require(sums == {"evidence.tar.gz": ARCHIVE_SHA, "ARCHIVE_CONTENTS.json": CONTENTS_SHA},
            "Outer SHA256SUMS entry set")
    entries = json.loads(contents_path.read_text())["members"]
    expected = {entry["path"]: entry for entry in entries}
    require(len(entries) == len(expected) == 205, "Expected archive member set")
    archived = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, "Unsafe archive name")
            if member.isdir():
                continue
            require(member.isfile() and member.name not in archived, "Archive member kind/uniqueness")
            require(member.name in expected, "Unexpected archive member")
            data = archive.extractfile(member).read()
            entry = expected[member.name]
            require(len(data) == entry["size_bytes"] and digest(data) == entry["sha256"],
                    f"Archive member mismatch: {member.name}")
            archived[member.name] = data
    require(set(archived) == set(expected), "Missing archive members")
    restored_paths = list(args.restore.rglob("*"))
    require(not any(path.is_symlink() for path in restored_paths), "Restore contains symbolic link")
    restored = {str(path.relative_to(args.restore)) for path in restored_paths if path.is_file()}
    require(restored == set(expected), "Restored regular-file set differs")
    for name, data in archived.items():
        require((args.restore / name).read_bytes() == data, f"Restored bytes differ: {name}")

    provenance = json.loads((args.inputs / "provenance.json").read_text())
    manifest_bytes = (args.inputs / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    require(provenance["archive_sha256"] == ARCHIVE_SHA and
            provenance["archive_contents_sha256"] == CONTENTS_SHA and
            provenance["verified_members"] == 205, "Provenance archive identities")
    require(provenance["derived_manifest_sha256"] == digest(manifest_bytes), "Derived manifest hash")
    require(provenance["script_sha256"] == digest((args.inputs / "restore_derive_inputs.py").read_bytes()),
            "Derivation helper hash")
    originals = provenance["originals"]
    require(len(originals) == 7 and len({row["archive_member"] for row in originals}) == 7,
            "Original input member count/uniqueness")
    for row in originals:
        data = archived[row["archive_member"]]
        require((args.inputs / row["path"]).read_bytes() == data and
                row["sha256"] == digest(data) and row["size_bytes"] == len(data),
                f"Original input copy mismatch: {row['path']}")
    source_co1 = json.loads(archived["inputs/co1_control_manifest.json"])
    full = {row["id"]: row for row in source_co1["candidates"]}
    require(set(full) == set(EXPECTED_DECOYS) and len(full) == len(source_co1["candidates"]),
            "Archived competitor identities")
    decoy_receipt = []
    for candidate_id, (gene, length) in EXPECTED_DECOYS.items():
        row = full[candidate_id]
        require(len(row["sequence"]) == length and digest(row["sequence"].encode()) == row["sequence_sha256"],
                f"Decoy sequence length/hash: {candidate_id}")
        require(all(member["gene"] == gene and member["length"] == length and
                    member["sequence_sha256"] == row["sequence_sha256"]
                    for member in row["source_candidates"]), "Competitor memberships")
        decoy_receipt.append({"id": candidate_id, "gene": gene, "length": length,
                              "sequence_sha256": row["sequence_sha256"],
                              "roles": sorted({m["role"] for m in row["source_candidates"]})})
    expected_decoys = [{key: row[key] for key in ("id", "sequence", "sequence_sha256")}
                       for row in source_co1["candidates"]]
    require(manifest["decoys"] == expected_decoys, "Derived competitor list differs from archive")
    historical = json.loads(archived["inputs/candidate.json"])["record"]
    require(historical["oriented_sequence"] == full[ITS_ID]["sequence"] and
            historical["sequence_sha256"] == full[ITS_ID]["sequence_sha256"], "Historical ITS identity")
    require([row["id"] for row in manifest["regions"]] == ["region_1", "region_2", "co1_control"],
            "Derived region identities/order")
    source_contexts = {row["class_id"]: row for row in provenance["classes"]}
    require(len(source_contexts) == len(provenance["classes"]) == 6, "Six unique source class records")
    class_receipt = []
    for region in manifest["regions"]:
        require(len(region["classes"]) == 2, "Two classes per region")
        is_its = region["id"] != "co1_control"
        source_member = (f"local_assays/{region['id']}.scanner_manifest.json"
                         if is_its else "inputs/co1_control_manifest.json")
        source = json.loads(archived[source_member])
        for index, row in enumerate(region["classes"]):
            if is_its:
                candidate = source["candidates"][index]
                event = next(e for e in source["events"] if e["candidate_id"] == candidate["id"])
                sequence, start, end = candidate["sequence"], 0, len(candidate["sequence"])
                require(row["id"] == candidate["id"] == f"{region['id']}_class_{index + 1}", "ITS class identity")
                require(event["start"] == 21 and event["end"] == end - 21, "ITS flank boundaries")
                observed_contexts = contexts(sequence, ITS_ID, historical["oriented_sequence"])
                require(observed_contexts == EXPECTED_CONTEXTS[row["id"]], "ITS full-context interval/strand")
                members = candidate["source"]["members"]
                require(sorted(m["stage"] for m in members) == ["post_prune", "pre_prune"], "Graph memberships")
                for member in members:
                    nodes = member["route_node_ids"]
                    require(all(len(node) == 18 for node in nodes), "Graph route node lengths")
                    require(all(a[1:] == b[:-1] for a, b in zip(nodes, nodes[1:])), "Graph route overlap")
                    route = nodes[0] + "".join(node[-1] for node in nodes[1:])
                    require(digest(route.encode()) == member["sequence_sha256"], "Graph route hash")
                    require((route if member["canonical_orientation"] == "+" else rc(route)) == sequence,
                            "Canonical graph route sequence")
                expected_kind = "historical_local_window" if index == 0 else "extra_cycle_local_window"
            else:
                event = source["events"][index]
                candidate = full[event["candidate_id"]]
                start, end = event["start"] - 21, event["end"] + 21
                sequence = candidate["sequence"][start:end]
                require((start, end) == ((289, 349) if index == 0 else (289, 370)), "CO1 window interval")
                require(row["id"] == event["id"], "CO1 class identity")
                observed_contexts = [(candidate["id"], "+", start, end)]
                expected_kind = event["role"]
            require(row["sequence"] == sequence and row["sequence_sha256"] == digest(sequence.encode()),
                    "Derived local class sequence/hash")
            require(row["candidate_id"] == candidate["id"] and row["window_start"] == start and
                    row["kind"] == expected_kind, "Derived class metadata")
            declared_contexts = [(c["decoy_id"], c["strand"], c["start"], c["end"])
                                 for c in row["self_context"]]
            require(declared_contexts == observed_contexts, "Derived self-context declarations")
            require(source_contexts[row["id"]] == {
                "class_id": row["id"], "source_member": source_member,
                "source_candidate_id": candidate["id"], "source_event_id": event["id"],
                "source_interval": [start, end], "self_context": row["self_context"]},
                "Class provenance differs")
            class_receipt.append({"id": row["id"], "length": len(sequence),
                                  "sequence_sha256": row["sequence_sha256"],
                                  "self_context": row["self_context"],
                                  "graph_memberships_verified": 2 if is_its else 0})
    return {"schema_version": 1, "status": "pass", "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "Independent archive, restore, literal class, graph membership, and finite-decoy checks; no biological reads inspected.",
            "argv": sys.argv, "verifier_sha256": digest(Path(__file__).read_bytes()),
            "derivation_helper_sha256": provenance["script_sha256"],
            "archive_sha256": ARCHIVE_SHA, "archive_contents_sha256": CONTENTS_SHA,
            "derived_manifest_sha256": digest(manifest_bytes),
            "provenance_sha256": digest((args.inputs / "provenance.json").read_bytes()),
            "archive_and_restore_members_verified": len(archived), "original_copies_verified": len(originals),
            "classes": class_receipt, "decoys": decoy_receipt}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--restore", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    require(not args.receipt.exists(), "Receipt already exists; choose a fresh output")
    receipt = verify(args)
    with args.receipt.open("x") as stream:
        json.dump(receipt, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps({"status": receipt["status"], "members": receipt["archive_and_restore_members_verified"],
                      "classes": len(receipt["classes"]), "decoys": len(receipt["decoys"])}))


if __name__ == "__main__":
    main()

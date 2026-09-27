#!/usr/bin/env python3
"""Verify archived graph lineage and derive unchanged ITS-only class inputs.

This preparation hashes archived files and sequences; it does not scan reads.
Outputs require a fresh directory and preserve byte-identical source members.
"""

import argparse
import hashlib
import json
from pathlib import Path
import tarfile


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
ARCHIVE = ROOT / "benchmarks/benchmark_results/v3.2-its-scc-audit-20260913/evidence.tar.gz"
CONTENTS = ARCHIVE.with_name("ARCHIVE_CONTENTS.json")
RESTORED = ROOT / "tmp/its-recruitment-20260925/prior"
PRIOR = ROOT / "benchmarks/benchmark_results/v3.2-its-recruitment-20260925/inputs/manifest.json"
ARCHIVE_SHA = "14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3"
CONTENTS_SHA = "dcf9821b93b1b993115890953b6c45a3ab82e7dda090e988100f2523a9deda10"
PRIOR_SHA = "8a61d8a1982fb5e4f2ad2cd40585a150f533066ad09949f9213cf66ebf1bc582"
SELECTED = {
    "graph_capture_v2/graphs/pre_prune.json": "graph_pre.json",
    "graph_capture_v2/graphs/post_prune.json": "graph_post.json",
    "graph_capture_v2/graphs/audit_manifest.json": "graph_archive_manifest.json",
    "local_assays/analysis.json": "scc_analysis.json",
    "local_assays/region_1.scanner_manifest.json": "originals/region_1.scanner_manifest.json",
    "local_assays/region_2.scanner_manifest.json": "originals/region_2.scanner_manifest.json",
    "inputs/expected_prefix.json": "expected_prefix.json",
}


def identity_bytes(raw):
    return {"size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def identity(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected regular nonsymlink input: {path}")
    before = path.stat()
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    after = path.stat()
    if any(getattr(before, key) != getattr(after, key) for key in
           ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")):
        raise ValueError(f"Input changed during hashing: {path}")
    return {"size_bytes": after.st_size, "sha256": hasher.hexdigest()}


def write_bytes(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as output:
        output.write(raw)


def write_json(path, value):
    write_bytes(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def prepare(output):
    if output.exists() or output.is_symlink():
        raise ValueError("Refusing to overwrite an input directory")
    pinned = {}
    for path, expected in ((ARCHIVE, ARCHIVE_SHA), (CONTENTS, CONTENTS_SHA), (PRIOR, PRIOR_SHA)):
        observed = identity(path)
        if observed["sha256"] != expected:
            raise ValueError(f"Pinned source identity mismatch: {path}")
        pinned[str(path.relative_to(ROOT))] = observed
    inventory_raw = CONTENTS.read_bytes()
    inventory = {item["path"]: {key: item[key] for key in ("sha256", "size_bytes")}
                 for item in json.loads(inventory_raw)["members"]}
    restored_bytes = {}
    for member in SELECTED:
        raw = (RESTORED / member).read_bytes()
        if identity_bytes(raw) != inventory[member]:
            raise ValueError(f"Restored member differs from archived inventory: {member}")
        restored_bytes[member] = raw
    found = set()
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        for member in archive:
            name = member.name.removeprefix("./")
            if name not in SELECTED:
                continue
            if name in found or not member.isfile():
                raise ValueError(f"Duplicate/nonregular selected archive member: {name}")
            if member.size != inventory[name]["size_bytes"]:
                raise ValueError(f"Archive member size differs: {name}")
            raw = archive.extractfile(member).read()
            if identity_bytes(raw) != inventory[name] or raw != restored_bytes[name]:
                raise ValueError(f"Archive/restored member mismatch: {name}")
            found.add(name)
    if found != set(SELECTED):
        raise ValueError("Missing selected archive members")
    original = PRIOR.read_bytes()
    if identity_bytes(original)["sha256"] != PRIOR_SHA:
        raise ValueError("Prior class manifest changed while reading")
    manifest = json.loads(original)
    regions = [item for item in manifest["regions"] if item["id"] in ("region_1", "region_2")]
    if len(regions) != 2 or sum(len(item["classes"]) for item in regions) != 4:
        raise ValueError("Expected exactly two regions and four unchanged ITS classes")
    for region in regions:
        for item in region["classes"]:
            if hashlib.sha256(item["sequence"].encode()).hexdigest() != item["sequence_sha256"]:
                raise ValueError("Local class sequence hash mismatch")
    derived = {"schema_version": 1, "regions": regions,
               "context_source": "prior_recruitment_manifest.json",
               "context_usage": "Archived self_context annotations and full candidates are lineage metadata only; no decoy matching is performed by this prototype."}
    # Bookend every external input before publishing the verified copies.
    for relative, observed in pinned.items():
        if identity(ROOT / relative) != observed:
            raise ValueError(f"Pinned source changed during preparation: {relative}")
    for member, raw in restored_bytes.items():
        if identity(RESTORED / member) != identity_bytes(raw):
            raise ValueError(f"Restored member changed during preparation: {member}")
    output.mkdir(parents=True)
    outputs = {}
    for member, destination in SELECTED.items():
        write_bytes(output / destination, restored_bytes[member])
        outputs[destination] = {"archive_member": member, **identity(output / destination)}
    write_bytes(output / "prior_recruitment_manifest.json", original)
    write_bytes(output / "archive_contents.json", inventory_raw)
    write_json(output / "its_manifest.json", derived)
    for name in ("prior_recruitment_manifest.json", "archive_contents.json", "its_manifest.json"):
        outputs[name] = identity(output / name)
    audit_manifest = json.loads(restored_bytes["graph_capture_v2/graphs/audit_manifest.json"])
    for phase in ("pre", "post"):
        name = f"graph_{phase}.json"
        if identity(output / name) != audit_manifest["files"][f"{phase}_prune.json"]:
            raise ValueError("Copied graph disagrees with graph export manifest")
    receipt = {"schema_version": 1, "status": "verified", "scope": "Archived input verification and literal ITS class selection only; no biological reads scanned",
               "prior_diagnostic_checkpoint": "69104d2", "sources": pinned,
               "restored_from": str(RESTORED.relative_to(ROOT)), "selected_archive_members_verified": len(found),
               "region_ids": [region["id"] for region in regions],
               "classes": [{"region_id": region["id"], "class_id": item["id"], "length": len(item["sequence"]),
                            "sequence_sha256": item["sequence_sha256"]} for region in regions for item in region["classes"]],
               "outputs": outputs, "preparation_source": identity(Path(__file__)),
               "decoy_limitation": derived["context_usage"]}
    write_json(output / "provenance.json", receipt)
    print(json.dumps({"status": "verified", "output": str(output), "members": len(found), "classes": 4}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "inputs")
    args = parser.parse_args()
    prepare(args.output_dir.resolve())

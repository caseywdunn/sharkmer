#!/usr/bin/env python3
"""Restore immutable historical metadata and derive the seven-loss comparison."""

import hashlib
import json
from pathlib import Path
import tarfile


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
PRIOR = ROOT / "benchmarks/benchmark_results/v3.2-withheld-diagnostics-20260913"
RELEASE = ROOT / "benchmarks/benchmark_results/v3.2-release-comparison-20260909"
SAMPLES = ("SRR27962769", "SRR31887760", "SRR1057608")
PINS = {
    PRIOR / "evidence.tar.gz": "04057bf295cf889c7c8db6fec59446923c2439e4be828d6d82b9286262574ae4",
    PRIOR / "ARCHIVE_CONTENTS.json": "4c947677b418d6fac4b98e4a45255e744e4fe1d82f9d99446ed2a2cc56dc5320",
    RELEASE / "inputs.json": "4763732cebacd43d3e295afe17ff903bdd1a6fc07449de24a03cd263b692cacf",
}
MEMBERS = ("protocol.json", "run_three_way.py", "test_run_three_way.py",
           "historical_source/candidate_manifest.json", "historical_lost_retained_map.json",
           "withheld_results.json", "changed_build/source_tree/panels/insecta.yaml")


def fingerprint(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected regular nonsymlink file: {path}")
    hasher = hashlib.sha256()
    before = path.stat()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    after = path.stat()
    if any(getattr(before, key) != getattr(after, key) for key in
           ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")):
        raise ValueError(f"File changed while hashing: {path}")
    return {"sha256": hasher.hexdigest(), "size_bytes": after.st_size}


def write_bytes(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as output:
        output.write(raw)


def write_json(path, value):
    write_bytes(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def main():
    output = HERE / "inputs"
    restored = ROOT / "tmp/uncertain-products-20260927/prior"
    if output.exists() or restored.exists():
        raise ValueError("Input or restore directory exists; refusing to overwrite")
    sources = {}
    for path, expected in PINS.items():
        observed = fingerprint(path)
        if observed["sha256"] != expected:
            raise ValueError(f"Pinned input differs: {path}")
        sources[str(path.relative_to(ROOT))] = observed
    inventory = json.loads((PRIOR / "ARCHIVE_CONTENTS.json").read_text())
    archive_bytes = {}
    with tarfile.open(PRIOR / "evidence.tar.gz", "r:gz") as archive:
        for member in archive:
            name = member.name.removeprefix("./")
            if name not in MEMBERS:
                continue
            if name in archive_bytes or not member.isfile():
                raise ValueError("Duplicate/nonregular selected archive member")
            raw = archive.extractfile(member).read()
            observed = {"size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
            if observed != inventory[name]:
                raise ValueError(f"Archive member differs: {name}")
            archive_bytes[name] = raw
    if set(archive_bytes) != set(MEMBERS):
        raise ValueError("Missing selected archive member")
    previous_protocol = json.loads(archive_bytes["protocol.json"])
    historical = json.loads(archive_bytes["historical_source/candidate_manifest.json"])
    losses = []
    for case in historical["cases"]:
        for item in case["lost_products"]:
            sequence = item["sequence"]
            if len(sequence) != item["length"] or hashlib.sha256(sequence.encode()).hexdigest() != item["sha256"]:
                raise ValueError("Historical lost sequence identity differs")
            losses.append({"accession": case["accession"], "taxon": case["taxon"], "gene": case["gene"],
                           "length": item["length"], "sequence": sequence, "sequence_sha256": item["sha256"],
                           "historical_header": item["header"],
                           "interpretation": "Historical assembled hypothesis; not an independent reference or validated haplotype"})
    if len(losses) != 7 or len({(item["accession"], item["sequence_sha256"]) for item in losses}) != 7:
        raise ValueError("Expected seven distinct historical loss identities")
    release_inputs = json.loads((RELEASE / "inputs.json").read_text())
    selected = {item["id"]: item for item in release_inputs["inputs"] if item["id"] in SAMPLES}
    inputs = []
    for accession in SAMPLES:
        sample = selected[accession]
        prefix = sample["historical_1m_prefix"]
        if prefix != previous_protocol["samples"][accession]["prefix"]:
            raise ValueError("Consumed prefix identities disagree between archived experiments")
        url = sample["source_urls_ordered"][0]
        if not url.endswith("_1.fastq.gz") or prefix["records"] != 1000000:
            raise ValueError("Expected the exact first million R1 records")
        local = ("data/its-recruitment-20260925/SRR27962769.first1M_R1.retry1.fastq" if accession == "SRR27962769"
                 else f"tmp/uncertain-products-20260927/reads/{accession}.first1M_R1.fastq")
        inputs.append({"accession": accession, "sample": previous_protocol["samples"][accession]["sample"],
                       "expected_prefix": prefix, "source_url": url, "path": local,
                       "archived_ena_metadata": sample["ena_metadata"],
                       "roles": previous_protocol["samples"][accession]["roles"]})
    panel_raw = (ROOT / "panels/insecta.yaml").read_bytes()
    if panel_raw != archive_bytes["changed_build/source_tree/panels/insecta.yaml"]:
        raise ValueError("Current full insect panel differs from prior frozen panel")
    for path in PINS:
        if fingerprint(path) != sources[str(path.relative_to(ROOT))]:
            raise ValueError("Pinned source changed during input preparation")
    output.mkdir(parents=True)
    restored.mkdir(parents=True)
    for name, raw in archive_bytes.items():
        write_bytes(restored / name, raw)
        write_bytes(output / "historical" / name, raw)
    write_bytes(output / "prior_archive_inventory.json", (PRIOR / "ARCHIVE_CONTENTS.json").read_bytes())
    write_bytes(output / "release_input_manifest.json", (RELEASE / "inputs.json").read_bytes())
    write_bytes(output / "insecta.yaml", panel_raw)
    write_json(output / "historical_losses.json", {"schema_version": 1, "rows": losses,
        "source_member": "historical_source/candidate_manifest.json",
        "scope": "All seven frozen loss identities; uncertain emission is separate from standard product recovery"})
    write_json(output / "read_inputs.json", {"schema_version": 1, "samples": inputs,
        "scope": "Exact consumed R1 prefixes only; original longer staged-file identities are not reproduced",
        "upstream_integrity_limit": "Prefix SHA256 verifies consumed bytes; partial download does not verify complete compressed ENA MD5"})
    files = {str(path.relative_to(output)): fingerprint(path) for path in sorted(output.rglob("*")) if path.is_file()}
    write_json(output / "provenance.json", {"schema_version": 1, "status": "verified", "sources": sources,
        "archive_members_verified": {name: inventory[name] for name in MEMBERS}, "outputs": files,
        "preparation_source": fingerprint(Path(__file__)), "biological_scans_performed": False})
    print(json.dumps({"status": "verified", "samples": len(inputs), "historical_losses": len(losses), "panel": fingerprint(output / "insecta.yaml")}))


if __name__ == "__main__":
    main()

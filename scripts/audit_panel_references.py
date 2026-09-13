#!/usr/bin/env python3
import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import yaml


SCHEMA_VERSION = 1
DEFAULT_MAX_SOURCE_BYTES = 5 * 1024 * 1024
DEFAULT_MAX_ARCHIVE_BYTES = 25 * 1024 * 1024
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_text(value):
    return sha256_bytes(value.encode())


def file_receipt(path):
    regular_path = Path(path).resolve()
    if not regular_path.is_file() or regular_path.is_symlink():
        raise ValueError(f"Expected regular file: {regular_path}")
    contents = regular_path.read_bytes()
    return {
        "path": str(regular_path),
        "size_bytes": len(contents),
        "sha256": sha256_bytes(contents),
    }


def accession_base(accession):
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def reverse_complement(sequence):
    if any(base not in "ACGTRYSWKMBDHVN" for base in sequence):
        raise ValueError("Reference has an unsupported non-IUPAC base")
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def all_positions(sequence, query):
    positions = []
    start_index = 0
    while True:
        position = sequence.find(query, start_index)
        if position < 0:
            return positions
        positions.append(position)
        start_index = position + 1


def panel_paths(panels_root):
    built_in_paths = sorted(
        path for path in panels_root.glob("*.yaml") if path.name != "reference.yaml"
    )
    example_path = panels_root / "examples/reference.yaml"
    if not example_path.is_file():
        raise ValueError(f"Missing reference example: {example_path}")
    return [("built_in", path) for path in built_in_paths] + [("example", example_path)]


def inventory_entries(panels_root):
    entries = []
    for panel_kind, panel_path in panel_paths(panels_root):
        panel = yaml.safe_load(panel_path.read_text())
        reference_groups = panel.get("references", []) if isinstance(panel, dict) else None
        if not isinstance(reference_groups, list):
            raise ValueError(f"Invalid references in panel: {panel_path}")
        panel_name = panel_path.stem
        for group_index, reference_group in enumerate(reference_groups):
            gene = reference_group.get("gene")
            sequences = reference_group.get("sequences")
            if not isinstance(gene, str) or not isinstance(sequences, list):
                raise ValueError(f"Invalid reference group in {panel_path}")
            for sequence_index, reference in enumerate(sequences):
                accession = reference.get("accession")
                taxon = reference.get("taxon")
                sequence = reference.get("sequence")
                if not all(isinstance(value, str) and value for value in (accession, taxon, sequence)):
                    raise ValueError(f"Reference must have nonempty accession, taxon, and sequence in {panel_path}")
                normalized_sequence = sequence.upper()
                entries.append({
                    "entry_id": f"{panel_name}:{group_index}:{sequence_index}",
                    "panel_kind": panel_kind,
                    "panel": panel_name,
                    "panel_path": str(panel_path),
                    "gene": gene,
                    "taxon": taxon,
                    "requested_accession": accession,
                    "requested_accession_base": accession_base(accession),
                    "panel_sequence_length": len(normalized_sequence),
                    "panel_sequence_sha256": sha256_text(normalized_sequence),
                    "panel_sequence": normalized_sequence,
                })
    return entries


def parse_metadata(path):
    metadata = json.loads(Path(path).read_text())
    result = metadata.get("result")
    if not isinstance(result, dict) or not isinstance(result.get("uids"), list):
        raise ValueError("Expected raw NCBI ESummary result object")
    records = {}
    for identifier in result["uids"]:
        record = result.get(str(identifier))
        if not isinstance(record, dict):
            raise ValueError(f"Metadata record missing for uid {identifier}")
        accession_version = record.get("accessionversion")
        if not isinstance(accession_version, str) or not accession_version:
            raise ValueError(f"Metadata record missing accession version for uid {identifier}")
        if accession_version in records:
            if records[accession_version] != record:
                raise ValueError(f"Conflicting metadata accession version: {accession_version}")
            continue
        records[accession_version] = record
    return records


def parse_fasta(path):
    records = {}
    header = None
    sequence_parts = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                accession_version = header.split(maxsplit=1)[0]
                if accession_version in records:
                    duplicate = {"header": header, "sequence": "".join(sequence_parts).upper()}
                    if records[accession_version] != duplicate:
                        raise ValueError(f"Conflicting FASTA accession version: {accession_version}")
                else:
                    records[accession_version] = {"header": header, "sequence": "".join(sequence_parts).upper()}
            header = line[1:]
            sequence_parts = []
        elif line.strip():
            sequence_parts.append(line.strip())
    if header is not None:
        accession_version = header.split(maxsplit=1)[0]
        if accession_version in records:
            duplicate = {"header": header, "sequence": "".join(sequence_parts).upper()}
            if records[accession_version] != duplicate:
                raise ValueError(f"Conflicting FASTA accession version: {accession_version}")
        else:
            records[accession_version] = {"header": header, "sequence": "".join(sequence_parts).upper()}
    return records


def resolve_metadata(requested_accession, metadata_records):
    if requested_accession in metadata_records:
        return metadata_records[requested_accession]
    if "." in requested_accession:
        return None
    matching = [
        record for accession_version, record in metadata_records.items()
        if accession_base(accession_version) == accession_base(requested_accession)
    ]
    if len(matching) == 1:
        return matching[0]
    if len(matching) > 1:
        return None
    return None


def match_candidates(panel_sequence, source_sequence, topology):
    if len(panel_sequence) > len(source_sequence):
        return []
    candidates = []
    orientations = (("forward", panel_sequence), ("reverse_complement", reverse_complement(panel_sequence)))
    for strand, query in orientations:
        for start_index in all_positions(source_sequence, query):
            candidates.append({
                "strand": "+" if strand == "forward" else "-",
                "start": start_index,
                "end": start_index + len(query),
                "wraps_origin": False,
                "kind": "exact_full" if len(query) == len(source_sequence) else "exact_subsequence",
                "region_sequence": query,
            })
        if topology.casefold() == "circular" and len(query) > 1:
            doubled_source = source_sequence + source_sequence[: len(query) - 1]
            for start_index in all_positions(doubled_source, query):
                if start_index >= len(source_sequence):
                    continue
                if start_index + len(query) <= len(source_sequence):
                    continue
                candidates.append({
                    "strand": "+" if strand == "forward" else "-",
                    "start": start_index,
                    "end": (start_index + len(query)) % len(source_sequence),
                    "wraps_origin": True,
                    "kind": "circular_wrap_subsequence",
                    "region_sequence": query,
                })
    unique = {}
    for candidate in candidates:
        key = (candidate["strand"], candidate["start"], candidate["end"], candidate["kind"])
        unique[key] = candidate
    return list(unique.values())


def taxon_relation(panel_taxon, source_organism):
    if not isinstance(source_organism, str) or not source_organism.strip():
        return "unavailable"
    if panel_taxon == source_organism:
        return "exact"
    if panel_taxon.casefold() == source_organism.casefold():
        return "casefold_equal"
    return "different"


def verify_entry(entry, metadata_records, fasta_records, max_source_bytes):
    record = resolve_metadata(entry["requested_accession"], metadata_records)
    output = {key: value for key, value in entry.items() if key != "panel_sequence"}
    if record is None:
        output["disposition"] = "unavailable_metadata"
        return output
    accession_version = record["accessionversion"]
    source_length = record.get("slen")
    output["resolved_accession_version"] = accession_version
    topology = "circular" if str(record.get("topology", "")).casefold() == "circular" else "linear"
    output["source_metadata"] = {
        "organism": record.get("organism"),
        "taxid": record.get("taxid"),
        "topology": topology,
        "declared_length": source_length,
    }
    output["taxon_relation"] = taxon_relation(entry["taxon"], record.get("organism"))
    if not isinstance(source_length, int) or source_length <= 0:
        output["disposition"] = "unavailable_source_length"
        return output
    if source_length > max_source_bytes:
        output["disposition"] = "skipped_oversize_source"
        return output
    source = fasta_records.get(accession_version)
    if source is None:
        output["disposition"] = "unavailable_source_bytes"
        return output
    source_sequence = source["sequence"]
    if len(source_sequence) != source_length:
        output["disposition"] = "unavailable_source_length_mismatch"
        return output
    output["source_fasta_header"] = source["header"]
    output["source_sequence_length"] = len(source_sequence)
    output["source_sequence_sha256"] = sha256_text(source_sequence)
    candidates = match_candidates(entry["panel_sequence"], source_sequence, topology)
    if candidates:
        candidates.sort(key=lambda candidate: (
            candidate["start"],
            candidate["end"],
            candidate["strand"],
            candidate["kind"],
        ))
        candidate = candidates[0]
        output["disposition"] = candidate["kind"]
        output["provenance"] = {
            "schema_version": 1,
            "kind": "public_record_region",
            "accession_version": accession_version,
            "source_sequence_sha256": sha256_text(source_sequence),
            "source_length": len(source_sequence),
            "start": candidate["start"],
            "end": candidate["end"],
            "strand": candidate["strand"],
            "wraps_origin": candidate["wraps_origin"],
            "sequence_sha256": sha256_text(entry["panel_sequence"]),
        }
        if len(candidates) > 1:
            output["nonunique_placement"] = True
            output["alternative_coordinates"] = [
                {key: value for key, value in alternative.items() if key != "region_sequence"}
                for alternative in candidates
            ]
        output["proposed_migration"] = {
            "accession": accession_version,
            "taxon": record.get("organism"),
            "sequence_sha256": sha256_text(entry["panel_sequence"]),
        }
    else:
        output["disposition"] = "no_exact_literal_match"
    return output


def source_records(metadata_records, fasta_records, max_source_bytes, source_url, retrieved_at):
    if not isinstance(source_url, str) or not source_url:
        raise ValueError("Source URL must be nonempty")
    if not isinstance(retrieved_at, str) or not retrieved_at:
        raise ValueError("Source retrieval timestamp must be nonempty")
    records = {}
    for accession_version, source in fasta_records.items():
        metadata = metadata_records.get(accession_version)
        if metadata is None:
            raise ValueError(f"FASTA source lacks metadata: {accession_version}")
        sequence = source["sequence"]
        declared_length = metadata.get("slen")
        if not isinstance(declared_length, int) or declared_length != len(sequence):
            raise ValueError(f"Source length receipt mismatch: {accession_version}")
        if len(sequence) > max_source_bytes:
            continue
        title = metadata.get("title")
        if not isinstance(title, str) or not title:
            raise ValueError(f"Source title is unavailable: {accession_version}")
        records[accession_version] = {
            "accession_version": accession_version,
            "sequence": sequence,
            "sequence_sha256": sha256_text(sequence),
            "organism": metadata.get("organism"),
            "taxid": metadata.get("taxid"),
            "topology": "circular" if str(metadata.get("topology", "")).casefold() == "circular" else "linear",
            "title": title,
            "url": source_url,
            "retrieved_at": retrieved_at,
        }
    return records


def build_catalog(panels_root, archive_dir, max_source_bytes, max_archive_bytes, source_url, retrieved_at):
    archive_path = Path(archive_dir)
    metadata_path = archive_path / "metadata.json"
    fasta_path = archive_path / "records.fasta"
    metadata_receipt = file_receipt(metadata_path)
    fasta_receipt = file_receipt(fasta_path)
    if fasta_receipt["size_bytes"] > max_archive_bytes:
        raise ValueError(f"FASTA archive exceeds {max_archive_bytes} byte cap")
    entries = inventory_entries(Path(panels_root))
    metadata_records = parse_metadata(metadata_path)
    fasta_records = parse_fasta(fasta_path)
    records = source_records(metadata_records, fasta_records, max_source_bytes, source_url, retrieved_at)
    verified = [
        verify_entry(entry, metadata_records, fasta_records, max_source_bytes)
        for entry in entries
    ]
    panel_receipts = {
        str(path.resolve()): file_receipt(path)
        for _, path in panel_paths(Path(panels_root))
    }
    dispositions = Counter(entry["disposition"] for entry in verified)
    return {
        "schema_version": SCHEMA_VERSION,
        "semantics": {
            "exact_match": "Literal panel IUPAC characters are compared directly; ambiguity codes are never wildcard-expanded.",
            "coordinate_system": "Coordinates are zero-based half-open on the downloaded source orientation. A circular-wrap region is source[start:] plus source[:end].",
            "scope": "A source sequence can rebuild only a recorded exact extracted region. Gene annotation and taxon labels are retained metadata, not inferred evidence.",
        },
        "source_archive": {
            "metadata": metadata_receipt,
            "fasta": fasta_receipt,
            "max_source_bytes": max_source_bytes,
            "max_archive_bytes": max_archive_bytes,
            "metadata_record_count": len(metadata_records),
            "fasta_record_count": len(fasta_records),
        },
        "records": records,
        "panel_receipts": panel_receipts,
        "summary": {
            "reference_entry_count": len(verified),
            "distinct_requested_accessions": len({entry["requested_accession"] for entry in entries}),
            "dispositions": dict(sorted(dispositions.items())),
            "taxon_relations": dict(sorted(Counter(entry.get("taxon_relation", "unavailable") for entry in verified).items())),
            "nonunique_placement_count": sum(bool(entry.get("nonunique_placement")) for entry in verified),
        },
        "entries": verified,
    }


def migration_document(catalog):
    verified = []
    quarantine = []
    for entry in catalog["entries"]:
        common = {
            "entry_id": entry["entry_id"],
            "panel": entry["panel"],
            "panel_path": entry["panel_path"],
            "gene": entry["gene"],
            "requested_accession": entry["requested_accession"],
            "legacy_taxon": entry["taxon"],
            "panel_sequence_sha256": entry["panel_sequence_sha256"],
            "disposition": entry["disposition"],
        }
        if "provenance" in entry:
            verified.append(common | {
                "public_taxon": entry["source_metadata"]["organism"],
                "public_taxid": entry["source_metadata"]["taxid"],
                "provenance": entry["provenance"],
                "proposed_migration": entry["proposed_migration"],
                "nonunique_placement": bool(entry.get("nonunique_placement")),
                "alternative_coordinates": entry.get("alternative_coordinates", []),
            })
        else:
            quarantine.append(common | {
                "source_metadata": entry.get("source_metadata"),
                "candidate_coordinates": entry.get("candidate_coordinates"),
            })
    return {
        "schema_version": 1,
        "purpose": "Proposed public-record reference migrations only; no panel edits are made by this audit.",
        "verified_migrations": verified,
        "quarantine": quarantine,
        "summary": {"verified_count": len(verified), "quarantine_count": len(quarantine)},
    }


def verify_existing(panels_root, catalog_path=None):
    try:
        from sharkmer_validate.reference_provenance import audit_references
    except ModuleNotFoundError:
        from scripts.sharkmer_validate.reference_provenance import audit_references

    results = []
    for panel_kind, panel_path in panel_paths(Path(panels_root)):
        panel = yaml.safe_load(panel_path.read_text())
        if not isinstance(panel, dict):
            raise ValueError(f"Invalid panel: {panel_path}")
        audit = audit_references(panel, catalog_path)
        results.append({
            "panel_kind": panel_kind,
            "panel_path": str(panel_path),
            "total_count": audit["total_count"],
            "verified_count": audit["verified_count"],
            "excluded_count": audit["excluded_count"],
            "catalog_status": audit["catalog"]["status"],
            "excluded": audit["excluded"],
        })
    failures = [
        result for result in results
        if result["catalog_status"] != "verified_local_snapshot" or result["excluded_count"]
    ]
    return {
        "panel_count": len(results),
        "total_references": sum(result["total_count"] for result in results),
        "verified_references": sum(result["verified_count"] for result in results),
        "excluded_references": sum(result["excluded_count"] for result in results),
        "failure_count": len(failures),
        "panels": results,
    }


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--panels-root", type=Path, default=Path("panels"))
    argument_parser.add_argument("--archive-dir", type=Path)
    argument_parser.add_argument("--output", type=Path)
    argument_parser.add_argument("--max-source-bytes", type=int, default=DEFAULT_MAX_SOURCE_BYTES)
    argument_parser.add_argument("--max-archive-bytes", type=int, default=DEFAULT_MAX_ARCHIVE_BYTES)
    argument_parser.add_argument("--source-url")
    argument_parser.add_argument("--retrieved-at")
    argument_parser.add_argument("--gzip-output", type=Path)
    argument_parser.add_argument("--migration-output", type=Path)
    argument_parser.add_argument("--verify-existing", action="store_true")
    argument_parser.add_argument("--reference-catalog", type=Path)
    arguments = argument_parser.parse_args()
    if arguments.verify_existing:
        verification = verify_existing(arguments.panels_root, arguments.reference_catalog)
        print(json.dumps(verification, sort_keys=True))
        if verification["failure_count"]:
            raise SystemExit(1)
        return
    if arguments.archive_dir is None or arguments.output is None:
        raise ValueError("--archive-dir and --output are required outside --verify-existing mode")
    if arguments.source_url is None or arguments.retrieved_at is None:
        raise ValueError("--source-url and --retrieved-at are required outside --verify-existing mode")
    if arguments.output.exists():
        raise ValueError(f"Output exists: {arguments.output}")
    catalog = build_catalog(
        arguments.panels_root,
        arguments.archive_dir,
        arguments.max_source_bytes,
        arguments.max_archive_bytes,
        arguments.source_url,
        arguments.retrieved_at,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")
    if arguments.gzip_output is not None:
        if arguments.gzip_output.exists():
            raise ValueError(f"Output exists: {arguments.gzip_output}")
        arguments.gzip_output.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(arguments.gzip_output, "wt") as output_stream:
            json.dump(catalog, output_stream, sort_keys=True, separators=(",", ":"))
    if arguments.migration_output is not None:
        if arguments.migration_output.exists():
            raise ValueError(f"Output exists: {arguments.migration_output}")
        arguments.migration_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.migration_output.write_text(json.dumps(migration_document(catalog), indent=2, sort_keys=True) + "\n")
    print(json.dumps(catalog["summary"], sort_keys=True))


if __name__ == "__main__":
    main()

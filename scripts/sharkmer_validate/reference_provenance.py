"""Verify reference extraction from a trusted, pinned public-record snapshot."""

import gzip
import hashlib
import json
import re
import zlib
from pathlib import Path

from . import runner


DEFAULT_CATALOG = runner.REPO_ROOT / "panels" / "reference_sources.json.gz"
MAX_CATALOG_BYTES = 25 * 1024 * 1024
MAX_SOURCE_BASES = 5 * 1024 * 1024
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
ACCESSION_VERSION = re.compile(r"[A-Z][A-Z0-9_]*\.[0-9]+\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")


def sequence_digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate catalog key: {key}")
        result[key] = value
    return result


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _sequence(value, name):
    if not isinstance(value, str) or not value or set(value) - DNA_ALPHABET:
        raise ValueError(f"{name} must contain literal uppercase DNA/IUPAC symbols")
    return value


def _digest(value, name):
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def load_catalog(catalog_path=None):
    path = Path(catalog_path) if catalog_path is not None else DEFAULT_CATALOG
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"reference catalog is not a regular file: {path}")
    with path.open("rb") as stream:
        compressed = stream.read(MAX_CATALOG_BYTES + 1)
    if len(compressed) > MAX_CATALOG_BYTES:
        raise ValueError("compressed reference catalog exceeds size limit")
    if compressed.startswith(b"\x1f\x8b"):
        import io

        with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
            content = stream.read(MAX_CATALOG_BYTES + 1)
    else:
        content = compressed
    if len(content) > MAX_CATALOG_BYTES:
        raise ValueError("decompressed reference catalog exceeds size limit")
    catalog = json.loads(content, object_pairs_hook=_unique_object)
    if not isinstance(catalog, dict) or type(catalog.get("schema_version")) is not int:
        raise ValueError("reference catalog schema_version must be integer 1")
    if catalog["schema_version"] != 1 or not isinstance(catalog.get("records"), dict):
        raise ValueError("unsupported reference catalog schema")
    for accession, record in catalog["records"].items():
        if not isinstance(accession, str) or ACCESSION_VERSION.fullmatch(accession) is None:
            raise ValueError("catalog accession must be versioned")
        if not isinstance(record, dict) or record.get("accession_version") != accession:
            raise ValueError(f"catalog accession binding differs: {accession}")
        sequence = _sequence(record.get("sequence"), "source sequence")
        if len(sequence) > MAX_SOURCE_BASES:
            raise ValueError(f"source exceeds size limit: {accession}")
        if _digest(record.get("sequence_sha256"), "source digest") != sequence_digest(sequence):
            raise ValueError(f"source checksum differs: {accession}")
        if not isinstance(record.get("topology"), str) or record["topology"] not in {"linear", "circular"}:
            raise ValueError(f"source topology is missing or unsupported: {accession}")
        if not isinstance(record.get("organism"), str) or not record["organism"].strip():
            raise ValueError(f"source organism is missing: {accession}")
        _integer(record.get("taxid"), "source taxid", 1)
        conflicts = record.get("annotation_conflicts", {})
        if not isinstance(conflicts, dict) or any(
            not isinstance(gene, str) or not gene or not isinstance(reason, str) or not reason
            for gene, reason in conflicts.items()
        ):
            raise ValueError(f"invalid reviewed gene annotation conflicts: {accession}")
        for name in ("url", "retrieved_at"):
            if not isinstance(record.get(name), str) or not record[name]:
                raise ValueError(f"source {name} is missing: {accession}")
    receipt = {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(compressed).hexdigest(),
        "status": "verified_local_snapshot",
        "record_count": len(catalog["records"]),
        "trust_basis": "pinned_local_catalog; not live public-record authentication",
    }
    return catalog, receipt


def extract_region(record, provenance):
    source = _sequence(record.get("sequence"), "source sequence")
    start = _integer(provenance.get("start"), "start")
    end = _integer(provenance.get("end"), "end")
    if start >= len(source) or end > len(source):
        raise ValueError("reference coordinates are outside the source")
    wraps = provenance.get("wraps_origin")
    if type(wraps) is not bool:
        raise ValueError("wraps_origin must be a boolean")
    if wraps:
        if record.get("topology") != "circular" or start < end:
            raise ValueError("origin wrapping requires a circular source and start >= end")
        region = source[start:] + source[:end]
    else:
        if end <= start:
            raise ValueError("linear extraction requires end > start")
        region = source[start:end]
    if not region or len(region) > len(source):
        raise ValueError("reference extraction must span at most one source traversal")
    strand = provenance.get("strand")
    if strand not in {"+", "-"}:
        raise ValueError("reference strand must be + or -")
    return reverse_complement(region) if strand == "-" else region


def verify_reference(reference, catalog):
    provenance = reference.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("missing public-record provenance")
    if type(provenance.get("schema_version")) is not int or provenance["schema_version"] != 1:
        raise ValueError("unsupported reference provenance schema")
    if provenance.get("kind") != "public_record_region":
        raise ValueError("reference source is not a public_record_region")
    accession = provenance.get("accession_version")
    if not isinstance(accession, str) or ACCESSION_VERSION.fullmatch(accession) is None:
        raise ValueError("reference accession must be versioned")
    if reference.get("accession") != accession:
        raise ValueError("reference accession differs from provenance accession")
    record = catalog["records"].get(accession)
    if record is None:
        raise ValueError("accession is absent from the trusted reference catalog")
    source_digest = _digest(provenance.get("source_sequence_sha256"), "source digest")
    if source_digest != record["sequence_sha256"]:
        raise ValueError("reference source checksum differs from catalog")
    source_length = _integer(provenance.get("source_length"), "source_length", 1)
    if source_length != len(record["sequence"]):
        raise ValueError("reference source length differs from catalog")
    sequence = _sequence(reference.get("sequence"), "panel reference sequence")
    reconstructed = extract_region(record, provenance)
    if reconstructed != sequence:
        raise ValueError("panel reference differs from the exact extracted public region")
    if _digest(provenance.get("sequence_sha256"), "region digest") != sequence_digest(reconstructed):
        raise ValueError("reference region checksum differs")
    if reference.get("taxon") != record["organism"]:
        raise ValueError("reference taxon differs from the public source organism")
    return {
        "taxon": record["organism"],
        "accession": accession,
        "sequence": reconstructed,
        "provenance": dict(provenance),
        "source_taxid": record["taxid"],
        "gene_assignment_basis": "panel_annotation; not established by sequence provenance alone",
        "contains_ambiguity": bool(set(reconstructed) - set("ACGT")),
    }


def audit_references(panel_data, catalog_path=None):
    path = Path(catalog_path) if catalog_path is not None else DEFAULT_CATALOG
    catalog_error = None
    try:
        catalog, receipt = load_catalog(path)
    except (OSError, ValueError, UnicodeError, EOFError, zlib.error) as error:
        catalog = None
        catalog_error = str(error)
        receipt = {"path": str(path.resolve()), "sha256": None, "status": "unavailable", "error": catalog_error}
    verified = []
    excluded = []
    for group_index, group in enumerate(panel_data.get("references") or []):
        gene_name = runner.derive_gene_name(group)
        for sequence_index, reference in enumerate(group.get("sequences") or []):
            identity = {
                "gene_name": gene_name,
                "taxon": reference.get("taxon"),
                "accession": reference.get("accession"),
                "group_index": group_index,
                "sequence_index": sequence_index,
            }
            try:
                if catalog_error is not None:
                    raise ValueError(catalog_error)
                result = verify_reference(reference, catalog)
                logical_gene = re.sub(r"_\d+$", "", gene_name).split("-", 1)[0].casefold()
                conflicts = catalog["records"][result["accession"]].get("annotation_conflicts", {})
                conflict = next((reason for gene, reason in conflicts.items() if gene.casefold() == logical_gene), None)
                if conflict:
                    raise ValueError(f"source region verified but gene annotation excluded: {conflict}")
            except (ValueError, TypeError, KeyError) as error:
                excluded.append({**identity, "status": "unverified_reference", "reason": str(error)})
            else:
                verified.append({**identity, **result})
    result = {
        "schema_version": 1,
        "verified": verified,
        "excluded": excluded,
        "verified_count": len(verified),
        "excluded_count": len(excluded),
        "total_count": len(verified) + len(excluded),
        "catalog": receipt,
        "coordinate_system": "zero_based_half_open; extract forward region then reverse complement",
        "biological_truth": "public-region provenance only; not sample haplotype truth",
    }
    result["audit_sha256"] = hashlib.sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return result


def verified_references(panel_data, catalog_path=None):
    return audit_references(panel_data, catalog_path)["verified"]

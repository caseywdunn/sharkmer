"""Assess recovered amplicons against verified external reference sequences."""

import hashlib
import json
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path

from . import reference_provenance
from .reference_targets import (
    logical_gene_name,
    reference_target_names,
    target_logical_genes,
)

MIN_QUERY_COVERAGE_PCT = 90.0
MIN_IDENTITY_PCT = 90.0
MIN_CONFLICT_COVERAGE_PCT = 20.0


@dataclass
class RefBlastResult:
    status: str
    expected_gene: str
    expected_taxon: str
    matched_gene: str | None = None
    matched_taxon: str | None = None
    matched_accession: str | None = None
    target_support: str = "unavailable"
    taxon_support: str = "unavailable"
    sequence_relationship: str = "unavailable"
    target_gene_absence: str = "not_established"
    haplotype_truth: str = "not_established"
    read_support: str = "not_evaluated"
    primer_region_support: str = "not_established"
    expected_logical_gene: str | None = None
    matched_logical_gene: str | None = None
    pct_identity: float | None = None
    identity_count: int | None = None
    align_length: int | None = None
    query_length: int | None = None
    query_coverage_pct: float | None = None
    aligned_fraction: float | None = None
    query_unaligned_bases: int | None = None
    unmatched_query_regions: list[list[int]] | None = None
    reference_length: int | None = None
    reference_coverage_pct: float | None = None
    unmatched_reference_regions: list[list[int]] | None = None
    gap_count: int | None = None
    query_gap_bases: int | None = None
    reference_gap_bases: int | None = None
    inter_hsp_query_gap_bases: int | None = None
    inter_hsp_reference_gap_bases: int | None = None
    hsp_count: int = 0
    coherent_hsp_count: int = 0
    split_alignment: bool | None = None
    chimeric_alignment: bool | None = None
    alignment_structure_status: str = "not_evaluated"
    alignment_evidence: list[dict] | None = None
    on_target: bool = False
    thresholds: dict | None = None
    reference_provenance: dict | None = None
    reference_database_provenance: dict | None = None
    error: str | None = None

    def __post_init__(self):
        if self.expected_logical_gene is None:
            self.expected_logical_gene = logical_gene_name(self.expected_gene)
        if self.matched_logical_gene is None and self.matched_gene:
            self.matched_logical_gene = logical_gene_name(self.matched_gene)


def check_blastn_available() -> bool:
    try:
        return subprocess.run(
            ["blastn", "-version"], capture_output=True, text=True
        ).returncode == 0
    except FileNotFoundError:
        return False


def check_makeblastdb_available() -> bool:
    try:
        return subprocess.run(
            ["makeblastdb", "-version"], capture_output=True, text=True
        ).returncode == 0
    except FileNotFoundError:
        return False


def extract_references(panel_data: dict, catalog_path: Path | None = None) -> list[dict]:
    return reference_provenance.verified_references(panel_data, catalog_path)


def _audit_summary(audit: dict, panel_data: dict | None = None) -> dict:
    target_mapping = target_logical_genes(panel_data or {})
    target_mapping_sha256 = hashlib.sha256(
        json.dumps(target_mapping, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    verified = []
    for reference in audit.get("verified", []):
        sequence = reference["sequence"].upper().encode()
        verified.append(
            {
                "gene": reference["gene_name"],
                "logical_gene": logical_gene_name(
                    reference["gene_name"], panel_data=panel_data
                ),
                "taxon": reference["taxon"],
                "accession": reference["accession"],
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence).hexdigest(),
                "provenance": reference.get("provenance"),
                "source_taxid": reference.get("source_taxid"),
                "gene_assignment_basis": reference.get("gene_assignment_basis"),
                "contains_ambiguity": reference.get("contains_ambiguity"),
            }
        )
    return {
        "verified": verified,
        "excluded": audit.get("excluded", []),
        "verified_count": len(verified),
        "excluded_count": len(audit.get("excluded", [])),
        "total_count": audit.get("total_count"),
        "catalog": audit.get("catalog"),
        "audit_sha256": audit.get("audit_sha256"),
        "coordinate_system": audit.get("coordinate_system"),
        "biological_truth": audit.get("biological_truth"),
        "target_logical_genes": target_mapping,
        "target_logical_genes_sha256": target_mapping_sha256,
        "primer_region_support": "not_established",
    }


def reference_checksums(panel_data: dict, catalog_path: Path | None = None) -> dict:
    return _audit_summary(
        reference_provenance.audit_references(panel_data, catalog_path), panel_data
    )


def available_reference_targets(
    panel_data: dict, catalog_path: Path | None = None
) -> set[str]:
    references = extract_references(panel_data, catalog_path)
    return reference_target_names(
        panel_data, {reference["gene_name"] for reference in references}
    )


def _metadata_path(db_prefix: Path) -> Path:
    return Path(f"{db_prefix}.reference_metadata.json")


def _metadata_checksum_path(db_prefix: Path) -> Path:
    return Path(f"{db_prefix}.reference_metadata.sha256")


def _database_manifest_path(db_prefix: Path) -> Path:
    return Path(f"{db_prefix}.database_manifest.json")


def _database_manifest_checksum_path(db_prefix: Path) -> Path:
    return Path(f"{db_prefix}.database_manifest.sha256")


def _file_receipt(path: Path) -> dict:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"BLAST database artifact is not a regular file: {path}")
    return {
        "path": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def build_reference_db(
    panel_data: dict, tmpdir: Path, catalog_path: Path | None = None
) -> Path | None:
    audit = reference_provenance.audit_references(panel_data, catalog_path)
    references = audit["verified"]
    if not references:
        return None
    if not check_blastn_available() or not check_makeblastdb_available():
        print("WARNING: blastn or makeblastdb not available; skipping reference BLAST.")
        return None
    fasta_path = tmpdir / "references.fasta"
    metadata = {}
    with fasta_path.open("w") as output:
        for reference_index, reference in enumerate(references):
            subject_id = f"reference_{reference_index:06d}"
            sequence = reference["sequence"].upper()
            metadata[subject_id] = {
                "gene": reference["gene_name"],
                "logical_gene": logical_gene_name(
                    reference["gene_name"], panel_data=panel_data
                ),
                "taxon": reference["taxon"],
                "accession": reference["accession"],
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "provenance": {
                    "source_region": reference.get("provenance"),
                    "source_taxid": reference.get("source_taxid"),
                    "gene_assignment_basis": reference.get("gene_assignment_basis"),
                    "contains_ambiguity": reference.get("contains_ambiguity"),
                },
            }
            output.write(f">{subject_id}\n")
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset : offset + 80] + "\n")
    db_prefix = tmpdir / "ref_db"
    completed = subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_prefix)],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        print(f"WARNING: makeblastdb failed: {completed.stderr}")
        return None
    Path(f"{db_prefix}.reference_count").write_text(str(len(references)))
    metadata_path = _metadata_path(db_prefix)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    _metadata_checksum_path(db_prefix).write_text(
        hashlib.sha256(metadata_path.read_bytes()).hexdigest() + "\n"
    )
    database_artifacts = sorted(
        path
        for path in tmpdir.glob(f"{db_prefix.name}.*")
        if path
        not in {
            _database_manifest_path(db_prefix),
            _database_manifest_checksum_path(db_prefix),
        }
    )
    if not database_artifacts:
        raise ValueError("makeblastdb did not produce database artifacts")
    database_manifest = {
        "schema_version": 1,
        "database_prefix": db_prefix.name,
        "reference_fasta": _file_receipt(fasta_path),
        "artifacts": [_file_receipt(path) for path in database_artifacts],
        "reference_audit": _audit_summary(audit, panel_data),
    }
    manifest_path = _database_manifest_path(db_prefix)
    manifest_path.write_text(json.dumps(database_manifest, indent=2, sort_keys=True) + "\n")
    _database_manifest_checksum_path(db_prefix).write_text(
        hashlib.sha256(manifest_path.read_bytes()).hexdigest() + "\n"
    )
    return db_prefix


def _load_metadata(db_path: Path) -> tuple[dict, dict]:
    metadata_path = _metadata_path(db_path)
    checksum_path = _metadata_checksum_path(db_path)
    if (
        not metadata_path.is_file()
        or metadata_path.is_symlink()
        or not checksum_path.is_file()
        or checksum_path.is_symlink()
    ):
        raise ValueError("Reference metadata provenance is missing")
    observed_sha256 = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
    if checksum_path.read_text().strip() != observed_sha256:
        raise ValueError("Reference metadata checksum differs")
    metadata = json.loads(metadata_path.read_text())
    if not isinstance(metadata, dict):
        raise ValueError("Reference metadata provenance is invalid")
    for subject_id, reference in metadata.items():
        if not isinstance(subject_id, str) or not isinstance(reference, dict):
            raise ValueError("Reference metadata entry is invalid")
        if any(
            not isinstance(reference.get(field), str) or not reference[field]
            for field in ("gene", "taxon", "accession", "sha256")
        ):
            raise ValueError("Reference metadata identity is invalid")
        if type(reference.get("length")) is not int or reference["length"] <= 0:
            raise ValueError("Reference metadata length is invalid")
        if not isinstance(reference.get("provenance"), dict):
            raise ValueError("Reference metadata public provenance is invalid")
    return metadata, {
        "metadata_path": str(metadata_path),
        "metadata_sha256": observed_sha256,
        "reference_count_path": str(Path(f"{db_path}.reference_count")),
    }


def _load_database_manifest(db_path: Path) -> tuple[dict, dict]:
    manifest_path = _database_manifest_path(db_path)
    checksum_path = _database_manifest_checksum_path(db_path)
    if (
        not manifest_path.is_file()
        or manifest_path.is_symlink()
        or not checksum_path.is_file()
        or checksum_path.is_symlink()
    ):
        raise ValueError("BLAST database manifest provenance is missing")
    manifest_sha256 = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if checksum_path.read_text().strip() != manifest_sha256:
        raise ValueError("BLAST database manifest checksum differs")
    manifest = json.loads(manifest_path.read_text())
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != 1
        or manifest.get("database_prefix") != db_path.name
        or not isinstance(manifest.get("artifacts"), list)
        or not manifest["artifacts"]
    ):
        raise ValueError("BLAST database manifest is invalid")
    expected_names = set()
    for receipt in [manifest.get("reference_fasta"), *manifest["artifacts"]]:
        name = receipt.get("path") if isinstance(receipt, dict) else None
        if (
            not isinstance(name, str)
            or Path(name).name != name
            or name in expected_names
        ):
            raise ValueError("BLAST database artifact receipt is invalid")
        expected_names.add(name)
        path = manifest_path.parent / name
        observed = _file_receipt(path)
        if observed != receipt:
            raise ValueError(f"BLAST database artifact checksum differs: {name}")
    observed_database_names = {
        path.name
        for path in manifest_path.parent.glob(f"{db_path.name}.*")
        if path
        not in {
            manifest_path,
            checksum_path,
        }
    }
    expected_database_names = {
        receipt["path"] for receipt in manifest["artifacts"]
    }
    if observed_database_names != expected_database_names:
        raise ValueError("BLAST database artifact set differs from manifest")
    return manifest, {
        "manifest_path": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "reference_audit_sha256": manifest.get("reference_audit", {}).get("audit_sha256"),
    }


def _validate_metadata_audit(metadata: dict, manifest: dict) -> None:
    reference_audit = manifest.get("reference_audit", {})
    verified = reference_audit.get("verified")
    target_mapping = reference_audit.get("target_logical_genes", {})
    if not isinstance(verified, list) or len(verified) != len(metadata):
        raise ValueError("BLAST metadata and verified-reference audit counts differ")
    for reference_index, expected in enumerate(verified):
        expected_logical_gene = logical_gene_name(
            expected.get("gene"), target_mapping=target_mapping
        )
        if expected.get("logical_gene") != expected_logical_gene:
            raise ValueError(
                "BLAST reference audit logical target differs from panel context"
            )
        observed = metadata.get(f"reference_{reference_index:06d}")
        if observed is None:
            raise ValueError("BLAST metadata subject ordering differs from reference audit")
        comparable = {
            "gene": observed.get("gene"),
            "logical_gene": observed.get("logical_gene"),
            "taxon": observed.get("taxon"),
            "accession": observed.get("accession"),
            "length": observed.get("length"),
            "sha256": observed.get("sha256"),
            "provenance": observed.get("provenance", {}).get("source_region"),
            "source_taxid": observed.get("provenance", {}).get("source_taxid"),
            "gene_assignment_basis": observed.get("provenance", {}).get("gene_assignment_basis"),
            "contains_ambiguity": observed.get("provenance", {}).get("contains_ambiguity"),
        }
        if comparable != expected:
            raise ValueError("BLAST metadata identity differs from verified-reference audit")


def blast_against_references(
    sequence: str,
    db_path: Path,
    expected_gene: str,
    expected_taxon: str,
    min_query_coverage_pct: float = MIN_QUERY_COVERAGE_PCT,
    min_identity_pct: float = MIN_IDENTITY_PCT,
) -> RefBlastResult:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as query:
        query.write(f">query\n{sequence}\n")
        query_path = query.name
    count_path = Path(f"{db_path}.reference_count")
    try:
        if not count_path.is_file() or count_path.is_symlink():
            raise ValueError("Reference-count provenance is missing")
        reference_count = int(count_path.read_text())
        metadata, metadata_provenance = _load_metadata(db_path)
        database_manifest, database_provenance = _load_database_manifest(db_path)
        _validate_metadata_audit(metadata, database_manifest)
        if reference_count <= 0 or reference_count != len(metadata):
            raise ValueError("Reference-count and metadata provenance disagree")
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
        Path(query_path).unlink(missing_ok=True)
        return RefBlastResult(
            "failed_run",
            expected_gene,
            expected_taxon,
            error=f"{error}; complete hit retrieval cannot be verified",
        )
    try:
        completed = subprocess.run(
            [
                "blastn", "-db", str(db_path), "-query", query_path, "-outfmt", "5",
                "-evalue", "1e-10", "-num_alignments", str(reference_count),
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        return RefBlastResult("failed_run", expected_gene, expected_taxon, error="BLAST timed out")
    finally:
        Path(query_path).unlink(missing_ok=True)
    if completed.returncode != 0:
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            error=f"blastn failed: {completed.stderr[:200]}",
        )
    try:
        metadata_after, metadata_provenance_after = _load_metadata(db_path)
        database_manifest_after, database_provenance_after = _load_database_manifest(db_path)
        _validate_metadata_audit(metadata_after, database_manifest_after)
        count_after = int(count_path.read_text())
        if (
            metadata_after != metadata
            or metadata_provenance_after != metadata_provenance
            or database_manifest_after != database_manifest
            or database_provenance_after != database_provenance
            or count_after != reference_count
        ):
            raise ValueError("BLAST database provenance changed during query")
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            error=f"{error}; BLAST output was not classified",
        )
    database_provenance = {**database_provenance, **metadata_provenance}
    target_mapping = database_manifest.get("reference_audit", {}).get(
        "target_logical_genes", {}
    )
    return _parse_blast_xml(
        completed.stdout,
        expected_gene,
        expected_taxon,
        min_query_coverage_pct,
        min_identity_pct,
        metadata,
        database_provenance,
        expected_logical_gene=logical_gene_name(
            expected_gene, target_mapping=target_mapping
        ),
    )


def _interval_union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged = []
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged


def _interval_union_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in _interval_union(intervals))


def _unmatched_regions(length: int, intervals: list[tuple[int, int]]) -> list[list[int]]:
    if length <= 0:
        return []
    regions = []
    cursor = 1
    for start, end in _interval_union(intervals):
        if cursor < start:
            regions.append([cursor, start - 1])
        cursor = max(cursor, end + 1)
    if cursor <= length:
        regions.append([cursor, length])
    return regions


def _legacy_metadata(hit_definition: str, hit_length: int) -> dict:
    parts = hit_definition.split("|", 2)
    metadata = {
        "gene": parts[0] if parts else None,
        "taxon": parts[1].replace("_", " ") if len(parts) >= 2 else None,
        "accession": parts[2] if len(parts) >= 3 else None,
        "length": hit_length or None,
        "provenance": {"status": "explicit_test_fixture"},
    }
    metadata["logical_gene"] = logical_gene_name(metadata["gene"])
    return metadata


def _hsp_orientation(hsp: dict) -> int:
    query_direction = 1 if hsp["query_to"] >= hsp["query_from"] else -1
    subject_direction = 1 if hsp["subject_to"] >= hsp["subject_from"] else -1
    return query_direction * subject_direction


def _chain_score(chain: list[dict]) -> tuple:
    return (
        sum(hsp["bit_score"] for hsp in chain),
        sum(hsp["identity_count"] for hsp in chain),
        _interval_union_length([hsp["query_interval"] for hsp in chain]),
        -sum(hsp["gap_count"] for hsp in chain),
    )


def _can_follow(previous: dict, current: dict) -> bool:
    if previous["orientation"] != current["orientation"]:
        return False
    if previous["query_interval"][1] >= current["query_interval"][0]:
        return False
    if current["orientation"] > 0:
        return previous["subject_interval"][1] < current["subject_interval"][0]
    return previous["subject_interval"][0] > current["subject_interval"][1]


def _best_chain(hsps: list[dict]) -> list[dict]:
    ordered = sorted(hsps, key=lambda hsp: (hsp["query_interval"], hsp["subject_interval"]))
    chains = []
    for current_index, current in enumerate(ordered):
        best = [current]
        for previous_index in range(current_index):
            candidate = chains[previous_index]
            if _can_follow(candidate[-1], current):
                extended = candidate + [current]
                if _chain_score(extended) > _chain_score(best):
                    best = extended
        chains.append(best)
    return max(chains, key=_chain_score, default=[])


def _parse_hit(hit: ET.Element, query_length: int, metadata_by_id: dict | None) -> dict:
    hit_definition = hit.findtext("Hit_def", "")
    hit_length = int(hit.findtext("Hit_len", "0"))
    metadata = (metadata_by_id or {}).get(hit_definition)
    if metadata is None:
        if metadata_by_id is not None:
            raise ValueError(f"BLAST subject is absent from reference metadata: {hit_definition}")
        metadata = _legacy_metadata(hit_definition, hit_length)
    elif hit_length != metadata.get("length"):
        raise ValueError(f"BLAST subject length differs from reference metadata: {hit_definition}")
    if not metadata.get("logical_gene"):
        metadata = {**metadata, "logical_gene": logical_gene_name(metadata.get("gene"))}
    hsps = []
    for hsp_index, hsp in enumerate(hit.findall(".//Hsp")):
        align_length = int(hsp.findtext("Hsp_align-len", "0"))
        identity_count = int(hsp.findtext("Hsp_identity", "0"))
        query_from = int(hsp.findtext("Hsp_query-from", "0"))
        query_to = int(hsp.findtext("Hsp_query-to", "0"))
        subject_from = int(hsp.findtext("Hsp_hit-from", "0"))
        subject_to = int(hsp.findtext("Hsp_hit-to", "0"))
        query_sequence = hsp.findtext("Hsp_qseq", "")
        subject_sequence = hsp.findtext("Hsp_hseq", "")
        if (
            align_length <= 0
            or identity_count < 0
            or identity_count > align_length
            or min(query_from, query_to) < 1
            or max(query_from, query_to) > query_length
            or min(subject_from, subject_to) < 1
            or (
                metadata.get("length") is not None
                and max(subject_from, subject_to) > metadata["length"]
            )
        ):
            raise ValueError(f"BLAST HSP coordinates or counts are invalid: {hit_definition}")
        parsed = {
            "hsp_index": hsp_index,
            "align_length": align_length,
            "pct_identity": 100.0 * identity_count / align_length if align_length else 0.0,
            "identity_count": identity_count,
            "bit_score": float(hsp.findtext("Hsp_bit-score", "0")),
            "query_interval": (min(query_from, query_to), max(query_from, query_to)),
            "subject_interval": (min(subject_from, subject_to), max(subject_from, subject_to)),
            "query_from": query_from,
            "query_to": query_to,
            "subject_from": subject_from,
            "subject_to": subject_to,
            "gap_count": int(hsp.findtext("Hsp_gaps", "0")),
            "query_gap_bases": query_sequence.count("-") if query_sequence else None,
            "reference_gap_bases": subject_sequence.count("-") if subject_sequence else None,
        }
        parsed["orientation"] = _hsp_orientation(parsed)
        hsps.append(parsed)
    chain = _best_chain(hsps)
    query_intervals = [hsp["query_interval"] for hsp in chain]
    subject_intervals = [hsp["subject_interval"] for hsp in chain]
    inter_hsp_query_gaps = 0
    inter_hsp_reference_gaps = 0
    for previous, current in zip(chain, chain[1:]):
        inter_hsp_query_gaps += max(
            current["query_interval"][0] - previous["query_interval"][1] - 1, 0
        )
        if current["orientation"] > 0:
            inter_hsp_reference_gaps += max(
                current["subject_interval"][0] - previous["subject_interval"][1] - 1, 0
            )
        else:
            inter_hsp_reference_gaps += max(
                previous["subject_interval"][0] - current["subject_interval"][1] - 1, 0
            )
    query_aligned = _interval_union_length(query_intervals)
    reference_length = metadata.get("length") or hit_length or None
    align_length = sum(hsp["align_length"] for hsp in chain)
    identity_count = sum(hsp["identity_count"] for hsp in chain)
    return {
        **metadata,
        "subject_id": hit_definition,
        "hsps": hsps,
        "chain": chain,
        "chain_hsp_indices": [hsp["hsp_index"] for hsp in chain],
        "bit_score": sum(hsp["bit_score"] for hsp in chain),
        "align_length": align_length,
        "identity_count": identity_count,
        "pct_identity": 100.0 * identity_count / align_length if align_length else 0.0,
        "query_coverage_pct": 100.0 * query_aligned / query_length if query_length else 0.0,
        "query_unaligned_bases": max(query_length - query_aligned, 0),
        "unmatched_query_regions": _unmatched_regions(query_length, query_intervals),
        "reference_coverage_pct": (
            100.0 * _interval_union_length(subject_intervals) / reference_length
            if reference_length
            else None
        ),
        "unmatched_reference_regions": (
            _unmatched_regions(reference_length, subject_intervals) if reference_length else None
        ),
        "gap_count": (
            sum(hsp["gap_count"] for hsp in chain)
            + inter_hsp_query_gaps
            + inter_hsp_reference_gaps
        ),
        "query_gap_bases": (
            sum(hsp["query_gap_bases"] for hsp in chain)
            if all(hsp["query_gap_bases"] is not None for hsp in chain)
            else None
        ),
        "reference_gap_bases": (
            sum(hsp["reference_gap_bases"] for hsp in chain)
            if all(hsp["reference_gap_bases"] is not None for hsp in chain)
            else None
        ),
        "inter_hsp_query_gap_bases": inter_hsp_query_gaps,
        "inter_hsp_reference_gap_bases": inter_hsp_reference_gaps,
        "all_query_union_coverage_pct": (
            100.0 * _interval_union_length([hsp["query_interval"] for hsp in hsps]) / query_length
            if query_length
            else 0.0
        ),
    }


def _hit_score(hit: dict) -> tuple:
    return (
        hit["bit_score"], hit["identity_count"], hit["query_coverage_pct"], -hit["gap_count"]
    )


def _complementary_conflict(
    hits: list[dict],
    query_length: int,
    min_query_coverage_pct: float,
    min_identity_pct: float,
) -> tuple[bool, bool]:
    if query_length <= 0:
        return False, False
    significant = []
    for hit in hits:
        intervals = [
            hsp["query_interval"]
            for hsp in hit["hsps"]
            if hsp["pct_identity"] >= min_identity_pct
            and 100.0 * _interval_union_length([hsp["query_interval"]]) / query_length
            >= MIN_CONFLICT_COVERAGE_PCT
        ]
        if intervals:
            significant.append((hit, intervals))
    cross_reference = False
    cross_gene = False
    for first_index, (first_hit, first_intervals) in enumerate(significant):
        for second_hit, second_intervals in significant[first_index + 1 :]:
            first_union = _interval_union_length(first_intervals)
            second_union = _interval_union_length(second_intervals)
            combined_union = _interval_union_length(first_intervals + second_intervals)
            added = combined_union - max(first_union, second_union)
            if (
                100.0 * combined_union / query_length >= min_query_coverage_pct
                and 100.0 * added / query_length >= MIN_CONFLICT_COVERAGE_PCT
            ):
                cross_reference = True
                if first_hit.get("logical_gene") != second_hit.get("logical_gene"):
                    cross_gene = True
    return cross_reference, cross_gene


def _relationship(hit: dict, structural: bool, qualifying: bool) -> str:
    if structural:
        return "structurally_conflicting"
    if not qualifying:
        return "partial_unresolved"
    if (
        hit["query_coverage_pct"] == 100.0
        and hit["reference_coverage_pct"] == 100.0
        and hit["pct_identity"] == 100.0
        and hit["gap_count"] == 0
        and hit["query_unaligned_bases"] == 0
    ):
        return "reference_identical"
    if hit["pct_identity"] == 100.0 and hit["gap_count"] == 0:
        return "partial_unresolved"
    return "aligned_differences"


def _parse_blast_xml(
    xml_text: str,
    expected_gene: str,
    expected_taxon: str,
    min_query_coverage_pct: float = MIN_QUERY_COVERAGE_PCT,
    min_identity_pct: float = MIN_IDENTITY_PCT,
    reference_metadata: dict | None = None,
    database_provenance: dict | None = None,
    fixture_provenance: str | None = None,
    expected_logical_gene: str | None = None,
) -> RefBlastResult:
    thresholds = {
        "min_query_coverage_pct": min_query_coverage_pct,
        "min_identity_pct": min_identity_pct,
        "min_complementary_conflict_coverage_pct": MIN_CONFLICT_COVERAGE_PCT,
        "coverage_rule": "best_collinear_nonoverlapping_hsp_chain",
        "criteria_scope": "recorded_marker_support_criterion_not_universal_biological_truth",
    }
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as error:
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            thresholds=thresholds, error=f"XML parse error: {error}",
        )
    iteration = root.find(".//Iteration")
    query_length = int(iteration.findtext("Iteration_query-len", "0")) if iteration is not None else 0
    if reference_metadata is None and fixture_provenance != "explicit_test_fixture":
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            query_length=query_length or None, thresholds=thresholds,
            error="Verified reference metadata is required for biological-reference classification",
        )
    try:
        hits = [
            _parse_hit(hit, query_length, reference_metadata)
            for hit in root.findall(".//Iteration/Iteration_hits/Hit")
        ]
    except (TypeError, ValueError) as error:
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            query_length=query_length or None, thresholds=thresholds,
            error=f"Invalid BLAST alignment evidence: {error}",
        )
    hits = [hit for hit in hits if hit["chain"]]
    if not hits:
        return RefBlastResult(
            "no_significant_hit", expected_gene, expected_taxon,
            target_support="insufficient", taxon_support="insufficient",
            sequence_relationship="partial_unresolved", query_length=query_length or None,
            thresholds=thresholds, alignment_structure_status="no_significant_alignment",
        )
    qualifying = [
        hit for hit in hits
        if hit["query_coverage_pct"] >= min_query_coverage_pct
        and hit["pct_identity"] >= min_identity_pct
    ]
    incoherent = any(
        len(hit["hsps"]) > 1
        and hit["query_coverage_pct"] < min_query_coverage_pct
        and hit["all_query_union_coverage_pct"] >= min_query_coverage_pct
        for hit in hits
    )
    cross_reference, cross_gene = _complementary_conflict(
        hits, query_length, min_query_coverage_pct, min_identity_pct
    )
    structural = not qualifying and incoherent
    multi_reference_unresolved = not qualifying and cross_reference
    selected_pool = qualifying or hits
    selected = max(selected_pool, key=_hit_score)
    selected_score = _hit_score(selected)
    tied = [hit for hit in selected_pool if _hit_score(hit) == selected_score]
    expected_logical_gene = expected_logical_gene or logical_gene_name(expected_gene)
    selected_gene_matches = selected.get("logical_gene") == expected_logical_gene
    tied_genes = {hit.get("logical_gene") for hit in tied}
    tied_taxa = {(hit.get("taxon") or "").lower() for hit in tied}
    if structural:
        status = "structurally_conflicting"
        target_support = "ambiguous"
        taxon_support = "ambiguous"
    elif multi_reference_unresolved and cross_gene:
        status = "ambiguous_gene"
        target_support = "ambiguous"
        taxon_support = "ambiguous"
    elif not qualifying:
        status = "insufficient_alignment"
        target_support = "insufficient"
        taxon_support = "insufficient"
    elif len(tied_genes) > 1:
        status = "ambiguous_gene"
        target_support = "ambiguous"
        taxon_support = "ambiguous"
    elif not selected_gene_matches:
        status = "wrong_gene"
        target_support = "conflicting"
        taxon_support = "unavailable"
    elif not expected_taxon:
        status = "gene_supported_taxon_not_evaluated"
        target_support = "supported"
        taxon_support = "not_evaluated"
    elif len(tied_taxa) > 1:
        status = "ambiguous_taxon"
        target_support = "supported"
        taxon_support = "ambiguous"
    elif (selected.get("taxon") or "").lower() != expected_taxon.lower():
        status = "gene_supported_other_taxon"
        target_support = "supported"
        taxon_support = "other_taxon"
    else:
        status = "gene_supported_expected_taxon"
        target_support = "supported"
        taxon_support = "supported"
    qualifying_selected = selected in qualifying
    evidence = []
    for hit in hits:
        for hsp in hit["hsps"]:
            evidence.append(
                {
                    "subject_id": hit["subject_id"],
                    "gene": hit.get("gene"),
                    "logical_gene": hit.get("logical_gene"),
                    "taxon": hit.get("taxon"),
                    "accession": hit.get("accession"),
                    "selected_in_coherent_chain": (
                        hit is selected and hsp["hsp_index"] in selected["chain_hsp_indices"]
                    ),
                    **{
                        key: value
                        for key, value in hsp.items()
                        if key not in {"query_interval", "subject_interval"}
                    },
                }
            )
    if structural:
        structure_status = "structurally_conflicting_alignment_evidence"
    elif multi_reference_unresolved:
        structure_status = "complementary_multi_reference_evidence_unresolved"
    elif len(selected["chain"]) > 1:
        structure_status = "coherent_collinear_multi_hsp"
    elif qualifying_selected:
        structure_status = "coherent_single_hsp"
    else:
        structure_status = "partial_alignment"
    reference_coverage = selected["reference_coverage_pct"]
    return RefBlastResult(
        status, expected_gene, expected_taxon,
        matched_gene=selected.get("gene"), matched_taxon=selected.get("taxon"),
        matched_accession=selected.get("accession"), target_support=target_support,
        taxon_support=taxon_support,
        expected_logical_gene=expected_logical_gene,
        matched_logical_gene=selected.get("logical_gene"),
        sequence_relationship=_relationship(selected, structural, qualifying_selected),
        pct_identity=round(selected["pct_identity"], 3),
        identity_count=selected["identity_count"], align_length=selected["align_length"],
        query_length=query_length or None,
        query_coverage_pct=round(selected["query_coverage_pct"], 3),
        aligned_fraction=round(selected["query_coverage_pct"] / 100.0, 6),
        query_unaligned_bases=selected["query_unaligned_bases"],
        unmatched_query_regions=selected["unmatched_query_regions"],
        reference_length=selected.get("length"),
        reference_coverage_pct=round(reference_coverage, 3) if reference_coverage is not None else None,
        unmatched_reference_regions=selected["unmatched_reference_regions"],
        gap_count=selected["gap_count"], query_gap_bases=selected["query_gap_bases"],
        reference_gap_bases=selected["reference_gap_bases"], hsp_count=len(selected["hsps"]),
        inter_hsp_query_gap_bases=selected["inter_hsp_query_gap_bases"],
        inter_hsp_reference_gap_bases=selected["inter_hsp_reference_gap_bases"],
        coherent_hsp_count=len(selected["chain"]),
        split_alignment=structural or multi_reference_unresolved or len(selected["chain"]) > 1,
        chimeric_alignment=None,
        alignment_structure_status=structure_status, alignment_evidence=evidence,
        on_target=status == "gene_supported_expected_taxon", thresholds=thresholds,
        reference_provenance=selected.get("provenance"),
        reference_database_provenance=database_provenance,
    )


def _unevaluated_result(
    status: str,
    gene: str,
    taxon: str,
    error: str | None = None,
    expected_logical_gene: str | None = None,
) -> dict:
    support = "unavailable" if status == "no_verified_reference" else "not_evaluated"
    return asdict(
        RefBlastResult(
            status,
            gene,
            taxon,
            target_support=support,
            taxon_support=support,
            expected_logical_gene=expected_logical_gene,
            error=error,
        )
    )


def _summarize_product_matches(products: list[dict]) -> dict | None:
    matches = [product.get("reference_match") for product in products]
    matches = [match for match in matches if match]
    if not matches:
        return None
    priority = {
        "failed_run": 0,
        "wrong_gene": 1,
        "ambiguous_gene": 2,
        "structurally_conflicting": 3,
        "ambiguous_taxon": 4,
        "insufficient_alignment": 5,
        "no_significant_hit": 6,
        "not_evaluated": 7,
        "no_verified_reference": 8,
        "gene_supported_other_taxon": 9,
        "gene_supported_taxon_not_evaluated": 10,
        "gene_supported_expected_taxon": 11,
    }
    summary = min(matches, key=lambda match: priority.get(match.get("status"), -1)).copy()
    summary["all_products_gene_supported"] = all(
        match.get("target_support") == "supported" for match in matches
    )
    summary["all_products_expected_taxon_supported"] = all(
        match.get("target_support") == "supported"
        and match.get("taxon_support") == "supported"
        for match in matches
    )
    summary["n_products_evaluated"] = len(matches)
    summary["on_target"] = summary["all_products_expected_taxon_supported"]
    summary["haplotype_truth"] = "not_established"
    summary["read_support"] = "not_evaluated"
    return summary


def blast_all_products(
    run_results: list,
    db_path: Path | None,
    sample_taxon: str,
    skip_blast: bool = False,
    reference_genes: set[str] | None = None,
    target_mapping: dict[str, str] | None = None,
):
    requested_target_mapping = dict(target_mapping or {})
    target_mapping = requested_target_mapping
    database_error = None
    if db_path is not None:
        try:
            manifest, _ = _load_database_manifest(db_path)
            reference_audit = manifest.get("reference_audit", {})
            stored_target_mapping = dict(
                reference_audit.get("target_logical_genes", {})
            )
            if (
                requested_target_mapping
                and requested_target_mapping != stored_target_mapping
            ):
                raise ValueError(
                    "BLAST database target mapping differs from current panel context"
                )
            target_mapping = stored_target_mapping
            reference_logical_genes = {
                reference.get("logical_gene")
                for reference in reference_audit.get("verified", [])
            }
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
            database_error = str(error)
            reference_logical_genes = set()
    else:
        reference_logical_genes = {
            logical_gene_name(gene, target_mapping=target_mapping)
            for gene in (reference_genes or set())
        }
    count = 0
    for run in run_results:
        if not run.get("success"):
            continue
        for gene_result in run.get("genes", []):
            gene = gene_result["gene"]
            products = gene_result.get("products", [])
            for product in products:
                count += 1
                expected_logical_gene = logical_gene_name(
                    gene, target_mapping=target_mapping
                )
                if database_error is not None:
                    match = _unevaluated_result(
                        "failed_run",
                        gene,
                        sample_taxon,
                        f"{database_error}; BLAST output was not classified",
                        expected_logical_gene,
                    )
                elif expected_logical_gene not in reference_logical_genes:
                    match = _unevaluated_result(
                        "no_verified_reference",
                        gene,
                        sample_taxon,
                        expected_logical_gene=expected_logical_gene,
                    )
                elif skip_blast:
                    match = _unevaluated_result(
                        "not_evaluated",
                        gene,
                        sample_taxon,
                        "BLAST disabled",
                        expected_logical_gene,
                    )
                elif db_path is None:
                    match = _unevaluated_result(
                        "not_evaluated", gene, sample_taxon,
                        "Verified reference database unavailable",
                        expected_logical_gene,
                    )
                else:
                    match = asdict(
                        blast_against_references(
                            product["sequence"], db_path,
                            expected_gene=gene, expected_taxon=sample_taxon,
                        )
                    )
                product["reference_match"] = match
                print(f"  {gene} product {product.get('product_index', '?')}: {match['status']}")
            gene_result["reference_match"] = _summarize_product_matches(products)
    if count == 0:
        print("No amplicons to validate.")

"""Validate every recovered amplicon against panel reference sequences."""

import hashlib
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from pathlib import Path

from . import runner

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
    pct_identity: float | None = None
    align_length: int | None = None
    query_length: int | None = None
    query_coverage_pct: float | None = None
    aligned_fraction: float | None = None
    hsp_count: int = 0
    split_alignment: bool | None = None
    chimeric_alignment: bool | None = None
    alignment_structure_status: str = "not_evaluated"
    alignment_evidence: list[dict] | None = None
    on_target: bool = False
    thresholds: dict | None = None
    error: str | None = None


def check_blastn_available() -> bool:
    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_makeblastdb_available() -> bool:
    try:
        result = subprocess.run(["makeblastdb", "-version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def extract_references(panel_data: dict) -> list:
    refs = []
    for gene_block in panel_data.get("references", []):
        gene_name = runner.derive_gene_name(gene_block)
        for seq_entry in gene_block.get("sequences", []):
            refs.append(
                {
                    "gene_name": gene_name,
                    "taxon": seq_entry["taxon"],
                    "accession": seq_entry.get("accession", "unknown"),
                    "sequence": seq_entry["sequence"],
                }
            )
    return refs


def reference_checksums(panel_data: dict) -> list[dict]:
    checksums = []
    for reference in extract_references(panel_data):
        sequence = reference["sequence"].upper().encode()
        checksums.append(
            {
                "gene": reference["gene_name"],
                "taxon": reference["taxon"],
                "accession": reference["accession"],
                "length": len(sequence),
                "sha256": hashlib.sha256(sequence).hexdigest(),
            }
        )
    return checksums


def _sanitize_for_fasta_header(value: str) -> str:
    return value.replace("|", "_").replace(" ", "_").replace("/", "_")


def build_reference_db(panel_data: dict, tmpdir: Path) -> Path | None:
    refs = extract_references(panel_data)
    if not refs:
        return None
    if not check_blastn_available() or not check_makeblastdb_available():
        print("WARNING: blastn or makeblastdb not available; skipping reference BLAST.")
        return None
    fasta_path = tmpdir / "references.fasta"
    with open(fasta_path, "w") as output:
        for reference in refs:
            gene = _sanitize_for_fasta_header(reference["gene_name"])
            taxon = _sanitize_for_fasta_header(reference["taxon"])
            accession = _sanitize_for_fasta_header(reference["accession"])
            output.write(f">{gene}|{taxon}|{accession}\n")
            sequence = reference["sequence"]
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset : offset + 80] + "\n")
    db_prefix = tmpdir / "ref_db"
    result = subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_prefix)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"WARNING: makeblastdb failed: {result.stderr}")
        return None
    Path(f"{db_prefix}.reference_count").write_text(str(len(refs)))
    return db_prefix


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
    reference_count_path = Path(f"{db_path}.reference_count")
    if not reference_count_path.exists():
        Path(query_path).unlink(missing_ok=True)
        return RefBlastResult(
            "failed_run",
            expected_gene,
            expected_taxon,
            error="Reference-count provenance is missing; complete hit retrieval cannot be verified",
        )
    reference_count = int(reference_count_path.read_text())
    try:
        result = subprocess.run(
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
    if result.returncode != 0:
        return RefBlastResult(
            "failed_run", expected_gene, expected_taxon,
            error=f"blastn failed: {result.stderr[:200]}",
        )
    return _parse_blast_xml(
        result.stdout, expected_gene, expected_taxon, min_query_coverage_pct, min_identity_pct
    )


def _interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    sorted_intervals = sorted(intervals)
    merged_length = 0
    current_start, current_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged_length += current_end - current_start + 1
            current_start, current_end = start, end
    return merged_length + current_end - current_start + 1


def _parse_hit(hit: ET.Element, query_length: int) -> dict:
    parts = hit.findtext("Hit_def", "").split("|")
    hsps = []
    for hsp in hit.findall(".//Hsp"):
        align_length = int(hsp.findtext("Hsp_align-len", "0"))
        identity_count = int(hsp.findtext("Hsp_identity", "0"))
        bit_score = float(hsp.findtext("Hsp_bit-score", "0"))
        query_from = int(hsp.findtext("Hsp_query-from", "0"))
        query_to = int(hsp.findtext("Hsp_query-to", "0"))
        subject_from = int(hsp.findtext("Hsp_hit-from", "0"))
        subject_to = int(hsp.findtext("Hsp_hit-to", "0"))
        query_span = abs(query_to - query_from) + 1
        hsps.append(
            {
                "align_length": align_length,
                "pct_identity": 100.0 * identity_count / align_length if align_length else 0.0,
                "identity_count": identity_count,
                "bit_score": bit_score,
                "query_coverage_pct": 100.0 * query_span / query_length if query_length else 0.0,
                "interval": (min(query_from, query_to), max(query_from, query_to)),
                "query_from": query_from,
                "query_to": query_to,
                "subject_from": subject_from,
                "subject_to": subject_to,
            }
        )
    best = max(
        hsps,
        key=lambda value: (value["query_coverage_pct"], value["pct_identity"]),
        default=None,
    )
    union_length = _interval_union_length([value["interval"] for value in hsps])
    return {
        "gene": parts[0] if parts else None,
        "taxon": parts[1].replace("_", " ") if len(parts) >= 2 else None,
        "accession": parts[2].replace("_", " ") if len(parts) >= 3 else None,
        "hsps": hsps,
        "best": best,
        "union_coverage_pct": 100.0 * union_length / query_length if query_length else 0.0,
    }


def _parse_blast_xml(
    xml_text: str,
    expected_gene: str,
    expected_taxon: str,
    min_query_coverage_pct: float = MIN_QUERY_COVERAGE_PCT,
    min_identity_pct: float = MIN_IDENTITY_PCT,
) -> RefBlastResult:
    thresholds = {
        "min_query_coverage_pct": min_query_coverage_pct,
        "min_identity_pct": min_identity_pct,
        "coverage_rule": "single_contiguous_hsp",
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
    hits = [_parse_hit(hit, query_length) for hit in root.findall(".//Iteration/Iteration_hits/Hit")]
    hits = [hit for hit in hits if hit["best"] is not None]
    if not hits:
        return RefBlastResult(
            "no_significant_hit", expected_gene, expected_taxon,
            query_length=query_length or None, thresholds=thresholds,
            alignment_structure_status="no_significant_alignment",
        )
    qualifying = [
        hit for hit in hits
        if hit["best"]["query_coverage_pct"] >= min_query_coverage_pct
        and hit["best"]["pct_identity"] >= min_identity_pct
    ]
    expected_gene_hits = [
        hit for hit in qualifying if (hit["gene"] or "").lower() == expected_gene.lower()
    ]
    expected_taxon_hits = [
        hit for hit in expected_gene_hits
        if (hit["taxon"] or "").lower() == expected_taxon.lower()
    ]
    significant_by_reference = []
    alignment_evidence = []
    for hit in hits:
        significant_intervals = []
        for hsp in hit["hsps"]:
            alignment_evidence.append(
                {
                    "gene": hit["gene"],
                    "taxon": hit["taxon"],
                    "accession": hit["accession"],
                    **{key: value for key, value in hsp.items() if key != "interval"},
                }
            )
            if (
                hsp["pct_identity"] >= min_identity_pct
                and hsp["query_coverage_pct"] >= MIN_CONFLICT_COVERAGE_PCT
            ):
                significant_intervals.append(hsp["interval"])
        if significant_intervals:
            significant_by_reference.append(
                {
                    "identity": (hit["gene"], hit["taxon"], hit["accession"]),
                    "gene": hit["gene"],
                    "intervals": significant_intervals,
                }
            )
    within_reference_split = any(
        len(hit["hsps"]) > 1
        and hit["best"]["query_coverage_pct"] < min_query_coverage_pct
        and hit["union_coverage_pct"] >= min_query_coverage_pct
        for hit in hits
    )
    cross_reference_split = False
    cross_gene_split = False
    for first_index, first in enumerate(significant_by_reference):
        for second in significant_by_reference[first_index + 1 :]:
            first_union = _interval_union_length(first["intervals"])
            second_union = _interval_union_length(second["intervals"])
            combined_union = _interval_union_length(first["intervals"] + second["intervals"])
            added_by_combining = combined_union - max(first_union, second_union)
            complementary = (
                query_length > 0
                and 100.0 * combined_union / query_length >= min_query_coverage_pct
                and 100.0 * added_by_combining / query_length >= MIN_CONFLICT_COVERAGE_PCT
            )
            if complementary:
                cross_reference_split = True
                if (first["gene"] or "").lower() != (second["gene"] or "").lower():
                    cross_gene_split = True
    split_alignment = within_reference_split or cross_reference_split
    chimeric_alignment = cross_gene_split

    def hit_score(hit):
        strongest_hsp = max(
            hit["hsps"],
            key=lambda hsp: (
                hsp["bit_score"], hsp["identity_count"], hsp["query_coverage_pct"]
            ),
        )
        return (
            strongest_hsp["bit_score"],
            strongest_hsp["identity_count"],
            strongest_hsp["query_coverage_pct"],
        )

    if chimeric_alignment or split_alignment:
        selected = max(hits, key=lambda hit: (hit["best"]["query_coverage_pct"], hit["best"]["pct_identity"]))
        status = "split_or_chimeric_alignment"
    elif qualifying:
        selected = max(qualifying, key=hit_score)
        selected_score = hit_score(selected)
        tied = [hit for hit in qualifying if hit_score(hit) == selected_score]
        selected_gene_matches = (selected["gene"] or "").lower() == expected_gene.lower()
        tied_other_gene = any(
            (hit["gene"] or "").lower() != expected_gene.lower() for hit in tied
        )
        tied_expected_gene = any(
            (hit["gene"] or "").lower() == expected_gene.lower() for hit in tied
        )
        tied_other_taxon = any(
            (hit["taxon"] or "").lower() != expected_taxon.lower() for hit in tied
        )
        tied_expected_taxon = any(
            (hit["taxon"] or "").lower() == expected_taxon.lower() for hit in tied
        )
        if tied_other_gene and tied_expected_gene:
            status = "ambiguous_gene"
        elif not selected_gene_matches:
            status = "wrong_gene"
        elif tied_other_taxon and tied_expected_taxon:
            status = "ambiguous_taxon"
        elif (selected["taxon"] or "").lower() != expected_taxon.lower():
            status = "confirmed_gene_other_taxon"
        else:
            status = "confirmed_product"
    else:
        selected = max(hits, key=lambda hit: (hit["best"]["query_coverage_pct"], hit["best"]["pct_identity"]))
        status = "insufficient_alignment"
    best = selected["best"]
    return RefBlastResult(
        status, expected_gene, expected_taxon,
        matched_gene=selected["gene"], matched_taxon=selected["taxon"],
        matched_accession=selected["accession"], pct_identity=round(best["pct_identity"], 3),
        align_length=best["align_length"], query_length=query_length or None,
        query_coverage_pct=round(best["query_coverage_pct"], 3),
        aligned_fraction=round(best["query_coverage_pct"] / 100.0, 6),
        hsp_count=len(selected["hsps"]), split_alignment=split_alignment,
        chimeric_alignment=chimeric_alignment, on_target=status == "confirmed_product",
        thresholds=thresholds,
        alignment_structure_status="split_evidence" if split_alignment else "single_hsp_gate",
        alignment_evidence=alignment_evidence,
    )


def _unevaluated_result(status: str, gene: str, taxon: str, error: str | None = None) -> dict:
    return asdict(
        RefBlastResult(
            status, gene, taxon, error=error,
        )
    )


def _summarize_product_matches(products: list[dict]) -> dict | None:
    matches = [product.get("reference_match") for product in products]
    matches = [match for match in matches if match]
    if not matches:
        return None
    status_priority = {
        "failed_run": 0,
        "wrong_gene": 1,
        "ambiguous_gene": 2,
        "ambiguous_taxon": 2,
        "split_or_chimeric_alignment": 3,
        "insufficient_alignment": 4,
        "no_significant_hit": 5,
        "not_evaluated": 6,
        "no_reference": 7,
        "confirmed_gene_other_taxon": 8,
        "confirmed_product": 9,
    }
    summary = min(matches, key=lambda match: status_priority.get(match.get("status"), -1)).copy()
    summary["all_products_confirmed"] = all(
        match.get("status") == "confirmed_product" for match in matches
    )
    summary["n_products_evaluated"] = len(matches)
    summary["on_target"] = summary["all_products_confirmed"]
    return summary


def blast_all_products(
    run_results: list,
    db_path: Path | None,
    sample_taxon: str,
    skip_blast: bool = False,
    reference_genes: set[str] | None = None,
):
    reference_genes = reference_genes or set()
    count = 0
    for run in run_results:
        if not run.get("success"):
            continue
        for gene_result in run.get("genes", []):
            gene = gene_result["gene"]
            products = gene_result.get("products", [])
            for product in products:
                count += 1
                if gene not in reference_genes:
                    match = _unevaluated_result("no_reference", gene, sample_taxon)
                elif skip_blast:
                    match = _unevaluated_result("not_evaluated", gene, sample_taxon, "BLAST disabled")
                elif db_path is None:
                    match = _unevaluated_result(
                        "not_evaluated", gene, sample_taxon, "Reference database unavailable"
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

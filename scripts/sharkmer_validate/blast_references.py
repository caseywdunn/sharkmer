"""BLAST amplicons against gold-standard reference sequences from panel YAML.

Replaces the old NCBI nt BLAST workflow with a fast, local, specific approach:
compile panel reference sequences into a temporary BLAST DB, then query amplicons
against it. Results report identity AND whether the match is on-target (correct
gene and taxon).
"""

import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from . import runner


@dataclass
class RefBlastResult:
    """Result of BLASTing an amplicon against the reference DB."""

    matched_gene: str | None = None
    matched_taxon: str | None = None
    matched_accession: str | None = None
    pct_identity: float | None = None
    align_length: int | None = None
    on_target: bool = False
    error: str | None = None


def check_blastn_available() -> bool:
    """Check if blastn is installed and accessible."""
    try:
        result = subprocess.run(
            ["blastn", "-version"], capture_output=True, text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def check_makeblastdb_available() -> bool:
    """Check if makeblastdb is installed and accessible."""
    try:
        result = subprocess.run(
            ["makeblastdb", "-version"], capture_output=True, text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def extract_references(panel_data: dict) -> list:
    """Extract all reference sequences from a panel's references section.

    Returns list of dicts with keys: gene_name, taxon, accession, sequence.
    """
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


def _sanitize_for_fasta_header(s: str) -> str:
    """Replace characters that could break FASTA header parsing."""
    return s.replace("|", "_").replace(" ", "_").replace("/", "_")


def build_reference_db(panel_data: dict, tmpdir: Path) -> Path | None:
    """Build a BLAST DB from panel reference sequences.

    Returns the DB path prefix, or None if no references or tools unavailable.
    Headers encode gene|taxon|accession for post-BLAST identification.
    """
    refs = extract_references(panel_data)
    if not refs:
        return None

    if not check_blastn_available() or not check_makeblastdb_available():
        print("WARNING: blastn or makeblastdb not available; skipping reference BLAST.")
        return None

    # Write multi-FASTA.
    fasta_path = tmpdir / "references.fasta"
    with open(fasta_path, "w") as f:
        for ref in refs:
            gene = _sanitize_for_fasta_header(ref["gene_name"])
            taxon = _sanitize_for_fasta_header(ref["taxon"])
            acc = _sanitize_for_fasta_header(ref["accession"])
            f.write(f">{gene}|{taxon}|{acc}\n")
            # Wrap sequence at 80 chars.
            seq = ref["sequence"]
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")

    # Run makeblastdb.
    db_prefix = tmpdir / "ref_db"
    result = subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta_path),
            "-dbtype",
            "nucl",
            "-out",
            str(db_prefix),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"WARNING: makeblastdb failed: {result.stderr}")
        return None

    return db_prefix


def blast_against_references(
    sequence: str,
    db_path: Path,
    expected_gene: str,
    expected_taxon: str,
) -> RefBlastResult:
    """BLAST a single amplicon against the reference DB.

    Returns a RefBlastResult indicating identity and on-target status.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as f:
        f.write(f">query\n{sequence}\n")
        query_path = f.name

    try:
        result = subprocess.run(
            [
                "blastn",
                "-db",
                str(db_path),
                "-query",
                query_path,
                "-outfmt",
                "5",  # XML output
                "-evalue",
                "1e-10",
                "-max_target_seqs",
                "5",
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        return RefBlastResult(error="BLAST timed out")
    finally:
        Path(query_path).unlink(missing_ok=True)

    if result.returncode != 0:
        return RefBlastResult(error=f"blastn failed: {result.stderr[:200]}")

    return _parse_blast_xml(result.stdout, expected_gene, expected_taxon)


def _parse_blast_xml(
    xml_text: str, expected_gene: str, expected_taxon: str
) -> RefBlastResult:
    """Parse BLAST XML output and extract the top hit."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        return RefBlastResult(error=f"XML parse error: {e}")

    iterations = root.findall(".//Iteration")
    if not iterations:
        return RefBlastResult()

    hits = iterations[0].findall(".//Hit")
    if not hits:
        return RefBlastResult()

    hit = hits[0]
    hit_def = hit.findtext("Hit_def", "")

    # Parse gene|taxon|accession from header.
    parts = hit_def.split("|")
    matched_gene = parts[0] if len(parts) >= 1 else None
    matched_taxon = parts[1].replace("_", " ") if len(parts) >= 2 else None
    matched_accession = parts[2].replace("_", " ") if len(parts) >= 3 else None

    # Extract alignment stats from the top HSP.
    hsp = hit.find(".//Hsp")
    if hsp is None:
        return RefBlastResult(
            matched_gene=matched_gene,
            matched_taxon=matched_taxon,
            matched_accession=matched_accession,
        )

    identity = int(hsp.findtext("Hsp_identity", "0"))
    align_len = int(hsp.findtext("Hsp_align-len", "1"))
    pct_identity = round(100.0 * identity / align_len, 1) if align_len > 0 else 0.0

    # Determine on-target status: matched gene and taxon match expectations.
    on_target = False
    if matched_gene and matched_taxon:
        gene_match = matched_gene.lower() == expected_gene.lower()
        taxon_match = matched_taxon.lower() == expected_taxon.lower()
        on_target = gene_match and taxon_match

    return RefBlastResult(
        matched_gene=matched_gene,
        matched_taxon=matched_taxon,
        matched_accession=matched_accession,
        pct_identity=pct_identity,
        align_length=align_len,
        on_target=on_target,
    )


def blast_all_products(
    run_results: list,
    db_path: Path | None,
    sample_taxon: str,
    skip_blast: bool = False,
):
    """BLAST all recovered amplicons against the reference DB.

    Mutates run_results in place, adding 'reference_match' to each product.
    Only the first (longest) amplicon per gene is validated.
    """
    if skip_blast:
        print("Skipping BLAST (--no-blast).")
        return

    if db_path is None:
        print("No reference DB available; skipping BLAST.")
        return

    count = 0
    for run in run_results:
        if not run["success"]:
            continue
        for prod in run["genes"]:
            seqs = prod.get("sequences", [])
            if not seqs:
                continue
            gene = prod["gene"]
            count += 1
            result = blast_against_references(
                seqs[0], db_path, expected_gene=gene, expected_taxon=sample_taxon
            )
            prod["reference_match"] = {
                "matched_gene": result.matched_gene,
                "matched_taxon": result.matched_taxon,
                "matched_accession": result.matched_accession,
                "pct_identity": result.pct_identity,
                "align_length": result.align_length,
                "on_target": result.on_target,
            }
            if result.error:
                prod["reference_match"]["error"] = result.error
                print(f"  {gene}: BLAST error: {result.error}")
            elif result.pct_identity is not None:
                target_str = "on-target" if result.on_target else "OFF-TARGET"
                print(
                    f"  {gene}: {result.pct_identity}% identity vs "
                    f"{result.matched_taxon} ({target_str})"
                )
            else:
                print(f"  {gene}: no significant hit in reference DB")

    if count == 0:
        print("No amplicons to BLAST.")

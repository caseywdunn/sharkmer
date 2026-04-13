#!/usr/bin/env python3
"""
Bootstrap gold-standard reference sequences for panel validation.

For each (gene, taxon) combination in a panel's validation samples, searches
NCBI Nucleotide for a reference sequence, extracts the expected amplicon by
in-silico primer matching, and outputs a references: YAML block for pasting
into the panel.

Usage:
    python scripts/bootstrap_references.py panels/cnidaria.yaml
    python scripts/bootstrap_references.py panels/cnidaria.yaml --genes 16S CO1
    python scripts/bootstrap_references.py panels/cnidaria.yaml --write

Requires network access for NCBI queries. Uses E-utilities (esearch + efetch).
"""

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib import error, parse, request

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


def derive_gene_name(primer: dict) -> str:
    """Derive output gene name from structured primer fields (mirrors Rust logic)."""
    gene = primer["gene"]
    region = primer.get("region")
    index = primer.get("index")
    name = f"{gene}-{region}" if region is not None else gene
    if index is not None:
        name = f"{name}_{index}"
    return name


# ---------------------------------------------------------------------------
# IUPAC primer matching
# ---------------------------------------------------------------------------

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
    "S": "[CG]", "W": "[AT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def iupac_to_regex(primer: str) -> str:
    """Convert an IUPAC primer sequence to a regex pattern."""
    return "".join(IUPAC.get(c.upper(), c) for c in primer)


def find_primer(seq: str, primer: str, mismatches: int = 2) -> list:
    """Find all positions where primer matches seq (forward strand).

    Uses exact IUPAC matching first, then falls back to allowing mismatches.
    Returns list of (start, end) tuples.
    """
    seq = seq.upper()
    primer = primer.upper()

    # Try exact IUPAC match first.
    pattern = iupac_to_regex(primer)
    matches = [(m.start(), m.end()) for m in re.finditer(f"(?={pattern})", seq)]
    if matches:
        return matches

    # Fall back to mismatch-tolerant search.
    plen = len(primer)
    results = []
    for i in range(len(seq) - plen + 1):
        subseq = seq[i : i + plen]
        mm = 0
        for a, b in zip(subseq, primer):
            allowed = set(IUPAC.get(b, b).strip("[]"))
            if a not in allowed:
                mm += 1
                if mm > mismatches:
                    break
        if mm <= mismatches:
            results.append((i, i + plen))
    return results


def extract_amplicon(
    seq: str, forward: str, reverse: str, min_length: int, max_length: int,
    mismatches: int = 2,
) -> str | None:
    """Extract amplicon from seq using forward and reverse primers.

    Returns the amplicon sequence (including primer binding sites) or None.
    """
    fwd_hits = find_primer(seq, forward, mismatches)
    rev_rc = revcomp(reverse)
    rev_hits = find_primer(seq, rev_rc, mismatches)

    # Find the best pair within length constraints.
    best = None
    for fs, _ in fwd_hits:
        for _, re_ in rev_hits:
            amp_len = re_ - fs
            if min_length <= amp_len <= max_length:
                if best is None or amp_len < best[1]:
                    best = (seq[fs:re_], amp_len)

    return best[0] if best else None


# ---------------------------------------------------------------------------
# NCBI E-utilities
# ---------------------------------------------------------------------------

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RATE_DELAY = 0.4  # seconds between requests (NCBI asks for <=3/sec)


def ncbi_get(url: str, params: dict) -> str:
    """Make a GET request to NCBI E-utilities."""
    query = parse.urlencode(params)
    full_url = f"{url}?{query}"
    req = request.Request(full_url, headers={"User-Agent": "sharkmer-bootstrap/1.0"})
    try:
        with request.urlopen(req, timeout=30) as resp:
            return resp.read().decode("utf-8")
    except error.HTTPError as e:
        print(f"  NCBI HTTP error {e.code}: {e.reason}")
        return ""
    except Exception as e:
        print(f"  NCBI request failed: {e}")
        return ""


def esearch(term: str, db: str = "nucleotide", retmax: int = 5) -> list:
    """Search NCBI and return list of GI/accession IDs."""
    xml = ncbi_get(ESEARCH_URL, {
        "db": db, "term": term, "retmax": retmax, "retmode": "xml",
    })
    if not xml:
        return []
    try:
        root = ET.fromstring(xml)
        return [id_el.text for id_el in root.findall(".//Id")]
    except ET.ParseError:
        return []


def efetch_fasta(ids: list, db: str = "nucleotide") -> str:
    """Fetch sequences in FASTA format from NCBI."""
    if not ids:
        return ""
    xml = ncbi_get(EFETCH_URL, {
        "db": db, "id": ",".join(ids), "rettype": "fasta", "retmode": "text",
    })
    return xml


def parse_fasta(text: str) -> list:
    """Parse multi-FASTA text into list of (header, sequence) tuples."""
    entries = []
    header = None
    seq_parts = []
    for line in text.strip().split("\n"):
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq_parts)))
    return entries


# ---------------------------------------------------------------------------
# Search strategies per gene type
# ---------------------------------------------------------------------------


def search_queries_for_gene(taxon: str, gene_name: str) -> list:
    """Return a list of NCBI search queries to try for a (taxon, gene) pair,
    ordered from most specific to most general."""
    genus = taxon.split()[0] if taxon else ""
    queries = []

    gene_upper = gene_name.upper()

    if gene_upper == "16S":
        queries.append(f'"{taxon}"[Organism] AND mitochondrion[Title] AND complete genome[Title]')
        queries.append(f'"{taxon}"[Organism] AND 16S ribosomal RNA[Title]')
        queries.append(f'"{genus}"[Organism] AND mitochondrion[Title] AND complete genome[Title]')
        queries.append(f'"{genus}"[Organism] AND 16S ribosomal RNA[Title]')
    elif gene_upper == "CO1" or gene_upper == "COI" or gene_upper == "COX1":
        queries.append(f'"{taxon}"[Organism] AND mitochondrion[Title] AND complete genome[Title]')
        queries.append(f'"{taxon}"[Organism] AND (cytochrome oxidase subunit 1[Title] OR COI[Title] OR CO1[Title])')
        queries.append(f'"{genus}"[Organism] AND mitochondrion[Title] AND complete genome[Title]')
    elif gene_upper == "18S":
        queries.append(f'"{taxon}"[Organism] AND 18S ribosomal RNA[Title] AND complete[Title]')
        queries.append(f'"{taxon}"[Organism] AND 18S ribosomal RNA gene[Title]')
        queries.append(f'"{taxon}"[Organism] AND small subunit ribosomal RNA[Title]')
        queries.append(f'"{genus}"[Organism] AND 18S ribosomal RNA[Title]')
    elif gene_upper == "28S" or gene_upper == "28S-V2":
        queries.append(f'"{taxon}"[Organism] AND 28S ribosomal RNA[Title]')
        queries.append(f'"{taxon}"[Organism] AND large subunit ribosomal RNA[Title]')
        queries.append(f'"{genus}"[Organism] AND 28S ribosomal RNA[Title]')
    elif gene_upper in ("ITS", "ITS-V2"):
        queries.append(f'"{taxon}"[Organism] AND internal transcribed spacer[Title]')
        queries.append(f'"{genus}"[Organism] AND internal transcribed spacer[Title]')
    elif gene_upper == "EF1A":
        queries.append(f'"{taxon}"[Organism] AND elongation factor 1 alpha[Title]')
        queries.append(f'"{genus}"[Organism] AND elongation factor 1[Title]')
    else:
        # Generic fallback.
        queries.append(f'"{taxon}"[Organism] AND {gene_name}[Title]')
        queries.append(f'"{genus}"[Organism] AND {gene_name}[Title]')

    return queries


def find_reference_for_gene(
    taxon: str,
    gene_name: str,
    forward: str,
    reverse: str,
    min_length: int,
    max_length: int,
    mismatches: int = 2,
) -> dict | None:
    """Search NCBI for a reference sequence and extract the amplicon.

    Returns {taxon, accession, sequence} or None.
    """
    queries = search_queries_for_gene(taxon, gene_name)

    for query in queries:
        time.sleep(RATE_DELAY)
        ids = esearch(query, retmax=5)
        if not ids:
            continue

        time.sleep(RATE_DELAY)
        fasta_text = efetch_fasta(ids[:3])  # Try top 3 hits.
        if not fasta_text:
            continue

        entries = parse_fasta(fasta_text)
        for header, seq in entries:
            if len(seq) < min_length:
                continue
            amplicon = extract_amplicon(
                seq, forward, reverse, min_length, max_length, mismatches
            )
            if amplicon:
                # Extract accession from header.
                acc = header.split()[0] if header else "unknown"
                # Clean accession (remove version for display).
                print(f"    FOUND: {acc} ({len(amplicon)} bp) from: {header[:80]}")
                return {
                    "taxon": taxon,
                    "accession": acc,
                    "sequence": amplicon,
                }

    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def bootstrap_panel(panel_path: Path, gene_filter: list | None = None):
    """Bootstrap references for all (gene, taxon) pairs in a panel."""
    with open(panel_path) as f:
        panel = yaml.safe_load(f)

    panel_name = panel.get("name", "unknown")
    primers = panel.get("primers", [])
    validation = panel.get("validation") or {}
    samples = validation.get("samples", [])

    if not samples:
        print(f"No validation samples in {panel_name}.")
        return {}

    # Collect unique taxa from validation samples.
    taxa = []
    for s in samples:
        taxon = s.get("taxon", "")
        if taxon and taxon not in taxa:
            taxa.append(taxon)

    if not taxa:
        print(f"No taxa specified in validation samples for {panel_name}.")
        return {}

    print(f"Panel: {panel_name}")
    print(f"Taxa: {', '.join(taxa)}")
    print(f"Genes: {', '.join(derive_gene_name(p) for p in primers)}")
    print()

    # Build references structure.
    # Keys are derived names (for dedup); values carry structured fields + sequences.
    references: dict[str, dict] = {}

    for primer in primers:
        gene = derive_gene_name(primer)
        if gene_filter and gene not in gene_filter:
            continue

        forward = primer["forward_seq"]
        reverse = primer["reverse_seq"]
        min_len = primer.get("min_length", 100)
        max_len = primer.get("max_length", 10000)
        mismatches = primer.get("mismatches", 2)

        print(f"=== {gene} ({forward[:15]}... / {reverse[:15]}...) ===")

        gene_refs = []
        for taxon in taxa:
            print(f"  Searching for {taxon}...")
            ref = find_reference_for_gene(
                taxon, gene, forward, reverse, min_len, max_len, mismatches
            )
            if ref:
                gene_refs.append(ref)
            else:
                print(f"    not found")

        if gene_refs:
            references[gene] = gene_refs
        print()

    return references


def format_references_yaml(references: dict) -> str:
    """Format references dict as YAML text for the references: block."""
    lines = ["references:"]
    for gene, refs in references.items():
        lines.append(f"  - gene: \"{gene}\"")
        lines.append(f"    sequences:")
        for ref in refs:
            lines.append(f"      - taxon: \"{ref['taxon']}\"")
            lines.append(f"        accession: \"{ref['accession']}\"")
            lines.append(f"        sequence: \"{ref['sequence']}\"")
    return "\n".join(lines)


def write_references_to_panel(panel_path: Path, references: dict):
    """Insert or replace references block in panel YAML using ruamel.yaml."""
    from ruamel.yaml import YAML
    from ruamel.yaml.comments import CommentedMap, CommentedSeq

    ryaml = YAML()
    ryaml.preserve_quotes = True
    ryaml.width = 4096

    with open(panel_path) as f:
        data = ryaml.load(f)

    # Build the references structure.
    refs_list = CommentedSeq()
    for gene, entries in references.items():
        gene_block = CommentedMap()
        gene_block["gene"] = gene
        seqs = CommentedSeq()
        for seq_entry in entries:
            seq_block = CommentedMap()
            seq_block["taxon"] = seq_entry["taxon"]
            seq_block["accession"] = seq_entry["accession"]
            seq_block["sequence"] = seq_entry["sequence"]
            seqs.append(seq_block)
        gene_block["sequences"] = seqs
        refs_list.append(gene_block)

    data["references"] = refs_list

    with open(panel_path, "w") as f:
        ryaml.dump(data, f)
    print(f"References written to {panel_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Bootstrap reference sequences from NCBI for panel validation."
    )
    parser.add_argument("panel", help="Path to the panel YAML file")
    parser.add_argument(
        "--genes", nargs="+",
        help="Only search for these genes",
    )
    parser.add_argument(
        "--write", action="store_true",
        help="Write references directly into the panel YAML file",
    )
    args = parser.parse_args()

    panel_path = Path(args.panel).resolve()
    if not panel_path.exists():
        print(f"Panel not found: {panel_path}")
        sys.exit(1)

    references = bootstrap_panel(panel_path, gene_filter=args.genes)

    if not references:
        print("No references found.")
        sys.exit(0)

    # Summary.
    total = sum(len(v) for v in references.values())
    print(f"\nFound {total} reference(s) across {len(references)} gene(s).")

    if args.write:
        write_references_to_panel(panel_path, references)
    else:
        print("\n" + format_references_yaml(references))
        print(f"\nRerun with --write to insert into {panel_path}")


if __name__ == "__main__":
    main()

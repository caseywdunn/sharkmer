"""Primer binding site analysis using IUPAC codes.

Compares observed amplicon ends against the IUPAC-coded primer sequences to
assess degeneracy utilisation. Extracted from scripts/validate_panel.py.
"""

from . import runner

IUPAC_SETS = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "M": frozenset("AC"),
    "K": frozenset("GT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_SET_TO_IUPAC = {v: k for k, v in IUPAC_SETS.items()}
_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"
)
DEFAULT_TRIM = 15


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def iupac_from_set(bases: frozenset) -> str:
    """Smallest IUPAC code covering exactly `bases`. Falls back to N."""
    return _SET_TO_IUPAC.get(bases, "N")


def _fmt_set(s: frozenset) -> str:
    return "{" + ",".join(sorted(s)) + "}"


def _best_sequence_for_gene(runs: list, gene: str) -> str | None:
    """Get the best (highest-depth) recovered sequence for a gene."""
    for run in runs:
        if not run.get("success"):
            continue
        for prod in run.get("genes", []):
            if prod["gene"] == gene:
                seqs = prod.get("sequences", [])
                if seqs:
                    return seqs[0]
    return None


def analyze_primer_bindings(
    panel_data: dict,
    sample_results: list,
    considered_genes: list,
) -> list:
    """Analyse primer bindings across samples for each gene.

    sample_results is a list of (sample_block, runs) tuples.
    Returns a list of per-gene analysis dicts.
    """
    primers = {p["gene_name"]: p for p in panel_data.get("primers", [])}
    analyses = []

    for gene in considered_genes:
        primer = primers.get(gene)
        if primer is None:
            continue
        forward_full = primer.get("forward_seq", "")
        reverse_full = primer.get("reverse_seq", "")
        trim = min(int(primer.get("trim", DEFAULT_TRIM)), runner.K)
        if not forward_full or not reverse_full:
            continue

        forward_spec = (
            forward_full[-trim:] if len(forward_full) >= trim else forward_full
        )
        reverse_spec = (
            reverse_full[-trim:] if len(reverse_full) >= trim else reverse_full
        )

        per_sample_fwd = []
        per_sample_rev = []
        missing = []

        for sample_block, runs in sample_results:
            accession = sample_block["accession"]
            seq = _best_sequence_for_gene(runs, gene)
            if seq is None or len(seq) < max(len(forward_spec), len(reverse_spec)):
                missing.append(accession)
                continue
            seq = seq.upper()
            fwd_obs = seq[: len(forward_spec)]
            rev_obs = revcomp(seq[-len(reverse_spec) :])
            per_sample_fwd.append({"accession": accession, "observed": fwd_obs})
            per_sample_rev.append({"accession": accession, "observed": rev_obs})

        analyses.append(
            {
                "gene_name": gene,
                "forward_full": forward_full,
                "reverse_full": reverse_full,
                "trim": trim,
                "missing_samples": missing,
                "forward": _analyze_one_primer(forward_spec, per_sample_fwd),
                "reverse": _analyze_one_primer(reverse_spec, per_sample_rev),
            }
        )
    return analyses


def _analyze_one_primer(spec: str, per_sample: list) -> dict:
    """Analyse one primer orientation across samples."""
    result = {
        "spec": spec,
        "per_sample": per_sample,
        "position_analysis": [],
        "verdict_lines": [],
    }
    if not per_sample:
        result["verdict_lines"].append(
            "No samples recovered this gene — cannot analyse."
        )
        return result

    spec_len = len(spec)
    per_position = []
    n_fully_utilised = 0
    n_could_reduce = 0
    n_widening = 0
    n_fixed = 0

    for i in range(spec_len):
        code = spec[i]
        allowed = IUPAC_SETS.get(code, frozenset())
        observed_bases = frozenset(
            s["observed"][i] for s in per_sample if i < len(s["observed"])
        )
        outside = observed_bases - allowed
        inside = observed_bases & allowed
        if len(allowed) == 1:
            status = "fixed"
            if outside:
                status = "mismatch_absorbed"
                n_widening += 1
            else:
                n_fixed += 1
        else:
            if outside:
                status = "mismatch_absorbed"
                n_widening += 1
            elif inside == allowed:
                status = "fully_utilised"
                n_fully_utilised += 1
            else:
                status = "could_reduce"
                n_could_reduce += 1

        per_position.append(
            {
                "pos": i + 1,
                "spec": code,
                "allowed": allowed,
                "observed": observed_bases,
                "status": status,
            }
        )

    result["position_analysis"] = per_position

    # Build human-readable verdict lines.
    for pos in per_position:
        if pos["status"] == "fully_utilised":
            result["verdict_lines"].append(
                f"  pos {pos['pos']:>2}: {pos['spec']} ({_fmt_set(pos['allowed'])}) "
                f"fully utilised — observed {_fmt_set(pos['observed'])}"
            )
        elif pos["status"] == "could_reduce":
            reduced = iupac_from_set(pos["observed"])
            result["verdict_lines"].append(
                f"  pos {pos['pos']:>2}: {pos['spec']} ({_fmt_set(pos['allowed'])}) "
                f"partially utilised — observed {_fmt_set(pos['observed'])}; "
                f"could reduce to {reduced}"
            )
        elif pos["status"] == "mismatch_absorbed":
            widened = iupac_from_set(pos["allowed"] | pos["observed"])
            result["verdict_lines"].append(
                f"  pos {pos['pos']:>2}: {pos['spec']} ({_fmt_set(pos['allowed'])}) "
                f"observed {_fmt_set(pos['observed'])} — off-code base(s) "
                f"absorbed by --mismatches; consider widening to {widened}"
            )

    if n_widening == 0 and n_could_reduce == 0 and n_fully_utilised == 0:
        result["verdict_lines"].append(
            "  all positions fixed (no degeneracy codes) and all observed bases "
            "match — nothing to tune."
        )
    elif n_widening == 0 and n_could_reduce == 0:
        result["verdict_lines"].insert(
            0,
            f"  {n_fully_utilised} degenerate position(s) fully utilised across samples.",
        )
    return result

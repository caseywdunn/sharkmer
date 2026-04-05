#!/usr/bin/env python3
"""
Validate a sharkmer panel against its declared ENA/SRA samples.

Reads the panel's `validation.samples` block, runs sharkmer against each
sample at each declared read depth, BLASTs the recovered amplicons for
identity, compares against the `expected` thresholds (if any), and emits
a markdown report to panels/reports/.

With --write, also updates the panel YAML in place: populates the
`expected` block for each gene with observed values (rounded, with a
small safety margin) and updates `validation.last_validated`.

With --genes, only runs the specified gene(s). Useful when iterating on
a single primer pair without re-running the whole panel.

Reuses infrastructure from benchmarks/ (sharkmer runner, FASTA parser,
local/remote BLAST). Run from the repo root with the sharkmer-bench
conda environment active.

Usage:
    python scripts/validate_panel.py panels/cnidaria.yaml
    python scripts/validate_panel.py panels/cnidaria.yaml --write
    python scripts/validate_panel.py panels/cnidaria.yaml --genes 16S CO1
    python scripts/validate_panel.py panels/cnidaria.yaml --no-blast
"""

import argparse
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

from ruamel.yaml import YAML

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
sys.path.insert(0, str(BENCHMARKS_DIR))

# Reuse benchmarks infrastructure.
from run_benchmark import (  # noqa: E402
    CACHE_DIR,
    DATA_DIR,
    K,
    SHARKMER_BIN,
    THREADS,
    get_sharkmer_version,
    parse_fasta_products,
)
from blast_validate import (  # noqa: E402
    blast_batch_remote,
    blast_sequence_local,
    check_blastn_available,
    find_local_blast_db,
    get_git_email,
    verify_local_blast_db,
)


REPORTS_DIR = REPO_ROOT / "panels" / "reports"
RUNS_DIR = REPO_ROOT / "panels" / "validation_runs"


def clean_sharkmer_version(raw: str) -> str:
    """Extract just the version number from sharkmer's --version output.

    Example: "sharkmer 3.0.0-dev (https://github.com/caseywdunn/sharkmer)"
             -> "3.0.0-dev"
    Falls back to the raw string if the expected format is not recognised.
    """
    # Drop anything in parentheses (repo URL), drop the "sharkmer" prefix.
    s = raw.split("(")[0].strip()
    parts = s.split()
    if len(parts) >= 2 and parts[0].lower() == "sharkmer":
        return parts[1]
    return s or raw

# Default safety margins used when writing observed values back as thresholds.
# See PANELS.md for rationale — the goal is to absorb normal run-to-run
# variation so the panel does not break on small shifts.
IDENTITY_SAFETY_MARGIN = 0.02
MIN_LENGTH_TOLERANCE = 20
LENGTH_TOLERANCE_FRACTION = 0.05  # 5% of observed length


# ---------------------------------------------------------------------------
# Panel loading
# ---------------------------------------------------------------------------


def load_panel(panel_path: Path):
    """Load a panel with ruamel.yaml for round-tripping."""
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.width = 4096
    with open(panel_path) as f:
        data = yaml.load(f)
    return yaml, data


def panel_gene_names(panel_data) -> list:
    return [p["gene_name"] for p in panel_data.get("primers", [])]


def build_filtered_panel(panel_path: Path, gene_filter):
    """If gene_filter is non-empty, write a temp panel YAML containing only
    those genes (plus all other panel metadata). Returns (path_to_use,
    temp_path_or_none). The temp file should be deleted by the caller.
    """
    if not gene_filter:
        return panel_path, None

    yaml, data = load_panel(panel_path)
    all_genes = panel_gene_names(data)
    missing = [g for g in gene_filter if g not in all_genes]
    if missing:
        raise SystemExit(
            f"--genes filter references genes not in panel: {missing}\n"
            f"Available: {all_genes}"
        )
    data["primers"] = [p for p in data["primers"] if p["gene_name"] in gene_filter]

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False, dir=str(RUNS_DIR)
    )
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    yaml.dump(data, tmp)
    tmp.close()
    return Path(tmp.name), Path(tmp.name)


# ---------------------------------------------------------------------------
# Sharkmer execution
# ---------------------------------------------------------------------------


def run_sharkmer_for_sample(
    panel_path_for_run: Path,
    panel_name: str,
    sample: dict,
    max_reads: int,
    run_dir: Path,
):
    """Run sharkmer once for (sample, max_reads). Returns a dict with
    sample_prefix, wall_time, success, genes (list of product dicts)."""
    accession = sample["accession"]
    k_reads = max_reads // 1000
    sample_prefix = f"{panel_name}_{accession}_{k_reads}k"
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(SHARKMER_BIN),
        "-k", str(K),
        "-t", str(THREADS),
        "--max-reads", str(max_reads),
        "-o", str(run_dir) + "/",
        "-s", sample_prefix,
        "--pcr-panel-file", str(panel_path_for_run),
    ]

    # Prefer a local FASTQ if one is present; otherwise stream via --ena.
    local_fq = DATA_DIR / f"{accession}.fastq"
    if local_fq.exists():
        cmd.append(str(local_fq))
        source = "local"
    else:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--ena", accession, "--cache-dir", str(CACHE_DIR)])
        source = "ena"

    print(f"  [{source}] running sharkmer for {accession} @ {k_reads}k reads...")
    start = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - start

    if result.returncode != 0:
        print(f"  ERROR: sharkmer failed ({wall:.1f}s)")
        print(f"  stderr tail: {result.stderr[-500:]}")
        return {
            "sample_prefix": sample_prefix,
            "accession": accession,
            "max_reads": max_reads,
            "wall_time_s": round(wall, 1),
            "success": False,
            "genes": [],
        }

    # Persist the sharkmer log for provenance.
    log_path = run_dir / f"{sample_prefix}.log"
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    products = parse_fasta_products(sample_prefix, run_dir)
    # parse_fasta_products returns gene names with the panel prefix (e.g.
    # "cnidaria_16S"). Strip it so we can match against the panel's primers.
    stripped = []
    prefix_to_strip = f"{panel_name}_"
    for p in products:
        gene = p["gene"]
        if gene.startswith(prefix_to_strip):
            gene = gene[len(prefix_to_strip):]
        stripped.append({**p, "gene": gene})

    print(
        f"  completed in {wall:.1f}s, {len(stripped)} genes amplified"
    )
    return {
        "sample_prefix": sample_prefix,
        "accession": accession,
        "max_reads": max_reads,
        "wall_time_s": round(wall, 1),
        "success": True,
        "genes": stripped,
    }


# ---------------------------------------------------------------------------
# BLAST identity
# ---------------------------------------------------------------------------


def blast_amplicons(run_results: list, skip_blast: bool):
    """Annotate each product in run_results with 'pct_identity' via BLAST.
    Only the longest (first) amplicon per gene is validated.

    Mutates run_results in place.
    """
    if skip_blast:
        print("Skipping BLAST (--no-blast).")
        return

    # Collect all (query_id, sequence, product_ref, label) entries to BLAST.
    entries = []
    for run in run_results:
        if not run["success"]:
            continue
        for prod in run["genes"]:
            seqs = prod.get("sequences", [])
            if not seqs:
                continue
            qid = f"q{len(entries)}"
            label = f"{run['accession']} {prod['gene']} ({len(seqs[0])} bp)"
            entries.append((qid, seqs[0], prod, label))

    if not entries:
        print("No amplicons to BLAST.")
        return

    local_db = None
    if check_blastn_available():
        local_db = find_local_blast_db()
        if local_db and not verify_local_blast_db(local_db):
            local_db = None

    if local_db:
        print(f"Using local BLAST database: {local_db}")
        for i, (_, seq, prod, label) in enumerate(entries, 1):
            print(f"  [{i}/{len(entries)}] {label}")
            try:
                hit = blast_sequence_local(seq, local_db)
            except Exception as e:
                print(f"    BLAST error: {e}")
                prod["blast_error"] = str(e)
                continue
            if hit:
                prod["blast_hit"] = hit
                prod["pct_identity"] = hit["pct_identity"]
                print(
                    f"    hit: {hit['accession']} {hit['pct_identity']}% "
                    f"{hit['hit_def'][:60]}"
                )
            else:
                prod["blast_hit"] = None
                print("    no significant hit")
    else:
        print("Using NCBI remote BLAST (slow; may fail on rate limits).")
        try:
            email = get_git_email()
        except Exception:
            email = "anonymous@example.com"

        BATCH = 50
        for start in range(0, len(entries), BATCH):
            batch = entries[start : start + BATCH]
            print(f"\nRemote batch {start // BATCH + 1}: {len(batch)} sequences")
            try:
                hits = blast_batch_remote(
                    [(qid, seq) for qid, seq, _, _ in batch], email
                )
            except Exception as e:
                print(f"  batch failed: {e}")
                for _, _, prod, _ in batch:
                    prod["blast_error"] = str(e)
                time.sleep(1)
                continue
            for qid, _, prod, label in batch:
                hit = hits.get(qid)
                if hit:
                    prod["blast_hit"] = hit
                    prod["pct_identity"] = hit["pct_identity"]
                    print(f"  {label} -> {hit['accession']} {hit['pct_identity']}%")
                else:
                    prod["blast_hit"] = None
                    print(f"  {label} -> no hit")
            time.sleep(1)


# ---------------------------------------------------------------------------
# Comparison to expected thresholds
# ---------------------------------------------------------------------------


def evaluate_gene(observed_length, observed_identity, expected):
    """Return ('pass'|'fail'|'unvalidated', reason_str)."""
    if expected is None or not expected:
        return "unvalidated", "no expected thresholds declared"

    reasons = []
    ok = True
    min_id = expected.get("min_identity")
    if min_id is not None:
        if observed_identity is None:
            ok = False
            reasons.append("identity unknown (BLAST missing)")
        elif observed_identity < min_id:
            ok = False
            reasons.append(
                f"identity {observed_identity:.3f} < min {min_id:.3f}"
            )

    expected_length = expected.get("length")
    tolerance = expected.get("length_tolerance", MIN_LENGTH_TOLERANCE)
    if expected_length is not None and observed_length is not None:
        if abs(observed_length - expected_length) > tolerance:
            ok = False
            reasons.append(
                f"length {observed_length} outside {expected_length}±{tolerance}"
            )

    return ("pass" if ok else "fail"), "; ".join(reasons) if reasons else "ok"


def suggest_expected(observed_length, observed_identity):
    """Build a suggested expected block from observed values."""
    suggestion = {}
    if observed_identity is not None:
        suggestion["min_identity"] = round(
            max(0.0, (observed_identity / 100.0) - IDENTITY_SAFETY_MARGIN), 3
        )
    if observed_length is not None:
        suggestion["length"] = int(observed_length)
        suggestion["length_tolerance"] = max(
            MIN_LENGTH_TOLERANCE,
            int(round(observed_length * LENGTH_TOLERANCE_FRACTION)),
        )
    return suggestion


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def build_per_sample_summary(sample, runs, panel_data):
    """Return a dict mapping gene_name -> {observed_length, observed_identity, ...}
    by selecting the best (highest read count) successful run per gene.
    """
    # runs is a list of run dicts for this sample, sorted highest max_reads first
    per_gene = {}
    for run in runs:
        if not run["success"]:
            continue
        for prod in run["genes"]:
            gene = prod["gene"]
            if gene in per_gene:
                continue
            seqs = prod.get("sequences", [])
            lengths = prod.get("lengths", [])
            length = lengths[0] if lengths else (len(seqs[0]) if seqs else None)
            per_gene[gene] = {
                "observed_length": length,
                "observed_identity": prod.get("pct_identity"),
                "max_reads_recovered": run["max_reads"],
                "kmer_count_median": prod.get("kmer_count_median"),
                "blast_hit": prod.get("blast_hit"),
                "sequence": seqs[0] if seqs else None,
            }
    return per_gene


# ---------------------------------------------------------------------------
# Primer binding analysis
# ---------------------------------------------------------------------------
#
# sharkmer matches the 3' `trim` bases of each primer (default 15) against
# k-length kmers in the read set, anchored at the start of the kmer. The
# product written to FASTA therefore contains, at its very ends, the genome
# sequence at the primer binding sites:
#
#   - product[0:trim] = forward primer binding, same orientation as user
#   - product[-trim:] = reverse complement of the reverse primer binding
#
# We compare each position against the user's IUPAC-coded primer to decide
# whether degeneracy was fully utilised, could be reduced (only a subset of
# the coded bases actually seen), or was exceeded (an off-code base was
# absorbed by sharkmer's --mismatches tolerance).

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
_COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")
DEFAULT_TRIM = 15


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def iupac_from_set(bases: frozenset) -> str:
    """Smallest IUPAC code covering exactly `bases`. Falls back to N."""
    return _SET_TO_IUPAC.get(bases, "N")


def panel_primer_map(panel_data):
    """Return dict: gene_name -> primer entry dict (forward_seq, reverse_seq, trim)."""
    return {p["gene_name"]: p for p in panel_data.get("primers", [])}


def analyze_primer_bindings(panel_data, sample_results, considered_genes):
    """For each gene in considered_genes, analyse forward and reverse primer
    bindings across all samples that recovered it. Returns a list of dicts
    suitable for rendering, one per gene.

    Each gene analysis has:
      - gene_name, forward_spec_full, reverse_spec_full, trim
      - forward: per_sample list + position_analysis list + verdict
      - reverse: same, already in user-facing orientation
      - missing_samples: list of accessions that did not recover this gene
    """
    primers = panel_primer_map(panel_data)
    analyses = []

    for gene in considered_genes:
        primer = primers.get(gene)
        if primer is None:
            continue
        forward_full = primer.get("forward_seq", "")
        reverse_full = primer.get("reverse_seq", "")
        trim = int(primer.get("trim", DEFAULT_TRIM))
        # Clamp trim to K to match sharkmer's own clamp (trim > k -> k).
        trim = min(trim, K)
        if not forward_full or not reverse_full:
            continue

        forward_spec = forward_full[-trim:] if len(forward_full) >= trim else forward_full
        reverse_spec = reverse_full[-trim:] if len(reverse_full) >= trim else reverse_full

        per_sample_fwd = []
        per_sample_rev = []
        missing = []

        for sample_block, runs in sample_results:
            accession = sample_block["accession"]
            per_gene = build_per_sample_summary(sample_block, runs, panel_data)
            info = per_gene.get(gene)
            if info is None or not info.get("sequence"):
                missing.append(accession)
                continue
            seq = info["sequence"].upper()
            if len(seq) < len(forward_spec) or len(seq) < len(reverse_spec):
                missing.append(accession)
                continue
            fwd_obs = seq[: len(forward_spec)]
            rev_obs = revcomp(seq[-len(reverse_spec):])
            per_sample_fwd.append({"accession": accession, "observed": fwd_obs})
            per_sample_rev.append({"accession": accession, "observed": rev_obs})

        analyses.append({
            "gene_name": gene,
            "forward_full": forward_full,
            "reverse_full": reverse_full,
            "trim": trim,
            "missing_samples": missing,
            "forward": _analyze_one_primer(forward_spec, per_sample_fwd),
            "reverse": _analyze_one_primer(reverse_spec, per_sample_rev),
        })
    return analyses


def _analyze_one_primer(spec, per_sample):
    """spec is the IUPAC-coded primer (trimmed to the match window).
    per_sample is a list of {accession, observed} dicts. Returns a dict with
    spec, per_sample, position_analysis, verdict.
    """
    result = {
        "spec": spec,
        "per_sample": per_sample,
        "position_analysis": [],
        "verdict_lines": [],
    }
    if not per_sample:
        result["verdict_lines"].append("No samples recovered this gene — cannot analyse.")
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
        observed_bases = frozenset(s["observed"][i] for s in per_sample if i < len(s["observed"]))
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

        per_position.append({
            "pos": i + 1,
            "spec": code,
            "allowed": allowed,
            "observed": observed_bases,
            "status": status,
        })

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
        # 'fixed' positions with exact match are not listed; they are the norm.

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


def _fmt_set(s: frozenset) -> str:
    return "{" + ",".join(sorted(s)) + "}"


def format_binding_section(analyses):
    """Render the primer binding analysis for the markdown report.
    Returns a list of lines."""
    lines = []
    lines.append("## Primer binding analysis")
    lines.append("")
    lines.append(
        "For each gene, the first and last `trim` bases of each recovered "
        "amplicon are compared against the user-specified primer sequence "
        "(3'-trimmed to the match window). Reverse primer bindings are shown "
        "reverse-complemented so they appear in the same orientation as the "
        "primer was written in the panel (i.e. facing the forward primer). "
        "Use this section to decide whether a primer's degeneracy should be "
        "reduced (only a subset of coded bases is actually seen) or widened "
        "(an off-code base was absorbed by sharkmer's `--mismatches` "
        "tolerance)."
    )
    lines.append("")

    for a in analyses:
        gene = a["gene_name"]
        trim = a["trim"]
        lines.append(f"### {gene}")
        lines.append("")
        if a["missing_samples"]:
            lines.append(
                f"**Not recovered** in: {', '.join(a['missing_samples'])}. "
                "This may indicate insufficient primer degeneracy, "
                "insufficient read coverage, or the target taxon lacks the "
                "locus."
            )
            lines.append("")
        if not a["forward"]["per_sample"] and not a["reverse"]["per_sample"]:
            lines.append("_(no samples recovered this gene; skipping alignment)_")
            lines.append("")
            continue

        for which, title in (("forward", "Forward primer"), ("reverse", "Reverse primer")):
            info = a[which]
            full = a[f"{which}_full"]
            lines.append(
                f"**{title}** — spec `{full}` (trim={trim}, match window "
                f"`{info['spec']}`)"
            )
            lines.append("")
            lines.append("```")
            lines.append(f"  {'spec':<14}{info['spec']}")
            # Match marker row: | = fixed match, . = degeneracy resolved,
            # x = mismatch absorbed
            marker = []
            for pos in info["position_analysis"]:
                if pos["status"] == "fixed":
                    marker.append("|")
                elif pos["status"] == "fully_utilised" or pos["status"] == "could_reduce":
                    marker.append(".")
                else:
                    marker.append("x")
            lines.append(f"  {'':<14}{''.join(marker)}")
            for s in info["per_sample"]:
                lines.append(f"  {s['accession']:<14}{s['observed']}")
            lines.append("```")
            lines.append("")
            for v in info["verdict_lines"]:
                lines.append(v)
            lines.append("")
    return lines


def format_identity(pct):
    if pct is None:
        return "—"
    return f"{pct:.1f}%"


def format_pass_fail(status):
    return {"pass": "PASS", "fail": "FAIL", "unvalidated": "—"}[status]


def write_report(
    panel_path,
    panel_data,
    sample_results,
    sharkmer_version,
    gene_filter,
    report_path,
):
    lines = []
    panel_name = panel_data.get("name", "unknown")
    panel_version = panel_data.get("version", "unversioned")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines.append(f"# Panel validation: {panel_name}")
    lines.append("")
    lines.append(f"- Panel file: `{panel_path}`")
    lines.append(f"- Panel version: `{panel_version}`")
    lines.append(f"- sharkmer version: `{sharkmer_version}`")
    lines.append(f"- Date: {now}")
    if gene_filter:
        lines.append(f"- Gene filter: {', '.join(gene_filter)}")
    lines.append("")

    declared_genes = panel_gene_names(panel_data)
    if gene_filter:
        considered_genes = [g for g in declared_genes if g in gene_filter]
    else:
        considered_genes = declared_genes

    for sample_block, runs in sample_results:
        accession = sample_block["accession"]
        taxonomy = sample_block.get("taxonomy", "")
        expected_block = sample_block.get("expected", {}) or {}
        lines.append(f"## Sample: {accession}")
        if taxonomy:
            lines.append(f"")
            lines.append(f"Taxonomy: {taxonomy}")
        lines.append("")

        if any(not r["success"] for r in runs):
            failed = [r for r in runs if not r["success"]]
            lines.append(
                f"Failed runs: {', '.join(str(r['max_reads']) for r in failed)}"
            )
            lines.append("")

        per_gene = build_per_sample_summary(sample_block, runs, panel_data)

        lines.append(
            "| Gene | Observed length | Observed identity | Expected | Status | Notes |"
        )
        lines.append("|------|----------------:|------------------:|----------|--------|-------|")
        for gene in considered_genes:
            info = per_gene.get(gene)
            expected = expected_block.get(gene)
            if info is None:
                status = "fail" if expected else "unvalidated"
                note = "not recovered"
                lines.append(
                    f"| {gene} | — | — | "
                    f"{_format_expected(expected)} | "
                    f"{format_pass_fail(status)} | {note} |"
                )
                continue

            obs_len = info["observed_length"]
            obs_id = info["observed_identity"]
            status, reason = evaluate_gene(obs_len, obs_id, expected)
            lines.append(
                f"| {gene} | {obs_len} | {format_identity(obs_id)} | "
                f"{_format_expected(expected)} | "
                f"{format_pass_fail(status)} | {reason} |"
            )
        lines.append("")

    # Primer binding analysis
    analyses = analyze_primer_bindings(panel_data, sample_results, considered_genes)
    if analyses:
        lines.extend(format_binding_section(analyses))

    # Suggested validation block
    lines.append("## Suggested validation block")
    lines.append("")
    lines.append(
        "Observed values from this run, formatted for pasting into the panel. "
        "Identity is set slightly below observed to absorb run-to-run variation."
    )
    lines.append("")
    lines.append("```yaml")
    lines.append("validation:")
    lines.append("  last_validated:")
    lines.append(f"    sharkmer_version: \"{sharkmer_version}\"")
    lines.append(f"    panel_version: \"{panel_version}\"")
    lines.append(f"    date: {datetime.now().strftime('%Y-%m-%d')}")
    lines.append("  samples:")
    for sample_block, runs in sample_results:
        accession = sample_block["accession"]
        taxonomy = sample_block.get("taxonomy", "")
        max_reads_list = sample_block.get("max_reads", [])
        lines.append(f"    - accession: {accession}")
        if taxonomy:
            lines.append(f"      taxonomy: \"{taxonomy}\"")
        lines.append(f"      max_reads: {list(max_reads_list)}")
        per_gene = build_per_sample_summary(sample_block, runs, panel_data)
        recovered_in_filter = [g for g in considered_genes if g in per_gene]
        if recovered_in_filter:
            lines.append("      expected:")
            for gene in recovered_in_filter:
                info = per_gene[gene]
                suggestion = suggest_expected(
                    info["observed_length"], info["observed_identity"]
                )
                lines.append(
                    f"        \"{gene}\": {_format_suggestion_inline(suggestion)}"
                )
    lines.append("```")
    lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"\nReport written to: {report_path}")


def _format_expected(expected):
    if not expected:
        return "—"
    parts = []
    if "min_identity" in expected and expected["min_identity"] is not None:
        parts.append(f"id≥{expected['min_identity']:.3f}")
    if "length" in expected and expected["length"] is not None:
        tol = expected.get("length_tolerance", MIN_LENGTH_TOLERANCE)
        parts.append(f"len {expected['length']}±{tol}")
    return ", ".join(parts) if parts else "—"


def _format_suggestion_inline(suggestion):
    if not suggestion:
        return "{}"
    parts = []
    if "min_identity" in suggestion:
        parts.append(f"min_identity: {suggestion['min_identity']:.3f}")
    if "length" in suggestion:
        parts.append(f"length: {suggestion['length']}")
    if "length_tolerance" in suggestion:
        parts.append(f"length_tolerance: {suggestion['length_tolerance']}")
    return "{ " + ", ".join(parts) + " }"


# ---------------------------------------------------------------------------
# --write: update panel YAML in place
# ---------------------------------------------------------------------------


def apply_write(
    panel_path: Path,
    panel_data,
    yaml_handle,
    sample_results,
    sharkmer_version,
    gene_filter,
):
    """Update validation.last_validated and validation.samples[*].expected
    for the samples that were run. Only touches genes that are in
    gene_filter (if set), and only samples already declared.
    """
    from ruamel.yaml.comments import CommentedMap

    validation = panel_data.get("validation")
    if validation is None:
        validation = CommentedMap()
        panel_data["validation"] = validation

    # last_validated
    last_validated = CommentedMap()
    last_validated["sharkmer_version"] = sharkmer_version
    last_validated["panel_version"] = panel_data.get("version", "unversioned")
    last_validated["date"] = datetime.now().strftime("%Y-%m-%d")
    validation["last_validated"] = last_validated

    samples_yaml = validation.get("samples", [])
    # Index by accession for lookup
    samples_by_acc = {s["accession"]: s for s in samples_yaml}

    panel_name = panel_data.get("name", "unknown")
    changed_any = False
    warnings = []
    for sample_block, runs in sample_results:
        accession = sample_block["accession"]
        target = samples_by_acc.get(accession)
        if target is None:
            continue  # --write refuses to add new samples
        per_gene = build_per_sample_summary(sample_block, runs, panel_data)

        expected = target.get("expected")
        if expected is None:
            expected = CommentedMap()
            target["expected"] = expected

        declared_genes = panel_gene_names(panel_data)
        genes_to_update = gene_filter if gene_filter else declared_genes
        for gene in genes_to_update:
            info = per_gene.get(gene)
            if info is None:
                # Not recovered; do not invent an expected block for it.
                continue
            new_suggestion = suggest_expected(
                info["observed_length"], info["observed_identity"]
            )
            if not new_suggestion:
                continue

            old = expected.get(gene)
            if old is not None:
                old_min_id = old.get("min_identity")
                new_min_id = new_suggestion.get("min_identity")
                if (
                    old_min_id is not None
                    and new_min_id is not None
                    and new_min_id < old_min_id
                ):
                    warnings.append(
                        f"{accession} {gene}: loosening min_identity "
                        f"{old_min_id:.3f} -> {new_min_id:.3f}"
                    )

            new_entry = CommentedMap()
            for k, v in new_suggestion.items():
                new_entry[k] = v
            expected[gene] = new_entry
            changed_any = True

    if warnings:
        print("\nWARNING: --write would loosen existing thresholds:")
        for w in warnings:
            print(f"  - {w}")
        print("  These may indicate regressions — inspect before committing.")

    with open(panel_path, "w") as f:
        yaml_handle.dump(panel_data, f)

    if changed_any:
        print(f"\nPanel updated in place: {panel_path}")
    else:
        print(f"\nPanel updated (last_validated only): {panel_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Validate a sharkmer panel against its declared samples."
    )
    parser.add_argument("panel", help="Path to the panel YAML file")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Rewrite the panel YAML in place, populating expected thresholds "
        "and last_validated from observed values.",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        help="Only run validation for these gene names (use panel-native names, "
        "without the panel prefix). Useful when iterating on a single primer.",
    )
    parser.add_argument(
        "--no-blast",
        action="store_true",
        help="Skip BLAST identity validation. Length-only comparisons still run.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPORTS_DIR,
        help=f"Directory for markdown reports (default: {REPORTS_DIR})",
    )
    args = parser.parse_args()

    panel_path = Path(args.panel).resolve()
    if not panel_path.exists():
        print(f"Panel file not found: {panel_path}")
        sys.exit(1)

    if not SHARKMER_BIN.exists():
        print("Building sharkmer (release)...")
        result = subprocess.run(
            ["cargo", "build", "--release"],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"Build failed: {result.stderr}")
            sys.exit(1)

    yaml_handle, panel_data = load_panel(panel_path)
    panel_name = panel_data.get("name")
    if not panel_name:
        print(f"Panel file missing required 'name' field: {panel_path}")
        sys.exit(1)

    validation = panel_data.get("validation") or {}
    samples = validation.get("samples") or []
    if not samples:
        print(
            f"Panel '{panel_name}' has no validation.samples declared.\n"
            "Add samples to the panel file before running the validator."
        )
        sys.exit(1)

    gene_filter = list(args.genes) if args.genes else None

    # If a gene filter is set, write a filtered temp panel to avoid running
    # primers we don't care about.
    panel_path_for_run, temp_panel_path = build_filtered_panel(
        panel_path, gene_filter
    )

    try:
        sharkmer_version = clean_sharkmer_version(get_sharkmer_version())
        print(f"sharkmer version: {sharkmer_version}")
        print(f"panel: {panel_name} v{panel_data.get('version', 'unversioned')}")
        print(f"samples: {len(samples)}")
        if gene_filter:
            print(f"gene filter: {', '.join(gene_filter)}")
        print()

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = RUNS_DIR / f"{panel_name}_{stamp}"

        # Run sharkmer for each sample × max_reads (descending, like benchmarks).
        sample_results = []
        for sample_block in samples:
            accession = sample_block["accession"]
            max_reads_list = sorted(
                sample_block.get("max_reads", [1_000_000]), reverse=True
            )
            print(f"=== {accession} ({len(max_reads_list)} depths) ===")
            runs = []
            for max_reads in max_reads_list:
                run = run_sharkmer_for_sample(
                    panel_path_for_run,
                    panel_name,
                    sample_block,
                    max_reads,
                    run_dir,
                )
                runs.append(run)
            sample_results.append((sample_block, runs))
            print()

        # BLAST the amplicons (just product 0 per gene).
        all_runs = [r for _, runs in sample_results for r in runs]
        blast_amplicons(all_runs, skip_blast=args.no_blast)

        # Write report.
        panel_version = panel_data.get("version", "unversioned")
        report_name = f"{panel_name}_{panel_version}_{sharkmer_version}_{stamp}.md"
        report_path = args.output_dir / report_name
        write_report(
            panel_path,
            panel_data,
            sample_results,
            sharkmer_version,
            gene_filter,
            report_path,
        )

        if args.write:
            apply_write(
                panel_path,
                panel_data,
                yaml_handle,
                sample_results,
                sharkmer_version,
                gene_filter,
            )
        else:
            print(
                "\nDefault run is read-only. Rerun with --write to update the "
                "panel file in place."
            )
    finally:
        if temp_panel_path is not None and temp_panel_path.exists():
            temp_panel_path.unlink()


if __name__ == "__main__":
    main()

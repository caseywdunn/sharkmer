# Independent pre-outcome diagnostic review

Date: 2026-09-25. Reviewer: independent `diagnostic_review` agent.
Scope: protocol, scanner, synthetic tests, and the frozen six-class/five-candidate
manifest. No biological FASTQ records or new biological outcomes were examined.

Scientific and implementation review: PASS for the source hashes below. The
pre-execution protocol corrections specified below must precede the freeze.

## Independent checks

`independent_synthetic_review.py` implements a full rectangular dynamic-programming
matrix, independently of the production banded distance and traceback code. It
checks distance and optimal-path multiplicity for 450 deterministic generated
sequence pairs at four distance cutoffs (1,800 comparisons). A separate exhaustive
substring enumerator checks every eligible position and length in 100 generated
subjects against the seeded scanner: 70 21-base/two-edit searches and 30
25-base/four-edit searches. Fixtures include repetitive strings, insertions,
deletions, substitutions, invalid internal bases, and cases beyond the cutoff.
These are exhaustive alignments of the generated subjects, not an exhaustive
enumeration of all DNA strings.

The script also uses all six frozen ITS/CO1 classes in nine synthetic variants:
exact, reverse complement, two/four interior substitutions, two flank
substitutions, flank insertion/deletion, and two-base interior insertion/deletion.
All 54 recover the intended unique local class at the independent whole-read
distance. These fixtures need at most 81 ordered pairs per read. Exact and
reverse-complement cases include the complete five-candidate decoy collection;
exact decoy placements, strands, physical coordinates, and declared self-context
intervals are checked independently by literal substring enumeration.

Three further synthetic reads concatenate two copies of each short historical
ITS or retained-CO1 class with a five-base linker. Their lengths are 127, 129, and
125 bases. Each produces 315 ordered pairs and correctly retains two exact best
physical intervals as ambiguous. This independently corroborates the pre-outcome
change of the proposed pair bound from 256 to 4096. The initial results under
256 and subsequent results under 4096 remain separate evidence; no biological
outcome informed the change. The bound is computational, and exceeding it must
fail the complete diagnostic rather than select a subset of pairs.

`independent_cli_review.py` executes the scanner as a subprocess on generated
FASTQ files. Two records with the same header count as two records and one
distinct ID. Successful output carries matching complete staged-file hashes
before/after and the actual argument vector. Reusing the same output is refused
without changing it. Wrong prefix identity, truncated input, malformed quality,
and the 1,001-base read-length overflow all fail with a failure receipt and no
successful output. No temporary synthetic FASTQ is retained as biological data.

## Scientific interpretation

The three-disjoint-seed/two-edit and five-disjoint-seed/four-edit arguments are
appropriate for unit-cost global edit distance. Seed matches alone never count
as evidence. Exact observed bases in both orientations, complete bounded flank
verification, ordered nonoverlapping pairs, full-span alignment, and explicit
class/orientation/physical-interval ties preserve the intended local question.
The unique/ambiguous local status is computed within each declared region's
class set; descriptive decoy placements do not change that status.

Full-candidate comparisons are finite contextual checks. Their absence cannot
establish genome-wide uniqueness. Self-context overlap flags appropriately avoid
equating expected endpoint shifts with independent off-target evidence. Read IDs
and offsets do not establish independent molecules. Q20 applies only to observed
bases, so it cannot validate bases represented as deletions. Representative
tracebacks and their edit-type summaries retain an explicit optimal-path tie
flag. Local alignments cannot resolve full-product phase or justify a repeat-copy
claim, admission decision, or removal of SCC protection.

The unchanged CO1 exact/one-substitution compatibility check remains a required
execution gate. Its expected counts and any new biological result were not used
to modify the independent synthetic controls. Independent held-out biological
evaluation, all-panel comparison, and production threading review remain open.

## Final source and execution clearance

Reviewed scanner SHA-256:
`5c6268e840a626eb980713da8b0a1992030c50bad0e258241dcce703b3f0a75d`.
Reviewed test suite SHA-256:
`0e3f8392659385c970fea136db716dfa52fc9f169976c8bc8738d2a4ce5c6a2e`.
Reviewed hypothesis manifest SHA-256:
`8a61d8a1982fb5e4f2ad2cd40585a150f533066ad09949f9213cf66ebf1bc582`.

The independent oracle's final rerun passed against this unchanged scanner
(`independent_synthetic_results.final.json`). The independent subprocess
checks passed the same scanner (`independent_cli_results.initial.json`). All
29 implementation tests were independently rerun and passed in 5.293 seconds
(`scanner_tests.final.log`). Initial/intermediate oracle receipts identify
their own earlier source hashes and remain preserved.

No remaining scientific or scanner implementation defect was found within this
diagnostic scope. Execution clearance applies to these sources and the following
reviewed pre-execution amendment: replace the draft protocol's 256 ordered-pair
cap with 4,096 and record the synthetic 315-pair rationale above; state that the
unique local status is per region's declared class set. This review explicitly
approves those clarifications before the first biological scan. The coordinator
must reconcile and hash the actual protocol before execution. The source freeze
must include this review and its receipts, exact input provenance, and actual
argv. Running unchanged CO1 compatibility controls and preserving failures
remain execution requirements, not evidence of biological validity supplied
by this review.

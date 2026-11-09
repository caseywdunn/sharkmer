# Independent pre-execution review

Date: 2026-09-26. Scope: the standalone exact ordered-history prototype,
its frozen method, prepared graph/class inputs, and execution wrapper.
No biological FASTQ was opened or scanned by this reviewer. Clearance applies
only to the two frozen graph scans; it does not authorize production admission,
repeat-count claims, full-product reconstruction, or release.

## Verdict

**Pass.** The prototype preserves the read-level information that the current
production edge/adjacent-branch summaries discard. It re-spells only observed
contiguous graph runs, retains repeated edge occurrences and positioned branch
and SCC events, and keeps both orientations and quality layers. Input gaps,
missing graph kmers, and quality boundaries do not become inferred sequence.

The exact literal-window witness screen is suitably narrower than the prior
variation-aware alignment diagnostic. A known indel-containing read may fail
this prototype even when its gapped alignment is credible. The two identical
ITS2 records with different quality endpoints cannot establish two independent
sequence patterns or positions; no terminal deletion or clipped-window rescue
is implemented. The old short-class zeros are known calibration observations,
not held-out targets or biological absence.

The current production representation cannot establish cycle multiplicity or
distant phase: it deduplicates repeated edges and local branch links. Separate
restrictions remain for SCC membership, omitted self-loops, retained collisions,
and links between distant regions. The prototype grants no exception to any of
them, and mapping existing graphs cannot restore their omitted edges.

## Independent algorithm and fixture checks

`independent_history_results.final.json` records **861 successful checks**.
The independent script uses a naive positioned-window mapper and SCCs obtained
by intersections of forward/backward reachability, rather than importing the
producer's run construction or Kosaraju traversal as its expected results.

- Both graph hashes equal the pinned archive members. Every local class string
  and checksum equals the original archived scanner manifests. Independent SCC
  membership agrees with graph annotations; represented anchors are outside
  those SCCs.
- All four archived classes pass positive controls in both graph stages,
  with both read orientations and fixed 150-base reads at starts 0, 1, and 2:
  **48 explicit synthetic fixtures**. This includes the reverse-complement
  representation required by region 1. These graph-derived fixtures test
  implementation sensitivity, not independent biological truth.
- Every checked run equals the independent window oracle and spells exactly
  its oriented read interval. Edge endpoints, repeated-occurrence counts,
  branch-event coordinates, and Q20 span bases are checked directly.
- Reverse-complement records, including IUPAC symbols outside literal windows,
  share complete-sequence groups and normalized footprints. Identical records
  with distinct IDs cannot satisfy the diversity screen. Endpoint/internal Q11
  bases cannot be dropped or hidden by shortened-window matching.
- All 512 three-by-three ID/footprint bipartite graphs agree with brute-force
  maximum matching. A Hall-condition counterexample demonstrates why unrelated
  three-ID/three-position cardinalities are insufficient. The implementation
  also checks that a canonical sequence group cannot have two unique footprints;
  consequently its three distinct matched footprints require three distinct
  sequence groups among those same witness records.
- Repeated anchors, competing low-quality placements, class ties, and
  palindromic strand ties remain ambiguous. Quality cannot select a unique
  occurrence by discarding another unfiltered placement.
- Three/four-copy repeat pools and AA+CC versus AC+CA variant pools have
  identical **summed per-read deduplicated edge and adjacent-branch Counter
  totals**, while their ordered histories differ. These are comparisons on a
  shared frozen synthetic graph, not claims of identical global kmer counts,
  graph construction, or production assembly outputs.
- Unequal repeat mixtures retain both observed classes. Missing self-loops,
  internal indels absent from the graph, and separate fragments of a repeat
  longer than a read remain unrepresented or unresolved. Completely
  unrepresented literal classes survive the input prefilter and cannot appear
  to be absent merely because they lack graph edges.
- Typed prefix limits, wrong hashes, malformed FASTQ, record/edge/read/line
  caps, final-newline-inclusive output caps, and no-clobber publication fail
  without exposing a successful partial JSON result.

The implementation's separate 25-test suite also covers source/input changes,
explicit cycle-node positions, and graph/manifest validation. Source fingerprints
in the final execution review bind the tested implementation and tests.

## Independent execution review

The durable `test_runner_independent.py` suite passes **20 tests**: the 13
previously reviewed execution-wrapper cases are applied to the new wrapper,
plus seven checks for required review entries, failed review status, exact
reviewed-source hashes, invalid review paths, snapshot mutation, freeze mutation,
and durable no-clobber failed-freeze receipts.

The wrapper enforces the reviewed source identities and pinned input hashes,
copies source snapshots, compares input/source/snapshot bookends, inherits the
40-GiB address-space bound, enforces command timeouts, preserves failure receipts,
and refuses to overwrite outputs. Input preparation verifies selected archive
members against the pinned archive and inventory and preserves unchanged class
strings; independent fixture checks verify their downstream graph/class use.

Review identified small pre-freeze consistency fixes: count the final JSON
newline toward its byte cap, retain literal placements with no graph edge,
and complement standard IUPAC symbols when computing full-read canonical groups.
These fixes are covered by the passing final checks. No evidence-driven
threshold or sensitivity change was requested or made.

## Limits and following work

Three matched IDs, sequence groups, and positions are an experimental screen,
not a calibrated confidence model or proof of independent molecules. Finite
class/graph uniqueness is not genome-wide specificity. Separate local witnesses
cannot establish full-product phase, and exact-run failure cannot adjudicate
the earlier gapped ITS alignments. Streaming reads alone does not bound evidence
memory; the frozen explicit caps and wrapper address-space limit fail closed.

Production integration still requires quality-preserving replay and bounded
graph/evidence batching (#145/#116), cause-specific reconstruction evidence,
separate review, and fixed-baseline plus released-version high-copy comparisons.
The 71-to-64 release comparison and issues #148/#153/#129 remain unresolved.
These checks do not cover hostile concurrent filesystem mutation or recovery
after power loss.

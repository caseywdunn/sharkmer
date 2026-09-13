# Gryllus ITS SCC and local read-evidence audit

Date: 2026-09-13. Development baseline: `ae9189a`. Issues #148/#153/#129.
This is a targeted follow-up to the
[withheld-path findings](WITHHELD_PATH_DIAGNOSTICS_RESULTS.md), not a held-out
validation, admission-policy change, or release. Prior marginal read outcomes
are already known. Freeze this method before new graph/bridge outcomes;
freeze derived event manifests before inspecting their read-support counts.

## Question

Why is the exact historical Gryllus ITS_2 978 bp hypothesis withheld, and do
complete, distinguishing local reads support its traversals of the responsible
cyclic components? The candidate SHA-256 is
`23163ca76f97f5d1685f0fcf4bda6ffc798bdf58c63953416833e78a480af68c`.
At k=19 and coverage threshold 4 its five pre-pruning SCC-node markers start
at 93, 94, 276, 277, and 278. These footprints alone do not establish complete
component membership, alternative routes, or a valid bridge interval.

All graph-derived strings are unsupported hypotheses. Same-input raw-read
evidence can corroborate local sequence structure without providing an
orthogonal accessioned reference or establishing the full 978 bp haplotype.
No diagnostic sequence becomes an active panel reference.

## Stage A: actual graph context

1. Build a scratch-only instrumented copy of `ae9189a`; archive the exact
   additive patch, source snapshots, compiler/default-feature build command,
   executable fingerprint, and controls. Do not edit production assembly code.
2. Export the actual complete stored graph immediately before pruning and
   after pruning for `insecta_ITS_2`, threshold 4, k=19 only. Include stable
   oriented 18-mer node identities, start/end flags, oriented 19-mer edges,
   edge counts, and repeat-cause memberships. Node/edge IDs alone are not
   evidence identities. Exports are diagnostic files outside FASTA products.
3. Hard bounds: 100,000 nodes and 500,000 edges per graph, with 256 MiB
   combined serialized bytes for both graphs and their completion manifest. Reject an
   oversize export rather than silently publishing a partial graph. Record
   omitted-self-loop and retained-collision markers explicitly. The stored
   extension graph is not all genome sequence or every count-supported edge;
   high-coverage exclusions and other extension restrictions remain relevant.
4. Use the frozen two-target Gryllus focused panel and exact original first
   1M R1 records: k=19, two threads, chunks=0, unpaired, no read threading.
   Preserve the 40 GiB address-space cap and 1,800-second timeout; verify
   full staged-file and consumed-prefix identities before/after. Later R2
   records are not consumed. This is not a performance benchmark.
5. Require equality with the preceding focused diagnostic-on result for
   all five count fields, ordered full FASTA records, complete gene outcomes,
   thresholds, and diagnostic payloads. Only executable/command/output-path
   and run-timing metadata differ. A mismatch blocks biological interpretation.
6. Independently recompute SCCs from exported edges and verify every candidate
   node/edge and all five markers. Identify complete components, internal
   and boundary edges, and their post-pruning survival. Group all candidate
   occurrences of each component; do not assume the two footprint runs are
   two independent components. Missing connectivity/identities is unresolved.

## Stage B: finite local alternative definitions

Cyclic components allow arbitrarily long walks. **Global alternative
exhaustiveness is not claimed.** The named assay universe is every directed
anchor-to-anchor spelling in the exported graph fitting a single 150 bp read,
including fixed 21 bp flanks on each side. Do not impose the production
node-visit limit as a biological repeat-copy bound.

- Determine candidate-coordinate spans from all traversed nodes belonging
  to each relevant SCC. Extend to literal 21 bp flanks and merge overlapping
  or coupled envelopes before treating regions independently. Lack of flanks,
  repeated anchor placements, or a window longer than 150 bp is explicit.
- Enumerate across the complete exported graph between the selected oriented
  anchors, not merely internal SCC edges: include alternate exits/re-entry
  routes. Analyze pre- and post-pruning graphs separately. Require the full
  fixed flank strings for focal alternatives, not just matching 18-mer anchor
  nodes. Keep walks with changed flank bases in the competing local sequence
  universe: dropping them can manufacture unique substitution-tolerant
  matches. Label them separately, not as absent biology.
- Deterministic traversal limits are 1,000,000 expansion states and 4,096
  distinct spellings per region/stage. Hitting a limit means incomplete
  enumeration. Retain overlength-frontier observations and distinguish
  dead ends from longer possibilities. Counts are graph/search observations,
  not molecule abundance or validated alternatives.
- Reaching the right anchor records a spelling but does not itself terminate
  traversal: longer returning walks must also be considered within the same
  bound. Require noncyclic flanking anchors for qualified local interpretation;
  retain repeated-anchor cases as unresolved rather than choosing a placement.
- Include the historical local spelling explicitly and preserve every
  pre/post route-to-sequence relationship. Deduplicate identical local
  spellings before read tests. Full-length haplotypes sharing a local window
  belong to the same local comparison class, not competing unique templates.
- Read-support comparisons use the union of complete pre/post local spelling
  classes for each region. If enumeration is incomplete, anchors fail, or
  a region's total assayed local alternatives exceed 128, mark that analysis
  unresolved; do not increase caps after observing desired strings/counts.
- Persist exact candidate/window sequences, hashes, graph lineage, coordinates,
  alternative class membership, scope, and completeness/limit status before
  running the read scanner. Independently review them before proceeding.

## Stage C: unchanged read assay

Reuse the reviewed `scripts/audit_read_bridges.py`, without weakening its
quality, orientation, placement, or distinct-read handling. Use one local
comparison manifest per region, with graph-derived window strings as local
hypotheses and events spanning the interior between the 21 bp flanks.

Scan the identical first 1M R1 records twice per qualified region: exact
matching as primary, and separately at most one substitution. Require Q20
throughout the full distinguishing window. The sensitivity result includes
exact matches; it is not an independent replicate. Competing placements or
classes remain ambiguous rather than arbitrarily assigned. Indels are not
silently accepted: enumerated graph-supported length alternatives are tested
as separate literal classes; arbitrary indel-tolerant alignment is outside
this assay and a remaining false-negative source.

Local corroboration requires a matching witness with at least three distinct
read IDs and three distinct projected footprints. These are observed R1
records, not paired fragments or UMI-certified molecules. Preserve complete
matching-record ledgers and all class counts, including zero/ambiguous cases.
Callable but unobserved is not biological absence; uncallable/incomplete is
not a negative test. Recheck the preceding Gryllus retained/lost CO1 structural
control on the same prefix without changing its event definitions.

Synthetic controls cover exact and substituted true alternatives, unequal
mixtures, ambiguous placements, local strings shared across full haplotypes,
repeat-copy changes, cyclic walks requiring more than two node visits,
overlength/cap exhaustion, alternate exits/re-entry, flank/boundary failures,
and distinct-ID/footprint requirements. Separate local positives must not be
concatenated into an unsupported multi-region or full-amplicon phase claim.

## Evidence and decision boundary

Use fresh outputs and durable per-command receipts with actual argv, exit
status, input/source/tool hashes, and complete output hashes. Preserve failed
attempts without overwriting. Source or specification drift invalidates the
associated attestation. Archive method snapshots, graphs, local hypotheses,
scripts/tests, read ledgers, independent checks, and all attempt lineage.

Report whether the actual old ITS local traversals are corroborated,
uncallable, ambiguous, or unobserved within the specified universe; which
alternative strings are supported; and exactly what pruning changed. A
positive result can motivate a narrowly specified read-supported policy
proposal, not automatic removal of SCC gates. Full-length phase, other
samples/targets, fresh high-copy regression comparisons, and user release
review remain separate requirements. Keep #148/#153/#129 open. No release,
tag, master merge, or active-reference change is authorized here.

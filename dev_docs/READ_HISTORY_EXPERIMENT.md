# Ordered read histories: experimental design and evaluation

Date: 2026-09-26. Branch: `experiment/its-read-threading`.
Production comparison baseline: `956e771437a2410a7bcbb1504b033c61dd1e11b6`.
Preceding diagnostic checkpoint: `69104d2`. Issues #148/#145/#116/#153/#129.
Freeze this method, inputs, sources, synthetic controls, and commands before
new biological execution. The preceding [ITS alignment findings](ITS_READ_RECRUITMENT_RESULTS.md)
are already known. This is a calibration experiment, not held-out validation.

## Concrete proposal

Retain complete ordered, positioned read-to-graph runs before reducing them
to edge counts and adjacent branch pairs. A run must preserve repeated edge
occurrences, read orientation, base coordinates, quality continuity, and the
identity of the graph it was mapped against. Spell a local sequence only from
the actual contiguous edge occurrences observed on that read. The first
implementation is a standalone prototype producing diagnostic JSON on this
branch. It does not modify Sharkmer's production assembly or emit amplicons.

This addresses a specific loss of information in the current implementation:

| Current seam | Information missing for repeat reconstruction | Proposed representation |
| --- | --- | --- |
| `src/io.rs`: `ReadRecord`, `reread_sequences` | Qualities and read identifiers; replay collects all sequences | Bounded borrowed read views with original identifiers, qualities, ordinals and mate metadata |
| `src/pcr/threading.rs`: `ReadRun`, `evidence_from_runs` | Runs lose coordinates; aggregate sets deduplicate repeated edges and local branch pairs | Separate positioned runs retaining every edge occurrence and discontinuity |
| `src/pcr/threading.rs`: `map_read` | Longest-run orientation selection does not establish all local placements | Preserve competing placements/orientations and explicit ambiguity |
| `src/pcr/bubble.rs`: `resolve_bubbles` | Local ranking does not establish cycle count or distant phase | Read-spanning ordered histories as separately evaluated evidence |
| `src/pcr/paths.rs`: `search_assembly_paths` | Every SCC, omitted-self-loop or collision cause remains a distinct unresolved obligation | Future cause-specific witness checks; no blanket SCC exception |

The stored Gryllus cycles survive pruning. The historical short path does not
itself revisit a node, so a node-revisit heuristic would not adjudicate the
current loss. Conversely, an edge set or adjacent-pair set cannot distinguish
walks differing only in repeated cycle traversals. Increasing a node-visit
budget is a search change, not evidence for the number of repeat copies.

## First implementation: exact graph-run prototype

1. Read the pinned actual pre- and post-pruning ITS graphs separately. Validate
   node/edge spellings, oriented k-mer uniqueness, connectivity, declared k,
   and graph identity. Use oriented sequence identities rather than unstable
   allocation indices. This experiment adds no missing graph edge.
2. Stream the identical first million Gryllus R1 records, retaining the original
   identifiers, ordinals, sequences, qualities, and byte-prefix hash. Map both
   read orientations to the graph. Reverse quality coordinates together with
   sequence orientation; keep all orientations rather than selecting a longest
   run as a unique assignment.
3. Produce an unfiltered exact layer and a Q20 exact layer. A contiguous run
   requires consecutive read k-mer starts and matching graph endpoints. Split
   at a non-ACGT base, non-graph k-mer, or (for the Q20 layer) any below-Q20
   base in the k-mer. Preserve positions and the reason for each discontinuity.
   Never concatenate runs across a mismatch, indel, invalid base, low-quality
   interval, read boundary, or mate boundary.
4. Retain repeated edge occurrences in order. Re-spelling a run must reproduce
   exactly its oriented raw-read interval; it must not insert, delete, extend,
   or patch observed bases. Record branch events and cycle-related context as
   positions within this run, not merely a set of supported edges or pairs.
5. Evaluate the four unchanged local ITS classes from the preceding SCC audit.
   Their complete fixed 21-base flanks remain part of each literal window.
   Check both class orientations because the stored first-region class is the
   reverse complement of the historical candidate interval. Mark classes not
   representable in the graph explicitly; absent edges do not imply absent
   biology. Verify the selected flanking anchor nodes are outside the relevant
   cyclic component before assigning local experimental eligibility.
6. A candidate local witness requires the entire literal class in a single
   contiguous Q20 run, including both complete flanks. Determine class,
   placement and strand ambiguity from all unfiltered exact class placements
   first; quality filtering cannot select a unique placement by discarding a
   competing low-quality occurrence. Preserve ambiguous observations but do
   not use them as uniquely placed witnesses.

Normalize a witness's full-read footprint to the class orientation:
`[-class_oriented_read_start, read_length - class_oriented_read_start]`.
Preserve the reverse-complement-canonical hash of the complete read sequence
before quality filtering. Identical sequences and reverse-complement copies
at the same class placement therefore cannot manufacture different footprints
through quality-dependent clipping. This exact prototype never clips a class
window or substitutes an end deletion to pass Q20.

The frozen **experimental eligibility screen** requires a matching witness
set with at least three distinct original read IDs and three distinct
normalized full-read footprints, representing at least three canonical
complete-read sequence groups. Merely having three of each in separate
sets is insufficient: require a coupled distinct-ID/footprint assignment.
The same selected records must also have distinct canonical sequence hashes;
do not combine unrelated cardinalities to pass the screen. For a uniquely
placed literal class, identical canonical complete sequences must normalize
to the same footprint; verify this invariant rather than allowing quality
filtering to create a new position.
Report complete-sequence groups and the actual witness records. Duplicate
records are descriptive counts, not additional positioned support. This
screen is not a probability model, proof of independent molecules, a repeat
copy-number guarantee, or an admission rule.

Eligibility is local and relative to the finite class/graph universe. Separate
eligible windows do not establish that they coexist on one haplotype. Preserve
the read identities and run boundaries needed to determine what is actually
linked; do not concatenate marginal support into a full-product claim.

## Frozen real-input evaluation

Use the archived k=19, threshold-4 pre/post graphs from
`v3.2-its-scc-audit-20260913`, whose archive SHA-256 is
`14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3`.
Pre-graph file SHA-256:
`0810b40673c4e4bdc7dcad0573bbe7f3f49e8a0b4542ff5789bd7ce2ec545bf9`;
post-graph file SHA-256:
`b23e89f11d6ed653c7d8113a3f2c07506932049c8ed2397d4083a2f5fd5dfd40`.
Freeze exact graph members, class manifests, and their lineage before scanning.
Any retained full-candidate/self-context annotations are provenance only.
This exact graph prototype does not repeat the preceding approximate decoy
alignment assay and makes no genome-wide uniqueness claim.

The read prefix is 1,000,000 records, 150,000,000 bases, 365,356,898 bytes,
SHA-256 `6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
The existing prefix-only restoration is allowed and must be re-verified.
Do not add R2, coverage, longer reads, new reference sequences, or tuned k.

Compare pre/post graph representability and complete literal witnesses with
the archived exact SCC assay, whose four Q20 class counts are zero. Agreement
checks assay compatibility; it does not prove biological absence. Report the
exact run breaks for the three previously identified ITS spanning records
without changing sensitivity to make their gapped alignments into graph paths.
In particular, the known identical ITS2 indel-bearing reads are outside the
prototype's exact-walk sensitivity. Their failure to supply a complete exact
run is not a newly discovered biological negative.

Report read/run/branch-history counts, repeated edge occurrences, graph or
quality discontinuities, represented/unrepresented classes, all local witness
IDs/placements, ambiguity and limit statuses. Keep unfiltered and Q20 results
separate. Archive every retained raw read and history needed to audit these
summaries. Streaming alone is not a memory guarantee: publish graph, retained
record, per-read history, and serialized-output limits in the source freeze.
Any exceeded bound invalidates a complete result; do not publish a truncated
zero or eligible verdict. Use one CPU, at most 40 GiB address space, and a
1,800-second limit per graph scan. Timing and RSS describe the prototype only.

The frozen implementation caps graph JSON at 64 MiB, each manifest/expected
prefix JSON at 1 MiB, graph nodes/edges at 100,000/200,000, regions/classes at
16/64, and the scanned prefix at 1,000,000 records. Reads are at most 1,000 bases,
FASTQ lines at most 4,096 bytes, and identifiers at most 256 bytes. Retained
history records are capped at 100,000, total edge occurrences across retained
layers/orientations at 10,000,000, literal placements per read at 4,096, and the
serialized report (including its final newline) at 1 GiB. The per-read run count
and event vectors are additionally bounded by read length and the four fixed
orientation/quality combinations. These are failure limits, not truncation
targets or guarantees that every allowed input fits below the address-space
limit.

## Required synthetic controls and acceptance

| Control | Required outcome |
| --- | --- |
| Exact spanning class at three distinct read starts/IDs | Local screen passes with a coupled witness and exact run spelling |
| Reverse-complement counterparts | Equivalent class-normalized footprints and eligibility |
| Repeated cycle visits, including more than two visits | Ordered occurrence vectors preserve the distinct walk lengths even when edge/adjacent-pair sets agree |
| Unequal mixtures of supported repeat lengths | Preserve each observed local walk; do not choose a majority copy count as the only product |
| N, low-quality interval, or non-graph k-mer | Split histories; disconnected fragments never form a witness |
| Internal substitution/indel absent from the graph | Explicit missing exact representation; no silent alignment or sequence repair |
| Low-quality endpoint with an otherwise matching span | No shortened-window or deletion-based Q20 rescue |
| Identical/RC reads with different qualities or IDs | Same sequence group/footprint; duplication cannot pass the screen |
| Reused IDs and uncoupled three-ID/three-footprint sets | Screen fails unless a valid joint matching witness exists |
| Repeated anchors, multiple placements, palindrome or class tie | Preserve competing observations; no unique-witness assignment |
| Two separated variable regions with matched edge/local-pair counts | Equal-sized pools have equal summed per-read deduplicated edge and adjacent-pair counts, while ordered histories preserve which choices each read spans; no unsupported full phase |
| Omitted self-loop / repeat longer than a read | Unrepresented or unresolved, even with abundant marginal edge coverage |
| Shuffled records | Same normalized class/history evidence and screen result; ordinals retain actual input provenance |
| Changed input/graph/source, malformed FASTQ, resource cap | Explicit failed/incomplete result, durable receipt, no successful partial report |

Freeze and independently review synthetic controls and source before real
execution. Independent checks must verify run spelling, repeated occurrences,
orientation/coordinates, quality boundaries, and the coupled witness screen.
Use no-clobber outputs and source/input bookends; preserve failures and all
attempt lineage. Earlier frozen experiments remain immutable.

## Production integration and comparison criteria

The prototype is the first bounded implementation step. Production integration
would need the quality-preserving replay/batching contracts of #145/#116:
bounded read views, target-graph batches, stable graph identities, and an
explicit budget or spill strategy for histories. The current all-read vector
does not satisfy that contract. Keep records and mate links distinct from
independent molecules or physically observed sequence across mate gaps.

A later path-search proposal must identify exactly which repeat restriction
each witness addresses. Evidence for a stored SCC traversal does not restore
an omitted self-loop, clear a collision marker, prove distant phase, or justify
an unrestricted visit-budget increase. Unrepresented sequences and unlinked
regions remain explicit uncertainty. Every proposed admission change needs
positive evidence for the actual restriction removed and separate review.

Before incorporating production changes into `dev`, freeze a comparison matrix:
same inputs/panels/settings/toolchain, experimental feature off and on against
`956e771`, then the required released-v3.1.0 high-copy regression comparison.
Compare all count fields, ordered product identities/sequences, diagnostics,
repeat/mixture correctness, runtime and peak memory. Feature-off compatibility
is required; a recovery increase is not success without sequence/phase evidence.
Protect clean abundant organelle/rRNA controls and document any retained
uncertainty. The existing 71-to-64 high-copy result remains the release issue.

This experiment changes neither admission nor production threading and is not
a new assembly benchmark. #148/#153/#129, held-out registration, all-panel
validation, dependency review, and user release review remain open. Keep the
experiment branch separate; no merge, tag, publication, or reference promotion
is included.

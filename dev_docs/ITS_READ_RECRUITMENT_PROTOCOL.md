# Gryllus ITS local read recruitment and alignment

Drafted 2026-09-25; reviewed and frozen 2026-09-26.
Baseline: `956e771437a2410a7bcbb1504b033c61dd1e11b6`
on `dev`; experiment branch: `experiment/its-read-threading`.
Issues #148/#153/#129 remain open. The user authorized proceeding with the
next diagnostic and evaluating threading separately before incorporation.

This is the next investigation from [the handoff](HANDOFF.md), following the
[SCC/read audit](ITS_SCC_AUDIT_RESULTS.md). Freeze this protocol, the derived
sequence manifest, source, and synthetic validation before examining new
biological results. Existing zero-match results are already known. This is
calibration evidence, not an independent held-out evaluation.

## Question and scope

Do reads from the original Gryllus input reach either flank of the two
withheld ITS junctions? Do reads spanning both flanks show substitutions or
small indels relative to the historical or extra-cycle local hypotheses?
Distinguish sparse observed coverage, low quality, local sequence variation,
and ambiguous placements without inferring absence from a negative result.

The branch first adds a standalone diagnostic and its evidence. Production
read threading, graph construction, admission gates, defaults, and panel
references stay at the baseline during this experiment. A later threading
implementation requires a concrete evidence-supported proposal, separate
commits, and a fresh comparison before incorporation into `dev`.

## Input identity and frozen hypotheses

Use exactly the first 1,000,000 SRR27962769 R1 records from the earlier audit:
150,000,000 bases, 365,356,898 uncompressed FASTQ bytes, SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
R2, additional coverage, or longer reads are separate experiments. The prior
staged file included unused later records; this new scanner may consume a
prefix-only restoration if its complete byte identity matches the recorded
prefix. This does not reproduce the older wrapper's stronger full-file
attestation. Record the new staged-file identity separately.

Verify the preceding evidence archive and every restored member against its
manifest before use. The SCC archive SHA-256 is
`14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3`.
Use its existing local classes without regenerating or selecting alternatives:

| Region | Historical candidate interval | Historical length | Extra-cycle length |
| --- | --- | ---: | ---: |
| ITS 1 | `[72,133)` | 61 | 124 |
| ITS 2 | `[255,317)` | 62 | 147 |

Intervals are zero-based and half-open. Each literal class supplies its first
and last 21 bases as flanks. Preserve all class sequences, hashes, source
members, and memberships. Keep the unchanged retained/lost CO1 assay and its
five competing full candidates as a separate positive/negative control.
Include all these finite local classes and archived full candidate sequences
in the documented comparison context. These are hypotheses, not references
or an exhaustive genome-wide set of competitors.

## Recruitment and local alignments

1. Stream the FASTQ prefix once, checking record structure and Phred+33
   encoding. Analyze each sequence in both orientations and reverse its
   quality string with the sequence. Do not bridge an `N` or another non-ACGT
   base within a tested flank or paired span; a valid independent flank
   elsewhere in the same record is still retained. Retain a record if either
   complete flank has an alignment with at
   most two unit edits (substitutions, inserted bases, or deleted bases).
   Recruitment precedes quality filtering so low-quality cases stay visible.
2. Seed each 21-base flank with its three disjoint 7-base pieces. For every
   exact seed occurrence, verify candidate flank starts within two bases of
   the seed-projected start and observed lengths 19 through 23. Global edit
   distance verifies the full flank. At most two edits leave at least one
   seed intact; the seed is only a filter, never evidence on its own. Keep
   all verified placements, including overlapping and competing placements.
3. Independently record left-only, right-only, both-flank, and no-ordered-pair
   recruitment. Form every ordered, nonoverlapping left/right pair in one
   oriented read. Preserve pairs even if their span does not fit a known
   class. Reads matching opposite flanks on separate records are never joined.
4. Globally align each complete paired span against the frozen local classes
   using unit-cost Levenshtein distance, at most four total edits. Do not
   penalize only gap opens: a gap of three bases costs three edits. Retain
   substitutions, insertion/deletion base counts, coordinates, edit script,
   and optimal-alignment ambiguity. The primary qualified span requires Q20
   at every observed base; report low-quality alignments separately. Deleted
   hypothesis bases have no observed quality score, so Q20 does not validate
   those bases or prove a biological deletion.
5. Within each declared region, compare all eligible classes and paired
   placements. Assign a unique local alignment within that region only when
   the best full-span distance has one class, orientation,
   and physical read interval. Ties across those identities remain ambiguous;
   ordering must not choose a winner. A tied edit traceback is recorded as
   alignment ambiguity even if class/interval agree. Retain other passing
   classes and finite decoy/context matches so a best match is not presented
   as the only biologically possible sequence.

For minimum-distance class/pair contenders, also search the observed full
span against both orientations of each archived full candidate (ITS, 16S,
12S, retained CO1, and lost CO1) at up to four edits. Use five disjoint seeds,
verify starts within four bases and substring lengths within four bases of
the span length, and retain every verified placement. Preserve original
candidate coordinates and orientation. Explicit manifest self-context
intervals label overlap and exact interval agreement; overlapping endpoint
variants at the expected locus are not automatically off-target evidence.
Other placements remain descriptive finite-context ambiguity flags. This
context test does not filter away recruitment or claim genome-wide uniqueness.

Summarize exact, substitution-only, indel-containing, ambiguous, low-quality,
one-flank-only, and paired-but-outside-distance cases separately. Report
distinct read IDs, literal observed spans, and oriented read placements;
records are not independent molecules. A gapped alignment does not establish
the historical repeat copy number. Neither a unique local best alignment nor
two separately supported regions establishes full-length phase. Do not turn
the diagnostic output into an admission decision or count extra-cycle zeros
as equally powered tests of shorter classes.

## Bounds, controls, and freeze

Use at most 100,000 recruited records and 4,096 ordered flank pairs per read.
The initial draft used 256 pairs. Before examining new biological results,
three synthetic repeated-class reads of 125, 127, and 129 bases each generated
315 valid pairs and two exact physical placements. Independent review approved
raising this resource cap to preserve these ambiguity controls; edit and
quality limits are unchanged. Individual frozen classes passed the original cap.

Additional fixed limits are: 1 MiB manifest; 16 regions; 16 classes per region
and 64 total; 32 decoys, each at most 10,000 bases and 100,000 bases total;
64 self-contexts per class; 256-byte identifiers; 1,000-base reads;
4,096-byte FASTQ lines; 4,096 flank placements per read; 4,096 decoy placements
per contender and 16,384 per read; 16 MiB serialized evidence per recruited
record and 1 GiB total serialized output. Hitting a bound fails the assay explicitly; no
successful result may silently represent truncated evidence. The execution
budget is 40 GiB address space and 1,800 seconds per primary scan. Preserve
timeouts and failures and do not increase biological sensitivity limits after
observing results. Synthetic performance checks may precede the freeze.

Synthetic controls must cover exact and reverse-complement reads; more than
one substitution; flank/interior insertions and deletions; both error limits;
left-only/right-only and disconnected support; low quality and invalid bases;
repeat-copy alternatives and unequal mixtures; repeated/palindromic placements;
class ties and shared cross-locus strings; duplicate read IDs; malformed or
changed input; evidence bounds; and refusal to overwrite existing outputs.
Independently compare seeded recruitment with exhaustive flank alignment on
small synthetic inputs and check edit distances against an independent DP.

Rerun the existing unchanged CO1 exact/one-substitution control on the same
prefix. Its archived exact/sensitivity counts are 27/28 retained records and
zero lost records. These are compatibility checks, not independent replicates;
the new four-edit diagnostic need not have those same counts. A mismatch in
the unchanged control blocks interpretation until explained.

Freeze the protocol, manifest, implementation, tests, independent review,
input identities, and actual argv in a durable hash manifest before the first
new biological scan. Use new output directories, no-clobber publication, and
before/after input and source checks. Record exit status, elapsed time, and
output hashes; archive failures and superseded attempts without rewriting
them. Keep raw reads outside Git. If exact inputs cannot be restored, report
the missing prerequisite rather than substituting a different subset.

## Decision and later threading comparison

Report observed recruitment and local sequence relationships within these
fixed limits, including empty or conflicting results. Positive local evidence
may motivate a narrowly scoped threading/reconstruction proposal; it does not
authorize removing SCC protection or establish all seven lost products.
Before integrating any later production change, compare the frozen `dev`
baseline and experiment on identical inputs, panels, settings, and toolchain:
product identities and biological evidence, all count fields, repeat/mixture
controls, timing, and peak memory. Keep diagnostics-on/off compatibility
distinct from parity with released v3.1.0 (still 71 to 64 high-copy products).

Independent held-out registration, all-panel validation, dependency review,
and user release review remain open. No merge, tag, publication, or active
reference promotion is part of this diagnostic.

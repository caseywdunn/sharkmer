# Read-backed audit of withheld high-copy products

Protocol date: 2026-09-13. Baseline source: `a3c4a34` on `dev`.
Issues: #129 and #153. This protocol precedes inspection of read-support
outcomes. No release, primer change, or gate relaxation is authorized by an
absent reference hit or by this protocol alone.
Revision 2, still before real support counts, adds uniformly tiled marginal
diagnostics for cases without serialized graph intervals and specifies the
three-distinct-ID corroboration guard. The original v1 snapshot is preserved.

## Question and frozen inputs

Audit seven historical high-copy losses: Gryllus CO1_1 373 bp and ITS_2
978 bp; Drosophila 12S 475/487/499 bp and 16S_2 590 bp; Heliconius ND1
263 bp. Include retained Gryllus CO1_1 352 bp and retained same-sample
products as controls. Extract candidates from receipt-verified old/new
outputs, deduplicate identical sequences, and pin hashes and sample identity.
They are hypotheses, not biological references.

Use the exact first 1,000,000 FASTQ records consumed by the frozen k=19
benchmark for SRR27962769, SRR31887760, and SRR1057608. The original staged
files contain additional records, but all three consumed prefixes are R1
only and were ingested as unpaired reads. Do not infer pairs from adjacent
records, append R2, or add later reads while calling the result same-input
evidence. Preserve full-file and consumed-prefix hashes, header, ordinal,
quality, and sequence receipts. A future R2/deeper-data audit must be separate.

## Freeze events before counting support

For each candidate, record the difference/repeat or ambiguous graph interval
whose identity must be resolved, its coordinates, and the competing templates
being compared. Candidate differences and graph diagnostics may define these
questions; the reads may not be used to choose favorable event boundaries.
Preserve alternative placements in repetitive alignments rather than selecting
one convenient edit alignment. Graph SCC membership and high edge coverage
are not themselves evidence that the candidate contains a collapsed repeat.

For a separate marginal diagnostic, tile every frozen candidate uniformly:
19 bp intervals starting at 21, 41, 61, and so forth while both 21 bp flanks
fit; also include the last eligible interval if not already present. Each
61 bp window measures local sequence occurrence, **not** resolution of a
graph SCC, structural event, or whole haplotype. Apply this same schedule to
lost and retained candidates before observing counts. Keep diagnostic roles
separate from structural events in every summary and any gate decision.

The primary assay uses a contiguous window containing the complete questioned
interval plus **21 literal A/C/G/T bases on each side**, with at least Q20
(Phred+33) at every assayed base. Coordinates are zero-based, half-open.
If both flanks are unavailable, duplicated, or a window cannot fit in an
observed read, report that limitation rather than silently reducing the assay.
Uniqueness is relative to the frozen candidate/alternative set in both
orientations, not a claim of uniqueness in an unknown whole genome.

Primary observations are exact contiguous window matches. A separately
reported sensitivity assay permits **one substitution, no indels**, while
requiring Q20 and retaining tied placements. This permits a population SNP or
read error to remain visible, without silently accepting an insertion/deletion
that changes the structural question. If that sensitivity assay is not run,
say so explicitly. Other mismatch/gapped methods require a separately recorded
extension of this protocol, with controls before their real-data outcomes.

## Evidence units and decisions

- Preserve record ordinals and read IDs. One observed R1 is at most one vote
  per event; its overlapping k-mers, repeated window occurrences, reverse
  complement, or duplicate representations are not additional votes.
- Report raw matching records, distinct IDs, distinct read sequences, and
  distinct projected start/end footprints separately. Without UMIs these are
  not certified independent molecules. Repeated identical read copies cannot
  manufacture distinct footprint corroboration.
- Retain shared/ambiguous placements and support for competing structures.
  Do not pick one best alignment and label it unique. N/low-quality gaps
  cannot be bridged by concatenating separated matching runs.
- A **corroborated local bridge** requires at least three distinct projected
  read footprints supporting the same complete interval, without a tied
  placement in the assayed alternative set. This is an audit corroboration
  rule, not a validated error rate or full-haplotype truth threshold. At least
  three distinct read IDs must also contribute; repeated identifiers cannot
  manufacture corroboration through altered records.
- A spanning read favoring an incompatible structure is a local conflict for
  the focal structure. Competing or mixed support does not establish absence
  of a lower-abundance allele. Zero/low support, ambiguous placement, and
  missing spanning observations remain **unresolved**, never automatically
  biologically wrong.
- Local bridges do not phase disconnected events into a whole amplicon.
  Nonoverlapping mates, even if separately acquired later, cannot establish
  the sequence or repeat copy count of their unsequenced gap.

Any proposed gate relaxation must identify and positively resolve every
ambiguous interval whose restriction it removes. Marginal edge counts,
branch-entry links, and reference equality do not substitute for that evidence.
Preserve the distinction between assay-level local evidence, candidate-wide
structural support, and unestablished sample haplotype truth.

## Controls and implementation review

Before real support counts, test:

1. True and collapsed homopolymer/tandem-repeat templates with unique flanks:
   spanning reads support the true length, not a collapsed alternative that
   shares strong marginal k-mer coverage. Include a clean short-repeat positive.
2. Recombinant candidates whose individual edges/one-sided junctions occur in
   real templates but whose complete entry-to-exit traversal is not sequenced.
3. Repeats longer than all reads and nonoverlapping mate gaps: remain unresolved.
4. Repeated records, duplicate IDs, reverse complements, palindromes, and
   multiple placements: no inflated support or arbitrary unique assignment.
5. N, low quality, and a single substitution within the bridge: exact and
   sensitivity observations stay separate; interrupted matches are not joined.
6. Unequal-depth mixed alleles: preserve evidence for both and never turn a
   low-support allele into an asserted error.
7. Malformed/truncated FASTQ and changed consumed-prefix or candidate hashes:
   fail closed, rather than publishing a partial successful audit.

Freeze protocol, scanner, control results, candidate/event manifests, and
their hashes before running on the real consumed prefixes. Any later change
gets a new receipt, an explanation, and a full rerun; preserve earlier results.
Sol investigates candidate differences; Terra verifies read provenance; Astra
independently reviews criteria, implementation, and conclusions.

If evidence justifies a targeted fix, add positive and negative regression
tests, review the patch, and rerun the frozen high-copy recovery/runtime/RSS
comparison. Otherwise document the unresolved release gate without changing
assembly policy. Nuclear/metagenomic optimization and #156 remain separate.

## Method boundaries

Read aligners distinguish paired concordance, multiple alignments, and
unpaired fallback; a returned alignment alone is not evidence of unique
fragment phase ([Bowtie 2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)).
BLAST task/word-size choices affect sensitivity, so a previous default-task
no-hit is not a demonstrated sequence absence ([NCBI BLAST+ features](https://www.ncbi.nlm.nih.gov/books/NBK569839/)).
The primary exact-window assay here is deliberately narrower and reports
its limitations rather than converting a no-hit into biological incorrectness.

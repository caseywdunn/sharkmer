# Read-backed audit of seven high-copy losses

Date: 2026-09-13. Issues: #129 and #153. **Release remains on hold.**
This follows the [public-reference replacement assessment](PUBLIC_REFERENCE_REPLACEMENTS.md)
and the [preregistered read-bridge protocol](READ_BACKED_AUDIT_PROTOCOL.md).
No assembly gate, primer, k-mer counter, or default parameter changes here.

## Question and scope

Do the exact reads consumed by the historical comparison corroborate the
structures of seven products emitted by v3.1.0 but withheld by development
code? Neither Sharkmer products nor matching historical panel strings are
independent truth. In particular, the 978 bp Gryllus ITS_2 panel string was a
bootstrap product, not an exact extraction of the 256 bp AK281180.1 record.
The audited candidate strings are explicitly **hypotheses**, not references.

The source comparison uses k=19, two threads, chunks=0, and read threading
off. Its inputs and results are preserved in the
[September 11 evidence](../benchmarks/benchmark_results/v3.2-high-copy-final-20260911/README.md).
This audit reuses exactly the consumed first 1,000,000 records per sample:

| Sample | Read accession | Consumed bases | Read length | Read evidence |
| --- | --- | ---: | --- | --- |
| Gryllus bimaculatus | SRR27962769 | 150,000,000 | 150 bp | R1 only |
| Drosophila melanogaster | SRR31887760 | 101,000,000 | 101 bp | R1 only |
| Heliconius pachinus | SRR1057608 | 89,994,628 | At most 90 bp | R1 only |

Although the staged files later append R2 data, **none of these consumed
prefixes reaches R2**. Adjacent records are not paired fragments. Input
manifests bind original public download URLs, full-file and consumed-prefix
hashes, exact record ordinals, and all six historical result receipts per
sample. No new reads are silently added to the same-input comparison.

## Assays and interpretation

The fixed structural assay requires the complete questioned interval plus
21 A/C/G/T bases on each side, contiguously observed within one read, with
Q20 or better throughout the window. Primary matching is exact. A separate
sensitivity run allows at most one substitution, **not indels**. This lets
a population SNP or sequencing substitution remain visible without erasing
the repeat-length question. Other variation can still produce false negatives;
a missing match is not a biological rejection.

All equivalent repeat placements and competing frozen candidate sequences
remain in the comparison. Uniqueness means unique within this candidate set
and both orientations, not unique in the entire genome. Corroboration requires
three distinct read-ID/footprint pairs, not three overlapping k-mers, copied
records, or reverse-complement representations. These are not UMI-certified
independent molecules. Each event retains its matching read ledger.

The separate marginal assay tiles **all** lost and retained candidates with
61 bp windows on the preregistered schedule. Such windows can reveal local
sequence occurrence when the structural question cannot be localized or fit
within a read. They do not resolve an SCC, repeat copy number, or full-length
haplotype; overlapping windows from different reads must not be concatenated
into a phase claim.

## Results

All six unchanged-assay replication commands exit zero, validate their
outputs, and preserve all 17 frozen source/interpreter hashes. The 273 frozen
events comprise 14 structural/supplementary events and 259 marginal windows.
Exact and sensitivity runs reuse the same three million R1 records; they are
not independent biological replicates. The sensitivity counts **include**
exact matches, rather than adding another independent set of reads.

### Complete structural windows

| Candidate | Full assay window | Exact / at most one substitution | Interpretation |
| --- | --- | --- | --- |
| Gryllus retained CO1_1, 352 bp | 60 bp | 27 / 28 unique records; 15 distinct footprints in each mode | Corroborates the shorter local structure |
| Gryllus lost CO1_1, 373 bp | 81 bp | 0 / 0 records | Callable but not observed; no corroboration of the extra-copy structure |
| Gryllus lost ITS_2, 978 bp | Graph interval unlocalized; whole-candidate declaration has no flanks and would need 1,020 bp | Not callable | Does not test the unknown SCC interval |
| Heliconius lost ND1, 263 bp | Graph interval unlocalized; whole-candidate declaration has no flanks and would need 305 bp | Not callable | Pairwise differences also lack a left flank and exceed the read length |
| Drosophila lost 12S, 475/487/499 bp | First repeat window 91 bp; second 123 or 146 bp | First window 0 / 0 for each; second cannot fit | The first interval is shared across alternatives; complete structures remain unresolved |
| Drosophila lost 16S_2, 590 bp | Repeat windows 133 and 310 bp | Neither fits a 101 bp read | Uncallable, not a negative structural test |

For CO1, the historical 373 bp candidate differs from retained 352 bp by a
21 bp repeat-copy segment with 19 equivalent deletion placements. The lost
candidate's complete placement union is `[310,349)`, assayed as `[289,370)`.
The corresponding retained interval is `[310,328)`, assayed as `[289,349)`.
Thus the comparison spans the complete alignment ambiguity rather than
selecting one convenient deletion boundary. A supplementary 63 bp window
around one difference placement also has zero observations in both modes.

The 27 exact retained observations have 27 distinct IDs and 15 projected
footprints. The sensitivity assay adds one record/ID, not another footprint.
These results **favor the retained shorter local CO1 structure**; they do not
prove absence of a rare longer allele or phase the entire 352 bp amplicon.
They provide no positive justification for rescuing the historical 373 bp
sequence by relaxing a gate.

### Marginal sequence occurrence

This is a different question from structural verification. The table reports
windows with at least one qualifying record and the union of their candidate
coordinate spans, **not** independent verified bases or whole-haplotype
coverage. A sensitivity-window span may itself contain one substitution.

| Lost candidate | Exact observed windows / all windows | Exact span union | Sensitivity observed windows / all windows | Sensitivity span union |
| --- | ---: | ---: | ---: | ---: |
| Gryllus ITS_2, 978 bp | 38 / 47 | 959 bp | 39 / 47 | 978 bp |
| Gryllus CO1_1, 373 bp | 14 / 17 | 321 bp | 14 / 17 | 321 bp |
| Heliconius ND1, 263 bp | 7 / 12 | 204 bp | 9 / 12 | 263 bp |
| Drosophila 12S, 475 bp | 0 / 22 | 0 bp | 0 / 22 | 0 bp |
| Drosophila 12S, 487 bp | 0 / 23 | 0 bp | 0 / 23 | 0 bp |
| Drosophila 12S, 499 bp | 0 / 23 | 0 bp | 0 / 23 | 0 bp |
| Drosophila 16S_2, 590 bp | 0 / 28 | 0 bp | 0 / 28 | 0 bp |

Every tiled window of all five retained controls is observed exactly:
Gryllus CO1_1 352 bp, 12S 423 bp, and 16S_2 534 bp; Heliconius ND1 240 bp;
and Drosophila CO1_1 352 bp. Shared candidate windows are still not unique
evidence distinguishing those alternatives.

**ITS and ND1 have substantial local read evidence.** In particular, the
overlapping sensitivity windows collectively span their full candidate
coordinates. These are not simply sequences with no corresponding read
content. However, the windows need not come from one full-length template,
and the actual withheld graph intervals remain unlocalized. This leaves open
the possibility of conservative false negatives; it does not establish them.
The Drosophila losses lack even these 61 bp observations under either fixed
assay, despite a positive retained CO1 control. That is a lack of corroboration
in this subset, not proof of biological absence: shorter k-mer support,
quality filtering, variation, and incomplete sampling are different criteria.

### Execution and validation

The [evidence archive](../benchmarks/benchmark_results/v3.2-read-backed-audit-20260913/README.md)
preserves manifests, scripts, source snapshots, complete matching-read ledgers,
summary, controls, independent checks, and execution lineage.

The first wrapper changed during its run, correctly failing its source-freeze
attestation despite six completed scanner outputs. Scanner, protocol, and
event definitions remained unchanged. A redundant recovery attempt was
interrupted without overwriting outputs. The final six-command replication
uses a reviewed immutable wrapper, fresh output directory, full-input preflight,
immediate per-job receipts, and before/after hashes including the interpreter.
It is explicitly **nonblind replication after provisional inspection**, and
all six output files are byte-identical to the first attempt. Both attempts
and the failed first attestation remain archived; no timing is retroactively
invented or presented as a Sharkmer performance benchmark.

Astra independently verifies the six reports, all matching-read ledger
observations and alternative placements, three original million-record
prefixes, and all marginal count/span summaries. A separate full-prefix exact
CO1 recount, without importing the scanner, reproduces the 27 retained-window
records and zero lost-window observations. The archived checker and review
receipt preserve this distinction between ledger verification and an
independent complete-prefix recount.

Validation passes: 112 Python regressions (15 new read-bridge tests), 223 Rust
unit tests, 22 integration tests, formatting, Clippy, and exact-source
verification of all 156 active panel references. Synthetic bridge controls
cover true/collapsed repeats, unequal mixtures, recombinant marginal evidence,
long repeats, duplicate IDs/footprints, orientation and placement ties,
quality/N/substitution limits, and malformed or checksum-mismatched input.

## What this does and does not establish

There are three different outcomes to distinguish:

1. **Positive local structural evidence:** a complete discriminating interval
   occurs in multiple distinct read footprints. This supports that local
   structure, not necessarily the entire amplicon or every graph ambiguity.
2. **Callable but not observed:** the structural window could fit, but no
   qualifying read was found in this fixed subset. This is not proof that a
   rare allele is absent, nor that every historical product is wrong.
3. **Uncallable or unlocalized:** reads are too short, flanks unavailable, or
   historical diagnostics do not identify the relevant candidate-local graph
   interval. These are unresolved questions, not negative biological results.

The historical SCC/collision counters describe completed searches at a
threshold. They do **not** serialize withheld path sequences or their marker
coordinates. An aggregate withheld count therefore cannot establish that a
particular old sequence was enumerated in the current search, or locate the
cause on that sequence. Intrinsic sequence repeats and pairwise differences
define useful diagnostic assays but are not interchangeable with those
missing graph-local intervals.

No result here establishes that all seven missing products were incorrect.
Nor does it demonstrate that simply relaxing the SCC gate would recover
known-correct sample haplotypes. Clean-template synthetic reconstruction and
independent public locus context remain distinct from real-sample structural
and haplotype truth.

## Decision and next implementation boundary

**Keep the existing assembly gates unchanged for now; do not release.**
This is a decision not to introduce an unsupported exception, not approval of
the seven-product recovery loss or a conclusion that every gate is optimal.

The next useful implementation is bounded **diagnostic retention of actual
current-run withheld paths and their graph-local ambiguous intervals**, with
explicit truncation and search-budget diagnostics. It must not publish these
paths as supported products. This makes the remaining structural questions
testable instead of attributing an aggregate SCC count to a historical string.

Existing follow-ups cover evidence collection rather than a global bypass:

- #148: read-spanning branch histories and candidate-local repeat evidence.
- #145: bounded, quality-aware replay/spooling over the exact consumed subset.
- #116: stream targeted evidence through bounded batches rather than retaining
  all reads for every gene in RAM.

Current hidden read threading stores sequences without qualities and reduces
contiguous mappings into marginal edge/adjacent-branch support. It does not
provide the complete, quality-qualified interval evidence used here. A safe
future verifier needs bounded current-run candidate requests, one targeted
streaming evidence pass, and final publication only after verification.
Uncached remote input, stdin/spooling, competing placements, and unresolved
long repeats must remain explicit. Nonoverlapping mates do not sequence their
gap. The archived Sol design note describes this boundary; it is not an
implemented gate exception or a release recommendation.

Any subsequent assembly-policy change needs positive and negative regression
controls, Astra review, and a fresh frozen high-copy recovery/runtime/RSS
comparison. This audit alone adds no recovered products and makes **no new
runtime or memory improvement claim**. The earlier 71-to-64 high-copy output
comparison remains historical evidence; no new Sharkmer benchmark ran here.
Held-out registration, all-panel validation, dependency review, the small
human runtime regression, and explicit user release approval remain open.

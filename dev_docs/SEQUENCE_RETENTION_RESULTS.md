# Sequence-first retention assessment

Date: 2026-10-01. Branch: `experiment/sequence-first-retention` from
`b5a6c7e`. This evaluates bounded delivery of completed repeat-uncertain
paths. It does not establish their biological repeat lengths or independent
read support.

Latest follow-up: [the complete target comparison and six new targets](#complete-target-comparison-and-six-new-targets-2026-10-10)
separate sequence availability from reference and local read support. The
original comparison and later diagnostic checkpoints below retain their own
baselines and dates.

## Implementation and Drosophila replay

Inclusive output now budgets previously retained detailed paths and extra
summary-only sequences separately. When an old diagnostic limit would reject
a path, the extra allowance can still emit its bounded sequence with
`assembly_status=uncertain`, `repeat_length=unresolved`,
`read_support=not_evaluated`, and `marker_detail=summary_only`. Its assessment
keeps the exact sequence/hash, threshold, repeat-cause counts, and node-visit
summary, but omits marker identities and coordinates. Complete-detail sources
remain labelled `marker_detail=complete`. The ordinary target FASTA remains a
single file. Standard admission and the graph search are unchanged.

On the fixed SRR31887760 first-million R1 input and full 21-target insect
panel, all 26 previously observed, distinct 12S completions are now retained
at threshold 2. The old inclusive build retained three. The new FASTA has 26
12S uncertain records, of which 23 are new; nine standard sequences across the
panel are identical. The 499 bp and 487 bp historical 12S sequences match by
full-length SHA-256, occur among those 23, and have summary-only marker detail.
No completed 12S path was dropped by sequence retention in this replay.

The 23 summary-only 12S sources omit 9,884 individual detailed marker entries;
the first three complete-detail sources retain 1,274. All 12,747 bases of the
26 completed 12S paths fit the bounded sequence budget. The assessment reports
`marker_detail_truncated=true`, `retention_truncated=false`, and 587
node-visit skips at that threshold. Search and node-visit bounds still prevent
an exhaustive enumeration claim. These 23 additional sequences are hypotheses,
not 23 validated biological products.

## Ten-sample paired comparison

The first candidate failed a non-loss requirement in the Liriodendron pair:
it added 39 uncertain identities but displaced one previously emitted 261 bp
`psbA-trnH` identity when a shared 32-path threshold cap filled. That run was
stopped, and its attempted outputs remain preserved. The corrected candidate
reserves old detailed paths separately and adds a bounded sequence allowance.

All **60 corrected invocations pass** on the exact ten input prefixes and five
panels, with three alternating old-inclusive/corrected pairs per sample. The
old build reproduces **64 standard + 62 uncertain = 126** records. The
corrected build retains every old sequence identity and emits **64 standard +
262 uncertain = 326** records, an additional **200 uncertain identities**.
The uncertain class is 262/326 (80.4%) of delivered records. No standard
identity or input/k-mer count changes. All three repetitions per build are
stable within each sample. Gains are:

| Sample | Target | Additional uncertain identities |
| --- | --- | ---: |
| Liriodendron | psbA-trnH, trnV-atpE | 32 + 32 |
| Drosophila | 12S | 23 |
| Heliconius | ND4 | 32 |
| Gryllus | 18S_1, ITS_2, Yp2 | 17 + 32 + 32 |

The other six samples add none. An independent raw-FASTA/YAML review passes
for all 60 runs. Matching against the previous released-identity inventory
finds all 71 high-copy candidate identities available (64 standard, seven
uncertain) and all 15 deferred identities available as uncertain. The 14
identities absent from the prior inclusive output are the two Drosophila 12S
and 12 Gryllus Yp2 sequences, all exact and labelled uncertain. The Yp2
sequences remain deferred, with no independent biological validation.

Across the ten cells, the sum of median wall times is **434.59 → 442.35 s**
(+1.8%). The largest positive per-sample median change is Drosophila +3.33%;
Gryllus is +2.82%, Liriodendron +1.06%. No sample exceeds the predeclared
5% wall-time or 2% peak-RSS review threshold. Peak RSS is effectively flat
for the large samples; the largest relative change, +0.53%, occurs in the
22 MiB Human cell and is about 0.12 MiB absolute. These are descriptive
two-core, fixed-prefix measurements, not portable performance guarantees.
The earlier v3.1-to-inclusive runtime flags in Rhopilema, Human,
Liriodendron, and Heliconius remain separate unresolved context; this
comparison does not erase them.

Five gene assessments still report sequence-retention truncation. They log
21,637 completed-path **observations** omitted by the extra path cap; these
are not a count of distinct missing sequences. All 12 targets with uncertain
products report search limits reached. The larger delivery therefore remains
bounded and non-exhaustive. It makes more hypotheses available, without
establishing which repeat structures are biological.

## Validation and evidence

The corrected frozen candidate passes release Rust tests on both hash backends (251 unit,
23 integration, and six documentation tests each), Clippy on both backends,
and formatting. The Python suite passes 198 tests with one existing skip;
all 156 active references pass the local audit. A second parser of raw FASTA
and YAML confirms the full comparison independently of its summary. See the
[protocol](SEQUENCE_RETENTION_PROTOCOL.md) and
[evidence bundle](../benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/).

A separate real-input `--no-emit-uncertain --diagnose-withheld-paths` replay
produces identical old/new standard FASTAs, counts, and per-target diagnostic
rows: nine standard records and 50 hidden threshold payloads. Those payloads
remain schema 1 with the original fields and limits.

No release, merge into `dev`, or biological correctness claim follows from
this branch assessment alone.

## Post-comparison reference check (2026-10-02)

The 200 newly emitted uncertain identities were compared with the pinned
121-record public-source catalog and the 61 verified reference regions in the
two relevant panels. This is an exploratory, local reference check, separate
from the frozen paired comparison. The candidate FASTA sequences were bound to
the gain table by full-length SHA-256 before query creation. An independent
review then checked all 200 query hashes against the gain table, all 121 source
sequences against the pinned catalog, saved file receipts, and every raw BLAST
row. See the [review files](../benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/blast-gains-v3/README.md).

With BLASTN 2.17.0+ (`dust=no`, E-value at most `1e-10`, at most 121 targets
and 20 HSPs per target), 168 queries have at least one public-source HSP and
88 have a panel-reference HSP. Across 5,289 public-source HSPs, **none**
combines at least 90% of the query in one HSP with at least 90% identity; none
of 1,103 panel HSPs does either. The five public-source HSPs that cover at
least 90% of a query all map Gryllus 18S hypotheses to a Drosophila source at
about 64% identity. Seventy-eight public-source HSPs reach 90% identity over
shorter segments. No full query is an exact forward or reverse-complement
substring of a catalog source assigned to the same taxon.

These results mean **no new identity has qualifying full-query reference
support in this bounded catalog under this search and criterion**. They do
not show that the hypotheses are false, that no relevant public record exists,
or that any repeat length is correct. Panel targets and public catalog records
have uneven organism and locus coverage; the catalog is a pinned local
snapshot, not a live exhaustive database. Multiple local HSPs may reflect
partial homology or repeats without establishing a coherent full-length
sequence. Full-product read support and independent same-sample truth are
still unevaluated. The earlier panel-only `blast-gains-v1`/`v2` summaries are
preserved as exploratory attempts; their `no_verified_reference` and
`no_significant_hit` statuses are not claims of biological absence.

## Bounded read follow-up (2026-10-09)

The follow-up asks whether nearby alternatives occur together in Gryllus ITS
and whether the recovered Drosophila 12S hypotheses have local read support.
It uses the saved candidates and original, hash-verified first-million R1
prefixes. Full repeat reconstruction remains a separate task.

### Gryllus ITS: some nearby combinations occur on individual reads

The [initial pilot](../benchmarks/benchmark_results/v3.2-its-gain-pilot-20261009/summary.json)
found all eight tested local strings among 25 saved ITS read records. The
[linkage audit](../benchmarks/benchmark_results/v3.2-its-linkage-20261009/run-v1/summary.json)
reuses those reads, preserving both orientations, competing placements,
intervening mismatches, and quality scores. It requires each linked interval
to occur on one read with every base at Q20 or higher. It never joins separate
reads into a linkage witness.

The original windows overlap and omit another variable stretch between them.
The audit therefore records three literal blocks, A/C/B, including that
intervening variation. These are descriptions of strings, not uniquely placed
indels or independently validated alleles.

| Evidence among the 32 gained ITS hypotheses | Hypotheses |
| --- | ---: |
| Complete nearby A/C/B span with both 21-base flanks observed on one Q20 read | 11 |
| Additional complete A/C/B core observed, but right outer flank incomplete | 4 |
| Nearby A/C/B combination unobserved in the selected saved reads | 17 |
| Entire product linked by these reads | 0 |

The fully flanked observations are `TTCGCC / TCTGT / C` (read 32796),
`TTCGCC / TCTGT / T` (671970), and `TTCGCCC / ACGGG / T` (295300).
Their 112–113-base spans match four, four, and three gained hypotheses,
respectively; the third also matches the historical ITS assembly. Reads
273451 and 427221 support `CTCGCCC / ACGGG / T` through the variable core,
with 21 bases of left flank and zero or two bases of right flank. That supports
four additional hypotheses locally, with the narrower flank limitation.
The [per-candidate table](../benchmarks/benchmark_results/v3.2-its-linkage-20261009/run-v1/candidate_disposition.tsv)
retains the individual dispositions.

One further Q20 flanked string, `CTCGCCC / TCGGG / T` in read 784113, is
outside the frozen candidate set. This does not establish a new biological
allele or a production failure: sequencing error, another locus, and bounded
candidate coverage remain possible. Lower-quality spans are recorded separately.

The complete variable extent across the gains is 424–426 bases with flanks,
longer than the 150-base reads. Compatible intervals in this saved ledger
begin at candidate position 722 (zero-based), leaving earlier differences
beginning at 518 unlinked. All 32 full hypotheses remain unresolved.
The 25 reads were selected by the earlier eight-window assay; an unobserved
combination here is not proven absent from the full prefix or sample.
Observed records also do not establish independent molecules or target identity.

Eleven synthetic controls pass, and a fresh-directory replay reproduces the
reports byte for byte. The [independent review](../benchmarks/benchmark_results/v3.2-its-linkage-20261009/independent_review.json)
recounts 231 segment assays, reconstructs 286 candidate placements, verifies
the quality/coordinate evidence, and reproduces the original ITS and CO1
controls. No new FASTQ scan was needed for this linkage analysis.

### Drosophila 12S: longer local strings remain uncorroborated

The [Drosophila assay](../benchmarks/benchmark_results/v3.2-dros12s-read-pilot-20261009/summary.json)
includes all 26 current 12S hypotheses: the three previously retained records
and all 23 gains. Full hashes bind them to the saved assembly output,
including the historical 475-, 487-, and 499-base sequences. Deterministic
61-base tiles, stepped by 21 bases with an end-aligned final tile, cover every
candidate. Nearest-difference windows provide additional local contrasts.
The 568 tile memberships and 52 nearest-window endpoints collapse to 117
distinct strings; shared strings are not independent evidence.

**None of the 117 strings occurs exactly in either orientation in the original
million-read prefix, with or without the Q20 filter.** All 26 hypotheses
therefore lack corroboration under this local assay. The three callable
historical 91-base structural windows also have zero matches. Three longer
windows (123, 123, and 146 bases) are explicitly uncallable with these 101-base
reads, rather than biological negative results. All 16 unchanged CO1 controls
reproduce their prior Q20 counts (9–24 per window) and footprint metrics.
There is no standard 12S positive control in this comparison.

The exact-only implementation uses literal 19-base prefixes to select which
events need the original scanner. Every exact full-window match necessarily
contains its oriented prefix. It retains overlapping seed matches and both
orientations, then delegates observations to the unchanged scanner. Controls
pass 1,835 optimized/direct observation comparisons and 835 supplementary
unfiltered-count comparisons; all 135 saved Q20 matching reads were also
checked directly. The [independent review](../benchmarks/benchmark_results/v3.2-dros12s-read-pilot-20261009/independent_review.json)
verifies source identities, candidate membership, windows, placements,
controls, and receipts. The scan completed in 96.64 seconds with 19.5 MiB peak
RSS inside the existing one-CPU allocation, below 600-second/2-GiB bounds.

Quality filtering does not explain the zeros. This exact-string assay excludes
mismatches and indels and samples only the fixed R1 prefix; it does not prove
biological absence, an incorrect target, or an incorrect full sequence.
Counts from separate windows cannot establish full-product linkage.

### Drosophila graph-word check and production disposition

The absence of longer exact matches prompted a specific invariant check:
does every canonical 19-mer in each emitted sequence occur in the same raw
input? The [bounded recount](../benchmarks/benchmark_results/v3.2-dros12s-read-pilot-20261009/k19-check/run-v1/results/report.json)
finds **every 19-mer start in all 26 hypotheses represented**, with a minimum
raw occurrence count of two in every candidate, matching the assembly's
threshold of two. The 26 hypotheses contain 815 distinct canonical words,
observed 2–392 times. The two historical control sequences are reported
separately. All 1,600 distinct canonical words across candidates and controls
occur; no missing-word coordinates were found. Occurrences include both
orientations and repeated positions within a read, without quality filtering;
they are not counts of independent molecules.

The scan verifies the same prefix hash while reading, completes in 60.60
seconds with 16 MiB peak RSS, and stays within 600-second/2-GiB bounds. Five
synthetic controls include direct randomized counting, reverse complements,
ambiguity resets, repeats, and preventing cross-read words. The
[independent check](../benchmarks/benchmark_results/v3.2-dros12s-read-pilot-20261009/k19-check/independent_review.json)
reviews the source, raw-count table, candidate coordinates, and receipts.

Source inspection also finds consistent sequence reconstruction: uncertain
and standard paths both start with the first `(k-1)`-base node and append the
last base of successive nodes (`src/pcr/paths.rs`, `reconstruct_path_sequence`
and `generate_sequences_from_paths`). Graph extensions require an observed
canonical k-mer meeting the count threshold (`src/pcr/graph.rs`), and primer
matching selects observed k-mers (`src/pcr/primers.rs`). Together with the
recount, this supports construction from observed short graph words; it
does not validate how those words are combined into longer sequences or
prove the absence of all assembly defects.

No specific production failure was demonstrated in this follow-up, so no
assembly, admission, retention, or output change is justified by these
results. Preserve inclusive FASTA delivery and the existing uncertainty,
unresolved repeat length, unevaluated full-product support, and bounded-search
labels. Neither local Gryllus support nor Drosophila 19-mer occurrence promotes
a hypothesis to a standard or biologically validated product.

The scoped linkage and recovered-sequence checks are complete. The next
release work is to resolve the four previously flagged runtime increases and
the outstanding validation/review gates in [PLAN.md](PLAN.md). Full repeat
reconstruction remains separate; further read studies need a concrete
question and explicit bounds. These findings do not authorize integration,
tagging, publication, or closure of the wider biological-validation issues.

## Complete target comparison and six new targets (2026-10-10)

The user-approved recommendations 1–3 freeze the implementation, assemble a
complete comparison table, and assess the six sample–target combinations that
had no output in v3.1.0. The
[implementation freeze](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/implementation_freeze.json)
records branch HEAD `1ea21c9`, 74 source/panel/test file hashes, inclusive
delivery, uncertainty labels, and current search/retention limits. This work
adds diagnostic evidence and documentation; production behavior is unchanged.

### Complete comparison

The [compact comparison CSV](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/final-v1/target_matrix_compact.csv)
contains **all 118 sample–target combinations**, including zero-output rows.
The [full table](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/final-v1/target_matrix.tsv)
adds source receipts, actual assembly status, failure reasons, search and
retention limits, and available reference/read evidence. The
[sequence inventory](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/final-v1/sequence_inventory.tsv)
retains individual sequence hashes and evidence links. Unassessed evidence
is explicit; a search with no qualifying match is distinct from no assessment.
Scoped read assessments are available for 228 of 326 current sequences and
reference assessments for 313; the other 98 and 13, respectively, remain
unassessed. These counts describe assessment coverage, not positive support.

| Matched ten-sample comparison against v3.1.0 | Count |
| --- | ---: |
| Released sequences retained exactly | 86 |
| Released sequences missing now | 0 |
| Current-only sequences | 240 |
| Current total | 326: 64 standard + 262 uncertain |
| Sample–targets with any output, released → current | 68 → 74 |
| Sample–targets still without output | 44 |
| New sequences on six previously empty sample–targets | 154 |
| New sequences on already represented sample–targets | 86 |

All 240 additions are uncertain. Of the 86 historical sequences, 64 remain
standard and 22 are now uncertain. The earlier **200 gains** used the
intermediate inclusive development build as baseline; **240 gains** use
v3.1.0. The older 13-sample, 179-product benchmark has a different scope and
is not the denominator here. The present matrix includes all five panels,
including the four insect single-copy targets excluded from the narrower
high-copy release comparison.

Actual per-target statuses are 64 success, 10 uncertain, and 44 fail, despite
successful CLI execution for every sample. Search bounds were reached on 63
targets; five report sequence-retention truncation and seven report marker
detail truncation. Most failures involve search bounds. Zero output therefore
does not establish biological absence.

### Disposition of the six newly represented targets

| Taxon / sample | Target | New hypotheses | Bounded evidence and disposition |
| --- | --- | ---: | --- |
| Liriodendron tulipifera / SRR25378184 | atpB-rbcL | 1 | Full 812 bp exact match in two annotated Liriodendron plastomes; all 37 local tiles observed at Q20. Strong support for the intended target sequence. |
| Liriodendron tulipifera / SRR25378184 | psbA-trnH | 49 | No qualifying full-query reference match or tested exact local read string. Recovery remains uncorroborated. |
| Liriodendron tulipifera / SRR25378184 | trnV-atpE | 41 | Same lack of corroboration; configured primer identity and length conflict resolved below. These records do not demonstrate recovery of the intended interval. |
| Acer monspessulanum / ERR14009273 | trnL-F | 1 | No qualifying full-query reference match or tested exact local read string. Recovery remains uncorroborated. |
| Heliconius pachinus / SRR1057608 | ND4 | 44 | Partial hits to other loci/taxa, including ND5; no qualifying ND4 reference support or tested exact local read string. Recovery remains uncorroborated. |
| Gryllus bimaculatus / SRR27962769 | 18S_1 | 18 | No qualifying full-query reference match or tested exact local read string. Recovery remains uncorroborated. |

The [reference assessment](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/reference/assessment-v2/report.json)
reuses 113 existing sequence searches and adds 41 bounded searches against
the same pinned 121 public records and 61 relevant panel references. Only
the Liriodendron atpB-rbcL hypothesis has qualifying full-query support:
100% identity over all 812 bases in NC_008326.1 and MN990625.1, with endpoints
in annotated atpB and rbcL. A separate panel-only metric uses a proxy taxon
and does not erase this direct annotated genomic evidence. The other 153
have no single public or expected-target panel HSP meeting both 90% query
coverage and 90% identity. The preserved `assessment-v1` is a draft;
`assessment-v2` is authoritative.

The [read assessment](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/reads/results-v1/summary.json)
tests 61-base tiles across all 154 hypotheses, plus complete-flank local
contrasts between nearest alternatives. All 37 tiles of the 812-base
atpB-rbcL sequence occur in the original Liriodendron reads at Q20 or higher.
For the other 153 hypotheses, **no tested exact local string occurs even
without quality filtering**. Every tile of all eight predefined unchanged
standard positive controls is observed at Q20. Four scans of the original,
hash-verified million-read prefixes take 18.97–37.73 seconds each, 105.50
seconds total, with at most 22.3 MiB peak RSS. Each stays within its one-CPU,
600-second, 2-GiB bound. A factored seed lookup replaces the slower preflight
implementation without changing exact-match semantics.

Counts cover each complete prefix in both orientations. Saved example reads
are selected and capped; they are not the basis for the zero-count results.
Incomplete-flank windows are explicitly uncallable. Overlapping tiles provide
local coverage, not full-product phase or evidence of independent molecules.
The fixed input prefixes, exact-only matching, and uneven reference catalog mean
that an uncorroborated hypothesis is not proven false or absent. Even the
strong atpB-rbcL result retains its production uncertainty classification.

Input-description correction (2026-10-11): the frozen Liriodendron million-record
input comprises 656,431 R1 records followed by 343,569 R2 records, not R1 alone.
Its SHA-256 remains `2a44558f81127ac7aae0f8b9789d3ac5d589986e172c5ba27784a5d84f287bd3`.
The saved manifests' R1 shorthand is inaccurate for this input; the analyzed
bytes, order, exact-match counts, and quality checks are unchanged. Reads were
analyzed individually, with no paired-fragment linkage claim. Subsequent
plant comparisons preserve this exact ordered input.

### trnV-atpE interpretation resolved; panel correction remains separate

The [annotation review](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/reference/annotation_review.json)
maps the exact configured primers and the effective final 15 bases used with
k=19, allowing up to two substitutions per primer. An exhaustive check on both orientations
of three circular plastomes finds exactly one effective oriented pair per
reference: **1,359 bp** in both Liriodendron records and **1,379 bp** in the
Acer opalus proxy. There is no alternative pair within the configured
400–1,000 bp range on these references. Full-primer spans are 1,376/1,396 bp.
The forward endpoint lies in trnV-UAC and the reverse endpoint in atpB;
atpE lies inside the interval. The exact reverse primer is also listed as
the atpB primer S1022 in
[Tank and Donoghue (2010), Appendix 2, p. 441](https://donoghuelab.yale.edu/sites/default/files/190_tank_systbot10.pdf).

This resolves the interpretation of the current configured assay: its label
and length range are inconsistent with the mapped reference interval. It
does not establish that every sampled individual must share that structure.
The 41 short Liriodendron hypotheses have no corroborating match to the mapped
locus. A focused correction under #156 must decide the intended assay,
update its definition and references together, and evaluate the affected
plant cases. This diagnostic does not silently rename the assay, enlarge
production limits, or promote an unreviewed reference.

### Verification and next work

Independent reviews cover the
[complete base matrix](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/review/matrix_review.json),
[read results](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/review/read_assay_review.json),
and [reference/annotation results](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/review/reference_review.json).
They reproduce sequence identities, frozen inputs, local exact/Q20 evidence,
full-reference placements, and exhaustive primer mapping. Source and result
receipts preserve the distinction between current evidence and earlier drafts.
The [final evidence-join review](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/review/final_matrix_review.json)
also passes, and the [freeze verification](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/implementation_freeze_verification.json)
confirms all 74 source hashes and an empty production diff from HEAD.

Recommendations 1–3 establish a stable baseline and an evidence-backed target
inventory. Keep inclusive, bounded delivery and the current uncertainty
labels. Prioritize the concrete plant assay correction, then investigate the
four earlier runtime flags and finish the release gates in [PLAN.md](PLAN.md).
Further biological reconstruction should address a named unresolved case
with a fixed computation budget. The 240 extra records are not evidence of
240 validated new products, and this work alone does not close release review.

# Sequence-first retention assessment

Date: 2026-10-01. Branch: `experiment/sequence-first-retention` from
`b5a6c7e`. This evaluates bounded delivery of completed repeat-uncertain
paths. It does not establish their biological repeat lengths or independent
read support.

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

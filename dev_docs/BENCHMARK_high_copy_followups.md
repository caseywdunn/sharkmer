# High-copy release follow-ups

## Scope frozen 2026-09-11

Protect abundant organelle/ribosomal recovery and runtime/RAM. Metagenomic
diversity and single-copy nuclear recovery are deferred, not additional v3.2
release blockers. Target annotations select candidates for review; they do
not establish measured abundance, read support, or exact repeat-copy truth.
No release is authorized by these measurements.

## #154: repeat-rejected path quota

Sol implementation, independent Astra review, implementation commit
`7b057dcde7594fd64cc0f23903f49c3fa1c23974`. Rejected complete paths no longer
consume the eligible-product quota. DFS work, node visits, and maximum length
remain bounded; unresolved repeat-copy claims remain withheld.

An isolated pre-fix regression enumerated the first 20 of 32 high-ranked
repeat-tainted routes, filtered all 20, and failed to return the later clean
route. The repaired test withholds 32 routes and returns the clean one.
Both hash backends pass 216 unit and 22 integration tests, formatting, and
Clippy with warnings denied. Existing short/long-repeat and pinned high-copy
fixture controls pass.

### Supplemental paired measurement

Compare pristine `4753620` development code with pristine `7b057dc`, using
the same locked release/default-ahash build configuration. This is **not**
a comparison against released v3.1.0. The pre-fix binary SHA-256 is
`387af7374e567828ed334ac488a032d79933f4f34a1da3823edc6c3ea7275b82`;
the #154 binary SHA-256 is
`bcbac0dd58940bbb381970abfd35a06dd9ab2b0f7933674e3fe242e02673a367`.

Three affected insect samples, one alternating paired run each, unchanged
whole panels and frozen one-million-record prefixes. Settings: k=19, two
physical cores, chunks=0, threading off, local uncompressed FASTQ, prefix
prewarming outside timing, 1,800-second timeout, 40 GiB address-space ceiling.
GNU time measures whole-command wall time and process peak RSS. BLAST runs
only after every timed invocation. This does not test network/gzip throughput
or establish statistical equivalence, held-out sensitivity, or laptop limits.

| Sample | Pre-fix wall s | #154 wall s | Pre-fix RSS bytes | #154 RSS bytes | Products |
| --- | ---: | ---: | ---: | ---: | ---: |
| Drosophila / SRR31887760 | 60.22 | 59.92 | 2,287,484,928 | 2,287,312,896 | 9 -> 9 |
| Heliconius / SRR1057608 | 44.18 | 44.23 | 2,287,435,776 | 2,287,296,512 | 6 -> 6 |
| Gryllus / SRR27962769 | 98.21 | 98.41 | 4,569,202,688 | 4,569,210,880 | 11 -> 11 |

All six invocations and per-product evaluations completed. Aggregate
read/base/k-mer occurrence counts, all 26 product sequences, and classifier
outcomes agree. Wall-time sum is 202.61 -> 202.56 seconds; differences are
small and do not demonstrate a speedup. No meaningful RSS change is observed.

The existing release-baseline Gryllus ITS_2 loss is **not repaired by #154**.
The regression proves the quota-ordering bug is fixed, not that every withheld
real-data product has a clean recoverable route. Repeat-marker scope and
threshold/search diagnostics remain #155; direct released-v3.1 comparisons
remain necessary for the release gate in #153.

### Evidence

The [archived evidence](../benchmarks/benchmark_results/v3.2-high-copy-issue154-20260911)
contains frozen protocol/build/validator/input receipts, raw logs, FASTA and
stats, all per-product classifications, paired comparisons, and checksums.
Inputs and binaries are not embedded; their exact identities and clean-build
source commits are recorded. Temporary absolute paths in receipts describe
the measured environment, not portable installation paths.

Review criteria flag every high-copy sequence loss or measured resource
increase; >5% per-cell median runtime or >2% median RSS increases receive
priority investigation. These are review triggers, not permitted regression
allowances or automatic release approval.

## #155: precise marker scope and threshold diagnostics

Sol implementation and independent Astra review are complete. Retained
high-coverage collisions mark the traversed edge, not both endpoints. Omitted
self-loops and nodes in pre-pruning cyclic SCCs remain conservative markers.
Pruning only removes graph nodes/edges, so surviving stable edge identifiers
are not reassigned before search. Tests cover actual collision-marker creation
for new and existing edges, clean acyclic shared endpoints, and repeat controls.

Stats now retain a bounded record for every attempted coverage threshold,
including successful attempts. Records distinguish connectivity, node/DFS/path
limits, candidate-local repeat causes, and unperformed SCC/path evaluation.
Cause-specific withheld-path counts can overlap; generated products are before
final cross-path deduplication. No graph identities or sequences are serialized
into diagnostics, and graph-level markers alone do not establish that every
candidate is ambiguous or the target is absent.

Both hash backends pass 221 unit and 22 integration tests, formatting and
Clippy; 39 Python validation regressions pass. Implementation commit:
`0c38d6a051082cc2a1eb961b4206837f82948fda`; clean-binary SHA-256:
`469065d8416d38f25993d537ee2956d88fdd2ba24a9461e356807c5bf4957a81`.

### Independent controls

All 13 invocations pass their applicable checks: five descriptive v3.1 probes,
five candidate synthetic expectations, and three candidate exact-fixture
checks. The 180 bp threshold target and A18 control remain exact. A19/A40
produce no product and report positive path-local self-loop withholding;
AC40 produces no product and reports positive collision-edge withholding.
This guards against simply exempting every acyclic collision.

ERR571460 has exactly 13,197,385 k-mer occurrences at k=19. At k=31,
both with and without read threading, it preserves exactly the pinned 1,783 bp
18S and 430 bp 28S_2 hashes. Current-run manifests validate. These are known
regression controls, not held-out biological sensitivity measurements.

### Fresh released-v3.1 comparison

**Release remains on hold: the high-copy recovery gate does not pass.**
The final matrix contains 60 fresh invocations: three alternating pairs for
each of ten historical samples across five unchanged non-bacterial panels.
Requested depth is 1M records; the human source exhausts at **108,518 actual
records**. The other nine inputs contain 1M records each. Algorithm, resource,
input-prefix, and timing boundaries match the supplemental protocol above.
The comparator is the attested pristine released-v3.1 binary, SHA-256
`4cb93c72c830e052bf9bf03638eec2a060a699c5741ea0630fd9225ecc5b5e18`.

All 60 CLI invocations and classifications complete. All 30 pairs have exact
aggregate read/base/k-mer occurrence parity. Each version's sequences and
classifier statuses are stable across all three repetitions. Seven non-insect
samples retain exactly the same sequences. Across the ten samples, the
candidate's 64 sequences are also unchanged from pre-fix development code;
#154/#155 fix the demonstrated mechanisms but do not recover additional
real-data products in this calibration.

| Sample | v3.1 median wall s | Candidate median wall s |
| --- | ---: | ---: |
| Xenia / SRR9278435 | 52.19 | 50.96 |
| Agalma / SRR25099394 | 48.39 | 47.03 |
| Rhopilema / SRR8617500 | 53.73 | 52.30 |
| Human / SRR17535371 | 0.39 | 0.44 |
| Nomeus / SRR22396603 | 46.10 | 44.54 |
| Liriodendron / SRR25378184 | 60.96 | 60.37 |
| Acer / ERR14009273 | 31.00 | 30.95 |
| Drosophila / SRR31887760 | 59.74 | 59.46 |
| Heliconius / SRR1057608 | 43.95 | 43.83 |
| Gryllus / SRR27962769 | 98.27 | 95.90 |

Sum of per-cell median wall times: **494.72 -> 485.78 seconds (-1.81%)**.
Nine cells have lower medians; the shortest human input is **0.05 seconds
slower (+12.82%)**, exceeding the percentage review trigger despite its small
absolute change. Do not call this universally regression-free or assign the
timing difference to a particular optimization without profiling. There is
no new counting backend in this release.

Process RSS is effectively unchanged, not improved: median increases are
12–156 KiB on these cells, including 88 KiB for the human sample. The largest
median is about 4.255 GiB. Full per-pair RSS, throughput, and available stage
metrics are in the evidence. Three same-host repeats remain descriptive;
deeper 2M/4M/8M cells, gzip/network ingestion, all nine panels, and a physical
16 GB laptop were not revalidated in this scoped follow-up.

### Remaining high-copy losses

Metadata-selected high-copy products are **71 -> 64**, with no gained
sequences. The seven missing sequences are not seven independently proven
biological false negatives, but they cannot all be dismissed as errors.

| Sample / gene | Lost lengths bp | Baseline median/min k-mer support | Frozen reference classifier | Candidate still recovers this gene? |
| --- | --- | --- | --- | --- |
| Drosophila / 12S | 475, 487, 499 | 6/2 each | No significant hit | No |
| Drosophila / 16S_2 | 590 | 5/2 | No significant hit | No |
| Heliconius / ND1 | 263 | 17/9 | No significant hit | Yes, another sequence |
| Gryllus / CO1_1 | 373 | 52/29 | Insufficient alignment | Yes, another sequence |
| Gryllus / ITS_2 | 978 | 17/4 | Confirmed; exact AK281180 match | No |

The additional CO1_1/ND1 products are alternatives, not established minor
alleles or low-abundance variants. Product index is not an abundance estimate.
No significant hit or insufficient alignment does not prove a sequence wrong.
The 15 lost Gryllus Yp2 products (median/min support 5/2) are single-copy
nuclear candidates and are explicitly outside the present recovery gate.
All-gene counts are 86 -> 64; repeated timing observations are not multiplied
into the unique product totals above.

Gryllus ITS_2 remains absent in every candidate replicate. Its baseline
978 bp sequence exactly equals frozen reference AK281180, SHA-256
`23163ca76f97f5d1685f0fcf4bda6ffc798bdf58c63953416833e78a480af68c`.
At threshold 4, the candidate evaluates **180 complete candidates**, withholding
all 180 because they touch pre-pruning SCC markers; none is withheld by a
collision-edge marker. Neither DFS nor path quota is reached there. At
threshold 2, graph extension reaches its node budget without connectivity.
The new diagnostic retains both outcomes instead of reporting only the final
node-budget failure. This rules out the repaired collision-marker/quota
mechanisms as sufficient rescue for this cell, not every possible algorithmic
cause or repeat-copy uncertainty.

The next review should focus on supported reconstruction of this abundant
ribosomal locus and disposition of the other high-copy losses. More data or
unconditionally removing SCC safeguards is not a demonstrated solution.
Read-spanning evidence is a candidate for a targeted follow-up; existing
`--read-threading` does not resolve repeat copy count. Do not expand this into
general metagenomic diversity or single-copy nuclear reconstruction without
separate review. #153 remains open; no release, tag, master merge, or package
publication is authorized.

### Final evidence

The [final archive](../benchmarks/benchmark_results/v3.2-high-copy-final-20260911)
contains all 60 raw runs, controls, validation logs, clean-build receipts,
per-product classifications, scoped losses, threshold diagnostics, and
checksums. Independent Astra review verifies artifact/source/reference hashes
and agrees that the high-copy release hold remains warranted. Frozen panel
references remain calibration evidence, not independent held-out truth.

The timing driver did not fingerprint the BLAST executables before its first
classification. A separate post-timing audit pins BLAST/makeblastdb 2.17.0+
paths, binary hashes and versions, validator/helper identities, and existing
reference-database file hashes. Reclassifying all 60 saved results reproduces
every full per-product `reference_match` record, with tool/DB/source-result hashes unchanged
before and after. This is supplementary verification, not retroactive
executable attestation for the original classification. Its first attempt
had 54 locator failures from attempting to open FASTAs for zero-product genes;
the corrected `final-classification-audit-r2` succeeds for all 60. Both audit
attempts are retained; no Sharkmer measurement or original result was replaced.

## Longer-k hypothesis and preregistered sweep

### What the existing evidence establishes

The release comparisons above use **k=19**, the CLI default. The separate
k=31 coral 18S/28S controls do not test the affected insect loci. Before this
sweep, no longer-k comparison had tested recovery of the seven missing sequences.

The repeat-withholding policy introduced by #132 (`b728f14`) marks cyclic
strongly connected components before pruning and rejects paths touching those
nodes, even on their first visit. Its synthetic controls prevent incorrect
shortened products. The ITS_2 diagnostics directly show SCC-based rejection,
but do not establish that its exact old 978 bp sequence is among the 180
complete candidates. Drosophila 12S/16S_2 additionally reach node/DFS limits;
Heliconius ND1 has one SCC-rejected candidate and one survivor with DFS
exhaustion; Gryllus CO1_1 rejects two SCC-touched candidates at threshold 18
and returns one clean candidate at threshold 13. These are observed rejection
mechanisms, not a one-change-at-a-time attribution for every missing sequence.

Longer k could distinguish contexts merged at k=19, reducing cyclic regions
and conservative rejection. Conversely, it changes seed context and usable
k-mer support and may break low-support paths. This is a hypothesis, not
measured rescue. General motivation is the repeat-resolution/connectivity
trade-off described in the [SPAdes paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/).
The current encoding permits k up to 31. All affected panel primers retain
15 bases at each k in this sweep; longer k adds observed inward seed context,
not an automatic change to primer trim or permission to ignore sequence changes.

### Protocol frozen before measurement

Status: preregistered 2026-09-11; all 72 discovery invocations and independent
evidence review complete. No longer k qualifies for conditional confirmation;
#153 remains open. The original prospective protocol follows unchanged.
Sol prepares the isolated driver, Terra the cross-k analysis, and Astra
independently reviews protocol, harness, and evidence. No production code,
default k, primer panel, repeat policy, or release state changes.

- Discovery: **72 invocations**, three affected insect samples, k=19/23/27/31,
  both released v3.1.0 and reviewed candidate `0c38d6a`, three paired repeats.
  The candidate production tree matches documentation-only dev `9105f29`;
  the previously attested clean binaries are reused after identity checks.
- Same frozen first 1M records, unchanged whole insect panel, two physical
  cores, chunks=0, threading off, 1,800 s timeout and 40 GiB address-space
  ceiling. Only k changes. Prewarming remains outside timing.
- Serial paired measurements rotate k order across pair/sample blocks and
  alternate version order. Every discovery timing finishes before BLAST.
  Record actual command/stats/manifest k, tool and database identities, input
  and binary hashes, and preserve failures rather than replacing measurements.
- Compare versions at each identical k and compare both against fixed k19
  released/candidate memberships. Read/base/prefix identities must agree
  across k; k-mer occurrence parity is required only within the same k.
- Track all seven missing high-copy sequences, exact-reference ITS_2, every
  previously retained high-copy sequence, per-product reference classification,
  threshold/repeat/search diagnostics, wall time and peak RSS. Same-k parity
  obtained by both versions losing products is not rescue. Full-sequence hashes
  are primary; any independently demonstrated boundary equivalence is separate,
  never counted as exact restoration. Deferred nuclear results stay separate.
- Select at most one longer k for confirmation. Eligibility requires valid,
  repeat-stable results, restoration of at least one of the seven exact missing
  sequences in every replicate, and preservation of every k19 candidate
  reference-confirmed sequence in the affected samples. Rank eligible settings
  by exact ITS_2 restoration, number of the seven restored, total released-k19
  high-copy sequence retention, then lower k. Do not select by noisy timing.
- If eligible, confirm the selected k on the seven remaining non-insect samples
  using three paired repeats (**42 additional invocations**). Otherwise stop
  at discovery and report the unresolved gate. Selection does not waive other
  high-copy losses or resource increases and does not change the default.

The machine-readable discovery protocol SHA-256 is
`464f70848477c5580482d67963774b9b0182a82df0ab45d09e996f61ff7646af`.
Its full file and execution/analysis receipts will accompany the resulting
archive. These are historical calibration inputs, not newly held-out data.
The user reviews the findings before any release or assembly-policy expansion.

### Longer-k results

**Longer k does not repair the seven targeted regressions.** All 72 invocations
complete and all 36 same-k pairs have matching aggregate counts. Sequences and
full per-product classifications are stable across all three repeats within
each sample/k/version. No k>19 restores any of the seven exact missing sequences,
including ITS_2; released v3.1.0 also loses all seven at those longer k values.
These are three-sample totals, not the earlier ten-sample 71 -> 64 aggregate.

| k | Released high-copy products | Dev high-copy products | Dev retains of its 26 k19 sequences | Targeted seven restored | Exact same-k version parity |
| --- | ---: | ---: | ---: | ---: | --- |
| 19 | 33 | 26 | 26 | 0 | No |
| 23 | 34 | 27 | 24 | 0 | No |
| 27 | 28 | 28 | 24 | 0 | Yes |
| 31 | 37 | 30 | 24 | 0 | No |

The equality at k=27 is not recovery: both versions have lost the targeted
k19 products. No longer k satisfies the preregistered eligibility rule, so
the conditional 42 confirmation invocations are **not run**, rather than
silently choosing a different success criterion after seeing the results.

There are real, separate benefits worth retaining as evidence:

- Drosophila CO2_1 at k=23/27/31 changes from a 291 bp product classified
  `wrong_gene` against the frozen indexed panel to a 325 bp product exactly
  matching the expected CO2_1 reference AC254620 (`confirmed_product`). This
  is changed sequence membership, not restoration of one of the seven losses.
- Gryllus 28S gains a 636 bp `confirmed_product` at k=27/31.
- Drosophila ND1 gains a 240 bp `confirmed_gene_other_taxon` product at every
  longer k; this is gene support, not expected-taxon confirmation.
- Additional Heliconius/Gryllus 12S or 18S products lack sufficient reference
  confirmation. Product-count growth alone does not establish biological gain.

All **12 k19 dev `confirmed_product` sequences** survive at longer k. That
does not mean all reference-supported products survive: Drosophila loses
its 241 bp ND4 `confirmed_gene_other_taxon` product at every longer k, as well
as the old CO2_1 sequence. Drosophila retains 7/9 old dev sequences; Heliconius
retains 6/6 and Gryllus 11/11. No changed sequence is silently accepted as a
primer-boundary equivalent. At identical k, dev has no sequences absent from
the released version: k23 additionally withholds one Heliconius 12S and six
Gryllus 28S baseline products; k31 withholds seven Gryllus 18S_1 products.
Those same-k withheld products have `no_significant_hit`, not proof of error.

### Longer-k runtime and memory

| k | Released sum of three median wall times (s) | Dev sum (s) | Dev vs released at same k | Dev vs its k19 |
| --- | ---: | ---: | ---: | ---: |
| 19 | 205.09 | 198.87 | -3.03% | baseline |
| 23 | 191.05 | 188.82 | -1.17% | -5.05% |
| 27 | 185.54 | 179.21 | -3.41% | -9.89% |
| 31 | 176.49 | 173.96 | -1.43% | -12.53% |

Higher k reduces runtime in this scoped experiment, but also changes accepted
k-mer workload, seed context, graph searches, and output membership; this is
not an isolated counter or repeat-resolution speedup. Median RSS stays around
2.13 GiB for Drosophila/Heliconius and 4.26 GiB for Gryllus, with only small
changes. There is no meaningful RAM improvement or demonstrated laptop limit.
All per-cell/per-pair increases and measurements remain in the archive.

### What the sweep adds to causal understanding

Gryllus ITS_2 diagnostics are identical across repetitions at each k:

| k | Threshold 4 | Threshold 2 |
| --- | --- | --- |
| 19 | 148 pre-pruning SCC nodes; 180 complete candidates, all SCC-rejected | Node budget reached without connectivity |
| 23 | No SCC nodes; connectivity found but no in-range candidate; maximum-length bound encountered | Node budget reached without connectivity |
| 27 | 1,912 SCC nodes; connectivity found but no in-range candidate; maximum-length bound encountered | No connectivity; node budget not reached |
| 31 | 1,915 SCC nodes; connectivity found but no in-range candidate; maximum-length bound encountered | No connectivity; node budget not reached |

Thus larger k can remove SCC evidence at one setting without restoring a
valid-length product; SCC size is not monotonic in k in these actual searches.
At k23/27/31 no complete in-range ITS_2 candidate is repeat-rejected. This
does not prove that changing length limits, coverage, or graph budgets would
recover the target, or that the exact baseline sequence exists in those graphs.
The lower-threshold graphs are independently constructed, not guaranteed
supersets of higher-threshold graphs. No primer-seed-specific count diagnostic
in these receipts supports attributing the failure to seed discovery alone.

A supplementary reference-only probe verifies that the exact 978 bp AK281180
sequence has **961 distinct oriented 18-mers with no duplicates**; its 22-,
26-, and 30-mers are likewise unique. The isolated reference path therefore
does not revisit an identical oriented node even at k19. This is not a read
graph: other templates, orientations, sequencing errors, and graph construction
can still produce cycles. It does not demonstrate that the exact old sequence
is a recoverable clean path in the current graph. The next targeted diagnostic
is to trace that sequence through seed discovery, extension, pruning, and
length/repeat checks, rather than assuming its own exact repeats cause failure.

### Sweep evidence and decision

The [sweep archive](../benchmarks/benchmark_results/v3.2-k-sweep-20260911)
preserves all 72 raw runs, immutable measurement results, per-product BLAST
outputs, raw stats/manifests/FASTAs, interleaved schedule, reviewed protocol,
analysis, summaries, reference probe, and tool/checksum receipts. Ten driver
tests and four analyzer tests pass; all eight real-binary actual-k adapter
smokes pass, including deliberate wrong-k rejection. Astra independently
validates measurement integrity, classification stability, and no-selection.

BLAST/makeblastdb 2.17.0+ binaries, reference databases, measurement-result
files, and source identities are checked before/after classification, which
starts after the last timed invocation. The analysis `protocol_sha256` field
is a canonical-JSON digest; the raw preregistered protocol digest is the
`464f7084...7646af` value above and in execution receipts. These are distinct
representations of the same protocol, not interchangeable hashes.

Driver bytes match the approved launch checkpoint and all five frozen copies.
An unused-helper working-file edit was reported during the approval handoff
and reverted; the archive's development-race note records this explicitly.
Do not claim continuous working-file immutability or independently attested
running Python bytecode. No measurements were replaced or rerun to hide it.

**Keep k=19 as the unchanged default and #153 open.** Longer k has useful
calibration gains and lower runtime, but neither targeted-regression rescue
nor universal high-copy retention. Do not release or broaden assembly policy
without user review. Metagenomic/nuclear recovery remains deferred.

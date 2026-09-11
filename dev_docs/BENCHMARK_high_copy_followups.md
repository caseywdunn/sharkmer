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

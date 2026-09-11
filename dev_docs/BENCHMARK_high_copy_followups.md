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

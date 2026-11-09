# Inclusive release comparison results

Date: 2026-09-30. Candidate: `ed04bdec5178d055b47cc8085135ee5ae2600998`.
Release baseline: v3.1.0, `5a664680c91ad59f32b8c2a847b8fef37f34a0ae`.

**All 60 invocations pass. Default inclusive output delivers 126 records:
64 standard and 62 uncertain.** Of the 71 released high-copy identities,
64 remain standard, five are available as uncertain hypotheses, and two
remain absent from the production output. The separate
[Drosophila diagnostic](DROS12S_RETENTION_RESULTS.md) finds both missing
identities in completed search paths, discarded by the retention budget.

Development integration is complete: `dev` includes `ed04bde`, and
`baseline/pre-inclusive-dev-20260929` preserves `956e771`. This comparison
does not authorize a release or establish biological repeat correctness.

## Delivery and historical identity availability

Counts below represent one deterministic output per sample; the three
repetitions agree. They are not totals inflated by repeated runs. Within
each sample, emitted records and unique target/oriented-sequence identities
agree. Targets are sample/target combinations, not inferred biological loci.

| Sample | Released records | Candidate standard | Candidate uncertain | Candidate total |
| --- | ---: | ---: | ---: | ---: |
| Xenia, SRR9278435 | 7 | 7 | 0 | 7 |
| Agalma, SRR25099394 | 5 | 5 | 0 | 5 |
| Rhopilema, SRR8617500 | 5 | 5 | 0 | 5 |
| Human, SRR17535371 | 5 | 5 | 0 | 5 |
| Nomeus, SRR22396603 | 4 | 4 | 0 | 4 |
| Liriodendron, SRR25378184 | 6 | 6 | 27 | 33 |
| Acer, ERR14009273 | 6 | 6 | 1 | 7 |
| Drosophila, SRR31887760 | 13 | 9 | 4 | 13 |
| Heliconius, SRR1057608 | 7 | 6 | 13 | 19 |
| Gryllus, SRR27962769 | 28 | 11 | 17 | 28 |
| **Total** | **86** | **64** | **62** | **126** |

Uncertain records are **49.2% of the complete default delivery** (62/126).
For the original high-copy scope alone, they are **48.0%** (59/123).
The prior three-insect interface evaluation's 34/60 fraction concerns a
different sample scope; it is not the denominator for this larger comparison.

| Released identity scope | Released | Standard exact | Uncertain exact | Absent |
| --- | ---: | ---: | ---: | ---: |
| High-copy | 71 | 64 | 5 | 2 |
| Deferred Gryllus Yp2 | 15 | 0 | 3 | 12 |
| **All targets** | **86** | **64** | **8** | **14** |

There are no reverse-complement-only recoveries. Candidate delivery also
contains **54 additional high-copy identities**, all uncertain, that were
absent from the released output in both orientations. Thus 126 products
means 72 exact historical identities plus 54 additional hypotheses; it
does not mean all 86 historical products were retained. The 15 Yp2 identities
remain explicitly outside the original high-copy comparison gate.

The seven previously lost high-copy sequences have the same disposition as
in the preceding inclusive-output assessment:

| Sample and target | Length | Current default delivery |
| --- | ---: | --- |
| Gryllus CO1_1 | 373 bp | Exact uncertain hypothesis |
| Gryllus ITS_2 | 978 bp | Exact uncertain hypothesis |
| Heliconius ND1 | 263 bp | Exact uncertain hypothesis |
| Drosophila 12S | 475 bp | Exact uncertain hypothesis |
| Drosophila 16S_2 | 590 bp | Exact uncertain hypothesis |
| Drosophila 12S | 487 bp | Absent; enumerated then discarded by retention |
| Drosophila 12S | 499 bp | Absent; enumerated then discarded by retention |

Standard sequences match the preceding September 11 candidate in all ten
samples. Standard target availability is 64, versus 68 in the release;
including uncertain hypotheses makes 74 sample/target combinations available.
Additional hypotheses do not establish additional validated haplotypes or
correct repeat lengths. Biological support is unassessed in this study.

## Runtime and memory

The frozen [protocol](INCLUSIVE_RELEASE_PROTOCOL.md) repeats the original
ten-sample, five-panel matrix with three alternating pairs on the same host,
two distinct physical cores, and exact historical inputs and panels. Both
pristine sources were compiled with the same Rust/Cargo 1.98.1 toolchain and
default hash backend. GNU time measures the complete CLI, with verification
and prewarming outside timing. Every invocation passes its 40 GiB address-space
limit, timeout, input/count checks, output validation and source identity checks.

| Sample | Release median seconds | Candidate median seconds | Change | Median RSS change |
| --- | ---: | ---: | ---: | ---: |
| Xenia | 42.26 | 42.83 | +1.35% | +0.0171% |
| Agalma | 37.89 | 37.95 | +0.16% | +0.0047% |
| Rhopilema | 43.87 | 46.53 | **+6.06%** | +0.0088% |
| Human | 1.55 | 4.29 | **+176.77%** | -0.6714% |
| Nomeus | 35.33 | 35.87 | +1.53% | +0.0071% |
| Liriodendron | 52.63 | 56.80 | **+7.92%** | +0.0002% |
| Acer | 24.58 | 25.30 | +2.93% | +0.0002% |
| Drosophila | 59.35 | 60.97 | +2.73% | +0.0097% |
| Heliconius | 41.59 | 44.07 | **+5.96%** | +0.0099% |
| Gryllus | 94.96 | 95.46 | +0.53% | +0.0125% |

All ten median wall times increase. Four exceed the protocol's 5% priority
review threshold; none exceeds its 2% peak-RSS threshold. The sum of sample
median wall times increases from 434.01 to 450.07 seconds (**+3.70%**), a
descriptive workload total that does not cancel individual regressions.
Median peak RSS is approximately 4.26 GiB for six samples, 2.13 GiB for
Acer, Drosophila and Heliconius, and 22 MiB for Human; exact per-sample
values are in the resource table.

Human's candidate times are 7.07, 3.75 and 4.29 seconds, compared with
1.55, 1.71 and 1.27 seconds for the release. The median difference is 2.74
seconds. Acer's first candidate time is 47.19 seconds, followed by 25.30 and
25.28 seconds; the outlier remains in the evidence. All repetitions and
paired deltas are retained. Three repetitions do not establish statistical
equivalence, and a shared host does not eliminate filesystem or scheduling
variation. These timings compare all intervening release changes together;
they cannot attribute overhead specifically to inclusive emission.

The [saved resource assessment](../benchmarks/benchmark_results/v3.2-inclusive-release-20260929/review/resource_assessment.json)
finds Human's summed user/system CPU median changes from 0.65 to 0.71 seconds,
while its wall median changes by 2.74 seconds; Human emits no uncertain
records. Rhopilema also emits none. Publication and filesystem waits are
possible contributors, but these receipts cannot establish a cause. The
four runtime flags remain open investigation priorities.

## Evidence and next decisions

The [evidence bundle](../benchmarks/benchmark_results/v3.2-inclusive-release-20260929/README.md)
contains the frozen protocol, exact input identities, source/build bindings,
all invocation receipts, class-aware tables, independent reviews and archive
verification. Input restoration preserves its initial failed attempt: two
checks had conflated a historical local gzip-cache SHA with the upstream
compressed encoding. A separately reviewed retry requires upstream MD5/size,
EOF where applicable, and the original exact uncompressed prefix identity;
all ten inputs pass. No input or benchmark scope was substituted.

The next production question is how to retain useful uncertain sequences
within explicit bounds when verbose per-path provenance consumes the budget.
The 12S evidence supports reviewing that retention policy; it does not select
new caps or justify changing standard admission. Runtime regressions also
need targeted investigation before a release decision. Wider panels, held-out
biological validation and the remaining dependency/release review remain
separate work.

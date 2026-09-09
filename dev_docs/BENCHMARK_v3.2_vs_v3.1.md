# v3.2 development versus released v3.1.0

2026-09-09. **Release remains on hold for user review.** Tracking:
[#153](https://github.com/caseywdunn/sharkmer/issues/153),
[#154](https://github.com/caseywdunn/sharkmer/issues/154), and
[#155](https://github.com/caseywdunn/sharkmer/issues/155).

## Bottom line

- No overall speedup or meaningful RAM reduction is demonstrated. Sum of
  primary per-sample median runtimes is **512.34 -> 515.91 seconds (+0.70%)**.
- There are important output regressions: **179 -> 71 primary products**,
  including **93 -> 7** across the three bacterial/metagenomic samples.
  All candidate sequences are exact baseline subsets, not replacements.
- Fewer products do not establish better precision or worse biological recall.
  Most losses lack reference truth, but one lost primary Gryllus ITS_2 product
  exactly matches the embedded 978 bp reference. Do not dismiss all losses as
  removal of incorrect repeat collapses.
- All **114 Sharkmer invocations exit zero**, without timeout or resource-limit
  failure. Read/base/accepted k-mer occurrence totals agree in all 57 pairs.
  This checks aggregate occurrences, not equality of every count-table entry.
- Peak RSS reaches **34.005 GiB in both versions**. The planned 16 GB hardware
  target still requires the v4.0 counting/storage work.

## What was compared

| | Released baseline | Candidate |
| --- | --- | --- |
| Source commit | `5a664680c91ad59f32b8c2a847b8fef37f34a0ae` (`v3.1.0`) | `ba64f573048b6c19a528278028810af5bd475b81` (`3.2.0-dev`) |
| Binary SHA-256 | `4cb93c72c830e052bf9bf03638eec2a060a699c5741ea0630fd9225ecc5b5e18` | `387af7374e567828ed334ac488a032d79933f4f34a1da3823edc6c3ea7275b82` |

Both binaries were built from clean `git archive` exports with
`cargo build --locked --offline --release`, default `ahashmap`, Rust/Cargo
1.98.1, and separate target directories. Source-tree hashes were checked before
and after building; archives, lockfiles, features, build environment, and binary
hashes are attested. Later evidence/docs commits do not change the tested code.

This replaces the *release-to-release comparison gap*, not the historical
record in [REVIEW_v3.2.md](REVIEW_v3.2.md). That earlier 71-product parity used
an intermediate repaired baseline, **not released v3.1.0**.

The protocol was frozen before timing: three alternating paired runs at a
1M-record cap for all 13 historical samples, then one paired run at 2M, 4M,
and 8M for each of the six configured cnidarian/insect samples. All primary
cells precede ascending deeper caps. This is 31 sample/depth cells, 57 pairs,
and 114 invocations, not 114 independent biological samples.

Each invocation uses k=19, two distinct physical cores (CPU 0/1), chunks=0,
unpaired sequential inputs, and read threading off. The same unchanged panel
files and frozen uncompressed four-line FASTQ inputs are used by both versions.
ENA URL order is retained; every staged 1M prefix matches the previously
verified calibration bytes. Four sources exhaust below 1M: human 108,518;
bacteria 525,982, 167,686, and 59,704 records. All deeper cells reach their cap.

Whole input SHA-256 is checked before each invocation, followed by prefix
prewarming; selected-prefix and binary/panel checks run after it. GNU time
measures wall time and peak process RSS. The host has 12 logical CPUs, six
physical cores, approximately 62.6 GiB RAM, and the `powersave` governor.
The 40 GiB `RLIMIT_AS` ceiling is virtual address space, **not an RSS cap**.
Each invocation has a 1,800-second timeout. BLAST database construction and
all per-product classification occur only after every timed invocation finishes.

## Runtime and memory

The primary sum of medians is approximately flat (+0.70%), not a demonstrated
speed improvement. The most notable larger-input primary slowdown is
Drosophila: **60.11 -> 64.10 seconds (+6.64%)**. Its candidate range is
60.86–64.14 seconds; candidate stage metrics place that within-version spread
in PCR rather than ingestion/finalization. That does not establish a specific
code-level cause. Liriodendron changes 60.76 -> 61.72 seconds (+1.58%).

Small jobs also slow: human 0.39 -> 0.43 seconds and coral metagenome
1.53 -> 1.67 seconds. Report absolute time as well as percentage; GNU time
resolution and fixed overhead matter at these durations. Several other cells
are slightly faster, but three repetitions on one host do not support strong
statistical or portable performance claims.

The deeper single-pair totals are **5,054.68 -> 5,118.93 seconds (+1.27%)**.
Drosophila 2M changes 131.29 -> 141.25 seconds (+7.59%); Agalma 4M changes
208.23 -> 219.17 seconds (+5.25%). These are descriptive observations requiring
replication before diagnosing small regressions. Peak RSS is about 4.256 GiB
at 1M, 8.505 GiB at 2M, 17.006 GiB at 4M, and 34.005 GiB at 8M across the
tested cells, with no meaningful version difference.

The [complete tables](../benchmarks/benchmark_results/v3.2-release-comparison-20260909/tables.md)
retain every cell's timings, ranges, RSS, actual records, and product changes.
Timings exclude downloading, application-cache management, gzip decompression,
prewarming/hash checks, and BLAST. They do not benchmark the gzip/cache fixes,
incremental counting, paired/threaded inference, or a real 16 GB laptop.

## Sequence and classifier outcomes

| Primary group | Baseline products | Candidate products | Lost |
| --- | ---: | ---: | ---: |
| Three bacterial/metagenomic samples | 93 | 7 | 86 |
| Three insect samples | 48 | 26 | 22 |
| Seven remaining samples | 38 | 38 | 0 |
| Total | 179 | 71 | 108 |

Sequence multisets and classifications are stable across all three primary
repetitions within each version. Comparison ignores output index for biological
gains/losses and records index changes separately. There are **no gained
candidate sequences in any of the 31 sample/depth cells**.

Primary lost-product classifications are: **86 no-reference, 20 no-significant
hit, one insufficient alignment, and one reference-confirmed product**.
Primary confirmed products change 40 -> 39. The bacterial loss has no reference
oracle and cannot be scored as either 86 corrected false positives or 86 proven
false negatives. Across all 31 representative cells, products change
470 -> 222 and confirmed products 117 -> 112; five confirmed-product observations
and one confirmed-gene/other-taxon observation are lost. Observations across
depths can represent the same locus and must not be counted as independent loci.

### Lost reference-supported products

All five expected-gene/taxon observations are **Gryllus ITS_2 (SRR27962769)**:

| Cap | Baseline sequence evidence against AK281180 | Candidate result |
| --- | --- | --- |
| 1M | 978 bp; 978/978 identities, exact full embedded reference | No ITS_2 product; final node-budget failure, earlier 20 repeat-touched paths withheld |
| 2M | 978 bp; 977/978 identities, one substitution | No product; repeat/path-cap and later node/DFS limits |
| 4M | 978 bp, 975/978 identities; second 966 bp product with a localized 12 nt deletion | No products; repeated path-cap withholding and node/DFS limits |
| 8M | 967 bp, 964 identities over 979 alignment columns; similar repeat-region deletion plus an insertion | No product; unresolved repeat plus exhausted search |

The remaining observation is **Drosophila 12S at 8M**: 419 bp, full-query
alignment at 95.487% identity to EU494495 from *D. tanythrix*, not the expected
*D. melanogaster*. It is other-taxon gene support, not expected-taxon confirmation.

The exact 1M ITS_2 match is a strong warning of reduced callable recovery.
However, repeat-rich reference agreement alone does not prove read-supported
repeat-copy reconstruction. Shorter 4M/8M variants could be collapsed products
or biological indels; this benchmark cannot distinguish them. We therefore
retain both the output loss and its uncertainty rather than declaring the
candidate biologically better or restoring every old product indiscriminately.

## Diagnosed follow-ups

Of 108 primary losses, five are in still-successful genes whose candidate
products are baseline subsets, consistent with standard first-valid-threshold
selection. The other 103 are in 14 genes that fail with repeat evidence somewhere
in stats/logs. Only 21 have the strongest terminal path-local statement that all
enumerated candidates touched repeat markers. Other failures combine graph-wide
repeat evidence with graph/node, DFS, or path limits.

1. [#154](https://github.com/caseywdunn/sharkmer/issues/154): enumeration can
   consume the 20-path cap **before repeat filtering**, potentially hiding a
   later clean route. Add a targeted graph regression; continue rejected paths
   within a bounded work budget without restoring unsafe repeat collapses.
2. [#155](https://github.com/caseywdunn/sharkmer/issues/155): preserve marker
   cause and threshold provenance, audit high-coverage existing-node vertex
   marking, and distinguish path-local repeat uncertainty from exhausted search.
   Gryllus final node-budget messages currently obscure earlier repeat withholding.

These are concerns requiring resolution or explicit review disposition, not
proof that every lost product is recoverable by either proposed change.

## Improvements demonstrated separately

Fresh known-synthetic probes against these exact binaries confirm that the
threshold fix recovers the 180 bp rare target missed by v3.1.0. The A18 control
remains exactly 138 bp. Old A19/A40/AC40 examples emit incorrect shortened
138 bp products; the candidate withholds them rather than claiming resolved
length. This demonstrates the specified safety/threshold behaviors, not a
general gain in real-sample sensitivity or metagenomic precision.

## Harness compatibility and evidence

The frozen isolated driver initially treated legacy abundance medians as
integers; v3.1.0 validly writes values such as `67.5`. This caused **22 parser
failures after successful CLI execution**. We did not change the active driver,
rerun timed commands, replace raw failures, or rewrite legacy headers. A separate
reviewed postprocessor revalidates the raw outputs, accepts integer/half medians,
preserves independently computed composite scores, and records original-result
SHA/failure lineage. All 114 normalized results validate and classify. The
original driver's nonzero suite exit therefore reflects this adapter defect,
not 22 Sharkmer failures. Both layers remain in the evidence archive.

Sol implemented the isolated harness/legacy correction; Terra independently
reviewed it and prepared analysis tooling; root reviewed artifacts and audited
real headers. Seventeen driver tests, eight compatibility-parser tests, tiny
synthetic analysis checks, and dual-version 100k smoke checks pass. Smoke timings
are excluded from the benchmark. No production code changed in this comparison.

See the [artifact index](../benchmarks/benchmark_results/v3.2-release-comparison-20260909/README.md)
for hashes, raw outputs, corrected lineage, scripts, limitations, and commands.
Inputs are historical calibration/regression data, not held-out truth. Leave
#129, dependency-alert triage, #153–155, and user release authorization open.
Incremental counting remains supported; no v4.0–v4.2 work is folded into v3.2.

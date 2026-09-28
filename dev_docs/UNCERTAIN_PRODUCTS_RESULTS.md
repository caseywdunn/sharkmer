# Uncertain-product emission: implementation and recovery assessment

Completed 2026-09-28 on `experiment/uncertain-products` for
[issue #157](https://github.com/caseywdunn/sharkmer/issues/157).
**Five of seven historical lost identities are now available as explicitly
uncertain FASTA records. Ordinary products remain exactly unchanged.**
Integration and release remain separate decisions.

## What was implemented

The public, opt-in `--emit-uncertain` flag publishes bounded completed,
in-range graph paths withheld for unresolved repeat structure in separate
`{sample}_{gene}.uncertain.fasta` files. Each header declares
`assembly_status=uncertain`, `repeat_length=unresolved`, and
`read_support=not_evaluated`, with sequence identity, causes, thresholds,
and search/retention flags. Stats preserve source paths and oriented marker
coordinates. Uncertain-only targets retain ordinary status `fail` and zero
ordinary products.

This is a production CLI/output implementation evaluated on an isolated
branch. It is distinct from the preceding standalone read-history diagnostic.
It makes hypotheses inspectable; it does not reconstruct repeat copy number
or provide new read support. Search order, thresholds, admission gates,
accepted-path quotas, and existing collection budgets are unchanged.
Ordinary-equivalent sequences are excluded; identical oriented uncertain
paths are merged, while reverse-complement alternatives remain separate.

FASTA, stats, and checksums share the existing publication transaction.
The validator verifies their agreement and excludes uncertainty-labelled
records from standard validation and reference bootstrap, including renamed
or mixed legacy files. See the frozen [output contract](UNCERTAIN_PRODUCTS.md).

## Exact historical recovery

All five matches below are literal, full-length, same-orientation sequence
matches. No historical sequence appears in ordinary output in any mode,
and there are no additional reverse-complement-only matches.

| Sample | Target | Historical length | New uncertain output | Source threshold | Retention truncated for target |
| --- | --- | ---: | --- | ---: | --- |
| Heliconius pachinus, SRR1057608 | ND1 | 263 | Exact match | 7 | No |
| Gryllus bimaculatus, SRR27962769 | CO1_1 | 373 | Exact match | 18 | No |
| Gryllus bimaculatus, SRR27962769 | ITS_2 | 978 | Exact match | 4 | Yes |
| Drosophila melanogaster, SRR31887760 | 12S | 475 | Exact match | 2 | Yes |
| Drosophila melanogaster, SRR31887760 | 12S | 487 | Not observed | — | Yes |
| Drosophila melanogaster, SRR31887760 | 12S | 499 | Not observed | — | Yes |
| Drosophila melanogaster, SRR31887760 | 16S_2 | 590 | Exact match | 2 | No |

All five matching hypotheses traverse pre-pruning cyclic-component markers;
the Drosophila matches also traverse retained collision-edge markers.
All five targets reached at least one search limit across attempted thresholds.
The header flags summarize the target, not just the matching source path.
Full hashes, identifiers, source records, and classifications are preserved
in [the seven-loss table](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/results/historical_losses.tsv)
and archived `comparison/runs/assessment.json`.

Drosophila 12S emits three hypotheses of 495, 489, and 475 bp. Neither of the
first two is an exact restoration of the historical 487 or 499 bp sequence.
At threshold 2, 26 withheld paths were observed, three retained, and 23
omitted by the marker budget. Graph node/DFS limits and node-visit restrictions
also apply. The two missing historical identities therefore remain
**unobserved under bounded search and retention**; these results do not locate
them among omitted paths or establish biological absence.

Gryllus ITS_2 emits twelve hypotheses: nine of 978 bp and three of 977 bp.
Equal length does not imply identical sequence. The exact historical 978 bp
identity is one of the twelve. The collection is not an exhaustive list of
repeat structures. The full panel's fixed run-budget shares allow only 12–13
path records per target before sequence/marker caps; unused shares are not
redistributed. No caps or panel membership were tuned after these results.

## Compatibility and resource observations

All nine preregistered runs completed. Baseline, new binary with flag off,
and new binary with flag on have identical read/base/k-mer counts, ordered
ordinary FASTA headers and sequences, ordinary gene outcomes, and existing
threshold diagnostics. New uncertainty metadata was assessed separately.

| Sample | Ordinary products in each mode | Additional uncertain hypotheses |
| --- | ---: | ---: |
| Gryllus | 11 | 17 |
| Drosophila | 9 | 4 |
| Heliconius | 6 | 13 |
| Total | 26 | 34 |

These totals cover three full-panel samples. They do not replace the earlier
71-to-64 high-copy release comparison, which had a different scope.
The 34 hypotheses are oriented graph sequences, not counts of validated loci,
molecules, or haplotypes. Per-target counts, lengths, and retention accounting
are in [uncertain_genes.tsv](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/results/uncertain_genes.tsv).

| Sample | Baseline seconds | Flag off seconds | Flag on seconds | Peak RSS GiB, baseline/off/on |
| --- | ---: | ---: | ---: | --- |
| Gryllus | 99.92 | 99.62 | 99.98 | 4.2552 / 4.2556 / 4.2552 |
| Drosophila | 61.26 | 60.91 | 63.91 | 2.1302 / 2.1306 / 2.1304 |
| Heliconius | 42.64 | 42.58 | 44.10 | 2.1306 / 2.1306 / 2.1304 |

Each sample has one rotated trio, so these are descriptive observations.
They do not establish statistical performance equivalence. Process peak RSS
is distinct from allocator peak and the enforced 40 GiB address-space cap.
The detailed [run table](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/results/sample_runs.tsv)
retains unrounded values and receipt identities.

## Methods and validation

The baseline is pristine `956e771437a2410a7bcbb1504b033c61dd1e11b6`.
Both release/default-hash binaries were built with project-local Rust/Cargo
1.98.1, identical locked offline build settings, and verified source snapshots.
Baseline binary SHA-256 is
`bbf418d8e3b83f7eb7bf94f294f906a666fcf92babede2c46772a5d31e3caec6`;
candidate is
`ef6074b4e817dcafce209b9ccce40ab33c70a6ebb1c47dd3b6240acc3f6c3ca2`.

The frozen full 21-target insect panel and exact previously consumed first
million R1 records per sample were used at k=19, two CPUs, chunks=0, and
without read threading. Inputs were prewarmed; invocations ran sequentially
with a 1,800-second timeout and observed CPU affinity/address-space enforcement.
Source, runtime, panel, binary, and read identities were checked before/after.
Slurm comparison job: `27727704`. The pre-execution freeze SHA-256 is
`ab139c02c2a14874fb66c41322268bdf5c690b837120c113f95f2ef2fe4cbb79`.

Formatting, Clippy with warnings denied, and release tests pass on both hash
backends: 246 unit and 26 integration tests per backend. Python ran 184 tests
with one existing skip and no failures; all 156 active references verified
across ten panels. Integration controls cover unchanged ordinary output,
flag combinations, metadata/coordinates, and publication/cleanup protection.
Independent parser controls, harness controls, execution clearance, and
[results review](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/results_review/final_results_review.json)
are preserved with the evidence. The first two compiled-check attempts exposed
test-fixture argument errors, which were corrected before the successful
full checks and candidate build; their failing logs remain archived.

The [evidence bundle](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/README.md)
describes source reconstruction, input restoration, complete outputs,
member checksums, and intentional archive omissions. Raw million-read inputs,
compiled binaries, and dependency/toolchain caches are not committed.

## Interpretation and next decision

The implementation meets the requested availability goal for five of seven
historical identities while preserving ordinary behavior in this comparison.
It is suitable to consider for integration as an **opt-in labelled output
feature**. Keep the default off and preserve the separate output/count contract.
This recommendation concerns output availability and compatibility, not
biological correctness of the sequences.

The five identities remain historical assembly hypotheses. In particular,
the Gryllus 978 bp ITS reference was a bootstrap assembly, not independent
full-length truth. Emission adds no read-linkage or repeat-length evidence.
The two missing Drosophila identities and read-supported reconstruction remain
separate questions under #148/#153. More coverage, paired/long reads, or
variation-aware mapping would need a separately specified assessment before
changing admission. The earlier read-history prototype remains diagnostic.

The branch is ready for the user's integration decision. No merge, release,
or closure of the held-out/all-panel/dependency review gates is implied.

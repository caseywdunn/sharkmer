# Benchmark dataset policy

Dataset roles are frozen before future assembly parameter tuning.

## Calibration and regression

- `benchmark.yaml` contains the historical real-data panel suite. These accessions and prior results have already informed implementation work, so they are calibration data.
- `known_truth.yaml` contains the included `ERR571460` 100k-read fixture, pinned by SHA-256. Its exact count and established amplicon lengths are executable offline regression checks, not estimates of biological sensitivity.
- Synthetic BLAST XML cases in `tests/fixtures/validation/` are validator correctness fixtures. They test identity, query coverage, wrong-gene, and split-alignment classification only.
- A panel accession label does not establish that its entire embedded sequence was fetched from that accession. Bootstrap amplicons annotated by partial BLAST hits are regression oracles, not independent full-product truth. Record deposited-region version/coordinates/strand separately from bootstrap provenance; see the [2026-09-13 reference audit](../dev_docs/BENCHMARK_reference_gates.md#reference-provenance-discovery).
- Active panel references now require exact extraction receipts against `panels/reference_sources.json.gz`; all 159 old entries are retained separately as regression-only history. The [complete provenance audit](../dev_docs/REFERENCE_PROVENANCE.md) finds 100 source-exact candidates and 59 source mismatches, then excludes two known gene-annotation conflicts. There are 99 active public regions (98 retained plus one new documentation example). Source mismatches alone do not establish biological incorrectness.
- Exactness applies to reference origin, not sample-to-reference agreement. Population SNPs/indels can preserve gene support; reference differences, unmatched regions, structural conflict, and missing reference evidence are separate outcomes. Haplotype truth and read support are not established by BLAST. Gene annotations remain curator-supplied, not certified by a sequence-source checksum.

Run the offline fixture in both modes:

```bash
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml -k 19 --scope counting-only --no-blast --executable target/release/sharkmer
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml -k 31 --scope end-to-end --no-blast --executable target/release/sharkmer
```

Omitting `--executable` rebuilds the release binary immediately before fingerprinting it. Normal timing excludes graph dumps; add `--diagnostic-graphs` to measure and label their overhead. `--cache-mode cold` gives sharkmer an empty per-run application cache for ENA inputs, while `warm` reuses its cache. Neither mode claims control of the operating-system page cache.

## Held-out evaluation

No existing biological dataset is genuinely held out: the historical accessions and included fixture have already been inspected or used during development. Biological held-out collection is therefore `reserved_not_yet_run`. Before any future assembly tuning begins, its accessions, immutable input checksums, read-subset rule, target truth, callable regions, and per-target tolerances must be added to a new versioned config and reviewed. Results must not be described as held-out until that registration occurs and the frozen inputs are run.

The deterministic XML fixtures are calibration regressions used while implementing the validator. They test variation-aware alignment, wrong-gene evidence, insufficient fragments, ambiguity, and conflicting arrangements; a split alignment alone is not a chimera. They are not held-out evidence and make no biological recovery claim.

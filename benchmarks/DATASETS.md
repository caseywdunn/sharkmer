# Benchmark dataset policy

Dataset roles are frozen before future assembly parameter tuning.

## Calibration and regression

- `benchmark.yaml` contains the historical real-data panel suite. These accessions and prior results have already informed implementation work, so they are calibration data.
- `known_truth.yaml` contains the included `ERR571460` 100k-read fixture, pinned by SHA-256. Its exact count and established amplicon lengths are executable offline regression checks, not estimates of biological sensitivity.
- Synthetic BLAST XML cases in `tests/fixtures/validation/` are validator correctness fixtures. They test identity, query coverage, wrong-gene, and split-alignment classification only.

Run the offline fixture in both modes:

```bash
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml -k 19 --scope counting-only --no-blast --executable target/release/sharkmer
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml -k 31 --scope end-to-end --no-blast --executable target/release/sharkmer
```

Omitting `--executable` rebuilds the release binary immediately before fingerprinting it. Normal timing excludes graph dumps; add `--diagnostic-graphs` to measure and label their overhead. `--cache-mode cold` gives sharkmer an empty per-run application cache for ENA inputs, while `warm` reuses its cache. Neither mode claims control of the operating-system page cache.

## Held-out evaluation

No existing biological dataset is genuinely held out: the historical accessions and included fixture have already been inspected or used during development. Biological held-out collection is therefore `reserved_not_yet_run`. Before any future assembly tuning begins, its accessions, immutable input checksums, read-subset rule, target truth, callable regions, and per-target tolerances must be added to a new versioned config and reviewed. Results must not be described as held-out until that registration occurs and the frozen inputs are run.

The deterministic XML fixtures are calibration regressions used while implementing the validator. They establish only that the harness rejects wrong genes, short fragments, and split products; they are not held-out evidence and make no biological recovery claim.

# v3.2 pre-release calibration evidence

These are **calibration/regression results**, not held-out biological validation
or a release announcement. See [the review](../../../dev_docs/REVIEW_v3.2.md)
and [dataset policy](../../DATASETS.md).

## Contents

- `baseline/`: six unmodified per-panel YAML reports, covering 13 warm-cache
  runs of the preserved executable historically labeled `b728f14`. This
  comparison point was used after #131/#132, not for untouched 3.1.0; the label
  is not verified build provenance.
- `candidate/`: six unmodified per-panel YAML reports covering the same 13
  samples at a one-million-record cap, built from clean `03c0fc6`.
- `known_truth_count/`: the pinned 100k-read fixture at k=19, counting only.
- `known_truth_pcr/`: the same fixture at k=31, with real reference BLAST.
- `comparison.json`: per-sample count/product/input-identity comparisons,
  timings/RSS, final classifier totals, and independently checked CLI probes.

Raw YAML paths preserve the original execution locations for provenance; they
are not portable links. Match their basenames to the files in these folders.
The final reports record actual Cargo artifact features, clean source identity,
executable SHA, selected cached-file SHA values, subset rules, and effective
commands. The explicitly supplied baseline binary has **unknown automated
source binding**: neither its top-level `git_commit` nor its workspace
observations prove the binary's build source. Compare its executable SHA, not
those workspace fields, when identifying the baseline.

All 13 comparisons preserve counts and all 71 product length/SHA values. The
five corrected human classifier labels result solely from previously missing
taxon metadata, not changed assemblies. All products are classified individually;
unreferenced, wrong-gene, ambiguous, split, insufficient, and no-hit outcomes
are not silently counted as validated products.

## Reproduction

Install the dependencies in `benchmarks/environment.yaml`, use the reviewed
source point, and run from the repository root. Omitting `--executable` lets the
harness build and fingerprint the actual release-profile artifact.

```bash
ulimit -v 12582912
python benchmarks/run_benchmark.py --config benchmarks/benchmark.yaml \
  -k 19 --scope end-to-end --threads 2 --max-reads 1000000 --cache-mode warm
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml \
  -k 19 --scope counting-only --no-blast --threads 2
python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml \
  -k 31 --scope end-to-end --threads 2
```

The limit shown is the Linux shell's 12 GiB virtual-address-space limit, not a
measurement of RSS. The original review redirected runs/reports to an isolated
`/tmp/sharkmer-v3.2-review/` workspace and reused its already-populated cache;
the normal commands use the repository's configured output/cache locations.
Populate the application cache before making a warm timing comparison and
check recorded cache state. Operating-system page-cache state is unknown.

These single-run timings do not establish a speedup. Four inputs contain fewer
than one million records. The configured 2/4/8M depth sweeps and biological
validation of all nine panels were **not** performed by this bounded matrix.
Register independent held-out truth under #129 before further policy tuning.

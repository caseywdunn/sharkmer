# Direct released-v3.1.0 comparison

See [the review report](../../../dev_docs/BENCHMARK_v3.2_vs_v3.1.md).
This is historical calibration/regression evidence, **not held-out biological
validation or release approval**. It is separate from `v3.2-review-20260909`,
whose baseline already contained some repairs.

## Contents

- `analysis.json`, `tables.md`: per-cell median/range timings, RSS, pairwise
  comparisons, index-independent sequence multisets, classifier changes, and
  repeat stability. Primary-only totals and all-depth totals are separate.
- `comparison.json`: normalized 114-invocation comparison and explicit list of
  22 legacy-reader corrections. All Sharkmer commands exited zero; these are
  parser compatibility corrections, not reruns or recovered biological products.
- `protocol.json`, `builds.json`, `inputs.json`: frozen settings, clean-source
  build attestations, and input/source/record-prefix receipts.
- `synthetic.json`: known threshold/repeat controls, separate from biological
  calibration and timing evidence.
- `evidence.tar.gz`: unchanged raw execution results, GNU time/stdout/stderr,
  stats, FASTA, current manifests, frozen harness receipts, original failed
  parser results, normalized results with lineage, ENA responses, and synthetic
  inputs/outputs. `ARCHIVE_CONTENTS.json` hashes every archived file.
- `tools/`: the exact isolated build/staging/harness/postprocessor/analysis
  scripts and focused tests. These are session-specific evidence tools, not a
  new supported Sharkmer CLI or a general portable benchmark framework.
- Logs retain preparation, original timing/classification, and post-processing
  activity. `SHA256SUMS` covers this directory's files except itself.

Large FASTQ inputs, compiled binaries, source exports, and Cargo target caches
are not committed. Their identities and reconstruction requirements are
recorded. The original staged data and full working artifacts remain at
`/tmp/sharkmer-release-comparison` on the measurement machine.

## Inspect and verify

From this directory:

```bash
sha256sum -c SHA256SUMS
INSPECTION=$(mktemp -d /tmp/sharkmer-evidence-inspection.XXXXXX)
tar -xzf evidence.tar.gz -C "$INSPECTION"
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tools -p 'test_*.py'
```

Use a fresh inspection directory; do not extract over the original measurement
workspace. Archive metadata timestamps are normalized; actual run timestamps
are in the JSON receipts. Python tests require PyYAML. The archive's raw records
retain their original absolute provenance paths instead of pretending they were
generated in the checkout or extraction directory.

## Measurement commands and reproduction boundaries

The original workspace used Python 3.12, BLAST 2.17.0, GNU time, taskset,
PyYAML/ruamel.yaml/jsonschema, Rust/Cargo 1.98.1, and cached locked Cargo
dependencies. Builds were clean exports of the two commits in `builds.json`.
The commands after build/input preparation were:

```bash
ROOT=/tmp/sharkmer-release-comparison
python "$ROOT/driver.py" --protocol "$ROOT/protocol.json" \
  --builds "$ROOT/builds.json" --inputs "$ROOT/inputs.json"
python "$ROOT/postprocess_legacy.py" --execution-root "$ROOT/execution" \
  --output-root "$ROOT/normalized"
python "$ROOT/analyze_comparison.py" --normalized-root "$ROOT/normalized" \
  --source-execution "$ROOT/execution" --output-dir "$ROOT/analysis"
```

The frozen driver returns nonzero because it rejects valid v3.1.0 half-integer
median headers. Run the postprocessor only after the complete source comparison
exists; it corrects that narrowly identified reader defect in a **separate**
output tree. It does not retry failed Sharkmer commands or rewrite original
results. The legacy composite score is independent of the abundance median.

These commands reproduce the original layout, not an automatic fresh-machine
setup. Scripts refuse to overwrite attested builds, incompatible schedules,
or existing normalized/analysis outputs. Preserve the original run and create
a separate workspace for new measurements. On another machine, adapt workspace
paths and CPU affinity, retain exact commits/panel bytes/algorithm parameters,
and generate new attestations rather than editing recorded evidence.

`stage_inputs.py` intentionally depends on the earlier verified calibration
receipts/cache to prove exact historical 1M-prefix continuity. Those large caches
are not in this bundle. Independently reconstructing data requires consuming
`inputs.json`'s ordered ENA source URLs, decompressing sequentially without
interleaving mates, retaining the first `requested_max_records` records or EOF,
and checking the full/prefix byte sizes, record/base counts, and SHA-256 values
against the receipts. The staged files are uncompressed four-line FASTQ. Partial
remote streams were not checked against full upstream compressed-file MD5s;
do not mistake recorded ENA MD5 metadata for a performed full-file verification.

## Interpretation limits

Three paired runs at 1M support descriptive medians/ranges; deeper cells have
only one pair. Four sources contain fewer than one million records. Inputs are
prewarmed, and downloading, gzip decoding, cache management, and BLAST are not
inside the measured command. A 40 GiB virtual-address-space ceiling is not an
RSS cap or validation on a 16 GB laptop. Aggregate k-mer occurrences agreeing
does not prove every abundance-table entry agrees. Reference classifier passes
are not proof of read-supported repeat-copy count or exact haplotype truth.

Review issues [#153](https://github.com/caseywdunn/sharkmer/issues/153),
[#154](https://github.com/caseywdunn/sharkmer/issues/154), and
[#155](https://github.com/caseywdunn/sharkmer/issues/155) before release.

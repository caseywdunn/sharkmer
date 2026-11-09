# High-copy follow-up evidence, 2026-09-11

This is a pre-release calibration comparison, **not release approval**.
See [the reviewed report](../../../dev_docs/BENCHMARK_high_copy_followups.md).

- Released v3.1.0: `5a664680c91ad59f32b8c2a847b8fef37f34a0ae`.
- Reviewed candidate: `0c38d6a051082cc2a1eb961b4206837f82948fda`.
- 60 fresh invocations: ten samples, three alternating pairs each, k=19,
  two physical cores, chunks=0, threading off, frozen local uncompressed
  prefixes. Human exhausts at 108,518 records; other samples use 1M.
- All invocations complete, all paired count totals agree, and sequence and
  classifier outcomes are stable within each version across repetitions.
- High-copy products: 71 -> 64, including the lost exact-reference Gryllus
  ITS_2. Fifteen nuclear Yp2 losses are deferred. No new sequences are gained.
- Median-runtime sum: 494.72 -> 485.78 s. Human is 0.05 s slower; RSS is
  effectively unchanged. These are descriptive measurements, not universal
  performance guarantees. The high-copy release gate remains open in #153.

## Contents

- `analysis.json` / `report.md`: per-cell timing/RSS, count parity, scoped
  sequence losses, support/reference evidence, and repetition stability.
- `diagnostic-audit.json`: cause/threshold audit of the five affected genes.
- `pre-fix-sequence-audit.json`: same input/panel-prefix identities and all
  64 candidate sequences unchanged from pre-fix development code.
- `comparison.json` / `provenance.json`: complete execution matrix, settings,
  build identities, frozen receipts, input hashes and actual prefix counts.
- `builds/`: clean candidate compilation output and manifest. The original
  release build is also attested in the prior September 9 evidence bundle.
- `evidence.tar.gz`: 1,000 original evidence files, including all timed
  runs, FASTA/stats/manifests, logs, reference databases, 13 control invocations,
  Rust validation logs, and supplementary classification audits.
- `tools/`: analysis/packaging/control tools. Executed measurement and build
  wrappers are frozen inside `final-execution/receipts/` in the archive.

The original classification did not pin BLAST executable identity. The
supplementary `final-classification-audit-r2` pins BLAST/makeblastdb 2.17.0+
binaries and records validator/helper identities. It verifies executable,
database, and source-result hashes before/after and reproduces every full
per-product `reference_match` record for all 60 results. The first
`final-classification-audit` is preserved with 54 zero-product FASTA-locator
failures; only the audit locator was corrected. No original result or Sharkmer
measurement was overwritten or rerun to obtain the successful audit.

## Integrity and reuse

From this directory:

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` records the size and SHA-256 of every archived file.
Large inputs, binaries, and source exports are not embedded; receipts identify
their checksums and source commits. Temporary absolute paths describe the
measured environment. Audit scripts depend on those paths; a new experiment
needs new protocol/build files and fresh output directories, not edited
historical receipts. Panel references are calibration evidence, not independent
held-out biological truth. Deeper inputs and all-nine-panel validation remain
outside this scoped follow-up.

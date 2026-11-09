# Bounded withheld-path diagnostics

Date: 2026-09-13. Issues #148/#153/#129. This is a diagnostic-only
development comparison, **not a release or a new v3.1.0 benchmark**.
See [results](../../../dev_docs/WITHHELD_PATH_DIAGNOSTICS_RESULTS.md) and
[protocol](../../../dev_docs/WITHHELD_PATH_DIAGNOSTICS.md).

Twelve commands compare development baseline `58ca2a6`, new default-off,
and new diagnostic-on behavior on three frozen first-million-record R1
prefixes at k=19, plus three focused diagnostic-on runs. All six full/focused
parity checks pass. Five of seven historical lost identities are retained
exactly in both diagnostic collections. Marker positions are not complete
ambiguity intervals, read-backed truth, or grounds for a gate exception.

## Archive contents

- `protocol.json`, `WITHHELD_PATH_DIAGNOSTICS.snapshot.md`: original frozen
  run plan and prose; the snapshot's incorrect 43-target statement is
  explicitly corrected in maintained documentation. The pinned panel has
  21 targets and all executable allocations used that actual count.
- `run_three_way.py`, `test_run_three_way.py`, `focused_panels/`: frozen
  sequential runner, controls, exact selected panel clones, and identity
  receipts. Input full-file/prefix hashes and actual resource-limit commands
  are bound in protocol and execution receipts.
- `baseline_build_receipt.json`, `changed_build/`, `build_changed.py`:
  baseline isolated-build verification, new default-feature build logs,
  exact compiler/Cargo versions, and all 34 compiled-source fingerprints
  and source snapshots. Frozen executable SHA-256 values are retained;
  compiled binaries and Cargo target trees are excluded.
- `execution/`: twelve fresh job directories containing complete raw stats,
  emitted FASTA, stdout/stderr, time/RSS reports, individual receipts,
  source bookends, aggregate execution, and comparison results.
- `historical_source/candidate_manifest.json`, `historical_lost_retained_map.json`:
  hypotheses from the preceding read audit, not orthogonal references.
  `historical_localization.json` is a historical-input-only summary;
  **it is not the new diagnostic outcome**.
- `withheld_results.json`, `analyze_withheld_results.py`,
  `test_analyze_withheld_results.py`: actual current-run identity matches,
  marker footprints, budgets/omissions, fail-closed raw-output/parity
  verification, and eleven analyzer tests.
- `smoke-execution*`, `execution-driver*`, `ATTEMPT_LINEAGE.md`:
  preserved failed schema/argparse attempts and successful fresh attempts.
  The rejected first launcher never starts a real benchmark job.
- `validation/`, `review/`: actual test logs, validation receipt, and
  independent pre-run/result checks. Tests/builds do not overlap timed jobs.
- `repository_snapshot/`, `repository_changes.patch`: documentation and
  source changes relative to the baseline. The documentation snapshot
  includes the transparent post-run target-count correction.

The multi-gigabyte staged FASTQ files are not duplicated. Restore them from
the public-input lineage in the
[preceding read audit](../v3.2-read-backed-audit-20260913/README.md) and verify
both full-file and consumed-prefix hashes. Later appended R2 records are
unused by these commands but are still bound by the full-file preflight.

## Verify and reproduce

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` binds each member's path, size, and SHA-256. Extract
only into a fresh directory. Absolute receipt paths document the original
environment; do not mistake an existing file at that name for verified input.

Build the baseline from `58ca2a6` and the changed code from its archived
compiled-source snapshot or the containing development commit. Use the
recorded default-feature release toolchain and compare binary fingerprints.
The runner requires PyYAML and the repository's validation package; scripts
also record the specific environment/tool executable hashes. Reproduction
on a different toolchain is a new comparison, not the original binary run.

The exact original execution command is preserved in `ATTEMPT_LINEAGE.md`.
The runner requires `--execute`,
`--diagnostic-arg=--diagnose-withheld-paths`, frozen source/build receipts,
and a fresh output directory. Do not overwrite archived outcomes. Restoring
only the consumed prefixes is enough to reproduce Sharkmer's read content,
but does not satisfy this wrapper's stronger full-staged-file identity check.

Timings/RSS describe single balanced trios, not statistical equivalence or
a demonstrated performance improvement. Full-length haplotype correctness,
read-supported repeat reconstruction, and release authorization remain open.

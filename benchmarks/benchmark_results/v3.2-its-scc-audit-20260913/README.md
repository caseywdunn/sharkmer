# Gryllus ITS SCC and local read audit

Date: 2026-09-13. Issues #148/#153/#129 remain open. **Not a release or
performance benchmark.** See [findings](../../../dev_docs/ITS_SCC_AUDIT_RESULTS.md)
and the [frozen protocol](../../../dev_docs/ITS_SCC_AUDIT_PROTOCOL.md).

Scratch instrumentation of development baseline `ae9189a` preserves focused
production outputs and exports both actual k19/threshold4 ITS graphs. The two
relevant cycles survive pruning. Four local sequence hypotheses, covering the
historical and extra-cycle traversals, have zero Q20 matches in both exact and
one-substitution modes. Independent exact recount confirms zero literal matches
even without quality filtering. Retained CO1 reproduces 27/28 records over 15
footprints. These results do not establish biological absence or justify a
gate exception; full-length haplotype and release decisions remain unresolved.

## Evidence

- `ITS_SCC_AUDIT_PROTOCOL.snapshot.md`, `inputs/`: pre-outcome protocol, exact
  historical hypothesis, two-target panel, prefix identity, unchanged CO1
  controls, and bindings to the preceding immutable evidence archives.
- `instrumentation.patch`, instrumented source copies, `build/`, `build-v2/`:
  scratch-only exporter, source fingerprints, compiler/Cargo/default-feature
  commands, executable hashes, and build/test logs. Binaries are excluded.
- `graph_capture/`, `graph_capture_v2/`: both attempts' complete graph files,
  raw Sharkmer outputs, commands, parity checks, and execution receipts. The
  first wrapper's phase-label validation failure is preserved, not a success.
  All three graph files are byte-identical to the successful second capture.
- `local_assays/`, derivation receipt, builder and tests: complete finite
  read-length alternatives, local-class/orientation mappings, and frozen
  scanner manifests. No claim of global repeat-copy exhaustiveness.
- `read_assay_plan.json`, `read_execution/`, `read_summary.json`: all six
  unchanged-scanner commands, source/input bookends, full reports and matching
  read ledgers, including zero and ambiguous counts. Exact and sensitivity
  scans reuse the same reads, not independent replicates.
- `review/`: independent SCC/walk checks, actual-manifest synthetic controls,
  source/ledger review, and the separately implemented full-prefix exact
  recount with its execution receipt and tests.
- `final_validation/`: 27 harness tests and 7 recount tests pass; scratch
  formatting/Clippy pass; repository Python suite runs 112 tests with one skip
  and no failures; all 156 active references pass local snapshot verification.
- Source mappings, attempt lineage, and repository documentation snapshots
  preserve reproduction context without rewriting original absolute paths.

The first build's invalid snapshot-hash field and the first graph wrapper's
failed phase-label check are superseded, not silently repaired in place.
Only the corrected build and successful capture feed the read assays.
The method and graph-derived manifests are frozen before the new read counts;
no post-outcome cap increase or relaxed scan is used for the findings.

## Verify and reproduce

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` binds every member's path, byte count, and SHA-256.
Extract into a fresh directory and verify member bytes before use. The archive
retains original absolute paths as provenance, not as a promise that files at
those paths are still the same. Source mappings identify archived copies.

Restore baseline `ae9189a7c8013c7160be616bdf3f6c668d5031f8`, apply the exact
instrumentation patch in a scratch checkout, and use the recorded toolchain
and default-feature commands. The original build wrapper also checks its
`git archive` source-tar hash; regenerate that tar from the pinned commit.
Use the per-stage wrappers and recorded argv with fresh output directories.
Python requires PyYAML and the repository validation package. Reproduction
with a different binary/toolchain is a new run, not the original attestation.

Raw FASTQ, executables, Cargo target trees, and nested older evidence archives
are intentionally excluded. Public-input lineage and original downloads are
bound by the [preceding read audit](../v3.2-read-backed-audit-20260913/README.md)
and [withheld diagnostic archive](../v3.2-withheld-diagnostics-20260913/README.md).
The assay consumes only the first million 150 bp R1 records, prefix SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
The wrappers additionally require the recorded full staged-file identity;
restoring only the prefix does not satisfy that stronger preflight.

Same-input read evidence is not an independent accessioned reference or
held-out biological validation. Production code, panel sequences, defaults,
and admission gates are unchanged. No release, tag, or master merge occurs.

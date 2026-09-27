# Ordered read-history experiment

Started and scanned 2026-09-26; reviewed and completed 2026-09-27 on
`experiment/its-read-threading`, fixed `dev` baseline `956e771`.
See [findings](../../../dev_docs/READ_HISTORY_RESULTS.md) and the frozen
[design/protocol](../../../dev_docs/READ_HISTORY_EXPERIMENT.md).

The standalone Python prototype preserves repeated edge visits, coordinates,
branch choices, quality breaks, and original read identity. Synthetic controls
show information lost by aggregate edge/local-pair counts. No production
assembly code or admission gate changes.

The same million Gryllus R1 records yield 773 graph-matching records before
pruning and 284 afterward. All four literal ITS classes are represented but
have zero complete exact/Q20 witnesses. The three previously identified
spanning records have identical histories across the graph stages; their
exact-walk breaks already occur before pruning. Known gapped ITS2 alignments
remain outside this prototype's sensitivity, not disproven biological evidence.

## Evidence layout

- `freeze.json` binds sources, input identities, exact commands, synthetic
  controls and independent clearance before scanning. SHA-256:
  `337e27452472e6bbe99371b0ee12e3d41709a0079c1677ead72ff249180c45e9`.
- `inputs/` includes byte-identical graph exports, four unchanged class
  windows, archived lineage, expected prefix identity, and preparation receipt.
  Retained self-context annotations are provenance only, not an alignment assay.
- `review/` preserves pre-execution review, 861 independent checks (including
  48 actual-class fixtures and 512 exhaustive matching cases), and 20 wrapper
  tests. The prototype's 25 implementation tests also pass.
- `results_review/` contains the independently tested full-prefix/ledger
  auditor, `provenance_ledger_review.retry1.json` pass receipt, known-read
  coordinate analysis, and final
  interpretation review. The auditor imports no prototype implementation.
- `evidence.tar.gz` packages this directory under `experiment/`, both full
  reports and retained raw-read/history ledgers under `run-v1/`, source
  snapshots, command receipts/timing, cluster logs, and repository document
  snapshots. Raw FASTQ, binaries and Python caches are excluded.

The original prefix is 365,356,898 bytes, SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
It is the same restored prefix used in the
[preceding recruitment experiment](../v3.2-its-recruitment-20260925/README.md).
No R2 or additional coverage was added. Absolute paths in receipts preserve
execution provenance; verify identities before trusting restored files.

## Verification and rerunning

Run `sha256sum -c SHA256SUMS` here. `ARCHIVE_CONTENTS.json` records every archive
member's size and hash; extract into a fresh directory and verify each member.
`package_evidence.py` verifies the archive by read-back before publishing the
inventory. Never overwrite original outcomes or mutate frozen sources.

Restore the prefix with the preceding experiment's `restore_prefix.py` into a
fresh location, verify its identity, and restore archived source versions.
Run `run_experiment.py --freeze <new.json> --fastq <restored.fastq> --run-dir
<fresh-directory>`; execute with `--freeze <new.json> --execute` on a compute
node with one CPU, 40 GiB and two hours for two individually bounded scans.
Then run `results_review/provenance_ledger_review.py --freeze </absolute/path/new.json>
--output <new-review.json>` with one CPU, 40 GiB and an 1,800-second limit.
The original run and independent audit have separate Slurm receipts.
The first audit failed in final artifact reporting when given a relative
freeze path; its failure receipt/log are preserved alongside the unchanged
checker's successful retry with an absolute path and fresh output name.

A moved checkout produces a new attestation, not a rewritten original freeze.
Do not tune eligibility or biological classes against these calibration
outcomes. Production replay/threading integration, variant-aware mapping,
admission review and a fresh high-copy comparison remain separate work.

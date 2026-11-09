# Inclusive-output verification, 2026-09-29

This bundle records the user-approved change from opt-in split uncertainty
files to default combined target FASTA delivery. See
[the contract](../../../dev_docs/INCLUSIVE_OUTPUT.md),
[results](../../../dev_docs/INCLUSIVE_OUTPUT_RESULTS.md), and
[product preference](../../../AGENTS.md).

Slurm job `27832506` completed the full checks and six default/opt-out runs in
`tmp/inclusive-output-20260929/attempt-v1`. Default output contains 26 standard
and 34 uncertain sequences, 60 total; opt-out output contains the same 26
standard sequences. Literal sequence parity with the preceding accepted
split-output evaluation passes. Five of seven historical identities appear
as uncertain hypotheses. Repeat correctness and performance equivalence are
not claimed.

Rust checks pass with 247 unit and 29 integration tests per hash backend,
formatting and Clippy. Python runs 197 tests with one existing skip; all 156
references verify across ten panels. The three comparison controls and seven
independent consumer controls also pass.

## Evidence map

- `verify_delivery.py`, `verify_delivery.sbatch`, `test_verify_delivery.py`:
  checked six-run harness and its controls. It reuses the immutable preceding
  bundle's compiled-check and input/receipt helpers.
- `review/source_review.json`, `review/test_combined_independent.py`, and adjacent
  log: independent source review and consumer controls.
- `review/outcome_review.json`: independent review of completed checks and
  sample outcomes, bound to the exact final verification receipt.
- `package_evidence.py`, `archive_receipt.json`, `evidence.tar.gz`,
  `SHA256SUMS`: exact source/output preservation and reopened archive hashes.

The archive's `MANIFEST.json` gives every payload's original path, size,
SHA-256, and inclusion reason. `payload/` preserves repository-relative paths,
including the complete successful attempt's `verification.json`, individual
run receipts, checks, logs, FASTAs, YAML outputs, and all 85 checked source
payloads. Directly compared prior outputs and their receipts are included;
the prior committed bundle remains unchanged.

External first-million FASTQ inputs, compiled binaries, Cargo targets, and
the local toolchain are excluded from the archive. Their recorded identities
remain in the original receipts. Small checked-source fixtures are retained.
The installed toolchain is Rust/Cargo 1.98.1; checks use Python 3.12.3 and
PyYAML 6.0.2. See the protocol and saved command receipts for exact arguments.

Verify `SHA256SUMS`, then extract into a fresh directory and check member hashes
against `MANIFEST.json`. Never extract over original outcomes. Archived docs
describe this checkpoint; later handoff edits can differ.

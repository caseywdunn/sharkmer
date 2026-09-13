# Read-backed high-copy audit

Date: 2026-09-13. Issues #129/#153. This is a **read-evidence audit**, not a
new Sharkmer throughput benchmark or a release. See the
[findings](../../../dev_docs/READ_BACKED_AUDIT_RESULTS.md) and
[preregistered protocol](../../../dev_docs/READ_BACKED_AUDIT_PROTOCOL.md).

The fixed candidate comparison contains seven historical losses and five
retained controls, deduplicated into 12 sequence identities. It uses the exact
first 1M R1 records previously consumed for each of three samples at k=19.
Six scanner invocations cover exact matching and the separate at-most-one-
substitution sensitivity assay for all samples. A read-bridge observation is
not a full-amplicon truth label; marginal diagnostics cannot relax a gate.

## Evidence layout

- `protocol/`: original v1 and pre-count v2 snapshots; v2 uniformly adds
  marginal diagnostics and strengthens distinct-read corroboration.
- `inputs/`: source download/staging lineage, full-file and consumed-prefix
  hashes, R1-only record semantics, and original result receipts.
- `candidates/`: frozen candidate sequences, all equivalent difference
  placements, intrinsic repeat intervals, source-repetition receipts,
  generator/tests, prior v1 snapshot, and Sol's future-design note.
- `scan_specs/`: 273 frozen event definitions, scanner/prefix bindings,
  adapter, and pre-execution command inventory. That inventory's
  `not_run_pending_independent_scanner_review` status is historical; the
  execution receipts below attest the completed commands.
- `execution/`: the first six outputs and failed wrapper-source attestation.
  The wrapper changed during execution; scanner/specs did not. Preserve this
  attempt as supplementary evidence, not a clean source-frozen execution.
- `ATTEMPT_LINEAGE.md`, `attempt_lineage.json`: the failed wrapper attestation
  and redundant, interrupted recovery attempt, without manufactured timings
  or overwritten scanner results.
- `execution_replication/`, `execution_replication_outputs/`: reviewed
  wrapper and clean, nonblind replication of all six unchanged assays in a
  fresh directory, with full-input preflight, per-job exit/command/hash/time
  receipts, complete matching-read ledgers, and before/after source hashes.
- `summary.json`, `summarize.py`: receipt-validated structural and marginal
  results. Union coverage of observed marginal windows is local occurrence,
  not their joint full-length haplotype.
- `review/`: independent Astra checks and findings; test logs and
  `source_tree/` preserve evaluated scanner/tests and final documentation.
- `repository_base.json`: starting commit and hash of the required
  historical source-results archive. No production assembly code changes.

The multi-gigabyte read files are not duplicated in this archive. Matching
read ledgers contain only the relevant records from these public inputs.
Original historical Sharkmer outputs remain in the
[September 11 archive](../v3.2-high-copy-final-20260911/README.md).

## Verify and reproduce

From this directory:

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` binds every member by path, size, and SHA-256. Extract
into a fresh directory. Do not overwrite earlier evidence or rerun the failed
first wrapper. Absolute paths in receipts describe the original environment,
not a requirement to trust an existing local file with the same name.

For one assay, restore the public R1 prefix and verify its exact bytes against
the archived expected-prefix manifest, then use the pinned scanner:

```bash
python3 source_tree/scripts/audit_read_bridges.py \
  --fastq /path/to/SRR27962769.first-million-R1.fastq \
  --events scan_specs/SRR27962769.scanner_manifest.json \
  --expected-prefix scan_specs/SRR27962769.expected_prefix.json \
  --max-substitutions 0 \
  --output /new/path/SRR27962769.exact.json
```

Repeat with `--max-substitutions 1` in another new output. The scanner accepts
uncompressed four-line FASTQ, requires a receipt-bound record limit, verifies
the consumed bytes, and refuses to overwrite an output. The replication
wrapper additionally verifies the entire historical staged file, including
unused later records; restoring only the consumed R1 prefix is sufficient
for the direct scanner command, but not that full-file preflight.

Run the repository regressions with
`python -m unittest discover -s tests -p 'test_*.py'` (the full suite needs
PyYAML and BLAST+; the read-bridge scanner/tests use only the standard library).
No gate change or new runtime/RSS comparison is implied by these audit runs.

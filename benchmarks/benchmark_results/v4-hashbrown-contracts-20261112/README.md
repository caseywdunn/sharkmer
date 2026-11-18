# Exact Hashbrown candidate conformance — 2026-11-12

Slurm job 33123876 passes 313 Rust tests per backend (281 unit and 32
integration), 217 Python tests, both Clippy configurations, formatting,
158 reference audits, ten schemas, nine exports, nineteen loader probes and
the pinned complete-sequence fixture. Before/after Rust/Cargo/checker hashes
agree with this source snapshot.

The replacement uses the already-present Hashbrown 0.16.1 dependency with an
explicit exact pin and the selected ahash/Fx hasher. Full u64 keys/counts retain
checked overflow, exact absence, bounded queries and move finalization. Before
growth, shared reservations cover the current allocation and the conservative
new bucket/control-byte bound. Allocation-size checks precede moving accepted
entries. Layout tests cover multiple growth sizes; budget failures preserve
counts, including updating an existing key when growth is unavailable.

Memory assembly access retains its concrete type and inlineable scan callback.
The shared fallible trait still covers opaque stores and injected access
failures. Full fixture records and count fingerprints agree at 1/2/4/8 threads
and through one-shot stdin in both feature configurations. Primer seeding and
graph thresholds retain their released distinction.

Owned byte measurements now use `owned_table_bytes`/`peak_table_bytes`, including
control bytes; the hidden table limit is `--count-table-budget-mib`. These
remain owned count-buffer limits, not process RSS or whole-run resource bounds.

The validated default binary is
`e72355620a4dda08827b69289742866a92fea17c4e973d8a956d6e7ed9c3ca41`.
Fresh paired measurements run separately in job 33124410 and the full ten-input,
two-hasher, 1/2/4/8-thread calibration matrix in job 33124697. Neither gate is
claimed passed here. Default selection remains unchanged pending those results;
the earlier custom-table attempts remain rejected.

Verify SHA256SUMS here and source-files.sha256 from the repository root at this
increment. Slurm's batch peak includes compilation/tests and is not application
RSS. The launcher retains local toolchain and dependency paths as provenance.

# Direct-counter compatibility refinement — 2026-11-12

Slurm job 33120256 passes 312 Rust tests per backend (280 unit and 32
integration), 217 Python tests, both Clippy configurations, formatting,
158 reference audits, ten schemas, nine embedded exports, nineteen loader dry
runs and the pinned full-sequence fixture. Source hashes match before/after.

Typed reverse complementation now uses constant work and is checked against a
scalar implementation at every supported k. Visiting immutable memory-store
keys no longer repeats canonicalization already established at insertion.
These changes preserve exact membership/overflow behavior and full fixture,
repeat-warning, supported-thread and one-shot-input equivalence.

The assembly adapter also restores the released distinction between primer and
graph thresholds: primer scanning sees the underlying entries and applies the
target's own minimum; graph lookups apply the global filter. A new regression
test covers a target whose minimum is lower than the global graph threshold
through both the legacy and exact-store adapters.

The comparison harness now normalizes historical uncertainty labels in both
observations and expectations, with regression tests against the retained
insect inventory. Tests also protect the unchanged aggregate/per-sample wall,
RSS and warmup-exclusion rules. This fixes a comparison error, without changing
the frozen protocol or biological claims.

These correctness results do not pass the performance gate. Both custom-table
attempts are rejected; see the independent audit and retained raw results in
the v4-direct-comparison-20261112 bundle. Default selection, external storage,
replay and complete resource ownership remain open.

Verify SHA256SUMS here and source-files.sha256 from the repository root at this
increment. The batch launcher retains its local toolchain and dependency paths.

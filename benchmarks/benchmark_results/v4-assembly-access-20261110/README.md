# Fallible exact-store assembly access — 2026-11-10

Slurm job 32789222 passed normal development checks: 300 Rust tests per hash
configuration (270 unit and 30 integration), 213 Python tests without skips,
both Clippy configurations, formatting, 158 reference audits, ten schemas,
nine embedded exports, nineteen loader dry runs and the pinned complete-sequence
fixture. The log records the same source hashes before and after execution.

The assembly boundary borrows a finalized exact store or an explicit legacy
adapter. Primer discovery visits entries without copying the global table;
membership/count queries propagate failures. Local directional primer matches
use a separate collection. Graph edges, coverage thresholds, path minima and
warning-source thresholds retain u64 values. Integer means use a wider sum to
avoid overflow before conversion to the existing floating-point statistic.

Tests compare actual exact-store PCR with the legacy path, including complete
record identities/descriptions/sequences, threshold diagnostics and warning
retention evidence across every existing repeat PCR fixture. A full 18S example
also retains its sequence after counts are scaled above u32::MAX; its emitted
minimum count and threshold remain larger than u32::MAX. Separate tests cover
directional strand identity, checked primer accumulation, threshold filtering,
and scan/lookup failures reaching the PCR caller as errors rather than no-data.

This increment does not select a new CLI counting backend. Normal ingestion,
replay, disk-store integrity and actual whole-run resource ownership remain
open under #139–#145. No candidate speedup or additional biological validation
is claimed. Slurm's batch peak includes compilation and checks; it is not
application RSS evidence. The frozen backend-promotion protocol remains in
the v4-performance-baseline-20261109 bundle.

`SHA256SUMS` binds this evidence bundle; `source-files.sha256` binds the source,
fixtures, scripts, panels and schema files at this increment. Verify the latter
from the repository root at the corresponding commit. The launcher retains its
repository-local toolchain/cache paths as execution provenance.

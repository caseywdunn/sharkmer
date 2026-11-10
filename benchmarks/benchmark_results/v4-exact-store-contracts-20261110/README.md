# Exact-store and shared-ledger conformance — 2026-11-10

The candidate passes normal development checks in Slurm job 32765482: 296 Rust
tests per backend (266 unit+30 integration), 213 Python tests without skips, both
Clippy checks, formatting, 158 reference audits, 10 schemas, 9 exports, 19 loader dry
runs, and the pinned complete sequence fixture. Both ahash and FxHash builds
now select the corresponding candidate table hasher and record it in metadata.

The new tests use an independent scalar count oracle, force hash collisions
and absent-key probes, include the zero-valued all-A key, reject mismatched k,
check per-key/total u64 overflow, exercise empty finalization, bounded batch
order/duplicates/errors, filtered visitation, and a moving finalize allocation.
Capacity failure preserves accepted counts. Resize accounting retains old and
new buffers simultaneously. Shared count/graph reservations reject combined
pressure and allow a retry after graph reservation release. Concurrent workers
cannot overcommit the ledger; guards release on normal exit and unwind. Disk
reservations have a separate ceiling and arithmetic/resize failure tests.

These are logical ownership reservations, not allocator hooks or measured RSS
limits. The exact in-memory candidate is not selected as the production
counter. Actual frontend, graph, replay and persisted-store producers still
need migration, and #139 remains open. See COUNT_STORE_CONTRACTS.md and PLAN.md
for those outstanding gates. The frozen calibration protocol precedes any
backend selection; no candidate speedup or biological accuracy claim is made.

The raw log binds changed source files, Cargo metadata and the shared checker
before and after execution. The source manifest binds all tracked code,
fixtures, panels, schemas and scripts plus new store/budget modules. The exact
batch launcher is retained; its repository-local toolchain/cache paths are
provenance. No binary or bulk input is copied into this small evidence bundle.
Verify SHA256SUMS here and source-files.sha256 from the repository root.

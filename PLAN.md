# v2.0 Development Plan

Execution order for v2.0 issues. Check off as completed. Run regression
benchmarks (#50) after each phase.

See ROADMAP.md for the full scope and rationale. See individual issues in the
[issue tracker](https://github.com/caseywdunn/sharkmer/issues?q=label%3Av2.0)
for detailed specifications.

## Phase 0 — Baseline and test infrastructure

Build testing and benchmarking infrastructure before any functional changes.

- [x] #18 Add CLAUDE.md (do first — provides context for all subsequent work)
- [x] #50 Build benchmark workflow and capture v1.0.1 baseline (Porites_lutea done; full suite across all samples pending — run on HPC)
- [x] #17 Fix broken test compilation
- [x] #20 Add workspace Cargo.toml at repo root
- [x] #19 Add GitHub Actions CI (depends on #17)
- [x] #30 Add test fixtures from ERR571460 (depends on #17, #19)
- [ ] #31 Deferred to Phase 6 — most tests would be invalidated by later phases (CLI #24, panels #33, error handling #21, FASTQ parsing #57). Add tests incrementally with each issue instead.
- [ ] Run benchmark, commit result

## Phase 1 — Quick independent fixes

Small, independent changes. Each should have its own branch from dev.

- [x] #38 Ensure deterministic output
- [x] #40 Add --cite flag
- [x] #42 Ensure FASTA line wrapping
- [x] #53 Free chunk hash tables after merging
- [ ] Run benchmark, commit result

## Phase 2 — Error handling

The biggest single item. Consider splitting into sub-branches:
kmer module → pcr module → main.rs.

- [x] #21 Replace panics/unwraps with proper error handling
- [x] #55 Fix unsigned integer underflow risks (depends on #21)
- [x] #60 Edge cases: empty inputs, degenerate primers (depends on #21)
- [x] Run benchmark, commit result

## Phase 3 — Core infrastructure

Depends on #21 (error handling). Order within phase by dependency, bug
fixes first.

- [x] #57 FASTQ parser validation
- [x] #56 Deduplication bug (make threshold configurable)
- [x] #22 Adopt log + env_logger, separate stdout/stderr
- [x] #23 Native gzip input support
- [x] #33 Refactor primer panels to YAML (do together with #25)
- [x] #25 Add --pcr-file for sideloading YAML panels (do together with #33)
- [x] #41 Validate input early
- [x] #45 Error when stdin is terminal
- [x] #52 Remove hardcoded singleton filtering
- [x] #26 Deduplicate code (do together with #58)
- [x] #58 Rename misleading get_canonical
- [x] Run benchmark, commit result

## Phase 4 — CLI, output format, performance

Depends on #22 (logging). Performance items (#27, #28, #29, #54) are
independent and can be interleaved.

- [x] #24 CLI flag improvements (do together with #35)
- [x] #35 Default to --chunks 0
- [x] #32 Structured YAML stats output
- [x] #34 Improve FASTA header format
- [x] #44 --color flag and NO_COLOR
- [x] #36 Parallelize sPCR across genes
- [x] #37 Progress indicators
- [x] #48 Warn on overwriting
- [x] #65 Add --sra flag for direct ENA download (depends on #24)
- [x] #27 HashMap for node lookup in graph construction
- [x] #28 Avoid storing RC in hash table (do together with #54, #52)
- [x] #54 Eliminate double hash lookups
- [x] #29 Pre-size hash maps, evaluate ahash
- [x] Run benchmark, commit result

## Phase 5 — Polish

Depends on Phase 4.

- [ ] #43 Summary line at completion (depends on #22, #32)
- [ ] #47 Peak memory reporting (depends on #43, #32)
- [ ] #39 Shell completions (depends on #24)
- [ ] #46 --dry-run mode (depends on #41, #24)
- [ ] #59 User warnings for common mistakes (depends on #22, #41, #23)
- [ ] #49 Polished CLI formatting (depends on #22, #43, #44)
- [ ] #67 Derive sample name from ENA metadata when --sra used without --sample
- [ ] #61 Update sharkmer_viewer (depends on #32)
- [ ] Run benchmark, commit result

## Phase 6 — Final validation

- [ ] #31 Update integration tests for final v2.0 CLI/output
- [ ] #50 Run full regression benchmark, compare to v1.0.1 baseline
- [ ] #51 Optimize benchmark downloads (optional, can defer)
- [ ] #64 sra_download.sh: default to R1 only
- [ ] #66 Update README examples for v2.0 CLI and test data
- [ ] #68 Document QIIME2/SILVA workflow and consider convenience flags
- [ ] Update CHANGELOG.md with all v2.0 changes
- [ ] Update bioconda recipe for new dependencies
- [ ] Tag v2.0.0-rc1, test bioconda build
- [ ] Tag v2.0.0 release

## Notes

- All work is on issue branches merged to dev. dev merges to master only
  for releases.
- Quality gates before merge: `cargo test`, `cargo clippy` (no warnings),
  `cargo fmt --check`.
- #31 appears in both Phase 0 (initial test coverage) and Phase 6 (update
  for final v2.0 behavior).
- #50 benchmark runs after every phase. Commit results to `benchmarks/`.
- Issues marked "do together" share code and should be in the same PR.

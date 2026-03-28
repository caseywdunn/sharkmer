# v3.0 Development Plan

Execution order for v3.0 issues. Check off as completed. Run regression
benchmarks after each phase.

See ROADMAP.md for the full scope and rationale. See individual issues in the
[issue tracker](https://github.com/caseywdunn/sharkmer/issues?q=label%3Av3.0)
for detailed specifications.

## Phase 0 — Benchmarks

Enhance benchmark infrastructure before any code changes. Capture v2.0.0
baseline at multiple coverage levels to measure the impact of later phases.

- [ ] Add multi-read-count sweep to benchmark suite (multiple subsampling
  levels per sample to measure coverage sensitivity — e.g., how many reads
  are required to recover single-copy nuclear genes). For insect and
  cnidarian samples, sweep max reads at 1M, 2M, 4M, 8M, 16M — these
  panels have single-copy nuclear genes (EF1A, EF1g, Fz4, Gpdh, Pgi, Yp2)
  that are the key targets for coverage sensitivity analysis. Run sweeps
  from high to low (16M, 8M, 4M, 2M, 1M) so the largest download populates
  the cache first and all smaller runs are cache hits.
- [ ] Add property-based integration tests that are algorithm-agnostic (e.g.,
  "recovered sequence aligns to reference with >99% identity") rather than
  golden-file tests that will break when algorithms change
- [ ] Capture v2.0.0 baseline benchmarks at multiple coverage levels
- [ ] Define measurable success targets for Phase 3 (e.g., recover 18S from
  ERR571460 with fewer reads than the current threshold)

## Phase 1 — Kmer pipeline optimizations

Internal changes deferred from v2.0. No result changes expected — benchmark
to confirm.

- [ ] Extract kmers directly from ASCII sequence in a single pass, eliminating
  the `Read` struct encoding/decoding round-trip
- [ ] Use `u32` for kmer counts instead of `u64` (~33% hash table memory
  savings). Saturate at u32::MAX.
- [ ] Run benchmarks, confirm identical results with performance improvement

## Phase 2 — Remote read caching

Cache reads downloaded from remote sources (currently ENA, other archives
may be added in the future). Reduces server load and enables reuse across
multiple runs on the same dataset. Internals should use generic terminology
(e.g., "URL read source") rather than ENA-specific naming, to support
additional archives later.

- [ ] Local cache for reads fetched from remote URLs (always cache by default).
  Default location: `dirs::cache_dir()/sharkmer/reads/` (e.g.,
  `~/.cache/sharkmer/reads/` on Linux, `~/Library/Caches/sharkmer/reads/`
  on macOS). Overridable with `--cache-dir`. Add `dirs` crate dependency.
- [ ] Cache metadata per entry (sidecar file): read count and whether the
  download was complete (unlimited) or partial (`--max-reads` limited).
  Cache hit logic:
  - Complete cached file: always a hit, regardless of `--max-reads` requested
    (there are no more reads to get)
  - Partial cached file: hit if requested `--max-reads` <= cached read count,
    re-download and replace otherwise
  - No cache entry: download
- [ ] Cache flags: `--no-cache` (skip reading from and writing to cache),
  `--clear-cache` (delete cache directory and exit)
- [ ] Log cache activity at info level: cache location, cache hit/miss per
  file, download vs reuse
- [ ] Run benchmarks, confirm no result changes
- [ ] Switch benchmark suite to use `--ena` with cached reads instead of
  pre-downloaded reads in `data/` directory. Benchmarks benefit from
  caching automatically — first run downloads, subsequent runs use cache.

## Phase 3 — Graph traversal

Replace ad-hoc graph heuristics with principled algorithms from the assembler
literature. All work in this phase is independent of read threading. Tackle
pruning improvements first (self-contained), then construction changes
(larger architectural impact, benefits from better pruning as a safety net).

### Pruning and path finding

- [ ] Replace heuristic ballooning detection with coverage-aware tip clipping
- [ ] Coverage-aware bubble popping (design scoring to be pluggable so
  Phase 6 can add read-support signal without restructuring)
- [ ] Replace backward-degree-based termination with principled traversal
- [ ] Improved handling of repeats (current cycle avoidance is too aggressive)
- [ ] Coverage-weighted best-path algorithm to replace `all_simple_paths`
  enumeration
- [ ] Better path scoring beyond `kmer_min_count` ordering
- [ ] #76 Fix O(N^2) dedup memory (compute distances on the fly instead of
  pre-allocating pairwise matrix)
- [ ] Evaluate #12 (duplicate product 0) — may be resolved by improved
  graph traversal and pruning

### Graph construction efficiency

- [ ] Build a single graph per gene seeded with all forward primer kmers
  simultaneously, instead of one graph per forward primer kmer (see
  [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#single-graph-seeded-with-all-forward-primer-kmers))
- [ ] Extend graphs incrementally across coverage threshold steps instead of
  rebuilding from scratch at each threshold
- [ ] Incremental histogram updates after each chunk merge

### Validation

- [ ] Run benchmarks, compare to Phase 2 baseline
- [ ] Evaluate against Phase 0 success targets

## Phase 4 — Read backend

Two-pass architecture for read threading (see
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#single-pass-vs-two-pass-read-ingestion-for-read-threading)
for analysis). Pass 1 counts kmers and builds graphs (unchanged). Pass 2
re-reads FASTQ to thread reads through the completed graphs. No graph
algorithm changes in this phase — just infrastructure. Remote reads are
already cached locally (Phase 2), so the second pass reads from cache.
Benchmark first to snapshot pre-threading state. The second pass is only
needed when sPCR is requested — if the run is histogram-only
(`--chunks >= 1` with no PCR), skip read retention entirely.

Read threading behavior by input source:
- **Local files**: two-pass, seek back to start for second pass
- **Remote URLs, cache enabled (default)**: download once to cache, read
  cache for both passes
- **Remote URLs, `--no-cache`**: download twice, log warning that reads
  will be fetched from server twice
- **stdin**: implies `--no-read-threading`, log info message explaining why

- [ ] Run benchmarks to snapshot state before backend changes
- [ ] Two-pass architecture: Pass 1 counts kmers (as now), Pass 2 re-reads
  FASTQ for threading
- [ ] `--no-read-threading` flag: skip second pass, single-pass kmer
  counting only. Automatically enabled for stdin input.
- [ ] File input: seek back to start for second pass
- [ ] Remote input: second pass reads from local cache (Phase 2). If
  `--no-cache`, re-download with warning.
- [ ] `--paired` flag: first file is R1, second is R2. Errors if not exactly
  2 input files (local or remote). Errors if used with stdin. When set,
  reads are ingested alternately from R1 and R2 so graph is built from the
  same balanced read set that will later be threaded and analyzed as pairs.
  `--max-reads` applies to the total (e.g., 1000 = 500 from each; if odd,
  round up to next even). Not implicit for `--ena` since some accessions
  are single-end. Without `--paired`, multiple files are read sequentially
  as in v2.0 (no pairing assumed).
- [ ] Run benchmarks, confirm no result changes from backend refactoring alone

## Phase 5 — Read threading

Thread reads through the assembled graph, annotating edges with read support.

- [ ] Map reads to graph paths, annotate edges with read-support counts.
  Store both total and unambiguous support per edge to allow strategy
  changes without re-running Pass 2 (see
  [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#ambiguous-read-mapping-during-read-threading))
- [ ] Per-edge coverage from actual read support (distinct from kmer frequency)
- [ ] Paired-end link annotation: when R1 maps to one branch and R2 to
  another, record the linkage
- [ ] Run benchmarks

## Phase 6 — Threading-dependent traversal

Use read-support signal to improve graph operations that were designed to be
pluggable in Phase 3.

- [ ] Read-supported pruning: remove edges/nodes not supported by any full
  read path
- [ ] Read-aware bubble resolution: use spanning reads to resolve bubbles
  where coverage alone is ambiguous
- [ ] Paired-end path constraints: use insert size distribution to constrain
  physically plausible paths
- [ ] Run benchmarks, compare to Phase 3 and Phase 5 results

## Phase 7 — Cleanup and validation

- [ ] Final benchmark comparison across all phases
- [ ] Update integration tests for v3.0 behavior
- [ ] Update README and documentation
- [ ] Update CHANGELOG.md
- [ ] Update bioconda recipe for new dependencies (e.g., `dirs`)
- [ ] Tag v3.0.0 release

## Open questions

- **Traversal/threading coupling**: How much will Phase 3 pruning and scoring
  designs need to change once read threading is available in Phase 6? The
  plan assumes pluggable scoring is sufficient — if it turns out the
  algorithm structure itself needs to change, Phase 6 scope could grow
  significantly.

- **Success criteria**: What measurable targets define "done" for Phase 3?
  E.g., recover 18S from ERR571460 at 50k reads? Recover specific nuclear
  genes that currently fail? These should be defined in Phase 0 using the
  multi-read-count benchmark sweep.

- **Benchmark sweep dataset selection**: Which datasets to include in the
  multi-read-count sweep? Should include at minimum all datasets for which
  EF1 was recovered in published v2.0.0 benchmarks, possibly others.

## Notes

- Run benchmarks after every phase. Commit results to `benchmarks/`.
- Issues not yet created in the tracker — create as work begins on each phase.
- See [CONTRIBUTING.md](../CONTRIBUTING.md) for branching model and quality gates.

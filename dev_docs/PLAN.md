# v3.0 Development Plan

Execution order for v3.0 issues. Check off as completed. Run regression
benchmarks after each phase.

See ROADMAP.md for the full scope and rationale. See individual issues in the
[issue tracker](https://github.com/caseywdunn/sharkmer/issues?q=label%3Av3.0)
for detailed specifications.

## Session log

Brief notes after each phase/session for cold-start context. Most recent
first.

(No sessions completed yet.)

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
- [ ] Ensure benchmark comparison detects changes in recovered amplicon
  sequences across runs (regression testing against previous results)
- [ ] #102 (optional) BLAST validation: batch-submit all amplicons to NCBI
  blastn against nt, parse per-amplicon e-value, identity, top hit
  accession, gene name, and taxon into benchmark results. E-value threshold
  1e-50. Uses `git config user.email` for NCBI API. Separate optional step
  since it requires network and takes a few minutes.
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
literature. All work in this phase is independent of read threading.

Key design change: adopt annotation-only model (see
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#graph-annotation-model-preserve-structure-defer-decisions)).
Replace destructive pruning with light structural cleanup + annotation-informed
path selection. Bubbles and variants are preserved until sequence emission.

### Light structural cleanup (replaces destructive pruning)

- [ ] #89 Replace heuristic ballooning detection with light tip removal:
  remove dead-end tips shorter than k with very low coverage (sequencing
  errors). Preserve all bubbles and meaningful branching.
- [ ] Remove disconnected components not reachable from any start node
- [ ] #91 Replace backward-degree-based termination with principled traversal

### Annotation-informed path finding

- [ ] #90 Pluggable scoring interface for path selection. Takes kmer
  coverage as input; Phase 6 adds read support and phasing signals
  without restructuring. Bubbles resolved at path selection, not by
  graph editing.
- [ ] #92 Improved handling of repeats (current cycle avoidance is too aggressive)
- [ ] #93 Coverage-weighted best-path algorithm to replace `all_simple_paths`
  enumeration, better path scoring beyond `kmer_min_count` ordering
- [ ] #76 Fix O(N^2) dedup memory (compute distances on the fly instead of
  pre-allocating pairwise matrix)
- [ ] Evaluate #12 (duplicate product 0) — may be resolved by improved
  graph traversal and path selection

### Graph construction efficiency

- [ ] #94 Build a single graph per gene seeded with all forward primer kmers
  simultaneously, instead of one graph per forward primer kmer (see
  [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#single-graph-seeded-with-all-forward-primer-kmers))
- [ ] #95 Extend graphs incrementally across coverage threshold steps instead
  of rebuilding from scratch at each threshold; incremental histogram updates
  after each chunk merge

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
- [ ] #96 Two-pass architecture: Pass 1 counts kmers (as now), Pass 2
  re-reads FASTQ for threading. `--no-read-threading` flag to skip second
  pass. File input seeks back; remote input reads from cache (Phase 2);
  `--no-cache` re-downloads with warning; stdin implies `--no-read-threading`.
- [ ] #97 `--paired` flag: first file is R1, second is R2. Errors if not
  exactly 2 input files (local or remote). Errors if used with stdin. When
  set, reads are ingested alternately from R1 and R2 so graph is built from
  the same balanced read set that will later be threaded and analyzed as
  pairs. `--max-reads` applies to the total (e.g., 1000 = 500 from each;
  if odd, round up to next even). Not implicit for `--ena` since some
  accessions are single-end. Without `--paired`, multiple files are read
  sequentially as in v2.0 (no pairing assumed).
- [ ] Run benchmarks, confirm no result changes from backend refactoring alone

## Phase 5 — Read threading

Thread reads through the assembled graph. Annotation only — no graph
structure changes. See
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#read-retention-and-threading-mechanics-phases-4-5)
for full mechanics and
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#graph-annotation-model-preserve-structure-defer-decisions)
for the annotation-only model.

- [ ] #98 Map reads to graph edges via maximal contiguous runs of adjacent
  graph kmers. Annotate per edge:
  - `read_support_total` — every read whose run includes this edge
  - `read_support_unambiguous` — only reads mapping to a single
    unbranched path
- [ ] Branch-point phasing: when a contiguous run passes through a branch
  point (in-degree > 1 or out-degree > 1), record (incoming_edge,
  outgoing_edge) link with count. This captures read-scale haplotype
  structure.
- [ ] Run benchmarks

## Phase 6 — Threading-dependent path selection

Use read-support and phasing annotations to improve path scoring via the
pluggable interface designed in Phase 3 (#90). No graph structure edits.

- [ ] #99 Add read-support signal to path scoring: penalize edges with
  zero read support, prefer edges with high unambiguous support
- [ ] #100 Read-aware bubble resolution: use branch-point phasing to
  determine which bubble arms connect to which — resolves cases where
  kmer coverage alone is ambiguous
- [ ] #101 Paired-end path constraints: use insert size distribution to
  link branches further apart than a single read can span (longer-range
  phasing)
- [ ] Run benchmarks, compare to Phase 3 and Phase 5 results

## Phase 7 — Cleanup and validation

- [ ] Final benchmark comparison across all phases
- [ ] Update integration tests for v3.0 behavior
- [ ] Update README and documentation
- [ ] Update CHANGELOG.md
- [ ] Update bioconda recipe for new dependencies (e.g., `dirs`)
- [ ] Tag v3.0.0 release

## Open questions

- **Read retention filtering**: During Pass 2, should all reads matching
  any graph edge kmer be retained, or should highly abundant (repetitive/
  low-complexity) kmers be excluded from the matching to avoid retaining
  irrelevant reads? The graph itself filters out most repeat kmers, but
  some may survive into the graph within amplicon regions. Assess with
  real results — if retained read counts are reasonable without filtering,
  no extra logic needed.

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

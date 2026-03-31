# v3.0 Development Plan

Execution order for v3.0 issues. Check off as completed. Run regression
benchmarks after each phase.

See ROADMAP.md for the full scope and rationale. See individual issues in the
[issue tracker](https://github.com/caseywdunn/sharkmer/issues?q=label%3Av3.0)
for detailed specifications.

## Session log

Brief notes after each phase/session for cold-start context. Most recent
first.

**2026-03-29 — Phase 3 implementation (graph traversal)**
Replaced ad-hoc graph heuristics with principled algorithms. Graph
construction: single graph per gene seeded with all forward primer kmers
(#94), incremental threshold extension across coverage steps (#95).
Structural cleanup: coverage-aware tip removal (#89), reachability pruning
(bidirectional BFS from start/end nodes), removed topology-based termination
heuristics (#91). Retained coverage-ratio filtering during extension (10×
median guard). Path finding: coverage-ratio edge annotations, pluggable
composite scoring (median count, coverage CV, coverage-ratio penalty) (#90),
coverage-weighted DFS replacing `all_simple_paths` (#93), bounded repeat
traversal with `MAX_NODE_VISITS=2` per path (#92), O(N) dedup memory (#76).
Removed: `pop_balloons`, `remove_side_branches`, backward-degree checks,
`would_form_cycle`, pairwise distance matrix. All 65 tests pass (43 unit +
22 integration). Added design decisions section on graph traversal analyzing
v1 shortcomings, assembler literature, and Phase 3 approach.

**2026-03-29 — Phase 2 implementation (remote read caching)**
Added persistent local cache for reads downloaded from remote URLs (ENA).
New `src/cache.rs` module with `CacheConfig`: lookup by SHA-256(URL) key,
download-then-read strategy (full gzipped file cached, `--max-reads` applied
at read time), SHA-256 checksum verification on cache hit, YAML sidecar
metadata (via `serde_yml`). Three new CLI flags: `--cache-dir` (override
cache location), `--no-cache` (stream directly as before), `--clear-cache`
(delete cache and exit). New deps: `dirs` (platform cache dir), `sha2`
(checksums). Benchmark `run_benchmark.py` switched entirely to `--ena --cache-dir
benchmarks/data/cache`, removing local file lookup (`find_sample_data` and
`run_sharkmer` removed; only `run_sharkmer_ena` remains).
All 66 tests pass, no code changes to existing logic.

**2026-03-29 — Phase 1 completion (kmer pipeline optimizations)**
Replaced Read struct encode/decode round-trip with `kmers_from_ascii()` —
single-pass kmer extraction directly from ASCII sequence bytes (commit
602236b). Switched kmer counts from u64 to u32 with `saturating_add`,
reducing hash table memory per entry. Benchmark output now uses timestamped
subdirectories to preserve results across runs. Added local BLAST support
to `blast_validate.py` (discovers databases in `/db/`, falls back to NCBI
remote API). Added `benchmarks/environment.yaml` conda env with blastn.
Added `dev_docs/overview.md` architecture overview. All 55 tests pass,
no regressions expected (internal-only changes).

**2026-03-28 — Phase 0 sweep benchmarks and BLAST validation**
Full coverage sweep completed (commit b9a9835). 34/38 runs succeeded; 4 OOM
failures at highest read counts (Porites 8M/16M, Agalma 8M/16M, Gryllus 16M
— all exceed 32GB RAM). 289 amplicons BLAST-validated against NCBI nt.
Results: `benchmarks/results/2026-03-28_sharkmer_3.0.0-dev_b9a9835.yaml`.

Key sweep findings for nuclear gene recovery:
- **cnidaria EF1A**: Porites recovered at 1M–2M (not at 4M — likely graph
  complexity at higher coverage); Agalma not recovered at any level; Rhopilema
  recovered at 4M and 16M (intermittent).
- **insecta Yp2**: Drosophila recovered at 16M only. No other insect nuclear
  genes (EF1g, Fz4, Gpdh, Pgi) recovered at any level.
- **More genes at higher coverage**: Drosophila goes 8→9 genes (1M→16M);
  Heliconius 7→10; Gryllus 13→15 (1M→4M). Diminishing returns above 4M
  for mitochondrial/rRNA genes.

**2026-03-28 — Phase 0 completion (infrastructure)**
Ran 1M baseline benchmark for all 14 samples with BLAST validation (commit
7944499). Results: `benchmarks/results/2026-03-28_sharkmer_3.0.0-dev_7944499.yaml`.
Restructured benchmark infrastructure: sweep levels are now per-sample in
config.yaml (cnidaria/insecta get 16M/8M/4M/2M/1M, others get 1M only). Data
files renamed to `{accession}.fastq` — one file per accession, `--max-reads`
controls coverage level. BLAST validation (product 0 per gene) is now a default
step in `run_benchmark.py` (skip with `--no-blast`). Coverage sweep benchmarks
(2M+) require >=16GB RAM — run on local machine with
`python benchmarks/run_benchmark.py`. Defined Phase 3 success targets
based on 1M baseline with sweep-dependent placeholders to fill after sweep.

**2026-03-28 — Phase 0 implementation (code + infrastructure)**
Implemented `--dump-graph` CLI flag with annotated DOT output (node shapes for
start/end/terminal, edge labels with kmer sequence and count). Added per-sample
multi-read-count sweeps to benchmark config (16M/8M/4M/2M/1M for insect and
cnidarian samples, high-to-low ordering). Updated `run_benchmark.py` to support
per-sample `max_reads` overrides and pass `--dump-graph` to all benchmark runs.
Enhanced `compare.py` to handle (sample, max_reads) pairs and report per-product
sequence fingerprint diffs. Added `blast_validate.py` for optional NCBI BLAST
validation of amplicons.

## Phase 0 — Benchmarks

Enhance benchmark infrastructure before any code changes. Capture v2.0.0
baseline at multiple coverage levels to measure the impact of later phases.

- [x] Add multi-read-count sweep to benchmark suite (per-sample `max_reads`
  in config.yaml). Insect and cnidarian samples sweep at 16M, 8M, 4M, 2M,
  1M — these panels have single-copy nuclear genes (EF1A, EF1g, Fz4,
  Gpdh, Pgi, Yp2) that are the key targets for coverage sensitivity
  analysis. All other samples run at 1M only. Sweeps run largest-first.
  Data files are `{accession}.fastq`; `--max-reads` controls coverage.
- [x] Ensure benchmark comparison detects changes in recovered amplicon
  sequences across runs (regression testing against previous results)
- [x] #102 BLAST validation: batch-submit all amplicons to NCBI blastn
  against nt, parse per-amplicon e-value, identity, top hit accession,
  gene name, and taxon into benchmark results. E-value threshold 1e-50.
  Uses `git config user.email` for NCBI API. Runs by default as the final
  step of `run_benchmark.py`; skip with `--no-blast`.
- [x] Implement `--dump-graph` flag: write per-gene annotated assembly graphs
  as DOT files (Graphviz) with kmer coverage, start/end status, terminal
  status per node/edge. Phase 5 adds read support and phasing annotations.
- [x] Enable `--dump-graph` in benchmark runs to capture baseline graphs
- [x] Capture baseline benchmarks at 1M reads for all 14 samples with BLAST
  validation. Full coverage sweeps (4M/2M/1M) completed for all sweep
  samples; 8M/16M completed for Rhopilema, Drosophila, Heliconius;
  8M/16M OOM-killed for Porites, Agalma (both levels) and Gryllus
  (16M only — 1.67GB genome). 289 amplicons BLAST-validated.
  Results: `benchmarks/results/2026-03-28_sharkmer_3.0.0-dev_b9a9835.yaml`
- [x] Define measurable success targets for Phase 3:
  - **No regressions at 1M** (baseline: 2026-03-28 commit 7944499):
    - Porites_lutea: 18S, 28S, 28S-v2, ITS, ITS-v2, EF1A (6 genes)
    - Agalma_elegans: 18S, 28S, 28S-v2, CO1, ITS, ITS-v2 (6 genes)
    - Rhopilema_esculentum: 16S, 18S, 28S, 28S-v2, CO1, ITS-v2 + 9 bacteria (15 genes)
    - Drosophila_melanogaster: 12S, 18S, 18S-v2, CO1-v2, CO2-v2, CO2, ND1, ND5 (8 genes)
    - Heliconius_pachinus: 12S, CO1-v2, CO1, CO2, CytB, ND1, ND5 (7 genes)
    - Gryllus_bimaculatus: 12S, 16S, 16S-v2, 18S-v2, 28S, CO1-v2, CO1, CO2,
      CytB, NADH, ND1, ND4, ND5 (13 genes)
    - All other samples: gene counts and identities must match baseline
  - **Single-copy nuclear gene recovery** (baseline from sweep):
    - cnidaria_EF1A: Porites recovered at 1M–2M (BLAST: no significant hit
      — may be divergent or misassembled); not recovered at 4M (graph
      complexity?). Agalma not recovered at any level (up to 4M tested).
      Rhopilema recovered intermittently (4M, 16M but not 8M).
    - insecta_Yp2: Drosophila recovered at 16M only (9 genes vs 8 at 1M).
      Not recovered for Heliconius or Gryllus at any level.
    - insecta nuclear genes (EF1g, Fz4, Gpdh, Pgi): not recovered at any
      level for any insect sample (up to 16M tested).
    - Phase 3 target: recover EF1A and Yp2 at the same or fewer reads than
      baseline; attempt recovery of currently-unrecovered nuclear genes
  - **Sequence stability**: product sequences for rRNA and mitochondrial
    genes should be identical (same MD5 fingerprints) after Phase 3 changes

## Phase 1 — Kmer pipeline optimizations

Internal changes deferred from v2.0. No result changes expected — benchmark
to confirm.

- [x] Extract kmers directly from ASCII sequence in a single pass, eliminating
  the `Read` struct encoding/decoding round-trip
- [x] Use `u32` for kmer counts instead of `u64` (~33% hash table memory
  savings). Saturate at u32::MAX.
- [x] Run benchmarks, confirm identical results with performance improvement

## Phase 2 — Remote read caching

Cache reads downloaded from remote sources (currently ENA, other archives
may be added in the future). Reduces server load and enables reuse across
multiple runs on the same dataset. Internals should use generic terminology
(e.g., "URL read source") rather than ENA-specific naming, to support
additional archives later.

- [x] Local cache for reads fetched from remote URLs (always cache by default).
  Default location: `dirs::cache_dir()/sharkmer/reads/` (e.g.,
  `~/.cache/sharkmer/reads/` on Linux, `~/Library/Caches/sharkmer/reads/`
  on macOS). Overridable with `--cache-dir`. Add `dirs` crate dependency.
- [x] Cache metadata per entry (sidecar file): whether the
  download was complete, and SHA-256 checksum of the cached file. Verify
  checksum on cache hit; treat mismatch as cache miss and re-download.
  Add `sha2` crate. Sidecar format: YAML (via existing `serde_yml`).
  Cache hit logic:
  - Complete cached file: always a hit, regardless of `--max-reads` requested
    (there are no more reads to get)
  - No cache entry: download
  - Checksum mismatch: re-download
- [x] Cache flags: `--no-cache` (skip reading from and writing to cache),
  `--clear-cache` (delete cache directory and exit)
- [x] Log cache activity at info level: cache location, cache hit/miss per
  file, download vs reuse
- [x] Run benchmarks, confirm no result changes (1M baseline: all 14
  samples match Phase 0 gene counts exactly)
- [x] Switch benchmark suite to support `--ena` with cached reads as
  fallback when local data files are absent. Local files in `data/` are
  preferred when present (faster). Cache stored in `benchmarks/data/cache/`.

## Phase 3 — Graph traversal

Replace ad-hoc graph heuristics with principled algorithms from the assembler
literature. All work in this phase is independent of read threading. See
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#graph-traversal-v1-shortcomings-and-phase-3-approach)
for detailed analysis of v1 shortcomings and the new approach.

Key design change: adopt annotation-only model (see
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#graph-annotation-model-preserve-structure-defer-decisions)).
Replace destructive pruning with light structural cleanup + annotation-informed
path selection. Bubbles and variants are preserved until sequence emission.

Implementation order: graph construction changes first (they simplify
everything downstream), then structural cleanup, then path finding.

### Step 1: Graph construction efficiency

- [x] #94 Build a single graph per gene seeded with all forward primer kmers
  simultaneously, instead of one graph per forward primer kmer (see
  [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#single-graph-seeded-with-all-forward-primer-kmers)).
  Refactor `do_pcr()` to pass the full `forward_primer_kmers` to
  `create_seed_graph()` instead of looping. Path finding runs from each
  start node to each end node within the single graph.
- [x] #95 Extend graphs incrementally across coverage threshold steps instead
  of rebuilding from scratch at each threshold. Build graph once at the
  highest threshold; at each subsequent (lower) threshold, add only newly
  qualifying edges to the existing graph. `break` on first product retained
  for now — will be removed when cross-threshold scoring (#90) matures.

### Step 2: Light structural cleanup (replaces destructive pruning)

- [x] #89 Replace heuristic ballooning detection with coverage-aware tip
  removal: remove dead-end tips shorter than k with coverage below a
  fraction of the local median (< 0.1× global median). Preserve
  all bubbles and meaningful branching. Removed `pop_balloons()`,
  backward-degree checks from `extend_graph()`. Retained
  `HIGH_COVERAGE_RATIO_THRESHOLD` (10× median) during extension as
  coverage-based guard against repeat regions.
- [x] Reachability pruning: after pruning, remove nodes and edges that
  cannot be part of any start-to-end path. Bidirectional BFS from start
  and end nodes; remove nodes not in the intersection. Subsumes
  `remove_orphan_nodes()`.
- [x] Remove disconnected components not reachable from any start node
  (subsumed by reachability pruning above).
- [x] #91 Replace backward-degree-based termination in `extend_graph()`
  with coverage-based decisions. Removed topology-based heuristics.
  Extension termination is now: end node reached, no qualifying
  successors, max-length exceeded, or MAX_NUM_NODES exceeded.
- [x] Coverage-ratio annotation: annotate each edge with count / global
  median. High ratios flag potentially repetitive edges; low ratios flag
  potential errors. Annotations feed into the scoring interface (#90).

### Step 3: Annotation-informed path finding

- [x] #90 Pluggable scoring interface for path selection. `PathScore`
  struct with kmer min/median count, coverage CV, max coverage ratio.
  `composite()` method combines signals. Phase 6 adds read support and
  phasing signals without restructuring.
- [x] #93 Coverage-weighted best-path algorithm to replace `all_simple_paths`
  enumeration. Iterative DFS exploring highest-count edges first. First
  paths found are highest-quality. Budget: `MAX_NUM_PATHS_PER_PAIR`
  paths per start node.
- [x] #92 Improved handling of repeats: replaced absolute cycle ban
  (`would_form_cycle()` BFS) with bounded repeat traversal during path
  finding. `MAX_NODE_VISITS = 2` per node per path. Cycles allowed in
  graph structure.
- [x] #76 Fix O(N^2) dedup memory: greedy clustering computes
  `bounded_levenshtein` on-the-fly against kept records only.
- [x] Evaluate #12 (duplicate product 0) — resolved by greedy
  Levenshtein dedup (#76) and sequential product re-numbering after
  dedup. No duplicate product IDs in any benchmark output.
- [x] #108 Unify forward and reverse primer handling. Remove the
  asymmetric reverse-complement step during preprocessing; both
  primers are processed identically (trim → expand → permute) and
  matched at the START of their respective kmers. Nodes are annotated
  as forward or reverse, paths emitted only from forward to reverse.
  Prerequisite refactor for #107 — makes reverse extension
  fall out naturally (forward seeds extend rightward, reverse seeds
  extend leftward).
- [x] #107 Reverse graph extension from forward and reverse primer
  seeds. Extend graph from both directions simultaneously; frontiers
  converge at the amplicon region and off-target seeds never meet.
  Naturally provides seed coherence, reduces path length through
  complex regions (exponential branching cut in half), and focuses
  the graph on the amplicon subgraph. Depends on #108. Benchmark
  results: gained 3 genes (Drosophila CO1, Gryllus 18S-v2, Covercrop
  16S-PRK341F) but lost 2 (Rhopilema CO1 and 16S-515F-Y-926R — both
  marginal cases near 50K node budget or DFS state limits). See
  diagnostic evidence in `tmp/porites_tests/ANALYSIS_16M.md`.
- [x] #105 Bounded seed evaluation before full graph extension. Two
  goals: (1) avoid wasted compute when no product exists; (2) avoid
  polluting the graph with spurious extensions that obscure the real
  product and exhaust node/DFS budgets. Approach: before full
  extension, give each seed a bounded local exploration (proportional
  to max_length), then evaluate structural signatures — local graph
  linearity, edge count consistency with primer kmer count. Seeds
  that fail are abandoned early, keeping their nodes out of the
  shared graph. If bounded exploration completes a full path to an
  opposite-direction seed, retain it but still evaluate remaining
  seeds (paralogs, alleles).

### Validation

- [x] Run benchmarks, compare to Phase 2 baseline (v2.0.0, bcb7818):
  14 new genes gained, 8 lost. Major gains: Drosophila 8→11, Gryllus
  13→15, Homo_sapiens 3→9. Losses: Agalma CO1 (pre-existing from DFS
  budget), Rhopilema -3 (16S pre-existing, CO1 and 16S-515F-Y from
  #107 shared node budget), Seawater 4→0 (pre-existing from Phase 3
  graph rewrite).
- [x] Evaluate against Phase 0 success targets: 4/6 target samples
  pass (Porites, Drosophila, Heliconius, Gryllus). Agalma misses CO1
  (pre-existing regression). Rhopilema misses 3 genes (1 pre-existing,
  2 from #107). Seawater regression is pre-existing. Nuclear gene
  recovery improved: Drosophila gains CO1+CytB+ITS at 1M reads.
  Sweep-level evaluation (2M/4M/8M/16M) not yet run on this branch.

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

**Design note for Phase 7 (runway) reuse:** The read retention
infrastructure should support filtering reads by primer kmer match (not
just graph edge match). Phase 7 needs primer-containing reads before
the full graph exists — during seed evaluation. Design the retention
API so it can be queried per-seed (e.g., "give me all reads containing
this primer kmer") rather than only per-graph-edge.

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

## Phase 5 — Read threading

Thread reads through the assembled graph. Annotation only — no graph
structure changes. See
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#read-retention-and-threading-mechanics-phases-4-5)
for full mechanics and
[DESIGN_DECISIONS.md](DESIGN_DECISIONS.md#graph-annotation-model-preserve-structure-defer-decisions)
for the annotation-only model.

**Design note for Phase 7 (runway) reuse:** The read-to-graph mapping
logic should work on arbitrary subgraphs, not just the full amplicon
graph. Phase 7 threads reads through bounded seed subgraphs using the
same mapping code. Keep the threading API graph-agnostic: accept any
`StableDiGraph<DBNode, DBEdge>` plus a set of reads, return annotations.

- [ ] #98 Map reads to graph edges via maximal contiguous runs of adjacent
  graph kmers. Annotate per edge:
  - `read_support_total` — every read whose run includes this edge
  - `read_support_unambiguous` — only reads mapping to a single
    unbranched path
- [ ] Branch-point phasing: when a contiguous run passes through a branch
  point (in-degree > 1 or out-degree > 1), record (incoming_edge,
  outgoing_edge) link with count. This captures read-scale haplotype
  structure.

## Phase 6 — Threading-dependent path selection

Use read-support and phasing annotations to improve path scoring via the
pluggable interface designed in Phase 3 (#90). No graph structure edits.

**Design note for Phase 7 (runway) reuse:** The read-support signals
used for path scoring (read count, consistency, unambiguous support)
are the same signals Phase 7 uses to evaluate seeds. Phase 7 asks
"does this seed have consistent read support?" rather than "which path
has the best read support?" — same data, different question. Keep the
signal computation separate from the scoring/decision logic so both
phases can reuse it.

- [ ] #99 Add read-support signal to path scoring: penalize edges with
  zero read support, prefer edges with high unambiguous support
- [ ] #100 Read-aware bubble resolution: use branch-point phasing to
  determine which bubble arms connect to which — resolves cases where
  kmer coverage alone is ambiguous
- [ ] #101 Paired-end phasing: when both mates of a pair map to the same
  amplicon graph in correct orientation (forward on one strand, reverse
  on the other), use the pair as a phasing link across the spanned
  region. Insert size is not known a priori and amplicons are only on
  the kb scale, so rather than estimating insert size distribution,
  simply require both mates to map to the same amplicon in valid
  orientation. This provides longer-range phasing than single reads
  without needing insert size calibration.
- [ ] Run benchmarks, compare to Phase 3 and Phase 5 results

## Phase 7 — Read-backed runway for seed evaluation

Reuse the read threading infrastructure (Phases 4-6) to replace the
kmer-table-only seed evaluation (#105) with read-backed evaluation.
For each seed, actual reads containing the primer kmer provide direct
evidence of whether the seed is real. See #110 for full design.

- [ ] #110 Read-backed runway: during Pass 2, collect reads matching
  each primer kmer. For each seed, thread its reads through a bounded
  local subgraph. Seeds with consistent read support (reads extending
  in the same direction, sharing overlapping kmers) are real. Seeds
  where reads diverge immediately are off-target. Seeds that pass
  get a pre-built read-backed subgraph ("runway") incorporated into
  the main graph, giving full extension a head start.
- [ ] Evaluate whether kmer-only seed evaluation (#105) should be
  retained as a fast pre-filter before read-backed evaluation, or
  replaced entirely.
- [ ] Run benchmarks, compare to Phase 3 (#105 kmer-only) results.

## Phase 8 — Performance optimizations

Hot-path performance improvements identified by code review. No behavioral
changes — benchmark to confirm identical results. See #103 for full analysis
and rationale.

### Easy wins (items 1–6 from #103)

- [ ] #103 Use entry API in `extend_with_histogram` to eliminate double
  hash lookup per kmer during chunk consolidation
- [ ] #103 Cache median edge count in `extend_graph` — recompute every N
  nodes instead of every iteration
- [ ] #103 Guard `summarize_extension` BFS descendants computation behind
  `log::log_enabled!(Level::Debug)` check
- [ ] #103 Replace `HashSet<u64>` with `[Option<u64>; 4]` stack array for
  candidate kmers in graph extension inner loop
- [ ] #103 Pre-allocate `KmerCounts` capacity from sum of chunk sizes
  before consolidation
- [ ] #103 Replace `get_path_length` per-call HashSet with bounded depth
  counter for cycle detection

### Larger refactors (items 7–9 from #103)

- [ ] #103 DFS path finding: replace path/visit_counts cloning with
  stack-based backtracking (push/pop instead of clone per state)
- [ ] #103 Avoid full graph clone for pruning — use in-place pruning on
  a copy-on-write structure or track removals separately
- [ ] #103 Implement byte-level lookup table for `revcomp_kmer` to reduce
  from O(k) to O(k/4) bit operations

### Memory reductions (#104)

- [ ] #104 Remove unused `DBEdge._kmer` field — reconstruct from node pair
  in `write_annotated_dot()` when needed
- [ ] #104 Drop `node_lookup` HashMap after graph extension completes
  (rebuilt between threshold steps anyway)
- [ ] #104 Add `kmer_last_base()` helper to avoid `kmer_to_seq()` String
  allocation per node during path assembly
- [ ] #104 Clone only histogram `Vec<u64>` instead of full `Histogram`
  struct (skip `FxHashMap` clone) in incremental counting
- [ ] #104 Stream histogram rows during output instead of materializing
  all histogram vectors simultaneously

### Parameter tuning (#109)

- [ ] #109 Revisit hard-coded graph parameters. Profile benchmark runs
  to understand how often each limit is hit and whether it causes
  regressions. Constants to review: MAX_NUM_NODES (50K shared budget
  — Rhopilema CO1 regression), MAX_DFS_STATES (100K — Rhopilema
  16S-515F-Y regression), HIGH_COVERAGE_RATIO_THRESHOLD, seed eval
  thresholds. Some may benefit from being derived from other
  parameters (e.g., node budget scaled by max_length).

### Validation

- [ ] Run benchmarks, confirm identical results to Phase 7
- [ ] Profile before/after to quantify gains

## Phase 9 — Cleanup and validation

- [ ] Final benchmark comparison across all phases (skip if already run
  at end of Phase 8)
- [ ] Update integration tests for v3.0 behavior
- [ ] Update documentation:
  - [ ] README.md (user-facing changes, new flags, updated examples)
  - [ ] CLAUDE.md (module descriptions, line counts, constants)
  - [ ] dev_docs/overview.md (architecture diagrams, data flow, design choices)
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

## Notes

- Run benchmarks after every phase. Commit results to `benchmarks/`.
- See [CONTRIBUTING.md](../CONTRIBUTING.md) for branching model and quality gates.

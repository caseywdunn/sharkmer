# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [3.0.0] - Unreleased

This release replaces the ad-hoc graph extension and pruning heuristics with
principled algorithms from the assembler literature. Results change but CLI
usage and output format are largely unchanged.

### Added

- **Bidirectional graph extension**: Both forward and reverse primer kmers
  seed the graph, and a single unified extension pass interleaves forward
  extension from start seeds with reverse extension from end seeds on a
  shared frontier. The two directions share one global node budget and one
  graph; a path is detected the moment forward and reverse extensions meet
  in the middle. Primer roles are symmetric — swapping the forward and
  reverse labels on a pair does not change the result, only the strand
  orientation of the product. This replaces the v2 forward-only extension
  that stalled in repetitive regions (#107).
- **Mismatch-aware primer kmer cap**: The `--max-primer-kmers` cap is now
  filled round by round, one mismatch level at a time. Exact-match kmers
  (0 mismatches) are added first; higher mismatch levels only contribute
  if capacity remains. Within each round, kmers are ranked by count
  (descending) then by kmer value (deterministic tiebreaker). Replaces the
  v2 count-threshold filter that filled the cap without regard to mismatch
  level (#118).
- **Dynamic global node budget**: The global node budget scales with data
  volume — 100K nodes at ≤150M bp, lerping up to 500K nodes at ≥750M bp —
  instead of the v2 fixed 50K default (#113).
- **Coverage-aware extension**: Edges with counts exceeding a configurable
  multiple of the local median (default 10×) are skipped during graph
  extension to avoid entering repetitive regions.
- **Frontier queue extension**: Graph extension uses a frontier queue
  instead of scanning all nodes each pass, eliminating O(n²) scaling.
- **Reachability pruning**: Forward BFS from start nodes and backward BFS
  from end nodes; any node not in both reachable sets is removed. Replaces
  the v2 orphan-removal and dead-branch heuristics with a single principled
  pass.
- **Coverage-aware tip clipping**: Dead-end tips are removed only if both
  short (fewer than k nodes) and low coverage (below 10% of graph median).
  High-coverage tips are preserved as potential real variants.
- **DFS state budget**: Explicit budget (default 100K states) and per-start
  path cap (default 20) prevent runaway exploration on complex graphs.
- **PCR failure reasons**: Per-gene stats YAML now reports why a gene
  failed: primer not found, node budget exceeded, or no valid path found.
- **Read caching**: Downloaded ENA reads are cached locally with SHA-256
  verification. Managed with `--cache-dir`, `--no-cache`, `--clear-cache`.
- **Panel versioning**: Primer panels now carry `name`, `version`,
  `description`, `maintainers`, and `changelog` metadata.
- **Panel validation blocks**: Panels include `validation` sections with
  sample accessions, read-depth sweeps, and reference sequences for
  systematic BLAST-based validation.
- **`expected_length` field**: Per-primer expected amplicon length for
  validation reporting.
- **New panel**: `c_elegans` primer panel.
- **Validation scripts**: `validate_panel.py` and SLURM batch scripts for
  systematic panel validation across taxa and read depths.
- Hidden CLI arguments for advanced tuning: `--node-budget-global`,
  `--high-coverage-ratio`, `--tip-coverage-fraction`, `--max-dfs-states`,
  `--max-paths-per-pair`, `--max-node-visits`, `--max-primer-kmers`.

### Changed

- **Default k changed from 21 to 19**. This allows primers up to 19 bp to be
  fully represented as single kmers, matching redesigned primer lengths.
- **Cnidaria CO1 primer redesigned** from 161 cnidarian mitogenomes:
  degeneracy reduced (48 → 32 kmer variants per pair), reverse primer
  extended to 19 bp. Recovery improved from 5/15 to 13/15 benchmark runs
  at 99.6%+ identity.
- **Panel cleanup**: Removed underperforming primer pairs from bacteria (8
  pairs removed) and insecta (2 pairs removed) panels based on validation
  results.
- **Path scoring**: Composite metric using median kmer count, coverage
  variance penalty, repeat-edge penalty, and (when available) read-support
  factor.
- **DFS edge ordering**: Branch points explored in descending coverage order
  so the best-coverage path is found first.
- **Deduplication scoring**: Paths sorted by composite score rather than
  raw kmer count.
- **Updated citation**: `--cite` and README now point to Dunn & Church
  (2026) Bioinformatics, btag163, the published sharkmer paper.
- Replaced `serde_yml` with `serde_yaml_ng` to resolve security advisory
  (#88).

### Performance

- **ahash for graph traversal**: Hot HashMap/HashSet data structures across
  graph extension, reachability pruning, DFS path finding, read threading,
  bubble resolution, and the histogram large-count tail are now ahash-backed
  instead of std SipHash. ahash is unconditionally pulled in; the
  `rustc-hash` dep is now optional, gated behind the `fxhashmap` feature.
- **Median computation in linear time**: The edge-count median recomputed
  every 1K nodes during extension now uses `select_nth_unstable` (expected
  O(E)) instead of `sort()` (O(E log E)). Same partition core is reused
  for all median callsites in `pcr/graph.rs` and `pcr/pruning.rs`.
- **EdgeIndex carried through DFS**: `get_assembly_paths` now records the
  `EdgeIndex` of each step alongside the node, eliminating three
  `find_edge` calls per edge per path in `generate_sequences_from_paths`.
- **SmallVec child frame**: The DFS child frame in `sorted_children` is a
  `SmallVec<[T; 4]>` instead of a `Vec`. De Bruijn out-degree is bounded
  by 4, so the SmallVec never spills and the per-call heap allocation
  is gone.
- **Single-pass reachability BFS**: Forward and backward BFS in
  `reachability_pruning` each push all start/end nodes onto one shared
  stack instead of restarting per start.
- **DFS uses explicit stack** with incremental push/pop instead of path
  cloning.
- **Pre-allocated data structures** throughout the PCR pipeline.

## [2.0.0] - 2026-03-22

This is a breaking release. CLI arguments, output file formats, and default
behaviors have changed to make the interface clearer, more consistent, and
to prepare for new features. See the Changed and Removed sections below for
details.

### Added

- **`--ena` flag**: Stream reads directly from ENA by accession — no SRA
  toolkit required (#65). Sample name is auto-derived from ENA metadata when
  `--sample` is omitted (#67).
- **Native gzip support**: `.fastq.gz` files are detected and decompressed
  automatically (#23).
- **`--pcr-panel`**: Select built-in primer panels by name (replaces `--pcr`
  for panel selection). Repeatable for multiple panels (#24).
- **`--pcr-panel-file`**: Sideload primer panels from user-supplied YAML files (#25).
- **`--pcr-primers`**: Specify primer pairs inline with key=value format
  (replaces `--pcr` for manual primers) (#24).
- **`--list-panels`**: List available built-in panels and exit (#24).
- **`--export-panel`**: Export a built-in panel as YAML for customization (#24).
- **`--help-pcr`**: Show detailed help for `--pcr-primers` format (#24).
- **`--dry-run`**: Validate inputs and show what would happen without
  processing (#46).
- **`--cite`**: Print citation information and exit (#40).
- **`--color auto|always|never`**: Control colored output; respects `NO_COLOR`
  environment variable (#44).
- **Shell completions**: `--completions` for bash, zsh, fish,
  elvish, powershell (#39).
- **Structured YAML stats**: `.stats.yaml` output with version, parameters,
  read/kmer counts, peak memory, and per-gene PCR results (#32).
- **Summary line**: Concise completion summary showing reads, kmers, genes,
  time, and memory (#43).
- **Peak memory reporting**: Track and report peak RSS via global allocator
  (#47).
- **Progress indicators**: Progress bars for read ingestion and PCR (#37).
- **Parallel sPCR**: Process genes concurrently with rayon (#36).
- **Overwrite warnings**: Warn before overwriting existing output files (#48).
- **User warnings**: Detect common mistakes — low read count for nuclear
  genes, missing panels, etc. (#59).
- **Input validation**: Validate file existence and FASTQ format before
  ingestion (#41, #57). Error when stdin is a terminal with no input (#45).
- **Gene name prefixing**: Output files use `{panel}_{gene}` naming for
  clarity when using multiple panels (#69).
- ROADMAP.md with plans for v2.0, v3.0, v4.0, and v5.0 releases
- CONTRIBUTING.md with development guidelines
- GitHub Actions CI (#19)
- Integration test suite with ERR571460 fixture (#30, #31)
- Regression benchmark suite with 14 SRA datasets (#50)

### Changed

- **CLI restructured**: `--pcr` split into `--pcr-panel`, `--pcr-panel-file`,
  `--pcr-primers`. `--sample` is required. `--verbosity N` replaced with
  `-v`/`-vv`/`-vvv`. `--chunks` defaults to 0 (skip histograms) (#24, #35).
- **Logging**: Adopted `log` + `env_logger`. Logs go to stderr, data to
  stdout (#22).
- **Primer panels**: Refactored from Rust code to YAML files in `panels/`,
  embedded via `include_str!()`. Shared parser for built-in and user panels
  (#33).
- **Hash backend**: Default changed from `fxhashmap` to `ahashmap` (AES-NI
  hardware acceleration). Removed `intmap` and `nohashmap` backends (#29).
- **FASTA headers**: Key=value metadata format with sample, gene, product
  index, length, and kmer count statistics (#34).
- **Histogram files**: Include header rows and comment lines with
  version/parameters.
- **Deduplication**: Made edit distance threshold configurable per primer
  pair via `dedup-edit-threshold` parameter (#56).
- Updated sharkmer_viewer for new histogram and stats formats (#61).
- Updated README with v2.0 CLI examples, `--ena` usage, and gzip support (#66).

### Fixed

- **Error handling**: Replaced all panics and unwraps with `anyhow::Result`
  propagation throughout (#21).
- **Deterministic output**: Sort kmer iterations for reproducible results
  across runs (#38).
- **Unsigned underflow**: Prevent subtraction underflow in coverage
  calculations (#55).
- **Edge cases**: Handle empty inputs, degenerate primers, division by
  zero (#60).
- **FASTA wrapping**: Enforce standard 80-character line wrapping (#42).
- **FASTQ validation**: Validate format on ingestion, reject malformed
  records (#57).
- **Deduplication bug**: Fix threshold comparison that could silently
  discard real variants (#56).

### Performance

- HashMap for node lookup in graph construction (#27)
- Avoid storing reverse complements in PCR kmer table (#28)
- Eliminate double hash lookups in graph extension (#54)
- Free chunk hash tables after merging (#53)
- Pre-size hash maps (#29)

### Removed

- `--pcr` flag (replaced by `--pcr-panel`, `--pcr-panel-file`, `--pcr-primers`)
- `--verbosity N` flag (replaced by `-v`/`-vv`/`-vvv`)
- `intmap` and `nohashmap` hash backend feature flags
- Hardcoded singleton kmer filtering before sPCR (#52)
- Misleading `get_canonical` method (renamed to accurate name) (#58)
- Duplicated code across kmer module (#26)

## [1.0.1] - 2025-01-23

### Fixed

- Added missing MIT LICENSE file
- Fixed build.sh path for bioconda build environment

## [1.0.0] - 2025-01-22

### Added

- Initial stable release
- Kmer counting with incremental histogram generation
- In silico PCR (sPCR) with seeded de Bruijn graph assembly
- Preconfigured primer panels: cnidaria, human, teleostei, angiospermae,
  insecta, bacteria, metazoa
- Manual primer pair specification via `--pcr` flag
- Support for reading FASTQ from files or stdin
- Configurable kmer length, chunk count, thread count
- Multiple hash map backends via feature flags (fxhashmap, intmap, nohashmap)
- Bioconda recipe for binary distribution
- sharkmer_viewer Python tool for histogram visualization

[3.0.0]: https://github.com/caseywdunn/sharkmer/compare/v2.0.0...v3.0.0
[2.0.0]: https://github.com/caseywdunn/sharkmer/compare/v1.0.1...v2.0.0
[1.0.1]: https://github.com/caseywdunn/sharkmer/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/caseywdunn/sharkmer/releases/tag/v1.0.0

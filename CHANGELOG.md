# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - v2.0.0

### Added

- **`--sra` flag**: Stream reads directly from ENA by accession — no SRA
  toolkit required (#65). Sample name is auto-derived from ENA metadata when
  `--sample` is omitted (#67).
- **Native gzip support**: `.fastq.gz` files are detected and decompressed
  automatically (#23).
- **`--pcr-panel`**: Select built-in primer panels by name (replaces `--pcr`
  for panel selection). Repeatable for multiple panels (#24).
- **`--pcr-file`**: Sideload primer panels from user-supplied YAML files (#25).
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
- **Shell completions**: `--generate-completions` for bash, zsh, fish,
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

- **CLI restructured**: `--pcr` split into `--pcr-panel`, `--pcr-file`,
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
- Updated README with v2.0 CLI examples, `--sra` usage, and gzip support (#66).

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

- `--pcr` flag (replaced by `--pcr-panel`, `--pcr-file`, `--pcr-primers`)
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

[Unreleased]: https://github.com/caseywdunn/sharkmer/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/caseywdunn/sharkmer/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/caseywdunn/sharkmer/releases/tag/v1.0.0

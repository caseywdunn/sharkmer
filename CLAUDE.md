# CLAUDE.md

Development context for AI-assisted work on sharkmer.

## What is sharkmer

A Rust CLI tool for kmer counting and in silico PCR (sPCR). Given whole genome
shotgun FASTQ reads, it can:

1. **sPCR**: Extract specific genomic regions using primer pairs without full
   genome assembly — produces FASTA amplicon files
2. **Incremental kmer counting**: Build kmer spectra across progressively
   larger data subsets for genome size estimation

## Repository structure

```
sharkmer/                  # Repo root
├── Cargo.toml             # Dependencies, feature flags, version
├── CLAUDE.md              # This file
├── PLAN.md                # v2.0 phased execution order with checkboxes
├── ROADMAP.md             # Multi-version release plan with issue references
├── CHANGELOG.md           # Release history
├── CONTRIBUTING.md        # Branching model, quality gates, bioconda recipe
├── readme.md              # User documentation
├── meta.yaml              # Bioconda recipe
├── build.sh               # Bioconda build script
├── src/
│   ├── main.rs            # CLI entry point, FASTQ ingestion, orchestration
│   ├── kmer/
│   │   └── mod.rs         # 2-bit encoding, kmer extraction, counting, histograms
│   └── pcr/
│       ├── mod.rs         # De Bruijn graph construction, extension, path finding
│       └── preconfigured.rs  # YAML panel loading (built-in via include_str!, user via --pcr-panel-file)
├── panels/                # Built-in primer panel YAML files (7 panels)
├── tests/
│   ├── fixtures/          # ERR571460 100k reads (gzipped) for integration tests
│   └── spcr_18s.rs        # Integration test: 18S recovery from ERR571460
├── sharkmer_viewer/       # Python tool for histogram visualization
└── benchmarks/            # Regression benchmark suite
```

## Build and test

```bash
cargo test --release     # Run tests (release mode for integration test speed)
cargo clippy             # Lint
cargo fmt --check        # Format check
cargo build --release    # Release build
```

All commands run from the repo root.

## Feature flags

Two hash map backends, selected at compile time:

- `ahashmap` (default): `ahash` — fast hash with AES-NI hardware acceleration
- `fxhashmap`: `rustc_hash::FxHashMap` — fast non-cryptographic hash

```bash
cargo build --features intmap --no-default-features
```

## Module architecture

### kmer/mod.rs (~940 lines)

- `Read`: 2-bit encoded DNA sequence (A=00, C=01, G=10, T=11)
- `Read::get_kmers()`: Sliding window extracts canonical kmers (min of
  forward/revcomp) as u64
- `KmerCounts`: Hash map wrapper (polymorphic via feature flags) storing
  kmer→count. Uses `KmerMap` trait to abstract over map implementations.
- `Chunk`: Groups reads for incremental processing
- `Histogram`: Kmer count frequency distribution

### pcr/mod.rs (~2300 lines)

Pipeline: preprocess primers → find primer kmers in data → seed graph →
extend graph → prune → find paths → generate sequences → deduplicate

Key data structures:
- `DBNode`: (k-1)-mer graph node with start/end/terminal flags
- `DBEdge`: Kmer edge with observed count
- `PCRParams`: Primer pair configuration
- Graph: `petgraph::StableDiGraph<DBNode, DBEdge>`

Key constants (may need tuning):
- `COVERAGE_MULTIPLIER = 2`: High coverage definition
- `COVERAGE_STEPS = 4`: Threshold reduction steps
- `MAX_NUM_PRIMER_KMERS = 100`: Primer variant filtering cap
- `MAX_NUM_NODES = 50_000`: Graph size limit
- `MAX_NUM_PATHS_PER_PAIR = 20`: Path enumeration limit
- `MAX_NUM_AMPLICONS = 20`: Output sequence limit
- `DEFAULT_DEDUP_EDIT_THRESHOLD = 10`: Levenshtein threshold for deduplication (configurable per primer pair)
- `BALLOONING_COUNT_THRESHOLD_MULTIPLIER = 10.0`: Spurious edge detection

### pcr/preconfigured.rs (~100 lines)

Loads primer panels from YAML files via `include_str!()` (built-in panels)
or from user-supplied files via `--pcr-panel-file`. Seven built-in panels:
cnidaria, human, teleostei, angiospermae, insecta, bacteria, metazoa.

### panels/ directory

YAML files defining built-in primer panels. Each file specifies panel name,
description, and a list of primers with PCRParams fields. Embedded at
compile time via `include_str!()`.

## Current known issues

None currently tracked.

## Git workflow

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full branching model, quality
gates, and patching workflow. Key points:

- Work on issue branches from `dev`, merge back to `dev` before starting
  the next issue.
- Quality gates before merge: `cargo test`, `cargo clippy` (no warnings),
  `cargo fmt --check`.
- Present a summary and any concerns to the user for review before committing.
  Check off the issue in PLAN.md when done.

## Current development

Version is `2.0.0` on `dev` branch. See PLAN.md for the phased
execution order. 45 issues across 7 phases. Run regression benchmarks
after each phase.

## Key design decisions for v2.0

- **CLI**: Split `--pcr` into `--pcr-panel`, `--pcr-panel-file`, `--pcr-primers`.
  Make `--sample` required. Default `--chunks 0` (skip histograms).
  Add `--list-panels`, `--export-panel`, `--help-pcr`.
  (implemented in Phase 4)
- **Output**: YAML for stats (`.stats.yaml`) with PCR results, histogram
  files get header rows and comment lines with version/params. Logs to
  stderr, data to stdout. (stats YAML and histogram headers implemented
  in Phase 4; FASTA header format pending #34)
- **Logging**: `log` + `env_logger`, `-v`/`-vv`/`-vvv` replaces `--verbosity N`
  (implemented in Phase 3)
- **Primer panels**: YAML files in `panels/` embedded via `include_str!()`.
  Shared parser for built-in and user-sideloaded panels via `--pcr-panel-file`.
  (implemented in Phase 3)
- **Error handling**: `anyhow` for error propagation (chosen over `thiserror`
  since sharkmer is a binary, not a library), all panics removed from library
  code, meaningful exit codes.
- **Performance**: Eliminate PCR hash table copy (#28, #52), HashMap for node
  lookup (#27), parallelize sPCR across genes (#36), free chunks after
  merge (#53).

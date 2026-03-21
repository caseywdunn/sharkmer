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
├── CLAUDE.md              # This file
├── PLAN.md                # v2.0 phased execution order with checkboxes
├── ROADMAP.md             # Multi-version release plan with issue references
├── CHANGELOG.md           # Release history
├── CONTRIBUTING.md        # Branching model, quality gates, bioconda recipe
├── readme.md              # User documentation
├── sharkmer/              # Rust crate (binary lives here)
│   ├── Cargo.toml         # Dependencies, feature flags, version
│   ├── meta.yaml          # Bioconda recipe
│   ├── build.sh           # Bioconda build script
│   ├── data/              # Test data (SRR5324768, gunzip before use)
│   └── src/
│       ├── main.rs        # CLI entry point, FASTQ ingestion, orchestration
│       ├── test.rs         # Integration test placeholder (mostly empty)
│       ├── kmer/
│       │   └── mod.rs     # 2-bit encoding, kmer extraction, counting, histograms
│       └── pcr/
│           ├── mod.rs     # De Bruijn graph construction, extension, path finding
│           └── preconfigured.rs  # Primer panel definitions (to be replaced with YAML)
├── sharkmer_viewer/       # Python tool for histogram visualization
├── docker/                # Dockerfile for development
├── tests/                 # Snakemake workflow for testing against SRA datasets
│   ├── Snakefile
│   ├── config.yaml        # 30+ test samples with SRA accessions
│   ├── sra.py             # SRA download utilities
│   └── aggregate_stats.py
└── genomescopemovie.sh    # Shell script for incremental visualization
```

## Build and test

```bash
cd sharkmer           # The Rust crate subdirectory
cargo test            # Run tests (currently broken — #17)
cargo clippy          # Lint
cargo fmt --check     # Format check
cargo build --release # Release build
```

Note: `Cargo.toml` is in `sharkmer/` not the repo root. Issue #20 will add a
workspace Cargo.toml at the root.

## Feature flags

Three hash map backends, selected at compile time:

- `fxhashmap` (default): `rustc_hash::FxHashMap` — fast non-cryptographic hash
- `intmap`: `intmap::IntMap` — optimized for integer keys
- `nohashmap`: Identity hash via std `HashMap` — for benchmarking

```bash
cargo build --features intmap --no-default-features
```

## Module architecture

### kmer/mod.rs (~940 lines)

- `Read`: 2-bit encoded DNA sequence (A=00, C=01, G=10, T=11)
- `Read::get_kmers()`: Sliding window extracts canonical kmers (min of
  forward/revcomp) as u64
- `KmerCounts`: Hash map wrapper (polymorphic via feature flags) storing
  kmer→count. Three nearly identical impl blocks — #26 will deduplicate.
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
- `DISTANCE_EDIT_THRESHOLD = 10`: Levenshtein threshold for deduplication
- `BALLOONING_COUNT_THRESHOLD_MULTIPLIER = 10.0`: Spurious edge detection

### pcr/preconfigured.rs (~1100 lines)

Hardcoded primer panels (cnidaria, human, teleostei, angiospermae, insecta,
bacteria, metazoa). Issue #33 will replace with YAML files via
`include_str!()`.

## Current known issues

- **Tests don't compile**: `pcr/mod.rs:2238` calls nonexistent
  `add_reverse_complements()`, line 2291 has wrong arity for
  `filter_primer_kmers()`. See #17.
- **Panics in library code**: ~30 `.unwrap()` calls and ~10 `panic!()` in
  library code. See #21.
- **No CI**: No `.github/workflows/`. See #19.
- **Unsigned underflow risks**: Several `usize` subtractions can wrap. See #55.
- **Deduplication may drop real variants**: Edit distance threshold is hardcoded
  and comparison is only against the first record. See #56.

## Branching model

- `master` — tagged releases only
- `dev` — active development for the next unreleased version
- `vN` (e.g. `v2`, `v3`) — release maintenance branches, created from the
  release tag when a major version ships. Used for patches only.
- Issue branches → `dev` (for new features) or `vN` (for patches)
- To patch a released version while working on the next: branch from `vN`,
  fix, merge to `vN`, tag, merge to `master`, cherry-pick into `dev`
- Quality gates: `cargo test` passes, `cargo clippy` no warnings,
  `cargo fmt --check` clean

## Current development

Version is `2.0.0-alpha` on `dev` branch. See PLAN.md for the phased
execution order. 45 issues across 7 phases. Run regression benchmarks
after each phase.

## Issue workflow

When working on an issue, complete a wrap-up phase before committing:

1. **Summarize changes**: List what was modified and why, including any
   non-obvious decisions made during implementation.
2. **Flag concerns**: Note anything that came up during the work that may need
   to be addressed in other issues — e.g. related code that looks fragile,
   assumptions that may not hold, or scope that was intentionally deferred.
   Open or reference issues for these as appropriate.
3. **User review**: Present the summary and concerns to the user for review
   before committing. Do not commit or close the issue until the user confirms.

## Key design decisions for v2.0

- **CLI**: Split `--pcr` into `--pcr-panel`, `--pcr-file`, `--pcr-primers`.
  Make `--sample` required. Default `--chunks 0` (skip histograms).
- **Output**: YAML for stats (`.stats.yaml`), FASTA headers with key=value
  metadata, histogram files get header rows. Logs to stderr, data to stdout.
- **Logging**: `log` + `env_logger`, `-v`/`-vv`/`-vvv` replaces `--verbosity N`
- **Primer panels**: YAML files in `panels/` embedded via `include_str!()`.
  Shared parser for built-in and user-sideloaded panels.
- **Error handling**: `anyhow` or `thiserror`, all panics removed from library
  code, meaningful exit codes.
- **Performance**: Eliminate PCR hash table copy (#28, #52), HashMap for node
  lookup (#27), parallelize sPCR across genes (#36), free chunks after
  merge (#53).

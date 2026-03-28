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
├── dev_docs/              # Development planning documents
│   ├── PLAN.md            # v3.0 phased execution order with checkboxes
│   └── DESIGN_DECISIONS.md # Background analysis for design decisions
├── ROADMAP.md             # Multi-version release plan with issue references
├── CHANGELOG.md           # Release history
├── CONTRIBUTING.md        # Branching model, quality gates, bioconda recipe
├── README.md              # User documentation
├── meta.yaml              # Bioconda recipe
├── build.sh               # Bioconda build script
├── src/
│   ├── main.rs            # Entry point, main() orchestration (~115 lines)
│   ├── cli.rs             # Args struct, CLI parsing, validation, early exits
│   ├── io.rs              # FASTQ reading, ENA streaming, FASTA writing
│   ├── format.rs          # Number/byte/duration formatting utilities
│   ├── stats.rs           # RunStats, PCR execution, YAML stats output
│   ├── kmer/
│   │   ├── mod.rs         # Re-exports and tests
│   │   ├── encoding.rs    # 2-bit encoding, Read struct, kmer extraction
│   │   ├── counting.rs    # KmerMap trait, KmerCounts, FilteredKmerCounts
│   │   ├── histogram.rs   # Histogram struct
│   │   └── chunk.rs       # Chunk struct
│   └── pcr/
│       ├── mod.rs         # do_pcr() orchestration, PCRParams, validation
│       ├── primers.rs     # Primer preprocessing, ambiguity, mismatch permutation
│       ├── graph.rs       # De Bruijn graph construction, extension, diagnostics
│       ├── pruning.rs     # Graph pruning (balloons, side branches, orphans)
│       ├── paths.rs       # Path finding, sequence extraction, deduplication
│       └── preconfigured.rs  # YAML panel loading (built-in + sideloaded)
├── panels/                # Built-in primer panel YAML files (7 panels)
├── tests/
│   ├── fixtures/          # ERR571460 100k reads (gzipped) for integration tests
│   └── spcr_18s.rs        # Integration test: 18S recovery from ERR571460
├── sharkmer_viewer/       # Python tool for histogram visualization
└── benchmarks/            # Regression benchmark suite
```

## Build and test

These commands match the CI quality gates exactly. All must pass before
committing:

```bash
cargo fmt --all --check        # Format check (CI-exact)
cargo clippy -- -D warnings    # Lint — warnings are errors (CI-exact)
cargo test --release           # Run tests (release mode for integration test speed)
cargo build --release          # Release build
```

**Important:** Use `cargo clippy -- -D warnings` (not bare `cargo clippy`).
CI treats warnings as errors. If clippy passes locally without `-D warnings`
but fails in CI, you have a warning that needs fixing.

All commands run from the repo root.

### Pre-commit hook

A git pre-commit hook enforces these same CI gates locally. Install it:

```bash
ln -sf ../../scripts/pre-commit .git/hooks/pre-commit
```

This prevents commits that would fail CI. The hook runs fmt, clippy
(`-D warnings`), and tests before every commit.

## Feature flags

Two hash map backends, selected at compile time:

- `ahashmap` (default): `ahash` — fast hash with AES-NI hardware acceleration
- `fxhashmap`: `rustc_hash::FxHashMap` — fast non-cryptographic hash

```bash
cargo build --features fxhashmap --no-default-features
```

## Module architecture

### main.rs (~135 lines)

Entry point. Parses CLI, delegates to helper modules, orchestrates the
pipeline: init logging → collect primers → validate → ingest reads →
consolidate → run PCR → write stats → print summary.

### cli.rs (~715 lines)

`Args` struct (clap), `ColorMode`, `parse_pcr_primers_string()`,
`init_logging()`, early exits, validation, `--dry-run`.

### io.rs (~640 lines)

`read_fastq()`, `validate_fastq_record()`, `get_ena_fastq_urls()`,
`write_fasta_record()`, `ingest_reads()`, `consolidate_and_histogram()`.

### format.rs (~45 lines)

`format_count()`, `format_bytes()`, `format_duration()`.

### stats.rs (~215 lines)

`RunStats`, `PcrGeneResult`, `run_pcr()`, `write_stats()`, `print_summary()`.

### kmer/ (4 submodules)

- `encoding.rs`: `Read` struct, 2-bit encoding (A=00, C=01, G=10, T=11),
  `get_kmers()`, `revcomp_kmer()`, `seq_to_reads()`, `kmer_to_seq()`
- `counting.rs`: `KmerMap` trait, feature-gated impls, `KmerCounts`,
  `FilteredKmerCounts`
- `histogram.rs`: `Histogram` struct
- `chunk.rs`: `Chunk` struct
- `mod.rs`: re-exports and tests

### pcr/ (5 submodules)

Pipeline: preprocess primers → find primer kmers → seed graph →
extend graph → prune → find paths → generate sequences → deduplicate

- `primers.rs`: Primer preprocessing, ambiguity resolution, mismatch
  permutation, kmer extraction
- `graph.rs`: `DBNode`, `DBEdge`, seed graph, `extend_graph()`,
  diagnostics. Graph: `petgraph::StableDiGraph<DBNode, DBEdge>`
- `pruning.rs`: `pop_balloons()`, `remove_side_branches()`,
  `remove_orphan_nodes()`
- `paths.rs`: `get_assembly_paths()`, sequence extraction, deduplication
- `mod.rs`: `do_pcr()` orchestration, `PCRParams`, validation, constants
- `preconfigured.rs`: YAML panel loading (built-in via `include_str!()`,
  user via `--pcr-panel-file`)

Key constants (may need tuning):
- `COVERAGE_MULTIPLIER = 2`: High coverage definition
- `COVERAGE_STEPS = 4`: Threshold reduction steps
- `MAX_NUM_PRIMER_KMERS = 100`: Primer variant filtering cap
- `MAX_NUM_NODES = 50_000`: Graph size limit
- `MAX_NUM_PATHS_PER_PAIR = 20`: Path enumeration limit
- `MAX_NUM_AMPLICONS = 20`: Output sequence limit
- `DEFAULT_DEDUP_EDIT_THRESHOLD = 10`: Levenshtein threshold for dedup
- `BALLOONING_COUNT_THRESHOLD_MULTIPLIER = 10.0`: Spurious edge detection

### panels/ directory

YAML files defining built-in primer panels (7 panels). Embedded at
compile time via `include_str!()`.

## Current known issues

None currently tracked.

## Git workflow

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full branching model, quality
gates, and patching workflow. Key points:

- Work on issue branches from `dev`, merge back to `dev` before starting
  the next issue.
- Quality gates before merge: `cargo fmt --all --check`,
  `cargo clippy -- -D warnings`, `cargo test --release`.
- The pre-commit hook (`scripts/pre-commit`) enforces these automatically.
- Present a summary and any concerns to the user for review before committing.
  Check off the issue in dev_docs/PLAN.md when done.

## Current development

Version is `3.0.0-dev` on `dev` branch. See dev_docs/PLAN.md for the phased
execution order and ROADMAP.md for scope and rationale.

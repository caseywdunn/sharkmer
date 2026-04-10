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
│   ├── DESIGN_DECISIONS.md # Background analysis for design decisions
│   └── OVERVIEW.md        # Architecture overview with data flow diagrams
├── ROADMAP.md             # Multi-version release plan with issue references
├── CHANGELOG.md           # Release history
├── CONTRIBUTING.md        # Branching model, quality gates, bioconda recipe
├── README.md              # User documentation
├── meta.yaml              # Bioconda recipe
├── build.sh               # Bioconda build script
├── src/
│   ├── main.rs            # Entry point, main() orchestration (~145 lines)
│   ├── cache.rs           # Read cache for remote downloads (SHA-256 verified)
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

### Running benchmarks

Benchmarks require the `sharkmer-bench` conda environment, which provides
`blastn` for BLAST validation of amplicons:

```bash
conda env create -f benchmarks/environment.yaml   # first time only
conda activate sharkmer-bench
python benchmarks/run_benchmark.py --samples Porites_lutea --max-reads 1000000 --no-blast
python benchmarks/run_benchmark.py                # all samples, with BLAST
```

Benchmarks use local FASTQ files from `benchmarks/data/` when present.
If local files are absent, they fall back to `--ena` with read caching
(cached in `benchmarks/data/cache/`). A local BLAST database in `/db/`
is used if `blastn` is available; otherwise falls back to NCBI remote API.

## Feature flags

Two hash map backends, selected at compile time:

- `ahashmap` (default): `ahash` — fast hash with AES-NI hardware acceleration
- `fxhashmap`: `rustc_hash::FxHashMap` — fast non-cryptographic hash

```bash
cargo build --features fxhashmap --no-default-features
```

## Module architecture

See [dev_docs/OVERVIEW.md](dev_docs/OVERVIEW.md) for a detailed
architecture overview with data flow diagrams.

### main.rs (~215 lines)

Entry point. Parses CLI, delegates to helper modules, orchestrates the
pipeline: init logging → collect primers → validate → pre-encode primer
Oligos (if `--read-eval`) → build Oligo filter → ingest reads (Pass 1)
→ consolidate → re-read sequences (if `--read-threading`, Pass 2) →
run PCR → write stats → print summary.

### cli.rs (~770 lines)

`Args` struct (clap), `ColorMode`, `parse_pcr_primers_string()`,
`init_logging()`, early exits, validation, `--dry-run`. New flags:
`--read-eval` (Pass 1 read retention for seed eval), `--read-threading`
(Pass 2 re-read for graph annotation), `--max-nodes` (hidden, graph
node budget). `--paired` exists but is hidden — paired-end R1/R2
reading works, but the paired-end phasing feature it is meant to
unlock is not yet wired into branch ranking; see issue #101.

### io.rs (~1260 lines)

`read_fastq()`, `validate_fastq_record()`, `get_ena_fastq_urls()`,
`write_fasta_record()`, `ingest_reads()`, `consolidate_and_histogram()`.
New: `OligoFilter` (bloom filter + AHashSet for Pass 1 read retention),
`RetainedRead`/`RetainedReads`, `ReadRecord`/`ReadPlan`/`ReadSourcePlan`
(Pass 2 re-reading), `reread_sequences()`, `read_fastq_paired()`.

### format.rs (~45 lines)

`format_count()`, `format_bytes()`, `format_duration()`.

### stats.rs (~235 lines)

`RunStats`, `PcrGeneResult`, `run_pcr()`, `write_stats()`, `print_summary()`.

### kmer/ (4 submodules)

- `encoding.rs`: `Read` struct, 2-bit encoding (A=00, C=01, G=10, T=11),
  `get_kmers()`, `revcomp_kmer()` (byte LUT), `seq_to_reads()`,
  `kmer_to_seq()`, `kmer_last_base()`
- `counting.rs`: `KmerMap` trait (with `insert_or_add_get_old`),
  feature-gated impls, `KmerCounts`, `FilteredKmerCounts`
- `histogram.rs`: `Histogram` struct
- `chunk.rs`: `Chunk` struct
- `mod.rs`: re-exports and tests

### pcr/ (9 submodules)

Pipeline: preprocess primers → find primer kmers → seed graph →
extend graph → prune → thread reads (if `--read-threading`) →
resolve bubbles → find paths → generate sequences → deduplicate

- `primers.rs`: Primer preprocessing, ambiguity resolution, mismatch
  permutation, kmer extraction.
- `graph.rs`: `DBNode`, `DBEdge`, seed graph, unified bidirectional
  `extend_graph()`, diagnostics.
  Graph: `petgraph::StableDiGraph<DBNode, DBEdge>`
- `pruning.rs`: `remove_low_coverage_tips()`, `reachability_pruning()`
- `paths.rs`: `get_assembly_paths()`, sequence extraction, deduplication.
  Edge ordering uses bubble resolution preferences when available.
- `threading.rs`: `thread_reads()`, `thread_reads_paired()`. Maps reads
  to graph edges via maximal contiguous runs. Graph-agnostic API.
  `EdgeReadSupport`, `BranchLink`, `ThreadingAnnotations`. Note:
  `PairedEndLink` and `thread_reads_paired` construct paired-end link
  data, but nothing downstream consumes it yet — see issue #101 for
  what needs to happen to complete paired-end phasing.
- `read_filter.rs`: `PrimerReadFilter` for per-gene read filtering
  during read threading.
- `bubble.rs`: `resolve_bubbles()`. Detects simple bubbles, ranks
  branches by read support + phasing, returns edge preferences.
- `mod.rs`: `do_pcr()` orchestration, `PCRParams`, `PathScore`
  (with read-support fields), validation, constants
- `preconfigured.rs`: YAML panel loading (built-in via `include_str!()`,
  user via `--pcr-panel-file`)

Key constants (all exposed as hidden CLI arguments via PCRParams unless
noted):
- `COVERAGE_MULTIPLIER = 2`: High coverage definition
- `COVERAGE_STEPS = 4`: Threshold reduction steps
- `DEFAULT_MAX_NUM_NODES = 500_000`: Graph size limit (computed from
  data volume by `compute_node_budget`, no CLI override)
- `DEFAULT_MAX_DFS_STATES = 100_000`: DFS state budget (`--max-dfs-states`)
- `DEFAULT_MAX_PATHS_PER_PAIR = 20`: Path enumeration limit (`--max-paths-per-pair`)
- `DEFAULT_MAX_NODE_VISITS = 2`: Cycle tolerance (`--max-node-visits`)
- `DEFAULT_MAX_NUM_PRIMER_KMERS = 20`: Primer variant cap (`--max-primer-kmers`)
- `DEFAULT_HIGH_COVERAGE_RATIO = 10.0`: Repeat edge filter (`--high-coverage-ratio`)
- `DEFAULT_TIP_COVERAGE_FRACTION = 0.1`: Tip pruning (`--tip-coverage-fraction`)
- `MAX_NUM_AMPLICONS = 20`: Output sequence limit
- `DEFAULT_DEDUP_EDIT_THRESHOLD = 10`: Levenshtein threshold for dedup

### panels/ directory

YAML files defining built-in primer panels (7 panels). Embedded at
compile time via `include_str!()`.

## Current known issues

- **Seed eval threshold (fixed in Phase 8)**: Changed from max to median
  primer kmer count. Degenerate primers with off-target matches previously
  inflated the threshold too high. Needs benchmark validation with
  Drosophila 16S/28S to confirm the fix works in practice.

## Git workflow

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full branching model, quality
gates, and patching workflow. Key points:

- Work on issue branches from `dev`, merge back to `dev` before starting
  the next issue.
- Quality gates before merge: `cargo fmt --all --check`,
  `cargo clippy -- -D warnings`, `cargo test --release`.
- The pre-commit hook (`scripts/pre-commit`) enforces these automatically.
- Present a summary and any concerns to the user for review before committing.
- **Before committing**, check off completed items in dev_docs/PLAN.md.
  This is easy to forget — do it as part of the commit, not as a follow-up.

## Current development

Version is `3.0.0-dev` on `dev` branch. See dev_docs/PLAN.md for the phased
execution order and ROADMAP.md for scope and rationale.

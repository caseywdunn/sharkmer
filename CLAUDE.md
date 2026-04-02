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
│   └── overview.md        # Architecture overview with data flow diagrams
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

See [dev_docs/overview.md](dev_docs/overview.md) for a detailed
architecture overview with data flow diagrams.

### main.rs (~180 lines)

Entry point. Parses CLI, delegates to helper modules, orchestrates the
pipeline: init logging → collect primers → validate → pre-encode primer
Oligos (if `--read-eval`) → build Oligo filter → ingest reads (Pass 1)
→ consolidate → re-read sequences (if `--read-threading`, Pass 2) →
run PCR → write stats → print summary.

### cli.rs (~770 lines)

`Args` struct (clap), `ColorMode`, `parse_pcr_primers_string()`,
`init_logging()`, early exits, validation, `--dry-run`. New flags:
`--read-eval` (Pass 1 read retention for seed eval), `--read-threading`
(Pass 2 re-read for graph annotation), `--paired` (paired-end R1/R2),
`--max-nodes` (hidden, graph node budget).

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
  `get_kmers()`, `revcomp_kmer()`, `seq_to_reads()`, `kmer_to_seq()`
- `counting.rs`: `KmerMap` trait, feature-gated impls, `KmerCounts`,
  `FilteredKmerCounts`
- `histogram.rs`: `Histogram` struct
- `chunk.rs`: `Chunk` struct
- `mod.rs`: re-exports and tests

### pcr/ (8 submodules)

Pipeline: preprocess primers → find primer kmers → seed graph →
evaluate seeds (with optional read divergence) → extend graph →
prune → thread reads (if `--read-threading`) → resolve bubbles →
find paths → generate sequences → deduplicate

- `primers.rs`: Primer preprocessing, ambiguity resolution, mismatch
  permutation, kmer extraction. New: `PrimerOligoSet`,
  `preprocess_primer_oligos()` for Pass 1 Oligo encoding
- `graph.rs`: `DBNode`, `DBEdge`, seed graph, `extend_graph()`,
  `extend_graph_reverse()`, diagnostics.
  Graph: `petgraph::StableDiGraph<DBNode, DBEdge>`
- `pruning.rs`: `remove_low_coverage_tips()`, `reachability_pruning()`
- `paths.rs`: `get_assembly_paths()`, sequence extraction, deduplication.
  Edge ordering uses bubble resolution preferences when available.
- `seed_eval.rs`: `evaluate_seeds()`, `bounded_extend()`,
  `check_read_divergence()`. Bounded local exploration to filter
  off-target seeds before full graph extension.
- `threading.rs`: `thread_reads()`, `thread_reads_paired()`. Maps reads
  to graph edges via maximal contiguous runs. Graph-agnostic API.
  `EdgeReadSupport`, `BranchLink`, `PairedEndLink`, `ThreadingAnnotations`.
- `read_filter.rs`: `PrimerReadFilter` for per-gene read filtering
  during Pass 2 threading
- `bubble.rs`: `resolve_bubbles()`. Detects simple bubbles, ranks
  branches by read support + phasing, returns edge preferences.
- `mod.rs`: `do_pcr()` orchestration, `PCRParams`, `PathScore`
  (with read-support fields), validation, constants
- `preconfigured.rs`: YAML panel loading (built-in via `include_str!()`,
  user via `--pcr-panel-file`)

Key constants (may need tuning — see #109):
- `COVERAGE_MULTIPLIER = 2`: High coverage definition
- `COVERAGE_STEPS = 4`: Threshold reduction steps
- `MAX_NUM_PRIMER_KMERS = 100`: Primer variant filtering cap
- `DEFAULT_MAX_NUM_NODES = 50_000`: Graph size limit (CLI: `--max-nodes`)
- `MAX_NUM_PATHS_PER_PAIR = 20`: Path enumeration limit
- `MAX_NUM_AMPLICONS = 20`: Output sequence limit
- `DEFAULT_DEDUP_EDIT_THRESHOLD = 10`: Levenshtein threshold for dedup
- `HIGH_COVERAGE_RATIO_THRESHOLD = 10.0`: Repeat edge filtering

### panels/ directory

YAML files defining built-in primer panels (7 panels). Embedded at
compile time via `include_str!()`.

## Current known issues

- **Seed eval threshold too stringent for degenerate primers (#109)**:
  Perfect-match failures (Drosophila 16S, 28S) are caused by the seed
  evaluation threshold being derived from the max primer kmer count.
  With degenerate primers, off-target matches inflate this count, making
  the threshold too high for real seeds to extend. All seeds are
  abandoned at "1 node". This is NOT a node budget issue — increasing
  `--max-nodes` to 500K has no effect. Fix requires #109 parameter
  tuning (use median primer count, or step down seed eval threshold).
- **~15% runtime overhead from Phase 4-6 code** (compared to pre-Phase-4):
  The `f64` edge scoring in DFS path finding (bubble resolution
  infrastructure) adds overhead even when `--read-threading` is off.
  This is the cost of the pluggable scoring architecture.

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

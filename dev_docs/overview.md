# Sharkmer Architecture Overview

High-level overview of how sharkmer works, from FASTQ input to amplicon
output.

## Pipeline

```
 FASTQ files / ENA stream / stdin
              |
              v
    +--------------------+
    |   Read ingestion   |   io.rs: read_fastq()
    |  Parse & validate  |   Batches of 1000 reads
    +--------------------+
              |
              v
    +--------------------+
    |  Chunk distribution |   Round-robin across N chunks
    |  (kmer extraction)  |   encoding.rs: kmers_from_ascii()
    +--------------------+   counting.rs: KmerCounts (u64 -> u32)
              |
              v
    +--------------------+
    |   Consolidation    |   io.rs: consolidate_and_histogram()
    |  Merge all chunks  |   Incremental histogram updates
    |  into one table    |   Write .histo files (if chunks > 0)
    +--------------------+
              |
              v
    +--------------------+
    |  Filtered view     |   counting.rs: FilteredKmerCounts
    | (lazy, no copy)    |   Threshold checked at lookup time
    +--------------------+
              |
              v
    +--------------------+
    |   In silico PCR    |   stats.rs: run_pcr() via rayon
    |  (per gene, parallel)   pcr/mod.rs: do_pcr()
    +--------------------+
              |
              v
       FASTA amplicons
       .stats.yaml
       .histo files
```

## Modules

### main.rs

Entry point. Orchestrates the pipeline in order:

1. Parse CLI, init logging
2. Collect primer params from panels, files, and inline `--pcr-primers`
3. Validate arguments, handle early exits (`--list-panels`, `--dry-run`, etc.)
4. Ingest reads into chunks
5. Consolidate chunks, write histograms
6. Run in silico PCR
7. Write stats YAML, print summary

### cli.rs

`Args` struct (clap-derived). Handles input sources (files, `--ena`,
stdin), PCR configuration, k-mer parameters, output options. Validates
arguments and clamps PCR min-count to `--min-kmer-count`.

### io.rs

Two main functions:

- **`ingest_reads()`** -- creates chunks, reads FASTQ from all sources
  (ENA/files/stdin), distributes sequences to chunks in round-robin
  batches of 1000. Auto-detects gzip.
- **`consolidate_and_histogram()`** -- merges chunk kmer tables into one
  `KmerCounts`. If `chunks > 0`, writes incremental histograms and
  validates count consistency.

### kmer/

```
encoding.rs    2-bit DNA encoding, kmer extraction
  |
  +-- kmers_from_ascii(seq, k) -> Vec<u64>   [hot path]
  |     Single-pass extraction, splits on N, returns canonical kmers
  |
  +-- revcomp_kmer(), kmer_to_seq()          [utility]
  +-- Read struct                             [legacy, used in tests]

counting.rs    Hash table: u64 kmer -> u32 count
  |
  +-- KmerCounts          Insert, merge, iterate, filtered view
  +-- FilteredKmerCounts  Lazy view with min_count threshold

chunk.rs       Groups reads for incremental counting
  |
  +-- Chunk { kmer_counts, n_reads, n_bases }

histogram.rs   Frequency distribution of kmer counts
  |
  +-- Histogram { histo: Vec<u64>, histo_large: FxHashMap }
  +-- move_count()  -- atomic bin transition for incremental merging
```

**Kmer encoding**: Each base is 2 bits (A=00, C=01, G=10, T=11). A kmer
of length k fits in a `u64` (max k=31). Canonical form = min(forward,
reverse complement).

**Count type**: `u32` with `saturating_add`. Saves ~25% memory per hash
table entry vs `u64`. Counts exceeding 4 billion are capped at
`u32::MAX`. The `get_n_kmers()` method returns `u64` (sum can overflow
u32).

### pcr/

In silico PCR reconstructs amplicon sequences from kmer counts using
De Bruijn graphs.

#### Data structures

```
DBNode {
    sub_kmer: u64,      (k-1)-mer, the overlap between adjacent kmers
    is_start: bool,     forward primer binding site
    is_end: bool,       reverse primer binding site
    is_terminal: bool,  no further extension possible
    visited: bool,      has been processed during extension
}

DBEdge {
    _kmer: u64,         the full kmer connecting two nodes
    count: u32,         observed frequency in the data
}

Graph: petgraph::StableDiGraph<DBNode, DBEdge>
```

#### Pipeline (per gene)

```
 Primer sequences (may contain IUPAC ambiguity codes)
              |
              v
 +---------------------------+
 | 1. Primer preprocessing   |   primers.rs
 |    Resolve ambiguities    |   Expand R->AG, Y->CT, N->ACGT, etc.
 |    Generate mismatches    |   Up to `mismatches` substitutions
 |    Trim to `trim` length  |   Shorter oligos for binding search
 +---------------------------+
              |
              v
 +---------------------------+
 | 2. Find primer kmers      |   primers.rs: find_oligos_in_kmers()
 |    Scan kmer table for    |   Check both orientations
 |    kmers containing       |   Filter to top 100 by count
 |    primer binding sites   |
 +---------------------------+
              |
              v
 +---------------------------+
 | 3. Seed graph             |   graph.rs: create_seed_graph()
 |    Forward primer kmers   |   Each forward kmer iterated separately
 |      -> start nodes       |   (one graph per forward primer kmer)
 |    Reverse primer kmers   |
 |      -> end nodes         |
 +---------------------------+
              |
              v
 +---------------------------+
 | 4. Extend graph           |   graph.rs: extend_graph()
 |    For each unvisited     |   Try 4 possible extensions per node
 |    node, find extending   |   Check kmer exists above threshold
 |    kmers in the data      |   Mark terminal if no extensions
 +---------------------------+
 |    Repeated with decreasing coverage thresholds
 |    (high -> min_count, in COVERAGE_STEPS steps)
 |
 |    Ballooning detection:
 |      - Skip edges with count > 10x median
 |      - Pop high-degree subgraphs periodically
 |      - Terminate nodes in rapidly branching regions
 |      - Abandon if graph exceeds 50,000 nodes
              |
              v
 +---------------------------+
 | 5. Prune graph            |   pruning.rs
 |    pop_balloons()         |   Remove high-degree explosions
 |    remove_side_branches() |   Trim dead-end paths
 |    remove_orphan_nodes()  |   Delete disconnected components
 +---------------------------+
              |
              v
 +---------------------------+
 | 6. Find paths             |   paths.rs: get_assembly_paths()
 |    All simple paths from  |   Max 20 paths per start-end pair
 |    start to end nodes     |   Respect min/max length constraints
 +---------------------------+
              |
              v
 +---------------------------+
 | 7. Generate sequences     |   paths.rs: generate_sequences_from_paths()
 |    Chain node sub_kmers   |   Record edge count statistics
 |    Deduplicate by edit    |   Levenshtein threshold (default 10)
 |    distance               |   Max 20 amplicons per gene
 +---------------------------+
              |
              v
       FASTA records with kmer count statistics in headers
```

#### Coverage threshold stepping

The graph is extended multiple times with decreasing minimum kmer count
thresholds. Starting from `primer_count / 2` and stepping down to
`min_count` in 4 steps. If a complete path (start -> end) is found at a
higher threshold, the search stops. This allows clean assembly at high
coverage before falling back to noisier low-coverage extension.

#### Primer panels

Preconfigured panels are YAML files in `panels/`, embedded at compile
time via `include_str!()`. Users can also load panels from local files
(`--pcr-panel-file`) or URLs. Each panel defines a set of primer pairs
with gene names, sequences, length constraints, and mismatch tolerance.

### stats.rs

- **`run_pcr()`** -- runs `do_pcr()` in parallel (rayon) across all
  primer pairs. Writes per-gene FASTA files sequentially for
  deterministic output.
- **`write_stats()`** -- serializes `RunStats` to YAML.
- **`print_summary()`** -- human-readable table of results.

### format.rs

Small utility module: `format_count()`, `format_bytes()`,
`format_duration()`.

## Key design choices

- **Lazy filtering**: `FilteredKmerCounts` wraps `KmerCounts` without
  copying. The min_count threshold is checked at each lookup. This
  allows the full kmer table to be shared across parallel PCR runs
  while each gene uses its own effective threshold.

- **Canonical kmers**: Only the lexicographically smaller of a kmer and
  its reverse complement is stored. Lookups check both orientations.
  This halves memory usage.

- **Incremental histograms**: When `chunks > 0`, histograms are updated
  incrementally during chunk merging using `move_count()` rather than
  rebuilding from scratch. This enables genome size estimation from
  progressively larger data subsets.

- **Deterministic output**: Forward primer kmers are sorted before
  iteration. FASTA files are written sequentially (not in parallel).
  Products are re-numbered after deduplication.

- **Separate graphs per forward primer kmer**: Currently, each forward
  primer kmer seeds its own independent graph. Phase 3 (#94) plans to
  unify these into a single graph per gene.

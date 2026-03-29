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
    _kmer: u64,           the full kmer connecting two nodes
    count: u32,           observed frequency in the data
    coverage_ratio: f64,  count / global median (annotated post-pruning)
}

PathScore {
    median_count: f64,        robust central tendency of edge counts
    coverage_cv: f64,         coefficient of variation across edge counts
    max_coverage_ratio: f64,  highest edge coverage_ratio on the path
}
    composite() = median * cv_penalty * repeat_penalty

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
 |    All forward primer     |   Single graph per gene seeded with
 |    kmers -> start nodes   |   all forward kmers simultaneously
 |    All reverse primer     |
 |    kmers -> end nodes     |
 +---------------------------+
              |
              v
 +---------------------------+
 | 4. Extend graph           |   graph.rs: extend_graph()
 |    For each unvisited     |   Try 4 possible extensions per node
 |    node, find extending   |   Check kmer exists above threshold
 |    kmers in the data      |   Mark terminal if no extensions
 +---------------------------+
 |    Extended incrementally across coverage threshold steps
 |    (high -> min_count, in COVERAGE_STEPS steps)
 |    prepare_for_lower_threshold() resets terminal flags
 |    between steps; only newly qualifying edges are added
 |
 |    Coverage-ratio filtering:
 |      - Skip edges with count > 10x median
 |      - Abandon if graph exceeds 50,000 nodes
 |      - Break on first threshold that produces a product
              |
              v
 +---------------------------+
 | 5. Prune graph            |   pruning.rs (post-extension, on clone)
 |  remove_low_coverage_tips |   Dead-end tips < k with count
 |                           |     < 0.1x global median
 |  reachability_pruning()   |   Bidirectional BFS from start/end
 |                           |     nodes; remove unreachable nodes
 +---------------------------+
              |
              v
 +---------------------------+
 | 5b. Annotate edges        |   graph.rs: annotate_coverage_ratios()
 |    coverage_ratio =       |   count / global median per edge
 |    count / global median  |   Feeds into path scoring
 +---------------------------+
              |
              v
 +---------------------------+
 | 6. Find paths             |   paths.rs: get_assembly_paths()
 |    Coverage-weighted DFS  |   Explores highest-count edges first
 |    from start to end      |   Max 20 paths per start node
 |    nodes                  |   Bounded revisitation: MAX_NODE_VISITS=2
 +---------------------------+
              |
              v
 +---------------------------+
 | 7. Generate sequences     |   paths.rs: generate_sequences_from_paths()
 |    Chain node sub_kmers   |   PathScore: median, CV, coverage ratio
 |    Deduplicate by edit    |   Greedy clustering (O(N) memory)
 |    distance               |   Max 20 amplicons per gene
 +---------------------------+
              |
              v
       FASTA records with kmer count statistics in headers
```

#### Coverage threshold stepping

The graph is extended incrementally with decreasing minimum kmer count
thresholds. Starting from `primer_count / 2` and stepping down to
`min_count` in 4 steps (`COVERAGE_STEPS`). Unlike earlier versions that
rebuilt the graph from scratch at each threshold, Phase 3 extends the
existing graph: `prepare_for_lower_threshold()` resets terminal flags on
nodes that may gain new successors, then `extend_graph()` adds only
newly qualifying edges. If a complete path (start -> end) is found at a
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

- **Single graph per gene**: All forward primer kmers seed a single
  graph simultaneously (Phase 3, #94). This avoids redundant graph
  construction and ensures shared structure is discovered once.

- **Annotation-only model**: Graph structure (bubbles, variants) is
  preserved through pruning. Only dead-end tips with low coverage and
  unreachable nodes are removed. Path selection uses coverage-based
  scoring (`PathScore` composite) rather than destructive pruning to
  choose the best sequences. This preserves information for future
  read-threading phases.

- **Coverage-weighted path finding**: DFS explores highest-count edges
  first, so the first paths found are highest quality. Bounded
  revisitation (`MAX_NODE_VISITS = 2`) handles tandem repeats without
  exponential blowup. Replaces exhaustive `all_simple_paths` enumeration.

- **O(N) deduplication**: Greedy clustering computes bounded Levenshtein
  distance on-the-fly against kept records only, replacing the previous
  O(N^2) pairwise distance matrix.

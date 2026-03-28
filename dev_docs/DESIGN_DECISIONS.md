# Design Decisions

Background analysis and reasoning behind architectural decisions for sharkmer
development. Each section records the considerations at the time a decision
was made, so future work can revisit assumptions if circumstances change.

## Single-pass vs two-pass read ingestion for read threading

### Context

Read threading (v3.0 Phases 4-6) requires reads to be available after the
de Bruijn graph is constructed. Two approaches:

- **Single-pass**: Retain reads in memory during kmer counting, then thread
  them through the graph immediately.
- **Two-pass**: Count kmers and build the graph in Pass 1, then re-read the
  FASTQ in Pass 2 to thread reads through the completed graph.

### Analysis: coverage requirements

To reliably recover a single-copy nuclear gene (e.g., EF1) as a ~3 kb
amplicon from WGS data, we need enough reads that the amplicon region has
adequate local coverage for de Bruijn graph construction.

For a 3 kb amplicon with k=21:
- ~2,980 distinct kmers in the amplicon
- At local coverage C, expected count per kmer ~ C x (L-k+1)/L ~ 0.87C
  (for L=150 bp reads)
- For a connected graph, most kmers must be observed at least once
- P(kmer observed at least once) = 1 - e^(-0.87C)
- For 95% of kmers observed: C ~ 3.4x local coverage
- For 99% of kmers observed: C ~ 5.3x local coverage
- Practical minimum for graph construction with pruning: ~5-10x local
  coverage

### Analysis: read count and RAM for different genomes

Total reads needed for Cx whole-genome coverage, assuming 150 bp reads:

N = (C x G) / L

| Organism              | Genome size | 5x reads   | 10x reads  | 20x reads  |
|-----------------------|-------------|------------|------------|------------|
| Drosophila            | 140 Mb      | 4.7M       | 9.3M       | 18.7M      |
| Porites (coral)       | 542 Mb      | 18.1M      | 36.1M      | 72.3M      |
| Human                 | 3.1 Gb      | 103M       | 207M       | 414M       |
| Liriodendron (tulip)  | 1.7 Gb      | 57M        | 113M       | 227M       |
| Wheat                 | 17 Gb       | 567M       | 1,133M     | 2,267M     |

RAM for retaining all reads in memory (compact 2-bit encoding, ~40 bytes
per 150 bp read including overhead):

| Organism              | 5x RAM  | 10x RAM | 20x RAM |
|-----------------------|---------|---------|---------|
| Drosophila            | 0.2 GB  | 0.4 GB  | 0.7 GB  |
| Porites (coral)       | 0.7 GB  | 1.4 GB  | 2.9 GB  |
| Human                 | 4.1 GB  | 8.3 GB  | 16.6 GB |
| Liriodendron (tulip)  | 2.3 GB  | 4.5 GB  | 9.1 GB  |
| Wheat                 | 22.7 GB | 45.3 GB | 90.7 GB |

A consumer machine with 16 GB RAM can hold reads for a human genome at ~10x
but not 20x, and this is before accounting for the kmer hash table which is
also in memory. Larger genomes (conifers, wheat, salamanders) are infeasible
for single-pass at useful coverage levels.

### Analysis: selective retention during single pass

Could we retain only reads relevant to the amplicons during Pass 1? A 3 kb
amplicon is a tiny fraction of any genome. Expected reads hitting a 3 kb
target at 10x coverage of a 3.1 Gb genome:

- Effective target size: ~3 kb + 150 bp (read length) = ~3.15 kb
- Fraction of genome: 3.15 kb / 3.1 Gb ~ 10^-6
- Reads hitting target: 207M x 10^-6 ~ 210 reads per amplicon
- Even 10 amplicons: ~2,100 reads to retain — negligible RAM

The problem is a chicken-and-egg: which reads to keep depends on the graph,
but the graph isn't built until all reads are counted. Filtering by primer
kmers alone misses reads in the interior of the amplicon that don't overlap
any primer binding site.

### Analysis: two-pass RAM cost

Two-pass avoids retaining reads entirely:
- Pass 1: count kmers, build graphs (kmer hash table in memory, reads
  discarded after counting — this is the current architecture)
- Pass 2: stream reads one at a time, check each against graph kmers,
  thread matches. RAM cost is just the graph + kmer table (already in
  memory) plus one read buffer.

The second pass is essentially free in RAM. The cost is I/O: re-reading
files from disk, or reading from cache for remote sources.

### Analysis: how other assemblers handle this

Most de Bruijn graph assemblers use post-processing for read threading:
- **Velvet** and **ABySS** build graphs from kmers only, then use
  paired-end information in a separate scaffolding stage.
- **SPAdes** is the notable exception — it tracks read paths during graph
  construction, requiring substantially more memory. This design is
  motivated by single-cell sequencing where coverage is extremely uneven
  and standard coverage-based heuristics fail.

For sharkmer's use case (WGS at moderate coverage), the post-processing
approach aligns with the majority of assemblers.

### Decision

**Two-pass architecture.** The single-pass approach is infeasible on consumer
hardware for large genomes at useful coverage levels. Selective retention
during a single pass has a chicken-and-egg problem. Two-pass adds I/O cost
but has negligible RAM overhead and matches how most assemblers handle read
threading. The remote read cache (Phase 2) ensures that two-pass does not
double network traffic for remote sources.

If future profiling shows the second-pass I/O is a bottleneck, selective
retention could be revisited — e.g., retaining reads that contain any kmer
above a coverage threshold during Pass 1, accepting that some relevant reads
will be missed. But this optimization is not worth the complexity for v3.0.

## Single graph seeded with all forward primer kmers

### Context

Currently (`do_pcr()` in `pcr/mod.rs`), each forward primer kmer seeds its
own graph. The code iterates over forward primer kmers (generated by
ambiguity expansion and mismatch permutation), creates a single-entry
`KmerCounts` for each, and calls `create_seed_graph()` + `extend_graph()`
separately. If a primer generates 10 variant kmers, 10 nearly-identical
graphs are built, extended, pruned, and path-found independently. This is
the largest performance waste in the tool per the ROADMAP.

The proposed optimization: seed one graph with all forward primer kmers
simultaneously. But the scope of "all" matters.

### Option A: single graph per gene, all forward primer kmers

Seed one graph per gene with all of that gene's forward primer kmers as
start nodes, plus all reverse primer kmers as end nodes. This is what
`create_seed_graph()` already supports structurally — it accepts
`KmerCounts` collections and creates start/end nodes for each. The current
bottleneck is that `do_pcr()` wraps each forward kmer in a single-entry
`KmerCounts` and calls the function in a loop.

**Advantages:**
- Straightforward refactor — pass the full `forward_primer_kmers` instead
  of looping. `create_seed_graph()` already handles multiple kmers.
- Extension from multiple start nodes explores the same genomic region,
  building one graph instead of N nearly-identical ones.
- Pruning and path finding run once instead of N times.
- No gene boundary issues — each gene still has its own graph.

**Disadvantages:**
- If different forward primer kmer variants reach different parts of the
  genome (e.g., one matches a paralog), the single graph includes both.
  Currently, separate graphs would produce separate results that could be
  compared. With a merged graph, paralog paths coexist and must be
  distinguished during path finding.
- Path enumeration between start and end nodes becomes more complex when
  there are multiple start nodes. `all_simple_paths` would need to run
  from each start node, or from a virtual source node.

**Assessment:** The paralog concern is minor — paralogs would also appear
in the kmer count table and affect individual-kmer graphs anyway. The
current approach just happens to separate them by accident when different
primer kmers match different paralogs, which is unreliable. Proper paralog
handling is a v4.0 metagenomics concern. For v3.0, the merged graph is
correct behavior.

### Option B: single graph across all genes

Seed one graph with forward primer kmers from all genes simultaneously.

**Advantages:**
- Maximum sharing of graph construction work. Genes on the same chromosome
  region or sharing repetitive flanking sequence would share subgraph
  structure.
- Single extension pass over the kmer table.

**Disadvantages:**
- Genes targeting different parts of the genome produce disjoint subgraphs
  within one large graph — no construction savings, just bookkeeping
  overhead.
- If two genes' amplicons overlap or are close together, their graphs
  merge. Distinguishing which paths belong to which gene requires tracking
  which start nodes correspond to which gene — significant added
  complexity.
- Graph size could exceed `MAX_NUM_NODES` more easily when multiple genes
  contribute nodes.
- Pruning heuristics (tip clipping, bubble popping) would operate on a
  merged graph where coverage statistics mix signals from different loci.
- Loses the natural parallelism of per-gene processing (currently
  parallelized with rayon in Phase 4 of v2.0).

**Assessment:** The disadvantages outweigh the benefits. Most primer panels
target unrelated genomic loci, so the graphs would be disjoint anyway —
no sharing benefit. The boundary/attribution complexity and loss of rayon
parallelism make this a net negative.

### Decision

**Option A: single graph per gene, seeded with all forward primer kmers.**
This is the straightforward refactor with large performance benefit and
minimal risk. The main code change is in `do_pcr()`: instead of looping
over forward primer kmers and building separate graphs, pass the full
`forward_primer_kmers` collection to `create_seed_graph()`. Path finding
runs from each start node to each end node within the single graph.

Cross-gene graph merging (Option B) is not worth the complexity for v3.0
and conflicts with the existing per-gene rayon parallelism.

## Ambiguous read mapping during read threading

### Context

During Pass 2 (Phase 5), reads are threaded through the completed graph by
matching read kmers to graph edges. A read "maps" to a path through the
graph. But a read may map ambiguously — its kmers are consistent with
multiple paths (e.g., the read falls entirely within a region shared by
two bubble arms, or within a repeat).

The question: how should ambiguous reads contribute to edge-support counts
used in Phase 6 (threading-dependent traversal)?

### Option A: count for all valid paths

Every edge the read traverses gets +1 support, regardless of ambiguity.

**Advantages:** Simple implementation. No read is wasted.

**Disadvantages:** Inflates support counts for shared regions. Edges in
shared segments (before a bubble diverges or after it reconverges) get
disproportionately high counts compared to edges within the bubble. This
is misleading — it doesn't mean those edges are "more real," just that
more reads happen to be consistent with them. Doesn't help resolve
bubbles.

### Option B: best path only

Choose the single best-scoring path for each read (by alignment score,
kmer match count, or coverage consistency) and count only that path's
edges.

**Advantages:** Each read contributes exactly once per edge.

**Disadvantages:** For truly ambiguous reads (identical score on multiple
paths), the choice is arbitrary. Tie-breaking could systematically favor
one bubble arm over another based on arbitrary ordering. Also more complex
to implement — requires scoring all candidate paths before committing.

### Option C: fractional counting

If a read maps to N paths, each edge gets +1/N support.

**Advantages:** Statistically principled — analogous to multi-mapper
handling in RNA-seq (RSEM, Salmon). Shared regions get proportional
credit. Bubble-specific edges get fractional credit that reflects
genuine ambiguity.

**Disadvantages:** More complex bookkeeping (floating-point edge counts
instead of integers). The fractional signal may be too weak to be useful
for pruning decisions in Phase 6, especially at low coverage.

### Option D: only count unambiguous reads

Discard reads that map to multiple paths. Only count reads with a unique
path through the graph.

**Advantages:** Cleanest signal — every counted read unambiguously supports
its path. No inflation, no arbitrary choices.

**Disadvantages:** Wastes reads. At low coverage, discarding any reads
hurts. More critically, reads in shared regions (which are the majority —
most of the amplicon is shared, bubbles are typically short) would all be
discarded since they're consistent with all paths through the bubble.
This could leave most edges with zero read support.

### Analysis

The key insight is that ambiguity arises from where the read falls relative
to graph structure, not from read quality:

- **Reads spanning a bubble divergence or convergence point**: These are
  unambiguous — they support one specific path through the bubble. These
  are the most valuable reads for Phase 6.
- **Reads entirely within a shared region** (before divergence or after
  convergence): These are consistent with all paths. They confirm the
  shared structure exists but don't help resolve bubbles.
- **Reads entirely within one arm of a bubble**: These are unambiguous
  for that arm. Also valuable for Phase 6.

The reads that matter most — those spanning branch points — are
unambiguous. The ambiguous reads (in shared regions) don't carry
useful information for bubble resolution regardless of counting method.

### Decision

**Not yet decided — evaluate during Phase 5 implementation.** The analysis
suggests that the choice matters less than it might seem, because the
reads that carry useful signal (spanning branch points) are unambiguous.
For initial implementation, **Option A (count all)** is simplest and
provides a baseline. If Phase 6 shows that inflated counts in shared
regions cause problems for pruning heuristics, switch to **Option C
(fractional)** or **Option D (unambiguous only)**.

The important design constraint for Phase 5: store enough information per
edge to support any of these strategies. At minimum, track both total
read support (Option A) and unambiguous read support (Option D) per edge.
This avoids re-running Pass 2 if the counting strategy changes during
Phase 6 development.

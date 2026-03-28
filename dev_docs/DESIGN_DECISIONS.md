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

## Read retention and threading mechanics (Phases 4-5)

### Context

Phase 4 (read backend) streams reads in a second pass and must decide which
reads to retain and how to store them. Phase 5 (read threading) maps
retained reads to graph edges. These decisions are tightly coupled — the
storage format determines what's efficient in Phase 5.

### Phase 4: identifying relevant reads

During Pass 2, each read's kmers are extracted and checked against the
kmer set present in any gene's graph. In the de Bruijn graph, edges *are*
kmers (`DBEdge._kmer`), so "kmer is in the graph" means "kmer corresponds
to an edge."

**Pre-filter**: retain a read if it contains *any* kmer present in any
graph. This is fast (hash lookup per kmer) and over-inclusive — some
retained reads may not map to a contiguous path. That's fine; Phase 5
does the precise mapping. The alternative (requiring a contiguous run of
graph kmers) would need adjacency checks during the streaming pass, which
is more complex for marginal benefit.

At typical coverage for a 3 kb amplicon in a 3 Gb genome, ~210 reads per
amplicon are relevant (see coverage analysis above). Even with 10
amplicons and false positives, retained reads number in the low thousands
— negligible memory.

### Phase 4: storage format

**Don't reuse the existing `Read` struct.** It stores 2-bit encoded
sequence, which then needs to be decoded back to extract kmers in
`get_kmers()`. Since Phase 1 eliminates this round-trip for kmer counting,
the threading path should follow the same pattern.

**Store retained reads as:**

```rust
struct ThreadableRead {
    /// Ordered kmers extracted from this read, as u64 values.
    /// These are the canonical kmers, same as stored in the graph.
    kmers: Vec<u64>,
    /// Index of the mate read in the retained read vector, if paired.
    /// None for unpaired reads or if mate was not retained.
    mate_index: Option<usize>,
}
```

Why store kmers rather than raw sequence:
- Kmers are already extracted during the pre-filter check (don't discard
  work already done)
- Phase 5 needs kmers, not sequence — storing sequence would require
  re-extracting kmers
- u64 per kmer is compact (8 bytes per kmer vs ~1 byte per base for
  raw sequence, but a 150 bp read has ~130 kmers at k=21 = 1040 bytes
  vs 150 bytes). The trade-off favors raw sequence on space, but kmers
  on avoiding redundant computation.

**Alternative: store raw ASCII sequence and re-extract during threading.**
This uses less memory (~7x less per read) and is simpler. Given that
retained read counts are in the low thousands, the memory difference is
negligible either way (~1 MB vs ~7 MB for 10,000 reads). The re-extraction
cost is also trivial for thousands of reads.

**Decision: store raw ASCII sequence.** The memory savings don't matter at
this scale, but the simplicity does — fewer data structures, easier
debugging, and the kmer extraction code path (which Phase 1 optimizes)
gets reused rather than duplicated.

```rust
struct RetainedRead {
    /// Raw ASCII sequence (ACGT bytes)
    sequence: Vec<u8>,
    /// Index of the mate read in the retained read vector, if paired.
    /// None for unpaired reads or if mate was not retained.
    mate_index: Option<usize>,
}
```

### Phase 4: pairing

Pairing is tracked via `mate_index` in the retained read struct. When
`--paired` is set and reads are ingested alternately from R1/R2, each
pair of consecutive reads are mates. During retention:

- If both mates are retained: set `mate_index` on each to point to
  the other
- If only one mate is retained: set `mate_index = None` on the
  retained read (the mate didn't hit any graph kmer)
- Pairing info is consumed in Phase 6 (#101) for insert-size
  constraints

No separate graph or structure needed — a flat `Vec<RetainedRead>` with
cross-references via indices is sufficient for thousands of reads.

### Phase 5: mapping reads to graph edges

A read maps to a graph path by matching its kmers to graph edges in order.
The graph edges are kmers; nodes are (k-1)-mers. Consecutive kmers in a
read that both exist as graph edges and share a node (suffix of first =
prefix of second) form a **contiguous run** through the graph.

**Mapping algorithm:**

1. Extract kmers from retained read in order
2. Look up each kmer in the graph's edge set
3. Identify maximal contiguous runs: consecutive kmers where each pair
   shares a graph node (i.e., they are adjacent edges in the graph)
4. Each contiguous run defines a path through the graph

**Why require adjacency, not just presence:**

- A read with kmers [A, B, C, D] where A and B are adjacent graph edges
  and C and D are adjacent graph edges, but B and C are not, produces two
  runs: [A, B] and [C, D]. This correctly handles reads that span a
  sequencing error (the error breaks the kmer chain) or a region pruned
  from the graph.
- Simply checking presence without adjacency would incorrectly suggest
  that disconnected graph regions are linked by the read.

**Edge annotation:**

For each contiguous run, annotate every edge in the run with +1 read
support. A run of length 1 (single kmer in graph, neighbors not in graph)
still counts — it confirms that edge exists in a real read.

Track per edge:

- `read_support_total: u32` — incremented for every read whose run
  includes this edge (Option A from ambiguous mapping analysis)
- `read_support_unambiguous: u32` — incremented only when the read's
  full set of graph kmers maps to exactly one path (no branching)

The unambiguous count is computed after all runs for a read are found:
if all graph-matching kmers form a single contiguous run with no branch
points, it's unambiguous.

### Phasing annotations

Per-edge read support counts answer "is this edge real?" but don't
answer "which edges go together?" — i.e., they don't capture phasing.
A read (or contiguous run) that spans a branch point links the specific
incoming and outgoing edges, which is information not captured by per-edge
counts alone.

**Branch-point phasing:** When a contiguous run passes through a node
with in-degree > 1 or out-degree > 1 (a branch point), record which
specific incoming edge and outgoing edge the read connects. Store as
a set of observed (edge_in, edge_out) pairs per branch-point node, with
counts.

```rust
/// Observed read-supported links through a branch point.
struct BranchPhasing {
    node: NodeIndex,
    /// (incoming_edge, outgoing_edge) -> count of reads supporting this link
    links: HashMap<(EdgeIndex, EdgeIndex), u32>,
}
```

This is distinct from paired-end constraints (#101), which operate at
insert-size scale (hundreds of bases apart). Branch-point phasing
operates at read scale (within a single read's contiguous run through
the graph). Both are useful:

- **Branch-point phasing** (read-scale): resolves which arm of a bubble
  connects to which arm of an adjacent bubble — local haplotype structure
- **Paired-end phasing** (insert-scale): links branches that are further
  apart than a single read can span — longer-range haplotype structure

## Graph annotation model: preserve structure, defer decisions

### Context

The current v2.0 approach destructively prunes the graph before path
finding: `pop_balloons()`, `remove_side_branches()`, `remove_orphan_nodes()`
all delete graph structure. Information lost during pruning cannot be
recovered.

### Decision

**Annotation-only model.** Read threading (Phase 5) annotates the graph
but never modifies its structure. All structural decisions — which paths
to report, which bubbles represent real variants vs artifacts — are
deferred to path finding and sequence emission. Bubbles and other
variants are retained all the way until sequences are emitted.

The pipeline becomes:

1. **Build graph** from kmers (Phase 3 construction)
2. **Light structural cleanup** — remove clearly spurious structure that
   would make downstream processing intractable:
   - Dead-end tips shorter than k with very low coverage (these are
     almost certainly sequencing errors, not real variants)
   - Disconnected components not reachable from any start node
   - This is the only step that edits graph structure
3. **Annotate with kmer coverage** (already available from construction)
4. **Annotate with read support** (Phase 5 — optional, skipped with
   `--no-read-threading`)
5. **Annotate with phasing** (Phase 5 — branch-point links from read
   runs; Phase 6 — paired-end links)
6. **Score and select paths** using all available annotations (Phase 6).
   Bubbles are resolved here by choosing the best-supported path(s),
   not by deleting the alternatives from the graph. Multiple paths can
   be emitted if they represent real variants (preparation for v4.0
   metagenomics).
7. **Emit sequences** from selected paths

This means Phase 3 "pruning" (tip clipping, bubble popping) becomes
Phase 3 "light cleanup + annotation" — the aggressive pruning is
replaced by annotation-informed path selection. The pluggable scoring
interface designed in #90 takes coverage, read support, and phasing as
inputs and returns path scores.

### Advantages

- No information loss — all evidence is available at decision time
- Same code path works with or without read threading (just fewer
  annotations available)
- Natural preparation for v4.0 metagenomics (multiple real variants
  coexist in the graph, scored by coverage profile)
- Easier to debug — the graph before and after annotation is the same
  structure, so you can inspect what annotations led to which decisions

### Concern: path-finding efficiency

Preserving all bubbles increases the number of paths. Mitigations:

- Light structural cleanup removes the worst offenders (low-coverage
  tips and disconnected components)
- Coverage-weighted path finding (Phase 3 #93) follows high-confidence
  edges first, naturally avoiding low-coverage artifacts without
  deleting them
- `MAX_NUM_PATHS_PER_PAIR` and `MAX_NUM_AMPLICONS` limits still bound
  output
- If path enumeration is still too expensive, annotations can be used
  to mask low-confidence edges during path finding (treat them as absent
  without deleting them — reversible)

### Phased implementation

- **Phase 3**: replace destructive pruning with light structural cleanup.
  Implement coverage-weighted path finding that uses kmer coverage
  annotations to score paths. Bubbles and variants are preserved.
- **Phase 4**: implement `RetainedRead`, pre-filtering, pairing. No
  graph annotation yet. Verify retained read counts match expectations.
- **Phase 5**: implement the mapping algorithm, edge annotation
  (total/unambiguous counts + branch-point phasing). Benchmark.
- **Phase 6**: add read-support and phasing signals to path scoring.
  Paired-end constraints for longer-range phasing. Second scoring pass
  with the richer annotation set.

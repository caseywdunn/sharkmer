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

## Graph traversal: v1 shortcomings and Phase 3 approach

### V1 shortcomings

The current graph construction and traversal has several interrelated
problems that limit gene recovery — particularly for single-copy nuclear
genes at moderate coverage.

**1. One graph per forward primer kmer.** `do_pcr()` iterates over every
forward primer kmer (produced by ambiguity expansion and mismatch
permutation) and builds a separate graph for each. If a primer produces 10
kmer variants, 10 nearly-identical graphs are built, extended, pruned, and
path-found independently. This is the single largest performance waste in
the tool. (Already discussed in the "Single graph seeded with all forward
primer kmers" section above — the fix is decided, this section focuses on
the traversal problems within each graph.)

**2. Greedy forward-only extension.** `extend_graph()` starts at seed
nodes and extends outward one node at a time, checking all four possible
successor kmers. This is conceptually simple but has a structural problem:
the extension frontier expands breadth-first across *all* growing tips
simultaneously. In repetitive regions of the genome, a single high-copy
kmer can inject hundreds of successors into the frontier, and each of those
can branch again. The graph balloons not because the amplicon region is
complex, but because the extension wandered into a repeat.

**3. Ad-hoc ballooning controls.** The current code has multiple overlapping
heuristics to contain graph explosion, none of which are principled:

- **Backward degree check** (graph.rs lines 471-500): if a node and its
  recent ancestors all have degree > 2, mark it terminal. This is a proxy
  for "we're in a repeat" but fires based on local topology, not coverage.
  It can falsely terminate real amplicon branches that happen to have a
  few consecutive branch points (e.g., heterozygous SNPs in close
  proximity).
- **Long-range degree check** (graph.rs lines 489-499): over a 15-node
  window, if 3+ ancestors have degree > 1, terminate. Same proxy problem
  with a wider window.
- **BALLOONING_COUNT_THRESHOLD_MULTIPLIER** (graph.rs lines 550-562):
  reject edges with count > 10× the median edge count. This prevents
  high-copy repeat kmers from being added, but the 10× multiplier is
  arbitrary and the median shifts as the graph grows, making the behavior
  unpredictable.
- **pop_balloons()** (pruning.rs): periodically count descendants within
  depth 4. If the count exceeds 4^3 = 64, delete the entire subtree.
  This is the most destructive heuristic — it removes real structure along
  with noise, and the threshold is topology-based (branching factor) rather
  than coverage-based.

These heuristics interact in hard-to-predict ways. A graph that happens to
balloon on one extension step may get pop_balloons'd, removing nodes that
would have been useful if the extension had proceeded in a different order.
The result is non-deterministic in practice (depending on HashMap iteration
order and when periodic checks fire) despite the sorted-kmer determinism
efforts.

**4. Destructive pruning before path finding.** After extension,
`remove_side_branches()` deletes all dead-end nodes that aren't end nodes,
regardless of their coverage. A dead end with 50× coverage (possibly a
real variant whose extension was terminated by one of the ballooning
heuristics) is treated the same as a dead end with 1× coverage (almost
certainly a sequencing error). `remove_orphan_nodes()` then deletes end
nodes with no incoming edges. Information lost during pruning cannot be
recovered.

**5. `all_simple_paths` enumeration is exponential.** Path finding uses
petgraph's `all_simple_paths`, which enumerates every distinct path between
start and end nodes. For a graph with B bubbles, this produces up to 2^B
paths. The `MAX_NUM_PATHS_PER_PAIR = 20` cap prevents runaway computation
but means that in a complex graph, only an arbitrary subset of paths is
considered. There is no guarantee that the best path (by any quality
metric) is among the first 20 enumerated.

**6. Cycle avoidance is absolute.** `would_form_cycle()` does a full BFS
from the child node to check if the parent is reachable. If so, the edge
is rejected entirely. This prevents any repeat from being traversed more
than once, which is correct for simple tandem repeats but wrong for
legitimate biological structures like tandem duplications or gene
conversion tracts where the same kmer sequence genuinely appears twice
in the amplicon.

**7. Path scoring is minimal.** Paths are ranked solely by `kmer_min_count`
— the lowest edge count along the path. This is a reasonable first-order
heuristic but ignores the distribution of counts along the path. A path
with one low edge and 200 high edges scores worse than a path with
uniformly mediocre edges.

**8. Threshold annealing rebuilds from scratch.** The stepping-down
threshold approach (high threshold → low threshold via
`compute_coverage_thresholds()`) is sound in principle: start with
high-confidence edges, then progressively include lower-coverage edges.
But each step rebuilds the graph from the seed, because `extend_graph()`
takes the immutable `seed_graph` and produces a new graph. Work done at
higher thresholds is discarded. And if a product is found at any threshold
step, the remaining steps are skipped (`break` at line 474), so lower
thresholds never run — even if they would have found a better product.

### How other de Bruijn assemblers handle traversal

The problems above are well-studied in the genome assembly literature.
Here is how the major assemblers address them:

**Tip clipping (Velvet, ABySS, SPAdes, MEGAHIT).** Dead-end paths shorter
than 2k with low coverage are removed as sequencing error artifacts. This
is the one form of structural editing that all assemblers agree on. The key
distinction from sharkmer's `remove_side_branches()` is the *coverage
criterion*: assemblers only clip tips that are both short AND low-coverage.
High-coverage dead ends are preserved because they may represent real
structural variants whose connection was broken by an under-sampled kmer.

**Bubble detection and resolution (Velvet, SPAdes).** Rather than deleting
bubbles during construction, assemblers identify them explicitly: a bubble
is a pair of paths that diverge from a node and reconverge at another
node. Velvet's Tour Bus algorithm does a modified Dijkstra traversal to
find bubbles and then merges them (keeping the higher-coverage arm).
SPAdes detects bubbles and marks them as alternative paths without deleting
either arm, deferring the choice to later stages with more information.

**Coverage-guided traversal (SPAdes, MEGAHIT).** Instead of enumerating all
paths and then scoring them, assemblers traverse the graph greedily,
following the highest-coverage outgoing edge at each branch point. When a
branch point is encountered, the traversal can:
- Follow the dominant edge immediately (greedy best-first)
- Record the branch for later exploration (depth-first with backtracking)
- Fork into parallel paths at the branch (breadth-first)

SPAdes uses coverage to distinguish between errors (low-coverage edges) and
real polymorphism (coverage proportional to the organism's allele
frequency). MEGAHIT's "mercy kmer" approach recovers low-coverage paths
that were filtered out of the initial graph.

**Repeat resolution.** Assemblers handle repeats through multiple
complementary strategies:
- **Read threading / read coherence** (SPAdes, Velvet): map original reads
  back to the graph and use the sequence of kmers within a single read to
  link edges across branch points. This is the most powerful technique and
  is what our Phase 5-6 implements.
- **Paired-end constraints** (Velvet scaffolding, SPAdes): use insert size
  to link branches that are further apart than a single read spans.
- **Coverage depth**: in a unique region, expected coverage is C; in a
  2-copy repeat, expected coverage is 2C. Edges with anomalously high
  coverage are flagged as repetitive.

**Iterative k-mer assembly (IDBA, SPAdes).** Rather than building one
graph at a fixed k, build graphs at multiple k values. Small k recovers
connections in low-coverage regions; large k resolves repeats. SPAdes
combines information across k values in its assembly graph. This is
analogous to sharkmer's threshold annealing but operates on the k
parameter rather than the count threshold. The key insight is that the
graphs at different parameters should *share information*, not be built
independently.

### Phase 3 approach

The new approach addresses each v1 shortcoming. Here are the key design
decisions and how specific challenges are handled:

**Single graph per gene, seeded with all primer kmers.** (Already decided
above.) Eliminates the N-graph-per-gene waste. All forward primer kmers
seed start nodes in one graph; all reverse primer kmers seed end nodes.

**Reverse extension with coverage awareness.** Instead of extending
only forward from start nodes, extend in reverse from end nodes as well.
The graph is complete when the two frontiers meet (start-rooted subgraph
connects to end-rooted subgraph) or when extension exhausts. This
halves the search space on average and ensures that the graph includes
the amplicon interior even if one primer site is in a complex region.

**Coverage-based edge filtering replaces topology-based ballooning
heuristics.** The backward-degree checks, descendant counting, and
pop_balloons are all removed. Instead:

- During extension, an edge is added only if its count ≥ `min_count`
  (this part is unchanged).
- After extension, edges are annotated with a "coverage ratio": the
  edge's count divided by the local median (median of counts in a
  neighborhood of depth ~5). Edges with a ratio far above 1.0 are flagged
  as potentially repetitive. Edges with a ratio far below 1.0 are flagged
  as potential errors.
- These annotations inform path scoring but do not cause edge deletion.

This addresses the "common kmer" problem directly: if the extension
encounters a kmer with count 10,000× the local median, that edge gets a
very high coverage ratio annotation. The path scorer can then strongly
penalize paths that traverse it, effectively routing around repeats without
deleting them from the graph. If a later phase (read threading) provides
evidence that the high-copy edge is genuinely part of the amplicon, the
evidence overrides the coverage penalty.

**Threshold annealing: incremental, not rebuild.** The stepping-down
threshold approach is retained because it is fundamentally sound — start
conservative, relax progressively. But the implementation changes:

- The graph is built once at the highest threshold.
- At each subsequent (lower) threshold, only newly qualifying edges are
  added to the existing graph. Edges and nodes from higher thresholds
  are preserved.
- Path finding runs after each threshold step. If a complete
  start-to-end path exists, it is recorded but extension continues at
  lower thresholds. All paths found across all thresholds are collected
  and scored together at the end.
- This means a high-quality path found at a high threshold is not
  discarded just because a lower threshold also produces paths. The
  scorer picks the best among all candidates.

The practical benefit: at high thresholds, only high-confidence edges
exist and the graph is small/fast. Most amplicons are found here. Lower
thresholds add edges incrementally, and the graph grows modestly because
most of the structure is already in place.

**Controlling graph size at low thresholds.** The concern with walking
thresholds down to `min_count` (which may be 2) is that at low thresholds,
many noise edges qualify and the graph can explode. The controls are:

1. **Reachability constraint.** After each threshold step, remove edges
   and nodes that are not on any path from a start node to an end node
   (or could not plausibly become part of such a path). Specifically,
   keep only the connected component that contains at least one start
   node and at least one end node. This is the "light structural cleanup"
   from the annotation model decision above.

2. **Extension frontier limiting.** Rather than extending all unvisited
   nodes simultaneously, extend from the current frontier (nodes added
   in the most recent threshold step). At each step, only nodes within
   `max_length` graph distance from a start node (or end node, for
   reverse extension) are eligible for extension. This prevents the
   graph from wandering arbitrarily far from the amplicon region.

3. **MAX_NUM_NODES limit is retained** as a hard safety bound.

**Coverage-weighted best-path traversal replaces `all_simple_paths`.**
Instead of enumerating all paths and taking the first 20, use a
coverage-weighted traversal:

1. From each start node, do a priority-queue traversal toward end nodes.
   The priority of an edge is its coverage (or a score combining coverage
   and coverage ratio). At each branch point, the traversal follows the
   highest-scoring edge first.
2. When an end node is reached, record the path. Continue exploring
   to find alternative paths (backtrack at branch points), but with a
   budget: stop after `MAX_NUM_PATHS_PER_PAIR` paths per start-end pair.
3. Because the traversal is coverage-ordered, the first paths found are
   the highest-quality paths. The arbitrary first-20 problem is eliminated.

This is essentially Dijkstra's algorithm with coverage as the edge weight,
run on a DAG (after cycle handling — see below). It naturally routes around
low-coverage noise and through high-coverage amplicon sequence.

**Cycle handling.** The absolute cycle ban is replaced with bounded repeat
traversal. During path finding, a node can be visited up to N times
(configurable, default 2). This allows the traversal to pass through a
tandem duplication twice but prevents infinite loops. For graph
construction, cycles are allowed in the graph structure — the acyclicity
constraint moves from construction time to path-finding time, where it
can be enforced per-path rather than globally.

**Pluggable path scoring.** Paths are scored by a function that takes
multiple signals as input:

- `kmer_min_count`: minimum edge count along the path (current metric,
  retained as one input)
- `kmer_mean_count` and `kmer_median_count`: distribution statistics
- `coverage_consistency`: variance of edge counts along the path (a good
  path through a unique region should have roughly uniform coverage)
- `coverage_ratio_penalty`: sum of log-ratios for edges flagged as
  potentially repetitive
- `path_length`: penalize paths far from the expected amplicon length,
  if known

Phase 3 uses only kmer-coverage-based signals. Phases 5-6 add read
support and phasing signals to the same interface without restructuring.

**Light structural cleanup.** The only structural edits to the graph are:

- Remove dead-end tips shorter than k with coverage below a fraction
  of the local median (e.g., < 0.1× local median). These are almost
  certainly sequencing errors.
- Remove connected components not reachable from any start or end node.
- Remove edges and nodes that cannot be part of any start-to-end path
  within the length bounds.

Everything else — bubbles, side branches with real coverage, repeat
edges — stays in the graph and is handled by the scorer.

### Threshold selection: how to avoid being too high or too low

This is the most delicate practical challenge. A threshold too high misses
real amplicon edges (especially in low-coverage regions of the amplicon,
which are common near GC-biased or repetitive primer binding sites). A
threshold too low admits noise and repeats, bloating the graph.

The annealing approach addresses this by trying thresholds from high to
low, but the v1 implementation gives up as soon as any product is found.
The Phase 3 approach collects products across all thresholds and scores
them together. This means:

- **Too-high thresholds** produce incomplete graphs (no start-to-end
  path) or paths with gaps. These are simply not scored — no harm done,
  and the lower thresholds fill in.
- **Too-low thresholds** produce noisy graphs with many paths. The scorer
  penalizes paths with high coverage variance and low minimum count. If a
  clean path was already found at a higher threshold, the noisy paths from
  lower thresholds will score worse.
- **The "just right" threshold** produces a graph with a clear, high-
  coverage path from start to end. This path scores best and is emitted.

The threshold sequence itself is unchanged: start at
`primer_count / COVERAGE_MULTIPLIER` and step down to `min_count` in
`COVERAGE_STEPS` steps. These constants may need tuning based on Phase 3
benchmarks, but the framework is robust to the exact values because all
threshold levels contribute candidates rather than short-circuiting.

**Why not use a single adaptive threshold?** Some assemblers (e.g.,
MEGAHIT) use iterative approaches that automatically find the right
threshold by monitoring graph quality metrics. This is appealing but adds
complexity. The multi-threshold approach is simpler to implement, reason
about, and debug: each threshold produces a concrete graph that can be
visualized with `--dump-graph`. If benchmarks show that the multi-threshold
approach misses cases that an adaptive threshold would catch, this can be
revisited.

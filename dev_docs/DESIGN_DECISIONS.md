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

The question is whether to build one graph per primer kmer variant, one
graph per gene, or one graph across all genes.

In v1, `do_pcr()` iterated over forward primer kmers (generated by
ambiguity expansion and mismatch permutation), created a single-entry
`KmerCounts` for each, and called `create_seed_graph()` + `extend_graph()`
separately. If a primer generated 10 variant kmers, 10 nearly-identical
graphs were built, extended, pruned, and path-found independently. This
was the largest performance waste in the tool.

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
minimal risk.

Cross-gene graph merging (Option B) is not worth the complexity and
conflicts with per-gene parallelism.

**Status (implemented in Phase 3):** `do_pcr()` now calls
`create_seed_graph()` once per gene with the full `forward_primer_kmers`
and `reverse_primer_kmers` collections. All forward primer kmers seed
start nodes and all reverse primer kmers seed end nodes in a single
graph. The graph is extended through multiple coverage threshold steps,
accumulating structure rather than rebuilding. Path finding runs from
each start node to each end node within the single graph.

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

**Option A (count all), with both total and unambiguous counts tracked.**

The analysis suggests that the choice matters less than expected, because
the reads that carry useful signal (spanning branch points) are
unambiguous.

**Status (implemented in Phase 5):** `thread_reads()` in `threading.rs`
implements Option A: every edge in a read's contiguous run gets +1
support. Both `total_read_support` and `unambiguous_read_support` are
tracked per edge via `EdgeReadSupport`, plus branch-point phasing via
`BranchLink` and paired-end links via `PairedEndLink`. This preserves
the option to switch counting strategies without re-running Pass 2.

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

### Implementation status

All phases are implemented:

- **Phase 3** (complete): `pruning.rs` has `remove_low_coverage_tips()`
  and `reachability_pruning()` — light structural cleanup only. The old
  destructive `pop_balloons()`, `remove_side_branches()`, and
  `remove_orphan_nodes()` are removed. Path finding in `paths.rs` uses
  coverage-weighted custom DFS (not petgraph `all_simple_paths`).
- **Phase 4** (complete): `RetainedRead` in `io.rs`, pre-filtering via
  `OligoFilter` (bloom filter + AHashSet), paired-end support.
- **Phase 5** (complete): `thread_reads()` and `thread_reads_paired()`
  in `threading.rs`. Edge annotation with `EdgeReadSupport`
  (total/unambiguous counts) and `BranchLink` (branch-point phasing).
- **Phase 6** (complete): `resolve_bubbles()` in `bubble.rs` detects
  simple bubbles, ranks branches by read support + phasing, returns
  edge preferences that feed into path scoring.

## Graph traversal: v1 shortcomings and v3 approach

### V1 shortcomings (all addressed in v3)

The v1 graph construction and traversal had several interrelated
problems that limited gene recovery — particularly for single-copy nuclear
genes at moderate coverage. All are addressed in v3 (Phases 3-7).

**1. One graph per forward primer kmer.** Fixed: single graph per gene
seeded with all primer kmers. See "Single graph seeded with all forward
primer kmers" section above.

**2. Greedy forward-only extension.** Fixed: v3 has both `extend_graph()`
(forward from start nodes) and `extend_graph_reverse()` (backward from
end nodes). The two frontiers meet in the amplicon interior.

**3. Ad-hoc ballooning controls.** Fixed: v1 had backward degree checks,
long-range degree checks, `BALLOONING_COUNT_THRESHOLD_MULTIPLIER`, and
`pop_balloons()` — all removed. V3 uses coverage-based edge filtering
(`high_coverage_ratio` to reject edges far above the median) and
`MAX_NUM_NODES` as a hard safety bound. Seed evaluation (`seed_eval.rs`)
filters off-target seeds before full extension based on bounded local
exploration and branching ratio.

**4. Destructive pruning before path finding.** Fixed: v1's
`remove_side_branches()` and `remove_orphan_nodes()` are replaced by
`remove_low_coverage_tips()` (coverage-aware, only clips short tips below
a fraction of local median) and `reachability_pruning()` (removes nodes
not on any start-to-end path). See "Graph annotation model" above.

**5. `all_simple_paths` enumeration is exponential.** Fixed: v3 uses a
custom coverage-weighted DFS in `paths.rs` (`get_assembly_paths()`).
Edges are sorted by coverage (and bubble-resolution preferences when
available), so the highest-quality paths are found first. Bounded by
`max_paths_per_pair` and `max_dfs_states`.

**6. Cycle avoidance is absolute.** Fixed: v3 allows cycles in the graph
structure. Cycle control is enforced during path finding via bounded
node revisitation (`max_node_visits`, default 2), allowing traversal
through tandem duplications while preventing infinite loops.

**7. Path scoring is minimal.** Fixed: v3 scores paths using multiple
signals including kmer min/mean/median counts, coverage consistency,
and (when available) read support and bubble-resolution preferences
from `resolve_bubbles()`.

**8. Threshold annealing rebuilds from scratch.** Fixed: v3 accumulates
the graph across threshold steps. `prepare_for_lower_threshold()` resets
terminal flags so previously-terminal nodes can extend at the lower
threshold. The graph grows incrementally rather than being rebuilt.

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

### V3 approach (implemented)

The v3 approach addresses each v1 shortcoming. Here are the key design
decisions and how specific challenges are handled:

**Single graph per gene, seeded with all primer kmers.** All forward
primer kmers seed start nodes in one graph; all reverse primer kmers
seed end nodes. Eliminates the v1 N-graph-per-gene waste.

**Bidirectional extension.** `extend_graph()` extends forward from start
nodes; `extend_graph_reverse()` extends backward from end nodes. Forward
extension runs first; reverse extension is skipped if forward already
reached an end node (to avoid wasting the node budget on off-target
reverse seeds).

**Coverage-based edge filtering.** During extension, edges with count
far above the median edge count (`high_coverage_ratio`, default 10x)
are skipped — they likely lead into repetitive regions. This replaces
the v1 topology-based ballooning heuristics.

**Seed evaluation.** Before full graph extension, `evaluate_seeds()`
performs bounded local exploration of each seed node. Seeds that show
excessive branching, terminate too quickly, or exhaust their node budget
with high branching ratios are marked terminal and excluded from full
extension. Seeds that reach an opposite-direction seed are kept for
early product recovery.

**Threshold annealing: incremental, not rebuild.** The stepping-down
threshold approach starts conservative and relaxes progressively:

- The graph is built once at the highest threshold.
- At each subsequent (lower) threshold, `prepare_for_lower_threshold()`
  resets terminal flags so previously-blocked nodes can extend further.
  Only newly qualifying edges are added; existing structure is preserved.
- Path finding runs after each threshold step. Products found at any
  threshold are collected and scored together.
- At high thresholds, only high-confidence edges exist and the graph is
  small/fast. Most amplicons are found here. Lower thresholds add edges
  incrementally.

**Graph size controls:**

1. **Seed evaluation** filters off-target seeds before full extension.
2. **High-coverage edge filtering** skips edges far above the median.
3. **Reachability pruning** removes nodes not on any start-to-end path.
4. **Low-coverage tip removal** clips short dead-end tips.
5. **`MAX_NUM_NODES`** (default 50,000) is a hard safety bound.

**Coverage-weighted DFS path finding.** `get_assembly_paths()` uses a
custom stack-based DFS instead of petgraph's `all_simple_paths`:

1. From each start node, traverse toward end nodes. Edges are sorted by
   coverage (and bubble-resolution preferences when available), so the
   highest-scoring edges are explored first.
2. When an end node is reached, record the path. Continue exploring
   alternative paths with backtracking, bounded by `max_paths_per_pair`
   and `max_dfs_states`.
3. Because the traversal is coverage-ordered, the first paths found are
   the highest-quality paths.

**Cycle handling.** Cycles are allowed in the graph structure. During path
finding, a node can be visited up to `max_node_visits` times (default 2).
This allows traversal through tandem duplications while preventing
infinite loops.

**Path scoring.** Paths are scored using `PathScore` which combines:

- `kmer_min_count`, `kmer_mean_count`, `kmer_median_count`
- `read_support_total`, `read_support_min` (when read threading is active)
- Bubble-resolution edge preferences from `resolve_bubbles()`

**Light structural cleanup.** The only structural edits to the graph:

- `remove_low_coverage_tips()`: tips shorter than k with coverage below
  `tip_coverage_fraction` of the local median.
- `reachability_pruning()`: removes nodes not on any start-to-end path.

Everything else — bubbles, side branches with real coverage, repeat
edges — stays in the graph and is handled by the scorer.

### Threshold selection: how to avoid being too high or too low

This is the most delicate practical challenge. A threshold too high misses
real amplicon edges (especially in low-coverage regions of the amplicon,
which are common near GC-biased or repetitive primer binding sites). A
threshold too low admits noise and repeats, bloating the graph.

The annealing approach addresses this by trying thresholds from high to
low. Products found at each threshold are collected and scored together:

- **Too-high thresholds** produce incomplete graphs (no start-to-end
  path) or paths with gaps. These are simply not scored — no harm done,
  and the lower thresholds fill in.
- **Too-low thresholds** produce noisy graphs with many paths. The scorer
  penalizes paths with high coverage variance and low minimum count. If a
  clean path was already found at a higher threshold, the noisy paths from
  lower thresholds will score worse.
- **The "just right" threshold** produces a graph with a clear, high-
  coverage path from start to end. This path scores best and is emitted.

The threshold sequence starts at
`primer_count / COVERAGE_MULTIPLIER` and steps down to `min_count` in
`COVERAGE_STEPS` steps. These constants may need tuning based on
benchmarks, but the framework is robust to the exact values because all
threshold levels contribute candidates rather than short-circuiting.

**Which primer count statistic to use matters.** Extension thresholds
and seed eval thresholds serve different purposes and need different
primer count statistics:

- **Extension thresholds use the max** primer kmer count (`max` of
  forward and reverse minimums). Starting high keeps the graph focused on
  confident kmers first, avoiding premature node budget exhaustion. If
  the initial threshold is too low, the graph extends aggressively into
  low-coverage regions and hits the 50K node budget before connecting
  forward and reverse primer sites.

- **Seed eval threshold uses the median** primer kmer count. Degenerate
  primers (ambiguity codes H, D, Y, N, R) generate off-target genomic
  matches that inflate the max far above real amplicon coverage. Using
  max for seed eval makes the threshold too stringent: real seeds cannot
  extend even a single node, and every seed is falsely abandoned. The
  median is robust to these outliers. See `failure_analysis.md` "Root
  cause: seed eval threshold" for the diagnostic evidence.

This split was validated by benchmark regression: using median for both
caused 41 gene regressions (graphs exhausting node budget), while using
max for both caused seed eval failures for degenerate primers. The split
approach recovers all genes from both failure modes.

**Why not use a single adaptive threshold?** Some assemblers (e.g.,
MEGAHIT) use iterative approaches that automatically find the right
threshold by monitoring graph quality metrics. This is appealing but adds
complexity. The multi-threshold approach is simpler to implement, reason
about, and debug: each threshold produces a concrete graph that can be
visualized with `--dump-graph`. If benchmarks show that the multi-threshold
approach misses cases that an adaptive threshold would catch, this can be
revisited.

### Node budget sizing and extension cost model

**Extension cost was O(n²), now O(n).** The original `extend_graph()`
loop iterated all node indices each pass to find unvisited nodes,
making the cost of adding n nodes proportional to n². Doubling the node
budget quadrupled runtime. Benchmark evidence (commit 5ed7445,
2026-04-03): raising the global budget from 50K to 200K (4×) increased
Agalma elegans runtime from 22s to 347s (~16×, consistent with 4²).

**Fix (commit bfea426):** Both `extend_graph()` and
`extend_graph_reverse()` now use a `VecDeque` frontier queue. Unvisited
nodes are pushed to the queue when created; the main loop pops from
the front. This makes extension O(n) in the number of nodes added.
Benchmark: Agalma with `--max-nodes 200000` dropped from 347s to 65s
(5.3× speedup).

**Current budget defaults.** The global `--node-budget-global` is 50K.
Per-component budgets (default 20K via `--node-budget-component`) bound
each connected component independently, preventing any single off-target
component from consuming the entire global budget.

**Why 50K global.** A sweep across 14 benchmark samples at 1M reads
(commit f7cbc90, 2026-04-03) shows diminishing returns above 50K:

| Budget | Total genes | Total time | Marginal genes | Marginal time |
| ---: | ---: | ---: | ---: | ---: |
| 50K | 106 | 159s | — | — |
| 100K | 109 | 214s | +3 | +55s |
| 200K | 110 | 296s | +1 | +82s |

The 50K→100K step recovers 3 additional genes (Gryllus) but costs 34%
more time. The 100K→200K step recovers only 1 more gene (Drosophila)
at 38% more time. Most of the extra time is spent on off-target
components that ultimately fail — Liriodendron goes from 9s to 67s
with no gene gain at 200K.

Users who need the marginal genes can increase the budget with
`--node-budget-global 100000`. The per-component budget (20K) provides
the real protection against runaway extension; the global budget is
a backstop.

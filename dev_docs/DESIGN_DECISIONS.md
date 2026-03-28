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

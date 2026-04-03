# sPCR failure analysis

A systematic analysis of why sharkmer's in silico PCR fails to recover
amplicons, organized by root cause. Used to guide development priorities,
parameter tuning, and user documentation.

## Failure categories

Failures fall into three groups: limitations of the target sequence,
limitations of the input data, and limitations of the method.

### Sequence limitations

These are properties of the target organism or locus that make recovery
inherently difficult.

#### Low-complexity / AT-rich regions

The amplicon contains a stretch where kmer diversity is low (typically
high AT content). In these regions, fewer distinct kmers exist, so even
with adequate read depth, individual kmer counts may fall below the
minimum threshold. The de Bruijn graph fragments into disconnected
forward and backward components with a gap in the low-complexity
region.

**Mechanism:** Shorter kmers have higher coverage per unique kmer
because fewer possible kmers exist in a low-complexity region. At k=31,
an 85% AT region has most 31-mers appearing only once; at k=21, the
same region may have kmers appearing 2+ times. However, below ~k=15,
kmers become too non-specific and the graph explodes with branching.

**Mitigation (user):** Decrease `-k` (e.g., from 31 to 21). This
bridges moderate AT-rich gaps but increases graph complexity elsewhere.

**Mitigation (developer):** Variable-k assembly — detect AT-rich gaps
and retry the gap region with smaller k. This is a significant
architectural change (see SPAdes/IDBA iterative k-mer assembly).

#### Repetitive amplicon flanking regions

The primer binding site or the region immediately adjacent is in or
near a repetitive element (e.g., retrotransposon insertions in rDNA,
tandem repeats). The de Bruijn graph includes edges from many genomic
copies, creating a dense tangle that exceeds the node budget or makes
path finding intractable.

**Mitigation (user):** Increase `--max-nodes` (hidden). May not help
if the repeat complexity is inherently too high.

**Mitigation (developer):** Read threading with repeat-aware scoring
(`--read-threading`) can help resolve which path through the repeat
is the real amplicon.

#### Amplicon too large for short-read graph assembly

The expected amplicon exceeds what de Bruijn graph assembly can
construct from short reads (~150 bp). For amplicons >2-3 kb, the
probability of spanning any given position with reads drops, and the
graph is more likely to fragment. Amplicons >5 kb are essentially
unrecoverable from 150 bp reads alone.

**Example:** Drosophila ITS amplicon is ~4835 bp.

**Mitigation (user):** Use longer reads (if available) or redesign
primers for a shorter amplicon.

**Mitigation (developer):** No practical fix with short reads. Document
the amplicon length limit.

#### rDNA copy-number variation

Ribosomal RNA genes exist in hundreds of tandem copies that are not
perfectly identical. Sequence variants across copies create bubbles and
branches in the de Bruijn graph that can exceed path-finding budgets.

**Mitigation (developer):** Bubble resolution with read threading
(`resolve_bubbles()`) already addresses simple cases. More aggressive
consensus calling across bubble arms could help.

### Data limitations

These are properties of the sequencing data that can be addressed by
providing different input.

#### Insufficient read coverage at the target locus

The most common failure cause. At low read counts, some kmers in the
amplicon region are never observed (count 0) or observed only once
(below min_count=2). The graph has gaps where coverage drops out.

The required read count depends on genome size. For a 150 bp read
library:

| Genome | Size | Reads for ~5x mt coverage | Reads for ~5x nuclear |
| --- | --- | --- | --- |
| Drosophila | 144 Mb | ~1M | ~5M |
| Coral (Porites) | 542 Mb | ~4M | ~18M |
| Human | 3.1 Gb | ~20M | ~103M |

Mitochondrial and plastid genes are easier because organellar copy
number provides ~100x enrichment over single-copy nuclear genes.

**Mitigation (user):** Increase `--max-reads`. More reads = higher
coverage per kmer = fewer gaps.

#### Uneven coverage across the amplicon

Even with adequate average coverage, specific regions may have low
coverage due to GC bias, library preparation artifacts, or secondary
structure. This creates the same gap problem as insufficient total
coverage but is harder to predict.

**Mitigation (user):** Increase `--max-reads` to raise the coverage
floor. Decrease `-k` to reduce the number of distinct kmers needed.

### Method limitations

These are failure modes caused by the algorithm's design or parameter
settings, addressable through code changes or parameter tuning.

#### Off-target seed nodes consuming the node budget

Degenerate primers (with ambiguity codes) match many genomic loci.
Each match creates a seed node that extends into an off-target
subgraph. When many off-target seeds survive seed evaluation, their
combined extension exhausts the 50K node budget before the on-target
amplicon is connected.

This failure is sample-specific even though primer degeneracy is fixed:
different genomes have different off-target match profiles.

**Mitigation (user):** Increase `--max-nodes` (hidden) to allow a
larger graph. Decrease primer `mismatches` to reduce off-target
matches. Increase primer `trim` to use more of the primer sequence
(more specific).

**Mitigation (developer):** Better seed evaluation — the current
`evaluate_seeds()` filters seeds by branching ratio and local
exploration size, but some off-target seeds pass because they extend
linearly into off-target regions. Read-backed seed evaluation
(`--read-eval`) helps by checking read divergence.

#### Seed evaluation threshold too stringent

The seed evaluation threshold (derived from median primer kmer count)
determines the minimum kmer count for bounded seed exploration. If the
median is inflated by off-target matches, real seeds fail to extend and
are falsely abandoned.

**Status:** Addressed in v3 by using median (not max) for seed eval
and max for extension thresholds. See DESIGN_DECISIONS.md "Which primer
count statistic to use matters."

#### Graph extension threshold too low

When extension thresholds start too low, the graph extends aggressively
from all seed nodes and hits the node budget before connecting forward
to reverse primer sites.

**Status:** Addressed in v3 by using max primer kmer count for
extension thresholds (starting high, stepping down).

#### Path finding budget exhausted

The DFS-based path finder has a state budget (`max_dfs_states`, default
100K). In complex graphs with many bubbles, the budget may be exhausted
before a valid start-to-end path is found.

**Mitigation (user):** Increase `--max-dfs-states` (hidden).

**Mitigation (developer):** Better edge ordering in DFS (coverage-
weighted, with bubble preferences) ensures the best paths are found
first. Already implemented.

## User-facing parameters and their effects

These are the knobs users can turn to improve recovery. Each has
trade-offs.

### Read count (`--max-reads`)

Controls how many reads from the input are used.

- **Increase:** More coverage, fewer gaps, recovers genes that fail at
  lower counts. Linear cost in runtime and memory for kmer counting.
- **Decrease:** Faster runtime, but genes with low coverage will fail.
  Useful for quick surveys or when only high-copy targets (mt, plastid)
  are needed.
- **Guidance:** Start with 1M for mt/plastid genes. Use 4-8M for
  nuclear genes in medium genomes. Use 16M+ for large genomes or
  single-copy targets.

### Kmer size (`-k`)

Controls the length of kmers used for counting and graph construction.

- **Increase (e.g., k=31):** More specific kmers, simpler graphs, but
  higher coverage needed per unique kmer. AT-rich regions may fragment.
- **Decrease (e.g., k=21):** Each kmer has higher coverage, bridging
  low-complexity gaps, but graphs have more branching and are larger.
  Below k=15, specificity is too low for most applications.
- **Default:** 31. Consider 21 for AT-rich organisms (Lepidoptera mt,
  some plastid genomes).

### Primer mismatches (`mismatches` in panel YAML)

Maximum allowed mismatches between primer and genome.

- **Increase:** Finds primers in more divergent species, but generates
  more primer kmer variants → more off-target seeds → larger graphs.
- **Decrease:** Fewer off-target matches, cleaner graphs, but may miss
  the target if the primer doesn't match exactly.
- **Default:** 2. Rarely needs changing.

### Primer trim (`trim` in panel YAML)

Number of 3' bases retained from the primer for kmer searching.

- **Increase:** More specific primer matching, fewer off-target seeds.
  But the trimmed primer must fit within a single kmer (trim ≤ k-1),
  and longer trim means fewer reads will contain the exact sequence.
- **Decrease:** More permissive matching, finds primers even with 3'
  divergence, but generates more off-target seeds.
- **Default:** 15. Increase for highly degenerate primers.

### Product length bounds (`min_length`, `max_length` in panel YAML)

Expected amplicon size range.

- **Too narrow:** May reject valid products that differ slightly from
  expected size (e.g., indels, species-specific length variation).
- **Too wide:** Allows spurious products from off-target amplification.
- **Guidance:** Set based on known amplicon sizes with ~20% margin.

### Read evaluation (`--read-eval`)

Enables read-backed seed evaluation (Pass 1 read retention).

- **On:** Filters off-target seeds more effectively by checking read
  divergence around seed nodes. Costs ~10% more memory (retained reads)
  and ~5% more time.
- **Off (default):** Faster, lower memory. Adequate for most cases.
- **Guidance:** Enable for degenerate primers or when seed explosion is
  suspected (many seeds abandoned in verbose output).

### Read threading (`--read-threading`)

Enables Pass 2 re-reading and read threading through the graph.

- **On:** Annotates edges with read support and phasing information.
  Enables bubble resolution. Costs a second read pass (I/O) and
  additional memory for retained reads.
- **Off (default):** Faster. Adequate when graph complexity is low.
- **Guidance:** Enable for genes with known heterozygosity or when
  multiple similar products are expected.

## Developer tuning parameters

These are hidden CLI arguments and hard-coded constants. They control
the algorithm's resource budgets and heuristic thresholds.

### Hidden CLI arguments

| Parameter | Default | Effect of increase | Effect of decrease |
| --- | ---: | --- | --- |
| `--max-nodes` | 50,000 | Allows larger graphs; may recover genes in complex regions but uses more memory and time | Smaller graphs; faster but may truncate before amplicon is connected |
| `--max-dfs-states` | 100,000 | Explores more paths; finds products in complex graphs | Faster path finding but may miss valid paths |
| `--max-paths-per-pair` | 20 | Reports more variant products | Fewer output sequences |
| `--max-node-visits` | 2 | Tolerates more cycles (tandem repeats) | Stricter cycle avoidance |
| `--max-primer-kmers` | 100 | Keeps more primer variants; better for degenerate primers | Fewer seeds; cleaner graphs |
| `--max-seed-nodes` | 500 | More thorough seed evaluation | Faster seed filtering |
| `--high-coverage-ratio` | 10.0 | Allows higher-coverage edges (less aggressive repeat filtering) | More aggressive repeat filtering |
| `--tip-coverage-fraction` | 0.1 | Prunes more tips (higher coverage threshold for keeping tips) | Preserves more tips |

### Hard-coded constants

| Constant | Value | Location | Controls |
| --- | ---: | --- | --- |
| `COVERAGE_MULTIPLIER` | 2 | `mod.rs` | Divisor for initial extension threshold: `primer_count / 2` |
| `COVERAGE_STEPS` | 4 | `mod.rs` | Number of threshold steps from initial to min_count |
| `MAX_NUM_AMPLICONS` | 20 | `paths.rs` | Hard limit on output FASTA records per gene |
| `EXTENSION_EVALUATION_FREQUENCY` | 1,000 | `graph.rs` | Graph size checked every N nodes |
| `MAX_BRANCHING_RATIO` | 0.4 | `seed_eval.rs` | Abandon seed if branching ratio exceeds this |
| `BUDGET_BRANCHING_THRESHOLD` | 0.2 | `seed_eval.rs` | Abandon seed if budget exhausted AND branching > this |
| `MIN_EXTENSION_FRACTION` | 0.1 | `seed_eval.rs` | Abandon seed if terminated with < 10% of expected nodes |
| `MAX_BUBBLE_DEPTH` | 50 | `bubble.rs` | Maximum depth for bubble detection |

## Real-world failure examples

Specific failures observed in benchmarks, with diagnosed root causes.
Full primer binding site data is in
[`primer_binding.yaml`](primer_binding.yaml).

### Insufficient read coverage

| Gene | Sample | Reads | Primer mm | Root cause |
| --- | --- | ---: | :---: | --- |
| 16S | Porites lutea | 1M | 0+0 | 542 Mb genome; ~0.4x mt coverage at 1M reads |
| CO1 | Xenia sp. | 1M | 0+0 | 223 Mb genome; insufficient mt coverage |
| ND4 | Drosophila | 1M | 0+0 | Primer kmer counts 3-11; **recovers at 4M** |
| ND4 | Heliconius | 1M | 0+0 | Reverse primer not found in reads at 1M |

### AT-rich kmer coverage gaps

| Gene | Sample | Reads | k | Root cause |
| --- | --- | ---: | ---: | --- |
| ND4 | Heliconius | 8M | 31 | AT-rich gap fragments graph; **recovers at k=21** |
| 16S | Drosophila | 8M | 31 | Severe AT-rich gap (83-90% AT); fails at k=21 and k=15 |

Graph evidence (k=31, 8M reads, exact primers):

- **Drosophila 16S:** 229 nodes. Forward component (146 nodes)
  terminates at `...CCCCAATAAAATATT` (83% AT). Backward component (83
  nodes) starts at `TTTTGACTAAAAAATAAAA...` (90% AT). Zero overlap.
- **Heliconius ND4:** 202 nodes. Forward component terminates after 1
  node — the primer's flanking region is deeply AT-rich and no 31-mer
  has count >= 2.

### Node budget exhaustion from off-target seeds

| Gene | Sample | Reads | Root cause |
| --- | --- | ---: | --- |
| ND4 | Drosophila | 1M (panel) | Degenerate Marquina primers; 50K node budget hit during reverse extension |

With exact (non-degenerate) species-specific primers, Drosophila ND4
recovers at 4M reads. The panel's degenerate primers create additional
off-target seeds that consume the node budget.

### Structural complexity

| Gene | Sample | Reads | Root cause |
| --- | --- | ---: | --- |
| ITS | Rhopilema | 1M | 0+0 primers; rDNA copy variation creates graph complexity |
| 28S | Drosophila | 8M | 1+0 mm; R1/R2 retrotransposon insertions in rDNA copies |
| ITS | Drosophila | 8M | 0+1 mm; amplicon ~4835 bp, too large for short-read graph |

### Primer mismatch failures

These are expected failures — the primers don't match the target
species well enough. Not a limitation of the method.

**Incompatible (>3 total mm or no binding site):**
- Drosophila: Fz4 (>5mm rev), 28S-v2 (5mm rev), ITS-v2/ITS-v3 (3mm
  shared rev site)
- Heliconius: CO2-v2 (7mm fwd), NADH (wrong gene entirely)
- Gryllus: CO2-v2 (4+4mm)
- Agalma: EF1A (3 consecutive 3' mm in fwd)

**Marginal (2-3 total mm, may work with tuning):**
- Porites CO1 (2+1mm), Drosophila NADH (2+2mm), Gpdh (2+0mm),
  Pgi (0+2mm)

**Gene absent in species:**
- Yp2 in Heliconius and Gryllus (Diptera-specific gene)

### Recovery statistics by gene type (1M reads, all benchmark samples)

| Gene type | Recovered | Total | Rate |
| --- | ---: | ---: | ---: |
| mitochondrial | 39 | 57 | 68% |
| rRNA | 45 | 100 | 45% |
| nuclear | 2 | 19 | 11% |
| plastid | 13 | 18 | 72% |
| **total** | **99** | **194** | **51%** |

Note: these rates include genes with known primer mismatches and
genes absent from the target species. Excluding those, recovery rates
for genes with compatible primers are substantially higher.

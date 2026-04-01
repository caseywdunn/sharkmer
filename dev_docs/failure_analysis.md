# Benchmark failure analysis

Primer binding site analysis for genes that fail to amplify in sPCR
benchmarks. Actual binding sequences in each species were identified from
publicly available sequences (mitochondrial genomes, rRNA, mRNA) and
compared to the panel primer sequences. Full binding site data is in
[`primer_binding.yaml`](primer_binding.yaml).

Failures at 16M reads (highest coverage tested) unless noted otherwise.

## Cnidaria panel

| Gene | Species | Fwd mm | Rev mm | Accession | Verdict |
| --- | --- | --- | --- | --- | --- |
| 16S | Porites lutea | 0 | 0 | KU159432.1 | Not primer mismatch — likely coverage (large 542 Mb genome) |
| CO1 | Porites lutea | 2 | 1 | KU159432.1 | 3 total mm; COX1 is a split gene with intron |
| CO1 | Agalma elegans | 0 | 1 | NC_080954.1 | 1mm at 3' critical position (pos 13/15); paradoxically works at lower coverages |
| EF1A | Agalma elegans | 3 | 0 | DQ157428.1 | 3 consecutive 3' mismatches in forward primer block extension |
| ITS | Rhopilema esculentum | 0 | 0 | JX845352.1 / XR_010517148.1 | Not primer mismatch — likely rDNA copy variation / graph complexity |

## Insecta panel — Drosophila melanogaster

| Gene | Fwd mm | Rev mm | Accession | Verdict |
| --- | --- | --- | --- | --- |
| 12S | 0 | 0 | KJ947872.2 | Not primer mismatch — works at 2M but fails at 16M (graph complexity?) |
| 16S | 0 | 0 | KJ947872.2 | Not primer mismatch |
| 16S-v2 | 1 | 0 | KJ947872.2 | 1mm internal, should not prevent amplification |
| ND1 | 0 | 0 | KJ947872.2 | Not primer mismatch — works at 2M-4M but fails at 16M |
| ND4 | 0 | 0 | KJ947872.2 | Not primer mismatch |
| NADH | 2 | 2 | KJ947872.2 | Actually targets ND2, not ND5; 4 total mm (all internal) |
| 28S | 1 | 0 | M21017.1 | 1mm marginal; R1/R2 retrotransposon insertions complicate graph |
| ITS | 0 | 1 | M21017.1 | Amplicon ~4835 bp — too large for de Bruijn graph from short reads |
| EF1g | 3 | 1 | NM_143743.3 | Hylaeus bee primers, too divergent for Diptera |
| Fz4 | 4 | >5 | NM_078513.3 | No reverse binding site — completely incompatible |
| Gpdh | 2 | 0 | NM_057219.4 | Marginal; may work with relaxed conditions |
| Pgi | 0 | 2 | NM_078939.3 | Marginal; forward is exact match |
| 28S-v2 | 3 | 5 | M21017.1 | Reverse non-binding; Drosophila 28S too diverged for Evans primers |
| ITS-v2 | 0 | 3 | M21017.1 | Reverse shares 28S-v2 binding site (3mm); amplicon would be ~1567 bp |
| ITS-v3 | 1 | 3 | M21017.1 | Same reverse primer problem as ITS-v2; amplicon ~1420 bp |

## Insecta panel — Heliconius pachinus

| Gene | Fwd mm | Rev mm | Accession | Verdict |
| --- | --- | --- | --- | --- |
| 16S-v2 | 1 | 3 | NC_024741.1 | Drosophila-specific; reverse too divergent for Lepidoptera |
| CO2-v2 | 7 | 0 | NC_024741.1 | Forward completely incompatible (only first 6 bases match) |
| ND4 | 0 | 0 | NC_024741.1 | Not primer mismatch — degenerate primer kmer explosion or AT-rich graph issues |
| NADH | 2 | 4 | NC_024741.1 | Primers bind tRNA-Met/ND2, not near ND5 at all |
| Yp2 | — | — | — | Gene does not exist in Lepidoptera (Diptera-specific) |

No public sequences available for: EF1g, Fz4, Gpdh, Pgi (no Heliconius
orthologs sequenced), 28S, 28S-v2, ITS, ITS-v2, ITS-v3 (no standalone
GenBank references).

## Insecta panel — Gryllus bimaculatus

| Gene | Fwd mm | Rev mm | Accession | Verdict |
| --- | --- | --- | --- | --- |
| CO2-v2 | 4 | 4 | MT993975.1 | Drosophila-specific; both primers exceed mismatch threshold |
| Yp2 | — | — | — | Diptera-specific gene; ortholog (vitellogenin) too divergent |

No public sequences available for: 28S-v2, ITS-v2, ITS-v3 (only partial
references), EF1g, Fz4, Gpdh, Pgi (no Gryllus gene sequences).

## Current status (post Phases 3-7, commit df4269e)

Updated analysis of perfect-match failures using current code (Phase 3
graph traversal improvements + Phase 7 read-backed seed evaluation).
Tested at multiple coverage levels with and without `--read-eval`.

### Recovery summary

| Gene | Species | Old status | 1M | 2M | 4M | 4M --read-eval | 8M | 16M | 16M --read-eval |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16S | Porites | Fail 4M | - | - | - | - | - | - | - |
| CO1 | Porites | Fail 4M | - | - | - | - | - | - | - |
| 12S | Drosophila | Fail 16M | Y | Y | Y | Y | Y | Y | Y |
| 16S | Drosophila | Fail 16M | - | - | - | - | - | - | - |
| 16S-v2 | Drosophila | Fail 16M | - | - | - | - | - | - | - |
| 28S | Drosophila | Fail 16M | - | - | - | - | - | - | - |
| ND1 | Drosophila | Fail 16M | Y | Y | Y | Y | - | Y | Y |
| ND4 | Drosophila | Fail 16M | - | - | Y | Y | - | - | - |
| ITS | Drosophila | Fail 16M | - | - | - | **Y** | - | - | - |
| ITS | Rhopilema | Fail 16M | - | - | - | - | - | - | - |
| ND4 | Heliconius | Fail 16M | - | - | Y | - | - | Y | Y |

Key: Y = recovered, - = still fails, **Y** = recovered only with --read-eval

### Analysis of changes since original failure analysis

**Newly recovered genes:**
- **Drosophila 12S**: Was failing at 16M in original analysis but now
  recovers at all coverage levels 1M-16M. Phase 3 graph traversal
  improvements (bidirectional extension, coverage-weighted DFS) resolved
  the graph complexity issue at high coverage.
- **Drosophila ND1**: Similar to 12S — was failing at 16M, now recovers
  at 1M-4M and 16M. Drops out at 8M, suggesting a coverage-dependent
  graph complexity sweet spot.
- **Drosophila ND4**: Recovers at 4M (both with and without --read-eval)
  and at some higher coverages. Previously failed at all levels.
- **Drosophila ITS**: Recovers at 4M **only with --read-eval**. The
  read-backed seed evaluation rejects off-target seeds that were consuming
  graph budget, allowing the correct seeds to extend successfully. This is
  the first gene where `--read-eval` makes a difference.
- **Heliconius ND4**: Now recovers at 4M and 16M. Previously failed at
  all levels.

**Still failing:**
- **Porites 16S, CO1 at 4M**: Still fail. Porites has a 542 Mb genome —
  at 4M 150bp reads, local coverage is only ~1x. The target regions
  likely don't have enough kmer coverage to build connected graphs.
- **Drosophila 16S, 16S-v2, 28S**: Still fail at all levels. These use
  degenerate Marquina primers that seed 100+ kmer nodes per primer due
  to ambiguity codes. The non-specific seeds create graphs with many
  false start/end nodes. ND1 recovery suggests some similar genes can
  be rescued by Phase 3 improvements, but 16S/28S remain too complex.
- **Rhopilema ITS at 16M**: Still fails. This is likely rDNA copy number
  variation creating graph complexity that exceeds node/DFS budgets.

**Coverage-dependent instability:**
- **Drosophila ND1**: Recovers at 1M-4M and 16M but not 8M. This
  suggests a narrow window where graph complexity is just right —
  too low and there's insufficient coverage, too high and repeat
  regions overwhelm the graph, but at very high coverage the
  threshold stepping finds a good balance.
- **Drosophila ND4**: Only at 4M. Similar coverage sensitivity.

### Effect of --read-eval

`--read-eval` had a clear positive effect in one case:
- **Drosophila ITS at 4M**: 11 → 12 genes. The ITS amplicon (~4835 bp)
  is at the edge of what de Bruijn graphs from short reads can handle.
  Read-eval rejects off-target primer seeds that waste the node budget,
  leaving more budget for the real ITS graph to extend fully.

No regressions from `--read-eval` were observed at any coverage level
for any species tested.

## Failure categories

### Perfect-match failures (8 entries)

Primers bind perfectly (0-1mm) but sPCR still fails. Diagnostic runs
with `-v` confirmed that **all 8 are graph traversal failures, not
primer-finding failures.** In every case, primer kmers are found in the
reads, a graph is seeded, but no path is found from forward to reverse
primer binding sites.

None of these failures are caused by insufficient read coverage to find
primers. The Drosophila 16S-v2 (1mm internal) and 28S (1mm) are
near-perfect matches that exhibit the same graph traversal failure
pattern and are included in the analysis below.

| Gene | Species | Reads | Fwd seeds | Rev seeds | Primer cov | Failure mode |
| --- | --- | --- | --- | --- | --- | --- |
| 16S | Porites lutea | 4M | 144 | 108 | 122 | Graph traversal, all 4 thresholds |
| CO1 | Porites lutea | 4M | 151 | 104 | 395 | Graph traversal, all 4 thresholds |
| ITS | Rhopilema esculentum | 16M | 45 | 102 | 966 | Graph traversal, all 4 thresholds |
| 12S | Drosophila melanogaster | 16M | 112 | 103 | 330 | Graph traversal, all 4 thresholds |
| 16S | Drosophila melanogaster | 16M | 118 | 125 | 162 | Graph traversal, all 4 thresholds |
| 16S-v2 | Drosophila melanogaster | 16M | 109 | 109 | 110 | End node found at min-count 2, but no valid path |
| 28S | Drosophila melanogaster | 16M | 100 | 113 | 980 | End node found at min-count 490, but no valid path |
| ND1 | Drosophila melanogaster | 16M | 106 | 102 | 724 | End node found at min-count 242, but no valid path |
| ND4 | Drosophila melanogaster | 16M | 104 | 100 | 957 | Graph traversal, all 4 thresholds |
| ND4 | Heliconius pachinus | 16M | 101 | 107 | 188 | Graph traversal, all 4 thresholds |

Key observations:

- **Non-specific primer seeding.** The degenerate Marquina primers
  (12S, 16S, ND1, ND4) seed 100-280 kmer nodes per primer due to
  ambiguity codes (H, D, Y, N, R). `MAX_NUM_PRIMER_KMERS` is 100, so
  many of these are at or above the cap. Most seed kmers are non-specific
  (from unrelated genomic regions), creating a graph with many false
  start/end nodes that don't connect.

- **End nodes found but no path.** For ND1, 16S-v2, and 28S in
  Drosophila, the log shows "End node incorporated into graph, complete
  PCR product found" — meaning forward and reverse primer regions are
  connected in the graph — but the path-finding step still fails to
  extract a valid product. This points to graph pruning or path
  enumeration issues rather than graph connectivity.

- **Coverage-dependent failures.** Drosophila 12S and ND1 amplify at
  2M reads but fail at 16M. At 2M, 12S seeds 278 forward + 150 reverse
  nodes (more non-specific seeds than at 16M) but has lower observed
  primer coverage (42 vs 330). The lower thresholds at 2M may allow the
  correct path to survive pruning, while at 16M the higher coverage
  creates a denser, more complex graph that obscures the target path.

- **Porites CO1 is also graph traversal.** Despite having 2+1 primer
  mismatches, Porites CO1 primers are still found (151 + 104 seed nodes,
  coverage 395) — the failure is in graph traversal, not primer binding.

### Primer mismatch failures (13 entries)

Range from marginal (1-2mm total) to completely incompatible:

- **Incompatible (>3mm or no binding site):** Drosophila Fz4 (>5mm rev),
  28S-v2 (5mm rev), ITS-v2/ITS-v3 (3mm shared rev site); Heliconius
  CO2-v2 (7mm fwd), NADH (wrong gene); Gryllus CO2-v2 (4mm each);
  Agalma EF1A (3 consecutive 3' mm)
- **Marginal (1-2mm, may work with tuning):** Porites CO1 (2+1mm),
  Drosophila NADH (2+2mm), Gpdh (2mm fwd), Pgi (2mm rev), 28S (1mm fwd)

Key pattern: the ITS-v2, ITS-v3, and 28S-v2 reverse primers all target
the same 5' region of 28S that has 3mm in Drosophila — a single point of
divergence causes three primer pairs to fail.

### Gene absent or no reference data (6 entries)

- **Gene absent:** Yp2 in Heliconius and Gryllus (Diptera-specific gene)
- **No public sequence:** Magnacca bee primers (EF1g, Fz4, Gpdh, Pgi) and
  rRNA markers in Heliconius/Gryllus lack standalone GenBank entries to
  check against

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

## Failure categories

### Perfect-match failures (8 entries)

Primers bind perfectly but sPCR still fails. These point to software or
coverage issues rather than biological incompatibility:

- Porites lutea 16S, Rhopilema esculentum ITS
- Drosophila melanogaster 12S, 16S, ND1, ND4
- Heliconius pachinus ND4

Several of these paradoxically work at lower read depths (1M-4M) but fail
at higher (8M-16M), suggesting coverage-dependent graph complexity. The
Drosophila 16S-v2 (1mm internal) and 28S (1mm) are near-perfect matches
that should also amplify.

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

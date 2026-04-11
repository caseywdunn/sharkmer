# Panel validation: angiospermae

- **Panel version**: `1.0.0`
- **sharkmer version**: `3.0.0-rc`
- **Date**: 2026-04-11 05:37:24
- **Machine**: Linux 6.12.76-linuxkit, 12 cores, 47.0 GB RAM

## Liriodendron tulipifera (SRR25378184)

| Gene | 1000k | Score |
|------|---:|:-----:|
| psbA-trnH | --- | `-*-` |
| rpl36-infA-rps8 | 552bp (100.0%) | `+**` |
| trnK-rps16 | --- | `---` |
| trnV-atpE | --- | `-*-` |
| trnC-ycf6 | 1099bp (100.0%) | `+**` |
| ycf6-psbM | 927bp (100.0%) | `+**` |
| psbM-trnD | 1287bp (100.0%) | `+**` |
| atpB-rbcL | 812bp (79.0%) | `+++` |
| trnL-F | 434bp (100.0%) | `+**` |
| ITS | --- | `-+-` |

Wall times: 1000k: 10.2s

## Acer monspessulanum (ERR14009273)

| Gene | 1000k | Score |
|------|---:|:-----:|
| psbA-trnH | 508bp (100.0%) | `+**` |
| rpl36-infA-rps8 | --- | `-+-` |
| trnK-rps16 | --- | `---` |
| trnV-atpE | --- | `-*-` |
| trnC-ycf6 | 593bp (100.0%) | `+**` |
| ycf6-psbM | 937bp (100.0%) | `+**` |
| psbM-trnD | 699bp (86.3%) | `+++` |
| atpB-rbcL | 812bp (100.0%) | `+**` |
| trnL-F | --- | `-*-` |
| ITS | 753bp (100.0%) | `+**` |

Wall times: 1000k: 9.0s

## Cross-sample summary (highest depth)

| Gene | Liriodendron tuli... | Acer monspessulanum |
|------|:---:|:---:|
| psbA-trnH | `-*-` | `+**` |
| rpl36-infA-rps8 | `+**` | `-+-` |
| trnK-rps16 | `---` | `---` |
| trnV-atpE | `-*-` | `-*-` |
| trnC-ycf6 | `+**` | `+**` |
| ycf6-psbM | `+**` | `+**` |
| psbM-trnD | `+**` | `+++` |
| atpB-rbcL | `+++` | `+**` |
| trnL-F | `+**` | `-*-` |
| ITS | `-+-` | `+**` |

**Scoring** — three positions: recovery / reference availability / BLAST result.

| Code | Meaning |
|------|---------|
| `+**` | Recovered, confirmed: same gene, same species |
| `+*+` | Recovered, ref for this species exists but hit different species |
| `+*-` | Recovered, ref for this species exists but no BLAST hit (suspicious) |
| `+++` | Recovered, hit same gene in a different species |
| `++-` | Recovered, refs for other species exist but no BLAST hit |
| `+--` | Recovered, no references for this gene |
| `-*-` | Not recovered, ref exists for this species |
| `-+-` | Not recovered, refs exist for other species |
| `---` | Not recovered, no references for this gene |

Position 1: `-` no product, `+` product recovered. Position 2: `-` no reference for any species for this gene, `+` reference for other species, `*` reference for this species. Position 3: `-` no BLAST hit to same gene, `+` hit same gene different species, `*` same gene same species.


## Primer binding analysis

For each gene, the first and last `trim` bases of each recovered amplicon are compared against the user-specified primer sequence (3'-trimmed to the match window). Reverse primer bindings are shown reverse-complemented so they appear in the same orientation as the primer was written in the panel. Use this section to decide whether a primer's degeneracy should be reduced (only a subset of coded bases is actually seen) or widened (an off-code base was absorbed by sharkmer's `--mismatches` tolerance).

### psbA-trnH

**Not recovered** in: SRR25378184. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `GTTATGCATGAACGTAATGCTC` (trim=15, match window `ATGAACGTAATGCTC`)

```
  spec          ATGAACGTAATGCTC
                |||||||||||||||
  ERR14009273   ATGAACGTAATGCTC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `CGCGCATGGTGGATTCACAATCC` (trim=15, match window `GTGGATTCACAATCC`)

```
  spec          GTGGATTCACAATCC
                |||||||||||||||
  ERR14009273   GTGGATTCACAATCC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### rpl36-infA-rps8

**Not recovered** in: ERR14009273. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `CACAAATTTTACGAACGAAG` (trim=15, match window `ATTTTACGAACGAAG`)

```
  spec          ATTTTACGAACGAAG
                |||||||||||xx||
  SRR25378184   ATTTTACGAACAGAG
```

  pos 12: G ({G}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos 13: A ({A}) observed {G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `TAATGACAGAYCGAGARGCTCGAC` (trim=15, match window `AYCGAGARGCTCGAC`)

```
  spec          AYCGAGARGCTCGAC
                |.||x||.|||||||
  SRR25378184   ATCGGGAAGCTCGAC
```

  pos  2: Y ({C,T}) partially utilised — observed {T}; could reduce to T
  pos  5: A ({A}) observed {G} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos  8: R ({A,G}) partially utilised — observed {A}; could reduce to A

### trnK-rps16

**Not recovered** in: SRR25378184, ERR14009273. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

### trnV-atpE

**Not recovered** in: SRR25378184, ERR14009273. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

### trnC-ycf6

**Forward primer** — spec `CCAGTTCAAATCTGGGTGTC` (trim=15, match window `TCAAATCTGGGTGTC`)

```
  spec          TCAAATCTGGGTGTC
                |||||||x|||||||
  SRR25378184   TCAAATCCGGGTGTC
  ERR14009273   TCAAATCTGGGTGTC
```

  pos  8: T ({T}) observed {C,T} — off-code base(s) absorbed by --mismatches; consider widening to Y

**Reverse primer** — spec `CCCAAGCAAGACTTACTATATCC` (trim=15, match window `AGACTTACTATATCC`)

```
  spec          AGACTTACTATATCC
                |||||||||||||||
  SRR25378184   AGACTTACTATATCC
  ERR14009273   AGACTTACTATATCC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### ycf6-psbM

**Forward primer** — spec `GGATATAGTAAGTCTTGCTTGGG` (trim=15, match window `TAAGTCTTGCTTGGG`)

```
  spec          TAAGTCTTGCTTGGG
                |||||||x|||||||
  SRR25378184   TAAGTCTCGCTTGGG
  ERR14009273   TAAGTCTCGCTTGGG
```

  pos  8: T ({T}) observed {C} — off-code base(s) absorbed by --mismatches; consider widening to Y

**Reverse primer** — spec `TTCTTGCATTTATTGCTACTGC` (trim=15, match window `ATTTATTGCTACTGC`)

```
  spec          ATTTATTGCTACTGC
                |||||||||||||||
  SRR25378184   ATTTATTGCTACTGC
  ERR14009273   ATTTATTGCTACTGC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### psbM-trnD

**Forward primer** — spec `GCGGTAGGAACTAGAATAAATAG` (trim=15, match window `AACTAGAATAAATAG`)

```
  spec          AACTAGAATAAATAG
                |||||||||x||x||
  SRR25378184   AACTAGAATGAACAG
  ERR14009273   AACTAGAATGAACAG
```

  pos 10: A ({A}) observed {G} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos 13: T ({T}) observed {C} — off-code base(s) absorbed by --mismatches; consider widening to Y

**Reverse primer** — spec `GGGATTGTAGTTCAATTGGT` (trim=15, match window `TGTAGTTCAATTGGT`)

```
  spec          TGTAGTTCAATTGGT
                |||||||||||x|||
  SRR25378184   TGTAGTTCAATCGGT
  ERR14009273   TGTAGTTCAATTGGT
```

  pos 12: T ({T}) observed {C,T} — off-code base(s) absorbed by --mismatches; consider widening to Y

### atpB-rbcL

**Forward primer** — spec `AGAAGTAGTAGGATTGATTCTCATA` (trim=15, match window `GGATTGATTCTCATA`)

```
  spec          GGATTGATTCTCATA
                |||||x|||||||||
  SRR25378184   GGATTGATTCTCATA
  ERR14009273   GGATTTATTCTCATA
```

  pos  6: G ({G}) observed {G,T} — off-code base(s) absorbed by --mismatches; consider widening to K

**Reverse primer** — spec `GAATCCAACACTTGCTTTAGTCTCT` (trim=15, match window `CTTGCTTTAGTCTCT`)

```
  spec          CTTGCTTTAGTCTCT
                |||||||||||||||
  SRR25378184   CTTGCTTTAGTCTCT
  ERR14009273   CTTGCTTTAGTCTCT
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### trnL-F

**Not recovered** in: ERR14009273. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `GGTTCAAGTCCCTCTATCCC` (trim=15, match window `AAGTCCCTCTATCCC`)

```
  spec          AAGTCCCTCTATCCC
                |||||||||||||||
  SRR25378184   AAGTCCCTCTATCCC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `ATTTGAACTGGTGACACGAG` (trim=15, match window `AACTGGTGACACGAG`)

```
  spec          AACTGGTGACACGAG
                |||||||||||||||
  SRR25378184   AACTGGTGACACGAG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### ITS

**Not recovered** in: SRR25378184. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `CCTTATCATTTAGAGGAAGGAG` (trim=15, match window `ATTTAGAGGAAGGAG`)

```
  spec          ATTTAGAGGAAGGAG
                |||||||||||||||
  ERR14009273   ATTTAGAGGAAGGAG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `TCCTCCGCTTATTGATATGC` (trim=15, match window `CGCTTATTGATATGC`)

```
  spec          CGCTTATTGATATGC
                |||||||||||||||
  ERR14009273   CGCTTATTGATATGC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

## Reference match details

| Sample | Gene | Sample taxon | Ref taxon | Ref accession | Identity | Align len |
|--------|------|-------------|-----------|---------------|----------|-----------|
| SRR25378184 | atpB-rbcL | Liriodendron tulipifera | Euphorbia nicaeensis subsp. nicaeensis | OR400606 | 79.0% | 861 |
| SRR25378184 | psbM-trnD | Liriodendron tulipifera | **same** | NC 008326 | 100.0% | 1287 |
| SRR25378184 | rpl36-infA-rps8 | Liriodendron tulipifera | **same** | NC 008326 | 100.0% | 552 |
| SRR25378184 | trnC-ycf6 | Liriodendron tulipifera | **same** | NC 008326 | 100.0% | 1099 |
| SRR25378184 | trnL-F | Liriodendron tulipifera | **same** | NC 008326 | 100.0% | 434 |
| SRR25378184 | ycf6-psbM | Liriodendron tulipifera | **same** | NC 008326 | 100.0% | 927 |
| ERR14009273 | ITS | Acer monspessulanum | **same** | MW070149 | 100.0% | 753 |
| ERR14009273 | atpB-rbcL | Acer monspessulanum | **same** | NC 056221 | 100.0% | 812 |
| ERR14009273 | psbA-trnH | Acer monspessulanum | **same** | NC 056221 | 100.0% | 508 |
| ERR14009273 | psbM-trnD | Acer monspessulanum | Liriodendron tulipifera | NC 008326 | 86.3% | 299 |
| ERR14009273 | trnC-ycf6 | Acer monspessulanum | **same** | NC 056221 | 100.0% | 593 |
| ERR14009273 | ycf6-psbM | Acer monspessulanum | **same** | NC 056221 | 100.0% | 937 |

## Performance

| Sample | Max reads | Wall time | Peak memory | Reads ingested | Bases ingested | Distinct kmers |
|--------|----------:|----------:|------------:|---------------:|---------------:|---------------:|
| Liriodendron tulipifera | 1000k | 10.2s | 2.1 GB | 656,431 | 99,121,081 | 87,286,422 |
| Acer monspessulanum | 1000k | 9.0s | 2.1 GB | 1,000,000 | 99,211,896 | 81,206,414 |

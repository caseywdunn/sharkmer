# Panel validation: cnidaria

- **Panel version**: `1.1.0`
- **sharkmer version**: `3.0.0-rc`
- **Date**: 2026-04-11 05:47:26
- **Machine**: Linux 6.12.76-linuxkit, 12 cores, 47.0 GB RAM

## Xenia sp. (SRR9278435)

| Gene | 1000k | 2000k | 4000k | 8000k | Score |
|------|---:|---:|---:|---:|:-----:|
| 16S | 621bp (100.0%) | 621bp (100.0%) | 621bp (100.0%) | 621bp (100.0%) | `+**` |
| CO1 | 692bp (100.0%) | 692bp (100.0%) | 692bp (100.0%) | 692bp (100.0%) | `+**` |
| 18S | 1807bp (100.0%) | 1807bp (100.0%) | 1807bp (100.0%) | 1807bp (100.0%) | `+**` |
| 28S | 3234bp (100.0%) | 3234bp (100.0%) | 3234bp (100.0%) | 3234bp (100.0%) | `+**` |
| ITS | 899bp (99.9%) | 899bp (99.8%) | 899bp (99.8%) | 899bp (99.8%) | `+**` |
| ITS-v2 | 748bp (99.9%) | 748bp (100.0%) | 748bp (100.0%) | 748bp (100.0%) | `+**` |
| 28S-v2 | 414bp (100.0%) | 414bp (100.0%) | 414bp (100.0%) | 414bp (100.0%) | `+**` |

Wall times: 1000k: 20.8s, 2000k: 35.6s, 4000k: 59.8s, 8000k: 106.1s

## Agalma elegans (SRR25099394)

| Gene | 1000k | 2000k | 4000k | 8000k | Score |
|------|---:|---:|---:|---:|:-----:|
| 16S | --- | --- | 542bp (100.0%) | 542bp (100.0%) | `+**` |
| CO1 | --- | 798bp | 692bp (100.0%) | 692bp (100.0%) | `+**` |
| 18S | 1785bp (100.0%) | 1785bp (100.0%) | 1785bp (100.0%) | 1785bp (100.0%) | `+**` |
| 28S | 3230bp (100.0%) | 3230bp (100.0%) | 3230bp (100.0%) | 3230bp (100.0%) | `+**` |
| ITS | 817bp (100.0%) | 817bp (100.0%) | 817bp (100.0%) | 817bp (100.0%) | `+**` |
| ITS-v2 | 671bp (100.0%) | 671bp (100.0%) | 671bp (100.0%) | 671bp (100.0%) | `+++` |
| 28S-v2 | 417bp (100.0%) | 417bp (100.0%) | 417bp (100.0%) | 417bp (100.0%) | `+*+` |

Wall times: 1000k: 13.3s, 2000k: 26.4s, 4000k: 51.2s, 8000k: 93.5s

## Rhopilema esculentum (SRR8617500)

| Gene | 1000k | 2000k | 4000k | 8000k | Score |
|------|---:|---:|---:|---:|:-----:|
| 16S | --- | 563bp (100.0%) | 563bp (100.0%) | 563bp (100.0%) | `+**` |
| CO1 | 692bp (100.0%) | 692bp (100.0%) | 692bp (100.0%) | 692bp (100.0%) | `+**` |
| 18S | 1796bp (100.0%) | 1796bp (100.0%) | 1796bp (100.0%) | 1796bp (100.0%) | `+**` |
| 28S | 3258bp (100.0%) | 3258bp (100.0%) | 3258bp (100.0%) | 3258bp (100.0%) | `+**` |
| ITS | --- | --- | --- | --- | `-+-` |
| ITS-v2 | 929bp (100.0%) | 929bp (100.0%) | 929bp (100.0%) | 929bp (100.0%) | `+**` |
| 28S-v2 | 426bp (100.0%) | 426bp (100.0%) | 426bp (100.0%) | 426bp (100.0%) | `+**` |

Wall times: 1000k: 14.1s, 2000k: 28.7s, 4000k: 52.1s, 8000k: 94.5s

## Cross-sample summary (highest depth)

| Gene | Xenia sp. | Agalma elegans | Rhopilema esculentum |
|------|:---:|:---:|:---:|
| 16S | `+**` | `+**` | `+**` |
| CO1 | `+**` | `+**` | `+**` |
| 18S | `+**` | `+**` | `+**` |
| 28S | `+**` | `+**` | `+**` |
| ITS | `+**` | `+**` | `-+-` |
| ITS-v2 | `+**` | `+++` | `+**` |
| 28S-v2 | `+**` | `+*+` | `+**` |

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

### 16S

**Forward primer** — spec `GRCTGTTTACCAAAAACATA` (trim=15, match window `TTTACCAAAAACATA`)

```
  spec          TTTACCAAAAACATA
                |||||||||||||||
  SRR9278435    TTTACCAAAAACATA
  SRR25099394   TTTACCAAAAACATA
  SRR8617500    TTTACCAAAAACATA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `AATTCAACATMGAGG` (trim=15, match window `AATTCAACATMGAGG`)

```
  spec          AATTCAACATMGAGG
                ||||||||||.||||
  SRR9278435    AATTCAACATCGAGG
  SRR25099394   AATTCAACATAGAGG
  SRR8617500    AATTCAACATCGAGG
```

  1 degenerate position(s) fully utilised across samples.
  pos 11: M ({A,C}) fully utilised — observed {A,C}

### CO1

**Forward primer** — spec `WAAYCATAAAGATAT` (trim=18, match window `WAAYCATAAAGATAT`)

```
  spec          WAAYCATAAAGATAT
                x||.||x||x|||||
  SRR9278435    TAATCATAAGGATAT
  SRR25099394   GAATCATAAAGATAT
  SRR8617500    TAACCACAAAGATAT
```

  pos  1: W ({A,T}) observed {G,T} — off-code base(s) absorbed by --mismatches; consider widening to D
  pos  4: Y ({C,T}) fully utilised — observed {C,T}
  pos  7: T ({T}) observed {C,T} — off-code base(s) absorbed by --mismatches; consider widening to Y
  pos 10: A ({A}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `GGRTGMCCAAAAAACCARA` (trim=18, match window `GRTGMCCAAAAAACCARA`)

```
  spec          GRTGMCCAAAAAACCARA
                |.||x||||||||x||.|
  SRR9278435    GGTGTCCAAAAAACCAAA
  SRR25099394   GATGACCAAAAAATCAAA
  SRR8617500    GATGTCCAAAAAATCAAA
```

  pos  2: R ({A,G}) fully utilised — observed {A,G}
  pos  5: M ({A,C}) observed {A,T} — off-code base(s) absorbed by --mismatches; consider widening to H
  pos 14: C ({C}) observed {C,T} — off-code base(s) absorbed by --mismatches; consider widening to Y
  pos 17: R ({A,G}) partially utilised — observed {A}; could reduce to A

### 18S

**Forward primer** — spec `AACCTGGTTGATCCTGCCAGT` (trim=15, match window `GTTGATCCTGCCAGT`)

```
  spec          GTTGATCCTGCCAGT
                |||||||||||||||
  SRR9278435    GTTGATCCTGCCAGT
  SRR25099394   GTTGATCCTGCCAGT
  SRR8617500    GTTGATCCTGCCAGT
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `TGATCCTTCTGCAGGTTCACCTAC` (trim=15, match window `TGCAGGTTCACCTAC`)

```
  spec          TGCAGGTTCACCTAC
                x||||||||||||||
  SRR9278435    CGCAGGTTCACCTAC
  SRR25099394   CGCAGGTTCACCTAC
  SRR8617500    CGCAGGTTCACCTAC
```

  pos  1: T ({T}) observed {C} — off-code base(s) absorbed by --mismatches; consider widening to Y

### 28S

**Forward primer** — spec `CCYYAGTAACGGCGAGT` (trim=15, match window `YYAGTAACGGCGAGT`)

```
  spec          YYAGTAACGGCGAGT
                ..|||||x|||||x|
  SRR9278435    CTAGTAATGGCGAAT
  SRR25099394   CTAGTAACGGCGAGT
  SRR8617500    CCAGTAACGGCGAGT
```

  pos  1: Y ({C,T}) partially utilised — observed {C}; could reduce to C
  pos  2: Y ({C,T}) fully utilised — observed {C,T}
  pos  8: C ({C}) observed {C,T} — off-code base(s) absorbed by --mismatches; consider widening to Y
  pos 14: G ({G}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `SWACAGATGGTAGCTTCG` (trim=15, match window `CAGATGGTAGCTTCG`)

```
  spec          CAGATGGTAGCTTCG
                |||||||||||||||
  SRR9278435    CAGATGGTAGCTTCG
  SRR25099394   CAGATGGTAGCTTCG
  SRR8617500    CAGATGGTAGCTTCG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### ITS

**Not recovered** in: SRR8617500. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `TACACACCGCCCGTCGCTACTA` (trim=15, match window `CGCCCGTCGCTACTA`)

```
  spec          CGCCCGTCGCTACTA
                |||||||||||||||
  SRR9278435    CGCCCGTCGCTACTA
  SRR25099394   CGCCCGTCGCTACTA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `ACTCGCCGTTACTRRGG` (trim=15, match window `TCGCCGTTACTRRGG`)

```
  spec          TCGCCGTTACTRRGG
                |||||x|||||..||
  SRR9278435    TCGCCATTACTAGGG
  SRR25099394   TCGCCGTTACTAGGG
```

  pos  6: G ({G}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos 12: R ({A,G}) partially utilised — observed {A}; could reduce to A
  pos 13: R ({A,G}) partially utilised — observed {G}; could reduce to G

### ITS-v2

**Forward primer** — spec `GTAGGTGAACCTGCAGAAGGATCA` (trim=15, match window `CCTGCAGAAGGATCA`)

```
  spec          CCTGCAGAAGGATCA
                |||||x|||||||||
  SRR9278435    CCTGCGGAAGGATCA
  SRR25099394   CCTGCGGAAGGATCA
  SRR8617500    CCTGCGGAAGGATCA
```

  pos  6: A ({A}) observed {G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `ACTCGCCGTTACTRRGG` (trim=15, match window `TCGCCGTTACTRRGG`)

```
  spec          TCGCCGTTACTRRGG
                |||||x|||||..||
  SRR9278435    TCGCCATTACTAGGG
  SRR25099394   TCGCCGTTACTAGGG
  SRR8617500    TCGCCGTTACTGGGG
```

  pos  6: G ({G}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos 12: R ({A,G}) fully utilised — observed {A,G}
  pos 13: R ({A,G}) partially utilised — observed {G}; could reduce to G

### 28S-v2

**Forward primer** — spec `CGTGAAACCGYTRRAAGGG` (trim=15, match window `AAACCGYTRRAAGGG`)

```
  spec          AAACCGYTRRAAGGG
                ||||||.|..xx|||
  SRR9278435    AAACCGTTGAAAGGG
  SRR25099394   AAACCGTTAGGAGGG
  SRR8617500    AAACCGTTGGAGGGG
```

  pos  7: Y ({C,T}) partially utilised — observed {T}; could reduce to T
  pos  9: R ({A,G}) fully utilised — observed {A,G}
  pos 10: R ({A,G}) fully utilised — observed {A,G}
  pos 11: A ({A}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R
  pos 12: A ({A}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `TTGGTCCGTGTTTCAAGACG` (trim=15, match window `CCGTGTTTCAAGACG`)

```
  spec          CCGTGTTTCAAGACG
                |||||||||||||||
  SRR9278435    CCGTGTTTCAAGACG
  SRR25099394   CCGTGTTTCAAGACG
  SRR8617500    CCGTGTTTCAAGACG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

## Reference match details

| Sample | Gene | Sample taxon | Ref taxon | Ref accession | Identity | Align len |
|--------|------|-------------|-----------|---------------|----------|-----------|
| SRR9278435 | 16S | Xenia sp. | **same** | XM 047005986 | 100.0% | 621 |
| SRR9278435 | 18S | Xenia sp. | **same** | XR 006952126 | 100.0% | 1807 |
| SRR9278435 | 28S-v2 | Xenia sp. | **same** | XR 006952139 | 100.0% | 414 |
| SRR9278435 | 28S | Xenia sp. | **same** | XR 006952137 | 100.0% | 3234 |
| SRR9278435 | CO1 | Xenia sp. | **same** | HG999760 | 100.0% | 692 |
| SRR9278435 | ITS-v2 | Xenia sp. | **same** | KY442629 | 100.0% | 748 |
| SRR9278435 | ITS | Xenia sp. | **same** | KC864852 | 99.8% | 899 |
| SRR25099394 | 16S | Agalma elegans | **same** | OQ957203 | 100.0% | 542 |
| SRR25099394 | 18S | Agalma elegans | **same** | AY937313 | 100.0% | 1785 |
| SRR25099394 | 28S-v2 | Agalma elegans | Agalma elegans | EU272542 | 100.0% | 417 |
| SRR25099394 | 28S | Agalma elegans | **same** | EU272542 | 100.0% | 3230 |
| SRR25099394 | CO1 | Agalma elegans | **same** | OQ957203 | 100.0% | 692 |
| SRR25099394 | ITS-v2 | Agalma elegans | Agalma elegans | AY937313 | 100.0% | 671 |
| SRR25099394 | ITS | Agalma elegans | **same** | AY937313 | 100.0% | 817 |
| SRR8617500 | 16S | Rhopilema esculentum | **same** | NC 035741 | 100.0% | 563 |
| SRR8617500 | 18S | Rhopilema esculentum | **same** | XR 010517292 | 100.0% | 1796 |
| SRR8617500 | 28S-v2 | Rhopilema esculentum | **same** | XR 010517322 | 100.0% | 426 |
| SRR8617500 | 28S | Rhopilema esculentum | **same** | XR 010517322 | 100.0% | 3258 |
| SRR8617500 | CO1 | Rhopilema esculentum | **same** | NC 035741 | 100.0% | 692 |
| SRR8617500 | ITS-v2 | Rhopilema esculentum | **same** | KR338966 | 100.0% | 929 |

## Performance

| Sample | Max reads | Wall time | Peak memory | Reads ingested | Bases ingested | Distinct kmers |
|--------|----------:|----------:|------------:|---------------:|---------------:|---------------:|
| Xenia sp. | 1000k | 20.8s | 4.3 GB | 1,000,000 | 150,000,000 | 131,779,117 |
| Xenia sp. | 2000k | 35.6s | 8.5 GB | 2,000,000 | 300,000,000 | 263,603,227 |
| Xenia sp. | 4000k | 59.8s | 8.5 GB | 4,000,000 | 600,000,000 | 527,231,543 |
| Xenia sp. | 8000k | 106.1s | 17.0 GB | 8,000,000 | 1,200,000,000 | 1,054,598,174 |
| Agalma elegans | 1000k | 13.3s | 4.3 GB | 1,000,000 | 151,000,000 | 132,986,604 |
| Agalma elegans | 2000k | 26.4s | 8.5 GB | 2,000,000 | 302,000,000 | 265,974,355 |
| Agalma elegans | 4000k | 51.2s | 17.0 GB | 4,000,000 | 604,000,000 | 531,942,944 |
| Agalma elegans | 8000k | 93.5s | 17.0 GB | 8,000,000 | 1,208,000,000 | 1,063,887,262 |
| Rhopilema esculentum | 1000k | 14.1s | 4.3 GB | 1,000,000 | 150,000,000 | 131,979,520 |
| Rhopilema esculentum | 2000k | 28.7s | 8.5 GB | 2,000,000 | 300,000,000 | 263,969,663 |
| Rhopilema esculentum | 4000k | 52.1s | 8.5 GB | 4,000,000 | 600,000,000 | 527,945,714 |
| Rhopilema esculentum | 8000k | 94.5s | 17.0 GB | 8,000,000 | 1,200,000,000 | 1,055,890,972 |

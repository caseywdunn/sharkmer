# Panel validation: bacteria

- **Panel version**: `1.0.0`
- **sharkmer version**: `3.0.0`
- **Date**: 2026-04-13 12:27:59
- **Machine**: Linux 6.12.76-linuxkit, 12 cores, 47.0 GB RAM

## rhizosphere metagenome (SRR19418213)

| Gene | 1000k | Score |
|------|---:|:-----:|
| 16S-341F-785R | --- | `---` |
| 16S-PRK341F-PRK806R | 268bp | `+--` |
| 16S-515F-806R | 283bp | `+--` |
| 16S-515F-806RB | 283bp | `+--` |
| 16S-68F-783Rabc | 208bp | `+--` |
| 16S-341F-783Rabc | --- | `---` |

Wall times: 1000k: 2.2s

## coral metagenome (SRR24806237)

| Gene | 1000k | Score |
|------|---:|:-----:|
| 16S-341F-785R | 435bp | `+--` |
| 16S-PRK341F-PRK806R | 437bp | `+--` |
| 16S-515F-806R | 283bp | `+--` |
| 16S-515F-806RB | 283bp | `+--` |
| 16S-68F-783Rabc | --- | `---` |
| 16S-341F-783Rabc | 265bp | `+--` |

Wall times: 1000k: 0.6s

## seawater metagenome (ERR2596344)

| Gene | 1000k | Score |
|------|---:|:-----:|
| 16S-341F-785R | --- | `---` |
| 16S-PRK341F-PRK806R | --- | `---` |
| 16S-515F-806R | --- | `---` |
| 16S-515F-806RB | --- | `---` |
| 16S-68F-783Rabc | --- | `---` |
| 16S-341F-783Rabc | --- | `---` |

Wall times: 1000k: 0.4s

## Cross-sample summary (highest depth)

| Gene | rhizosphere metag... | coral metagenome | seawater metagenome |
|------|:---:|:---:|:---:|
| 16S-341F-785R | `---` | `+--` | `---` |
| 16S-PRK341F-PRK806R | `+--` | `+--` | `---` |
| 16S-515F-806R | `+--` | `+--` | `---` |
| 16S-515F-806RB | `+--` | `+--` | `---` |
| 16S-68F-783Rabc | `+--` | `---` | `---` |
| 16S-341F-783Rabc | `---` | `+--` | `---` |

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

### 16S-341F-785R

**Not recovered** in: SRR19418213, ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `CCTACGGGNGGCWGCAG` (trim=15, match window `TACGGGNGGCWGCAG`)

```
  spec          TACGGGNGGCWGCAG
                ||||||.|||.|||x
  SRR24806237   TACGGGAGGCAGCAA
```

  pos  7: N ({A,C,G,T}) partially utilised — observed {A}; could reduce to A
  pos 11: W ({A,T}) partially utilised — observed {A}; could reduce to A
  pos 15: G ({G}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `GACTACHVGGGTATCTAATCC` (trim=15, match window `HVGGGTATCTAATCC`)

```
  spec          HVGGGTATCTAATCC
                ..|||||||||||||
  SRR24806237   TGGGGTATCTAATCC
```

  pos  1: H ({A,C,T}) partially utilised — observed {T}; could reduce to T
  pos  2: V ({A,C,G}) partially utilised — observed {G}; could reduce to G

### 16S-PRK341F-PRK806R

**Not recovered** in: ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `CCTACGGGRBGCASCAG` (trim=15, match window `TACGGGRBGCASCAG`)

```
  spec          TACGGGRBGCASCAG
                ||||||..|x|.x|x
  SRR19418213   TACGGGGGGGACAAG
  SRR24806237   TACGGGAGGCAGCAA
```

  pos  7: R ({A,G}) fully utilised — observed {A,G}
  pos  8: B ({C,G,T}) partially utilised — observed {G}; could reduce to G
  pos 10: C ({C}) observed {C,G} — off-code base(s) absorbed by --mismatches; consider widening to S
  pos 12: S ({C,G}) fully utilised — observed {C,G}
  pos 13: C ({C}) observed {A,C} — off-code base(s) absorbed by --mismatches; consider widening to M
  pos 15: G ({G}) observed {A,G} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `GGACTACYVGGGTATCTAAT` (trim=15, match window `ACYVGGGTATCTAAT`)

```
  spec          ACYVGGGTATCTAAT
                ||..|||||||||||
  SRR19418213   ACCAGGGTATCTAAT
  SRR24806237   ACCGGGGTATCTAAT
```

  pos  3: Y ({C,T}) partially utilised — observed {C}; could reduce to C
  pos  4: V ({A,C,G}) partially utilised — observed {A,G}; could reduce to R

### 16S-515F-806R

**Not recovered** in: ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `GTGCCAGCMGCCGCGGTAA` (trim=15, match window `CAGCMGCCGCGGTAA`)

```
  spec          CAGCMGCCGCGGTAA
                ||||.||||||||||
  SRR19418213   CAGCAGCCGCGGTAA
  SRR24806237   CAGCAGCCGCGGTAA
```

  pos  5: M ({A,C}) partially utilised — observed {A}; could reduce to A

**Reverse primer** — spec `GGACTACHVGGGTWTCTAAT` (trim=15, match window `ACHVGGGTWTCTAAT`)

```
  spec          ACHVGGGTWTCTAAT
                ||..||||.||||||
  SRR19418213   ACTCGGGTTTCTAAT
  SRR24806237   ACCCGGGTTTCTAAT
```

  pos  3: H ({A,C,T}) partially utilised — observed {C,T}; could reduce to Y
  pos  4: V ({A,C,G}) partially utilised — observed {C}; could reduce to C
  pos  9: W ({A,T}) partially utilised — observed {T}; could reduce to T

### 16S-515F-806RB

**Not recovered** in: ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `GTGCCAGCMGCCGCGGTAA` (trim=15, match window `CAGCMGCCGCGGTAA`)

```
  spec          CAGCMGCCGCGGTAA
                ||||.||||||||||
  SRR19418213   CAGCAGCCGCGGTAA
  SRR24806237   CAGCAGCCGCGGTAA
```

  pos  5: M ({A,C}) partially utilised — observed {A}; could reduce to A

**Reverse primer** — spec `GGACTACNVGGGTWTCTAAT` (trim=15, match window `ACNVGGGTWTCTAAT`)

```
  spec          ACNVGGGTWTCTAAT
                ||..||||.||||||
  SRR19418213   ACTCGGGTTTCTAAT
  SRR24806237   ACCCGGGTTTCTAAT
```

  pos  3: N ({A,C,G,T}) partially utilised — observed {C,T}; could reduce to Y
  pos  4: V ({A,C,G}) partially utilised — observed {C}; could reduce to C
  pos  9: W ({A,T}) partially utilised — observed {T}; could reduce to T

### 16S-68F-783Rabc

**Not recovered** in: SRR24806237, ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `TNANACATGCAAGTCGRRCG` (trim=15, match window `CATGCAAGTCGRRCG`)

```
  spec          CATGCAAGTCGRRCG
                x|x||||||||..||
  SRR19418213   TAGGCAAGTCGGGCG
```

  pos  1: C ({C}) observed {T} — off-code base(s) absorbed by --mismatches; consider widening to Y
  pos  3: T ({T}) observed {G} — off-code base(s) absorbed by --mismatches; consider widening to K
  pos 12: R ({A,G}) partially utilised — observed {G}; could reduce to G
  pos 13: R ({A,G}) partially utilised — observed {G}; could reduce to G

**Reverse primer** — spec `CTACCAGGGTATCTAATCCTG` (trim=15, match window `GGGTATCTAATCCTG`)

```
  spec          GGGTATCTAATCCTG
                ||||x||||||||||
  SRR19418213   GGGTTTCTAATCCTG
```

  pos  5: A ({A}) observed {T} — off-code base(s) absorbed by --mismatches; consider widening to W

### 16S-341F-783Rabc

**Not recovered** in: SRR19418213, ERR2596344. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

**Forward primer** — spec `CCTACGGGNGGCWGCAG` (trim=15, match window `TACGGGNGGCWGCAG`)

```
  spec          TACGGGNGGCWGCAG
                ||||||.||x.|||x
  SRR24806237   TACGGGGGGATGCAA
```

  pos  7: N ({A,C,G,T}) partially utilised — observed {G}; could reduce to G
  pos 10: C ({C}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to M
  pos 11: W ({A,T}) partially utilised — observed {T}; could reduce to T
  pos 15: G ({G}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `CTACCAGGGTATCTAATCCTG` (trim=15, match window `GGGTATCTAATCCTG`)

```
  spec          GGGTATCTAATCCTG
                ||||x||||||||||
  SRR24806237   GGGTTTCTAATCCTG
```

  pos  5: A ({A}) observed {T} — off-code base(s) absorbed by --mismatches; consider widening to W

## Performance

| Sample | Max reads | Wall time | Peak memory | Reads ingested | Bases ingested | Distinct kmers |
|--------|----------:|----------:|------------:|---------------:|---------------:|---------------:|
| rhizosphere metagenome | 1000k | 2.2s | 272 MB | 262,991 | 79,160,291 | 74,411,239 |
| coral metagenome | 1000k | 0.6s | 67 MB | 83,843 | 20,964,921 | 19,455,591 |
| seawater metagenome | 1000k | 0.4s | 68 MB | 29,852 | 8,836,192 | 8,298,856 |

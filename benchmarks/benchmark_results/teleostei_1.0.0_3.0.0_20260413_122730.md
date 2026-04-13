# Panel validation: teleostei

- **Panel version**: `1.0.0`
- **sharkmer version**: `3.0.0`
- **Date**: 2026-04-13 12:56:23
- **Machine**: Linux 6.12.76-linuxkit, 12 cores, 47.0 GB RAM

## Nomeus gronovii (SRR22396603)

| Gene | 1000k | Score |
|------|---:|:-----:|
| 18S | 472bp (98.5%) | `+++` |
| 16S | 605bp (100.0%) | `+**` |
| CO1 | 685bp (100.0%) | `+**` |
| CytB | --- | `-+-` |
| 12S | 419bp (100.0%) | `+**` |

Wall times: 1000k: 13.5s

## Cross-sample summary (highest depth)

| Gene | Nomeus gronovii |
|------|:---:|
| 18S | `+++` |
| 16S | `+**` |
| CO1 | `+**` |
| CytB | `-+-` |
| 12S | `+**` |

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

### 18S

**Forward primer** — spec `TAACATATGCTTGTCTCAAAG` (trim=15, match window `ATGCTTGTCTCAAAG`)

```
  spec          ATGCTTGTCTCAAAG
                |||||||||||||||
  SRR22396603   ATGCTTGTCTCAAAG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `CCTGTATTGTTATTTTTCGTCAC` (trim=15, match window `GTTATTTTTCGTCAC`)

```
  spec          GTTATTTTTCGTCAC
                |||||||||||||||
  SRR22396603   GTTATTTTTCGTCAC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### 16S

**Forward primer** — spec `CGCCTGTTTATCAAAAACAT` (trim=15, match window `GTTTATCAAAAACAT`)

```
  spec          GTTTATCAAAAACAT
                |||||x|||||||||
  SRR22396603   GTTTACCAAAAACAT
```

  pos  6: T ({T}) observed {C} — off-code base(s) absorbed by --mismatches; consider widening to Y

**Reverse primer** — spec `CCGGTCTGAACTCAGATCACGT` (trim=15, match window `GAACTCAGATCACGT`)

```
  spec          GAACTCAGATCACGT
                |||||||||||||||
  SRR22396603   GAACTCAGATCACGT
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### CO1

**Forward primer** — spec `TCAACCAACCACAAAGACATTGGCAC` (trim=15, match window `CAAAGACATTGGCAC`)

```
  spec          CAAAGACATTGGCAC
                x||||||||x|||||
  SRR22396603   TAAAGACATCGGCAC
```

  pos  1: C ({C}) observed {T} — off-code base(s) absorbed by --mismatches; consider widening to Y
  pos 10: T ({T}) observed {C} — off-code base(s) absorbed by --mismatches; consider widening to Y

**Reverse primer** — spec `TAGACTTCTGGGTGGCCAAAGAATCA` (trim=15, match window `GTGGCCAAAGAATCA`)

```
  spec          GTGGCCAAAGAATCA
                x||||||||||||||
  SRR22396603   ATGGCCAAAGAATCA
```

  pos  1: G ({G}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to R

### CytB

**Not recovered** in: SRR22396603. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

### 12S

**Forward primer** — spec `ACTGGGATTAGATACCCCACTATG` (trim=15, match window `AGATACCCCACTATG`)

```
  spec          AGATACCCCACTATG
                |||||||||||||||
  SRR22396603   AGATACCCCACTATG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `GAGAGTGACGGGCGGTGT` (trim=15, match window `AGTGACGGGCGGTGT`)

```
  spec          AGTGACGGGCGGTGT
                |||||||||||||||
  SRR22396603   AGTGACGGGCGGTGT
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

## Reference match details

| Sample | Gene | Sample taxon | Ref taxon | Ref accession | Identity | Align len |
|--------|------|-------------|-----------|---------------|----------|-----------|
| SRR22396603 | 12S | Nomeus gronovii | **same** | NC 082547 | 100.0% | 419 |
| SRR22396603 | 16S | Nomeus gronovii | **same** | NC 082547 | 100.0% | 605 |
| SRR22396603 | 18S | Nomeus gronovii | Hirundichthys speculiger | OP088600 | 98.5% | 474 |
| SRR22396603 | CO1 | Nomeus gronovii | **same** | NC 082547 | 100.0% | 685 |

## Performance

| Sample | Max reads | Wall time | Peak memory | Reads ingested | Bases ingested | Distinct kmers |
|--------|----------:|----------:|------------:|---------------:|---------------:|---------------:|
| Nomeus gronovii | 1000k | 13.5s | 4.3 GB | 1,000,000 | 151,000,000 | 132,972,251 |

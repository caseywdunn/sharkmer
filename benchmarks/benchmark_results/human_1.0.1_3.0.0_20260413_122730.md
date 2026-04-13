# Panel validation: human

- **Panel version**: `1.0.1`
- **sharkmer version**: `3.0.0`
- **Date**: 2026-04-13 12:40:41
- **Machine**: Linux 6.12.76-linuxkit, 12 cores, 47.0 GB RAM

## SRR17535371

| Gene | 1000k | Score |
|------|---:|:-----:|
| mt-1404-3947 | --- | `-+-` |
| mt-3734-6739 | 2990bp (99.9%) | `+++` |
| mt-6511-9220 | 2694bp (99.9%) | `+++` |
| mt-8910-10648 | 1725bp (99.7%) | `+++` |
| mt-10360-12226 | 1853bp (99.6%) | `+++` |
| mt-11977-13830 | 1837bp (99.9%) | `+++` |
| mt-13477-15349 | --- | `-+-` |
| mt-14898-151 | 1808bp (99.3%) | `+++` |
| mt-16488-1677 | --- | `-+-` |

Wall times: 1000k: 0.1s

## Cross-sample summary (highest depth)

| Gene | SRR17535371 |
|------|:---:|
| mt-1404-3947 | `-+-` |
| mt-3734-6739 | `+++` |
| mt-6511-9220 | `+++` |
| mt-8910-10648 | `+++` |
| mt-10360-12226 | `+++` |
| mt-11977-13830 | `+++` |
| mt-13477-15349 | `-+-` |
| mt-14898-151 | `+++` |
| mt-16488-1677 | `-+-` |

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

### mt-1404-3947

**Not recovered** in: SRR17535371. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

### mt-3734-6739

**Forward primer** — spec `AAGTCACCCTAGCCATCATTCTA` (trim=15, match window `CTAGCCATCATTCTA`)

```
  spec          CTAGCCATCATTCTA
                |||||||||||||||
  SRR17535371   CTAGCCATCATTCTA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `GATATCATAGCTCAGACCATACC` (trim=15, match window `AGCTCAGACCATACC`)

```
  spec          AGCTCAGACCATACC
                |||||||||||||||
  SRR17535371   AGCTCAGACCATACC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-6511-9220

**Forward primer** — spec `CTGCTGGCATCACTATACTACTA` (trim=15, match window `ATCACTATACTACTA`)

```
  spec          ATCACTATACTACTA
                |||||||||||||||
  SRR17535371   ATCACTATACTACTA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `GATTGGTGGGTCATTATGTGTTG` (trim=15, match window `GGTCATTATGTGTTG`)

```
  spec          GGTCATTATGTGTTG
                |||||||||||||||
  SRR17535371   GGTCATTATGTGTTG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-8910-10648

**Forward primer** — spec `CTTACCACAAGGCACACCTACA` (trim=15, match window `CAAGGCACACCTACA`)

```
  spec          CAAGGCACACCTACA
                |||||||||||||||
  SRR17535371   CAAGGCACACCTACA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `GGCACAATATTGGCTAAGAGGG` (trim=15, match window `TATTGGCTAAGAGGG`)

```
  spec          TATTGGCTAAGAGGG
                |||||||||||||||
  SRR17535371   TATTGGCTAAGAGGG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-10360-12226

**Forward primer** — spec `GTCTGGCCTATGAGTGACTACA` (trim=15, match window `CTATGAGTGACTACA`)

```
  spec          CTATGAGTGACTACA
                ||||||x||||||||
  SRR17535371   CTATGAATGACTACA
```

  pos  7: G ({G}) observed {A} — off-code base(s) absorbed by --mismatches; consider widening to R

**Reverse primer** — spec `CAGTTCTTGTGAGCTTTCTCGG` (trim=15, match window `TGTGAGCTTTCTCGG`)

```
  spec          TGTGAGCTTTCTCGG
                |||||||||||||||
  SRR17535371   TGTGAGCTTTCTCGG
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-11977-13830

**Forward primer** — spec `CTCCCTCTACATATTTACCACAAC` (trim=15, match window `CATATTTACCACAAC`)

```
  spec          CATATTTACCACAAC
                |||||||||||||||
  SRR17535371   CATATTTACCACAAC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `AAGTCCTAGGAAAGTGACAGCGA` (trim=15, match window `GGAAAGTGACAGCGA`)

```
  spec          GGAAAGTGACAGCGA
                |||||||||||||||
  SRR17535371   GGAAAGTGACAGCGA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-13477-15349

**Not recovered** in: SRR17535371. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

### mt-14898-151

**Forward primer** — spec `TAGCCATGCACTACTCACCAGA` (trim=15, match window `GCACTACTCACCAGA`)

```
  spec          GCACTACTCACCAGA
                |||||||||||||||
  SRR17535371   GCACTACTCACCAGA
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

**Reverse primer** — spec `GGATGAGGCAGGAATCAAAGAC` (trim=15, match window `GCAGGAATCAAAGAC`)

```
  spec          GCAGGAATCAAAGAC
                |||||||||||||||
  SRR17535371   GCAGGAATCAAAGAC
```

  all positions fixed (no degeneracy codes) and all observed bases match — nothing to tune.

### mt-16488-1677

**Not recovered** in: SRR17535371. This may indicate insufficient primer degeneracy, insufficient read coverage, or the target taxon lacks the locus.

_(no samples recovered this gene; skipping alignment)_

## Reference match details

| Sample | Gene | Sample taxon | Ref taxon | Ref accession | Identity | Align len |
|--------|------|-------------|-----------|---------------|----------|-----------|
| SRR17535371 | mt-10360-12226 |  | Homo sapiens | NC 012920.1 | 99.6% | 1838 |
| SRR17535371 | mt-11977-13830 |  | Homo sapiens | NC 012920.1 | 99.9% | 1822 |
| SRR17535371 | mt-14898-151 |  | Homo sapiens | NC 012920.1 | 99.3% | 1794 |
| SRR17535371 | mt-3734-6739 |  | Homo sapiens | NC 012920.1 | 99.9% | 2975 |
| SRR17535371 | mt-6511-9220 |  | Homo sapiens | NC 012920.1 | 99.9% | 2679 |
| SRR17535371 | mt-8910-10648 |  | Homo sapiens | NC 012920.1 | 99.7% | 1710 |

## Performance

| Sample | Max reads | Wall time | Peak memory | Reads ingested | Bases ingested | Distinct kmers |
|--------|----------:|----------:|------------:|---------------:|---------------:|---------------:|
| SRR17535371 | 1000k | 0.1s | 5 MB | 30,024 | 6,395,186 | 5,852,957 |

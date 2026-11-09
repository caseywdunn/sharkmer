# Released 3.1.0 vs Candidate Analysis

This is calibration/regression analysis, not held-out biological validation.

## Primary 1M Timing

Median [min, max] wall seconds and RSS are across three paired runs. `n_kmers` means accepted k-mer occurrences.

| Cell | Actual records | Baseline wall s | Candidate wall s | Baseline RSS bytes | Candidate RSS bytes | Paired speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| angiospermae/ERR14009273/1000000 | 1000000 | 31.120 [30.980, 31.200] | 31.150 [31.140, 31.390] | 2287230976.000 [2287030272.000, 2287300608.000] | 2287312896.000 [2287222784.000, 2287316992.000] | 0.995 [0.994, 0.999] |
| angiospermae/SRR25378184/1000000 | 1000000 | 60.760 [60.720, 60.860] | 61.720 [61.480, 61.940] | 4569047040.000 [4568961024.000, 4569096192.000] | 4569096192.000 [4569010176.000, 4569169920.000] | 0.984 [0.980, 0.990] |
| bacteria/ERR2596344/1000000 | 59704 | 0.710 [0.710, 0.720] | 0.740 [0.740, 0.750] | 76865536.000 [76742656.000, 76967936.000] | 77090816.000 [77070336.000, 77090816.000] | 0.959 [0.959, 0.960] |
| bacteria/SRR19418213/1000000 | 525982 | 12.950 [12.900, 13.130] | 12.850 [12.840, 12.870] | 1146527744.000 [1146404864.000, 1146527744.000] | 1146613760.000 [1146585088.000, 1146802176.000] | 1.009 [1.004, 1.020] |
| bacteria/SRR24806237/1000000 | 167686 | 1.530 [1.520, 1.530] | 1.670 [1.660, 1.670] | 148267008.000 [148209664.000, 148279296.000] | 148353024.000 [148320256.000, 148414464.000] | 0.916 [0.910, 0.922] |
| cnidaria/SRR25099394/1000000 | 1000000 | 48.560 [48.540, 49.210] | 48.540 [48.360, 48.690] | 4568944640.000 [4568924160.000, 4568961024.000] | 4569063424.000 [4568825856.000, 4569223168.000] | 1.004 [0.997, 1.014] |
| cnidaria/SRR8617500/1000000 | 1000000 | 54.160 [53.970, 54.220] | 53.710 [53.490, 53.790] | 4568965120.000 [4568895488.000, 4569104384.000] | 4569227264.000 [4569137152.000, 4569391104.000] | 1.008 [1.005, 1.013] |
| cnidaria/SRR9278435/1000000 | 1000000 | 52.590 [52.520, 52.700] | 52.220 [51.880, 52.340] | 4569030656.000 [4569014272.000, 4569051136.000] | 4569276416.000 [4569104384.000, 4569358336.000] | 1.007 [1.003, 1.016] |
| human/SRR17535371/1000000 | 108518 | 0.390 [0.390, 0.390] | 0.430 [0.430, 0.440] | 23347200.000 [23293952.000, 23556096.000] | 23658496.000 [23506944.000, 23805952.000] | 0.907 [0.886, 0.907] |
| insecta/SRR1057608/1000000 | 1000000 | 44.440 [44.390, 44.820] | 44.310 [44.070, 44.780] | 2287226880.000 [2287128576.000, 2287292416.000] | 2287378432.000 [2287349760.000, 2287452160.000] | 1.002 [0.992, 1.017] |
| insecta/SRR27962769/1000000 | 1000000 | 98.700 [98.560, 99.070] | 98.510 [98.140, 98.520] | 4569018368.000 [4568944640.000, 4569108480.000] | 4569178112.000 [4569153536.000, 4569219072.000] | 1.004 [1.002, 1.006] |
| insecta/SRR31887760/1000000 | 1000000 | 60.110 [60.060, 60.140] | 64.100 [60.860, 64.140] | 2287316992.000 [2287190016.000, 2287349760.000] | 2287427584.000 [2287280128.000, 2287611904.000] | 0.937 [0.937, 0.988] |
| teleostei/SRR22396603/1000000 | 1000000 | 46.320 [46.290, 46.390] | 45.960 [45.870, 46.060] | 4569010176.000 [4568928256.000, 4569092096.000] | 4568920064.000 [4568903680.000, 4568997888.000] | 1.007 [1.006, 1.011] |

| Version | Sum of available per-cell median wall seconds | Cells with medians | Missing cells |
| --- | ---: | ---: | ---: |
| baseline | 512.340 | 13 | 0 |
| candidate | 515.910 | 13 | 0 |

## Deeper Descriptive Runs

Each 2M/4M/8M cell has one paired run; values are descriptive, not replicated timing estimates.

| Cell | Actual records | Baseline status/time/RSS | Candidate status/time/RSS | Retained | Lost | Gained |
| --- | ---: | --- | --- | ---: | ---: | ---: |
| cnidaria/SRR25099394/2000000 | 2000000 | complete/103.610/9132384256 | complete/102.840/9132335104 | 5 | 2 | 0 |
| cnidaria/SRR25099394/4000000 | 4000000 | complete/208.230/18259046400 | complete/219.170/18196353024 | 7 | 0 | 0 |
| cnidaria/SRR25099394/8000000 | 8000000 | complete/396.140/18259046400 | complete/395.110/18259312640 | 7 | 0 | 0 |
| cnidaria/SRR8617500/2000000 | 2000000 | complete/114.960/9132388352 | complete/114.100/9132539904 | 6 | 0 | 0 |
| cnidaria/SRR8617500/4000000 | 4000000 | complete/213.970/9132437504 | complete/214.670/9132380160 | 6 | 0 | 0 |
| cnidaria/SRR8617500/8000000 | 8000000 | complete/390.550/18259111936 | complete/392.960/18259304448 | 6 | 0 | 0 |
| cnidaria/SRR9278435/2000000 | 2000000 | complete/112.950/9132277760 | complete/111.990/9132490752 | 7 | 0 | 0 |
| cnidaria/SRR9278435/4000000 | 4000000 | complete/205.460/9132240896 | complete/205.970/9132363776 | 7 | 0 | 0 |
| cnidaria/SRR9278435/8000000 | 8000000 | complete/367.590/18259148800 | complete/371.050/18259304448 | 7 | 0 | 0 |
| insecta/SRR1057608/2000000 | 2000000 | complete/101.410/4568911872 | complete/100.830/4569034752 | 7 | 1 | 0 |
| insecta/SRR1057608/4000000 | 4000000 | complete/229.080/9132552192 | complete/236.420/9132703744 | 8 | 1 | 0 |
| insecta/SRR1057608/8000000 | 8000000 | complete/458.820/18259087360 | complete/469.260/18259124224 | 6 | 8 | 0 |
| insecta/SRR27962769/2000000 | 2000000 | complete/204.840/9132437504 | complete/205.580/9132429312 | 13 | 11 | 0 |
| insecta/SRR27962769/4000000 | 4000000 | complete/417.830/18259046400 | complete/420.830/18259640320 | 14 | 15 | 0 |
| insecta/SRR27962769/8000000 | 8000000 | complete/815.110/36512829440 | complete/819.730/36513075200 | 14 | 1 | 0 |
| insecta/SRR31887760/2000000 | 2000000 | complete/131.290/4568928256 | complete/141.250/4569120768 | 10 | 29 | 0 |
| insecta/SRR31887760/4000000 | 4000000 | complete/224.260/4568834048 | complete/225.920/4568981504 | 11 | 35 | 0 |
| insecta/SRR31887760/8000000 | 8000000 | complete/358.580/9132457984 | complete/371.250/9132527616 | 10 | 37 | 0 |

## Per-Cell Product and Status Changes

Product comparisons use the `(gene, length, SHA-256)` multiset, so biological gains/losses do not depend on output index. The JSON includes per-gene recovery/failure reasons and classification status counts.

| Cell | Baseline status | Candidate status | Retained | Lost | Gained | Product multiset identical | Index ordering identical |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| angiospermae/ERR14009273/1000000 | complete | complete | 6 | 0 | 0 | True | True |
| angiospermae/SRR25378184/1000000 | complete | complete | 6 | 0 | 0 | True | True |
| bacteria/ERR2596344/1000000 | complete | complete | 5 | 3 | 0 | False | False |
| bacteria/SRR19418213/1000000 | complete | complete | 2 | 47 | 0 | False | False |
| bacteria/SRR24806237/1000000 | complete | complete | 0 | 36 | 0 | False | False |
| cnidaria/SRR25099394/1000000 | complete | complete | 5 | 0 | 0 | True | True |
| cnidaria/SRR8617500/1000000 | complete | complete | 5 | 0 | 0 | True | True |
| cnidaria/SRR9278435/1000000 | complete | complete | 7 | 0 | 0 | True | True |
| human/SRR17535371/1000000 | complete | complete | 5 | 0 | 0 | True | True |
| insecta/SRR1057608/1000000 | complete | complete | 6 | 1 | 0 | False | False |
| insecta/SRR27962769/1000000 | complete | complete | 11 | 17 | 0 | False | False |
| insecta/SRR31887760/1000000 | complete | complete | 9 | 4 | 0 | False | False |
| teleostei/SRR22396603/1000000 | complete | complete | 4 | 0 | 0 | True | True |
| cnidaria/SRR25099394/2000000 | complete | complete | 5 | 2 | 0 | False | False |
| cnidaria/SRR25099394/4000000 | complete | complete | 7 | 0 | 0 | True | True |
| cnidaria/SRR25099394/8000000 | complete | complete | 7 | 0 | 0 | True | True |
| cnidaria/SRR8617500/2000000 | complete | complete | 6 | 0 | 0 | True | True |
| cnidaria/SRR8617500/4000000 | complete | complete | 6 | 0 | 0 | True | True |
| cnidaria/SRR8617500/8000000 | complete | complete | 6 | 0 | 0 | True | True |
| cnidaria/SRR9278435/2000000 | complete | complete | 7 | 0 | 0 | True | True |
| cnidaria/SRR9278435/4000000 | complete | complete | 7 | 0 | 0 | True | True |
| cnidaria/SRR9278435/8000000 | complete | complete | 7 | 0 | 0 | True | True |
| insecta/SRR1057608/2000000 | complete | complete | 7 | 1 | 0 | False | False |
| insecta/SRR1057608/4000000 | complete | complete | 8 | 1 | 0 | False | False |
| insecta/SRR1057608/8000000 | complete | complete | 6 | 8 | 0 | False | False |
| insecta/SRR27962769/2000000 | complete | complete | 13 | 11 | 0 | False | False |
| insecta/SRR27962769/4000000 | complete | complete | 14 | 15 | 0 | False | False |
| insecta/SRR27962769/8000000 | complete | complete | 14 | 1 | 0 | False | False |
| insecta/SRR31887760/2000000 | complete | complete | 10 | 29 | 0 | False | False |
| insecta/SRR31887760/4000000 | complete | complete | 11 | 35 | 0 | False | False |
| insecta/SRR31887760/8000000 | complete | complete | 10 | 37 | 0 | False | False |

## Primary Repeat Stability

Sequence and classification stability compare product multisets within each version; failed observations are not interpreted as zero-product runs.

| Cell | Version | Assessable | Sequence stable | Classification stable | Observed statuses |
| --- | --- | --- | --- | --- | --- |
| angiospermae/ERR14009273/1000000 | baseline | True | True | True |  |
| angiospermae/ERR14009273/1000000 | candidate | True | True | True |  |
| angiospermae/SRR25378184/1000000 | baseline | True | True | True |  |
| angiospermae/SRR25378184/1000000 | candidate | True | True | True |  |
| bacteria/ERR2596344/1000000 | baseline | True | True | True |  |
| bacteria/ERR2596344/1000000 | candidate | True | True | True |  |
| bacteria/SRR19418213/1000000 | baseline | True | True | True |  |
| bacteria/SRR19418213/1000000 | candidate | True | True | True |  |
| bacteria/SRR24806237/1000000 | baseline | True | True | True |  |
| bacteria/SRR24806237/1000000 | candidate | True | True | True |  |
| cnidaria/SRR25099394/1000000 | baseline | True | True | True |  |
| cnidaria/SRR25099394/1000000 | candidate | True | True | True |  |
| cnidaria/SRR8617500/1000000 | baseline | True | True | True |  |
| cnidaria/SRR8617500/1000000 | candidate | True | True | True |  |
| cnidaria/SRR9278435/1000000 | baseline | True | True | True |  |
| cnidaria/SRR9278435/1000000 | candidate | True | True | True |  |
| human/SRR17535371/1000000 | baseline | True | True | True |  |
| human/SRR17535371/1000000 | candidate | True | True | True |  |
| insecta/SRR1057608/1000000 | baseline | True | True | True |  |
| insecta/SRR1057608/1000000 | candidate | True | True | True |  |
| insecta/SRR27962769/1000000 | baseline | True | True | True |  |
| insecta/SRR27962769/1000000 | candidate | True | True | True |  |
| insecta/SRR31887760/1000000 | baseline | True | True | True |  |
| insecta/SRR31887760/1000000 | candidate | True | True | True |  |
| teleostei/SRR22396603/1000000 | baseline | True | True | True |  |
| teleostei/SRR22396603/1000000 | candidate | True | True | True |  |

## Observed Failures

Failures remain failures; no missing result is converted to a zero-product observation.

| Invocation | Status | Failure | Error |
| --- | --- | --- | --- |
| None | — | — | — |

## Representative Totals

Totals use primary pair 1 once per primary cell plus each deeper cell once; three primary repetitions are not tripled.

| Version | Complete cells | Products | Classification counts |
| --- | ---: | ---: | --- |
| baseline | 31 | 470 | {"ambiguous_gene": 17, "confirmed_gene_other_taxon": 10, "confirmed_product": 117, "insufficient_alignment": 38, "no_reference": 160, "no_significant_hit": 104, "split_or_chimeric_alignment": 12, "wrong_gene": 12} |
| candidate | 31 | 222 | {"ambiguous_gene": 17, "confirmed_gene_other_taxon": 9, "confirmed_product": 112, "insufficient_alignment": 31, "no_reference": 8, "no_significant_hit": 21, "split_or_chimeric_alignment": 12, "wrong_gene": 12} |

## Primary-Only Product Totals

One representative run per primary cell; no deeper cells or repeated runs included.

| Version | Complete cells | Products | Classification counts |
| --- | ---: | ---: | --- |
| baseline | 13 | 179 | {"ambiguous_gene": 5, "confirmed_gene_other_taxon": 2, "confirmed_product": 40, "insufficient_alignment": 7, "no_reference": 94, "no_significant_hit": 25, "split_or_chimeric_alignment": 3, "wrong_gene": 3} |
| candidate | 13 | 71 | {"ambiguous_gene": 5, "confirmed_gene_other_taxon": 2, "confirmed_product": 39, "insufficient_alignment": 6, "no_reference": 8, "no_significant_hit": 5, "split_or_chimeric_alignment": 3, "wrong_gene": 3} |

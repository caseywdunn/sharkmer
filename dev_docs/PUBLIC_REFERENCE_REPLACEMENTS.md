# Public replacements for quarantined references

Date: 2026-09-13. Follow-up to the [initial provenance audit](REFERENCE_PROVENANCE.md)
under [#129](https://github.com/caseywdunn/sharkmer/issues/129) and
[#153](https://github.com/caseywdunn/sharkmer/issues/153). **Not a release.**

## What the 61 entries mean

The 61 were **not all proven Sharkmer products**: 59 lacked an exact extraction
from their cited public source, and two matched a public sequence annotated
as mutS despite being labeled 16S. Origins of many mismatches remain unknown.
The Gryllus 978 bp ITS_2 reference is a demonstrated bootstrap product.
None of these old strings is reinstated as independent evidence.

**58 of the 61 now have public-reference alternatives:**

- 57 legacy entries are covered by 56 newly added reference identities; the
  two old C. elegans H3 entries share one independently supported genomic locus,
  not two asserted sample haplotypes.
- The documentation example was already replaced with independently sourced
  Hydra 16S in the initial audit. It stays Hydra, consistent with its Hydrozoa
  clade; the newly acquired Xenia 16S region is used in the cnidaria panel.
- Three remain unresolved: two plant trnV-atpE entries and Agalma ITS.
- Separately, the short public Gryllus CO1 reference expands from 304 to
  1,531 bp, and a 939 bp public Heliconius ND1 reference is added.

Among those 58 legacy-entry decisions, 41 have same-species/strain naming
(including one archive spelling difference), 16 use explicitly related-taxon
or distinct-isolate evidence, and one retains the Hydra documentation example.
These are entry counts, not 58 unique independent samples or complete haplotypes.

Active references increase **99 → 156**, backed by **121 versioned public
records**. Six panel patch versions advance. The complete original 159-entry
inventory and prior evidence archive remain unchanged and regression-only.

## Acquisition and interpretation

Sources were fetched independently from NCBI nucleotide records with raw
GenBank feature tables, requests, retrieval receipts, and versioned accessions.
Exact source coordinates, strand, topology, organism, and hashes are preserved.
No Sharkmer product was used to select bases, fill gaps, trim a replacement,
or infer source provenance. Reviewed annotations establish the stated gene or
locus scope, not infallible public-assembly accuracy.

Use the archive organism as the reference label. Related species, undescribed
isolates, predicted rRNAs, partial features, and ambiguity codes remain explicit.
The table's same-species labels describe source/sample naming, not a claim that
the archive individual shares the read set's haplotype. Undescribed Morbakka
and Xenia samples are not equated to named species or other isolates.

Important curation decisions:

- Rhopilema COX1 crosses the circular source origin on the reverse strand;
  the extraction follows the annotated join, not the whole mitochondrial genome.
- HG999760.1 contributes its contiguous partial COI gene feature, not its
  anomalous joined CDS annotation.
- Haliclystus 16S, 18S, and ITS now use same-species public sequences rather
  than silently relabeling congeners.
- C. elegans H3 uses the complete public YAC Z98866.1 with an independent
  his-72 annotation chain from the curated genomic locus. Its 433 bp region
  retains the 59 bp intron. EF1A uses complete cosmid FO080392.2; independent
  exon annotation and transcript identity show the 160 bp span does not cross
  a splice junction. Neither a cropped chromosome tile nor spliced mRNA is
  mislabeled as complete genomic truth.
- Five plant spacer references include independently mapped public primer
  flanks, not just intergenic bases. Each has one compatible full-primer pair
  in the searched source; original spacer coordinates remain in the annotation
  receipt. This does not guarantee primer binding in another individual.
- The legacy insecta NADH target is retained without renaming. Independently
  mapped configured primers span tRNA-Met into ND2 in public Gryllus and
  D. sechellia records. New ND2 references supply homologous gene evidence;
  partial D. tanythrix ND2 does not contain confirmed primer sites.
  An initial analysis accidentally selected CO1_2 by array index; the corrected
  receipt selects the primer by gene name and supersedes that analysis.

## Remaining gaps

**trnV-atpE:** in both inspected plant plastomes, the independently mapped
full-primer endpoints lie near trnV-UAC and in atpB, not atpE. Their spans are
1,396 and 1,376 bp, exceeding the configured 1,000 bp maximum. This is an
assay-label/length issue requiring separate primer review, not permission to
invent an atpE reference. Full-primer mapping alone does not establish
uncallability under Sharkmer's trimmed/mismatch-tolerant search. Primers and
bounds remain unchanged.
Follow-up: [#156](https://github.com/caseywdunn/sharkmer/issues/156).

**Agalma ITS:** AY937313.1 is annotated 18S, not ITS. The bounded targeted
public-record searches found no independently annotated Agalma ITS replacement.
This is missing evidence, not proof that the gene or a historical product is
absent or incorrect.

Some other references are partial or proxies. In particular, short ITS flanks
and partial Drosophila 12S/16S records can still limit coverage of a valid product.
An unsupported alignment is not automatically a wrong sequence.

## Evaluation safeguards

Reference **provenance** remains exact; product **comparison** remains
variation-tolerant. SNPs, indels, and coherent split alignments can support
gene identity without exact equality. Full sample haplotypes, repeat lengths,
variant phasing, and read support remain separate questions.

Reviewed indexed 16S/18S/28S/CO1/CO2 labels can share gene-level evidence while
retaining their original target IDs. In the reviewed insecta, cnidaria, and
c_elegans panels, ITS_1/ITS_2 are alternative primer pairs spanning the same
rDNA cluster, not names for biological ITS1 versus ITS2. Context-specific
equivalence must remain bound to the reviewed primer definitions; arbitrary
external ITS labels must not be collapsed. Gene support does not claim that
a product has the intended primer-bounded endpoints:
`primer_region_support: not_established`.

Counting, assembly, repeat gates, sample definitions, and default **k=19**
are unchanged. This is a source/annotation/evaluation update, not a new
throughput comparison or newly held-out validation cohort. Release approval
and the high-copy recovery investigation remain open.

## Reassessment of saved products

Reclassified the same 60 immutable prior invocations, after verifying their
result receipts and all 516 raw output-file receipts. The 20 three-replicate
groups have identical product fingerprints; 150 representative product
assessments are retained. No new Sharkmer runs or reads were used here.

Both versions now retain the **same 51 gene-supported product identities**
(38 expected-taxon and 13 other-taxon), versus the same 30 in the earlier
limited-reference evaluation. This is improved evaluation coverage from the
expanded references and logical-target handling, **not improved assembly
recovery**, and is not an ablation attributing the increase to either change
alone.

| Historical high-copy product | Expanded-reference assessment |
| --- | --- |
| Gryllus ITS_2, 978 bp | No significant hit. |
| Gryllus CO1_1, 373 bp | Same-taxon public CO1: 328/373 bp aligned at 100% identity; 87.94% query coverage, below the unchanged 90% gate. |
| Drosophila 12S, 475/487/499 bp | No significant hit. |
| Drosophila 16S_2, 590 bp | No significant hit. |
| Heliconius ND1, 263 bp | Same-taxon public ND1: 99.26% identity over 135 aligned query bases; 51.33% query coverage. |
| Retained Gryllus CO1_1, 352 bp (both versions) | Gene-supported: 100% identity and 99.72% query coverage against the expanded public CO1 reference. |

The seven losses still are **not demonstrated to be biologically wrong**.
Two have substantial gene-matching portions but insufficient full-query
coverage; five lack significant hits in this reference set. Missing evidence,
reference divergence/scope, population variation, or assembly artifacts are
not distinguished by these results alone. Do not relax thresholds just to
recover the old classifications, or claim that conservative assembly gates
are vindicated by absent reference support. Read-backed and independently
sample-linked evidence remain the next discriminating checks.

Validation: 97 Python tests (including real BLAST), 223 Rust unit tests,
22 Rust integration tests, all 10 panel schemas, source-extraction audit,
formatting, and Clippy pass. Sol implemented the scoped logical-target change;
Terra acquired plant/nematode sources and reassessed saved outputs; Astra
independently reviewed source annotations, integration, and focused BLAST
outcomes. Prior timing measurements are not rerun or relabeled here.

## Entry-level source map

Coordinates are zero-based, half-open on the forward public source; reverse
complement after extraction for `-`. A wrap traverses the circular origin.
Legacy IDs refer to the frozen pre-audit inventory, **not current array indices**.
Feature qualifiers, partial/predicted scope, source receipts, and per-entry
decisions are in the [replacement evidence archive](../benchmarks/benchmark_results/v3.2-public-replacements-20260913/README.md).

| Legacy ID / extra addition | Target | Public accession | Source organism | Coordinates; length | Source relation |
| --- | --- | --- | --- | --- | --- |
| `angiospermae:0:1` | psbA-trnH | [MN990625.1](https://www.ncbi.nlm.nih.gov/nuccore/MN990625.1) | Liriodendron tulipifera | 128–663 -; 535 bp | same species |
| `angiospermae:0:2` | psbA-trnH | [OQ613386.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ613386.1) | Opuntia engelmannii var. cuija | 88190–88684 +; 494 bp | related taxon |
| `angiospermae:2:0` | trnV-atpE | Unresolved | — | — | same species or unresolved |
| `angiospermae:2:1` | trnV-atpE | Unresolved | — | — | same species or unresolved |
| `angiospermae:3:1` | trnC-ycf6 | [NC_008326.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_008326.1) | Liriodendron tulipifera | 29375–30487 +; 1112 bp | same species |
| `angiospermae:6:1` | atpB-rbcL | [OR400606.1](https://www.ncbi.nlm.nih.gov/nuccore/OR400606.1) | Euphorbia jolkinii | 58344–59196 +; 852 bp | related taxon |
| `angiospermae:7:0` | trnL-F | [NC_049166.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_049166.1) | Acer sutchuenense subsp. tienchuanense | 49682–50111 +; 429 bp | related taxon |
| `angiospermae:8:1` | ITS | [LN680643.1](https://www.ncbi.nlm.nih.gov/nuccore/LN680643.1) | Euphorbia glareosa | 0–737 +; 737 bp | related taxon |
| `angiospermae:8:2` | ITS | [MZ366774.1](https://www.ncbi.nlm.nih.gov/nuccore/MZ366774.1) | Opuntia bonaerensis | 2406–2980 +; 574 bp | related taxon |
| `c_elegans:6:0` | EF1A | [FO080392.2](https://www.ncbi.nlm.nih.gov/nuccore/FO080392.2) | Caenorhabditis elegans | 8149–8309 -; 160 bp | same species |
| `c_elegans:7:0` | H3 | [Z98866.1](https://www.ncbi.nlm.nih.gov/nuccore/Z98866.1) | Caenorhabditis elegans | 320–753 -; 433 bp | same species |
| `c_elegans:7:1` | H3 | [Z98866.1](https://www.ncbi.nlm.nih.gov/nuccore/Z98866.1) | Caenorhabditis elegans | 320–753 -; 433 bp | same species; shared H3 |
| `cnidaria:0:1` | 16S | [KU257501.1](https://www.ncbi.nlm.nih.gov/nuccore/KU257501.1) | Haliclystus octoradiatus | 0–581 +; 581 bp | same species |
| `cnidaria:0:3` | 16S | [NC_035741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_035741.1) | Rhopilema esculentum | 14008–15687 -; 1679 bp | same species |
| `cnidaria:0:4` | 16S | [LC467070.1](https://www.ncbi.nlm.nih.gov/nuccore/LC467070.1) | Xenia sp. 2 TK-2019 | 0–698 +; 698 bp | related taxon |
| `cnidaria:1:3` | CO1 | [NC_035741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_035741.1) | Rhopilema esculentum | 15807–1542 - (wrap); 1590 bp | same species |
| `cnidaria:1:4` | CO1 | [HG999760.1](https://www.ncbi.nlm.nih.gov/nuccore/HG999760.1) | Xenia umbellata | 0–880 +; 880 bp | related taxon |
| `cnidaria:2:0` | 18S | [AY937313.1](https://www.ncbi.nlm.nih.gov/nuccore/AY937313.1) | Agalma elegans | 0–1755 +; 1755 bp | same species |
| `cnidaria:2:1` | 18S | [AY845346.1](https://www.ncbi.nlm.nih.gov/nuccore/AY845346.1) | Haliclystus octoradiatus | 0–1755 +; 1755 bp | same species |
| `cnidaria:2:2` | 18S | [GQ849083.1](https://www.ncbi.nlm.nih.gov/nuccore/GQ849083.1) | Morbakka virulenta | 0–1733 +; 1733 bp | related taxon |
| `cnidaria:2:4` | 18S | [XR_006952126.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_006952126.1) | Xenia sp. Carnegie-2017 | 0–1823 +; 1823 bp | same named strain |
| `cnidaria:3:0` | 28S_1 | [EU272542.1](https://www.ncbi.nlm.nih.gov/nuccore/EU272542.1) | Agalma elegans | 0–2403 +; 2403 bp | same species |
| `cnidaria:3:1` | 28S_1 | [KU308592.1](https://www.ncbi.nlm.nih.gov/nuccore/KU308592.1) | Haliclystus octoradiatus | 0–3166 +; 3166 bp | same species |
| `cnidaria:3:2` | 28S_1 | [GQ849060.1](https://www.ncbi.nlm.nih.gov/nuccore/GQ849060.1) | Morbakka virulenta | 0–2859 +; 2859 bp | related taxon |
| `cnidaria:3:4` | 28S_1 | [XR_006952137.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_006952137.1) | Xenia sp. Carnegie-2017 | 0–3589 +; 3589 bp | same named strain |
| `cnidaria:4:0` | 28S_2 | [EU272542.1](https://www.ncbi.nlm.nih.gov/nuccore/EU272542.1) | Agalma elegans | 0–2403 +; 2403 bp | same species |
| `cnidaria:4:1` | 28S_2 | [KU308592.1](https://www.ncbi.nlm.nih.gov/nuccore/KU308592.1) | Haliclystus octoradiatus | 0–3166 +; 3166 bp | same species |
| `cnidaria:4:2` | 28S_2 | [GQ849060.1](https://www.ncbi.nlm.nih.gov/nuccore/GQ849060.1) | Morbakka virulenta | 0–2859 +; 2859 bp | related taxon |
| `cnidaria:5:0` | ITS_1 | Unresolved | — | — | unresolved |
| `cnidaria:5:1` | ITS_1 | [KU308625.1](https://www.ncbi.nlm.nih.gov/nuccore/KU308625.1) | Haliclystus octoradiatus | 0–653 +; 653 bp | same species |
| `cnidaria:5:2` | ITS_1 | [KC864852.1](https://www.ncbi.nlm.nih.gov/nuccore/KC864852.1) | Xenia sp. 1202010 | 0–951 +; 951 bp | related taxon |
| `cnidaria:6:0` | ITS_2 | [KU308625.1](https://www.ncbi.nlm.nih.gov/nuccore/KU308625.1) | Haliclystus octoradiatus | 0–653 +; 653 bp | same species |
| `cnidaria:6:1` | ITS_2 | [KR338966.1](https://www.ncbi.nlm.nih.gov/nuccore/KR338966.1) | Rhopilema esculentum | 0–813 +; 813 bp | same species |
| `cnidaria:6:2` | ITS_2 | [KY442629.1](https://www.ncbi.nlm.nih.gov/nuccore/KY442629.1) | Xenia sp. | 0–935 +; 935 bp | related taxon |
| `extra:Heliconius_ND1` | ND1 | [NC_024741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024741.1) | Heliconius pachinus | 11661–12600 -; 939 bp | same species |
| `hydrozoa:0:1` | 16S | [EU293971.1](https://www.ncbi.nlm.nih.gov/nuccore/EU293971.1) | Craspedacusta sowerbii | 0–539 +; 539 bp | same species archive spelling |
| `insecta:0:0` | 12S | [EU494495.1](https://www.ncbi.nlm.nih.gov/nuccore/EU494495.1) | Drosophila adunca | 256–697 -; 441 bp | related taxon |
| `insecta:10:0` | Yp2 | [L14423.1](https://www.ncbi.nlm.nih.gov/nuccore/L14423.1) | Drosophila melanogaster | 0–1114 +; 1114 bp | same species |
| `insecta:12:0` | 28S | [NR_133562.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_133562.1) | Drosophila melanogaster | 0–3970 +; 3970 bp | same species |
| `insecta:12:1` | 28S | [KM508879.1](https://www.ncbi.nlm.nih.gov/nuccore/KM508879.1) | Gryllus bimaculatus | 0–1266 +; 1266 bp | same species |
| `insecta:13:2` | 18S_2 | [XR_006246778.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_006246778.1) | Drosophila grimshawi | 0–1994 +; 1994 bp | related taxon |
| `insecta:14:0` | ITS_1 | [KR270069.1](https://www.ncbi.nlm.nih.gov/nuccore/KR270069.1) | Drosophila diamphidiopoda | 0–558 +; 558 bp | related taxon |
| `insecta:15:0` | CytB | [PP230540.1](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1) | Gryllus bimaculatus | 10361–11498 +; 1137 bp | same species |
| `insecta:15:1` | CytB | [NC_024741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024741.1) | Heliconius pachinus | 10432–11581 +; 1149 bp | same species |
| `insecta:16:0` | ND5 | [NC_024741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024741.1) | Heliconius pachinus | 6332–8066 -; 1734 bp | same species |
| `insecta:17:0` | ITS_2 | [MK441842.1](https://www.ncbi.nlm.nih.gov/nuccore/MK441842.1) | Gryllus bimaculatus | 0–664 +; 664 bp | same species |
| `insecta:1:3` | 16S_1 | [NC_024741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024741.1) | Heliconius pachinus | 12668–14048 -; 1380 bp | same species |
| `insecta:2:0` | CO1_1 | [OM754561.1](https://www.ncbi.nlm.nih.gov/nuccore/OM754561.1) | Drosophila tanythrix | 0–650 +; 650 bp | same species |
| `insecta:2:1` | CO1_1 | [PP230540.1](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1) | Gryllus bimaculatus | 1408–2939 +; 1531 bp | same species |
| `insecta:3:1` | CO2_1 | [HQ170726.1](https://www.ncbi.nlm.nih.gov/nuccore/HQ170726.1) | Drosophila tanythrix | 0–688 +; 688 bp | same species |
| `insecta:5:0` | ND4 | [MK659833.1](https://www.ncbi.nlm.nih.gov/nuccore/MK659833.1) | Drosophila neonasuta | 8267–9606 -; 1339 bp | related taxon |
| `insecta:5:1` | ND4 | [PP230540.1](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1) | Gryllus bimaculatus | 8066–9410 -; 1344 bp | same species |
| `insecta:6:1` | 16S_2 | [HQ171031.1](https://www.ncbi.nlm.nih.gov/nuccore/HQ171031.1) | Drosophila tanythrix | 0–383 +; 383 bp | same species |
| `insecta:6:2` | 16S_2 | [PP230540.1](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1) | Gryllus bimaculatus | 12590–13903 -; 1313 bp | same species |
| `insecta:6:3` | 16S_2 | [NC_024741.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_024741.1) | Heliconius pachinus | 12668–14048 -; 1380 bp | same species |
| `insecta:7:1` | CO1_2 | [NC_005780.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_005780.1) | Drosophila sechellia | 1475–3012 +; 1537 bp | same species |
| `insecta:7:2` | CO1_2 | [OM754561.1](https://www.ncbi.nlm.nih.gov/nuccore/OM754561.1) | Drosophila tanythrix | 0–650 +; 650 bp | same species |
| `insecta:7:3` | CO1_2 | [KP074829.1](https://www.ncbi.nlm.nih.gov/nuccore/KP074829.1) | Heliconius pachinus | 0–1521 +; 1521 bp | same species |
| `insecta:8:2` | CO2_2 | [HQ170726.1](https://www.ncbi.nlm.nih.gov/nuccore/HQ170726.1) | Drosophila tanythrix | 0–688 +; 688 bp | same species |
| `insecta:9:1` | NADH | [KM252368.1](https://www.ncbi.nlm.nih.gov/nuccore/KM252368.1) | Drosophila tanythrix | 31–520 +; 489 bp | same species |
| `insecta:9:2` | NADH | [PP230540.1](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1) | Gryllus bimaculatus | 204–1221 +; 1017 bp | same species |
| `reference:0:0` | 16S | [NC_011220.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_011220.1) | Hydra vulgaris | 1286–1804 +; 518 bp | different taxon documentation example |
| `teleostei:0:0` | 18S | [XR_010033285.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_010033285.1) | Engraulis encrasicolus | 0–1857 +; 1857 bp | same species |

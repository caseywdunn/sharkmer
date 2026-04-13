---
name: primer-panel
description: >
  Guides creation and optimization of sharkmer primer panel YAML files for
  in silico PCR (sPCR). Use this skill whenever the user wants to build a new
  primer panel for a taxonomic group, add primers to an existing panel, find
  SRA validation samples, search for reference sequences, run the panel
  validator, or troubleshoot why amplicons aren't matching references. Trigger
  on any mention of: primer panel, panel YAML, sharkmer panel, sPCR primers,
  validate_panel.py, SRA samples for validation, barcoding primers, or
  "my primers aren't working". The argument is the path to the panel file to
  create or edit (e.g. `panels/mollusca.yaml`).
---

# Primer Panel Builder

You are guiding the user through the creation or optimization of a sharkmer
primer panel YAML. Work interactively — confirm with the user at key decision
points rather than charging ahead.

The panel file format and primer design criteria are documented in PCR.md at
the repo root. Read it before starting if you haven't yet.

---

## Step 1 — Panel file scaffolding

Open the panel file the user specified (create it if it doesn't exist).

Check for these required fields and fill any that are missing, asking the user
if you need information you don't have:

**`name`** — short identifier, usually the clade name in lowercase
(e.g. `mollusca`). Ask the user if absent.

**`description`** — one sentence describing what the panel targets. Ask the
user if absent.

**`schema_version`** — set to `"2"`.

**`panel_version`** — set to `1.0.0` if absent. This is the panel's own
version, independent of sharkmer's version.

**`clade`** — NCBI-preferred name for the target taxon (e.g. `"Cnidaria"`).
Ask the user if absent.

**`taxon_id`** — NCBI Taxonomy ID for the clade. Look it up at
https://www.ncbi.nlm.nih.gov/taxonomy if needed.

**`maintainers`** — at least one entry with `name`. Ask for the user's name
and ORCID (ORCID is optional).

**`changelog`** — initial entry:
```yaml
changelog:
  - panel_version: 1.0.0
    date: <today>
    sharkmer_version: "3.0.0"
    changes: "Initial panel."
```

---

## Step 2 — Primers

### 2a. If the panel already has primers — assess them

For each primer, calculate its effective specificity and flag problems:

**Specificity per position:**
- Fully resolved (A/C/G/T): 2.0 bits
- 2-fold degenerate (R/Y/S/W/K/M): 1.0 bit
- 3-fold degenerate (B/D/H/V): 0.42 bits
- N: 0 bits

Sum over all positions (both forward and reverse). Then subtract 2 bits per
mismatch allowed (`mismatches` field, default 2) × 2 primers = 4 × 2 = 8 bits.

Flag and explain to the user:
- **Total effective specificity < 40 bits**: low — likely to produce many
  off-target amplicons, especially in large or complex genomes.
- **Either primer < 15 bp**: shorter than the default `trim=15`; set `trim`
  to `len(primer) - 1` at most, or better, note that this primer may underperform.
- **High degeneracy relative to length**: if one primer has more than 3
  degenerate bases in a 15 bp stretch, suggest finding a less degenerate
  15 bp sub-region (see §2b for how).
- **Specific or conserved bases located in the 5' overhang**: sharkmer seeds
  the de Bruijn graph from the **3' `trim` bases** of each primer only — by
  default the 3'-most 15 bp. Bases upstream of that window (the 5' tail) do
  not contribute to graph seeding and provide no specificity benefit during
  assembly. If the most diagnostic or conserved positions in a primer fall in
  the 5' tail, those positions are wasted. To exploit them, either increase
  `trim` (up to k−1 = 18 at default k=19) so the window reaches them, or
  redesign the primer so the specific positions are within the 3' trim window.

### 2b. If the panel has no primers — mine the literature

Search for common phylogenetics and barcoding primers used for the clade.
Good sources: published papers, the Primer3, Barcode of Life Data System
(BOLD), and GenBank primer notes. Use web search.

For each candidate primer pair:
1. **Record the published sequence verbatim** as a YAML comment (`#`) alongside
   the trimmed entry — this is provenance for future edits.
2. **Assess length and degeneracy.** If a primer is long (> 20 bp) and
   degenerate, find the least-degenerate contiguous 15 bp sub-region to use
   as the working sequence, maximising specificity bits. Show your reasoning.
   If the best 15 bp window still has > 3 degenerate bases, note this and ask
   the user whether to proceed or look for alternatives.
3. **Set trim** to `min(len(trimmed_seq), 18)` — the maximum useful trim at
   default k=19 is k−1 = 18. Never set trim > 18 without changing k.
4. **Set initial length window** conservatively wide (e.g. ±30% of the
   published expected amplicon size) so the first validation run can show you
   what sizes actually come out.

### 2c. Verify citations before writing the panel

For every primer pair — whether sourced from the literature or already in the
panel — verify the citation before writing it into the file. Hallucinated DOIs
are common and silently wrong citations undermine the panel's provenance.

**For each citation with a DOI:**

1. Fetch the DOI URL to confirm it resolves:
   ```
   WebFetch https://doi.org/<doi>
   ```
   If the DOI returns a 404 or does not resolve, the DOI is wrong. Do not
   write it into the panel. Instead tell the user what you found and ask them
   to supply the correct reference.

2. If the paper is accessible (open access or abstract available), read it and
   confirm the primer sequence appears in the text, supplementary table, or
   methods section. A citation is only trustworthy if the primer sequence can
   actually be found in the cited source.

3. If the paper is paywalled and the full text is not available, note this
   explicitly alongside the citation in the panel's `notes:` field, e.g.:
   `notes: "... Citation verified DOI resolves; full text paywalled, primer sequence not confirmed in paper."`

4. If you cannot verify a primer in its cited source, do not silently accept
   it. Flag it to the user: "I could not confirm this primer sequence appears
   in [citation] — please verify manually before using."

**For citations without a DOI** (books, conference proceedings, etc.), do a
web search for the title + author to confirm the work exists and note any
concerns.

Offer a table of candidate primer pairs to the user before writing the panel,
showing: gene name, forward seq, reverse seq, total specificity bits, citation
verified (yes / paywalled / failed), any flags.
Get confirmation before proceeding.

---

## Step 3 — Validation samples

### 3a. Map the major phylogenetic groups before searching

Before searching SRA, establish what sub-clades need to be represented.
**Do not rely solely on NCBI taxonomy** — it reflects administrative
classification and contains non-phylogenetic nodes that do not need sampling:

- Nodes named "environmental samples", "unclassified X", "X sp.",
  "metagenomes", or "other X" are not phylogenetic groups. Skip them.
- NCBI sometimes splits or lumps clades in ways that do not match current
  phylogenetic consensus.

Instead, do a brief web search for the phylogenetic literature on the focal
clade (e.g. "Anthozoa phylogeny", "Mollusca higher-level phylogeny") to
identify the major crown-group lineages that span the root of the clade.
These are the groups that must each have at least one sample for the panel
to be considered broadly validated.

**Example — Anthozoa:** Phylogenetics recognises two major clades:
Hexacorallia (stony corals, sea anemones, zoanthids) and Octocorallia
(soft corals, sea fans, sea pens). NCBI Taxonomy also lists "unclassified
Anthozoa" and "environmental samples Anthozoa" — those are not phylogenetic
groups and do not need representation. A valid Anthozoa validation set
needs at least one Hexacorallia and one Octocorallia sample.

Build a target group list before searching, e.g.:

```
Target groups for Anthozoa panel:
  [x] Hexacorallia
  [ ] Octocorallia   ← must find a sample
```

### 3b. Search for samples

Search NCBI SRA for 4–10 samples, aiming for one or more per target group.

**Required criteria for each sample:**
- Platform: Illumina (preferably NovaSeq 6000 or NovaSeq X)
- Layout: PAIRED or SINGLE, read length 150 bp
- Spots: ≥ 20 million (search filter: `spots >= 20000000`)
- Library attributes: Strategy = WGS, Source = GENOMIC, Selection = RANDOM
  (these filters exclude RNA-seq, amplicon, and targeted sequencing)

Search SRA with queries like:
```
("Mollusca"[Organism] OR "Gastropoda"[Organism]) AND "WGS"[Strategy] AND
"GENOMIC"[Source] AND "RANDOM"[Selection] AND "Illumina"[Platform] AND
"PAIRED"[Layout]
```

For each candidate accession, check: the taxonomy (does it represent a
distinct part of the clade?), the number of spots, and whether the library
metadata matches WGS/GENOMIC/RANDOM.

### 3c. Check breadth and alert on gaps

After assembling the candidate list, map each sample back to the target
groups from §3a. If any major phylogenetic group has no sample:

1. **Search specifically for that group.** Use the group name directly as
   the Organism filter. Accept a genus or family-level accession if a
   species-level one is unavailable.

2. **If no WGS data exists for a group**, tell the user explicitly:

   > "No WGS data is available in SRA for `<group>`. The panel cannot be
   > validated for this lineage and its primers may not perform there.
   > Consider either: (a) narrowing the stated scope of the panel to exclude
   > `<group>` until data become available, or (b) keeping the broad scope
   > but documenting in the panel `notes:` that `<group>` is unvalidated."

   Ask the user which they prefer before proceeding.

3. **Annotate the proposed sample table** with which target group each
   sample represents, so the breadth is visible at a glance.

Set `max_reads: [1000000, 2000000, 4000000]` as the default depth sweep for
each sample. If a sample has a very large genome (> 2 Gb), note that 4 million
reads may be insufficient coverage — suggest increasing to 8 million.

Present the proposed sample list to the user before writing it into the panel,
including taxon, accession, spots, target group represented, and why it was
chosen. Get confirmation.

---

## Step 4 — Reference sequences

For each primer pair and each validation sample taxon, search NCBI for a
reference sequence to use for BLAST validation. A reference is a known-good
sequence for that gene from that species (or a close relative).

Good sources: NCBI Nucleotide, RefSeq (mitochondrial genomes for mt genes,
rRNA databases for 18S/28S/ITS), BOLD for CO1.

Search strategy per gene:
```
"<taxon>" [Organism] AND "<gene>" [Gene Name]
```
Or for mitochondrial genes:
```
"<taxon>" [Organism] AND "mitochondrion" [Title] AND complete genome
```

For each found sequence, extract the relevant amplicon region (the stretch
between the forward and reverse primer binding sites). Write the references
directly into the panel YAML file under a top-level `references:` block —
**not in `notes:` fields**. The `gene:` key in each reference entry must
match the derived output name for that primer:
- `gene` only → `gene` (e.g., `"CO1"`)
- `gene` + `index` → `gene_index` (e.g., `"CO1_1"`)
- `gene` + `region` → `gene-region` (e.g., `"16S-V3"`)
- `gene` + `region` + `index` → `gene-region_index` (e.g., `"18S-V9_2"`)

The format is:

```yaml
references:
  - gene: "CO1"
    sequences:
      - taxon: Lymnaea stagnalis
        accession: AY382548
        sequence: GGTCAACAAATCATAAAGATATTGG...TAAACTTCAGGGTGACCAAAAAATCA
      - taxon: Biomphalaria glabrata
        accession: AF317857
        sequence: GGTCAACAAATCATAAAGATATTGG...TAAACTTCAGGGTGACCAAAAAATCA
  - gene: "16S"
    sequences:
      - taxon: Lymnaea stagnalis
        accession: AJ390977
        sequence: CGCCTGTTTAYCAAAAACAT...CCGGTCTGAACTCAGATCACGT
```

It is OK if not all taxa have a
reference for all genes — aim for at least one reference per gene, and at
least partial coverage across taxa.

Be selective: only include sequences where you are confident the gene
annotation is correct and the sequence is from the intended gene. GenBank
records with vague annotations or misidentified taxonomy should be skipped.

Show the user the reference list before writing it and flag any gene/taxon
combinations where you could not find a reference.

---

## Step 5 — Run the validator

```bash
conda activate sharkmer-bench
python scripts/validate_panel.py panels/<panel>.yaml
```

If `conda activate` fails with *"Run 'conda init' before 'conda activate'"*:
```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate sharkmer-bench
```

Read the markdown report written to `panels/validation_reports/`. For each
gene × sample combination, you will see one of:

| Outcome | What it means |
|---------|---------------|
| `+ / * / *` | Product recovered, BLAST matches same gene same species — ideal |
| `+ / + / +` | Product recovered, BLAST matches same gene different species — good |
| `+ / - / -` | Product recovered, no reference — confirm on-target manually |
| `+ / * / -` | Product recovered but BLAST *misses* own-species reference — suspicious |
| `- / * / -` | No product, reference exists — primer failure |
| `- / - / -` | No product, no reference — ambiguous |

---

## Step 6 — Diagnose and iterate

Work through failures systematically. For each gene with problems:

### Products present but no or wrong BLAST match

Run remote BLAST on the product sequence:
```bash
blastn -query <product.fasta> -db nt -remote -outfmt "6 qseqid sseqid stitle pident length evalue" | head -10
```

If the top hits are all from the expected gene/clade → the product is on-target;
the reference may be from the wrong region or have annotation issues. Pull a
better reference.

If top hits are from unrelated genes or taxa → off-target amplification. The
primer pair needs more specificity (see §6b).

### No products from some samples but products from others

The primers may not match binding sites in the non-recovering taxa. To diagnose:

1. Get the reference sequences for the non-recovering taxa (from §4 or from
   NCBI BLAST of a related species).
2. Manually examine the binding site region in those references against the
   primer. Note every position where the primer and reference differ.
3. If differences are at positions where the primer is already degenerate →
   the `mismatches` parameter is absorbing them; the real issue may be coverage
   (try increasing `max_reads`) or a graph size problem (see genome size below).
4. If differences are at resolved positions in the primer → consider adding
   degeneracy at those positions (R/Y/M/K etc.). Calculate the new specificity
   to make sure you are not dropping below ~40 bits total. Keep the original
   primer sequence as a comment.
5. If the primer sequence consistently does not match the clade at all →
   discuss with the user whether to redesign the primer from alignment of
   multiple reference sequences, or narrow the panel scope to a sub-clade.

### No products from any sample for a gene

Two possible causes:

**A. Primer design problem.** Align available reference sequences for the
target gene across the clade. Look at the binding site region for both primers.
If the published primer has many mismatches to all references, redesign. Use
the most conserved 15 bp stretch in each binding site as the new working
sequence.

**B. Coverage / genome size problem.** Check the approximate genome size for
the clade. The reads needed depend critically on whether the target is
high-copy or single-copy:

- **Mitochondrial / high-copy genes** (16S, CO1, CytB, etc.): mitochondria
  are present at hundreds to thousands of copies per cell, so effective
  coverage is far higher than the genome-wide average. 1–2 M reads is
  usually sufficient even for genomes up to several Gb.

- **Single-copy nuclear genes** (H3, ITS, 28S, etc.): need at minimum ~3×
  average genomic coverage for reliable assembly. With 150 bp reads:
  `reads ≈ 3 × genome_size_bp / 150`. Rough targets:
  - 1 Gb genome → ~20 M reads
  - 2 Gb genome → ~40 M reads
  - 5 Gb genome → ~100 M reads

  The default 1–4 M sweep is almost never sufficient for single-copy nuclear
  genes in animals. Tell the user and suggest either higher `max_reads` tiers
  (8–32 M) or switching to a WGS dataset with more spots.

### Poor-performing samples

After reviewing the validation report, look across all genes for each sample
and compute a rough per-sample gene recovery rate (number of genes with a
product ÷ total genes in the panel). If a sample recovers considerably fewer
genes than the others — a useful threshold is recovering fewer than half the
genes that pass for the median sample, or failing on every gene that succeeds
for all other samples — treat it as a potential technical outlier before
concluding it is a biological signal.

**Step 1 — Find an alternative sample.**
Search SRA for a second accession from the same species, or from the most
closely related taxon with available WGS data, using the same criteria as
Step 3 (Illumina, PAIRED, ≥ 20 M spots, WGS/GENOMIC/RANDOM). Add it to the
validation block temporarily alongside the original.

**Step 2 — Re-run the validator.**
```bash
python scripts/validate_panel.py panels/<panel>.yaml
```

**Step 3 — Compare recovery between the two accessions.**

- **Alternative performs markedly better** (recovers most genes that the
  original missed): the original dataset likely has a technical problem —
  low base quality, heavy contamination, high duplicate rate, or an unusual
  library prep. Tell the user:

  > "Sample `<accession>` (`<taxon>`) recovered `N` of `M` genes, while
  > `<alternative>` from the same taxon recovered `N'` of `M`. The original
  > sample may have a technical issue. Should I remove it and keep the
  > alternative, or keep both? Keeping the original is reasonable if it
  > represents a biologically distinct population or a divergent specimen —
  > in that case add a `notes:` field explaining the known limitation."

- **Alternative also performs poorly**: the low recovery is probably
  biological, not a data quality problem. The primers may not match this
  lineage well. Continue with §6b (No products from some samples) for the
  affected genes. Remove the temporary alternative sample if it does not add
  taxonomic breadth.

- **Both perform similarly and reasonably well**: the original sample is
  fine; the low recovery for specific genes is a primer issue, not a data
  issue. Remove the temporary alternative and proceed with §6b.

### Primer changes

When modifying a primer sequence:
- Keep the original sequence as a YAML comment: `# original: ORIGINAL_SEQ`
- Bump the `panel_version` patch digit (e.g. 1.0.0 → 1.0.1)
- Add a `changelog` entry describing what changed and why
- Re-run the validator to confirm improvement

### Splitting the panel

If no single primer pair can robustly cover the full clade due to too much
sequence variation at the binding sites, tell the user explicitly. Suggest
either:
- Creating sub-clade panels (e.g. `cnidaria_anthozoa.yaml`,
  `cnidaria_medusozoa.yaml`)
- Keeping the broad primers but documenting expected failure for certain
  sub-groups in the `notes` field

---

## Step 7 — Finalize

Once the panel validates well (most targets recovering with correct BLAST
matches at the depths you care about):

1. Update `validation.last_validated` using `--write`:
   ```bash
   python scripts/validate_panel.py panels/<panel>.yaml --write
   ```
2. Confirm `panel_version`, `changelog`, and `maintainers` are up to date.
3. Summarize what the panel covers and any known limitations to the user.

---

## Quick reference: primer quality checklist

Before finalising any primer pair, confirm:
- [ ] Both primers ≥ 15 bp
- [ ] `trim` ≤ min(len(primer), k−1) — default k=19 means trim ≤ 18
- [ ] Total effective specificity ≥ 40 bits (after mismatch deductions)
- [ ] Original published sequence recorded as a comment if the sequence was modified
- [ ] `min_length`/`max_length` window is wide enough to capture biological variation (≥ ±20% of expected)
- [ ] `citation` field populated and doi valid and correct

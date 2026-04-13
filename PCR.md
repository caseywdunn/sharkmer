# PCR primers and panels

This document covers primer design considerations for *in silico* PCR (sPCR)
and the panel file format used to bundle primer pairs for sharkmer. For
general tool usage, see [README.md](README.md).

## Primer specificity

The tradeoffs in designing primers for *in silico* PCR are the same as for
benchtop PCR. More specific primers bind fewer off-target sites, producing
faster and more reliable results on the taxa they were designed for — at the
extreme, species-specific primers. But more specific primers work across a
narrower phylogenetic range. A panel targeting a broad clade like Cnidaria
needs enough non-specificity to accommodate sequence variation across the
clade, while remaining specific enough that off-target hits do not overwhelm
real ones. Ideally, primers are as specific as possible for the clade they
will be used in, with just enough degeneracy to span its variation.

### Measuring specificity in bits

Primer optimization is largely about tuning specificity. Increase specificity for a given taxon, and it works faster, more reliably, and with lower read coverage, but it will work less well for other taxa. Decrease specificity, and you can increase the phylogenetic breadth at which the primers return a product, but it will take longer (since there are more off target hits to evaluate), may drop out for some species, and will become more sensitive to coverage (often requiring more data).

Primer specificity can be quantified in bits of information. Each position
in a primer constrains the set of sequences it can bind to, and the amount
of constraint is the specificity:

| Position type | States | Bits of specificity |
|---|---|---|
| Fully resolved (A, C, G, T) | 1 of 4 | 2.00 |
| 2-state degenerate (R, Y, S, W, K, M) | 2 of 4 | 1.00 |
| 3-state degenerate (B, D, H, V) | 3 of 4 | 0.42 |
| Unresolved (N) | 4 of 4 | 0.00 |

A 15 bp fully resolved primer with no mismatches has 30 bits of specificity.
Since sPCR requires both a forward and reverse primer to match, a fully
resolved primer pair provides 60 bits of specificity total.

By default sharkmer allows up to 2 mismatches at any position in each primer
(see `sharkmer --help-pcr`). Each mismatch makes one position effectively
unresolved, removing 2 bits of specificity per mismatch. Under default
settings, a fully resolved 15 bp primer pair therefore has 60 − 2×4 = 52
bits of effective specificity.

### Increasing specificity

There are three ways to add specificity to a primer:

1. **Reduce mismatches.** Lower the `mismatches` parameter (see
   `sharkmer --help-pcr`). Each mismatch removed adds 2 bits per primer.
2. **Reduce degeneracy.** Replace a more degenerate base with a less
   degenerate one — e.g. D (0.42 bits) → R (1 bit) → G (2 bits). 
3. **Lengthen the primer.** Increase the `trim` value (see
   `sharkmer --help-pcr`; default is 15) and add more resolved positions.
   Each fully resolved position adds 2 bits. Note that the effective primer
   length cannot exceed *k*−1 (default 18): the seed node for a primer of
   length *k* would only capture the first *k*−1 bases, silently dropping
   the 3' terminal base from the assembled amplicon. `trim` is automatically
   clamped to *k*−1 with a warning if set higher.

These approaches are complementary. If you increase degeneracy to
accommodate a broader clade (reducing specificity), consider lengthening
the primer to recover some of the lost bits. Whether this is possible
depends entirely on the sequence variation within the primer binding region
— the same classic primer design challenge as in wet-lab PCR.

### How much specificity is enough?

The number of bits of specificity required depends on:

- **Genome size.** Larger genomes have more potential off-target binding
  sites, requiring more bits to discriminate.
- **Genome complexity.** Repetitive regions can match primers at many
  locations, requiring more bits to resolve.
- **Target copy number.** High-copy targets (rRNA genes, mitochondrial
  genes) are easier to recover because their signal stands out above
  off-target noise. Single-copy nuclear genes require more specificity
  (or more reads) to reliably assemble.

The [node budget](README.md#node-budget) interacts with specificity: when
primers are less specific, more off-target seeds survive and consume graph
nodes. Raising the node budget can compensate to a point, but improving
primer specificity is generally more effective than increasing the budget.

## Panels

A **panel** is a YAML file under `panels/` that declares a named set of primer
pairs for *in silico* PCR. Panels are embedded into the sharkmer binary at
compile time via `include_str!()`, and can also be loaded at runtime from a
file via `--pcr-panel-file`.

## Panel file format

Panel files use **schema version 2**. A fully annotated example covering
every supported field is at [`panels/examples/reference.yaml`](panels/examples/reference.yaml).
A minimal working example:

```yaml
name: my_panel
schema_version: "2"
panel_version: 1.0.0
description: "Example panel."
clade: "Cnidaria"
taxon_id: 6073

maintainers:
  - name: "Samuel Church"
    orcid: "https://orcid.org/0000-0002-8451-103X"

changelog:
  - panel_version: 1.0.0
    date: 2026-04-05
    sharkmer_version: "3.0.0"
    changes: "Initial panel."

primers:
  - gene: "16S"
    compartment: "mitochondrion"
    gene_type: "rRNA_LSU"
    copy_number: "high_copy"
    forward_seq: "GRCTGTTTACCAAAAACATA"
    reverse_seq: "AATTCAACATMGAGG"
    min_length: 500
    max_length: 700
    expected_length: 580
    min_count: 2
    mismatches: 2
    trim: 15
    citation: "Cunningham and Buss 1993 https://doi.org/10.1016/0305-1978(93)90009-G"
    notes: "Amplifies portions of domains IV and V."

validation:
  samples:
    - accession: SRR9278435
      taxon: "Xenia sp."
      taxonomy: "Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Cnidaria; Anthozoa; Octocorallia"
      max_reads: [1000000]
      expected:
        "my_panel_16S": { min_identity: 0.97, length: 578, length_tolerance: 30 }
```

### Required fields

**Panel-level:** `name`, `description`, `primers`.

**Per-primer:** `gene`, `forward_seq`, `reverse_seq`.

### Optional panel-level fields

| Field | Description |
|---|---|
| `$schema` | JSON Schema URL for editor tooling. Set to `https://raw.githubusercontent.com/caseywdunn/sharkmer/main/schemas/panel/v2.json` to enable inline validation and autocomplete in VS Code (YAML extension) and JetBrains IDEs. |
| `schema_version` | Must be `"2"` for schema v2 panels. |
| `panel_version` | Semver string (e.g. `1.0.0`). Independent of sharkmer's version. |
| `clade` | NCBI-preferred name for the target taxon (e.g. `"Cnidaria"`). |
| `taxon_id` | NCBI Taxonomy ID for the clade. |
| `gene_prefix` | Overrides `name` as the output filename prefix. |
| `status` | `"experimental"` \| `"stable"` \| `"deprecated"`. |
| `source_url` | URL where this panel file can be obtained. |
| `license` | SPDX license identifier (e.g. `"CC-BY-4.0"`). |
| `citation` | Panel-level citation. |
| `notes` | Free-text notes about the panel. |
| `maintainers` | List of `{name, orcid, contact, notes}`. |
| `changelog` | List of `{panel_version, date, sharkmer_version, changes, notes}`. |
| `validation` | Validation block (see below). |
| `references` | Reference sequences for BLAST-identity validation. |

### Optional per-primer fields

| Field | Description |
|---|---|
| `region` | Sub-region of the gene (no `_`). Combined with `gene` to form `gene-region` in output names. |
| `index` | Integer (≥ 1) distinguishing multiple primer pairs for the same `(gene, region)`. Required when two entries share the same gene+region. Output name becomes `gene_index` or `gene-region_index`. |
| `compartment` | INSDC /organelle value. Absent = nuclear. `"mitochondrion"` \| `"plastid:chloroplast"` \| etc. |
| `gene_type` | `"protein_coding"` \| `"rRNA"` \| `"tRNA"` \| `"rRNA_SSU"` \| `"rRNA_LSU"` \| `"rRNA_5S"` \| `"ITS"`. |
| `copy_number` | `"single_copy"` \| `"low_copy"` \| `"high_copy"`. |
| `deprecated` | `true` to soft-deprecate (runs, warns). Default `false`. |
| `deprecated_by` | Name of the replacement entry. |
| `deprecated_reason` | Human-readable reason for deprecation. |
| `min_length`, `max_length` | Amplicon length search window (bp). |
| `min_count` | Minimum kmer count to seed the graph. |
| `mismatches` | Allowed mismatches per primer during kmer expansion. |
| `trim` | Bases to trim from each primer end (≤ k−1, default k=19). |
| `expected_length` | Expected amplicon length for reporting. |
| `citation` | Primer-level citation. |
| `notes` | Free-text notes about this primer pair. |
| `dedup_edit_threshold` | Levenshtein distance below which two products are merged (default 10). |

### Output naming

The output filename for each amplicon is `{prefix}_{gene_name}` where:

- `prefix` = `gene_prefix` if set, otherwise `name`.
- `gene_name` is derived from structured fields:
  - `gene` only → `gene` (e.g., `18S`)
  - `gene` + `region` → `gene-region` (e.g., `16S-V3`)
  - `gene` + `index` → `gene_index` (e.g., `CO1_1`)
  - `gene` + `region` + `index` → `gene-region_index` (e.g., `18S-V9_2`)

Character constraints: `gene` must not contain `_` (index delimiter); `gene`
must not contain `-` when `region` is also set (ambiguity); `region` must not
contain `_`.

### Field semantics

- `panel_version` is a semver string and is **independent of sharkmer's version**.
  Bump the patch digit for notes/citation/validation edits, the minor digit
  for new primers or changed expectations, and the major digit for breaking
  changes (renamed or removed genes).
- `expected_length` is the canonical expected amplicon length. It is distinct
  from the `min_length`/`max_length` *search* window, which should be wider
  to accommodate biological variation.
- `dedup_edit_threshold` controls the Levenshtein distance below which two
  output products are collapsed into one (default 10). Lower it (e.g. 0–2)
  for panels targeting complex samples where distinct but closely related
  products should be retained separately.
- `validation.samples[*].taxon` is the species (or lowest-rank identification)
  of the sample organism. This is a short human-readable label, distinct from
  the full lineage in `taxonomy`.
- `validation.samples[*].taxonomy` should be the NCBI lineage string for the
  sample organism (when the sample is not a metagenome). This makes it easy
  to see at a glance how much taxonomic diversity a panel has been tested
  against.
- `validation.samples[*].max_reads` is a list of read counts to run sharkmer
  at, processed in descending order. A single-element list (e.g.
  `[1000000]`) is fine; multiple values provide a sensitivity sweep.
- `validation.samples[*].expected` keys must match the **prefixed** gene
  names as they appear in output files (i.e., `{prefix}_{gene_name}`).
- Unknown fields are rejected at load time. If the Rust loader complains
  about a field, check for typos.

## Panel development cycle

The typical workflow for creating or modifying a panel is:

### 1. Start from an existing panel

Copy the closest existing panel in `panels/` to a new filename, or dump
a built-in panel with `--export-panel`:

```bash
sharkmer --export-panel cnidaria > panels/my_panel.yaml
```

Edit in place rather than starting from a blank file — the existing panels
carry conventions (mismatch counts, trim values, length windows) that are
tuned for sharkmer's graph assembly and are a reasonable starting point.

### 2. Edit primers

Add, remove, or tweak primers. Keep `gene` (and `region`/`index` if set)
stable across edits — the validation block's `expected:` keys derive from
these. If you rename a gene, update the `expected:` keys to match.

Bump `panel_version` in the panel header. Add a `changelog` entry describing
what you changed and why. Dates are in `YYYY-MM-DD`.

### 3. Declare validation samples

In the `validation.samples` list, declare one or more ENA/SRA accessions that
the panel should be tested against, with `max_reads` values and optional
`taxon` (species label) and `taxonomy` (NCBI lineage) strings. Leave
`expected` thresholds empty on first pass — the validator will populate them
on the next step.

```yaml
validation:
  samples:
    - accession: SRR9278435
      taxon: "Xenia sp."
      taxonomy: "Eukaryota; Opisthokonta; Metazoa; ... Xenia"
      max_reads: [1000000]
```

Choose samples that exercise the primers under realistic conditions: the
taxon group the panel targets, at a read depth where you would expect
recovery to succeed. One sample per major clade within the panel's scope
is usually enough. If you are actively tuning the panel, include a sweep
of read depths so you can see the sensitivity floor.

### 4. Run the validator

The validator reuses the `sharkmer-bench` conda environment, which
provides `blastn` (for amplicon identity validation), `ruamel.yaml` (for
round-tripping panel YAML in `--write` mode), and the other Python
dependencies shared with `benchmarks/`. If you have not created it yet:

```bash
conda env create -f benchmarks/environment.yaml
```

Then activate and run:

```bash
conda activate sharkmer-bench
python scripts/validate_panel.py panels/my_panel.yaml

# Or validate all panels at once:
scripts/validate_all_panels.sh

# Skip BLAST (faster, just checks recovery):
scripts/validate_all_panels.sh --no-blast
```

If `conda activate` errors with *"Run 'conda init' before 'conda
activate'"*, your shell has not been initialised for conda. Either run
`conda init zsh` (or `bash` — whichever shell you use), open a new
terminal, and try again; or, for a one-off without touching your shell
rc, source the conda shell script first:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sharkmer-bench
```

This runs sharkmer against each declared sample at each declared read depth
(using the shared cache under `benchmarks/data/cache/`), BLASTs the
recovered amplicons for identity (using a local `/db` BLAST database if
present, otherwise NCBI remote), compares recovery to any `expected`
thresholds, and emits a markdown report to `panels/validation_reports/`.

On first run with empty `expected` blocks, every gene will be reported as
unvalidated — that is the expected state. At the bottom of the report is a
**Suggested validation block**: a ready-to-paste YAML fragment containing
the observed identity, length, and tolerances.

Read the report. If the numbers look reasonable for your target taxon,
either paste the suggested block into the panel file by hand, or rerun
with `--write` to have the validator do it for you:

```bash
python scripts/validate_panel.py panels/my_panel.yaml --write
git diff panels/my_panel.yaml
```

`--write` only updates `expected` entries for samples already declared in
the panel — it will not add new samples. It rounds identities to 3 decimals
and sets thresholds slightly below observed values so normal run-to-run
variation does not cause false failures. It also updates
`validation.last_validated` with the current sharkmer version, panel
version, and date.

The validator writes reports to `panels/validation_reports/` by default. If
your panel lives outside the repo (e.g. a private panel under development),
use `--output-dir` to redirect reports to a directory of your choice:

```bash
python scripts/validate_panel.py ~/my_panels/arachnida.yaml \
    --output-dir ~/my_panels/reports/
```

This also works with `--write`, which edits the panel file in place regardless
of where it lives:

```bash
python scripts/validate_panel.py ~/my_panels/arachnida.yaml \
    --output-dir ~/my_panels/reports/ --write
```

### 5. Iterate on a single primer

When you are tuning one primer pair and do not want to rerun the whole
panel, use `--genes`:

```bash
python scripts/validate_panel.py panels/my_panel.yaml --genes 16S
python scripts/validate_panel.py panels/my_panel.yaml --genes 16S --write
```

The validator writes a temporary panel file containing only the requested
genes, runs sharkmer against that, and then updates only those genes in the
main panel file (for `--write`). Other genes' existing `expected` entries
are left untouched.

### 6. Commit

Once the panel validates at the depths you care about, commit the panel
file and its validation report together. The report is a small artifact
but it documents what "passing" meant on a specific sharkmer version,
which is useful provenance when debugging regressions later.

## Sharkmer version and panel stability

Sharkmer is under active development. Panel recovery may change slightly
across sharkmer versions as the PCR engine evolves — usually as
improvements, occasionally as tradeoffs where one gene gains at another's
expense.

To keep panels robust across these changes:

- `validate_panel.py --write` deliberately sets thresholds **below**
  observed values, absorbing small shifts.
- Validation reports record the exact sharkmer version they were produced
  with. When a panel that used to pass starts failing on a newer sharkmer,
  comparing reports pins down what moved.
- The `changelog` field in the panel is the place to record tradeoffs
  ("CO1 identity threshold lowered 0.99 → 0.97 after sharkmer 3.1.0 graph
  pruning changes") so future readers can tell a regression from a
  deliberate re-baselining.

If a panel that worked on an older sharkmer stops working on a newer one,
open an issue. Sharkmer aims not to silently reduce panel recovery, and
when it does, that is information we want.

## Contributing a panel

External panel contributions are welcome. 

Contributed panels must include a `validation:` block with samples for validating the primers, and a `references:` block with expected results. All primers in the panel should return at least some correct targets. If during development you find that some primers return no products, off-target products (according to BLAST), or are highly unreliable, remove them from the panel before submitting. Dead primers slow down runs (often by a lot) and confuse the user about what their expected results should be.

Open a PR that adds a new file
under `panels/`. The PR should include:

- The panel file with `version`, at least one `maintainers` entry, an
  initial `changelog` entry, and a populated `validation` block.
- The validation report generated by `validate_panel.py` for that panel.

A reviewer will re-run the validator against the declared samples and
confirm that the expected thresholds are reasonable.

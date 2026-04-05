# Panels

A **panel** is a YAML file under `panels/` that declares a named set of primer
pairs for in silico PCR. Panels are embedded into the sharkmer binary at
compile time via `include_str!()`, and can also be loaded at runtime from a
file via `--pcr-panel-file`.

This document describes the panel file format and the development cycle for
creating or modifying a panel. For tool usage, see [README.md](README.md).

## Panel file format

```yaml
name: cnidaria
version: 1.0.0
description: "Universal and custom primers targeting mitochondrial and nuclear products."

maintainers:
  - name: "Samuel Church"
    orcid: "https://orcid.org/0000-0002-8451-103X"

changelog:
  - version: 1.0.0
    date: 2026-04-05
    sharkmer_version: "3.0.0-dev"
    changes: "Initial versioned schema."

primers:
  - gene_name: "16S"
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
      taxonomy: "Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Cnidaria; Anthozoa; Octocorallia; Malacalcyonacea; Xeniidae; Xenia; unclassified Xenia"
      max_reads: [1000000, 2000000, 4000000, 8000000, 16000000]
      expected:
        "16S": { min_identity: 0.97, length: 578, length_tolerance: 30 }
```

### Required fields

**Panel-level:** `name`, `description`, `primers`.

**Per-primer:** `gene_name`, `forward_seq`, `reverse_seq`.

### Optional fields

**Panel-level:** `version`, `maintainers`, `changelog`, `validation`.

**Per-primer:** `min_length`, `max_length`, `min_count`, `mismatches`,
`trim`, `expected_length`, `dedup_edit_threshold`, `citation`, `notes`.

### Field semantics

- `version` is a semver string and is **independent of sharkmer's version**.
  Bump the patch digit for notes/citation/validation edits, the minor digit
  for new primers or changed expectations, and the major digit for breaking
  changes (renamed or removed genes).
- `expected_length` is the canonical expected amplicon length. It is distinct
  from the `min_length`/`max_length` *search* window, which should be wider
  to accommodate biological variation.
- `min_length`, `max_length`, `mismatches`, `trim`, and `min_count` are passed
  through to sharkmer's PCR engine. See the main README for their meanings.
- `dedup_edit_threshold` controls the Levenshtein distance below which two
  output products are collapsed into one (default 10). Lower it (e.g. 0–2)
  for panels targeting complex samples where distinct but closely related
  products should be retained separately.
- `validation.samples[*].taxon` is the species (or lowest-rank identification)
  of the sample organism, e.g. `"Xenia sp."` or `"Homo sapiens"`. This is a
  short human-readable label, distinct from the full lineage in `taxonomy`.
- `validation.samples[*].taxonomy` should be the NCBI lineage string for the
  sample organism (when the sample is not a metagenome). This makes it easy
  to see at a glance how much taxonomic diversity a panel has been tested
  against.
- `validation.samples[*].max_reads` is a list of read counts to run sharkmer
  at, processed in descending order. A single-element list (e.g.
  `[1000000]`) is fine; multiple values provide a sensitivity sweep.
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

Add, remove, or tweak primers. Keep `gene_name` stable across edits if you
can — the validation block keys on it.

Bump `version` in the panel header. Add a `changelog` entry describing what
you changed and why. Dates are in `YYYY-MM-DD`.

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
thresholds, and emits a markdown report to `panels/reports/`.

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

External panel contributions are welcome. Open a PR that adds a new file
under `panels/`. The PR should include:

- The panel file with `version`, at least one `maintainers` entry, an
  initial `changelog` entry, and a populated `validation` block.
- The validation report generated by `validate_panel.py` for that panel.

A reviewer will re-run the validator against the declared samples and
confirm that the expected thresholds are reasonable.

# Inclusive sequence delivery

User-approved direction, 2026-09-29: default behavior should deliver standard
products and bounded repeat-uncertain hypotheses together in the normal
per-target FASTA. [AGENTS.md](../AGENTS.md) records this product preference.
This supersedes the opt-in, separate-file interface evaluated at `f0b7949`;
the earlier [protocol](UNCERTAIN_PRODUCTS.md), [results](UNCERTAIN_PRODUCTS_RESULTS.md),
and evidence bundle remain historical checkpoints.

## Output contract

For every target with at least one emitted sequence, write one
`{sample}_{gene}.fasta`. Standard records come first in their existing order,
followed by deterministically ordered uncertain hypotheses. Every header has
exactly one explicit `assembly_status=standard` or `assembly_status=uncertain`
field. Standard means the existing assembly gates accepted that sequence;
it does not mean independent biological validation.

Uncertain records retain the detailed repeat-structure labels, sequence hash,
thresholds, source-path count, and search/retention flags. Their repeat length
remains unresolved and read support remains unevaluated by this output feature.
They use identifiers ending in `.uncertain`, distinct from ordinary identifiers
ending in a numeric product index. There is no separate uncertainty FASTA.

Emission is enabled by default. `--no-emit-uncertain` requests standard records
only. The existing `--emit-uncertain` spelling remains accepted as an explicit
request for inclusive output; specifying both flags is an error.

To select standard records with a shell command while preserving complete
FASTA records, including wrapped sequences:

```sh
awk '/^>/ {keep = / assembly_status=standard( |$)/} keep' sample_gene.fasta > sample_gene.standard.fasta
```

Replace `standard` with `uncertain` to select the other class. Output files
remain owned, checksummed members of the sample's publication transaction.
Reruns remove only verified, unchanged files from a previous owned generation,
including the earlier separate uncertainty files; modified or foreign files
remain protected.

## Counts, statuses, and consumers

Output manifests use schema version 2. Stats identify
`product_output_schema_version: 2` and
`product_output_layout: combined_classified`. Within each target:

- `n_products` and `product_lengths` describe all delivered FASTA records.
- `n_standard_products` and `n_uncertain_products` partition that total.
- `status: success` means at least one standard record was emitted.
- `status: uncertain` means only uncertain hypotheses were emitted; the
  ordinary assembly failure reason remains available.
- `status: fail` means no sequence was emitted and no target FASTA is created.

The optional `uncertain_products` assessment retains its own count, lengths,
source metadata, and limits. Its `output_file` points to the shared target
FASTA when uncertain records exist. The file has one transaction receipt.
Terminal summaries describe actual inclusive output and distinguish classes.

Updated validation first verifies the complete mixed FASTA against stats,
class labels, uncertainty metadata, and the transaction manifest. Standard
validation and reference bootstrap then use only the standard subset.
They do not count uncertain-only targets as standard recovery. Previous
schema-1 separate-file output and legacy files remain readable. Earlier
current-run consumers that only support manifest schema 1 reject the new
schema instead of silently counting all hypotheses as standard products.

## Assembly and verification scope

This changes default delivery and output accounting. It does not change graph
construction, search order, thresholds, repeat admission, accepted-path quota,
or existing retention limits. Inclusive output remains bounded and is not an
exhaustive enumeration of possible repeat structures.

The previous evaluation emitted 26 standard and 34 uncertain records across
three selected insect samples. Those are the fixed identity expectations for
this interface regression: 60 records by default, 26 with the opt-out, with
the same five of seven historical identities available as uncertain records.
They are predictions from the previous completed evaluation, not a claim that
the revised interface has already passed its tests.

Verify both Rust hash backends, formatting/Clippy, Python consumers, and the
reference audit. Regression controls must cover mixed, uncertain-only, empty,
and standard-only outputs; conflicting flags; custom target names; safe
migration from split schema-1 output; and uncertainty exclusion from standard
validation/bootstrap. Run six invocations (default and opt-out on each frozen
first-million R1 prefix) against the previous accepted outputs, comparing
literal sequences, original standard header fields, counts, threshold
diagnostics, and bounded uncertainty assessments. Header class fields and
uncertainty ID/output paths are intentional format changes. This is an output
compatibility regression, not a new repeat-correctness or performance study.

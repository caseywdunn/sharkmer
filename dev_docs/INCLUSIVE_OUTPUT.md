# Inclusive sequence delivery

The [sequence annotation conventions](SEQUENCE_ANNOTATIONS.md) define the
current schema-3 header names, values, and omission rules. The focused v3.2
revision preserves assembly, retention, and released ordinary-header fields;
breaking interface cleanup is deferred to v4.0 under #159.

User-approved direction, 2026-09-29: default behavior should deliver standard
products and bounded repeat-uncertain hypotheses together in the normal
per-target FASTA. [AGENTS.md](../AGENTS.md) records this product preference.
This supersedes the opt-in, separate-file interface evaluated at `f0b7949`;
the earlier [protocol](UNCERTAIN_PRODUCTS.md), [results](UNCERTAIN_PRODUCTS_RESULTS.md),
and evidence bundle remain historical checkpoints.

## Output contract

For every target with at least one emitted sequence, write one
`{sample}_{gene}.fasta`. Records without assembly warnings come first in their
existing order, followed by deterministically ordered warning-bearing
hypotheses. The former omit `assembly_status`; the latter have exactly one
`assembly_status=warning` and a nonempty `assembly_status_causes` list.
No warning means no warning was reported by the checks performed; it does
not mean independent biological validation.

Warning headers also carry `sample`, `gene`, `product`, and `length`, using
the existing ordinary-header field names. Causes are a sorted unique list of
`homopolymer_length_unresolved`, `repeat_structure_unresolved`, and/or
`high_coverage_repeat_suspected`. Detailed markers, hashes, source thresholds,
path counts, and search/retention observations remain in stats. No read-support
assessment is performed by this feature, so no `read_support` field is emitted.
Identifiers retain the `.uncertain` suffix for compatibility and remain
distinct from ordinary identifiers ending in a numeric product index.
There is no separate uncertainty FASTA.

Emission is enabled by default. `--no-emit-uncertain` requests records without
assembly warnings only. The existing `--emit-uncertain` spelling remains accepted as an explicit
request for inclusive output; specifying both flags is an error.

To select records without warnings in schema-3 output while preserving
complete FASTA records, including wrapped sequences:

```sh
awk '/^>/ {keep = !/ assembly_status=warning( |$)/} keep' sample_gene.fasta > sample_gene.no_warnings.fasta
```

Remove `!` to select records with warnings. For historical schema-2 output,
select `assembly_status=standard` or `assembly_status=uncertain` explicitly;
the old contract requires those fields. These shell filters select annotations,
not biologically validated sequences. Output files
remain owned, checksummed members of the sample's publication transaction.
Reruns remove only verified, unchanged files from a previous owned generation,
including the earlier separate uncertainty files; modified or foreign files
remain protected.

## Counts, statuses, and consumers

Output manifests use schema version 3. Stats identify
`product_output_schema_version: 3` and
`product_output_layout: combined_warnings`. Within each target, existing stats
names and status values remain for compatibility:

- `n_products` and `product_lengths` describe all delivered FASTA records.
- `n_standard_products` and `n_uncertain_products` partition that total into
  records without and with assembly warnings, respectively.
- `status: success` means at least one standard record was emitted.
- `status: uncertain` means only uncertain hypotheses were emitted; the
  ordinary assembly failure reason remains available.
- `status: fail` means no sequence was emitted and no target FASTA is created.

The optional `uncertain_products` assessment uses schema 3 and retains its own
count, lengths, source metadata, and limits. Each record's structured
`assembly_status_causes` must agree with its internal source evidence and FASTA
description. Its `output_file` points to the shared target FASTA when warning
records exist. The file has one transaction receipt. Terminal summaries
distinguish records with and without warnings and report target-wide search,
retention, and diagnostic-detail truncation.

Updated validation first verifies the complete mixed FASTA against stats,
warning annotations, uncertainty evidence, and the transaction manifest. Standard
validation and reference bootstrap then use only the standard subset.
They do not count uncertain-only targets as standard recovery. Previous
schema-1 separate-file output, schema-2 combined output, and legacy files
remain readable according to their own contracts. Missing status in schema 2
is invalid. Older current-run consumers reject unsupported schema 3 rather
than silently counting all hypotheses as standard products. Standalone and
renamed FASTA readers exclude both legacy uncertain labels and new warnings
from standard/reference-bootstrap products; malformed labels fail closed.

## Assembly and verification scope

The schema-3 revision changes annotations and format dispatch. It does not change graph
construction, search order, thresholds, repeat admission, accepted-path quota,
or existing retention limits. Inclusive output remains bounded and is not an
exhaustive enumeration of possible repeat structures.

The [schema-3 verification](../benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/README.md)
passes exact before/after agreement of sequence bytes, IDs, order, counts,
warning membership, and source evidence on the bounded repeat/clean fixture,
plus schema migration and consumer checks. The v3.2 subset of
[#159](https://github.com/caseywdunn/sharkmer/issues/159) is implemented;
the issue remains open for the v4.0 cleanup of released interface fields.

### Historical schema-2 verification protocol (2026-09-29)

The following preserves the original inclusive-interface test scope; it is
not a new acceptance requirement or a claim about schema-3 results.

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

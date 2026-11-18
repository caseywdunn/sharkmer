# Semantic conventions for sequence annotations

**Status: implemented and locally verified in v4.0 development; not released.**
This is the canonical design guide for sequence annotations in FASTA headers
and their associated structured evidence. Apply it when adding or revising
annotations. The [README](../README.md) describes implemented behavior;
[inclusive output](INCLUSIVE_OUTPUT.md) records the existing output contract.
The integrated verification (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
passes; Git integration and CI remain open under #159.

The v3.2 revision introduced warning annotations while preserving previously
released identity fields. The v4.0 revision renames those identity fields and
moves the ranking score into structured per-record evidence. It preserves
sequence bytes, identifiers, ordering, warning delivery and biological policy.
Implementation and the outstanding verification gate are tracked in
[#159](https://github.com/caseywdunn/sharkmer/issues/159) and the
[active plan](PLAN.md#v40--scalable-exact-counting-and-replay).
The [v3.2 verification](../benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/README.md)
retains the earlier warning-format checkpoint; it does not verify this migration.

## Principles

1. **Keep ordinary headers simple.** Emit annotations that help interpret or
   select a sequence. Omit routine status fields and assessments that were
   not performed, with an explicit meaning of omission for each field.
2. **Make relationships visible.** Related fields share a meaningful prefix,
   such as `assembly_status` and `assembly_status_causes`. Avoid unqualified
   names such as `causes` whose subject depends on context or field order.
3. **Use language meaningful to biologists.** Values describe the concern
   and its consequence for interpreting the sequence. Internal algorithms,
   graph identifiers, and pruning stages belong in detailed evidence.
4. **Use stable, documented values.** A reason code has one defined meaning.
   Keep codes readable and suitable for shell filtering; do not make users
   parse changing prose. Use lowercase words separated by underscores.
5. **Separate warning, reason, evidence, and scope.** A warning does not prove
   a sequence incorrect. An assessment being performed does not imply a
   favorable result. Search limits describe how completely alternatives were
   explored; they are separate from warnings about an individual sequence.
6. **Keep concise annotations backed by evidence.** FASTA carries the useful
   summary; the stats retain methods, measurements, provenance, and bounded
   diagnostic detail, associated unambiguously with the emitted record.
7. **Preserve inclusive delivery.** Warning-bearing sequences remain in the
   normal per-target FASTA by default, alongside sequences without warnings.
   Preserve an explicit opt-out and report search and retention limits.

These are semantic conventions for structured metadata, rather than an
event-log format. Shared namespaces and descriptive names have a precedent
in [OpenTelemetry's naming conventions](https://opentelemetry.io/docs/specs/semconv/general/naming/).
Record-level quality annotations also have a biological precedent in
[VCF's FILTER field](https://samtools.github.io/hts-specs/VCFv4.3.pdf), which
distinguishes passing filters from filters not being applied. These are
design references, not formats or vocabularies that Sharkmer implements.

## Agreed field behavior

The following meanings apply to output explicitly identified as implementing
this convention. They must not be inferred from arbitrary or older FASTA.

| Field | Value and meaning | When emitted | Meaning when absent |
| --- | --- | --- | --- |
| `assembly_status` | `warning`: one or more assembly warnings were detected for this sequence. | Only when a warning is present. | No assembly warning was reported by the checks performed. |
| `assembly_status_causes` | A nonempty comma-separated list of documented reasons for the assembly warning. | Required whenever `assembly_status=warning`; absent otherwise. | No warning reasons are reported. It must not be missing from a warning-bearing record. |
| `read_support` | `evaluated`: a defined read-support assessment was performed for this particular sequence. It says nothing about whether support was found. | Only after that assessment has actually been performed and its method and results recorded. | No read-support assessment is reported. |

`read_support=evaluated` is reserved for a future assessment feature. Schemas 3
and 4 do not emit or accept a read-support header field, since that assessment is
not implemented. Introducing it requires a defined assessment and compatible
format handling, not simply enabling the existing threading option.

Keep related fields adjacent for readability, but parsers must use their
names rather than their positions. Reason lists have deterministic ordering
and no duplicate values. A warning can have more than one reason.

Omission is field-specific. It is not a universal synonym for passed,
false, unknown, or not applicable. In particular, no warning does not mean
biological validation, and an absent read-support field does not mean
that reads contradict the sequence. The run metadata identifies the
software, output format, and checks performed.

The current optional `--read-threading` uses local graph evidence during
assembly; enabling it alone must not set `read_support=evaluated` on every
record. A future assessment needs a defined scope and results for the actual
emitted sequence. Related fields should use the `read_support_` prefix if
they are needed in the header; detailed methods and results belong in stats.

Illustrative headers (ordinary measurements omitted here for clarity):

```text
>sample_gene_0 sample=sample target=gene sequence_index=0 sequence_length=160
>sample_gene_0.uncertain sample=sample target=gene sequence_index=0 sequence_length=137 assembly_status=warning assembly_status_causes=repeat_structure_unresolved
```

## Warning causes

These stable, user-facing cause names summarize the exact internal evidence
retained in stats. They are conservative warnings, not proof of an incorrect
sequence or a biologically confirmed repeat.

| Internal evidence cause | FASTA and structured cause | Meaning |
| --- | --- | --- |
| `pre_prune_omitted_self_loop_node` | `homopolymer_length_unresolved` | The path encounters a run of one repeated base whose biological length is unresolved. |
| `pre_prune_cyclic_scc_node` | `repeat_structure_unresolved` | The path intersects a cyclic graph region, leaving repeat structure unresolved; the emitted path need not itself complete a loop. |
| `retained_collision_edge` | `high_coverage_repeat_suspected` | An unusually high-coverage connection into existing graph sequence suggests a repeat, without proving one. |

Choose final codes by their meaning to a biologist, and preserve the exact
internal trigger separately in the evidence. A simpler name must not make
a stronger biological claim. Do not invent a warning category for an
assessment that the software does not perform.

## Historical v3.2 scope and compatibility

Product output schema 3 uses `product_output_layout: combined_warnings`,
manifest schema 3, and warning-assessment schema 3. Warning records contain
exactly `sample`, `gene`, `product`, `length`, `assembly_status`, and
`assembly_status_causes`. Identity and length fields use the released ordinary
header names; `length` is emitted bases, and `product` is an ordinal within
its record class. Identifiers remain unique, including the `.uncertain` suffix.

Ordinary headers preserve `sample`, `gene`, `product`, `length`, `score`, and
all `kmer_count_*` fields. Only the unreleased `assembly_status=standard`
annotation is removed. K-mer counts are not independent read depth; `score`
remains an assembly ranking heuristic, not a confidence probability.

The new warning header replaces `hypothesis_length` with the existing `length`
field. Redundant `uncertainty` and `repeat_length` annotations are expressed by
the causes. The sequence hash, source thresholds/path counts, marker detail,
and target-wide search/retention limits remain in structured evidence rather
than FASTA descriptions. `assembly_status_causes` is also a structured list
on each warning record and must agree with its internal source evidence.
The existing stats names `n_standard_products`, `n_uncertain_products`, and
`uncertain_products` remain; they mean records without/with warnings and the
latter's evidence. A target's `status: uncertain` still means warning-only
output. Terminal summaries use warning terminology and report reached limits.

Readers preserve legacy/schema-1 split-output and schema-2 combined-output
contracts. Schema 2 requires explicit standard/uncertain labels; absence there
is invalid. Absence in schema 3 means no warning reported. Standalone readers
recognize both old uncertain labels and new warnings and exclude them from
standard validation/bootstrap. Malformed or unsupported labels fail closed.
Run metadata, receipts, and evidence bind the format to the particular output;
do not interpret arbitrary unlabeled FASTA as biologically validated.

## v4.0 interface and migration

Product output schema 4 uses `combined_warnings`, manifest schema 4 and
warning-assessment schema 4. These declarations distinguish the renamed fields
from historical output. Header order is stable for readability; consumers use
field names. `sample`, `target`, `sequence_index` and `sequence_length` are
required on every emitted record. Warning records add the two warning fields;
ordinary records retain the existing `kmer_count_*` measurements.

| Annotation | Meaning/type | Emission and omission |
| --- | --- | --- |
| `sample` | Sample name from this run; nonempty string. | Required on every record. |
| `target` | Target identifier, including an existing panel prefix when present; nonempty string. The target can be a spacer or span several genes. | Required on every record; it does not rename panel-definition fields or output filenames. |
| `sequence_index` | Nonnegative, zero-based ordinal within the target's ordinary or warning record class, after final selection. | Required; it is neither a locus count nor a haplotype identifier. Each class starts at zero. Use the FASTA identifier for joins. |
| `sequence_length` | Positive number of bases actually emitted. | Required; it does not resolve biological repeat length. |
| `kmer_count_mean`, `kmer_count_median`, `kmer_count_min`, `kmer_count_max` | Existing observed k-mer count summaries along the selected assembly path; same definitions and formatting as v3.2. | Retained on ordinary records. Absent on warning records, where these measurements are not fabricated. They are not independent read depth. |
| `assembly_rank_score` | Nonnegative, unrounded assembly ranking heuristic; not calibrated confidence or biological support. | In structured stats for ordinary sequences only, not in FASTA. Omission means no such score is reported for this sequence. |

### Structured per-record annotations

Each target result has `sequence_annotations` with one entry per emitted
record, in the same order as the combined FASTA. Empty output omits the list.
Entries have `sequence_id` (the exact FASTA identifier), `sequence_index` and
`sequence_length`. An ordinary sequence also has `assembly_rank_score`; warning
sequences omit that score. The exact score follows the retained typed assembly
record through sorting, deduplication, truncation and final numbering. It is
not reconstructed from a rounded header string.

The run's `sample` and target's existing `gene_name` locate the result; the
record ID joins the annotation to its FASTA. Warning sources and causes remain
in `uncertain_products.products`, joined by the same ID. Warning cause names,
source identities, bounded diagnostic detail, search and retention limits
retain their meaning. No new confidence or read-support assessment is implied.

### Mapping from released fields

| v3.1–v3.2 header field | v4.0 home | Compatibility |
| --- | --- | --- |
| `sample` | `sample` | Same value. |
| `gene` | `target` | Same identifier; panel definitions are unchanged. |
| `product` | `sequence_index` | Same ordinal and class scope. |
| `length` | `sequence_length` | Same emitted base count. |
| `score` | `pcr_results[].sequence_annotations[].assembly_rank_score` | Same ranking calculation, retained at full precision. Historical headers rounded to two decimals. |
| `kmer_count_*` | Same FASTA fields | Same values and formatting. |

Sequence identifiers, including existing `.uncertain` suffixes, output filenames,
`--no-emit-uncertain`, and existing stats count/status names are retained.
Changing those is a separate compatibility decision. New consumers read saved
schemas 1, 2 and 3 under their original rules and schema 4 under this dictionary.
They reject conflicting old/new aliases in one header. Schema 4 requires
complete, ordered annotation joins, valid score types and matching record
identities/lengths. Unsupported versions do not fall back to an older format.
Standalone FASTA readers recognize both index spellings and continue excluding
warning hypotheses from standard validation and reference bootstrap.

Update scripts to select fields by name and join stats by the FASTA identifier.
For example, use `sequence_length` for emitted bases and read ranking scores
from `sequence_annotations`, rather than expecting `score` in a header. The
warning filter in the [README](../README.md) applies unchanged to schema 4;
keep the schema-specific interpretation when reading historical output.

A default read-support assessment remains separate future inference work.

## Maintaining the contract

For each added or changed annotation, record its name, meaning, type or
allowed values, emission condition, omission meaning, scope, and evidence
source here. Distinguish agreed designs, proposals, and implemented behavior.
Changing the meaning of an existing value is an output compatibility change.

Implementation of this design must coordinate the FASTA writer, stats,
validators, output-format versioning, filtering examples, and relevant tests.
Existing schema-2 output requires explicit standard/uncertain labels; a
missing label in that format must not silently acquire the new meaning.
Preserve the interpretation of saved historical outputs and evidence.

Keep the user-facing README and CLI help synchronized with the emitted fields
and update this document's status as verification completes. Header
simplification must preserve sequence delivery, warning
provenance, and the distinction between an emitted hypothesis and a
biologically validated sequence.

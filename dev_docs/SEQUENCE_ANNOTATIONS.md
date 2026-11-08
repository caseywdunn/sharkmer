# Semantic conventions for sequence annotations

**Status: v3.2 implementation verified, 2026-11-08; v4.0 cleanup remains open.**
This is the canonical design guide for sequence annotations in FASTA headers
and their associated structured evidence. Apply it when adding or revising
annotations. The [README](../README.md) describes implemented behavior;
[inclusive output](INCLUSIVE_OUTPUT.md) records the existing output contract.

The v3.2 revision applies these conventions to the unreleased warning
annotations. Released ordinary-record fields and all sequence identifiers
remain unchanged. Breaking cleanup of released fields is explicitly scoped
to v4.0 below. These changes do not alter assembly, sequence delivery, or
release authorization.

Implementation is tracked in [#159](https://github.com/caseywdunn/sharkmer/issues/159)
and [PLAN.md](PLAN.md#sequence-annotation-revision-159).
The [verification evidence](../benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/README.md)
records exact compatibility on the bounded fixture and passing local checks.

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

`read_support=evaluated` is reserved for a future assessment feature. Schema 3
does not emit or accept a read-support header field, since that assessment is
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
>sample_gene_0 sample=sample gene=gene product=0 length=160
>sample_gene_0.uncertain sample=sample gene=gene product=0 length=137 assembly_status=warning assembly_status_causes=repeat_structure_unresolved
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

## v3.2 scope and compatibility

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

## v4.0 interface cleanup (deferred)

The remaining recommendations change fields already released in v3.1 and
could break scripts. They belong at the v4.0 major-version boundary, with an
explicit migration guide and preserved reading of historical outputs.
Exact replacement names remain proposals for that work.

| Current annotation | Recommended treatment |
| --- | --- |
| `length` | Consider `sequence_length` for the number of bases actually emitted. This is not a claim about biological repeat length. |
| `gene` | Prefer `target`, since assays can target spacers or span several genes. This proposal concerns output annotations, not panel definitions. |
| `product` | If the numeric index is retained, call it `sequence_index`. It does not count biological loci or haplotypes. |
| `score` | Preserve in detailed stats as `assembly_rank_score`, an internal ranking heuristic without a calibrated confidence interpretation. |
| Existing stats names and record identifiers | Review separately if user-facing names should align with warnings. Any change needs explicit compatibility treatment. |

Moving a released field out of a header must preserve its value in associated evidence.
Some ordinary-record measurements currently exist only in the header; provide
a structured home before removing them. Keep identifiers stable through this
migration, including existing `.uncertain` suffixes; changing identifiers is a
separate compatibility decision. Preserve the existing opt-out spelling for
compatibility. A default read-support assessment is separate future inference
work, not implied by a metadata cleanup.

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

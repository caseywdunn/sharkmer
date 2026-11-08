# Semantic conventions for sequence annotations

**Status: agreed design, 2026-11-08; output changes are not yet implemented.**
This is the canonical design guide for sequence annotations in FASTA headers
and their associated structured evidence. Apply it when adding or revising
annotations. The [README](../README.md) describes implemented behavior;
[inclusive output](INCLUSIVE_OUTPUT.md) records the existing output contract.

Current output still uses `assembly_status=standard|uncertain`, `causes`,
and `read_support=not_evaluated` on uncertain records. The convention below
specifies their intended replacement. It does not change which sequences
are assembled or delivered, and does not authorize a release.

Implementation is tracked in [#159](https://github.com/caseywdunn/sharkmer/issues/159)
and [PLAN.md](PLAN.md#sequence-annotation-revision-159).

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

Illustrative header descriptions under the new convention (other metadata
omitted here for clarity; the example reason remains provisional):

```text
>sample_gene_0
>sample_gene_1 assembly_status=warning assembly_status_causes=cycle_encountered
```

## Cause vocabulary still to finalize

The current labels expose implementation details. The following alternatives
were discussed; they are proposals, not a finalized enumeration:

| Current internal cause | Proposed wording | Naming question to resolve |
| --- | --- | --- |
| `pre_prune_omitted_self_loop_node` | `repeat_count_unresolved` or `homopolymer_length_unresolved` | The second identifies the homopolymer case more precisely; the first also overlaps the concern behind other repeat warnings. |
| `pre_prune_cyclic_scc_node` | `cycle_encountered` | Still requires graph knowledge. Explain the unresolved repeat structure; the emitted path need not itself complete a cycle. |
| `retained_collision_edge` | `high_coverage_repeat_suspected` | Communicate a suspected repeat from an unusually high-coverage connection into existing graph sequence, without asserting a proven repeat. |

Choose final codes by their meaning to a biologist, and preserve the exact
internal trigger separately in the evidence. A simpler name must not make
a stronger biological claim. Do not invent a warning category for an
assessment that the software does not perform.

## Other annotations to review

The header audit identified the following recommended changes. These remain
proposals to settle before implementation, rather than additions to the
agreed field contract above. Ordinary and warning-bearing records should
share the same basic identity and length fields.

| Current annotation | Recommended treatment |
| --- | --- |
| `length` and `hypothesis_length` | Use `sequence_length` for the number of bases actually emitted. This is not a claim about biological repeat length. |
| `gene` | Prefer `target`, since assays can target spacers or span several genes. This proposal concerns output annotations, not panel definitions. |
| `product` | If the numeric index is retained, call it `sequence_index`. It does not count biological loci or haplotypes. |
| `score` | Preserve in detailed stats as `assembly_rank_score`, an internal ranking heuristic without a calibrated confidence interpretation. |
| `kmer_count_mean`, `kmer_count_median`, `kmer_count_min`, `kmer_count_max` | Keep these measurements and their shared prefix. Define them as k-mer counts, not independent read depth. Do not invent values for records where these measurements were not made. |
| `uncertainty=repeat_structure`, `repeat_length=unresolved` | Consolidate into warning reasons once the final vocabulary fully expresses their meaning. |
| `marker_detail`, `thresholds`, `source_path_records` | Keep in associated diagnostic evidence rather than FASTA headers. Preserve exactly what was retained and what was omitted. |
| `enumeration`, `limit_scope`, `search_limits_reached`, `retention_truncated` | Keep target-wide search and retention information in target stats and summaries, without presenting it as a property of each sequence. |
| `sequence_sha256` | The name is already clear. Prefer structured provenance as its home, retaining an unambiguous link between each emitted record and its evidence. |

Moving a field out of a header must preserve its value in associated evidence.
Some ordinary-record measurements currently exist only in the header; provide
a structured home before removing them. Keep identifiers stable through this
migration, including existing `.uncertain` suffixes; changing identifiers is a
separate compatibility decision. Update user-facing summaries to explain
warnings, and preserve the existing opt-out spelling for compatibility.

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

When implementation is complete, update the user-facing README and CLI help
to describe the emitted fields in plain language, and update this document's
status. Header simplification must preserve sequence delivery, warning
provenance, and the distinction between an emitted hypothesis and a
biologically validated sequence.

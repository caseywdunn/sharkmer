# Default inclusive delivery results

Completed 2026-09-29 on `experiment/uncertain-products`, following the user's
request for simple, inclusive delivery recorded in [AGENTS.md](../AGENTS.md).
Standard products and bounded repeat-uncertain hypotheses now appear together
in each target's normal FASTA by default. Every header identifies its class;
`--no-emit-uncertain` requests only standard products. The existing
`--emit-uncertain` spelling remains accepted.

## Sequence availability

All six default/opt-out invocations passed against the exact accepted outputs
of the preceding [split-output evaluation](UNCERTAIN_PRODUCTS_RESULTS.md).
The full 21-target insect panel, k=19, two CPUs, first million R1 records,
unchanged bounded search/retention, and no read threading were used again.

| Sample | Standard | Uncertain | Default total | Uncertain fraction | Opt-out total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gryllus, SRR27962769 | 11 | 17 | 28 | 60.7% | 11 |
| Drosophila, SRR31887760 | 9 | 4 | 13 | 30.8% | 9 |
| Heliconius, SRR1057608 | 6 | 13 | 19 | 68.4% | 6 |
| **Total** | **26** | **34** | **60** | **56.7%** | **26** |

All standard sequences and original standard header fields match exactly;
their new class annotation is the intended addition. All uncertain sequences,
ordering, descriptions, and bounded assessments match the previous evaluation,
apart from their intentional ID and output-path changes. Input/count fields
and threshold diagnostics also match. The revised output format does not
introduce additional sequence hypotheses beyond those already evaluated.

Five of seven historical identities are available as uncertain records:
Gryllus CO1_1 (373 bp) and ITS_2 (978 bp), Heliconius ND1 (263 bp), and
Drosophila 12S (475 bp) and 16S_2 (590 bp). Drosophila 12S (487 and 499 bp)
remain unobserved under the unchanged search and retention bounds.
Availability does not validate repeat length or biological correctness.
These three selected samples do not estimate the uncertain fraction in all
datasets, and this interface regression makes no performance-equivalence claim.

## Interface and validation

One normal target FASTA contains standard records first and uncertain records
after them. Stable `assembly_status=standard` and `assembly_status=uncertain`
fields support the shell filters in [the output contract](INCLUSIVE_OUTPUT.md).
Both filters were checked against every generated target FASTA, including
wrapped sequences. The uncertain IDs end in `.uncertain` to avoid collisions
with standard IDs from custom targets.

Stats and transaction manifests identify the combined schema-2 format.
`n_products` counts all delivered records, with explicit standard and uncertain
subtotals. Uncertain-only targets have `status: uncertain` and retain the
ordinary admission failure reason. Updated validation verifies the full file
and then uses only standard records for validation, recovery reporting, and
reference bootstrap. Schema-1 split output remains readable and owned older
outputs migrate safely; modified or foreign files stay protected.

The compiled checks passed on the same source snapshot used for all six runs:

- Rust: 247 unit and 29 integration tests for each hash backend; formatting and
  Clippy also pass for both applicable configurations.
- Python: 197 tests, one existing skip, no failures.
- Reference provenance: 156/156 references across ten panels.
- Regression harness: three comparison controls, including wrapped shell
  filtering and rejection of metadata/sequence changes.
- Independent source review: seven additional synthetic consumer controls.

The run receipts record source/input/compiler identities, process constraints,
before/after checks, transaction receipts, and exact sequence comparisons.
Slurm job `27832506`, attempt `attempt-v1`, completed all checks and runs without
a failed verification attempt. See the [evidence bundle](../benchmarks/benchmark_results/v3.2-inclusive-output-20260929/README.md)
for the final independent outcome review and archive verification.

This change is a local experiment-branch checkpoint. It does not integrate
read threading, change standard admission, close reconstruction/validation
issues, or authorize a merge or release. `dev` remains at `956e771`.

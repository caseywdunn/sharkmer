# Developer documentation

Start here when joining development or resuming work. Current requirements,
implementation descriptions, and dated evidence have different roles.

The [benchmark companion](https://github.com/caseywdunn/sharkmer-benchmarks)
holds released reports, protocols and historical result archives. References
marked "local evidence" identify paths relative to the repository root in the
development workspace. Archive that evidence before deleting working data.

## Current guidance

| Document | Use |
| --- | --- |
| [PLAN.md](PLAN.md) | Active v4+ work, dependencies, resource contracts, acceptance gates, and outstanding distribution follow-up. |
| [COUNT_STORE_CONTRACTS.md](COUNT_STORE_CONTRACTS.md) | v4 identity, exact-count, replay, and resource contracts; implementation status is explicit. |
| [OVERVIEW.md](OVERVIEW.md) | Architecture of the released default and current v4 development paths; update with code changes. |
| [SEQUENCE_ANNOTATIONS.md](SEQUENCE_ANNOTATIONS.md) | Canonical annotation meanings, omission rules, and the schema-4 interface and migration. |
| [YAML_WORKSPACE.md](YAML_WORKSPACE.md) | Pinned YAML serializer/loader workspace, historical input/model owners, conservative bounds and audit requirements. |
| [HTTP_WORKSPACE.md](HTTP_WORKSPACE.md) | Current HTTP/library scope, ordinary failure handling, the minimal chunk-framing fix, and historical accounting evidence. |
| [CONTRIBUTING.md](../CONTRIBUTING.md) | Build, quality checks, branches, and release procedure. |
| [ROADMAP.md](../ROADMAP.md) | Release goals and longer-term scope. |
| [HANDOFF.md](HANDOFF.md) | Short resumption entry point and access to the former session history. |

The user-facing [README](../README.md) and [PCR guide](../PCR.md) describe usage.

## Release baseline and evidence

The v4 starting point is [v3.2.1](https://github.com/caseywdunn/sharkmer/releases/tag/v3.2.1).
The [maintenance checklist](RELEASE_v3.2.1.md) and
[publication receipts](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v3.2.1-release-20261109/PUBLICATION.md)
record its scope and verification. [REVIEW_v3.2.md](REVIEW_v3.2.md) preserves the
biological/performance assessment and the explicitly deferred studies under
#129. The [final comparison](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v3.2-final-check-20261108/RESULTS.md)
and [sequence-retention assessment](SEQUENCE_RETENTION_RESULTS.md) are the main
recovery evidence. Additional emitted hypotheses do not establish additional
validated loci or haplotypes.

The following reports retain the inputs, methods, limitations, and decisions
from their dated investigation. Read them as evidence for those checkpoints;
use the active plan for outstanding work and the current overview for code.
Downloaded archives restore the original paths for protocol and result files.

| Investigation | Documents |
| --- | --- |
| Released-baseline benchmark and follow-ups | [v3.2 versus v3.1](BENCHMARK_v3.2_vs_v3.1.md), [high-copy follow-ups](BENCHMARK_high_copy_followups.md), [reference gates](BENCHMARK_reference_gates.md), [failure analysis](FAILURE_ANALYSIS.md) |
| Reference origin and replacements | [Provenance](REFERENCE_PROVENANCE.md), [public replacements](PUBLIC_REFERENCE_REPLACEMENTS.md) |
| Inclusive delivery | [Output contract](INCLUSIVE_OUTPUT.md), [output results](INCLUSIVE_OUTPUT_RESULTS.md), [release protocol](INCLUSIVE_RELEASE_PROTOCOL.md), [release results](INCLUSIVE_RELEASE_RESULTS.md) |
| Repeat uncertainty | [Design](UNCERTAIN_PRODUCTS.md), [results](UNCERTAIN_PRODUCTS_RESULTS.md) |
| Sequence retention | [Protocol](SEQUENCE_RETENTION_PROTOCOL.md), [results](SEQUENCE_RETENTION_RESULTS.md), [Drosophila 12S assessment](DROS12S_RETENTION_RESULTS.md) |
| Read-backed audit | [Protocol](READ_BACKED_AUDIT_PROTOCOL.md), [results](READ_BACKED_AUDIT_RESULTS.md) |
| Withheld-path diagnostics | [Design](WITHHELD_PATH_DIAGNOSTICS.md), [results](WITHHELD_PATH_DIAGNOSTICS_RESULTS.md) |
| ITS cyclic-region audit | [Protocol](ITS_SCC_AUDIT_PROTOCOL.md), [results](ITS_SCC_AUDIT_RESULTS.md) |
| ITS read recruitment | [Protocol](ITS_READ_RECRUITMENT_PROTOCOL.md), [results](ITS_READ_RECRUITMENT_RESULTS.md) |
| Ordered read histories | [Experiment](READ_HISTORY_EXPERIMENT.md), [results](READ_HISTORY_RESULTS.md) |

[primer_binding.yaml](primer_binding.yaml) retains historical primer-binding analysis for failed targets.
The [benchmark artifact directories](https://github.com/caseywdunn/sharkmer-benchmarks/tree/main/benchmarks/benchmark_results) hold
frozen protocols, raw results, manifests, hashes, and execution/review receipts.
Preserve their contents when changing navigation or current guidance.

The v4 annotation prototype (local evidence: `benchmarks/benchmark_results/v4-sequence-annotations-20261112/README.md`)
retains a tested schema-4 migration patch, source bindings and native comparisons.
The migration is now applied to the working tree alongside the counter changes;
combined local verification (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
passes. Git integration and CI remain open under #159. Use the active plan for
completion status and the canonical annotation document for the current
development interface.

## Historical planning and design

- [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md) preserves earlier design exploration
  and release-era implementation descriptions. Check current claims against
  the overview and code before carrying them into a new design.
- The [v3.0 plan/session log](https://github.com/caseywdunn/sharkmer/blob/v3.2.1/dev_docs/PLAN_v3.md)
  is preserved at the release tag; its superseded implementation checklists
  are outside the active plan.
- The [v3.2 release plan](https://github.com/caseywdunn/sharkmer/blob/v3.2.0/dev_docs/PLAN.md)
  retains the release gates and decisions removed from the active plan.
- The [pre-cleanup handoff](https://github.com/caseywdunn/sharkmer/blob/17842be085821969981585b24dca1e540b0b6009/dev_docs/HANDOFF.md)
  retains the former session history and its original evidence links.

# Gryllus ITS cyclic-graph and read-evidence findings

Date: 2026-09-13. Baseline: `ae9189a`. Issues #148/#153/#129 remain open.
**Release remains on hold; production code and assembly gates are unchanged.**
This completes the [frozen SCC/read protocol](ITS_SCC_AUDIT_PROTOCOL.md),
following the [withheld-path localization](WITHHELD_PATH_DIAGNOSTICS_RESULTS.md).

## Answer

The two critical local traversals in the historical 978 bp Gryllus ITS_2
sequence are **not corroborated by this read assay**. Neither the historical
local strings nor their enumerated extra-cycle alternatives has a matching
Q20 read, even allowing one substitution. An independent exact recount also
finds zero literal matches before quality filtering. The retained CO1 positive
control reproduces its previous support.

This does **not** establish that the old ITS product is biologically wrong,
but it provides no positive evidence that its rejection is too conservative.
The relevant cycles remain in the graph after pruning: this is not simply
a path retaining obsolete markers from cycles that pruning removed.
Keep the SCC gate unchanged. Do not restore the historical product on the
basis of broad marginal coverage or its former bootstrap-panel match.

## Input and graph findings

The target is the exact historical ITS hypothesis with SHA-256
`23163ca76f97f5d1685f0fcf4bda6ffc798bdf58c63953416833e78a480af68c`.
It is an assembled hypothesis, **not an orthogonal reference**. No graph-derived
string is added to the active reference catalog.

A scratch-only instrumented build exports the actual stored pre- and
post-pruning graphs at **k=19, coverage threshold 4**. It uses the frozen
two-target Gryllus panel, two threads, chunks=0, no read threading, and the
same first 1,000,000 SRR27962769 R1 records: 150,000,000 bases in 150 bp reads.
The consumed prefix is 365,356,898 bytes, SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
Later appended R2 records are not consumed. Full-file and prefix hashes are
checked before/after execution; the 40 GiB address-space and 1,800-second
limits remain fixed. These diagnostic runs are not performance benchmarks.

The successful capture matches the prior focused diagnostic-on run for all
five count fields, ordered full FASTA records, and complete gene, threshold,
and diagnostic payloads. Instrumentation does not change product recovery.

| Stored graph | Nodes | Edges | Relevant cyclic component sizes |
| --- | ---: | ---: | --- |
| Before pruning | 3,112 | 3,058 | 63 and 85 nodes |
| After pruning | 1,313 | 1,323 | 63 and 85 nodes |

An independent SCC calculation confirms the full candidate path in both
graphs and its five cyclic-node starts at 93, 94, 276, 277, and 278. The
historical path itself does not revisit a node; membership in these cyclic
components, rather than an observed node revisit, causes withholding.
The stored graph is already restricted by extension and coverage rules;
it is not every sequence or count-supported edge in the sample.

## Local alternatives and read results

All coordinates below are zero-based, half-open, in the historical candidate.
The complete component occurrences plus fixed 21 bp flanks define two
nonoverlapping windows. Enumerating all directed anchor-to-anchor spellings
of at most 150 bp in each complete stored graph gives two local classes per
region. The pre/post sets agree exactly. Each enumeration visits 175 states,
hits no cap, and records an overlength frontier. Thus it is complete within
this finite read-length universe, **not over arbitrarily many repeat copies**.

| Region / local hypothesis | Historical coordinates | Window length | Exact Q20 records | Q20 records, at most one substitution |
| --- | --- | ---: | ---: | ---: |
| ITS region 1, historical traversal | `[72,133)` | 61 bp | 0 | 0 |
| ITS region 1, extra 63 bp cycle | Same anchors | 124 bp | 0 | 0 |
| ITS region 2, historical traversal | `[255,317)` | 62 bp | 0 | 0 |
| ITS region 2, extra 85 bp cycle | Same anchors | 147 bp | 0 | 0 |
| Retained 352 bp CO1, original control | `[289,349)` in CO1 | 60 bp | 27 | 28 |
| Lost 373 bp CO1, original control | `[289,370)` in CO1 | 81 bp | 0 | 0 |

The ITS windows fit the read length and have valid flanks. Their zeros are
zero matches, not matches suppressed solely by competition or a failed
three-read threshold. The CO1 positive has 27/28 distinct IDs and 15 distinct
projected footprints in each mode, with a three-ID/three-footprint witness.
Its original event definitions and all five competing candidate sequences
are preserved. Sensitivity counts include exact matches, not a second
independent replicate. Local uniqueness is relative to the frozen classes,
not the entire genome; records are not UMI-certified molecules.

An independently implemented exact-only scanner, frozen before its author
inspected the primary counts, rereads the same full prefix once. All exact
event counts, read IDs, matches, and projected footprints agree. Both ITS
classes in both regions have zero literal matches on either strand even
without Q20 filtering. Retained CO1 has 33 literal matches before quality
filtering and 27 after; lost CO1 has zero in both cases. Quality filtering
alone therefore does not explain the ITS **exact-match** zeros. This check
does not independently rescan substitution-tolerant matching without Q20.

## Why earlier broad ITS coverage is insufficient

The [earlier marginal assay](READ_BACKED_AUDIT_RESULTS.md#marginal-sequence-occurrence)
found 38/47 exact and 39/47 sensitivity windows, whose coordinate unions
cover 959/978 bp respectively. Those unions are not continuous read bridges
or evidence of full-length phase. In fact, the earlier windows neighboring
these critical regions, `[60,121)`, `[80,141)`, `[240,301)`, and `[260,321)`,
already had **zero records in both modes**. Coverage from other overlapping
windows cannot establish that the questioned junctions coexist on one read.
The current results resolve a missing structural assay, not a contradiction
with previously demonstrated junction support.

## Validation and attempt lineage

Sol implements the scratch Rust exporter and independent exact recount;
Terra implements the graph-to-local-assay builder; Astra independently reviews
the stages and evidence. Root freezes inputs/protocol, runs the bounded
capture/read pipeline, and integrates the report. Validation includes:

- Scratch instrumented default-backend Rust suite: 236 unit and 23 integration
  tests pass; formatting and Clippy pass. No production Rust changes here.
- Frozen harness suite: 27 tests pass across graph capture, local derivation,
  and read execution; independent exact recount: 7 tests pass.
- Independent SCC and bounded-walk recomputation agrees for both graphs and
  all four region/stage enumerations, without importing the builder.
- All six primary read commands complete with frozen source/input identities.
  Review checks their report bindings and all 55 CO1 ledger observations
  against the earlier raw-verified ledgers. Those 55 include repeated records
  across exact/sensitivity modes, not 55 independent molecules.
- Twenty-four direct synthetic checks using the actual ITS manifests cover
  both orientations and exact/substitution matching; no orientation or
  manifest-callability defect is found. The separate full-prefix exact
  recount verifies the reported exact zeros rather than trusting empty ledgers.

The archive preserves two superseded attempts. The first build receipt used
a logical source-map hash where a source-snapshot file-byte hash was required;
`build-v2` corrects the binding and rebuilds before the accepted capture.
The first graph wrapper expected abbreviated phase names rather than the
exporter's actual phase labels. Its Sharkmer command completed with unchanged
outputs, but wrapper validation failed. The corrected wrapper has a regression
test and runs into fresh `graph_capture_v2` outputs. All three exported graph
files are byte-identical between attempts. Failed receipts are not rewritten
as successes; no outcome-driven assay relaxation or cap increase occurs.

The [evidence archive](../benchmarks/benchmark_results/v3.2-its-scc-audit-20260913/README.md)
contains frozen methods, source/build lineage, complete graphs, alternatives,
raw outputs and read ledgers, receipts, and independent checks. Public input
lineage remains in the preceding archives; raw FASTQ and compiled binaries
are not duplicated.

## Decision and remaining uncertainty

This local test does not corroborate either ITS structure. It does not show
that a rare allele is absent, adjudicate all seven historical losses, or
establish a correct full 978 bp haplotype. More than one substitution,
unmodeled indels, different anchor contexts, sparse sampling, and walks longer
than one read remain outside its sensitivity. The 124/147 bp alternatives
also have fewer possible spanning read starts than the shorter windows;
their zeros are not equally powered negative tests. Same-input reads provide
local evidence, not an independent sample or an accessioned truth sequence.

The next useful experiment is a **separately frozen read-recruitment/alignment
diagnostic** around both SCCs: inspect coverage and sequence differences using
both flanks, explicit quality/orientation rules, and controls for repeat-copy
and cross-locus ambiguity. Evaluate substitution/indel variation without
silently treating a length-changing alignment as proof of the historical
copy number. Report an empty or conflicting recruitment result honestly;
additional coverage or longer/paired reads would be a separately labeled
new-input experiment. Do not tune production admission against these zeros.

Any subsequent read-supported reconstruction proposal still needs positive
evidence for the restrictions removed, true/collapsed-repeat and mixture
controls, a fresh high-copy regression comparison, and user review. #148's
reconstruction work, #153's loss disposition, and #129's independent validation
remain open. No release, tag, master merge, or gate exception is made here.

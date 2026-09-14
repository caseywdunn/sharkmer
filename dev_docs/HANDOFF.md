# Pause/resume handoff

Updated 2026-09-14. **Project paused at the user's request. No release approval.**
Start here in a new session; this handoff does not require the previous chat
or live coding agents. [PLAN.md](PLAN.md) remains the detailed execution plan.

## Saved state

- Work branch: `dev`. Latest scientific evidence commit:
  `539d2d81f954aff54cc57be3e6f7d10f6ba19700` (committed and pushed).
  The containing handoff commit adds documentation only.
- Package remains `3.2.0-dev`. The latest production-code change is `ae9189a`,
  adding bounded, hidden, opt-in withheld-path diagnostics without changing
  admission gates or emitted FASTA recovery. The subsequent ITS graph exporter
  is scratch instrumentation archived as evidence, not a production feature.
- Sol/Terra implementation and Astra independent review are complete for the
  latest audit. No agent implementation or measurement is awaiting completion.
- [CI for the evidence commit](https://github.com/caseywdunn/sharkmer/actions/runs/34779935567)
  passes Linux/macOS Rust tests, Clippy, formatting, and Python/reference
  validation. Local production checks pass 234 unit and 23 integration tests;
  the local Python suite runs 112 tests with one skip and no failures.

## What we know

The scoped released-v3.1.0 versus development comparison still has **71 -> 64
high-copy products**, not demonstrated recovery parity. Counts are unchanged;
runtime is broadly similar and RAM is not meaningfully reduced. Do not confuse
later diagnostic-on/off parity with parity against the released version.
The seven lost identities are not seven proven incorrect sequences.

The latest [ITS SCC/read findings](ITS_SCC_AUDIT_RESULTS.md) use **k=19,
threshold 4, the same first million 150 bp Gryllus R1 reads**. Both relevant
cycles survive pruning. Historical local windows of 61/62 bp and their
124/147 bp extra-cycle alternatives all have zero Q20 matches in exact and
one-substitution assays. Independent exact recount also finds zero before
quality filtering. Retained CO1 reproduces 27 exact / 28 sensitivity records
over 15 footprints. These results do not corroborate the historical ITS
junctions, but do not prove biological absence or full-product incorrectness.
Earlier broad marginal coverage did not bridge these junctions. **There is
no positive basis here for relaxing the SCC gate.**

Reference provenance is separate from product correctness. There are 156 active
source-verified references; assembled hypotheses are not independent truth.
The former 978 bp Gryllus ITS panel match was a bootstrap-product match, not
full-length accession support. Cross-individual SNPs/indels need not mean a
wrong product. Two plant trnV-atpE annotation/assay conflicts and an independent
Agalma ITS source remain unresolved; do not fill them with Sharkmer products.

## Where the evidence lives

Read these in order as needed; older reports preserve historical checkpoints:

1. [PLAN.md](PLAN.md): release sequence, completed work, dependencies, and gates.
2. [High-copy follow-ups](BENCHMARK_high_copy_followups.md): seven-loss comparison
   and k19/23/27/31 sweep. Longer k did not restore any of the seven exact losses;
   k19 remains the default.
3. [Reference provenance](REFERENCE_PROVENANCE.md) and
   [public replacements](PUBLIC_REFERENCE_REPLACEMENTS.md): independent evidence,
   quarantine, variation-aware assessment, and remaining reference gaps.
4. [Read-backed audit](READ_BACKED_AUDIT_RESULTS.md),
   [withheld-path diagnostics](WITHHELD_PATH_DIAGNOSTICS_RESULTS.md), then
   [ITS SCC/read audit](ITS_SCC_AUDIT_RESULTS.md): successive localization and
   structural tests, their frozen protocols, limitations, and archive links.
5. [Pre-release review](REVIEW_v3.2.md): remaining all-panel, held-out, dependency,
   and release gates; implementation/test totals in older sections are historical.

The [latest evidence archive](../benchmarks/benchmark_results/v3.2-its-scc-audit-20260913/README.md)
contains 205 members, including methods, scratch source/patch, graphs, read
ledgers, tests, receipts, and independent reviews. Its SHA-256 is
`14799e1b8da1c1efc8675647a8c5a26076edde9894e55107f37db0ec0ba6abc3`.
Issue comments also record the final independent archive verification and CI:
[#148](https://github.com/caseywdunn/sharkmer/issues/148#issuecomment-5655808445),
[#153](https://github.com/caseywdunn/sharkmer/issues/153#issuecomment-5655808502),
[#129](https://github.com/caseywdunn/sharkmer/issues/129#issuecomment-5655808567).

**Do not rely on `/tmp` surviving.** Raw multi-gigabyte FASTQ, compiled binaries,
and Cargo target trees are intentionally excluded; archived manifests record
public input lineage, full-file/prefix hashes, toolchain, and rebuild commands.
Restore and verify them before rerunning. Original absolute paths are provenance,
not trustworthy input identities. Three unavailable superseded helper-source
versions are explicitly documented in the archive; they are not the accepted
scientific producer versions. Extract into a fresh directory, verify
`SHA256SUMS` and member hashes, and never overwrite original outcomes. Archived
documentation snapshots describe their pinned commit, not later handoff edits.

## Resume from here

1. Check the current `dev` branch, local changes, remote updates, and open issue
   status. Do not assume this dated handoff is newer than subsequent work.
2. Propose and freeze a **local read-recruitment/alignment diagnostic** around
   both ITS SCCs before inspecting new results. Test coverage and substitution/
   indel variation using both flanks, quality/orientation rules, and repeat-copy/
   cross-locus controls. The next protocol and implementation do not yet exist.
3. Keep the same-input experiment distinct from adding coverage, paired reads,
   or longer reads. Missing matches remain unresolved; alignments must not
   silently establish repeat count or concatenate local support into full phase.
4. Only positive evidence for the actual restrictions removed can motivate a
   separately reviewed reconstruction/admission change. It needs appropriate
   repeat/mixture controls and a fresh high-copy regression comparison.

#148 read-supported reconstruction, #153 high-copy loss disposition, and #129
independent validation remain **open**. Independent held-out registration,
all-panel biological validation, dependency/security review, and user release
review are not waived. Do not merge to `master`, tag, or publish without the user.

The user prioritizes abundant organelle/rRNA recovery for v3.2; nuclear and
metagenomic optimization is deferred. v4.0 is exact bounded-memory counting,
v4.1 low-coverage nuclear recovery, and v4.2 metagenomic diversity. Incremental
counting is **retained as an isolated legacy path**, not deleted; removing its
machinery from standard sPCR belongs to the planned counting redesign.
Use Sol/Terra for bounded implementation/harness tasks and Astra for independent
review as previously requested. A new session may start fresh agents; no old
agent memory is required.

# Pause/resume handoff

Updated 2026-09-30. **Development integration and all 60 release-comparison
invocations are complete. `dev` includes inclusive output at `ed04bde`.
The focused Drosophila diagnostic found both missing 12S identities in search
completions, discarded by retention. Independent raw outcome review passes;
archive identity and verification are recorded in the evidence bundle.
No release or publication is authorized.**
Start here in a new session; this handoff does not require the previous chat
or live coding agents. [PLAN.md](PLAN.md) remains the detailed execution plan.

## Saved state

- Completed execution: [inclusive release results](INCLUSIVE_RELEASE_RESULTS.md)
  follow the [frozen protocol](INCLUSIVE_RELEASE_PROTOCOL.md). All 60 runs in
  the original ten-sample/five-panel matrix pass against pristine v3.1.0.
  Default delivery has 64 standard and 62 uncertain records (126 total);
  69/71 released high-copy identities are available (64 standard, five
  uncertain). Two high-copy and 12 deferred Yp2 identities remain absent.
  All ten median runtimes increase; four exceed the 5% review threshold,
  while no RSS increase exceeds 2%. Biological support remains unassessed.
  Evidence: `benchmarks/benchmark_results/v3.2-inclusive-release-20260929/`;
  scratch: `tmp/inclusive-release-20260929/`. Clean paired builds `build-v1`
  pass (Slurm `27910574`); restored inputs are verified in `availability.json`;
  `comparison-v1/execution.json` passes all 60 runs (Slurm `27936421`).
  Failed restoration/report attempts remain preserved. The bundle's final
  outcome review binds the report and documentation; `archive_receipt.json`
  records the archive identity and independent member verification.
- [Drosophila retention diagnosis](DROS12S_RETENTION_RESULTS.md) is complete
  and independently reviewed: 26 unique completed 12S paths, three retained
  and 23 marker-cap omissions. Exact 499/487 bp historical paths occur at
  ordinals 13/19 and exceed the remaining 287 slots in the fixed 1,561-marker
  run share. Normal output remains nine standard plus four uncertain; neither
  missing sequence is emitted in production. Review the bounded retention
  policy next, alongside the resource regressions; do not infer validated
  repeat structure or silently increase caps from this diagnostic.
- Integration is complete: `dev` fast-forwarded from `956e771` to `ed04bde`.
  `baseline/pre-inclusive-dev-20260929` preserves the former dev commit;
  `experiment/uncertain-products` preserves the evaluated checkpoint. The
  earlier instruction to hold `dev` fixed is superseded by this authorization.

### Completed inclusive-output checkpoint (`ed04bde`)

- Completed revision after `f0b7949`: standard and repeat-uncertain records belong
  together in the normal target FASTA by default, with explicit status headers
  and `--no-emit-uncertain` for standard-only output. [AGENTS.md](../AGENTS.md)
  records the user's simple, inclusive default preference. The
  [new output contract](INCLUSIVE_OUTPUT.md) versions the combined format and
  makes total stats counts match all delivered records, with separate class
  counts. All six default/opt-out runs pass: 60 records by default (26 standard,
  34 uncertain), 26 with the opt-out, and five of seven historical identities
  available as uncertain hypotheses. Both Rust backends pass 247 unit and 29
  integration tests, formatting and Clippy; Python runs 197 tests with one
  existing skip; 156/156 references verify. See [the results](INCLUSIVE_OUTPUT_RESULTS.md)
  and [evidence bundle](../benchmarks/benchmark_results/v3.2-inclusive-output-20260929/README.md)
  for source review, completed outcome review and archive receipts. The
  preceding split-output evidence below remains immutable history.

### Preceding completed checkpoint (`f0b7949`)

- Work branch: `experiment/uncertain-products`, created from the completed
  diagnostic branch at `3cc48a3`. The user prioritized emitting explicitly
  uncertain sequences over further production read-threading integration.
  [Issue #157](https://github.com/caseywdunn/sharkmer/issues/157) tracks the
  implementation and assessment; [the protocol](UNCERTAIN_PRODUCTS.md)
  defines the output and frozen comparison. `--emit-uncertain` now writes
  separate labelled FASTA and structured provenance for bounded actual graph
  paths; standard admission and threshold stopping remain unchanged. Both
  hash backends pass 246 unit and 26 integration tests, formatting and Clippy.
  Full Python suite: 184 tests, one existing skip; reference audit: 156/156
  verified across ten panels. Independent parser review: 18 methods,
  69 subcases. See the new comparison bundle's validation and review receipts.
- The [comparison checkpoint](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/README.md)
  contains the fixed first-million R1 inputs' identities, reviewed harness,
  successful baseline-v3/candidate-v1 build receipts, passing checks-v3, and
  all nine completed baseline/off/on runs. Standard sequences and counts match:
  11 Gryllus, 9 Drosophila and 6 Heliconius products. Opt-in emission adds
  17, 4 and 13 separately labelled uncertain hypotheses, respectively.
  Five of seven exact historical hypotheses are re-emitted; see
  [the results report](UNCERTAIN_PRODUCTS_RESULTS.md). The two earlier failing
  test-fixture attempts and passing focused probes are preserved. The bundle's
  final results review and archive receipts document evidence verification.
  Changes remain unmerged on the experiment branch.
- Toolchain blocker resolved at the user's request: project-local Rust/Cargo
  1.98.1, rustfmt 1.9.0, and Clippy 0.1.98 are installed. From this repository,
  run `source tmp/uncertain-products-20260927/toolchain/activate.sh` to set the
  isolated Cargo/Rustup homes and PATH. Compiled checks and both comparison
  builds completed with the same effective toolchain and configuration;
  `validation/toolchain_install.json` in the new evidence bundle records
  installer integrity and actual tool versions. Rust 1.86 remains unsuitable
  for existing baseline let-chain syntax; its failed probe is retained.
  Python checks use modules `Python/3.12.3-GCCcore-13.3.0` and
  `PyYAML/6.0.2-GCCcore-13.3.0`.
- Preceding branch: `experiment/its-read-threading`, created from clean `dev`
  `956e771437a2410a7bcbb1504b033c61dd1e11b6`. Keep `dev` as the fixed comparison
  baseline while evaluating the experiment. The preceding SCC evidence commit
  is `539d2d81f954aff54cc57be3e6f7d10f6ba19700` (committed and pushed).
- The preceding [read-recruitment diagnostic](ITS_READ_RECRUITMENT_RESULTS.md),
  committed as `69104d2`, is
  standalone Python; production threading and admission remain unchanged.
  Its frozen scanner, 29 tests, independent algorithm checks, and 13 wrapper
  tests pass. All three biological/control scans completed with matching
  input/source identities; see the new evidence bundle for receipts.
- The subsequent [ordered-history prototype](READ_HISTORY_RESULTS.md) is
  implemented and evaluated on the same branch. Its 25 tests, 861 independent
  checks and 20 wrapper tests pass; both frozen graph scans complete. The
  independent full-prefix/ledger review and complete evidence are archived.
  This is standalone Python, not a production threading/admission change.
- Package remains `3.2.0-dev`. Before the #157 change, the latest
  production-code change was `ae9189a`,
  adding bounded, hidden, opt-in withheld-path diagnostics without changing
  admission gates or emitted FASTA recovery. The subsequent ITS graph exporter
  is scratch instrumentation archived as evidence, not a production feature.
- The previous SCC implementation and independent review were complete before
  this resumption. Current experiment sources, review, and execution are saved
  without requiring live agent memory.
- [CI for the preceding evidence commit](https://github.com/caseywdunn/sharkmer/actions/runs/34779935567)
  passes Linux/macOS Rust tests, Clippy, formatting, and Python/reference
  validation. That checkpoint's local production checks passed 234 unit and
  23 integration tests; its Python suite ran 112 tests with one skip and no
  failures. The new isolated diagnostic's checks are recorded separately above.

## What we know

The [current release comparison](INCLUSIVE_RELEASE_RESULTS.md) has **64/71
released high-copy identities retained as standard and 69/71 available with
uncertain delivery**. All 60 runs pass with unchanged aggregate input counts.
The 126 delivered records contain 64 standard and 62 uncertain hypotheses;
62/126 (49.2%) are uncertain. No biological correctness assessment is implied.
All ten median runtimes increase; Rhopilema, Human, Liriodendron and Heliconius
exceed the 5% investigation threshold. No median peak-RSS increase exceeds 2%.

Five historical high-copy sequences are available as exact uncertain
hypotheses: Gryllus CO1_1 (373 bp) and ITS_2 (978 bp), Heliconius ND1 (263 bp),
and Drosophila 12S (475 bp) and 16S_2 (590 bp). The separate pre-retention
observer now establishes that Drosophila 12S 487/499 bp are also generated,
then discarded by the marker-retention budget. They remain absent from
production FASTA output. These are unresolved repeat hypotheses with
unevaluated full-product read support, not newly validated products.
Standard admission and biological reconstruction questions remain open.

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

The [new variation-aware recruitment](ITS_READ_RECRUITMENT_RESULTS.md)
finds 27 ITS1 and 28 ITS2 recruited records. ITS1 has one low-quality
extra-cycle alignment. ITS2 has two Q20 best alignments to the longer class,
but the two complete read sequences and their unfiltered placements are
identical. One Q20 alignment shortens its span with a terminal deletion to
exclude a Q11 base. This is one observed sequence pattern, not two independent
footprints or proof of repeat count. Both historical short classes still have
zero qualifying alignments. The unchanged CO1 controls reproduce 27/28
retained records and zero lost records. Keep the repeat gate unchanged.

The [ordered-history experiment](READ_HISTORY_RESULTS.md) preserves repeat
visits and linked branch choices that can disappear in aggregate edge/pair
counts. On the fixed prefix it retains 773 graph-matching records before
pruning and 284 afterward. All four classes remain represented, but all have
zero complete exact/Q20 witnesses. The three previously identified spanning
reads have identical histories before/after pruning: their specific exact-walk
breaks already exist in the earlier graph. The known gapped ITS2 pattern is
outside exact-walk sensitivity; this is not evidence of biological absence.

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
   [ITS SCC/read audit](ITS_SCC_AUDIT_RESULTS.md), then the
   [variation-aware recruitment](ITS_READ_RECRUITMENT_RESULTS.md), then
   [ordered read histories](READ_HISTORY_RESULTS.md): successive
   localization, structural, and alignment tests, with frozen protocols,
   limitations, and archive links.
5. [Uncertainty-output results](UNCERTAIN_PRODUCTS_RESULTS.md) and
   [comparison bundle](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/README.md):
   implemented opt-in output, compiled checks, nine-run standard parity,
   separately reported hypotheses, and exact historical-identity assessment.
   Consult the bundle's final results review and archive verification receipts.
6. [Pre-release review](REVIEW_v3.2.md): remaining all-panel, held-out, dependency,
   and release gates; implementation/test totals in older sections are historical.

The [new recruitment bundle](../benchmarks/benchmark_results/v3.2-its-recruitment-20260925/README.md)
contains the frozen protocol/input/source plan, source snapshots, complete
recruited-read ledger, unchanged CO1 reruns, and independent reviews. The
directory retains its 2026-09-25 start date; scanning completed 2026-09-26.
The prefix was restored from ENA and hash-verified, rather than assuming
original temporary files survived. Its source freeze SHA-256 is
`cb68a200f410a89980380cf7d4c290c3d825c7b501f43acf16736ec954abe600`.

The [ordered-history bundle](../benchmarks/benchmark_results/v3.2-read-histories-20260926/README.md)
adds the frozen design, standalone prototype, actual graph inputs, synthetic
controls, complete retained-read/run ledgers, and independent reviews.
Its source freeze SHA-256 is
`337e27452472e6bbe99371b0ee12e3d41709a0079c1677ead72ff249180c45e9`.
The directory retains its 2026-09-26 start date; final review completed
2026-09-27. Both bundles include per-member inventories and checksums.

The [preceding SCC evidence archive](../benchmarks/benchmark_results/v3.2-its-scc-audit-20260913/README.md)
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

The user selected default combined delivery, then authorized integration and
the [release comparison](INCLUSIVE_RELEASE_RESULTS.md). Integration, the
earlier [six-run regression](INCLUSIVE_OUTPUT_RESULTS.md), all 60 new comparison
runs and the independently reviewed Drosophila diagnostic are complete.
Use the new bundle's final review and archive receipt for the saved evidence.
Preserve old evidence and keep production assembly/repeat gates unchanged.

1. Check `dev`, local changes, and the new inclusive-release benchmark
   bundle. #157 default combined delivery, compiled validation and all six
   default/opt-out regressions are complete. Check the bundle's independent
   outcome review, archive receipt, and checksums when reconstructing
   evidence. Preserve failed attempts and old evidence; do not rerun
   successful frozen comparisons without cause.
2. Default/output and development-integration decisions are already made;
   do not ask the user to decide them again. Report the new complete release
   comparison separately from the three-sample interface regression. Standard
   recovery, uncertain availability, bounded non-observations and biological
   correctness remain distinct. Keep the saved `956e771` baseline ref intact.
   The [completed read-history prototype](READ_HISTORY_RESULTS.md) remains a
   diagnostic; production integration under #145/#116 is a separate future
   evidence/reconstruction task, not required for inclusive delivery.
3. The next production decision concerns bounded retention: completed useful
   uncertain sequences can be lost because repeated marker provenance consumes
   a fixed per-target share. Assess compact provenance or another explicit
   bounded policy before changing caps. The observed 12S omission is not a
   search failure or biological validation result. Any new implementation
   needs its own controls and comparison.
4. Investigate the four current runtime flags. Preserve the original three
   pairs and distinguish CPU, filesystem/publication and other waiting costs.
   Compare suitably frozen versions to isolate causes; this release-baseline
   matrix cannot assign the cost specifically to uncertain emission.
5. For later reconstruction, specify variant-aware graph mapping: preserve
   competing indel placements and missing graph sequence as uncertainty.
   More independently positioned support and linkage between regions remain
   missing. Freeze any additional-coverage, paired-read, or longer-read study
   separately. Local alignments must not silently establish repeat count or
   concatenate support into full phase.
6. Only positive evidence for the actual restrictions removed can motivate a
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

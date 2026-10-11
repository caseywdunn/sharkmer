# Pause/resume handoff

Updated 2026-10-11. **The authorized plant correction, runtime investigation
and local validation are complete. No release was initiated.** The user will
have collaborators assess performance before release. The user has authorized
committing the reviewed candidate and integrating it into **`dev`**. Collaborators
should use `dev` and record `git rev-parse HEAD` with their results. The tested
source originated on `experiment/sequence-first-retention`, base `1ea21c9`;
the earlier freezes and completion receipts describe that pre-integration state.

Plant panel 2.0.0 uses published 6R and two exact public references. In the
eight reviewed before/after runs, Liriodendron trnV-atpE changes from 41
uncertain hypotheses to one exact-reference 909 bp standard product; Acer
changes from no product to one 929 bp standard product matching the related
Acer opalus reference. Both synthetic positives pass; all nine other targets
retain their sequences, statuses and diagnostics. Assembly, admission,
retention and inclusive output policy are unchanged.

The derived inventory combining these two corrected plant runs with eight
saved samples contains 287 records (66 standard, 221 uncertain), with all 86
release identities retained, zero lost and 201 new. It represents 75/118
targets. This combines different run dates and a changed plant assay; it is
not a fresh ten-sample final-build benchmark. The original unchanged-assay
comparison remains 326 records and 240 new identities.

The four-flag runtime investigation and bounded Human output-location check
are complete and independently reviewed. The Human diagnostic shows a
release-to-retention wall increase of 0.43 s with NFS output and 0.04 s locally;
historical causality and final-build performance acceptance remain open.
Minimum rand/webpki security patches are applied. All 11 final local gates
pass: 280 Rust tests per backend, 200 Python tests with no skips, 158 public
references, panel/schema/fixture checks and all eight plant runs. The
[independent terminal review](../benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/independent_results_review.json)
confirms results and unchanged source bindings.

Start with [the current review](REVIEW_v3.2.md#current-assessment-2026-10-11)
and its [tested source archive](../benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/candidate-source-v2.tar.gz)
for collaborators. Performance acceptance, wider/held-out biological
validation and an independent Agalma ITS reference remain outstanding;
the review distinguishes these scopes. The authorized integration into `dev`
does not authorize a merge to `master`, release tag or package publication.
No release is authorized.

## Earlier checkpoints

Completed prior checkpoint (2026-10-10): **Recommendations 1–3 are complete: the implementation is
frozen, the full 118-row comparison is recorded, and all 154 hypotheses on
six newly represented targets have bounded reference/read assessments.**
The work remains on `experiment/sequence-first-retention` at `1ea21c9`;
`dev` remains the inclusive baseline at `b5a6c7e` (including `ed04bde`).
Start here in a new session; this handoff does not require the previous chat
or live coding agents. [PLAN.md](PLAN.md) remains the detailed execution plan.

At that checkpoint, [findings](SEQUENCE_RETENTION_RESULTS.md#complete-target-comparison-and-six-new-targets-2026-10-10)
and [complete comparison CSV](../benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/final-v1/target_matrix_compact.csv):
within the matched ten-sample scope, all **86 v3.1.0 sequences are retained,
zero are lost, and 240 are new**. Output was 64 standard plus 262
uncertain records. Targets with any output increase 68 → 74 of 118; the
remaining 44 have no output, often after reaching search bounds. The older
200-gain comparison uses the intermediate inclusive build as baseline.

One newly represented target has strong corroboration: the 812 bp
Liriodendron atpB-rbcL hypothesis exactly matches two annotated public
plastomes, and all 37 tested local tiles occur at Q20 in the original reads.
The other 153 hypotheses on the six new targets have no qualifying full-query
reference match or tested exact local read string. All tiles of all eight
positive controls reproduce. These bounded negatives do not prove the
sequences false; local coverage also does not phase an entire product.
Four scans took 18.97–37.73 seconds each within fixed limits.

The trnV-atpE interpretation is resolved: the effective primers map from
trnV into atpB, including atpE, over 1,359–1,379 bp on the inspected plastomes;
no alternative pair falls within the configured 400–1,000 bp range.
Primer identity is corroborated by the primary literature. Correcting the
intended assay and panel under #156 is a separate next task. Keep the current
inclusive delivery and uncertainty labels. Then investigate the four earlier
runtime flags and finish release review. Production code and panels did not
change in this diagnostic. Evidence and documentation remain uncommitted;
the branch is not merged or released.

Preceding read follow-up (2026-10-09): [the findings](SEQUENCE_RETENTION_RESULTS.md#bounded-read-follow-up-2026-10-09)
resolve the scoped Gryllus linkage and Drosophila recovery checks. Among 32
Gryllus ITS gains, 11 share a complete nearby Q20 flanked span observed on
one saved read; four more share the complete variable core with incomplete
outer flank. The remaining 17 hypotheses share five nearby combinations
unobserved in the selected 25-read ledger, not proven absent from the sample. All full products remain
unphased. For all 26 Drosophila 12S hypotheses, including historical
475/487/499 bp sequences, none of 117 tested longer local strings occurs
exactly in the original million-read prefix, even without quality filtering;
all 16 CO1 controls reproduce. A separate recount finds every candidate
19-mer in that same prefix at least twice. Thus observed short graph pieces
do not establish their longer combinations. Independent reviews pass.
No concrete production defect was demonstrated; no production code changed.
Keep inclusive uncertainty labels and proceed to remaining runtime and release
validation work. Full repeat reconstruction remains separate. New evidence
and documentation are currently uncommitted.

Post-comparison update (2026-10-02): an independent audit of the exploratory
local BLAST check of all 200 gained uncertain identities passes. None has a
single HSP at both 90% query coverage and 90% identity in the pinned 121-source
catalog or 61 relevant panel references. This limited absence of qualifying
reference support is not a biological rejection or read-support assessment.
See [the results](SEQUENCE_RETENTION_RESULTS.md) and the
[BLAST review](../benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/blast-gains-v3/README.md).

Resumption update (2026-10-09): the remaining question is the biological
support of the additional hypotheses. The loss/recovery history is confirmed;
the available logs do not substantiate a later analysis stopped for excessive
computation. Slurm accounting confirms the corrected 60-invocation comparison
(`28037442`) completed in 1:03:27. The preceding comparison (`28034535`) was
cancelled after 35:10 because it displaced an old plant identity. The earlier
ITS recruitment and ordered-history jobs completed in 3:56 and 2:56,
respectively. These are whole-job elapsed times, separate from the measured
1.8% assembly overhead. Do not infer that an unrecorded expensive attempt
could not have occurred.

A new [saved-sequence triage](../benchmarks/benchmark_results/v3.2-gain-triage-20261009/report.json)
binds all 200 gains to their saved FASTA hashes and compares sequences within
each sample/target without assembly or read scans. Every gain has another gain
with identical sequence outside one window of at most 40 bases; 143 have such
a neighbor within five bases. All 32 Gryllus ITS gains meet the latter bound.
These are pairwise neighborhoods, not 200 independent variant sites or a
shared 40-base interval across all products. Exact canonical 19-mer sets form
196 groups; only four Yp2 pairs share a set, and all 200 sets differ at k=31
and k=51. Thus identical short-sequence content alone does not explain most
of the extra output. This describes sequence relationships, not biological
correctness, target identity, or independent loci.

The [target summary](../benchmarks/benchmark_results/v3.2-gain-triage-20261009/target_summary.tsv)
and [nearest variant windows](../benchmarks/benchmark_results/v3.2-gain-triage-20261009/nearest_windows.tsv)
provide a bounded starting point for read-support work. First inspect flank
uniqueness and alternate indel placements, select a few informative contrasts
(Gryllus ITS is a compact first case), and compare their quality-aware local
support in the existing fixed read prefixes with positive controls. Retain
ambiguous placements and distinguish local support from full-product linkage.
Some windows reach sequence ends and lack two complete flanks. This initial
triage made no read-support or production change; unsupported does not mean
proven spurious. The subsequent bounded pilot is recorded below.
Seven focused diagnostic controls pass. A separate calculation reproduces
the 200 windows, the 40-base maximum, and the 143 five-base neighborhoods;
all window substitutions reconstruct their paired sequences. Source and
result checksums were verified. The new files remain uncommitted.

### Bounded Gryllus ITS gain pilot (2026-10-09)

The user authorized the next read-support investigation. The
[pilot summary](../benchmarks/benchmark_results/v3.2-its-gain-pilot-20261009/summary.json)
reports one exact/Q20 scan of the original first million Gryllus R1 reads,
with matching input SHA-256. The 32 gained ITS hypotheses' nearest-neighbor
contrasts reduce to eight literal windows with complete 21-base flanks;
two unchanged historical CO1 controls accompany them. All 32 gained full
sequences and five historical control sequences remain in the finite
comparison context, with every original window membership preserved.

All eight ITS strings are observed, with 11, 8, 2, 3, 4, 4, 3, and 5 Q20
exact matching records respectively. Their union contains 25 distinct read
IDs; overlapping-window counts cannot be summed as independent observations.
The CO1 controls reproduce 27 records for the retained local structure and
zero for the longer alternative. Each ITS string occurs in 10–12 frozen
candidate templates, so none uniquely supports a full candidate. Two reads
contain strings from both surveyed neighborhoods, providing limited local
co-occurrence. The full variable extent across all 32 gains requires
424–426 bases including flanks, longer than these 150-base reads. These
observations do not establish biological polymorphism, target specificity,
independent molecules, full-product phase, or correct repeat structure.

Execution completed in 28.66 seconds with 16,896 KiB peak RSS inside the
existing one-CPU Slurm allocation, below the fixed 300-second/2-GiB address
space limits. Input/source/output receipts and complete matching-read evidence
are saved under the pilot's `run-v1/`. Independent review reproduces all
event matches, quality tests, coordinates, template placements, and counts
from the retained ledger without repeating the FASTQ scan. No production
code or admission policy changed. The next diagnostic question concerns
which combinations are linked by observed reads and which remain unresolved;
the tested local strings cannot all be dismissed as absent from the reads.

## Saved state

- [Sequence-first retention results](SEQUENCE_RETENTION_RESULTS.md) follow the
  [protocol](SEQUENCE_RETENTION_PROTOCOL.md). The corrected branch delivers
  64 standard and 262 uncertain records across the ten samples (326 total),
  versus 64 + 62 on `dev`. All 71 released high-copy candidate identities
  are now available (64 standard, seven uncertain); all 15 deferred identities
  are available as uncertain, including 12 previously absent Gryllus Yp2.
  The first candidate was rejected after it displaced one old plant identity;
  its partial run is preserved. The corrected design reserves old detailed
  paths independently of extra summary-only sequences. Focused plant/Drosophila
  replays, an identical real-input standard-only/hidden-diagnostic opt-out
  control, both hash backends, Python consumers, the reference audit, and the
  completed three-pair matrix pass. Sum of median wall times rises 1.8% from
  the old inclusive build; no sample crosses 5% wall or 2% RSS review flags.
  Five genes still report retention truncation; these are bounded hypotheses,
  not validated repeat structures. See the
  [evidence bundle](../benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/README.md)
  and scratch `tmp/sequence-retention-20261001/` for raw receipts. The branch
  is not merged into `dev` or released.
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
  and 23 marker-cap omissions in the fixed inclusive baseline. Exact 499/487 bp historical paths occur at
  ordinals 13/19 and exceed the remaining 287 slots in the fixed 1,561-marker
  run share. At that checkpoint normal output was nine standard plus four
  uncertain, with neither missing sequence emitted. The subsequent bounded
  retention implementation above recovers both on the experiment branch;
  the new read follow-up keeps their biological structure unresolved.
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

The fixed `dev` [release comparison](INCLUSIVE_RELEASE_RESULTS.md) has **64/71
released high-copy identities retained as standard and 69/71 available with
uncertain delivery**. All 60 runs pass with unchanged aggregate input counts.
The 126 delivered records contain 64 standard and 62 uncertain hypotheses;
62/126 (49.2%) are uncertain. No biological correctness assessment is implied.
All ten median runtimes increase; Rhopilema, Human, Liriodendron and Heliconius
exceed the 5% investigation threshold. No median peak-RSS increase exceeds 2%.

Five historical high-copy sequences are available as exact uncertain
hypotheses: Gryllus CO1_1 (373 bp) and ITS_2 (978 bp), Heliconius ND1 (263 bp),
and Drosophila 12S (475 bp) and 16S_2 (590 bp). The separate pre-retention
observer established that Drosophila 12S 487/499 bp are also generated,
then discarded by `dev`'s marker-retention budget. The isolated sequence-first
branch now emits both, along with the previously absent Yp2 identities, with
explicit uncertainty labels. These are unresolved repeat hypotheses with
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
3. Bounded retention and the scoped Gryllus/Drosophila read follow-up are
   complete; do not repeat these checks without a concrete remaining question.
   Read [the final findings](SEQUENCE_RETENTION_RESULTS.md#bounded-read-follow-up-2026-10-09)
   and the linked independent reviews. Some nearby Gryllus combinations occur
   on single reads. Drosophila candidates contain observed 19-mers but lack
   the tested longer exact strings. Full-product correctness remains open,
   and no demonstrated production defect calls for a targeted fix. Preserve
   inclusive, labelled delivery. Integration into `dev` remains a separate
   decision; the next development work is the outstanding runtime and
   validation review below.
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

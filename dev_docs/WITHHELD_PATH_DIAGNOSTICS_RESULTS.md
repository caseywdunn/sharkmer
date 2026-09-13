# Withheld-path diagnostic results

Date: 2026-09-13. Issues #148, #153, and #129. **Release remains on hold.**
This diagnostic-only change follows the [read-backed audit](READ_BACKED_AUDIT_RESULTS.md)
and the [preregistered twelve-run protocol](WITHHELD_PATH_DIAGNOSTICS.md).
No assembly gate, emitted-product policy, reference, or counting algorithm
changes. The comparison baseline is development commit `58ca2a6`, **not
released v3.1.0**; this is not a resolution of all previously reported losses.

## What changed

The hidden `--diagnose-withheld-paths` flag retains bounded actual completed,
in-range, repeat-withheld candidate sequences in transaction-owned YAML stats.
Each retained record includes its exact oriented sequence/SHA-256 and every
encountered self-loop, pre-pruning cyclic-SCC, and retained-collision marker,
with typed oriented sequence identities and zero-based half-open spans.

Default-off output omits this payload and performs no per-path diagnostic
sequence reconstruction, hashing, or marker-list allocation. Diagnostic
retention never admits a path, consumes an eligible-product quota, or stops
the existing search. Whole records are retained or explicitly omitted; caps
apply independently to paths, sequence bases, and marker occurrences per
threshold, gene, and run. Static panel-index allocations avoid scheduling
dependent retention budgets. These are payload caps, not process-RSS limits.

**Diagnostic sequences remain unsupported hypotheses, not public references
or recovered products.** Marker footprints explain current withholding;
their unions are not necessarily complete ambiguous intervals with
distinguishing flanks, and are not automatically valid read-bridge assays.

## Compatibility and descriptive performance

All twelve commands exit zero, validate their complete output transactions,
and preserve frozen sources, executables, panels, and input identities.
Each sample's baseline-off, new-off, and new-on full-panel runs have identical:

- All five read/base/k-mer count fields.
- Ordered FASTA records, full headers, and exact sequence bytes/hashes.
- Complete gene outcomes, threshold ordering, and pre-existing aggregate
  diagnostics, after removing only the new optional payload.

All three focused-panel comparisons also preserve these properties for their
selected targets. Full-panel recovered products remain **11 Gryllus, 9
Drosophila, and 6 Heliconius**, each from the same number of successful targets.
The new diagnostic mode neither recovers nor removes an emitted sequence.

Inputs are the exact previously consumed first 1M R1 records per accession;
settings are k=19, two threads, chunks=0, unpaired input, and no read threading.
The staged files' later R2 data is not consumed. Each command uses prefix-only
prewarming outside the timed interval, CPU affinity `0,1`, a 40 GiB address-space
limit, and a 1,800-second timeout. No builds/tests run during timed invocations.
The observed host is x86-64 Linux with an Intel Core i7-8700, not Apple Silicon.

| Sample | Baseline off, seconds | New off, seconds | New on, seconds | Baseline / new off / new on peak RSS, GiB |
| --- | ---: | ---: | ---: | --- |
| Gryllus, SRR27962769 | 98.614 | 98.756 | 98.057 | 4.255856 / 4.255733 / 4.255661 |
| Drosophila, SRR31887760 | 59.747 | 59.343 | 59.279 | 2.130791 / 2.130669 / 2.130905 |
| Heliconius, SRR1057608 | 44.175 | 43.760 | 43.693 | 2.130920 / 2.130653 / 2.130646 |

New-off wall times differ from baseline by +0.14%, -0.68%, and -0.94%,
respectively. These are **single balanced trios, not replicated estimates**:
they show no conspicuous observed overhead here, but establish neither
statistical performance equivalence nor a speed/RAM improvement. Focused
run times (39.038, 23.139, and 19.705 seconds in the same sample order) are
diagnostic-only observations, not unchanged-panel performance comparisons.

## Historical identity localization

Five of seven historical lost sequence identities are retained exactly in
both their full and focused current runs: ten matching rows across fourteen
historical-identity/run comparisons. Each matching identity appears at one
threshold in each run. Coordinates below refer to the candidate orientation;
the archived records preserve complete identities and every occurrence.

| Historical candidate | Threshold | Marker occurrences by cause | Marker-covered union | Full / focused exact identity |
| --- | ---: | --- | --- | --- |
| Gryllus CO1_1, 373 bp | 18 | 22 SCC | `[310,349)` | Yes / yes |
| Gryllus ITS_2, 978 bp | 4 | 5 SCC | `[93,112)`, `[276,296)` | Yes / yes |
| Drosophila 12S, 475 bp | 2 | 403 SCC, 11 collision | `[24,458)` | Yes / yes |
| Drosophila 12S, 487 bp | — | Not retained | Unresolved | No / no |
| Drosophila 12S, 499 bp | — | Not retained | Unresolved | No / no |
| Drosophila 16S_2, 590 bp | 2 | 459 SCC, 19 collision | `[94,570)` | Yes / yes |
| Heliconius ND1, 263 bp | 7 | 128 SCC | `[0,145)` | Yes / yes |

None of the five matching records carries an omitted-self-loop marker.
The two unobserved 12S identities are **not demonstrated absent from the
graph, enumerated candidate set, or biological sample**: both retention
truncation and search limits apply. The focused collection does not remove
that uncertainty.

### What this adds to the previous read audit

- **Gryllus ITS is now a concrete conservatism test case.** Its exact 978 bp
  current path visits each traversed graph node only once and is withheld
  solely for five pre-pruning SCC-node occurrences, at starts 93, 94, 276,
  277, and 278. It has no self-loop or collision-edge cause. Thus actual node
  revisiting by this path is not the reason for its rejection. Membership
  in a cyclic component elsewhere can still leave alternatives unresolved;
  this does not establish that the whole ITS haplotype is correct. Combined
  with the prior substantial marginal read evidence, it prioritizes focused
  evaluation of whether this restriction is too conservative.
- **CO1 localization agrees with the prior structural question.** The actual
  withheld 373 bp path's SCC-footprint union equals the previous complete
  alignment-placement union `[310,349)`. The earlier read assay favors the
  retained shorter local structure (27 exact bridge records versus zero for
  the longer structure). This provides no positive basis for rescuing the
  longer product; it is not proof of the absence of a rare longer allele.
- **ND1 is localized but not resolved.** The marker footprint reaches the
  candidate's left boundary and extends through base 145. Prior marginal
  support cannot phase this broad region, and the footprint is not itself
  a complete flanked SCC alternative. No justified rescue follows.
- **Drosophila remains weakly corroborated.** Two of its four lost identities
  are now linked to actual heavily SCC/collision-marked paths. The prior
  absence of qualifying 61 bp marginal windows remains a lack of evidence,
  not proof of biological absence. The other two identities remain unresolved.

## Retention and search limits

The pinned full panel contains **21 targets**, giving 12–13 paths and
1,560–1,561 markers per gene before threshold/gene caps. An incorrect
43-target/5–6-path statement in the preregistered prose is corrected in the
maintained protocol. The original snapshot is preserved. All code used the
actual target count; panel bytes, run settings, caps, and allocation rules
were frozen and unchanged throughout execution.

| Run | Observed withheld completions | Retained records | Retained bases | Retained marker occurrences |
| --- | ---: | ---: | ---: | ---: |
| Gryllus, full | 528 | 18 | 15,215 | 2,489 |
| Drosophila, full | 27 | 4 | 2,049 | 1,752 |
| Heliconius, full | 14,450 | 13 | 2,113 | 615 |
| Gryllus, focused | 182 | 34 | 32,019 | 183 |
| Drosophila, focused | 27 | 10 | 4,992 | 4,327 |
| Heliconius, focused | 1 | 1 | 263 | 128 |

These are search completions/retained diagnostic records, **not unique
biological molecules, lost products, or a sensitivity denominator**.
Full and focused totals cover different target sets. Focused Drosophila
exceeds 4,096 total markers legitimately: that is the per-threshold cap,
whereas the run cap is 32,768 and the two targets have separate thresholds.

Gryllus ITS retains 12/180 completions in the full panel and 32/180 in its
focused run, with explicit path-cap omissions. Drosophila 12S retains 3/26
and 9/26, respectively, with marker-cap omissions; its full-panel marker
share and focused per-threshold marker cap bind before the path cap.
The exact historical 475 bp sequence is retained; increasing the allocation
does not capture the historical 487/499 bp identities in this bounded study.

ITS hits the existing maximum-length boundary and, at a later threshold,
the node budget. Drosophila 12S/16S_2 hit
node-budget, DFS-state, and maximum-length limits; ND1 hits DFS-state and
maximum-length limits. Node-visit skips remain recorded separately. None
of these observations proves a valid alternative exists beyond a bound.

## Validation and evidence

Sol implements the Rust slice; Terra prepares and runs the comparison;
Astra independently reviews implementation, pre-run controls, and results.
Validation passes for both hash backends: **234 unit and 23 integration
tests each**, Clippy with warnings denied, and formatting. All 112 repository
Python regressions pass. The archived harness has 9 tests and the analyzer
11, including a full synthetic twelve-job analysis and adversarial identity,
accounting, matrix, source-bookend, and parity checks.

The frozen analyzer revalidates raw output manifests, exact commands,
source bindings, sequence/marker hashes and coordinates, complete marker
unions, exclusive omission accounting, and all retention budgets. It
recomputes full/focused parity instead of trusting receipt success flags.
A separate Astra checker, without importing production or analyzer code,
verifies all 80 retained records and 9,494 marker occurrences, raw outputs,
six parity comparisons, and frozen source bindings. It checks recorded input
identities rather than independently rescanning the multi-gigabyte inputs.
The real-binary smoke validates plumbing, not biological localization;
its synthetic run has no emitted or retained paths.

See the [evidence archive](../benchmarks/benchmark_results/v3.2-withheld-diagnostics-20260913/README.md)
for the original protocol, all twelve outputs/commands/receipts, immutable
source/build snapshots, historical hypothesis lineage, analyzer/tests,
validation logs, and independent review. Failed pre-benchmark attempts
are preserved: a smoke fixture omitted required schema-2 metadata, and the first
real launcher spelling was rejected by argparse before starting a job.
Fresh successful smoke v2/v3 and the correctly spelled twelve-run command
do not overwrite those failures. No real-data rerun or post-outcome cap
tuning is used for the reported comparison.

## Next decision

The diagnostic slice is complete, but **#148's read-supported reconstruction,
#153's biological disposition, and #129's release gates remain open**.

Next prioritize the two localized ITS regions: reconstruct the relevant
SCC entry/exit context and all bounded competing repeat/sequence alternatives,
then preregister complete distinguishing intervals and variation-aware
read-evidence assays. Use exact primary evidence and explicitly separate
substitution/indel sensitivity without collapsing the copy-number question.
Marginal windows must not be combined into a full-length phase claim.

Only positive evidence covering the actual restrictions proposed for removal,
appropriate true/collapsed-repeat and mixture controls, and a fresh high-copy
comparison can justify a targeted policy change. Do not merely admit every
non-revisiting path through a cyclic SCC, increase caps until desired
historical strings appear, or promote diagnostic candidates to references.
No release, tag, master merge, or gate exception is made here.

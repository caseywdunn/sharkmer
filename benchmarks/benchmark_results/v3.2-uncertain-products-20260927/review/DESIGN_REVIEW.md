# Independent uncertainty-output review

Reviewer: `threading_review_requirements`. Date: 2026-09-27.

The proposed opt-in feature is suitable for exposing bounded graph-path
hypotheses. It does not establish repeat copy number, read support, or a true
haplotype. Standard output remains the result of the existing admission rules;
its records likewise are not certified biological truth.

## Required behavior

- Enable collection through the existing withheld-path collector and its
  unchanged threshold, gene, and statically allocated run budgets. Only actual
  completed in-range paths qualify. Do not extend DFS, node visit limits,
  accepted quotas, or threshold stopping to obtain uncertain output.
- Finalize uncertain hypotheses after standard deduplication. Suppress only
  exact sequence or reverse-complement matches to final standard records.
  Merge exact oriented uncertain sequences while preserving every bounded
  source record. Preserve reverse-complement alternatives as separate oriented
  records so source coordinates require no transformation.
- Verify both accounting identities: retained path records equal excluded
  standard equivalents plus exact duplicates plus unique oriented hypotheses;
  total source records equal retained records minus excluded equivalents.
- Keep standard status, count, lengths, and output paths unchanged. An
  uncertain-only target still has standard status `fail` and zero standard
  products. Disabled output must omit the new optional stats field.
- Every uncertain header explicitly marks assembly status, repeat uncertainty,
  unresolved biological repeat length, and unevaluated read support. Sequence
  length is the emitted hypothesis length. Cause/threshold/source metadata is
  bounded by the collector; headers do not imply exhaustive enumeration.
  Header search/retention flags summarize the gene's evaluated thresholds.
- Stats retain source coordinates, hashes, marker identities, collection caps,
  drop reasons, and search limits even when the collector emits zero records.
  Combining diagnostic and emission flags may copy the same bounded payload,
  but must not allocate a second collector budget.
- Publish uncertain FASTA through the same checksummed output transaction.
  On-to-off reruns remove owned unchanged uncertain artifacts; altered or
  unowned files are refused and preserved. Detect uncertain-versus-standard
  output basename collisions before expensive work.
- Validation/recovery/bootstrap must never count uncertain hypotheses as
  standard products. The record header guard must still work if a file is
  renamed or mixed with standard records. Manifest receipt sets must include
  declared uncertain outputs without adding them to standard recovery.

## Source inspection so far

`src/pcr/uncertain.rs` consumes retained diagnostic records after admission,
uses exact oriented keys, applies exact/RC final-standard suppression, and
checks both accounting identities. `src/pcr/mod.rs` enables the existing
collector for either flag and leaves the standard search loop unchanged.
Early no-search returns include an explicit empty assessment when enabled.
Publication in `src/stats.rs` keeps standard status/count separate.

One namespace issue was reported to the authors: custom gene names permit
dots, so uncertain output for gene `A` collides with standard output for gene
`A.uncertain`. A blanket legacy `.uncertain.fasta` suffix filter also changes
handling of historical standard gene names ending in `.uncertain`; manifested
standard paths and explicit header markers can distinguish these cases.
The implementation now checks collisions before input work when emission is
enabled, preserves default-off names, and filters legacy records by header.

The independent validator review found a raw-token versus shell-token parsing
mismatch for quoted uncertainty labels. The author aligned the guard with the
standard header parser and retained detection of every duplicate label. At the
author's request, this reviewer also implemented source-local checks for cause
versus identity kind, typed node-visit bounds, and the exact merged union of
marker-covered intervals in `scripts/sharkmer_validate/uncertain.py`.

The resulting independent suite passed 18 test methods with 69 table-driven
subcases. Its receipt binds the inspected source and test hashes before and
after execution. `tests/test_uncertain_output.py` is a regression copy for CI;
it is not a second independent control set. These checks bind FASTA and sidecar
records, accounting, headers, sequence hashes, and source-local geometry. They
do not reconstruct graphs, replay threshold counters, establish biological
correctness, or clear the pending biological comparison.

## Required evidence before biological comparison

Use independent controls for flag-off/on standard invariance, actual
uncertain-only completion, multiple thresholds, unchanged accepted quotas,
exact/RC exclusions, near-variant retention, orientation/marker coordinates,
empty/truncated collection, and both diagnostic flags together. Exercise
publication ownership/cleanup and parser/header isolation. Bind reviewed
source, binary, fixture, command, and input hashes before the comparison.

Compare the same archived consumed R1 prefixes with baseline `956e771`, new
feature off, and new feature on. Report standard identity changes separately
from uncertain hypotheses. For each of the seven historical lost products,
distinguish exact standard match, exact uncertain match, and non-observation
within the bounded search/collector. Reverse-complement equivalence may be a
separate match class; approximate matches must not silently become exact
restorations. Seven visible sequences would not demonstrate seven corrected
products or seven biologically valid haplotypes. Two historical Drosophila
paths were not observed in previous bounded diagnostics, so this design does
not promise all seven will appear.

This file records design requirements and preliminary inspection. It is not
an execution clearance or a claim that pending controls passed.

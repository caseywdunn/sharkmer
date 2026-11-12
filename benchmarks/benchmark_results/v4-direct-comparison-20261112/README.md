# Rejected direct-table candidates — 2026-11-12

The independent audit passed for all 100 v2 invocations: complete sequence
inventories, warning causes, ingestion totals, all-entry count fingerprints,
completion manifests and metrics agree with the frozen protocol. The ten
inputs yield 287 sample-associated records, including warning hypotheses.
These are engineering equivalence checks, not additional biological validation.

**The candidate fails promotion.** The geometric mean wall ratio is
1.2846253579 (about 28.5% slower than v3.2.1), every sample violates its
unchanged wall allowance and two samples also violate their RSS allowance.
The computation completed; its result does not authorize default selection.
The custom linear-probing table is rejected as the production design.

| Input | Candidate/baseline median wall | Candidate/baseline median RSS |
| --- | ---: | ---: |
| SRR9278435 | 1.210 | 0.706 |
| SRR25099394 | 1.145 | 0.706 |
| SRR8617500 | 1.188 | 0.706 |
| SRR17535371 | 1.428 | 0.773 |
| SRR22396603 | 1.197 | 0.706 |
| SRR31887760 | 1.353 | 1.411 |
| SRR1057608 | 1.548 | 1.411 |
| SRR27962769 | 1.161 | 0.706 |
| SRR25378184 | 1.191 | 0.706 |
| ERR14009273 | 1.501 | 0.707 |

The v1 run stopped after 50 recorded invocations when the harness omitted the
historical `uncertain`→`warning` normalization used by the frozen baseline
producer. The sixth input's sequence identities/hashes match; it was a label
comparison bug. All attempted raw artifacts are retained. The first five
completed sample medians already reject v1 on runtime. V2 fixes that comparison
and avoids repeated canonicalization of validated memory-store keys, with
constant-work typed reverse complementation. Full v2 correctness passes, but
its measured runtime/RSS still rejects promotion. Neither attempt is discarded
or described as a successful speedup.

The baseline executable is the independently rebuilt v3.2.1 binary from the
frozen baseline bundle. Fresh alternating pairs use the same two-CPU allocation,
inputs, panels and unmodified numerical criteria. Count fingerprints run in
separate untimed invocations. Every selected workload has fewer than 2^32 total
windows, so its legacy per-key counts cannot have saturated; comparing its full
key/count fingerprint with the exact store is meaningful.

The independent verifier uses a separate line-oriented FASTA parser, reads raw
GNU time/stats/manifests, reconciles all 100 scheduled invocations, rechecks all
full input/panel hashes and binary fingerprints, and recomputes the medians and
promotion decision. `verification-v2.json` distinguishes a passed audit from a
failed candidate. Slurm 33016469 is the completed v2 measurement, 33120716 the
independent audit and 33120728 the checked archival job. V1 is terminal failed
job 32947879. Full source bindings are verified against commits 27ea351 and
724db89 rather than against a subsequently edited worktree.

Raw tar archives retain every file, with member hashes and sizes in the
artifact manifests. The archiver re-opened both archives and checked every
member. `archives.json` records byte counts and archive checksums. The exact
harness versions and launch/audit scripts are retained. Slurm batch RSS includes
helpers and verification; use per-invocation GNU time RSS for comparisons.
The next table design must pass the same gates; whole-run resource ownership,
external storage and replay remain open.

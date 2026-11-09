# v3.2 focused final check: results and disposition

2026-11-08. **Accept the measured costs for the scoped v3.2 release.** All
32 invocations passed the frozen sequence and output checks. The remaining Human
runtime flag is acceptable at a 0.04-second median increase; the other three
median increases are below 1.5%, and process memory is essentially unchanged.
This recommendation concerns the tested candidate and scope; final release
metadata and packaging checks remain separate.

The [summary](results-v1/summary.json), [complete result receipt](results-v1/result.json),
and [independent review](results-v1/independent_review.json) retain the findings.
The [raw archive](results-v1.tar.gz) preserves every execution receipt, output,
and log; the [artifact manifest](artifact_manifest.json) records their hashes.
The [frozen protocol](protocol-v1.json) and [protocol review](protocol_review.json)
define the comparisons and provenance. The measured candidate is the previously
validated `3.2.0-dev` binary with SHA-256 beginning `7ef94fedafec`, compared with
the pristine v3.1 default-build binary beginning `e751e5eeca91`.

The matched historical-panel comparisons retain all **86 v3.1 sequence identities**;
22 now carry assembly warnings, and none are lost. The historical candidate
inventory exactly reproduces 326 records: 64 ordinary and 262 warning records.
The two separate corrected-plant replays produce the expected current inventory
of **287 records: 66 ordinary and 221 warning records**. These are exact matches
to the frozen expectations, not new evidence of biological correctness. The
corrected-panel runs use a changed assay and are not performance comparisons.

There were three alternating baseline/candidate pairs for each timing sample,
followed by six remaining historical inventory replays and two corrected-plant
replays. The following are medians computed from the raw GNU `time` receipts.
Runtime values and RSS values are retained at their recorded precision;
percentage changes are rounded.

| Sample | v3.1 wall (s) | Candidate wall (s) | Wall change | v3.1 peak RSS (KiB) | Candidate peak RSS (KiB) | RSS change (KiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Human, SRR17535371 | 0.34 | 0.38 | +0.04 s (+11.76%) | 22,224 | 22,256 | +32 |
| Rhopilema, SRR8617500 | 28.01 | 28.07 | +0.06 s (+0.21%) | 4,461,704 | 4,461,724 | +20 |
| Liriodendron, SRR25378184, historical panel | 34.68 | 35.11 | +0.43 s (+1.24%) | 4,461,208 | 4,461,728 | +520 |
| Heliconius, SRR1057608 | 27.61 | 28.01 | +0.40 s (+1.45%) | 2,232,960 | 2,233,492 | +532 |

All median RSS changes are below 0.15%. Human exceeds the predeclared 5% runtime
review trigger because its baseline is very short. Its baseline range is
0.32–0.42 seconds and candidate range 0.37–0.40 seconds; the three paired changes
are +0.05, −0.04, and +0.06 seconds. Its median CPU time changes from 0.40 to
0.41 seconds. Accepting this flag is an explicit release-cost judgment, not a
claim that the effect is zero or statistically insignificant. No other runtime
or RSS review trigger fired.

Raw execution-artifact hashes, all 32 timing receipts, and the four sets of
performance medians/ranges reconcile with the result summary. All 32 constraint
receipts confirm CPU affinity `[39, 48]` and a 40 GiB address-space limit. These
are distinct physical cores 7 and 16 on socket 1 of
`a1130u05n01.mghpcc.ycrc.yale.edu`; all runs shared this allocation and wrote to
the same NFSv3 mount at `/nfs/roberts/project`. Slurm job `32602794` completed
successfully in **14 minutes 49 seconds**, with two allocated CPUs and 16 GiB
requested memory. Its batch-step peak RSS was **7,431,480 KiB**; that scheduler
measurement covers the batch workload and is distinct from Sharkmer process
RSS in the table.

GNU `time` directly measured Sharkmer. Exact input hashing and warming occurred
outside that timing; no cold-cache condition was established. Candidate stage
timings also reconcile with the raw YAML. v3.1 does not expose comparable stage
timings, so these results do not locate the cause of the runtime difference or
establish a storage-related explanation. Before/after bound-source fingerprints
are identical.

Three pairs on one node, two cores, these bounded input prefixes, and this NFS
environment do not establish statistical equivalence or calibrate performance
on other storage, larger inputs, thread counts, or hardware. This run does not
reproduce the large historical timing flags, but does not explain those earlier
measurements. The held-out biological dataset and full nine-panel biological
campaign remain explicitly deferred under the
[approved v3.2 scope](../../../dev_docs/PLAN.md), not passed. They remain required
before future assembly tuning or general accuracy claims.

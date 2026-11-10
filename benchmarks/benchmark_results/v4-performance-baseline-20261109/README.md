# Frozen v4 engineering calibration baseline

The protocol was frozen before measurement and before selecting a new default
counting backend. It rebuilds the immutable v3.2.1 source, uses ten retained
calibration inputs with their verified full-prefix hashes, and uses the released
panels including the corrected plant assay. These are existing calibration
samples; this work does not add biological truth or discharge #129.

`protocol.json` and `protocol.sha256` bind input/assay identities, expected
sequence inventories and ingestion counts, commands, repetitions, and acceptance
criteria. Every sample gets one untimed warmup and three timed runs at k19,
two physical CPUs, and `--chunks 0`. Each invocation must match the accepted
sequence inventory and counts and have a valid completion manifest. Human's
short source is explicitly preserved; the requested cap is not actual EOF.

Baseline job 32725595 runs on Bouchet's day partition with two CPUs, 16GiB,
and a two-hour job limit. A separate target directory isolates compilation.
The rebuilt v3.2.1 binary SHA is
`0d57521d5c6060df3e954d366b72fdccb0b3e40f9532162b400d62b14a59109c`,
matching the retained v3.2.1 engineering-check executable. Runtime outputs and
partial/terminal result.json are under tmp/v4-performance-baseline-20261109/run-v1.
The result is not accepted until all runs and artifact/source checks complete.

## Frozen promotion criteria

Default in-memory backend selection requires a paired baseline/candidate run
on the same allocation and inputs, with three alternating pairs initially.
Resolve noise with up to seven pairs when necessary; retain every attempt.
For the normal two-thread calibration matrix:

- Exact ingestion counts, target/sequence identities, warning membership, and
  completed-output receipts agree. Intentional header-name migration has its
  own schema/compatibility checks and cannot excuse different sequence output.
- The geometric mean of per-sample median candidate/baseline wall ratios is
  at most 0.95 (at least 5% aggregate improvement).
- No sample's median wall exceeds the larger of 1.05x baseline or baseline+0.05s.
- No sample's median process RSS exceeds the larger of 1.02x baseline or
  baseline+16MiB. Report allocator heap independently.
- Validate 1/2/4/8-thread agreement and both supported build configurations.
  Two-thread performance alone does not establish scaling or all-platform speed.

The external backend has separate requirements: exact count and full-sPCR
agreement, one pass over the original source, and measured RSS/temp-disk bounds
when the global store exceeds the configured RAM allowance. It need not outrun
an unlimited-memory backend. Graph, cache, queue, merge and replay reservations
all count toward the run limit; a bounded counter alone does not pass this gate.

GNU time reports full-process wall/user/system time and peak RSS. The runner
also retains application stage timings and exact product inventories. Cached
input and NFS output policy remain fixed within paired comparisons; hardware,
affinity, toolchain and filesystem evidence accompany each comparison. The
initial baseline receipts characterize the starting point, while default
promotion requires fresh paired measurements rather than a cross-job timing claim.

## Accepted baseline results

All 40 invocations passed. A separate verification job 32734740 reconciled
all 40 raw outputs with the frozen cases, rechecked input and panel hashes,
verified 37 source blobs against v3.2.1, and reproduced all completion receipts,
GNU time metrics and sequence inventories. Across the ten inputs there are 287
distinct sample-associated output records, including warning-bearing hypotheses.
The raw bundle contains 541 files/15,102,497 bytes, archived in
[raw-results-v1.tar.gz](raw-results-v1.tar.gz) (1,935,043 bytes). Its member hashes
are in [artifact-manifest.json](artifact-manifest.json); the archive checksum and
verification scope are in [verification.json](verification.json).

| Input | Median wall (s) | Median process RSS (MiB) | Output records |
| --- | ---: | ---: | ---: |
| SRR9278435 | 40.98 | 4357.29 | 7 |
| SRR25099394 | 36.45 | 4357.29 | 5 |
| SRR8617500 | 42.86 | 4357.29 | 5 |
| SRR17535371 | 0.63 | 22.09 | 5 |
| SRR22396603 | 33.41 | 4357.12 | 4 |
| SRR31887760 | 60.58 | 2181.54 | 36 |
| SRR1057608 | 42.26 | 2181.54 | 51 |
| SRR27962769 | 97.05 | 4357.56 | 109 |
| SRR25378184 | 52.93 | 4357.14 | 57 |
| ERR14009273 | 25.14 | 2181.11 | 8 |

The Slurm batch peak includes compilation and helper processes; use GNU time
per-invocation RSS for application comparisons. These measurements establish
the starting point only. No candidate has passed the promotion gates yet.

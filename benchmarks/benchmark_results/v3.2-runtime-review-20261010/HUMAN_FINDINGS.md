# Human output-location diagnostic

Slurm job **29012047** completed successfully in 77 seconds on
`c1106u05n03`, using two distinct physical CPUs `[32,36]`. All twelve normal
measurements and four separate traced diagnostics passed. Every invocation
consumed the frozen 108,518-read input and reproduced all five full product
identities and all nine target outcomes. All outputs and traces are preserved
under `results-v1/`, including the originally node-local files.

## Normal measurements

These measurements exclude tracing. GNU time directly wraps Sharkmer; hashing,
input validation, output inspection and artifact copying occur outside the
timer. Both locations use the same NFS input and frozen Human panel, with a
complete pre-read before each invocation. Product outputs vary between NFS and
node-local ext4; logs and time receipts are always node-local.

| Output location | Binary | Wall times, seconds | Median wall | Median CPU |
|---|---|---:|---:|---:|
| NFS | v3.1 release | 0.52, 1.06, 0.71 | 0.71 | 0.62 |
| NFS | Frozen retention | 0.95, 1.14, 1.90 | 1.14 | 0.69 |
| Node-local ext4 | v3.1 release | 0.48, 0.47, 0.47 | 0.47 | 0.62 |
| Node-local ext4 | Frozen retention | 0.51, 0.51, 0.50 | 0.51 | 0.67 |

The release-to-retention median wall difference is **+0.43 seconds (+60.6%) on
NFS**, compared with **+0.04 seconds (+8.5%) on node-local ext4**. The paired
within-block differences are NFS +0.43, +0.08 and +1.19 seconds, and local +0.03,
+0.04 and +0.03 seconds. Percentages are large partly because Human is a short
workload; absolute differences and the substantial NFS variation should guide
interpretation.

The same five sequences appear in every cell. The measured location effect
therefore does not require additional product records. Both binary lineages are
slower with NFS output in this batch, with a larger difference for the retained
candidate. This supports investigating output publication and filesystem
interaction. It does not isolate every contribution of the source changes.

The candidate's median ingestion time was 0.3121 seconds on NFS versus 0.2937
seconds with local output. Its median assembly-plus-FASTA-staging time was
0.2376 versus 0.1933 seconds. Median wall time after the pipeline timer was
0.3954 versus 0.0041 seconds. These independently calculated medians are not
additive. One local residual is −0.000051 seconds because GNU wall time is
rounded to hundredths of a second while stage timers have finer precision.

## Separate syscall diagnostics

Four extra invocations trace publication and metadata calls, one per binary and
location. Their wall times are **excluded** from normal summaries; the tracer
adds substantial overhead, especially on NFS with decoded file descriptors.
The following values are summed durations of calls whose decoded paths contain
the exact product output directory. They describe these traced runs only.

| Candidate output call | NFS calls | NFS summed seconds | Local calls | Local summed seconds |
|---|---:|---:|---:|---:|
| `openat` | 22 | 0.536715 | 22 | 0.000332 |
| `renameat` | 4 | 0.185964 | 4 | 0.000115 |
| `unlink` | 6 | 0.140992 | 6 | 0.000049 |
| `fsync` | 15 | 0.000792 | 15 | 0.004357 |

The narrower hypothesis that **fsync itself dominated the current delay is not
supported by these traces**. NFS metadata operations show much larger measured
durations here. The released binary also incurred NFS output-open latency:
six `openat` calls total 0.241406 seconds, versus 0.000102 seconds locally.

These syscall sums are not elapsed blocked time or a decomposition of normal
runtime. Calls can overlap across threads, tracing perturbs execution, and the
declared filter does not observe every possible syscall. Calls without a
resolved output-directory path remain explicitly unassigned. No syscall line
recognized by the parser was left unparsed. Expected `mkdir` failures reflect
already-created directories; all CLI runs and output checks passed.

## Meaning and next step

This bounded experiment reproduces a current **output-location sensitivity**
with matched input, panel, product identities, CPU allocation and frozen
binaries. It gives concrete metadata/publication calls to inspect if an
optimization is requested. It provides no basis for disabling synchronization:
the measured NFS synchronization calls were short in this run.

It does not prove what caused the original September Human wait, establish
acceptability, or close the Liriodendron, Rhopilema and Heliconius findings.
Those remain documented in `SAVED_FINDINGS.md`. It compares the frozen v3.1 and
retention binaries, not the subsequently corrected plant/dependency build.
Collaborators can now assess the absolute costs and intended deployment
filesystem before choosing an optimization or final-build performance check.
No additional assembly experiment was launched after this diagnostic.

## Reproducible evidence

- `protocol-v1.json` and `independent_protocol_review.json`: fixed sources,
  workload, schedule, constraints and approved interpretation boundaries.
- `results-v1/`: complete per-invocation commands, process constraints,
  statistics, FASTAs, normal GNU time receipts and raw traces.
- `analysis-v1/runs.tsv`: every timing and available stage measurement, with an
  explicit `traced` field.
- `analysis-v1/normal_medians.tsv`: normal median comparisons.
- `analysis-v1/trace_summary.tsv` and `analysis.json`: syscall summaries and raw
  line references, full artifact hashes, mounts and interpretation limits.
- `slurm-29012047-accounting.json`: terminal job state and Slurm log digest.

`analyze_human_results.py` independently authenticates preserved artifacts and
recomputes the analysis. To reproduce without overwriting the preserved result,
pass a fresh `--output` directory.

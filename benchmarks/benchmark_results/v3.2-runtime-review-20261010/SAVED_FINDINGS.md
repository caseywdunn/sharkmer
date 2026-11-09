# Runtime investigation: saved evidence

This bundle investigates the four runtime flags raised in the original v3.1
versus inclusive comparison. It preserves those historical findings. It is not
a release decision and does not establish that an emitted sequence is valid.

## Recomputed evidence

`analyze_saved_runs.py` reads the raw run receipts, GNU time files, actual process
constraints, statistics and FASTAs for all three pairs of Human, Liriodendron,
Rhopilema and Heliconius in two existing batches: 48 invocations. It verifies
their saved fingerprints, input/ingestion counts, stable full sequence identities
and frozen binaries. `saved_runs.tsv` contains each invocation;
`saved_summary.tsv` contains medians and all three values;
`paired_deltas.tsv` preserves the within-batch pairs.

`saved_analysis.json` records the sources and limitations. Two historical
context bindings now differ: the live uncertain-output validator and its test.
They supply no runtime values or compiled code to this analysis. Their old and
current digests are recorded explicitly. The saved timing harnesses, binary
snapshots, binaries and raw measurements remain checked against their receipts.

`summarize_saved_identities.py` additionally verifies equal archived prefix
digests across both batches and writes every retained, gained and lost literal
identity to `saved_identity_deltas.json`. These prefix bindings describe the
historical pre-run validation; this saved-data analysis does not reread all four
FASTQ files. The new Human diagnostic does reread its exact input.

### Original release → inclusive comparison

Medians in seconds; CPU is GNU user + system across all threads.

| Sample | Wall, release → inclusive | Absolute wall change | Wall change | CPU, release → inclusive | Products, release → inclusive |
|---|---:|---:|---:|---:|---:|
| Human | 1.55 → 4.29 | +2.74 | +176.77% | 0.65 → 0.71 | 5 → 5 |
| Liriodendron | 52.63 → 56.80 | +4.17 | +7.92% | 87.43 → 88.79 | 6 → 33 |
| Rhopilema | 43.87 → 46.53 | +2.66 | +6.06% | 64.36 → 64.81 | 5 → 5 |
| Heliconius | 41.59 → 44.07 | +2.48 | +5.96% | 71.66 → 72.28 | 7 → 19 |

All four exceeded the original 5% wall review threshold. Their CPU increases
were much smaller in absolute time than their wall increases. This difference
does not by itself measure I/O wait: CPU is summed over two worker threads.
Human consumed only 108,518 reads to EOF, so its short compute workload makes
fixed costs especially prominent in percentages.

Every released literal sequence remains present in these four cells. One
Heliconius sequence changed from standard to uncertain; the inclusive counts
are Human 5 standard, Liriodendron 6 standard + 27 uncertain, Rhopilema 5 standard,
and Heliconius 6 standard + 13 uncertain.

### Later inclusive → sequence-retention comparison

This is a separate batch on a different node, with the exact earlier inclusive
binary as its baseline. It is not a direct controlled release-to-current test.

| Sample | Wall, inclusive → retention | Absolute wall change | Wall change | CPU, inclusive → retention | Products, inclusive → retention |
|---|---:|---:|---:|---:|---:|
| Human | 0.66 → 0.61 | −0.05 | −7.58% | 0.70 → 0.70 | 5 → 5 |
| Liriodendron | 53.92 → 54.49 | +0.57 | +1.06% | 88.91 → 90.86 | 33 → 97 |
| Rhopilema | 43.05 → 43.71 | +0.66 | +1.53% | 63.30 → 64.41 | 5 → 5 |
| Heliconius | 42.42 → 42.88 | +0.46 | +1.08% | 73.23 → 74.10 | 19 → 51 |

All earlier inclusive identities remain present. The added records are 64 in
Liriodendron and 32 in Heliconius. Human and Rhopilema produce the same five full
sequences. The latest retained totals in these cells are 5, 97, 5 and 51.

The modest CPU increases in the three longer workloads are real measurements
in the saved pairs. They cannot all be explained by emitting more records,
because Rhopilema's record set is unchanged. They do not isolate graph search,
retention scoring, allocation, serialization or publication as a cause.

## What the stage measurements establish

The candidate's `pcr_s` wraps parallel assembly **and sequential FASTA staging**.
It is not a timer for assembly alone. `pipeline_total_s` is recorded before
statistics serialization, final output transaction commit/reporting and process
cleanup. `wall_minus_pipeline_seconds` is a residual, not an fsync measurement.

For later inclusive → retention, the median `pcr_s` changes are Human −0.0479 s,
Liriodendron +1.2651 s, Rhopilema +0.6613 s and Heliconius +0.4049 s. Other stages
also move: median read-ingestion time falls in all three longer workloads while
count-finalization time rises. Stage medians are independently calculated and
must not be added to reconstruct the median total.

The candidate statistics payload is also larger than the release payload. For
example, Liriodendron's median statistics file grows from 1,856 bytes in v3.1 to
1,859,068 bytes inclusive, then 1,964,150 bytes with retention; its FASTA payload
grows from 6,465 to 29,030 to 96,562 bytes. This is additional work, but the saved
data do not attribute the observed time to it.

The archived code shows candidate staging flushes/synchronizes files and final
publication synchronizes directories/manifests. The released implementation
writes FASTA files directly. A storage-dependent publication cost is therefore
a concrete hypothesis; no syscall trace was collected during the original runs.

## Same-binary bridge: why investigate Human first

The exact old inclusive binary, input-prefix digest, counts and products recur
in the second batch. Human's wall median changes from 4.29 s to 0.66 s while CPU
changes from 0.71 s to 0.70 s. Its median `pcr_s` changes from 1.2110 s to 0.2668 s,
and the median residual after pipeline timing changes from 1.9154 s to 0.0495 s.
The original Human inclusive values were 7.07, 3.75 and 4.29 seconds; the later
same-binary values were 0.66, 0.66 and 0.77 seconds.

This is strong evidence that the large wall effect is not stable across the
two execution environments. It does not identify the old cause. The original
batch used node `c1102u09n02` and affinity `[0,2]`; the later batch used
`c1102u05n03` and affinity `[20,24]`. Both selected two distinct physical cores.
Storage load, scheduler delays and syscall profiles were not recorded. The same
binary bridges for the other samples likewise change wall time; see
`same_binary_across_batches.tsv`. Do not multiply those ratios into a synthetic
release-to-current speedup.

## Cache and workload limits

Before every historical timed invocation, the harness validates and fingerprints
the entire staged input. That can warm input pages; there was no cache-drop or
verified cold-cache design. All timed inputs were uncompressed files and outputs
were beneath the shared repository. Human is a complete 108,518-record input;
the other cells use their frozen million-read prefixes. In particular, the
historical Liriodendron prefix combines 656,431 R1 and 343,569 R2 records. A new
R1-only prefix would change the workload.

Both comparisons use k=19, two threads, chunks=0, max-reads=1,000,000, frozen full
panels and no read threading. The frozen panel files must remain in the runtime
protocol while the live plant panel is corrected independently.

## Bounded next experiment

The reviewed `protocol-v1.json` and `run_human_diagnostic.sbatch` fix twelve
untraced Human measurements: two frozen binaries × NFS/node-local product output
× three blocks, with binary and location order alternated. Four separate traced
runs then inspect publication/metadata calls. Every product, log and trace is
preserved. The job requests two distinct physical CPUs, 16 GiB and 15 minutes;
each normal invocation has a 60-second timeout, each traced run 120 seconds, and
the child address-space limit is 12 GiB. Expected useful execution time is a few
minutes including validation and evidence copying. This small sample is far
below that memory limit in the saved receipts.

The diagnostic compares v3.1 to the **frozen retention executable**, which
predates the current dependency and plant-panel correction. It does not measure
the final corrected build. It tests whether output location changes current
wall/publication behavior; it cannot prove what caused the old waits, establish
biological validity, or close the three larger flags. There are no assembly
warm-ups; input, executable and panel pre-reads occur outside each timer. Logs
and GNU time files are node-local in both location cells.

After this diagnostic, further runs should address a named remaining question.
If collaborators need to localize the small retention CPU increase, Rhopilema is
the smallest useful next contrast: the old inclusive and retained builds have
the same five output identities, while the saved later pairs show +1.11 s CPU
and +0.6613 s in assembly-plus-staging. Three alternating pairs on one node with
one fixed output location would cost roughly 4.5 minutes of timed work. Profiling
or finer stage instrumentation would still be required to attribute a reproduced
increase to a particular function. This is a proposed follow-up, not an automatic
additional run or a claim that the original performance flags are resolved.

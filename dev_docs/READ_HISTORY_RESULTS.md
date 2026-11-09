# Ordered read histories: prototype findings

Completed 2026-09-27 on `experiment/its-read-threading`, starting from unchanged
`dev` baseline `956e771437a2410a7bcbb1504b033c61dd1e11b6`. The experiment
directory retains its 2026-09-26 start date. See the frozen
[design and evaluation protocol](READ_HISTORY_EXPERIMENT.md) and
[evidence bundle](../benchmarks/benchmark_results/v3.2-read-histories-20260926/README.md).

The standalone prototype preserves information that current aggregate
threading evidence discards: repeated edge occurrences, positioned branch
choices, exact run boundaries, qualities, and read identity. Synthetic controls
demonstrate the benefit with equal-sized read pools that have identical summed
per-read deduplicated edge and adjacent-pair counts but different repeat
traversals or distant linkage. This is an evidence-representation result;
production threading, path admission, and FASTA output are unchanged.

On the fixed Gryllus input, all four ITS classes remain represented in both
graphs, but none has a complete literal read witness. The three known spanning
records have exactly the same graph histories before and after pruning.
Their missing exact walks therefore were not created by the pruning step
between these snapshots. This does not determine why the missing k-mers were
absent from the earlier graph or establish biological absence.

## Frozen comparison and checks

The two scans use the original k=19, threshold-4 graph exports and the same
first 1,000,000 R1 records: 150,000,000 bases, 365,356,898 bytes, SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
No R2, extra coverage, new reference, or changed class window was introduced.
The prior approximate alignment/decoy assay was not repeated here.

The source/input/command freeze SHA-256 is
`337e27452472e6bbe99371b0ee12e3d41709a0079c1677ead72ff249180c45e9`.
Scanner SHA-256 is
`feb61bdb227a14fa51ff5a3b5c2faa789ae4ff92e15a4eedb9f5413befee0c27`.
All 25 implementation tests pass. Pre-execution independent review passed
861 checks, including 48 actual-class synthetic 150 bp placements across
both graphs/strands and 512 exhaustive ID/footprint matching cases. Twenty
execution-wrapper tests pass. Controls cover repeats, unequal mixtures,
unresolved distant phase, ambiguity before quality filtering, duplicate/RC
records, absent graph sequence, quality boundaries, and explicit failures.

Both biological commands completed with exit 0 and matching source/input
identities under Slurm job `27634080` (one CPU, 40 GiB allocation, 1,800-second
per-command limit). All prototype caps remained satisfied. Frozen sources,
raw retained records, full histories, timing, receipts, and source snapshots
are archived. Post-run auditing independently reconstructs the full-prefix
selection, raw-record bindings, runs and summaries; its checker also passes
two synthetic graph executions, 200 matching cases, and nine deliberate
corruption rejections. The final pass receipt is in `results_review/`.

The first real audit attempt (job `27649979`) failed while constructing its
final artifact list because the command supplied a relative freeze path to a
checker expecting an absolute path. Its failure receipt and log are retained.
The unchanged checker was rerun with an absolute freeze path under job
`27650157`, using the fresh output `provenance_ledger_review.retry1.json`.
This reporting retry changes neither the biological results nor the audit
criteria; reproduction commands must pass an absolute freeze path.

## Graph evidence

| Measure | Before pruning | After pruning |
| --- | ---: | ---: |
| Graph nodes / edges | 3,112 / 3,058 | 1,313 / 1,323 |
| Records with at least one exact graph edge | 773 | 284 |
| Unfiltered runs, both orientations | 917 | 331 |
| Unfiltered edge occurrences | 45,734 | 28,233 |
| Unfiltered cyclic-node visits | 2,234 | 2,234 |
| Repeated edge occurrences within unfiltered runs | 14 | 14 |
| Q20 runs, both orientations | 1,058 | 465 |
| Q20 edge occurrences | 33,874 | 20,789 |
| Q20 cyclic-node visits | 1,669 | 1,669 |
| Repeated edge occurrences within Q20 runs | 6 | 6 |
| Scanner wall time | 81.81 s | 83.80 s |
| Scanner peak RSS | 40,908 KiB | 30,908 KiB |

Run counts can increase after quality filtering because a run splits. Counts
sum reported orientations and events; they are not independent molecules or
copy-number estimates. The 773/284 selection requires any graph edge and is
broader than the earlier flank-based recruitment of 27 ITS1/28 ITS2 records.
These counts must not be interpreted as a direct change in that assay's
support. Timings are single prototype scans, not assembly performance or a
released-version benchmark; Slurm allocation/job RSS is a separate metric.

| Frozen local class | Window length | Exact records, pre/post | Q20 exact records, pre/post |
| --- | ---: | ---: | ---: |
| ITS1 historical | 61 bp | 0 / 0 | 0 / 0 |
| ITS1 extra cycle | 124 bp | 0 / 0 | 0 / 0 |
| ITS2 historical | 62 bp | 0 / 0 | 0 / 0 |
| ITS2 extra cycle | 147 bp | 0 / 0 | 0 / 0 |

All classes are graph-representable with noncyclic outer anchor nodes in
their appropriate orientation: ITS1 reverse, ITS2 forward. None passes the
preregistered local eligibility screen. The four Q20 zeros reproduce the
archived exact SCC assay. This is assay compatibility, not proof of absence,
independent truth, a full-product phase result, or a reason to relax gates.

## Known spanning records

Coordinates below are zero-based, half-open base intervals on the specified
oriented read. Break intervals instead index failing k-mer **starts**; they
must not be read as the length or location of a biological deletion.
The complete coordinate projections are saved in
`results_review/known_spanning_reads.json`.

- **ITS1 `SRR27962769.700435`, graph orientation `+`:** the unfiltered run is
  `[0,130)`, while the earlier approximate class alignment projects to
  `[9,133)`. Missing graph k-mers have starts `[112,132)`. Q20 splits the run
  into `[0,65)`, `[66,104)`, and `[105,124)`. The prior one-substitution,
  low-quality alignment therefore supplies no complete exact class witness.
- **ITS2 `SRR27962769.442978` and `.443018`, graph orientation `-`:** both
  have unfiltered runs `[0,73)` and `[87,150)`, with missing graph k-mer starts
  `[55,87)`. Neither run spans the prior gapped alignment's observed
  `[1,148)` interval. The two complete sequences and their unfiltered
  histories are identical.
- **ITS2 quality difference:** `.443018` retains both runs at Q20.
  `.442978` has Q11 at raw position 2 (reverse-oriented position 147),
  shortening only its second Q20 run to `[87,147)`. This quality boundary
  does not create another full-read footprint or an independent sequence
  pattern. The prototype does not repair this endpoint through a deletion.

Every one of these histories is identical between the two graph stages.
The earlier [gapped-alignment evidence](ITS_READ_RECRUITMENT_RESULTS.md)
remains valid within its stated scope: the shared ITS2 sequence pattern is
outside this exact graph-walk prototype's sensitivity. Zero complete exact
witnesses must not be presented as a new negative biological finding.

## Decision and next step

Keep the experiment isolated. Ordered histories retain useful evidence that
aggregate threading counts cannot encode, and the fixed-data run identifies
the exact continuity limits without changing reconstruction policy. There is
still no positive basis here for admitting the withheld Gryllus product.

The next production implementation step is bounded, quality-preserving replay
and positioned history collection under #145/#116, initially behind an opt-in
experimental path with feature-off count/product compatibility. Before a
reconstruction change, separately specify and test how variant-bearing reads
map to graph walks while retaining competing indel placements and unavailable
edges as uncertainty. More independent positions or linkage would require a
separately frozen evidence study. Do not fill gaps or stitch regions together
from marginal support.

The frozen protocol describes the later feature-off/on comparison against
`956e771` and the required released-v3.1.0 high-copy comparison. No new
production benchmark was run here; the historical 71-to-64 product loss and
#148/#153/#129 remain open. Held-out registration, all-panel biological
validation, dependency review, and user release review remain outstanding.
No merge, push, tag, release, or reference promotion is part of this result.

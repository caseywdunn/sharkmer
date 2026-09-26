# Gryllus ITS read-recruitment findings

Date: 2026-09-26. Branch: `experiment/its-read-threading`, based on
`956e771437a2410a7bcbb1504b033c61dd1e11b6`. Issues #148/#153/#129 remain open.
Production assembly, read threading, repeat gates, defaults, and references
are unchanged. This completes the standalone diagnostic in the
[frozen protocol](ITS_READ_RECRUITMENT_PROTOCOL.md).

## Findings

The variation-aware diagnostic recruits reads matching flanks of both ITS
cycle hypotheses, but
does not corroborate either historical short junction. ITS1 has one
low-quality read aligning to its extra-cycle alternative. ITS2 has two
records with Q20 local alignments to its extra-cycle alternative, but their
complete read sequences and placement are identical. Their optimal edit
tracebacks are ambiguous. This is one observed sequence pattern, not two
independent footprints or a reconstructed full-length haplotype.

Keep the repeat gate unchanged. These findings explain part of the earlier
zero-match result through quality and indel sensitivity; they do not show
that the historical 978 bp product should be restored. The released-v3.1.0
comparison remains 71 to 64 high-copy products because no assembly change
was made or newly benchmarked here.

| Region | Recruited records | Records with ordered flank pairs | Q20 unique local best alignments |
| --- | ---: | ---: | --- |
| ITS1 | 27 | 1 | 0 to either class |
| ITS2 | 28 | 2 | 2 to extra-cycle class; 0 to historical class |
| CO1 control | 62 | 36 | 29 to retained class; 0 to lost class |

There are 117 recruited FASTQ records in total. ITS1 has 11 left-only and
15 right-only records; ITS2 has 13 of each. Separate flank observations are
not spanning evidence. Quality, placement, and edit-type categories in the
machine report may overlap and must not be summed as disjoint populations.
"Unique" means best class/orientation/read interval within one declared
region, not a unique molecule, unique edit traceback, or genome-wide locus.

## The three ITS records spanning both flanks

All read coordinates below refer to the oriented read, zero-based half-open.
`I` consumes an observed base; `D` consumes a hypothesis base.

- ITS1 record `SRR27962769.700435` aligns to the 124 bp extra-cycle class
  at `[17,141)` on the reverse orientation, with one substitution
  (`2=1X121=`). Its minimum observed-span quality is Q11. There is no
  Q20-qualified alignment under the frozen limits.
- ITS2 record `SRR27962769.443018` aligns to the 147 bp extra-cycle class
  at `[1,148)` on the reverse orientation. A representative optimal edit
  script is `71=1D3=1I10=1X61=`: one deletion, one insertion, and one
  substitution, minimum Q25. Multiple optimal tracebacks exist.
- ITS2 record `SRR27962769.442978` has the identical 150 bp sequence and
  identical unfiltered 147 bp alignment, but that span includes a Q11 base.
  Its Q20 best match is `[1,147)`, a 146 bp observed span with an additional
  terminal deletion (`71=1D3=1I10=1X60=1D`), minimum Q25. Deleted hypothesis
  bases have no observed quality, so this is not evidence for a biological
  terminal deletion.

The two ITS2 complete read sequences share SHA-256
`5a3b30477ccc2c7e0bccd250ca150511e5e00a3e9ed400f6486b6f4d5f6a0e8e`.
Their two reported Q20 span lengths/offsets arise from the quality boundary;
they must not be presented as independently positioned biological support.
The differing instrument coordinates do not establish independent molecules.
No minimum-distance ITS contender matches the five full archived decoy
hypotheses within four edits. That finite-context result does not establish
genome-wide uniqueness or validate the extra repeat copy number.

## Controls, identity, and validation

The new diagnostic uses 21 bp flanks with at most two unit edits, and full
paired spans with at most four. Q20 applies to every observed base of a
qualified span. All passing placements, lower-quality matches, and raw
recruited records are retained. No bound was reached.

Separate commands using the unchanged CO1 scanner rerun the original exact and
one-substitution assays on the same prefix: retained CO1 reproduces
27 exact and 28 sensitivity records over 15 footprints, with zero for lost
CO1 in both modes. The new four-edit diagnostic has 29 retained-class
alignments: 27 exact, one substitution-only, and one indel-containing.
These nested sensitivities are not independent experiments.

The restored input is byte-identical to the original first million
SRR27962769 R1 records: 150,000,000 bases, 365,356,898 FASTQ bytes, SHA-256
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`.
Only this prefix was restored; the unused remainder of the older 8M-record
R1/R2 staged file was not required by the new protocol. Raw reads remain
outside Git. The initial sandbox DNS failure and successful ENA restoration
are both preserved as receipts.

The pre-outcome freeze SHA-256 is
`cb68a200f410a89980380cf7d4c290c3d825c7b501f43acf16736ec954abe600`.
It binds 36 source/review/input files and exact commands. Scanner SHA-256 is
`5c6268e840a626eb980713da8b0a1992030c50bad0e258241dcce703b3f0a75d`.
All three command receipts report success and matching before/after input
and source identities. Slurm job `27584622` used one CPU and a 40 GiB memory
allocation, with 1,800 seconds allowed per scan. The primary scan took about
194 seconds and 40,448 KiB peak RSS; these are diagnostic execution figures,
not a Sharkmer assembly performance comparison.

Pre-outcome validation includes 29 implementation tests, 1,800 independent
distance/traceback comparisons, 100 exhaustive-versus-seeded searches,
54 frozen-class variants, repeated-placement and context controls, and
13 independent execution-wrapper tests. All pass. The protocol's initial
256-pair resource cap was raised to 4,096 before biological scanning because
three synthetic repeated-class reads required 315 pairs; the evidence
preserves that rationale and the independent review.

See the [evidence bundle](../benchmarks/benchmark_results/v3.2-its-recruitment-20260925/README.md)
for complete ledgers, receipts, source snapshots, input verification, and
independent result review.

## Next decision

Keep development on the experiment branch. The next useful proposal should
specify how read-supported longer local alternatives would be reconstructed
and evaluated, with explicit uncertainty and repeat/mixture controls. More
independently positioned support and evidence linking the two regions are
needed before a full-product claim. Additional coverage or paired/longer
reads require a separately frozen new-input experiment; these calibration
results must not silently tune an admission exception.

Any production threading change remains a separate implementation and review,
followed by an identical-input comparison against the fixed `dev` baseline
and the required released-version high-copy regression comparison. No merge,
release, active-reference promotion, or closure of #148/#153/#129 follows
from this diagnostic. Held-out, all-panel, dependency, and user release
reviews remain open.

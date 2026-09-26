# Independent post-outcome review

Date: 2026-09-26. Scope: the completed frozen diagnostic and its two unchanged
CO1 compatibility scans. This review does not change the protocol, hypotheses,
source, thresholds, or biological input. Review files are separate from the
pre-outcome freeze.

## Validation and provenance

`local_alignment_review.json` records a passing independent unbanded
dynamic-programming review of all reported evidence: 1,910 flank matches,
6,058 paired-span/class checks, 5,916 representative edit scripts, 234
unfiltered/qualified best-contender sets, 700 finite-context minimum-distance
checks, and 2,572 reported decoy placements. Caching reduces duplicate numerical
work without merging records or physical placements. Reported coordinates,
observed spans, quality decisions, minimum distances, optimal-traceback tie
flags, complete pairing of reported flanks, and principal region counts all
agree. Independent semiglobal query-to-candidate DP checks finite-context
match absence/presence without using the production seeds or distance bands.

The separate `provenance_ledger_review.json` records verification of the freeze,
all frozen source files and snapshots, command receipts, output hashes,
input/source bookends, raw FASTQ record bindings, and unchanged CO1 controls.
Its accompanying verifier reads the declared first million records to check
ledger headers, sequences, qualities, and ordinals. This is evidence validation
on the same input, not a new discovery scan or independent biological replicate.

Freeze SHA-256:
`cb68a200f410a89980380cf7d4c290c3d825c7b501f43acf16736ec954abe600`.
Primary result SHA-256:
`d2d3792eed0574f2d3113dd7f23b8ea94e9e4a5be36a5491d3f02e3f253c9669`.

## Claim-level interpretation

| Region | Recruited records | Left-only / right-only | Paired records | Qualified local best |
| --- | ---: | --- | ---: | --- |
| ITS1 | 27 | 11 / 15 | 1 | None |
| ITS2 | 28 | 13 / 13 | 2 | Extra-cycle class in two records sharing one full-read sequence |
| CO1 | 62 | 16 / 10 | 36 | Retained class in 29 records; lost class in none |

The 117 recruited records are the sum across these regions in this run. Separate
left/right flank observations cannot establish a same-read junction, complete
product, or linkage between the two ITS regions. Flank similarity itself does
not prove a unique genomic locus. The absence of qualified historical-class
matches under these bounds is not proof that the historical sequence is absent
or that every historical product was erroneous.

ITS1 record `SRR27962769.700435` has a 124-base observed span matching the
extra-cycle class at edit distance one, with the unique representative script
`2=1X121=`. Four observed span positions are Q11; it has no Q20-qualified local
match. This is a low-quality sequence relationship, not qualified bridge support.

ITS2 records `SRR27962769.442978` and `SRR27962769.443018` have **identical complete
150-base sequences**, SHA-256
`5a3b30477ccc2c7e0bccd250ca150511e5e00a3e9ed400f6486b6f4d5f6a0e8e`.
They share the reverse-oriented 147-base interval `[1,148)` and unfiltered edit
distance three to the longer class. The representative script is
`71=1D3=1I10=1X61=`; optimal tracebacks are tied. Record 443018 has minimum Q25
over this span. The other record has Q11 at the final observed span base, so
its best Q20 interval is `[1,147)`: 146 observed bases, distance four, and an
additional terminal `D`. This terminal deletion is an alignment consequence
of excluding an observed low-quality endpoint, not evidence for a biological
deletion. Q20 does not assess any deleted hypothesis base.

The two reported qualified span lengths and offsets therefore do **not** supply
two independent footprints. There is one observed read-sequence pattern in two
records with different quality strings. Different instrument coordinates and
read IDs do not establish independent molecules. Both optimal traceback flags
must remain visible even though the per-region class/orientation/interval best
is unique for each qualified record.

No minimum-distance ITS contender matches any of the five full archived
candidates within four edits. Independent DP confirms that finite comparison
result. It cannot establish genome-wide uniqueness, repeat copy number, or
phase between ITS regions. The longer hypotheses also have different
read-length/callability constraints from the shorter hypotheses; negative
counts are not equally powered biological tests.

The original CO1 compatibility modes reproduce 27 exact and 28 one-substitution
retained records, with zero lost records. The broader diagnostic yields 29
qualified retained-class records: 27 exact, one substitution-only, and one
indel-containing representative alignment. These nested modes use the same
input and are not independent experiments. Descriptive edit-type categories
and quality categories may overlap; a representative gapped traceback is not
a biological variant call.

## Decision

The reviewed evidence supports reporting local similarity to longer ITS
hypotheses, with explicit quality and independence limits. It does not justify
restoring the historical short product, relaxing repeat gates, inferring an
extra repeat copy, reconstructing a phased full product, or incorporating a
production threading change. The 71-to-64 released-version discrepancy remains
unresolved by this diagnostic. The result can inform a separately specified
experimental reconstruction/threading proposal on the isolated branch.

Any stronger claim needs additional evidence, including independent placements
and linkage across the relevant regions, plus a separately frozen protocol for
new input and the required baseline/product/mixture/performance comparisons.
The draft `dev_docs/ITS_READ_RECRUITMENT_RESULTS.md` was reviewed and accurately
states these limitations; a minor opening-language refinement was suggested to
describe flank similarity rather than imply unique locus placement.

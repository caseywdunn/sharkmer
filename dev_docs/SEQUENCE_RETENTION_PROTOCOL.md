# Sequence-first uncertain retention assessment

Authorized 2026-10-01. Work on `experiment/sequence-first-retention`, from
`b5a6c7e`; keep the completed inclusive comparison and `dev` baseline fixed.
This evaluates delivery of completed repeat-uncertain sequences, not standard
admission or biological repeat correctness.

## Change under test

Keep standard graph construction, search, threshold selection, repeat
withholding, and accepted-product quota unchanged. The inclusive collector
retains bounded complete path sequences even when its detailed marker-entry
budget fills. Sequence path/base and detailed marker budgets are independent.
For each source path, preserve its literal sequence, hash, threshold, cause
counts, and node-visit summary. When detailed markers fit, preserve their full
identities and covered runs. Otherwise explicitly label the source and FASTA
record `marker_detail=summary_only`, omit detailed markers and covered runs,
and record the omitted detail count. Do not describe absent coordinates as
complete. Continue to label every such product `assembly_status=uncertain`,
`repeat_length=unresolved`, and `read_support=not_evaluated`.

The prior hidden diagnostic without inclusive emission retains its original
limits and format. Inclusive emission reserves the original 32-path threshold,
64-path gene and 256-path run budget for the exact full-detail paths that the
previous build would retain. An independent summary-only allowance adds up to
32 paths per threshold, 64 per gene, and a deterministic threefold per-gene
share of the former run path allocation (at most 768 across a run). The
displayed combined maxima are therefore 64/128/1,024 paths, with doubled
sequence-base bounds; the 4,096/8,192/32,768 detailed marker limits remain.
Extra sequences cannot consume the original reserve or displace old output.
Once any sequence cap binds, report truncation. The default output remains one
combined FASTA per target; `--no-emit-uncertain` remains standard-only.

Version the uncertainty assessment schema and make consumers accept the older
full-detail schema as well. Validate complete and summary sources, marker
counts, declared detail status, header labels, stats counts, and manifest
receipts. No historical target sequence or hash may appear in the production
retention rule.

## Verification

Run Rust formatting, Clippy and release tests on both hash backends, Python
consumer tests, and reference audit. Cover the detail boundary, path and base
limits, duplicate and standard-equivalent handling, source accounting,
malformed summary records, and opt-out parity. Keep old diagnostic-only
behavior under a direct regression control.

Replay SRR31887760 with the same first-million R1 bytes, full 21-target insect
panel, k=19, two CPUs and no read threading. Compare the candidate with
`ed04bde` using exact input and panel identities. Require the 487 bp and 499 bp
historical 12S sequences to appear as exact labelled uncertain products, while
all existing standard sequences and input counts agree. Account for every
additional 12S hypothesis and every detailed-provenance omission. Check
literal and reverse-complement identities; do not count length alone.

Then assess the original ten samples/five panels with three alternating old
inclusive/new candidate pairs, preserving the original input prefixes and
panels. Report record and unique-identity changes by class, old product
availability, per-sample runtime and peak RSS, all three repetitions, and
resource regressions. The old comparison's four runtime priority flags remain
visible context; this new comparison must not use a favorable subset to hide
them. Tests and replay use a frozen candidate source/build identity. Preserve
failed attempts, input/build/output receipts, the comparison protocol, and
independent review in a versioned evidence bundle.

The first candidate was rejected during the paired matrix: an early plant
path consumed a shared 32-path threshold allowance and displaced a previously
emitted 261 bp `psbA-trnH` identity. Preserve that failed run and require old
uncertain identities to remain a subset of new output in every cell. The
reserve design above is the corrected protocol for the next build.

Successful retention means completed bounded paths are delivered with honest
labels. It does not imply correct repeat counts, independently supported
haplotypes, or complete enumeration beyond existing graph and DFS limits.

# Final independent result and interpretation review

Date: 2026-09-27. **Pass for the standalone experiment and its stated limits.**
This review does not approve production integration, repeat admission, or release.

The accepted full-prefix audit is `provenance_ledger_review.retry1.json`, SHA-256
`058fa13963811e29d97b8427bc809734d15b2918b0965c5760cf9da8cefa77f0`.
It binds the original freeze
`337e27452472e6bbe99371b0ee12e3d41709a0079c1677ead72ff249180c45e9`
and unchanged checker
`f22bfb7ce204318d60bde8368a533228bd2165e2503fe942bb4861034c645a29`.
The first audit failed only when formatting a relative freeze path in its final
artifact list. Its failed receipt is retained. The unchanged checker succeeded
with an absolute path and fresh output; no producer, graph, class, threshold,
or audit criterion was changed. Reproduction instructions document this path
requirement rather than presenting the failed attempt as successful.

The accepted audit independently verifies all 1,000,000 prefix ordinals,
complete retained-ledger selection and raw-record binding, every reported
positioned history, quality boundary, branch/cycle event, class placement,
counter, witness screen, source snapshot, receipt, and output identity.
The reviewer also inspected its source and independently reran two synthetic
graph executions, 200 matching cases, and nine deliberate corruption checks.

The findings document, handoff, plan, release-review update, and evidence README
correctly distinguish an evidence representation prototype from production
threading and reconstruction. Their graph/run/event numbers match the actual
reports. Scanner timings of 81.81/83.80 seconds and peak RSS of
40,908/30,908 KiB match GNU time receipts and are identified as diagnostic
execution figures, not an assembly benchmark.

The 773 pre-pruning and 284 post-pruning retained records are graph-matching
records, not independent molecules or comparable counts from the earlier
27/28 flank-recruitment assay. All four local classes remain representable,
with ITS1 in reverse and ITS2 in forward graph orientation. Both graph stages
have zero complete literal class matches before and after quality filtering;
this reproduces the archived exact assay without adjudicating biological
absence or the earlier indel-aware alignments.

The three known records' complete histories are identical between graph stages.
ITS1's exact unfiltered interval `[0,130)` does not cover the earlier observed
alignment span projected to `[9,133)` in graph orientation. ITS2's two exact
runs `[0,73)` and `[87,150)` do not span its observed alignment interval.
The report correctly labels missing-kmer intervals as **kmer-start ranges**,
not locations or lengths of biological deletions. The duplicate ITS2 sequences
remain one sequence pattern; Q11 at raw position 2 changes one Q20 endpoint
without creating independent positioned support.

Synthetic controls establish a limited but useful representation benefit:
equal-sized read pools can have identical summed per-read deduplicated edge
and adjacent-branch counts while differing in repeat multiplicity or distant
phase. Ordered histories preserve that difference on the tested graph. This
does not establish identical global kmer abundance, correct reconstruction from
missing edges, genome-wide specificity, or full-product phase.

The proposed next step—bounded, quality-preserving replay and positioned
histories under #145/#116 with experimental opt-in and feature-off compatibility—
is consistent with the implementation gap. Variant-aware mapping still needs
an explicit ambiguity-preserving design. Additional coverage or linkage requires
a separate frozen evidence study. No current result supplies positive evidence
for removing repeat restrictions or restoring the withheld Gryllus product.
The 71-to-64 released-version comparison and existing release gates remain open.

`final_results_review.json` binds the reviewed document and evidence bytes.
Archive creation/read-back verification is the following packaging step; this
review does not attest to an archive that has not yet been published.

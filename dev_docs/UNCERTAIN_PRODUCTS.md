# Uncertain repeat-structure candidates

Issue [#157](https://github.com/caseywdunn/sharkmer/issues/157). Implemented for
evaluation on `experiment/uncertain-products`; integration and release are
separate decisions. This protocol was written before the new real-data runs.

## Using the output

Add `--emit-uncertain` to an ordinary sPCR invocation:

```sh
sharkmer --pcr-panel insecta --emit-uncertain --max-reads 1000000 \
  -s sample -o output reads.fastq.gz
```

Ordinary products remain in `{sample}_{gene}.fasta`. Additional completed
in-range paths rejected for unresolved repeat structure are saved separately
in `{sample}_{gene}.uncertain.fasta`, when any are retained. Their identifiers
use `{sample}_{gene}.uncertain_{index}`. Headers include:

```text
assembly_status=uncertain uncertainty=repeat_structure repeat_length=unresolved read_support=not_evaluated
```

`hypothesis_length` is the length of the emitted sequence, not an established
biological repeat length. The sequence can contain a collapsed repeat or an
unsupported repeat traversal. This option makes such sequences available for
inspection without declaring their structure resolved. It does not evaluate
read support, correct repeats, enumerate all possible repeat counts, or
establish that one molecule spans the candidate.

The option is off by default. It preserves ordinary gene status, product
counts, output, search order, thresholds, accepted-path quota, and graph
admission rules. A target with uncertain output alone retains standard status
`fail` and `n_products: 0`. The new optional `uncertain_products` mapping in
each gene's stats records a separate count, lengths, output file, candidates,
and source evidence. When enabled with no candidates, the mapping still
reports search and retention observations; no empty uncertain FASTA is made.

## Provenance and limits

Only actual completed, in-range paths retained by the existing bounded
withheld-path collector qualify. The flag never extends graph search or
reconstructs historical products by reference matching. It shares the
collector with `--diagnose-withheld-paths`; using both does not double the
budget. Hidden diagnostic sequences are included only when that hidden flag
is also requested.

| Scope | Path records | Sequence bases | Marker occurrences |
| --- | ---: | ---: | ---: |
| Threshold | 32 | 524,288 | 4,096 |
| Gene | 64 | 1,048,576 | 8,192 |
| Run | 256 | 4,194,304 | 32,768 |

Run capacity is divided deterministically among selected target indices.
Unused shares are not redistributed. The 21-target insect panel consequently
has only 12–13 retained path records per target before other limits apply.
These are collection budgets, not process-memory limits. Oversize sequences
or incomplete marker sets are omitted as whole records, with reasons counted.
The output is a bounded sample, and a missing sequence is not evidence of
biological absence.

After ordinary product deduplication, exact matches or reverse complements
of final ordinary products are excluded from uncertain output. Remaining
uncertain paths are merged only for identical oriented sequences, preserving
all bounded source records. Near variants and reverse-complement alternatives
remain separate hypotheses. Counts therefore do not estimate loci, molecules,
or haplotypes. Product order is deterministic by oriented sequence.

Each source identifies its coverage threshold, node-visit limit, maximum
observed visits, repeat-marker causes, and marker occurrences. Marker identity
sequences and SHA-256 hashes bind zero-based half-open coordinates to the
emitted orientation. The causes are pre-pruning omitted self-loop nodes,
pre-pruning cyclic-component nodes, and retained collision edges. Their
covered runs describe marker footprints, not necessarily whole ambiguous
regions or intervals ready for a read-bridge test.

The assessment preserves each threshold's search counters and retention
accounting. Header `search_limits_reached` summarizes node/DFS/accepted-path
quota or node-visit restrictions across all attempted thresholds for the
gene; `retention_truncated` summarizes collection omissions across the same
gene. The configured product-length range remains a separate search bound.
Neither flag describes exhaustive sequence recovery, nor implies that an
additional candidate must exist beyond a limit.

## Publication and consumers

Uncertain FASTA is staged, checksummed, and published in the same transaction
as ordinary FASTA and stats. Rerunning with the option disabled removes
previous unchanged owned uncertain files. Modified, unowned, or symlink
collisions remain protected by the ordinary publication contract. Custom
target pairs such as `A` and `A.uncertain` are rejected when this option is
enabled because their output filenames would collide.

The validation parser verifies uncertain metadata, source-local coordinates,
sequence identities, FASTA/sidecar agreement, and transaction receipts. It
does not replay the graph to validate biological structure or all search
counters. Standard validation and reference bootstrap exclude uncertain
records. The header guard also applies to renamed or mixed legacy files;
explicitly declaring an uncertain record as standard is an error. External
consumers should honor `assembly_status=uncertain` and preserve this metadata.

## Preregistered assessment

Compare baseline `956e771437a2410a7bcbb1504b033c61dd1e11b6`, the new binary with
the flag off, and the same new binary with the flag on. Use the frozen full
21-target insect panel and exact previously consumed first-million R1
prefixes of Gryllus SRR27962769, Drosophila SRR31887760, and Heliconius
SRR1057608. Prefix hashes, the seven historical lost sequences, and source
lineage are saved in the
[comparison bundle](../benchmarks/benchmark_results/v3.2-uncertain-products-20260927/).

Run nine sequential invocations at k=19, two CPUs, chunks=0, and without read
threading. Rotate order: Gryllus baseline/off/on; Drosophila on/baseline/off;
Heliconius off/on/baseline. Prewarm inputs, use a 40 GiB address-space cap and
1,800-second timeout, and record effective affinity and enforcement. Avoid
concurrent builds/tests during timed invocations. Keep fresh attempt
directories and receipts, including failures. Freeze source, harness,
executable, interpreter, panel, and input identities before execution, after
independent review and controls pass.

Require exact equality of read/base/k-mer count fields, ordered standard
FASTA records including headers, gene outcomes, and pre-existing threshold
diagnostics across all three modes. Remove only the new optional assessment
for standard parity. Validate uncertain publication and retained source
accounting separately. An unexplained parity mismatch blocks acceptance.

For each historical loss, report exact oriented matches in standard and
uncertain outputs separately. Reverse-complement matches are a separate
equivalence category; approximate matches do not count as exact restoration.
Non-observation remains conditional on the unchanged search and retention
limits. Do not tune caps or target membership after seeing recovery. The
prior diagnostic found five of the seven identities; this is a prediction
to check, not a new result.

Single trios per sample provide descriptive runtime/RSS observations only.
This is an availability and compatibility assessment on calibration samples,
not held-out biological validation, a new released-v3.1 recovery benchmark,
or evidence that re-emitted sequences have correct repeat lengths. The
read-supported reconstruction and biological-loss questions in #148/#153,
and the remaining release gates in #129, stay separate.

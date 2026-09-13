# Bounded withheld-path diagnostics

Date: 2026-09-13. Baseline: `58ca2a6` on `dev`. Issues #148, #153, and #129.
The user approved diagnostic implementation after the
[read-backed audit](READ_BACKED_AUDIT_RESULTS.md), **not** a gate relaxation or
release. This protocol is recorded before the new real-data outcomes.

## Scope and interpretation

The preceding audit could identify historical sequences and aggregate SCC
withholding counts, but could not link them to the exact paths withheld by
the current search. Add a hidden, opt-in `--diagnose-withheld-paths` mode that
records a bounded sample of actual completed, in-range, repeat-withheld paths
in the existing transaction-owned YAML stats. These sequences are diagnostic
hypotheses, never FASTA products, public references, or validated haplotypes.

For each retained path, preserve its exact oriented sequence and SHA-256 plus
every encountered marker occurrence and its sequence-relative coordinates.
Use zero-based, half-open base spans and stable oriented sequence identities,
not transient graph node/edge indexes. The enclosing gene/threshold identifies
the search context. Node-marker spans cover their `(k-1)`-mer; collision-edge
spans cover their k-mer. Repeated occurrences have separate coordinates.

Pre-pruning SCC membership, omitted-self-loop markers, and retained collision
edges explain **why the path was withheld**. Their footprint spans are not
automatically the complete ambiguous interval between distinguishing entry
and exit flanks. Do not feed them directly into a bridge verifier as though
all repeat-copy alternatives have been enumerated or localized. This change
does not add SCC reconstruction, read evidence, replay, or path admission.

## Resource and compatibility contract

- Default off: omit the new payload and perform no per-path diagnostic
  sequence reconstruction, hashing, or marker-list allocation.
- Keep existing search order, thresholds, extension/pruning, eligible quotas,
  withholding counters, FASTA generation/ranking, and counting unchanged.
- Retain whole sequences and complete marker-occurrence sets only. Test
  budgets before materialization; an oversize record is omitted with an
  explicit count, never published as a partial complete candidate.
- Independently cap diagnostic paths, sequence bases, and marker occurrences
  per threshold, per gene, and for the run. Assign run quotas from stable
  panel indices independently of parallel scheduling. Unused shares are not
  redistributed.
- Expose allocated/effective budgets, observed versus retained path counts,
  omissions and their reasons, and search/visit-limit observations. A zero
  allocation or saturated cap must not look like no withheld candidates.
- Diagnostic truncation never consumes an eligible path slot or terminates
  assembly search. Absence of a historical hash from a capped collection is
  not proof that its path was absent from the graph or never enumerated.

The reviewed design caps each threshold at 32 paths, 512 KiB of sequence,
and 4,096 marker occurrences; each gene at 64 paths, 1 MiB, and 8,192 markers;
the run at 256 paths, 4 MiB, and 32,768 markers. These are payload budgets,
not process-RSS ceilings. With 21 insect targets, static sharing permits only
12–13 retained paths per gene before other limits; report this sampling limit.
Focused target panels below increase the share without increasing any cap.

Post-run prose correction: the archived preregistration snapshot incorrectly
said 43 targets and 5–6 paths per gene. The pinned panel actually contains 21
targets. Implementation, harness, and analyzer use the actual target count;
no panel, cap, allocation algorithm, or invocation changed to make this
correction. The original snapshot remains preserved in the evidence archive.

## Preregistered verification

### Controls and independent review

Sol implements; Astra independently reviews the design, implementation,
coordinate/identity checks, bounds, and interpretation. Controls must cover:

1. Off/on modes preserve ordered eligible paths/products, count totals,
   threshold order, existing aggregate counters, and search-limit flags.
2. Positive withheld paths retain exact sequence/hash and all self-loop, SCC,
   and collision-edge occurrences, including repeated nodes/orientations.
3. Pruned graphs retain correct pre-pruning marker semantics without graph
   index leakage or marking an unrelated clean path.
4. Zero, tiny, saturated, and oversize retention budgets remain bounded and
   explicit; no partial sequences, hidden truncation, or eligibility changes.
5. Default serialization omits the payload; opt-in serialization is tied to
   the existing stats transaction and clearly labels unsupported hypotheses.

Run focused tests first, then the full Python/Rust suites, formatting, and
Clippy. Freeze reviewed sources, executables, harness, panels, inputs, and
interpreter before real-data execution. Do not edit those files during runs.

### Full-panel three-way comparison

Perform nine sequential invocations: one three-way comparison for each of
SRR27962769 (Gryllus), SRR31887760 (Drosophila), and SRR1057608 (Heliconius).

| Sample | First | Second | Third |
| --- | --- | --- | --- |
| Gryllus | Baseline off | New off | New on |
| Drosophila | New on | Baseline off | New off |
| Heliconius | New off | New on | Baseline off |

Use the same frozen current insect panel and the exact consumed first 1M R1
records from the preceding audit, k=19, two threads, chunks=0, and no read
threading. Preserve prefix prewarming, CPU affinity, 40 GiB address-space
limit, and 1,800-second timeout from the scoped comparison where available;
record actual enforcement rather than assuming it. No builds/tests run
concurrently with timed invocations. Capture source/binary/interpreter hashes,
commands, exit status, timings/RSS, full-file and consumed-prefix identity,
output manifests, and artifact hashes in fresh per-job directories/receipts.
No silent reuse or overwritten attempts.

Require equality of all five read/base/k-mer count fields, ordered emitted
FASTA records including headers, gene outcomes and threshold ordering, and
all pre-existing aggregate threshold diagnostics. Strip only the explicitly
new optional diagnostic payload for this comparison. Validate its sequence
hashes, marker identities/coordinates, budgets, and observed/omitted counts
separately. A mismatch fails the compatibility gate and requires investigation.

### Focused diagnostic localization

Preregister three additional diagnostic-on invocations using the same read
prefixes/settings and exact copies of selected current panel entries:

- Gryllus: CO1_1 and ITS_2.
- Drosophila: 12S and 16S_2.
- Heliconius: ND1.

Preserve panel identity metadata and every selected target parameter/name;
freeze the selected panels and their distinct hashes before observing results.
Removing other targets changes the diagnostic budget allocation, not the
selected gene's assembly settings. Compare each selected gene's products,
existing aggregate diagnostics, and input/count totals to its full-panel
diagnostic-on result. These focused runs are **not** unchanged-panel timing
comparisons. Do not tune caps or target membership after seeing support.

Compare the seven frozen historical lost-sequence hashes to retained actual
withheld sequences. An exact identity allows marker localization for that
current-run path; report every threshold where retained. No match under
retention or search bounds stays unresolved. A marked region still requires
independent graph/alternative analysis before becoming a complete read-bridge
question; no biological verdict or gate exception follows from its presence.

### Reporting and release boundary

These are twelve diagnostic/compatibility invocations, not a new released-v3.1
benchmark or held-out biological validation. Single three-way comparisons
give descriptive timing/RSS observations only, not statistical equivalence,
universal absence of regressions, or a demonstrated speed/RAM improvement.
Archive successful and failed attempts with exact lineage; stop on unexplained
parity failures. Document what was localized and what remains truncated.

Keep #148's read-supported reconstruction, #145/#116's bounded quality-aware
replay/batching, and #153's biological release disposition open. No merge to
master, release tag, publication, or gate relaxation is authorized here.

The completed comparison and its limitations are recorded in
[withheld-path diagnostic results](WITHHELD_PATH_DIAGNOSTICS_RESULTS.md).

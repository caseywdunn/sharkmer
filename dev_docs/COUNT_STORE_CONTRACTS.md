# v4 count-store, replay, and resource contracts

Implementation under [#139](https://github.com/caseywdunn/sharkmer/issues/139).
This document defines the interfaces required by #140–#145. The current branch
implements typed identities and positioned rolling extraction. An exact in-memory
store candidate and shared-ledger conformance API now pass their focused tests.
Assembly consumers now use fallible borrowed count access and retain u64 values;
normal exact ingestion is now available through a hidden development selector,
while default selection and the replay vector are still being migrated. The
[active plan](PLAN.md) records verification and completion.

## Identity and preprocessing

`kmer::identity` defines `KmerSize`, `OrientedKey`, `CanonicalKey`, and
`KmerWindows`. The supported representation has k in 1–31. Raw values with
bits outside their declared k are errors. Canonical count identity is
`min(forward, reverse_complement)`; graph and primer identities preserve
orientation. Keys crossing a storage boundary carry k, and a query with the
wrong k must fail rather than return absence or alias another sequence.
Internally a store may pack just the bits after validating this boundary.

Rolling windows retain original zero-based byte positions. Uppercase ACGT
contribute bases; N breaks a run. Every other byte is an error under the
released ingestion policy. Extraction does not silently uppercase, discard
quality, splice runs, or collapse gaps. The iterator emits one error and then
terminates; a streaming ingestion consumer must propagate that error and fail
publication. The legacy collecting wrapper returns an error without exposing
its partial vector. This maintains its existing count-before-insertion behavior.

## Exact store boundary

The new normal counter must implement the following semantics, independent of
whether counts live in memory or on disk:

| Operation | Contract |
| --- | --- |
| Add `(CanonicalKey, NonZeroU64)` to a builder | Checked addition; maximum exact count is `u64::MAX`. Overflow is an error, never saturation or wrapping. A failed update cannot replace a valid count with a partial sum. |
| Finalize builder | Consume or move its state into an immutable store. Flush and verify persisted state before publishing its completion marker. No mandatory full second in-memory table. |
| Lookup | `Result<Option<NonZeroU64>>`: `None` means the exact key is absent, a positive integer is present, and I/O/corruption/wrong-k errors remain errors. Zero is never stored as an observed count. |
| Bounded batch lookup | Same results and input order as individual queries, including duplicates. Enforce configured batch capacity; do not buffer all caller queries. |
| Visit entries | Yield every exact key/count once, with bounded buffering and propagated visitor/I/O errors. Hash-map iteration order is unspecified; callers requiring stable order must arrange it explicitly. |
| Filtered view | Check thresholds on exact lookup/iteration results without copying the entire store or changing stored counts. |

Fingerprints and minimal-perfect hashes cannot establish membership by
themselves. Any such index must verify the complete queried key. An absent
key must never acquire another key's count. Finalization records distinct-key
count and checked occurrence totals; totals that exceed the supported integer
range must fail explicitly as well.

The read-only interface must be safe to share across target workers. A disk
implementation accounts for its file descriptors, lookup buffers, and cache;
it cannot materialize the full count table as an implicit assembly prerequisite.
Legacy saturating `u32` tables do not automatically satisfy this exact interface.

## Replay and input identity

Replay is optional for counting-only operation. When requested, it must be
produced during the original ingestion pass, including for stdin and remote
streams, without a second original-source download. Its logical records retain:

- Global selected-record index, source identity, original read identifier,
  fragment/pair identity, and mate designation.
- Original sequence bytes, complete quality bytes and encoding, and original
  orientation. Quality length must equal sequence length.
- Original positions and invalid-base discontinuities: store complete reads,
  rather than concatenating valid runs. Any transformed graph mapping must
  retain its relation to these original coordinates.

Selection metadata records the effective read cap, paired/interleaved policy,
source order, and a deterministic digest of the selected records. Preserve
released mate-selection behavior explicitly, including any paired cap rounding;
test that count and replay consume exactly the same selection. A requested cap
is distinct from actual records available at EOF.

Define a versioned, length-framed replay format with a fixed byte order,
per-record validation, and a completed-stream checksum/count footer. Validate
record and frame lengths against the configured byte budget before allocation.
Truncated, oversized, checksum-invalid, or incompatible frames are errors.
The same rules apply to uncompressed and compressed source streams, including
concatenated gzip members. Opening replay must not imply loading all records.

## Resource ownership and persisted identity

Reserve independent byte budgets for mutable counting, ingestion/worker queues,
external partition/sort buffers, immutable-store lookup caches, active target
graphs, evidence batches, replay buffers, and temporary disk. Their combined
reservations plus accounted process overhead must fit the configured run budget.
Graph concurrency is bounded by its own reservation and is not simply the
number of available threads. Temporary storage has a separate limit and path.

Use a 16 GB laptop as the initial hardware class and reserve OS headroom.
Concrete defaults, overhead factors, and quantitative performance gates must
be measured and frozen before backend selection. Record both allocator heap
and process RSS; a logical queue capacity alone does not prove a process bound.
Resizing must account for simultaneous old/new storage. Strongly skewed
partitions and long reads must fit or fail visibly within the same bounds.

Persist format version, k, canonicalization, count width/overflow policy,
preprocessing, effective input selection and digest, backend parameters, and
completion/integrity state. Include replay identity when evidence uses it.
Reject incompatible stores before lookup. Write new state to owned staging
paths; errors, interruption, and disk-full conditions preserve existing complete
stores and original input. Incomplete work remains identifiable as incomplete;
it cannot replace a completed count store or output manifest. Resume semantics
must be specified and tested before claiming restart support.

## Assembly/evidence boundary and legacy adapter

Assembly consumes immutable exact membership/count lookup plus optional bounded
replay/evidence access. Storage changes preserve standard assembly, inclusive
warning delivery, sequence identities, and existing admission limits. No-data,
failed I/O, search truncation, unresolved repeats/phase, and partial sequence are
separate outcomes. Observed k-mers and local read support do not certify a whole
haplotype. New evidence-driven admission belongs to the later work in PLAN.

Dispatch explicit `--chunks > 0` into the retained incremental implementation.
Its snapshot format, viewer, and combined histogram/PCR behavior remain covered
by compatibility tests. The exact-counter interface does not require histogram
snapshots. Label legacy count precision honestly and reject unsupported backend
combinations; never cast a possibly saturated legacy table into an allegedly
exact store. A final count histogram is a separate operation from incremental
snapshots.

## Verification before closing #139

- Typed identity and independent scalar window-oracle tests cover k boundaries,
  strand normalization, wrong-width keys, N gaps, positions, and invalid bytes.
- Actual store/builder implementations test absent-key non-aliasing, exact
  overflow, finalization, bounded batch order, and failure/corruption handling.
- Replay tests bind sequences, qualities, IDs, mate selection, gaps, and EOF to
  the same counted subset, including one-shot input and malformed frames.
- Resource-budget checks exercise combined reservations and allocation/resize
  failures, rather than testing nominal capacities alone.
- The existing full correctness suite and pinned products agree with the released
  baseline. The broader biological studies remain deferred under #129; this
  interface work does not satisfy or waive those studies.

The `kmer::store` candidate implements exact lookup/iteration, bounded batches,
non-copying filtered views, checked per-key/total overflow, move finalization,
and slot-buffer limits including resize peaks. These are owned table-buffer
measurements, not complete process budgets. The shared `kmer::budget` ledger
now enforces simultaneous stage reservations plus declared process overhead,
with separate temporary-disk guards. It uses explicit `reserved_*` names to
distinguish logical reservations from allocator heap and RSS. The candidate
store holds reservation guards through allocation, resizing, finalization and
destruction; actual graph/input/replay resource ownership remains open. Production dispatch, global resource
reservations, persistence, and replay integration remain open. The candidate is
not selected as the default backend before the
[frozen performance gates](../benchmarks/benchmark_results/v4-performance-baseline-20261109/README.md)
are measured.

The [identity/extraction verification](../benchmarks/benchmark_results/v4-identity-contracts-20261109/README.md) establishes the first item. Store/replay and
budget integration remain open until exercised by their actual implementations.

[Exact-store/ledger verification](../benchmarks/benchmark_results/v4-exact-store-contracts-20261110/README.md)
records the tested candidate. Production dispatch, persisted corruption/I/O
checks, and replay/producer integration remain required before closing #139.

`kmer::access::FilteredKmerCounts` is the assembly boundary for a borrowed exact
store or an explicit legacy adapter. Primer scans use a fallible visit operation;
graph membership/count queries propagate errors to the PCR caller. Neither
access path materializes a second global table. The legacy adapter widens its
already-saturated u32 values and does not claim exactness. `kmer::primers::PrimerKmers`
holds local directional candidates separately, with checked u64 addition.
Thresholds, graph edges, path minima and warning-source thresholds retain u64
counts; floating-point ranking statistics retain their existing interpretation.
The CLI defaults to the existing ingestion implementation pending #140 and
the measured promotion gates. Persisted I/O integrity and whole-run resource
enforcement are still required when those producers are integrated.

[Assembly-access verification](../benchmarks/benchmark_results/v4-assembly-access-20261110/README.md)
records exact/legacy product and repeat-warning agreement, large-count PCR and
injected scan/lookup failures. This does not close the remaining frontend,
replay, persisted-integrity or combined-resource gates.

## Direct-ingestion comparison under #140

`--count-backend memory` is a hidden development control until the frozen
performance gates have been evaluated. It requires `--chunks 0`; explicit
incremental counting stays on the released implementation. The direct path
owns one logical builder, emits rolling canonical windows without a per-read
vector, and finalizes by moving its allocation into the borrowed assembly path.
Valid-base totals are gathered during rolling extraction, avoiding another
sequence scan or a post-ingestion global count scan.

`fastq::FastqReader` reuses one byte arena per source reader. Original header,
sequence, separator and quality fields remain available until the next record.
LF and CRLF delimiters are removed equivalently; other original bytes are
retained. The parser preserves the released UTF-8 rejection and sampled
structural-validation policy. Length/quality consistency is therefore only
guaranteed at the configured validation points in this counting path; the
quality-preserving replay producer still requires its own per-record validation.
The record limit applies before allocation, including line delimiters. Shared
reservations charge old and new arenas during growth and are released on EOF
or failure. The selected-record digest uses a versioned domain, global index,
mate tag, then length-framed source/header/sequence/separator/quality bytes.
Effective caps remain in run metadata. Paired cap rounding and the released
inclusion of the first unmatched record at uneven EOF remain explicit.

The hidden `--count-slot-budget-mib` bounds owned count slots, including resize
peaks. `--fastq-record-bytes` bounds each complete FASTQ arena. These controls
are not whole-process RAM limits: allocator overhead, decompression/transport
buffers, graphs and evidence still require migration. Exact runs report count
precision, input identity and owned slot/record-buffer measurements under
`count_store`. `peak_count_input_reserved_bytes` is logical shared reservation
accounting; use allocator/RSS metrics separately. Default/legacy runs do not
emit an unperformed exact-store assessment.

The hidden `--count-fingerprint` performs an opt-in verification scan over all
full key/count pairs, with k included in each pair hash. Its versioned SHA-256
XOR digest, distinct cardinality and checked occurrence total support
order-independent backend comparisons. This scan is excluded from normal
operation and timed promotion runs. The paired comparison harness checks it
separately before timing each input against the immutable v3.2.1 binary.

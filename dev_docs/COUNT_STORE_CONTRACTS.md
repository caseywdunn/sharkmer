# v4 count-store, replay, and resource contracts

Implementation under [#139](https://github.com/caseywdunn/sharkmer/issues/139).
This document defines the interfaces required by #140–#145. The current branch
implements typed identities and positioned rolling extraction. Exact in-memory
and streamed sorted-file candidates now exercise the shared-ledger API.
Assembly consumers now use fallible borrowed count access and retain u64 values;
normal sPCR now selects exact hybrid counting, bounded replay and shared graph
limits. Integrated default acceptance remains open under the
[current resource policy](#normal-run-resource-policy-under-validation). Dated
candidate checkpoints below retain their original evidence scope; statements
about prototype defaults describe those checkpoints. The
[active plan](PLAN.md) records verification and completion.

## V4 resource boundary

**User decision, 2026-11-16:** keep v4 focused on exact counting and bounded
application buffers for ingestion, lookup, replay, assembly and retained
results. Combined application reservations and declared process headroom must
fit the planned envelope. Existing application body/record limits and owned
disk/descriptor controls remain; exact allocator overhead, physical filesystem
metadata and every process descriptor are not promised by those counters.

HTTP parsing, URL handling, connection pools, resolver work, sockets, TLS and
automatic content decoding belong to the client library. CLI/runtime/error
formatting and allocator internals likewise use headroom and platform limits.
Do not add per-allocation formulas, custom library admission hooks, replacement
subsystems or a portable RSS guarantee for them. The previous exhaustive
networking/whole-process accounting gate is **waived by the user 2026-11-16**
for maintainability and release scope. Existing application owners stay unless
a concrete simplification removes them with regression checks.

A network/status/decoder failure propagates as an operation failure. Memory
exhaustion in an infallible dependency allocation may abort; a kernel limit is
a last safeguard, not a guarantee of a typed error or immediate cleanup.
Failed work must not be accepted as complete or overwrite a previously
accepted generation. Retry/recovery tests verify that boundary.

Final-source product/count agreement, application-budget spilling/refusal,
remote-input and publication/recovery failures under existing process limits,
measured RSS/headroom, and the original performance comparisons remain
required. Historical checkpoints below describe their own stricter experiment;
references to still-unowned dependency internals are not new work requirements.
See [PLAN.md](PLAN.md) for current release gates and [HTTP scope](HTTP_WORKSPACE.md).

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

## Shared startup plan

`io::prepare_ingestion` builds a numeric table/parallel/spill layout and one
shared parent ledger after configuration validation and dry-run handling,
before output-directory and transaction creation. It allocates no count tables,
starts no workers and creates no temporary paths. Ingestion consumes the same
borrowed arguments, plan and ledger. Table, queue, sort, lookup, replay, graph
and combined memory ceilings retain their previous formulas; the result class
still shares the existing combined ceiling instead of adding an allowance.
The explicit legacy and incremental branches are unchanged.

The output transaction borrows that parent ledger from startup under the
existing development result-allowance control. Manifest YAML uses the audited
preflight before temporary-file creation and retains its reservation through
emitter destruction. New-generation transaction labels, OS paths, staged receipt
arrays and digest text now reserve under that same parent before allocation.
Manifest and generated recovery-journal views borrow their receipt owners;
journal emitter workspace is also reserved before temporary-file creation.
Verification cleanup paths and duplicate tables, recovery receipt indexes and
the complete restore plan now use separate reservations on the same parent.
All managed restore paths are admitted before public mutation; the sorted
borrowed index preserves lexical recovery order without copying receipt text.
Historical input buffers and retained typed models now use separate parent
reservations. Serde seeds admit scalar copies and receipt-vector growth before
allocation, including aliases; each parsed model holds its aggregate guard
through data destruction. The loader's scanner/event/anchor workspace also has an audited precharge
before deserializer construction; returned errors and process headroom remain
separate requirements.
The ledger now begins before panel collection and ENA/sample resolution.
Raw CLI/clap processing and complete HTTP/TLS transport remain separate owners.

A scoped memory observation starts before actual counter construction and ends
after finalization. It measures simultaneous run reservations during that
phase, including live baseline owners, while excluding completed earlier
startup workspace. The lifetime peak remains intact. This preserves the meaning
of `count_store.peak_count_input_reserved_bytes` when the shared ledger also
serves publication. Nested measurement is rejected; refusal, unwinding and
invalid accounting prevent a stale successful measurement.

Full startup-ledger verification (local evidence: `benchmarks/benchmark_results/v4-startup-ledger-20261113/README.md`)
passes 479 Rust tests per hasher and 224 Python tests, including 128 layouts
compared with source-frozen resource formulas, parent-identity/lifetime checks,
early-exit behavior and manifest refusal preserving completed outputs. All
1,024 native candidates and 16 frozen references pass under actual hard/soft
64-descriptor limits, with complete FASTA/PCR, count/input identity and stable
metadata agreement. Independent auditing verifies both complete sequence/PCR output
equivalence and unchanged temporary pool with the explicit planned lease.
**Scoped shared-ledger/manifest-workspace gate: met 2026-11-13.** Whole-run
resource enforcement remains open.

New-generation publication metadata verification (local evidence: `benchmarks/benchmark_results/v4-publication-metadata-20261113/README.md`)
passes 485 Rust tests per hasher and 224 Python tests. Borrowed manifest and
recovery-journal models preserve the prior YAML bytes; isolated allocator
checks cover owned paths and 2,000 staged receipts. Metadata refusal precedes
staged-file creation and can retry after competing owners release. Independent
fault/restart auditing passes 160 EIO/ENOSPC/SIGTERM cases (480 processes).
All 1,024 native candidates and 16 frozen references pass complete FASTA/PCR,
count/input identity and stable metadata comparisons under actual hard/soft
64-descriptor limits; independent auditing and retained-artifact verification
also pass. **Scoped new-generation publication metadata ownership gate: met
2026-11-13.**

Verification/recovery workspace (local evidence: `benchmarks/benchmark_results/v4-recovery-workspaces-20261113/README.md`)
passes full quality checks (489 Rust tests per hasher, 224 Python tests).
Prior-output cleanup vectors and paths, borrowed duplicate tables, sorted
receipt indexes and the complete recovery plan reserve under the parent.
Every managed destination and restore-source path is admitted before public
mutation. Tests preserve all interrupted public/backup files on typed refusal,
retry after competing memory releases and verify borrowed receipt identities
against a BTreeMap oracle. Eight isolated index allocator observations cover
0/1/13/2,000 receipts under both hashers. All 1,024 native candidates and 16
frozen references pass complete FASTA/PCR, count/input identity and stable
metadata checks under actual hard/soft-64-descriptor limits. All 160 fault cases
(480 processes), independent auditing and retained-artifact verification pass.
**Scoped recovery workspace ownership gate: met 2026-11-13.**

Historical input and typed-model ownership (local evidence: `benchmarks/benchmark_results/v4-parsed-metadata-20261113/README.md`)
passes full quality (493 Rust tests per hasher, 224 Python tests). An admitted
original input buffer feeds the borrowed-slice deserializer. Value seeds reserve
each scalar copy and receipt-vector growth before allocation, including aliases.
An aggregate owner retains actual backing capacities until its model data drops.
Historical format comparisons, alias-copy refusal/retry, destruction order and
real completed-output preservation pass. All 1,024 native candidates and 16
frozen references agree on complete FASTA/PCR, identities and stable metadata
under actual hard/soft-64-descriptor limits. All 160 fault cases (480 processes),
independent audits and retained-artifact checks pass. **Scoped input/typed-metadata
ownership gate: met 2026-11-13.**

Historical loader workspace (local evidence: `benchmarks/benchmark_results/v4-loader-workspace-20261113/README.md`)
now reserves scanner buffers/stacks, event arrays, scalar boxes and anchor/alias
dictionaries before decoding. The pinned source model covers array replacement,
ignored nesting, second-document overlap, Unicode escape expansion, repeated
tag prefixes, allocation headers and audited Rust BTree/Vec layouts. The guard
remains live through loader destruction and unwinding. Full checks pass 495
Rust tests per hasher and 224 Python tests; all 56 isolated loader observations
fit the reserved bytes and release their measured heap. All 1,024 native
candidates and 16 frozen references agree; 160 fault cases (480 processes),
independent audits and retained-artifact verification pass. **Scoped loader
workspace ownership gate: met 2026-11-13.** Boundary counting can
overreserve scalar text. Returned error reporting and other process owners
remain separate; these tests do not establish a whole-process RAM guarantee.

Combined-budget native integration (local evidence: `benchmarks/benchmark_results/v4-combined-budget-20261113/README.md`)
passes one million stdin reads, exact spill, original-pass replay and complete
reference sPCR with graph/result owners enabled. The original hard 256 MiB
address-space cap, 32 MiB table, 8 MiB sort, 64 GiB logical disk and timeout
remain unchanged; this trial adds 8 MiB shared/per-target graph ceilings and
actual hard/soft 64-descriptor limits. All seven FASTA outputs and PCR evidence
agree, every target fits its ceiling and owned temporary files are removed.
Independent input framing, count fingerprint, reference and receipt audits pass.
This is scoped integration evidence; the remaining owner table and production
resource/default gates are not closed by this one trial.

Publication copy workspace (local evidence: `benchmarks/benchmark_results/v4-publication-copy-20261113/README.md`)
uses one reusable 64 KiB parent-owned byte buffer for all backup or restore files.
Backup admission precedes backup-directory creation; recovery admission follows
complete validation and plan construction but precedes any public file deletion
or restoration. Reads retry interruption, writes handle short progress, and I/O
failure or unwinding retains the guard until the buffer drops. The buffer is
released before journal emission or final recovery-directory cleanup. Unbudgeted
runs preserve their existing standard-library copy path. Targeted checks pass
38 output tests and Clippy per hasher, including byte equality, refusal/retry,
short writes and panic cleanup; 12 isolated allocator observations exactly match
64 KiB and release all measured heap. Full quality passes 498 Rust tests per
hasher and 224 Python tests. All 1,040 native comparisons, 160 fault/restart
cases (480 processes), independent audits and archive checks pass. The original
fault-run limits remain unchanged; a corrected tracer regression requires actual
terminal wait status after a disappearing tracee and preserves the initial
harness failure. **Scoped copy workspace ownership gate: met 2026-11-13.**
Tempfile metadata, error reporting and broader process/filesystem ownership
remain open.

## Remote cache-miss record buffers

Exact ingestion now passes its record ceiling and shared input ledger to the
remote cache downloader. Cache misses use the same reusable bounded FASTQ
arena as original-source counting, reserving both buffers during growth before
allocation. The raw four-line byte ceiling applies before cache newline
normalization. EOF/completion flags, selected-record caps, gzip-member handling,
UTF-8 rejection and verified cache publication retain their meanings. The
explicit legacy downloader keeps its existing implementation.

Cache-record checks (local evidence: `benchmarks/benchmark_results/v4-cache-records-20261113/README.md`)
pass all 32 cache tests and Clippy under both hashers. Real HTTP tests verify
successful publication and preservation of prior cache data, metadata and
foreign files on actual parent-budget or record-limit refusal. Input guards
release on success and error. Full checks pass 507 Rust tests per hasher and
225 Python tests, including cache miss/hit, literal count/input-frame oracles
and original-pass replay. All 36 native comparisons, 160 fault cases and
retained source/archive checks pass. **Scoped cache-record buffer gate: met
2026-11-13.** This
scopes record arenas only. The later
cache-metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-cache-metadata-20261114/README.md`)
bounds sidecar input/output at 10 MiB, owns input capacity, audited YAML
workspace and surviving strings on the shared retained-results ledger, and
keeps closed temporary sidecar bytes charged through cache-directory sync.
Existing metadata/default/alias semantics and typed resource causes are
preserved. Full checks pass 560 Rust tests per hasher and 225 Python tests;
all 80 cache/HTTP component tests pass under actual 64-descriptor/512 MiB
limits, with independent auditing and retained source/raw evidence.
**Scoped cache-metadata ownership gate: met 2026-11-14.**

Compressed-cache extent ownership (local evidence: `benchmarks/benchmark_results/v4-cache-data-20261114/README.md`)
also uses the common temporary writer: admit growth before I/O, preserve typed
refusal through gzip and buffered output, then verify, sync and close before
publication. Both data and sidecar charges remain live through cache-directory
sync. Confirmed live rollback syncs before releasing the new data charge.
Real HTTP quota tests cover fresh/replacement refusal before a complete gzip
file and after data but before metadata, preserving old/foreign files without
temporary residue. Success measures the exact combined data/metadata extent.
Full checks pass 562 Rust tests per hasher and 225 Python tests; all 84 native
cache/HTTP component tests pass with independent auditing and retained evidence.
**Scoped compressed-cache extent gate: met 2026-11-14.** That earlier checkpoint
scopes new compressed-file extents; registered pair recovery is verified below.

The registered cache-recovery checkpoint (local evidence: `benchmarks/benchmark_results/v4-cache-recovery-20261114/README.md`)
now stages on the cache filesystem under an explicit ownership marker and
immutable preparation/journal receipts. Backups survive failed restoration;
public names are restored through separate links. Completed acceptance and
partial cleanup have distinct restart paths. Recovery preadmits the common
inode inventory and preserves unregistered/modified evidence. Full checks pass
570 Rust tests per hasher and 225 Python tests; 100 hard-limit cache/HTTP
component tests pass independent auditing. The full native matrix and independent
audit pass 544 publication/restoration cases (1,488 invocations, 212 actual
SIGTERMs) under actual 64-descriptor/512 MiB limits. Both old backups survive
failed restoration; cold retry preserves the accepted pair, foreign files and
complete cache metadata/FASTQ without registered leftovers. Three exact source
archives and all 16,196 raw payloads pass retained verification, including the
original hard-link restoration failure and failed harness versions.
**Scoped registered-cache publication/restoration I/O/SIGTERM gate: met
2026-11-14.** This is no whole-system power-loss proof. Broader process/library/
physical-filesystem requirements and final-source integration remain open.

The application FASTQ codec/buffer checkpoint (local evidence: `benchmarks/benchmark_results/v4-codec-20261114/README.md`)
now owns the source buffer, gzip state, optional member headers and resize
overlap, cache compressor and final flushing. Pinned-source formulas use public
backend state sizes plus audited nested allocations; constructors reserve before
allocation or eager header parsing. These owners borrow the existing shared
retained-results allowance without increasing the parent or record-buffer limit.
Targeted checks and isolated heap observations pass under both hashers. The
largest decoder peak is 313,655 bytes within its 313,927-byte envelope; the
encoder peak is 360,296 within 360,456, with complete allocation release.
Full checks pass 576 Rust tests per hasher and 225 Python tests. All 112
hard-limit component tests and 544 native cache fault/recovery cases (1,488
invocations, 212 actual SIGTERMs) pass independent auditing. Retention verifies
two exact source archives, 118 live joins, 45 dependency source/manifest files,
the dependency archive and all 15,447 raw native payloads. **Scoped application
FASTQ-buffer/gzip backing gate: met 2026-11-14.** These are logical backing
reservations and requested-heap observations. HTTP/TLS/content-encoding and
preprocessing, allocator/OS overhead and complete descriptor/physical ownership
remain separate requirements.

The ENA ownership checkpoint (local evidence: `benchmarks/benchmark_results/v4-ena-20261114/README.md`)
preadmits raw response growth, exact lossy-UTF-8 expansion, URL-vector/string
models and sample spelling on the original count ledger. Owners survive move
into the uncached source plan. Reading preserves the original 10 MiB logical
body limit, complete-body error checking and first-row/URL/sample semantics.
Dry runs use a separate bounded metadata preview without evaluating unused
numeric count settings; an initial full-test regression and its unchanged
passing early-exit test are retained. Full checks pass 583 Rust tests per hasher
and 225 Python tests. All 126 hard-limit component tests and 544 native workspace
I/O/retry cases pass independent auditing, as do both source archives, 76 audited
dependency files and all 13,342 raw payloads. Isolated observed and reserved
requested-heap bytes agree and release to baseline. **Scoped ENA body/model/sample
ownership gate: met 2026-11-14.** Request/header/TLS/content-decoding allocations,
earlier panel/CLI preprocessing and broader process/filesystem ownership remain
open. The preview allowance never increases an ordinary run's ledger.

The automatic HTTP gzip checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-codec-20261114/README.md`)
reserves the pinned client's eager content decoder before constructing ENA/read
requests. Its guard survives redirects, response-to-reader transfer, late read
failure and retained HTTP-status errors that own a response. The transport
reader or error drops before the guard. Nested HTTP/application gzip uses two
simultaneous owners on the unchanged original ledger. Source-audited bounds and
ten isolated observations per phase fit a 305,671-byte decoder envelope and
release completely. Full checks pass 589 Rust tests per hasher and 225 Python
tests. All 138 hard-limit component tests and 544 native cache I/O/SIGTERM cases
(1,488 invocations, 212 actual SIGTERMs) pass independent auditing. The exact
122-file source archive, 76 dependency files and all 15,447 raw payloads pass
retention. **Scoped automatic HTTP gzip gate: met 2026-11-14.** This does not
cover request/URL/header, resolver/socket/TLS, earlier panel/CLI preprocessing,
whole-process descriptor/physical resources or final-source promotion. Remote
panel fetching remains part of its earlier preprocessing path.

The later response-metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-metadata-20261114/README.md`)
admits status/header storage before constructing requests, including redirects
and formatted parser errors. Success consumes Response before releasing this
owner; retained errors destroy their library payload before the guard. See
[HTTP workspace](HTTP_WORKSPACE.md) for the pinned source-derived bound and
fixture composition. Full checks, original-kernel-cap component/native audits
and source/dependency/raw retention pass. **Scoped response status/header and
error-lifetime gate: met 2026-11-14.** Request/URL/Agent/pool and
resolver/socket/TLS remain separate owners; this does not close fresh
integrated cache/publication I/O/SIGTERM or whole-run resource gates.

The chunk-framing checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-chunk-framing-20261114/README.md`)
replaces the pinned client's unbounded size-token Vec with a constant-size
streamed parser, preserving grammar and source-error ordering. The licensed
client source is retained locally with its original version/features and
documented patch. Full checks, zero-allocation observations, native low-RSS
replay, fresh cache I/O/SIGTERM auditing and complete evidence retention pass.
**Scoped chunk-framing and current-cache compatibility gate: met 2026-11-14.**
Request/URL/Agent/pool, resolver/socket/TLS, raw CLI/clap, other process owners,
complete physical/descriptor accounting and final-source promotion remain open.

The original request-text checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-request-text-20261115/README.md`)
admits the raw URL and GET method copies before construction on the same
ledger. Their exact reservation survives an unused builder or consuming call
until those original Strings drop. Parsed/redirect URL storage remains a
different owner. Full checks pass 634 Rust tests per hasher and 225 Python
tests. Independent HTTP auditing passes 26 native invocations and 84 isolated
allocation observations under unchanged 64-descriptor/512-MiB caps. Fresh cache
auditing passes 544 cases, 1,488 invocations and 212 actual SIGTERMs. Source,
dependency, raw retention and the 51-file seal pass. **Scoped original
request-text and current HTTP/cache integration gates: met 2026-11-15.**
Parsed/redirect URLs, request headers, Agent/pool, resolver/socket/TLS and the
broader process, storage, control and final-source gates remain open.

The URL/redirect workspace checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-url-workspace-20261115/README.md`)
now admits each initial parse/redirect join and retains source-derived
URL/IDNA workspace, clones and history until response/error destruction.
Failed admission returns its original typed error. The pinned client hook
leaves URL syntax, redirects, configuration and diagnostics unchanged. Full
checks pass 639 Rust tests per hasher and 225 Python tests. Native HTTP auditing
passes 36 invocations and 372 allocation observations at unchanged
64-descriptor/512-MiB caps; fresh cache auditing passes 544 cases, 1,488
invocations and 212 actual SIGTERMs. Three exact source archives, 749
dependency/compiler bindings, licensed dependency and raw retention, and the
73-file seal pass. **Scoped URL/parser/clone and current HTTP/cache integration
gates: met 2026-11-15.** Derived request method/header/prelude copies,
Agent/pool/resolver/socket/TLS and broader process/storage/control/final-source
gates remain open. See [HTTP workspace](HTTP_WORKSPACE.md) for the bound and
positive-fixture composition.

The fixed Agent checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-agent-backing-20261116/README.md`)
admits default/panel client container backing before construction, using one
ledger and preserving timeouts. Reader/error leases follow possible Weak state
backing without keeping the actual pool alive. Isolated warm constructor
observations fit the source/layout-derived allowance; global TLS initialization
is explicitly excluded. Full checks pass 647 Rust tests per hasher and 225
Python tests. Independent auditing passes 52 HTTP invocations, 388 allocation
observations, 16 wire pairs and the unchanged 544-case/1,488-invocation cache
matrix with 212 actual SIGTERMs at 64 descriptors/512 MiB. Source/dependency/raw
retention and the 57-file seal pass. **Scoped fixed Agent/Weak backing and
current HTTP/cache integration gates: met 2026-11-16.** Pool allocations,
resolver/socket/TLS and broader resource/final-source gates remain open.
See [HTTP workspace](HTTP_WORKSPACE.md) for the bound and lifetime contract.

The derived GET request checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-request-buffers-20261116/README.md`)
now charges automatic headers before parsing and each Unit's method/header/
authentication/prelude workspace before construction. Actual header capacities
and User-Agent length drive admission; a separate guard releases after the
call while URL/response owners remain with retained errors. Full checks pass
643 Rust tests per hasher and 225 Python tests. Native auditing passes 44 HTTP
invocations and 16 wire pairs with exact live/released reservation checks at
unchanged 64-descriptor/512-MiB caps. Fresh cache auditing passes 544 cases,
1,488 invocations and 212 actual SIGTERMs. Source/dependency/raw retention and
the 61-file seal pass. **Scoped derived-GET request and current HTTP/cache
integration gates: met 2026-11-16.** These wire/ledger checks do not claim an
isolated request-only heap measurement. Agent/pool/resolver/socket/TLS and the
broader resource/control/final-source gates remain open.

## Staged output byte ownership

Budgeted output transactions now charge staged FASTA and stats bytes to the
same temporary-disk ledger as counts and replay. The write-only owner reserves
the attempted extent before each write and accounts for short progress. A
configured-limit failure remains typed even if a callback catches its I/O
error. Permanent I/O failures poison finalization; changed extents are rejected.
The file is flushed, synced and closed before receipt hashing.

Private bytes remain charged through confirmed removal. Charges follow public
links until durable acceptance or completed rollback. A failed recovery keeps
private payloads and their charges; failed cleanup conservatively retains the
charge for the rest of the process. Accepted final products become persistent
output and leave the temporary allowance. Interrupted writes remain retryable
under the [Write contract](https://doc.rust-lang.org/std/io/trait.Write.html#method.write_all),
holding the attempted extent throughout retry.

Staged-output verification (local evidence: `benchmarks/benchmark_results/v4-output-disk-20261113/README.md`)
passes 514 Rust tests per hasher and 225 Python tests. Independent native
audits pass 36 default/reference checks, all 160 original publication fault
cases, 32 added FASTA/stats write/fsync error cases and 18 interrupted-write
invocations. The latter preserve a two-case pre-fix failure and verify 16
actual EINTR retries against complete reference output. Retained-artifact
verification passes 170 retained files, nine source archives and 107 source
bindings. **Scoped staged-output disk gate: met 2026-11-13.**
Memory-only unit fixtures now allow 1 MiB of
staged disk; their original memory criteria and every native protocol ceiling
remain unchanged. Dedicated disk tests use exact 3/6-byte ceilings.

This scopes staged logical extents. Backup/cache disk ownership, manifests and
filesystem overhead, complete descriptor accounting and pre-journal orphan
cleanup remain open.

### Prior-generation backup copies

The next development checkpoint also admits the complete prior-generation
backup extent before creating its directory. It sums regular-file extents
with checked arithmetic; each copy is bounded by the remaining admitted bytes
and rejects source-size changes. The transaction keeps the aggregate charge
until confirmed backup-directory removal, including failed publication and
recovery. A directory-creation refusal releases the unused allowance.

Both corrected targeted suites pass 50 output/recovery tests and Clippy. The
source-bound checkpoint (local evidence: `benchmarks/benchmark_results/v4-backup-disk-20261113/README.md`)
records the first full checks (517 Rust tests per hasher and 225 Python tests),
passing native/fault checks and a new native backup-disk refusal/retry trial.
The latter
uses 1,024 synthetic targets on both hashers at threads 1/8, with space for
one staged generation but insufficient space for its independent backup.
All accepted bytes and foreign files must survive refusal; a retry with the
original disk allowance must reproduce complete products and evidence. The
first trial preserved all accepted bytes on disk refusal and reproduced the
complete output on retry, but left backups because decoding its own recovery
journal exceeded the unchanged memory allowance. Corrected accepted cleanup
uses captured hashes and existing typed receipts, verifies the complete backup
inventory and accepted generation, then removes recovery data. Missing or
modified evidence preserves recovery files and their charge. Historical
restart recovery keeps its full on-disk parsing and validation. The corrected
full/native retry and independent audits pass with the original criteria and
limits. Combined retained verification passes 153 evidence files, seven source
archives, 108 source joins and 188 local links. **Scoped new-backup-copy and
accepted-cleanup gate: met 2026-11-13.** This does not yet account for
fresh-process adoption of recovery files, restoration temporary copies,
new manifests/journals or the other storage requirements above.

### Restoration copy admission

Rollback now admits the total additional copy extent, including restoration
of the prior manifest, after its paths and shared copy buffer are admitted and
before any public deletion. Each copy is limited to the planned extent and
verified against its receipt before publication. Final checks reuse owned
paths and stack digest buffers, avoiding new application metadata reservations
after mutation. The allowance releases only after the complete prior
generation is verified and its manifest durable. A failed recovery retains
the charge conservatively; configured-limit refusal leaves interrupted state
unchanged and retryable.

The restoration checkpoint (local evidence: `benchmarks/benchmark_results/v4-restore-disk-20261113/README.md`)
records full initial checks (521 Rust tests per hasher, 225 Python), native
compatibility, fault and large quota-refusal auditing. Its strict raw-state
audit exposed transient NFS inodes from unlinking open temporary files. The
corrected owner closes its file before dropping its path on errors, and hands
only closed paths to publication. Both corrected targeted suites pass 53
output/recovery tests and Clippy. Corrected full checks pass 522 Rust tests per
hasher and 225 Python tests; all native compatibility, fault and large quota
audits pass. Combined retained verification passes 213 evidence files, six
source archives and 108 source joins. The original failed audits and exact
source versions remain preserved. **Scoped restoration-copy and close-order
gate: met 2026-11-13.** Fresh-process ownership of existing
recovery data, other publication/cache resources, physical filesystem overhead
and pre-journal orphan cleanup retain their existing requirements.

The interrupted-restoration checkpoint (local evidence: `benchmarks/benchmark_results/v4-restore-signals-20261113/README.md`)
now creates new restore copies beside their verified backups, within the
receipted staging directory. A native before/after study and independent audits
verify 24 actual SIGTERM cases per version, covering copy writes, syncs and
publication across two hashers and threads 1/2/4/8. The preceding frozen
version restores the prior generation but leaves 24 root copies; the new
version restores every prior output byte and leaves none. Full checks pass
523 Rust tests per hasher and 225 Python tests. Compatibility and existing
I/O-fault/quota audits also pass. Corrected retained
verification passes 146 evidence files, three source archives and 108
source joins. **Scoped restoration-signal cleanup gate: met 2026-11-13.**
Historical receipts and failure criteria remain preserved; unrelated
unreceipted old files are not reclaimed.

### Existing recovery storage adoption

Fresh-process recovery now inventories the validated journal's private staging
tree, its backups and abandoned restore copies, plus newly published public
payloads that have not been accepted. It admits those existing bytes together
with all planned restoration copies before changing public files. Unix
device/inode identities count a hard-linked private/public payload once;
other platforms conservatively count each path. The fixed-depth inventory
rejects unlisted entries and nonregular payloads rather than following links
or removing foreign subdirectories.

In-process rollback credits the same transaction's live staged/backup guards
on the same ledger; the previously unowned journal is always added. All
memory/copy owners precede the final disk admission. A configured refusal
releases unused owners, preserves interrupted files and allows retry with the
same budget. Once mutation can begin, uncertain failure conservatively keeps
the allowance; confirmed cleanup releases it after output-directory sync.
Accepted persistent public outputs are outside this temporary allowance.

The adoption checkpoint (local evidence: `benchmarks/benchmark_results/v4-recovery-adoption-20261113/README.md`)
passes 57 output/recovery tests and Clippy per hasher, and full development
checks with 526 Rust tests per hasher and 225 Python tests. Native compatibility,
publication/staged/restoration faults and actual restoration SIGTERM audits
pass. The larger quota trial preserved interrupted bytes on refusal and
restored prior outputs on retry, but directory removal left 17 private files.
The failed source and raw evidence are being retained. The correction owns all
private removal paths before mutation instead of enumerating while deleting,
and includes already-restored public files while their manifest is incomplete.
The corrected checkpoint (local evidence: `benchmarks/benchmark_results/v4-recovery-completion-20261113/README.md`)
passes 58 output/recovery tests and Clippy per hasher, full checks with 527 Rust
tests per hasher and 225 Python tests, and native compatibility/I/O/SIGTERM
audits. Large quota, restored-public and completed-replacement trials and
independent audits also pass with the original resource ceilings and no
unreclaimed private files. Combined retained verification passes 127 evidence
files in nine groups, three source archives and 109 current-source joins. It
reopens and checks every payload in eight new raw archives and the retained
historical quota archive. **Scoped receipted-recovery storage and planned
cleanup gate: met 2026-11-13.**
This does not complete pre-journal/count/replay orphan ownership, new
manifest/journal write accounting, other publication/cache/library resources,
physical filesystem accounting or process descriptor ownership.

### New manifest and journal byte ownership

The metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-metadata-disk-20261113/README.md`)
uses the staged writer for temporary manifest/journal bytes: admit attempted
growth before I/O, account accepted short writes, retry interrupted writes and
recover the original typed refusal through YAML serialization. The writer
closes before hashing and path publication. A manifest's allowance releases
after durable publication as run-status metadata; the transaction retains a
journal guard until removal. Recovery credits the matching journal guard
separately from staged/backup payloads; a fresh process supplies no credit.

Counting summaries now use a separate non-nested disk observation from counter
startup through store finalization. It includes live baseline owners and any
overlapping replay, while excluding completed startup and later publication.
The whole-run ledger and its lifetime peak retain all charges. Native before/
after measurement reproduces 52,536 bytes attributed to counting by the old
lifetime snapshot and 8,512 bytes for the same count under the corrected phase
observation. This corrects measurement scope, not physical storage usage.

Local full checks pass 534 Rust tests per hasher and 225 Python tests. Targeted
checks pass 63 output/recovery cases, two disk-phase tests and the real-startup
production-layout test per hasher, plus Clippy. All native runtime and
independent audits pass, including 96 metadata write/sync fault/retry cases
and both measurement comparisons. The two stale audit provenance references
were corrected in separately frozen verifiers; failed originals remain
preserved. Combined retained verification passes 187 evidence files in 11
groups, six source archives and 110 current-source joins, with every payload
in ten new raw archives and the reused historical archive checked.
**Scoped temporary-metadata and count-disk measurement gate: met 2026-11-13.**
The staging preparation checkpoint (local evidence: `benchmarks/benchmark_results/v4-preparation-20261113/README.md`)
passes its frozen-source registered-staging SIGTERM cleanup gate and retained
verification. It preserves the actual earlier-restore-copy failure and the
corrected native before/after proof. The
registered count/replay candidate (local evidence: `benchmarks/benchmark_results/v4-registered-workspace-20261113/README.md`)
passes its frozen-source registered-workspace SIGTERM/quota gate with full
checks, 240 independently audited native cases and retained verification.
The later missing-receipt guard is included in the cache-metadata/data full
checks. The registered-workspace I/O matrix (local evidence: `benchmarks/benchmark_results/v4-cache-data-20261114/README.md`)
passes 544 cases and 1,360 invocations under actual 64-descriptor/3 GiB limits,
both hashers and threads 1/2/4/8. It injects 176 EIO, 176 ENOSPC, 176 EINTR and
16 short-write cases across 11 registry/root/owner boundaries. Retry preserves
the accepted generation and foreign files without registered leftovers.
Independent auditing and verification of all 13,342 archived raw payloads pass.
**Scoped registered-workspace I/O gate: met 2026-11-14.** Broader final integration
remains open. This does not complete
whole-system durability proof, remaining cache/library resources, complete
process descriptor/physical-filesystem ownership, public controls, final-source
performance or Git/CI integration.

## Optional process address-space safeguard

Normal sPCR commands now request a 12 GiB Linux address-space ceiling after CLI
parsing and before panel/ENA decoding. Integrated default acceptance remains
open. Explicit incremental counting and developer backend selectors preserve
their separate behavior.

`SHARKMER_PROCESS_ADDRESS_SPACE_LIMIT_MIB` overrides this policy. A positive
integer is installed before CLI decoding; `0` preserves inherited limits and
disables automatic installation. Absence uses the normal-run default or leaves
inherited limits unchanged on an explicit developer/incremental path.
A positive integer requests a Linux virtual-address-space ceiling;
positive requests on other platforms fail explicitly. The implementation uses
the existing safe Rustix interface, takes the minimum of requested and inherited
soft/hard ceilings, rounds down to a page, rejects a ceiling smaller than the
existing image, lowers both kernel limits and verifies their readback.

Stats omit `process_address_space` unless the safeguard was installed. When
present, its byte-valued fields mean:

| Field | Meaning |
| --- | --- |
| `requested_limit_bytes` | Requested limit from the environment, `--memory-budget`, or the normal-run default. |
| `inherited_soft_limit_bytes`, `inherited_hard_limit_bytes` | Limits observed before installation; omission inside this object means unlimited. |
| `effective_limit_bytes` | Actual installed, page-aligned soft and hard ceiling. |
| `observed_address_space_bytes` | Virtual image size observed before installation, not a peak or an RSS measurement. |

Linux [RLIMIT_AS](https://man7.org/linux/man-pages/man2/getrlimit.2.html) bounds
virtual mappings, including library allocations and thread stacks. It excludes
kernel caches and other processes. Exhaustion can terminate the process during
stack growth or an infallible allocation; Rust's default
[allocation-error handler](https://doc.rust-lang.org/std/alloc/fn.handle_alloc_error.html)
aborts. Owned admission, controlled spilling and publication recovery remain
necessary. Disabling this safeguard does not disable owned class budgets. The
safeguard does not prove a portable whole-run RAM guarantee.

Process-cap verification (local evidence: `benchmarks/benchmark_results/v4-process-memory-20261113/README.md`)
passes 501 Rust tests per hasher, 224 Python tests, 1,040 compatibility
invocations and 160 fault cases (480 processes). A separate 78-invocation
matrix observes actual kernel limits, disabled-field omission, inherited caps
and early refusal without changing completed outputs. A million-read native
trial installs its own 256 MiB cap from an unlimited parent and preserves all
seven reference sequences and complete PCR results. Its count and replay
files each exceed that ceiling. Scoped runtime checks pass; retained-evidence
verification is recorded alongside these results. Production resource choices,
remaining owned allocations and physical filesystem accounting remain open.

## Normal-run resource policy under validation

The frozen parameter pilot (local evidence: `benchmarks/benchmark_results/v4-default-resources-20261113/README.md`)
verified all ten calibration inputs and 287 records with a 12 GiB process
ceiling, 8 GiB table, 64 MiB sort, 1 GiB shared graph allowance and 256 MiB per
target. Peak native RSS was 3,348,764 KiB, allocator heap 3,422,623,873 bytes
and target buffers 14,723,232 bytes. These are observed peaks, not maximum
requirements for arbitrary inputs or biological validation.

Normal commands now use this policy with the exact hybrid counter. Omitted
graph controls select these allowances; explicit zero disables graph
accounting. Explicit developer backend selectors preserve their opt-in controls,
and `--chunks > 0` keeps the legacy path. Before creating output, normal-run
preflight requires the entire shared-ledger ceiling plus 2 GiB planning
headroom to fit the effective process ceiling. The public
`--memory-budget` option uses the smaller of its value and any inherited
Linux limit. After reserving headroom and fixed buffers, it divides the rest
among the count table, sort buffer, and assembly graphs. The per-target
allowance follows the shared graph budget. An explicit `12G` budget uses
the established defaults when they fit. `--temp-disk-budget` sets the
existing shared temporary disk allowance. An oversized plan fails before it
can change completed output.

The frozen integrated checkpoint (local evidence: `benchmarks/benchmark_results/v4-production-defaults-20261113/README.md`)
passes both complete suites, actual normal-command cap/incremental/refusal
checks, the original 100-invocation speed/RSS comparison and its conditional
80-case matrix, with independent audits and retention. **Scoped normal-default
snapshot gate: met 2026-11-13.** The elapsed geometric-mean ratio is
0.7760642119, with every original per-input wall/RSS criterion met. These
results bind its 106-file source, preceding cache and staged/backup ownership.
Final-source default promotion remains not met: application-resource and
failure/recovery integration, measured process headroom, performance revalidation
and Git/CI remain required under the resource boundary above. Exhaustive
physical-filesystem and dependency allocation accounting are outside v4 scope.

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
and owned-table limits including resize peaks. These are owned table-buffer
measurements, not complete process budgets. The shared `kmer::budget` ledger
now enforces simultaneous stage reservations plus declared process overhead,
with separate temporary-disk guards. It uses explicit `reserved_*` names to
distinguish logical reservations from allocator heap and RSS. The candidate
store holds reservation guards through allocation, resizing, finalization and
destruction; actual graph/input/replay resource ownership remains open. Production dispatch, global resource
reservations, persistence, and replay integration remain open. The candidate is
not selected as the default backend before the
[frozen performance gates](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-performance-baseline-20261109/README.md)
are measured.

Configured-limit failures now carry a downcastable `ResourceExhausted` error
through contextual error propagation. `ResourceKind` distinguishes a particular
stage's memory ceiling, combined run memory, and temporary disk. The error
records required and allowed bytes; projected usage uses u128 so a request
beyond u64 is reported without wrapping accounting. A rejected reservation or
resize leaves live charges and previous peaks unchanged. The count-table growth
precheck uses the same type and preserves accepted counts. Invalid accounting,
wrong-k queries, count overflow, allocation failures and physical I/O errors
remain separate errors. Future spill/concurrency logic must inspect the type
and resource, rather than matching messages or treating every failure as a
reason to spill. External counting now exercises this foundation as described
below; whole-run resource enforcement remains open.

Typed-exhaustion verification (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
records both-hasher context/downcast, failure-preservation and retry tests.
`kmer::temporary::OwnedDiskFile` now applies the disk ledger to actual temporary
files. It reserves logical extents before writes, rejects out-of-range reads,
poisons failed reads/writes, and releases a charge only after confirmed removal.
If removal fails, the charge stays live conservatively for the rest of the run.
The current descriptor-budget revision also reserves one shared temporary-file
slot before each count/replay file is created. Its `TemporaryFileReservation`
lasts until the actual descriptor closes, including failed-open rollback and
unlink errors; an unlink failure retains disk bytes but releases the closed
descriptor. `TemporaryFilesExhausted` reports required/allowed file counts and
is distinct from byte-limit and physical I/O errors. Run/merge files, count
owners and replay use the same ledger. Exact frontends request at most 128 such
files and perform a process-aware startup preflight. On Unix it reads the OS
soft limit; Linux also counts `/proc/self/fd` descriptors, excluding the census
directory. The exact run now prepares that ledger before opening its output
lease. It subtracts one explicitly planned output-lease descriptor and keeps
the existing 16 additional slots for other stages/libraries, saturating at
zero. Where an open-file census is available, the resulting temporary pool
is unchanged from the former census that already included the output lease.
Without a census, the known future lease is conservatively deducted while the
reported observation remains absent. Generic
conformance budgets keep explicit caller-supplied limits. A zero pool permits
a run opening no count/replay files; its first required temporary file fails
before creation.

`count_store.process_file_limits` reports requested/effective temporary-file
limits, `other_file_allowance` and performed observations as `process_soft_limit`
and `process_open_at_start`. Prepared runs also report
`output_lease_allowance=1`; standalone counting preflights omit it when no future
output lease is planned. The observed census stays literal, rather than adding
a descriptor that is not open yet. An absent OS value is unlimited or unavailable;
an absent census means no census was reported. No global descriptor peak is
implied. Input/output/cache/library opens are not individually controlled yet,
and later OS-limit changes remain physical I/O failures. The
process-aware revision (local evidence: `benchmarks/benchmark_results/v4-process-file-limits-20261113/README.md`)
passes 392 Rust tests per hasher and 224 Python tests, including actual finite
child limits, inherited-descriptor census, completed spill/replay/PCR and
pre-creation refusal without partial publication. Its 1,512-case native input
matrix and independent audit also pass under actual hard/soft 64-descriptor
limits. This establishes startup preflight, not complete per-stage
process-descriptor ownership.
Full descriptor checks (local evidence: `benchmarks/benchmark_results/v4-file-descriptors-20261112/README.md`)
pass 383 Rust tests per hasher and 224 Python tests, along with native
I/O-fault/interruption auditing. Complete production FD/process accounting
remains open.
File-owner verification (local evidence: `benchmarks/benchmark_results/v4-owned-disk-file-20261112/README.md`)
includes actual truncation and failed-removal cases. Filesystem allocation
rounding/metadata, descriptor/concurrency limits and measured headroom still
need production accounting; count/replay framing, integrity and integration
remain required.

Primary graph-buffer ownership and target admission now have
source-bound verification (local evidence: `benchmarks/benchmark_results/v4-graph-buffers-20261113/README.md`):
401 Rust tests per hasher and 224 Python tests. Hidden development controls
`--graph-buffer-budget-mib` (0 disables) and `--graph-buffer-target-mib` (256)
reserve each active target's allowance on the parent count/graph ledger.
A child ceiling charges actual primary node/edge buffers before growth,
including overlapping replacements and independently owned seed clones. The
last child graph owner releases the parent allowance and target admission slot.
The backing representation is pinned to petgraph 0.6.5. Optional
`assembly_graph_buffers` and per-target `graph_buffers` stats report these
performed engineering measurements; disabled assessment fields are omitted.
Eight exact candidates at threads 1/2/4/8 preserve counts, FASTA and full PCR
evidence. A real 25 kb target refuses a 1 MiB buffer allowance and completes
at 8 MiB. Lookup maps, traversal/pruning workspaces, primer/evidence allocation,
allocator overhead and whole-process memory remain outside this control.

The node-lookup migration (local evidence: `benchmarks/benchmark_results/v4-graph-lookups-20261113/README.md`)
uses pinned hashbrown buffers under the same child ceiling, with precharged
replacement growth and unchanged stable IDs. Full integrated verification
passes 412 Rust tests per hasher and 224 Python tests, including eight exact
candidates at threads 1/2/4/8, actual growth refusal and completed-generation
preservation. Optional graph-buffer measurements now include node-map buffers;
they still do not describe all target allocations. The subsequent
workspace migration (local evidence: `benchmarks/benchmark_results/v4-graph-workspaces-20261113/README.md`)
passes full quality checks (421 Rust tests per hasher and 224 Python tests).
It charges frontiers, compact directional/visited flags, repeat-marker sets,
iterative component traversal, pruning/median buffers and debug summaries
before growth. Independent graph oracles preserve component membership and
reachability policy. Native comparison and independent auditing pass all 128
candidate cases plus two frozen previous-build references, preserving full
FASTA/PCR evidence across eight backends, threads 1/2/4/8 and both hashers with
the allowance enabled/disabled. Independent raw-record/count oracles also
pass. **Scoped workspace conformance gate: met 2026-11-13.** Target peaks
include these workspaces. The subsequent
primer-buffer migration (local evidence: `benchmarks/benchmark_results/v4-primer-buffers-20261113/README.md`)
passes 428 Rust tests per hasher and 224 Python tests. It charges encoded
variant levels, reverse-orientation sets, capped candidate heaps, exact
directional maps and seed/count/median buffers before allocation. Preflight
borrows the bounded primer and uses fixed stack arrays. Streaming top-candidate
selection preserves the old mismatch/abundance/key priority while avoiding a
full matching table before capping; a broad scan keeps two candidates per
direction within 4 KiB of logical buffers. Native comparison and independent
count/record-hash auditing pass 512 candidate cases plus eight references,
including exact/IUPAC primers, caps 1/2/8, eight backends, both hashers and
threads 1/2/4/8 with target ownership enabled/disabled. All FASTA/PCR evidence
agrees with the frozen prior build. Existing expansion/work limits and
biological policy remain unchanged. **Scoped primer buffer/selection gate:
met 2026-11-13.** The subsequent
path-search migration (local evidence: `benchmarks/benchmark_results/v4-path-buffers-20261113/README.md`)
passes full quality checks (434 Rust tests per hasher and 224 Python tests).
End membership, visit counts, working path/child stacks, completed path copies
and edge-count/index/read-support metric scratch share the target ceiling.
Consuming iterators retain backing reservations until storage is destroyed;
nested path owners release independently. Ordinary DFS frames hold four
children inline; larger manually constructed degrees use charged spill and
stable-order sort scratch. Traversal preserves score ties, path/state/visit
limits and repeat rejection. Path statistics preserve their original floating
summation order before an in-place median. The extension lookup is released
before pruning and traversal, avoiding unnecessary overlap. The existing
25 kb test still completes at 8 MiB and refuses 1 MiB without publication or
temporary leaks. Native differential comparison and independent auditing pass
512 candidate cases and eight frozen prior-build references, preserving all
FASTA/PCR evidence and exact count/selected-read identities across the same
exact/IUPAC primer, cap, backend, thread, hasher and allowance matrix.
**Scoped path-buffer conformance gate: met 2026-11-13.** The subsequent
read-evidence/bubble migration (local evidence: `benchmarks/benchmark_results/v4-evidence-buffers-20261113/README.md`)
passes full quality checks (440 Rust tests per hasher and 224 Python tests).
Threading charges lookup/support/link maps, positioned and reverse windows,
nested contiguous runs, evidence sets, fragment grouping and paired-link
copies under the same target ceiling. Replay still holds one pending fragment;
the legacy slice adapter preserves arbitrary fragment order and duplicate-mate
merging. Bubble detection, nested branch buffers, ranking and edge preferences
also retain owned reservations. Two-edge iterator checks avoid temporary
per-step branch lists, and in-place score/key/ordinal sorting preserves stable
ties without hidden sort scratch. Existing u32 support counters and bubble
sums now reject overflow instead of wrapping. Frozen pre-migration algorithms
agree across 72 threading cases and 48 bubble cases; successive allocation
limits test cleanup, and replay tests check refusal followed by retry. Expanded
native comparison and independent auditing pass 768 candidate cases and 12
frozen path-buffer references. They preserve all FASTA/PCR evidence and exact
count/read identities across eight backends, threads 1/2/4/8, both hashers and
enabled/disabled target ceilings, adding paired reverse-mate and unmatched-tail
cases to the exact/IUPAC cap matrix. Independent oracles distinguish raw bases
from valid ingested bases across N gaps and bind actual paired source/mate
identity. **Scoped read-evidence/bubble conformance gate: met 2026-11-13.**
The ordinary sequence-result migration (local evidence: `benchmarks/benchmark_results/v4-result-records-20261113/README.md`)
passes full quality checks (447 Rust tests per hasher and 224 Python tests).
Native comparison and independent auditing pass all 768 candidate cases and
12 frozen evidence-buffer references, preserving complete FASTA/PCR evidence,
count fingerprints and selected-read identities across the same six-scenario
exact/IUPAC/paired matrix. **Scoped ordinary-sequence ownership gate: met
2026-11-13.**
`RetainedResults` is a separate run-level
memory class on the same parent ledger. Its stage ceiling shares the existing
combined run limit rather than increasing it. UTF-8 text reserves replacement
growth before allocation and streams formatting without an intermediate
String, preserving typed reservation errors. Ordinary candidate/record buffers,
sequence text, identifiers and descriptions, and per-record annotation lists
and identifiers retain parent guards through their actual destruction. Stable
score/sequence ties use charged ordinal decorations and in-place sorting;
renumbering moves sequence storage and publication borrows FASTA fields.
Tests retain three real PCR outcomes while admitting subsequent targets into
one slot, reject a result-stage limit independently of graph space, and verify
class-preserving copies, UTF-8 retry, unchanged formatted bytes and moved
sequence pointers. The subsequent
diagnostic/warning migration (local evidence: `benchmarks/benchmark_results/v4-warning-payloads-20261113/README.md`)
passes full quality checks (451 Rust tests per hasher and 224 Python tests).
Native comparison and independent auditing pass 1,024 candidate cases and 16
frozen ordinary-result references, preserving complete FASTA/PCR and read/count
identities across eight exact/IUPAC/paired/hidden-diagnostic scenarios.
**Scoped diagnostic/warning ownership gate: met 2026-11-13.**
Retained path text/hashes, marker identities and
hashes, marker/run/path/threshold vectors, warning sequences, descriptions,
source/product buffers and full sequence lookup keys now use parent-run guards.
Hidden diagnostic copies are fallible and independently owned. Fixed schema
labels and the three cause labels borrow static text, avoiding repeated String
allocations. Canonical standard-equivalence keys use a borrowed reverse-complement
iterator followed by owned bytes. Warning lookup retains full oriented sequence
keys; lexical sorting restores the previous BTreeMap numbering. Digests remain
provenance rather than sequence identity. Typed refusal propagates through DFS
and finalization instead of returning a partial warning collection.
Tests check independent-copy refusal/cleanup, successive allocation limits,
both hidden modes, canonical byte/hash equivalence and deliberately identical
digests on distinct sequences. Existing retention limits and selection policy
are unchanged. The subsequent
stats metadata migration (local evidence: `benchmarks/benchmark_results/v4-stats-metadata-20261113/README.md`)
passes full quality checks (457 Rust tests per hasher and 224 Python tests).
Native comparison and independent auditing pass 1,024 candidate cases and 16
frozen warning-payload references, preserving complete FASTA/PCR, read/count
identities, run/input metadata and recorded commands across the same
eight-scenario matrix. **Scoped stats-metadata ownership gate: met 2026-11-13.**
Target/run/input strings and vectors, failure
reasons, filenames, product-length arrays and digest summaries retain parent
result guards. Fixed statuses, version/layout and source-kind labels borrow
static text. Failure summaries preserve their previous cause/threshold order
without raw string joins. Indexed parallel collection reserves its exact
backing array before callbacks and retains that charge through element and
storage destruction; panic cleanup also releases partially initialized results.
Stats/result lists continue owning their data through streaming YAML output.
The [pinned YAML workspace model](YAML_WORKSPACE.md) now also reserves emitter
buffers, queued scalar copies, default directives, boxed owners, depth stacks
and replacement overlap before the staged stats file is created. An
allocation-free borrowed Serde preflight computes scalar/depth bounds; the
prepared owner retains the same immutable input and guard through emitter
destruction. Unsupported tagged/formatting forms fail preflight.
Full quality checks (local evidence: `benchmarks/benchmark_results/v4-yaml-workspace-20261113/README.md`)
pass 464 Rust tests per hasher and 224 Python tests, including 80 isolated
allocator observations and unchanged YAML bytes. All 1,024 native candidates
and 16 frozen reference runs agree on complete FASTA/PCR, counts and metadata;
independent count and input-frame auditing also passes. **Scoped YAML workspace
ownership gate: met 2026-11-13.**

Deduplication calls pinned `editdistancek 1.0.2` with the same clamp and distance
calculation previously supplied by bio's adapter. After an allocation-free
length-difference rejection, its entire heap workspace is two `Vec<isize>`
arrays of `2*k+3` elements each, where `k` is the threshold clamped to the longer
input. The prefix-matching routines only borrow the input. Reserve both arrays
before the call and retain the guard until they drop, under the result owner's
existing parent ledger and class. Checked arithmetic includes each array's
addressable layout limit. A fallible comparison loop propagates refusal while
dropping kept and unconsumed records; sorting, early duplicate acceptance and
final truncation preserve their prior order.
Full quality verification (local evidence: `benchmarks/benchmark_results/v4-dedup-workspace-20261113/README.md`)
passes 470 Rust tests per hasher and 224 Python tests. Independent full-matrix
and previous-adapter comparisons agree; 48 isolated allocator observations
exactly match the reserved bytes and release both arrays on return. All 1,024
native candidates and 16 frozen reference runs match complete FASTA/PCR, counts
and stable metadata; independent input-frame auditing also passes. **Scoped
deduplication workspace ownership gate: met 2026-11-13.** The former remaining
owner inventory is superseded by the [v4 resource boundary](#v4-resource-boundary).
Library/allocator/filesystem internals use measured headroom and existing
platform safeguards; owned publication data and accepted-output integrity
retain their existing checks.

The workspace namespace checkpoint (local evidence: `benchmarks/benchmark_results/v4-workspace-identity-20261114/README.md`)
validates parent storage before accepting even a missing registered child. Live
creation, cleanup planning and final cleanup repeat verification. NFS binds an
opaque kernel object handle rather than endpoint addresses or client-local
mount/device IDs; local storage also binds boot/device/inode and available birth
time. The fixed-buffer libc adapter isolates the necessary syscall missing
from rustix 1.1.4, with checked C layout, capacity and returned lengths. Path
conversion is preadmitted through verification. libc was already locked; no
locked package version or checksum changed.

The original missing-parent bug is reproduced and fixed. The first address-based
candidate passed full/I/O checks but failed a real shared-NFS retry across two
hosts with alternate server endpoints. Its source and preserved failure remain
in the evidence. The kernel-object correction passes 595 Rust tests per hasher,
225 Python tests, 168 hard-limit component tests and 544 native I/O/retry cases
(1,360 invocations). The corrected two-node matrix and independent audit pass
32 cases, 96 application invocations and 32 actual SIGTERMs. Shared recovery
succeeds remotely, a different local namespace preserves output/registration,
and the live origin reclaims the original private bytes. Kernel probes
independently reproduce the stored NFS fingerprint on both hosts.

Retention verifies five exact application archives, 125 live joins, 94 audited
dependency files and their licensed archive, original helper bindings and all
28,664 payloads in four raw archives. Isolated pathname-allocation observations
fit their precharges and release completely. **Scoped workspace-parent namespace
identity and cross-host recovery gate: met 2026-11-14.** Unsupported identity
queries preserve registration/data; Mac/Windows native checks are not claimed.

The markerless development-journal checkpoint (local evidence: `benchmarks/benchmark_results/v4-legacy-journal-20261114/README.md`)
covers schema-1 journals from earlier unreleased v4 builds, before preparation
controls existed. Released v3.2.1 did not emit recovery journals. Cleanup now
preplans both marker paths and admits all existing payload/journal extents
before mutation. While the journal still grants the validated scope, it creates
and syncs a zero-byte owner, then seals the accepted public manifest. Only then
may it remove the journal, payloads and final owner. Pending zero-byte paths are
not mistaken for existing inventory entries and do not raise the original disk
ceiling. Producer-marked hexadecimal run IDs can be reclaimed; unmarked
nonstandard directories remain untouched.

The original native SIGTERM regression leaves four orphans; the candidate
pilot leaves none. Full checks pass 597 Rust tests per hasher and 225 Python
tests. Independent audits pass 346 component tests under 64 descriptors/512 MiB
and 1,088 native fault cases under 64 descriptors/3 GiB, covering both hashers,
threads 1/2/4/8, rollback and accepted replacement, EIO/ENOSPC/EINTR and actual
SIGTERM. The 2,184 invocations include 276 actual SIGTERMs. Accepted
sequence/count/input evidence and foreign scopes survive; original-ceiling disk
refusal occurs before creating grants and retries at the exact allowance.
Retention verifies three source archives, 126 current-source joins, 94
dependency files and all 30,344 payloads plus 10,301 directories in four raw
archives. **Scoped markerless development-journal interruption gate: met
2026-11-14.** Broader physical/descriptor ownership, other library/preprocessing
resources and final-source promotion remain open. Process interruption does not
establish whole-system power-loss behavior.

The panel configuration snapshot (local evidence: `benchmarks/benchmark_results/v4-panel-preprocessing-20261114/README.md`)
admits raw local/URL input and the audited loader before decoding existing
Serde models. Visitors admit strings, sequences/maps and alias expansions;
derived labels, inline settings, uniqueness tables and validation reports use
the original run ledger. Discarded root metadata and intermediate arrays drop
before a capacity inventory reduces the owner to retained primer data.
Standalone Linux previews observe installed/inherited process limits while
ignoring unused counting controls. Original URL limits/timeouts and strict
UTF-8 remain. The shared decoder-error follow-up (local evidence: `benchmarks/benchmark_results/v4-yaml-error-lifetime-20261114/README.md`)
transfers the loader allowance with heap-owning library errors, dropping the
payload before its guard. Typed admission failures instead discard the temporary
Serde error while the guard lives and return their original downcastable type.

Both hashers pass 609 Rust tests and 225 Python tests. Independent checks pass
428 component executions under 64 descriptors/512 MiB and 26 native cases/48
invocations, including complete FASTA/annotation/count/input parity at threads
1/2/4/8, inherited-preview refusal and paired malformed panel/manifest/journal
diagnostics and preservation. Seven baseline plus two follow-up source archives
and their current-source/compiler/dependency joins pass retention, along with
the sealed baseline and all new raw captures. **Scoped configuration
input/model/preview and YAML decoder-error lifetime gates: met 2026-11-14.**
Raw CLI/clap and other validation/library/error/transport owners, complete
physical/descriptor accounting, public controls and final-source promotion
remain required. These component gates do not establish a whole-run guarantee.

Retained results cannot keep the target admission lease: parallel result
collection holds earlier outcomes while later targets wait for admission, so a
single-slot pool would deadlock. Ordinary records now use separate parent
reservations; diagnostic and warning buffers now follow the same separation.
Sequence RAM can release after its data is staged; annotation owners continue
through stats serialization and cleanup. Streaming alone does not bound
serializer workspace; the separate source-audited reservation covers those
internal allocations. The YAML gate does not cover other libraries or process
headroom.
Configured limits and logical allocation charges still require measured
allocator/process headroom before claiming a whole-run RAM guarantee.

Pre-publication rerun preservation (local evidence: `benchmarks/benchmark_results/v4-rerun-preservation-20261113/README.md`)
passes 404 Rust tests per hasher and 224 Python tests. Completed FASTA/stats and
manifest bytes survive an actual graph-budget failure during a rerun, with
failed-attempt metadata kept separately. This checkpoint does not cover
rollback or recovery during publication. The durable recovery implementation (local evidence: `benchmarks/benchmark_results/v4-publication-recovery-20261113/README.md`)
passes full quality checks (409 Rust tests per hasher, 224 Python tests). It
copies and verifies the prior generation, publishes a synced recovery receipt
before public changes, rolls back incomplete publication and accepts a fully
published complete manifest on restart. Conflicting edits preserve backups and
stop recovery. Native kernel-call fault/interruption auditing and its independent audit now
pass all 160 cases (480 process invocations) under both hashers and threads
1/2/4/8. They verify EIO/ENOSPC, SIGTERM at each publication position, exact
manifest/output restoration and accepted replacements after completion.
**Scoped publication failure/restart gate: met 2026-11-13.** Pre-journal orphan
reclamation and resource accounting for output backups remain open.

The [identity/extraction verification](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-identity-contracts-20261109/README.md) establishes the first item. Store/replay and
budget integration remain open until exercised by their actual implementations.

[Exact-store/ledger verification](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-exact-store-contracts-20261110/README.md)
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

[Assembly-access verification](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-assembly-access-20261110/README.md)
records exact/legacy product and repeat-warning agreement, large-count PCR and
injected scan/lookup failures. This does not close the remaining frontend,
replay, persisted-integrity or combined-resource gates.

## External counting and lookup under #143–#144

`--count-backend external` is a hidden development path, requires `--chunks 0`,
and preserves its independent developer controls. It uses the same bounded
FASTQ reader, rolling windows, selection/cap/paired rules and selected-record
digest as the memory candidate. No original-source replay is needed for counting
or ordinary PCR. Requested read threading now captures original-pass replay in
the exact candidates as described below. Legacy ingestion still uses its
existing source-reopen/read-vector adapter pending production dispatch.

`ExternalCountBuilder` reserves one reusable sort buffer before allocation.
It sorts and reduces positive u64 key/count pairs in place, writes a verified
sorted run, and compacts equal-level runs with a streaming merge. A fixed vector
of 64 run slots bounds retained descriptors and metadata for a u64-sized stream.
Each run owns one charged page; at most 66 run/output pages are allowed by the
current cache ceiling. Inputs remain live until the merged output is complete
and verified. The final reduction is streamed into a source-bound store; it
never becomes a second global RAM table. Invalid identity/overflow preserve the
accepted prefix. Spill I/O/resource failure invalidates the producer so it
cannot publish a partial selected dataset.

The development controls are `--count-sort-buffer-mib` (64 by default),
`--count-disk-budget-mib` (65,536 logical MiB) and `--count-temp-dir` (OS temporary
directory by default). A private owned run directory contains derived count
files; confirmed file and directory cleanup precedes output publication.
These controls cover owned count buffers/extents, not whole-process RAM or
filesystem allocation overhead. They are not promoted production defaults.

### Native sorted-file format 1

| Region | Layout and validation |
| --- | --- |
| 64-byte header | `SHKCTP01`, little-endian format version 1, k, count width 64, canonicalization/preprocessing codes, 254 entries/page, page size 4,096, a fresh 128-bit file identity, and zero reserved bytes. |
| Each 4,096-byte page | 254 little-endian `(u64 canonical key, u64 positive count)` entries and a 32-byte SHA-256 binding the header, page index and payload. The final page has zero padding. File identity rejects foreign-page substitution. |
| 96-byte completion footer | `SHKEND01`, distinct-key count, occurrence total, page count, selected-record SHA-256 and a SHA-256 binding the complete header/pages/footer metadata. |

Opening validates the exact extent, format, full checksum, canonical keys,
strict ordering/uniqueness, nonzero counts, padding and checked totals before
exposing an immutable store. Every later page read verifies its checksum.
Known access/integrity failure invalidates later access; it never becomes key
absence. Intermediate run footers complete their own files with an unbound
zero input digest; only the final whole-input reduction is supplied to PCR,
with the actual selected-record digest and checked total.

`DiskCountStore` binary-searches full keys through one shared 4 KiB lookup cache.
Each active visit uses a separate charged 4 KiB scan buffer and invokes callbacks
outside the file mutex. This permits reentrant lookup and prevents concurrent
scans from evicting one another's page on every entry. No fingerprint is used
as membership evidence. Batch lookup preserves order/duplicates and configured
capacity checks. The initial implementation is an ephemeral owned store; it
does not claim cross-run reuse, restart or resume.

External runs emit `count_store.external` measurements only when performed:
count-file bytes, permanent lookup-cache bytes, peak charged count-page bytes,
peak sort-buffer and logical temporary-disk bytes, sorted-run count and merge
count. Existing table capacity/owned-table bytes are zero because this backend
has no global hash table. Logical reservations remain distinct from allocator
heap, process RSS and filesystem allocation. Production selection, whole-run
resource ownership and complete quantitative gates remain open in the plan.

## Direct-ingestion comparison under #140

`--count-backend memory` is a hidden development control until the frozen
performance gates have been evaluated. It requires `--chunks 0`; explicit
incremental counting stays on the released implementation. The direct path
owns one logical builder, stages rolling canonical windows in a bounded reusable
buffer without per-read allocation, and finalizes by moving its table allocation
into the borrowed assembly path.
Valid-base totals are gathered during rolling extraction, avoiding another
sequence scan or a post-ingestion global count scan.

The serial emission buffer reserves 1,024 typed keys before allocation and
reuses that capacity across records (16 KiB on tested 64-bit targets). Encoding
and table updates run in separate loops. Full blocks and each record's final
partial block are flushed; a failed flush retains the first unaccepted key's
cursor so an explicit retry cannot double-count its accepted prefix. Finalization
rejects unfinished batches. The shared ledger charges this allocation under
`WorkerQueues` in addition to table and input arenas, reports
`peak_kmer_buffer_bytes`, and includes it in `peak_count_input_reserved_bytes`.
This is owned-buffer accounting; the parallel candidate below extends it to
worker key buffers. Whole-process resource enforcement remains open.

## Original-pass replay under #145

The memory, parallel and external development paths capture replay during
ingestion when `--read-threading` and PCR are requested. They do not capture it
for ordinary counting/PCR. One-shot stdin no longer needs to fall back to
k-mer-only scoring on these paths. Replay stores the same selected subset as
counting, including rounded paired caps and the first unmatched record selected
at uneven paired EOF. The original global record index defines fragment identity
as `index / 2` for paired input; unpaired fragment identity is the record index.
This preserves the released pairing policy rather than inferring pairing from
header text.

`ReplayWriter` stages one bounded frame at a time in reusable charged storage.
Its footer binds the actual selected-record digest calculated by counting;
finalization rejects any mismatch and validates the complete persisted stream.
Each `ReplayStore` visit owns only its current charged frame and performs
callbacks outside the file lock. Threading maps records directly from this
visitor, retaining one pending fragment rather than a vector of all reads or
a map of all fragment evidence. Existing graph annotations and paired-link
storage still need their separate evidence/resource integration.

### Native replay format 1

| Region | Layout and validation |
| --- | --- |
| 64-byte header | `SHKRPL01`, little-endian version 1, pairing policy, original-field encoding code, maximum frame bytes, effective selected-read cap, fresh 128-bit file identity and zero reserved bytes. |
| Each variable frame | 64-byte prefix: record index, fragment index, mate, reserved zeros and five u64 lengths. Payload: source, complete original FASTQ header, sequence, separator and quality, followed by a SHA-256 binding the file header and complete frame. |
| 96-byte completion footer | `SHKREN01`, actual record/base totals, completed frame extent, selected-record SHA-256 and complete-stream SHA-256. |

Frame lengths are checked against the configured limit and completed extent
before allocation. Visits check record order, fragment/mate policy, UTF-8,
header/separator structure and equal sequence/quality lengths. They verify
each frame before exposing it and reconcile totals, selected identity and the
completed-stream checksum at EOF. N gaps, original sequence positions, full
identifiers, source identities, quality and orientation survive unchanged.
The format preserves quality bytes without assuming an unrecorded quality-score
encoding. Corruption, incompatible formats and incomplete frames are errors.

Replay buffers use the same `ResourceBudget` as counts, input and external
storage. The current development stage ceiling is 64 MiB; frame limits derive
from the FASTQ record limit plus bounded source/framing space. Logical replay
file extents share `--count-disk-budget-mib` with count files, and owned replay
directories use `--count-temp-dir`. Explicit replay removal precedes publication.
Optional `read_replay` stats record actual format, selected records/bases/digest,
file bytes, frame ceiling and peak owned replay buffers. The field is absent
when this replay path was not used; it does not establish per-sequence read
support and adds no FASTA read-support annotation.

Local integration checks (local evidence: `benchmarks/benchmark_results/v4-replay-20261112/README.md`)
pass 348 Rust tests per hasher and 224 Python tests, including format/integrity,
selected-pair and one-shot threaded PCR. The broader native matrix is in progress.
The explicit legacy backend retains its read vector. Subsequent normal-run
dispatch uses bounded replay; integrated graph/process acceptance and broader
measurements remain open under the current resource policy.

## Packed-table evaluation under #142

The hidden `--count-backend packed` candidate requires `--chunks 0` and shares
the exact reader, original-pass replay, assembly access and resource ledger.
It stores complete canonical keys in `ceil(2*k/8)` bytes per bucket and exact
u64 counts in a separate column. Count zero identifies an empty bucket; key
zero remains valid, and positive abundance is never compressed or saturated.
Linear probing compares the entire key, including when hashes collide.

Seven trailing key-column bytes permit safe eight-byte reads at every bucket,
followed by a width mask. The implementation uses no unsafe loads. Load stays
below 7/8; growth reserves both new columns while old columns remain charged.
Rejected identity, arithmetic and allocation/resource operations preserve
accepted counts. Finalization moves the existing columns into a borrowed exact
store. Concrete lookup/scan dispatch avoids building a second global table.

Optional `count_store.packed` stats describe actual key/count bytes per bucket,
bucket count and padding. Owned table bytes and resize peaks remain distinct
from allocator heap and process RSS. This storage experiment does not change
FASTA annotations, sequence admission or the default backend.

Local integration passes 352 Rust tests per hasher and 224 Python tests.
All supported internal k values, full-key collision/wraparound behavior,
abundance above u32, wrong-k/batch/overflow failures, resize failure/retry and
native PCR/replay equality are covered. A separate synthetic experiment compares
insertion, finalization, mixed present/absent lookup and table storage with
Hashbrown under both hashers. Early measurements exchange smaller storage for
slower operations. The completed real-data paired protocol independently
verifies all 100 invocations, exact counts and 287 records, but its 1.3477940354
geometric-mean wall ratio and all ten input wall allowances fail. All RSS
allowances pass. The packed candidate is not selected for the default.

## Automatic memory-to-disk transition under #143

The hidden `--count-backend adaptive` candidate starts with the exact Hashbrown
builder. Typed count-slot or combined-memory exhaustion during table growth
initiates a one-way transfer. Invalid identity, count overflow, physical allocation
failure and I/O errors do not select a fallback. The rejected update has not
changed the accepted table. Its entries are streamed with their full u64 counts
into the existing bounded external producer, verified against the accepted
occurrence total, and the memory table is freed before retrying that update.
No source is reopened, and no second global count table is built. Original-pass
replay and selected-input hashing continue through the same ingestion path.

Count, sort, run-slot and page reservations share one ledger throughout the
transition, including the time the original table and external buffers coexist.
External buffers are allocated only when spilling. A fitting run keeps its
memory store and removes the unused owned temporary directory explicitly.
A spilled run reports its actual disk backend and retains the earlier table
peak in `count_store.peak_table_bytes`. Temporary count and replay files share
the disk limit. Any failed transfer poisons finalization, cleans owned resources
and cannot publish a partial selected input.

Local spill conformance (local evidence: `benchmarks/benchmark_results/v4-adaptive-spill-20261112/README.md`)
passes 361 Rust tests per hasher and 224 Python tests, including native one-shot
spill/PCR/replay agreement and combined-run pressure with another stage's
reservation held live. The independent input-form audit passes all 1,008
invocations, including 224 actual spills. This serial candidate does not recover asynchronous worker state and does not change
the default selector, whole-run budget status or quantitative promotion gates.

## Parallel automatic spilling under #141 and #143

The hidden `--count-backend adaptive-parallel` candidate uses the same fixed
owner queues as the memory-only comparison. Each owner starts in memory and
can transfer its own accepted weighted prefix to the serial adaptive producer.
A spill runs inside that owner; queued batches retain their buffers and the
coordinator observes bounded backpressure. Memory-only and adaptive owners use
statically selected builder types, keeping the existing memory comparison free
of virtual calls in its counting loop. One requested thread uses serial adaptive
counting.

All owners share the count-table and disk budgets. The one configured sort
buffer allowance is divided into positive pair-aligned shares; it is not granted
in full to every worker. Fixed run descriptors, retained pages and key pools
reserve their aggregate ownership. Finalization joins every owner and verifies
the checked total before exposing a disjoint memory/disk shard store. A lookup
routes to one owner and verifies complete keys; a scan visits each shard once.
Every shard binds the same complete selected-input digest. Explicit cleanup tries
every disk close before returning an error, then removes the owned directory
before output commit. Original-pass replay shares the existing disk ledger.

`count_store.adaptive_parallel` reports actual memory/disk shard counts, final
count-file/cache bytes, sort/merge work and shared sort/disk peaks. Existing
backends omit this field. `count_store.parallel` retains queue and owner
measurements; FASTA annotations and sequence admission are unchanged. Shared
resize scheduling can change which owners spill, while counts, source identity
and configured limits must remain invariant.

Parallel-spill local checks (local evidence: `benchmarks/benchmark_results/v4-parallel-spill-20261112/README.md`)
pass 365 Rust tests per hasher and 224 Python tests, and the native input matrix
passes independent auditing of all 1,008 invocations, including 214 actual
spills and reconstructed storage composition. The hard 256 MiB parallel trial
failed on allocation with 47,132 KiB peak RSS. Deferred assembly-pool startup
passes full 366-Rust/224-Python checks and an actual thread-lifetime probe.
The unchanged-cap retry and independent audit (local evidence: `benchmarks/benchmark_results/v4-deferred-pool-20261112/bounded/summary.json`)
complete one million stdin reads with all three owners spilling, exact counts,
seven byte-identical FASTA outputs and complete reference PCR evidence under
the observed hard 256 MiB address-space limit. Native peak RSS is 55,332 KiB;
aggregate count files are 1,190,621,664 bytes and replay is 428,777,952 bytes.
This establishes the scoped parallel beyond-budget check. Whole-run production
accounting and the original performance promotion gate remain open. Default dispatch, complete process/graph/
evidence accounting and the original performance promotion gates remain open.

## Ordered input-identity pipeline under #140 and #141

Selected-input hashing now uses `ring` 0.17.14, already locked for HTTPS, behind
a fallible cloneable SHA-256 wrapper. The domain, framing and digest bytes are
unchanged. A checked byte census rejects streams exceeding SHA-256's 64-bit
bit-length representation before calling the provider; rejection permanently
poisons that state, so finalization cannot expose a successful prefix digest.
Count fingerprints and persisted count/replay checksums retain their existing
implementation and bytes. The
isolated native provider trial (local evidence: `benchmarks/benchmark_results/v4-native-sha-20261113/README.md`)
passes 390 Rust tests per hasher and 49 independently audited Human invocations,
including complete input/count/FASTA/PCR agreement. Its one-input throughput
improvement does not establish full default promotion. The
combined worktree checks (local evidence: `benchmarks/benchmark_results/v4-selected-sha-20261113/README.md`)
also pass 390 Rust tests per hasher and 224 Python tests. Both exported binaries
match the native provider pilot exactly. Native phase-failure auditing and all
1,512 input-matrix invocations pass. The actual-default-handoff million-read
replay/PCR trial and independent audit pass under a hard 256 MiB address-space
cap, with 31,824 KiB peak RSS and exact full reference agreement. The full
performance comparison passes all 100 correctness checks and RSS allowances
but misses aggregate speed (0.9632352084 versus 0.95) and Human wall time
(0.67 versus 0.66 seconds allowed). Its conditional matrix is explicitly skipped;
default promotion remains not met.

The hidden `--count-backend pipeline` candidate retains the exact
`sharkmer-selected-fastq-v1` byte framing and SHA-256 digest. The coordinator
serializes record index, mate, source and all four original FASTQ fields to
one FIFO hashing worker. Three reusable 64 KiB buffers stream arbitrarily long
accepted fields without a whole-record clone. The shared `WorkerQueues`
reservation covers their storage, producer descriptor and fixed channel
allowance before allocation or worker startup; it does not account for OS
thread stacks or complete process overhead.

At one requested thread, hashing stays inline with serial adaptive counting.
At two, the coordinator counts while one worker hashes. Higher requests use
`threads - 2` adaptive count owners, the identity worker and coordinator.
This shares the requested ingestion thread allowance rather than adding an
extra worker to an existing owner pool. Assembly workers start after ingestion
workers join, as for the other exact candidates.

Record indices must be contiguous. Checked frame lengths and byte totals are
validated before appending a record. Partial append or worker failures poison
the producer; finalization closes channels, drains the tail, joins and verifies
the consumed byte total before exposing a digest. The digest seals further
ingestion and is reused for count-store and replay finalization. Drop joins
before releasing buffer reservations. `count_store.input_identity_worker`
reports the selected-record count, framed bytes and fixed pool ownership only
when background hashing actually ran; inline paths omit it.

Local checks and the independent native audit (local evidence: `benchmarks/benchmark_results/v4-identity-pipeline-20261112/README.md`)
pass both hashers, malformed/failed-worker handling, actual thread lifetimes,
and all 1,008 input-form cases including 223 spills. The independent oracle
reconstructs the complete frame and count pairs from original FASTQ bytes.
Expanded full checks pass 375 Rust tests per hasher and 224 Python tests,
including 48 native failure invocations per hasher across adaptive backends
and threads 1/2/4/8. Later gzip-member checksum/trailer failures, invalid bases
after an accepted prefix and configured disk exhaustion produce failed
manifests without sequence/stats publication; original/foreign files survive
and owned temporary files are removed. Physical disk-full and interruption
matrices remain separate requirements.
The native storage-fault matrix (local evidence: `benchmarks/benchmark_results/v4-storage-faults-20261112/README.md`)
also verifies controlled partial-write/ENOSPC and sync/EIO propagation, failed
publication and fresh restart after SIGTERM. Interrupted count/replay files
remain uniquely named and unchanged; automatic reclamation/resume is not
implemented. The Human attribution pilot improves on the earlier parallel
producer but does not establish default promotion. The full unchanged
comparison preserves counts, 287 records, source identities and RSS allowances,
but fails promotion with a 1.0139565284 wall ratio and the Human per-input
allowance. Its conditional matrix is skipped because that prerequisite failed.
This candidate's own hard-limit replay/PCR trial failed
on allocation in assembly; a coordinator-participation correction and
unchanged-cap retry follow. The correction now passes 376 Rust tests per
hasher and 224 Python tests, including actual total OS-thread probes at
1/2/4/8. The coordinator occupies one Rayon participant slot, leaving
`threads - 1` background assembly workers and preserving the requested number
of computing participants. Its unchanged-cap retry and independent audit now
pass: one million selected records, exact count/replay identity, seven
byte-identical FASTA outputs and complete reference PCR evidence, with 12,828
KiB peak RSS and at most four observed OS threads. The subsequent descriptor
revision also passes its own unchanged-cap trial and audit: all seven FASTA
outputs/full PCR agree, with 18,876 KiB peak RSS and at most 21 sampled process
descriptors. Those samples are observations, not a true FD peak. Whole-run
production resource accounting remains open.

## Hybrid worker roles under #140 and #141

The hidden `--count-backend hybrid` candidate initially follows the ordered
pipeline. At two threads, a successful complete record whose current count
table allocation reaches 16 MiB triggers a change in ownership. The coordinator
flushes its key tail, drains and joins the identity worker, recovers its
unfinished SHA-256 state, and moves the existing `AdaptiveCountBuilder` into
one counting owner. It continues the same hash stream inline while the owner
accepts subsequent keys. No global table is copied, no source is reopened and
no prefix digest is mistaken for the completed input identity.

The owner inherits k, the same resource ledger and the exact accepted occurrence
total from its prefix. Startup/finalization errors cannot publish a subset.
The old identity buffers release their reservations before the counting worker
starts; both bounded queue configurations are allowed for in the initial ledger.
There are at most two OS threads through the handoff. Threads 1 and requests
above two retain the corresponding pipeline dispatch. The default threshold
is fixed before performance evaluation; the hidden threshold override exists
to exercise the transition in small conformance fixtures.

`count_store.count_owner_transfer` is emitted only after an actual successful
handoff. `selected_records_before_transfer` and `occurrences_before_transfer`
describe the completed prefix; `current_table_bytes` is its allocation at the
handoff, and `threshold_bytes` is the trigger. When this field exists,
`input_identity_worker` describes only that worker's prefix work. The final
`selected_input_sha256` still covers every selected record, including the
suffix hashed inline. A run without transfer omits the phase field.

Full hybrid conformance and independent native auditing (local evidence: `benchmarks/benchmark_results/v4-hybrid-owners-20261112/README.md`)
pass 388 Rust tests per hasher, 224 Python tests and 1,512 native invocations,
including 56 actual handoffs and 202 spills. The first native harness's invalid
no-spill resize prediction is preserved with its corrected verification; actual
counts, identity and resource ceilings remain required throughout. Native
phase faults also verify owner-startup EAGAIN after identity joins and partial
ENOSPC/sync EIO after the new owner starts. The
hard-limit trial reaching the actual default handoff (local evidence: `benchmarks/benchmark_results/v4-hybrid-owners-20261112/bounded/summary.json`)
passes at two threads under a 256 MiB address-space cap, with 32,028 KiB maximum
RSS and exact full input/count/PCR/FASTA agreement. The identity worker joins
after 3,889 complete records; one count owner subsequently spills. Both the
1,190,613,152-byte final count file and 428,777,952-byte replay exceed the
address-space cap. The original full paired comparison passes all 100
correctness checks and per-input wall/RSS allowances but misses aggregate speed
(geometric-mean wall ratio 0.9656945923 versus the required 0.95). Its conditional
matrix is skipped explicitly. No default promotion follows from these
correctness and hard-limit proofs; the newer provider comparison is separate.

## Parallel ingestion under #141

`--count-backend parallel` is a hidden development selector requiring
`--chunks 0`. With one requested thread it runs the serial memory candidate.
Otherwise one coordinator parses/encodes input and routes canonical keys to
`threads - 1` counting owners. A deterministic mix assigns each key to exactly
one owner. Owners retain disjoint tables, sharing the count-slot reservation
ledger; finalization joins owners and moves their tables into a sharded immutable
store without consolidation or duplicate global maps. Lookup routes to one
owner's completed table, while primer scans visit every shard.

Each owner has one queue slot and three preallocated buffers of 1,024 canonical
keys. The compact-queue development revision stores eight-byte encodings after
the coordinator checks the typed key's k, and reconstructs trusted typed keys
inside the owner. Earlier retained comparisons used sixteen-byte key objects.
The coordinator moves full buffers through the queue and reuses cleared ones.
Partial tails are drained before finalization. The shared ledger charges the
actual key storage and shard/worker descriptors, with a conservative fixed
channel allowance. A development ceiling prevents thread requests from making
the key pool unbounded; OS thread stacks and allocator/process overhead still
need production accounting.

Channel failure joins the affected owner and preserves its typed underlying
error. Finalization and producer destruction close channels and join every
started owner before releasing key-pool reservations. No failed owner can
publish a partial selected-input count store. `count_store.parallel` records
worker, batch, queue and reusable-buffer measurements only when background
counting actually ran; the single-thread fallback omits it.

The candidate passes both-hasher local conformance and clean/warning one-shot
PCR at threads 1/2/4/8. Native input-form coverage, the unchanged quantitative
promotion gates and whole-run resource integration remain open in the plan.

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

The hidden `--count-table-budget-mib` bounds owned count storage, including resize
peaks. `--fastq-record-bytes` bounds each complete FASTQ arena. These controls
are not whole-process RAM limits: allocator overhead, decompression/transport
buffers, graphs and evidence still require migration. Exact runs report count
precision, input identity and owned table/record-buffer measurements under
`count_store`. `peak_count_input_reserved_bytes` is logical shared reservation
accounting; use allocator/RSS metrics separately. Default/legacy runs do not
emit an unperformed exact-store assessment.

The hidden `--count-fingerprint` performs an opt-in verification scan over all
full key/count pairs, with k included in each pair hash. Its versioned SHA-256
XOR digest, distinct cardinality and checked occurrence total support
order-independent backend comparisons. This scan is excluded from normal
operation and timed promotion runs. The paired comparison harness checks it
separately before timing each input against the immutable v3.2.1 binary.

The [independent paired audit](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-direct-comparison-20261112/README.md)
verifies all 100 revised-candidate invocations, but rejects the custom table on
the unchanged performance gates. The [compatibility refinement](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-direct-refinement-20261112/README.md)
preserves the released separation of primer and graph thresholds: primer
scans use the target's minimum on unfiltered entries; graph lookups use the
global view threshold. A production table must pass these same semantics and
the frozen quantitative gates before default selection.

The replacement candidate uses Hashbrown 0.16.1 with the selected ahash/Fx
hasher and full u64 keys/counts. The exact dependency pin protects the audited
pre-allocation bound: `(u64,u64)` buckets require 16 bytes plus a control byte
per bucket and at most sixteen trailing control bytes. Capacity starts at
fourteen entries and doubles at the library's 7/8 load factor. Reserve the
old allocation plus the new conservative bound before calling `try_reserve`.
The library reserves replacement storage before moving entries and preserves
accepted contents on allocation failure. Keep both guards live through that
operation, then replace the freed allocation's guard and check the public
`allocation_size()`. This avoids a second key-insertion pass during growth.
An existing key can still be updated at full capacity without reserving growth.
Layout tests exercise the bound over multiple growth sizes; a dependency update
must recheck this layout and the actual resource/performance gates.

`count_store.owned_table_bytes` and `peak_table_bytes` include bucket/control
allocation, unlike allocator overhead or process RSS. Earlier rejected
prototype receipts retain their original `*_slot_bytes` names and meanings.
Concrete memory access and visitors retain their types through assembly so
the inner scan can inline its callback. Opaque stores retain the fallible
trait interface, including injected I/O-error tests. This candidate still
requires fresh measured promotion; the earlier table's correctness evidence
does not establish its performance.

The first Hashbrown comparison (local evidence: `benchmarks/benchmark_results/v4-hashbrown-comparison-20261112/README.md`)
passes independent auditing of all 100 paired invocations and all 80 two-hasher,
1/2/4/8-thread matrix invocations. Full count fingerprints, complete sequence
inventories, warning causes and selected-input digests agree. Its geometric
mean wall ratio of 1.0909445088 fails the unchanged performance requirement,
despite lower RSS. This is a rejected default candidate. Subsequent
native-growth and inline-emission conformance (local evidence: `benchmarks/benchmark_results/v4-native-emission-checks-20261112/README.md`)
is separate evidence and requires its own measured promotion decision.

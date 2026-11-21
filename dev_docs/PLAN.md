# Development plan: scalable, reliable sPCR

**Active version: 4.0.0-dev. Released baseline: [v3.2.1](https://github.com/caseywdunn/sharkmer/releases/tag/v3.2.1) (2026-11-09).**

**Status:** implementation and final local verification are complete; Git/CI
integration remains open. The user narrowed
v4 resource scope on 2026-11-16: bound Sharkmer-owned count, ingestion, lookup,
replay, assembly and retained-output buffers; use established libraries and
process headroom for networking/runtime internals. Do not implement pool, DNS,
socket or TLS accounting, or replace those subsystems to prove a memory bound.

**Previous whole-process/per-allocation networking gate: waived by the user
2026-11-16**, because maintaining client-internal formulas and hooks adds
complexity outside v4's counting/ingestion goal. This is a scope decision, not
proof that every library allocation fits. Ordinary network/HTTP/decoder errors
fail the operation; an OS-enforced limit can terminate the process. Previously
accepted output must remain intact, and failed work must not become complete.
Existing logical disk/descriptor controls, process safeguards, application body
limits and timeout behavior remain. No new cross-platform RSS guarantee is made.

The networking simplification removes the experimental Agent/request/URL/header
allocation owners and restores upstream client code except the demonstrated
chunk-framing fix and a warning-only middleware change. HTTP content decoding
also returns to the client. Application FASTQ decoding/buffering and ENA/panel
body limits retain their existing owners. The
simplification checks (local evidence: `benchmarks/benchmark_results/v4-http-scope-20261116/README.md`)
pass 625 Rust tests per hasher, 225 Python tests, all quality checks, eight native
HTTP invocations and the unchanged cache matrix (544 cases, 1,488 invocations,
212 SIGTERMs) at 64 descriptors/512 MiB. Five client files match upstream;
only chunk framing and warning cleanup remain patched. **Scoped HTTP
simplification gate: met 2026-11-16.** Earlier accounting checkpoints retain
their dated evidence. See [HTTP scope](HTTP_WORKSPACE.md) and the
[resource boundary](COUNT_STORE_CONTRACTS.md#v4-resource-boundary).

**Final local verification: met 2026-11-18.** The existing resource controls
are visible in `--help` and documented; defaults and calculations are unchanged.
The [final integration report](https://github.com/caseywdunn/sharkmer-benchmarks/tree/main/benchmarks/benchmark_results/v4-final-integration-20261117)
records 625 Rust tests per hasher, 225 Python tests, development metadata,
36 default/refusal checks, 1,024 compatibility cases plus 16 references,
160 publication/recovery fault cases, the 256-MiB beyond-budget count/replay/PCR
trial, 1,024-target rollback/backup refusals and 64-GiB retries, the original 100-invocation paired
speed/RSS gates and 80-case thread/hasher matrix. Source hashes join those runs.

**Next:** integrate the issue branch into `dev` after CI passes. No release
ritual has started, and explicit approval is still required. New analysis
results remain outside version control to keep production downloads small.
Exhaustive
dependency allocation proofs remain outside scope. Biological #129 stays in v4.1.

## Repository size and benchmark storage

The [public companion repository](https://github.com/caseywdunn/sharkmer-benchmarks)
holds browseable reports and protocols, including the final v4 verification.
Its [historical archive](https://github.com/caseywdunn/sharkmer-benchmarks/releases/tag/legacy-2026-11)
preserves all 1,538 result files at `7bd1ac3` and a benchmark-only Git bundle
covering both historical results directories. The downloaded archives matched
their checksums. The original Git bundle and unpublished result files are
backed up separately in the development workspace.

The cleanup keeps the benchmark tools in production and gives their unit tests
a small synthetic fixture. Removing the committed results reduces the tracked
checkout from 252 MiB to about 21 MiB. The 32-ref rehearsal removed
`benchmarks/benchmark_results/` and its earlier name, `benchmarks/results/`.
Every other file matched at each branch and tag tip. The Git pack fell from
187 MiB to 32.8 MiB. This measures branch and tag history in the rehearsal.

The companion holds [the original source tarballs](https://github.com/caseywdunn/sharkmer-benchmarks/releases/tag/original-sources-2026-11)
for all seven existing release tags. Before applying the rewrite on GitHub,
protect package source URLs. Bioconda
3.2.1 shipped via [PR #71402](https://github.com/bioconda/bioconda-recipes/pull/71402).
It pins the v3.2.1 generated tarball at SHA-256
`e8031a7b346303f05e553886d473781b01a1a9c0dda0752d5b4573af6ad6061c`.
The current recipe needs to point to a stable copy of those exact bytes
before the v3.2.1 tag can move. Open PR #71407 is a separate recipe update.
Review affected release tags, downstream pins and the ref map together before
cutover. Developers will need to re-clone or transplant their work onto the
rewritten history. The original Git bundle provides a rollback point.

**History cutover gate: not met.** The rehearsal and archive checks are
complete. The published tags and live package source pin still need the
compatibility step before a coordinated remote rewrite.

## Previous resource checkpoints

These dated results preserve their original source and evidence scope. Their
lists of unfinished whole-process/networking owners are superseded by the
user-approved scope above; they do not create new v4 release gates.

The previous
fixed Agent checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-agent-backing-20261116/README.md`)
admits default/panel client containers before construction. Shared accounting
leases cover body/error Weak backing without extending the connection pool's
strong lifetime; panel timeouts and the original ledger are preserved. Full
checks pass **647 Rust tests per hasher and 225 Python tests**. Native HTTP
auditing passes 52 invocations at unchanged 64-descriptor/512-MiB limits and
388 allocation observations, including 16 warm constructor measurements.
The independent 16-pair wire audit and unchanged cache matrix pass: 544 cases,
1,488 invocations and 212 actual SIGTERMs. Retention reopens two exact 197-file
source archives, 789 dependency/compiler bindings, licensed sources and
15,609 raw files plus 2,330 directories. The independent seal verifies 57
files. **Scoped fixed Agent/Weak backing and current HTTP/cache integration
gates: met 2026-11-16.** Pool tables/keys/streams, resolver/socket/global and
connection TLS, remaining process/preprocessing/library/error owners, complete
physical/descriptor accounting, public controls, final-source promotion and
Git/CI remain open. No release ritual has started.

The
derived-request checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-request-buffers-20261116/README.md`)
admits method/header/authentication/serialization storage before Unit creation,
using actual header capacities and configured User-Agent length. A separate
request reservation releases after the consuming call; retained HTTP errors
keep only their still-live URL/response/decoder owners. Full checks pass
**643 Rust tests per hasher and 225 Python tests**. Independent auditing passes
44 native HTTP invocations and 16 direct/budgeted wire pairs at unchanged
64-descriptor/512-MiB limits, preserving 372 existing allocation observations.
Fresh cache auditing passes 544 cases, 1,488 invocations and 212 actual SIGTERMs.
Retention reopens three 195-file source archives, 785 dependency/compiler
bindings, the new 36-member licensed dependency archive and both reused
licensed archives, plus 15,585 raw files and 2,322 directories. The final seal
verifies 61 files. **Scoped derived-GET request and current HTTP/cache
integration gates: met 2026-11-16.** Agent/pool/resolver/socket/TLS, remaining
process/preprocessing/library/error owners, complete physical/descriptor
accounting, public controls, final-source promotion and Git/CI remain open.
No release ritual has started.

The
URL/redirect workspace checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-url-workspace-20261115/README.md`)
admits initial parsing, each redirect join, URL/IDNA workspace, retained clones
and history before allocation. Successful responses release the owner before
body delivery; retained errors keep it through payload destruction, and
admission failures preserve their original typed error. Full checks pass
**639 Rust tests per hasher and 225 Python tests**. Independent HTTP auditing
passes 36 native invocations and 288 URL/four-clone allocation observations
plus 84 other HTTP observations at unchanged 64-descriptor/512-MiB caps. Fresh
cache auditing passes 544 cases, 1,488 invocations and 212 actual SIGTERMs.
Retention verifies three 194-file source archives, 749 dependency/compiler
bindings, both licensed dependency archives and 15,560 raw files plus 2,314
directories. The seal verifies 73 files, preserving the initial source and
the ignored-control/Punycode-prefix correction. **Scoped URL/parser/clone and
current HTTP/cache integration gates: met 2026-11-15.** Derived request
method/header/prelude copies, Agent/pool/resolver/socket/TLS, other process
owners, complete physical/descriptor accounting, public controls, final-source
promotion and Git/CI remain open. No release ritual has started.

The
original HTTP request-text checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-request-text-20261115/README.md`)
admits the raw URL and method copies before either request builder, retaining
their reservation until the original builder drops. Full checks pass **634
Rust tests per hasher and 225 Python tests**. Independent native auditing passes
all 13 HTTP tests separately under both hashers (26 invocations), including
28 exact raw-text allocation observations, at unchanged 64-descriptor/512-MiB
limits. Fresh cache auditing passes 544 cases, 1,488 invocations and 212 actual
SIGTERMs. Retention reopens two 193-file source archives, both raw archives
(15,530 files and 2,304 directories), and the original licensed dependency
archive; the final seal verifies 51 files. **Scoped original request-text and
current HTTP/cache integration gates: met 2026-11-15.** Parsed/redirect URLs,
request headers, Agent/pool/transport, other resource owners, final-source
performance and Git/CI remain open. No release ritual has started.

The
markerless-journal checkpoint (local evidence: `benchmarks/benchmark_results/v4-legacy-journal-20261114/README.md`)
now gives schema-1 journals from earlier unreleased v4 builds a durable owner and
accepted-manifest anchor before journal removal. The original four-case native
SIGTERM regression left four orphans; the corrected pilot leaves none. Full
checks pass 597 Rust tests per hasher and 225 Python tests, and independent
hard-limit auditing passes 346 component tests under 64 descriptors/512 MiB.
The native matrix and independent audit pass 1,088 cases, 2,184 invocations and
276 actual SIGTERMs, preserving accepted generations and foreign scopes.
Retention verifies three source archives, 126 current-source joins, 94
dependency files and all 30,344 payloads plus 10,301 directories in four raw
archives. **Scoped markerless development-journal interruption gate: met
2026-11-14.** Released v3.2.1 did not emit recovery journals. No release ritual
has started.

The panel configuration v7 snapshot (local evidence: `benchmarks/benchmark_results/v4-panel-preprocessing-20261114/README.md`)
now admits raw local/URL input, YAML workspace, typed models, derived labels,
inline settings and validation reports on the original ledger before collection.
Discarded root metadata and arrays release before counting. Linux previews
honor the inherited cap while ignoring unused counting settings. Full checks
pass 608 Rust tests per hasher and 225 Python tests; independent auditing passes
426 component executions at 64 descriptors/512 MiB and 20 native cases/36
invocations. Headers, sequences, structured annotations, counts and input
identity match the pre-change v4 executable at threads 1/2/4/8 under both
hashers; large-panel refusal preserves accepted/foreign files. Seven source
archives, 130 snapshot source files, 245 dependency/compiler bindings and all
314 raw payloads plus 94 directories pass retention. The seal verifies 265
checkpoint files. Original ceilings are unchanged.

The shared decoder-error lifetime follow-up (local evidence: `benchmarks/benchmark_results/v4-yaml-error-lifetime-20261114/README.md`)
reproduces retained YAML errors losing admission under both hashers, then keeps
the loader allowance with library errors in configuration and historical
manifest/journal readers. Original diagnostics and typed resource errors remain
intact. That snapshot's full checks pass 609 Rust tests per hasher and 225 Python tests;
independent audits pass 428 hard-limit component executions and 26 native
cases/48 invocations, including paired malformed panel/manifest/journal behavior.
Two exact source archives, 130 current files, 245 dependency/compiler bindings,
the sealed 265-file baseline and all 298 new raw payloads plus 94 directories
pass retention. **Scoped configuration input/model/preview and YAML decoder-error
lifetime gates: met 2026-11-14.** Raw CLI/clap, other validation/library/error and
transport owners, complete process/filesystem accounting, public controls,
final-source performance and Git/CI remain open. The broader release gate is
not met; no release ritual has started.

The response-metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-metadata-20261114/README.md`)
now admits status/header parsing before connecting and retains its owner with
library errors, releasing before successful body buffering. The
[source-derived allowance](HTTP_WORKSPACE.md) is separate from gzip. Full
checks pass 612 Rust tests per hasher and 225 Python tests; independent auditing
passes 434 component executions at 64 descriptors/512 MiB, plus 20 paired native
cases/40 invocations preserving HTTP behavior, diagnostics and unrelated files.
Retention verifies three source archives, 131 current files, 245 dependency
bindings and all 439 raw payloads plus 90 directories. The final seal verifies
75 evidence files. **Scoped response status/header and error-lifetime gate:
met 2026-11-14.** Request/URL/resolver/socket/TLS ownership and fresh integrated
cache/publication I/O/SIGTERM checks remain open. No release ritual has started.

The native-cache follow-up (local evidence: `benchmarks/benchmark_results/v4-http-cache-integration-20261114/README.md`)
corrects its ignored driver's admission while preserving original non-response,
disk, descriptor and kernel ceilings. Full checks pass 612 Rust tests per hasher
and 225 Python tests; production executable hashes are unchanged. Independent
auditing passes all 544 cases, 1,488 invocations and 212 actual SIGTERMs.
Retention verifies two source archives, 131 current files, 245 dependency
bindings and all 15,454 raw files plus 2,278 directories. The seal verifies
61 evidence files. **Scoped response-owner/cache I/O/SIGTERM integration gate:
met 2026-11-14.** This snapshot precedes the chunk-framing patch below.

The chunk-framing checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-chunk-framing-20261114/README.md`)
fixes the reproduced allocator abort in the locally pinned client. Permanent
parser tests stream 768 MiB with zero heap allocation under both hashers. Full
checks pass 630 Rust tests per hasher and 225 Python tests. Independent native
replay succeeds at about 6.4 MiB RSS, versus the original 262 MiB/abort case,
at unchanged 64-FD/512-MiB caps. Fresh cache auditing passes 544 cases,
1,488 invocations and 212 actual SIGTERMs. Source/upstream/dependency/raw
retention passes; the seal verifies 98 evidence files. See
[HTTP workspace](HTTP_WORKSPACE.md). **Scoped chunk-framing and current-cache
compatibility gate: met 2026-11-14.** Broader resource and final-source promotion
gates remain open. No release ritual has started.

The
workspace namespace checkpoint (local evidence: `benchmarks/benchmark_results/v4-workspace-identity-20261114/README.md`)
now validates original temporary storage before accepting a missing child or
cleaning private data. Kernel NFS object handles support alternate server
addresses; local storage binds boot/device/inode. Full checks pass 595 Rust
tests per hasher and 225 Python tests, all 168 hard-limit component tests and
544 native workspace I/O/retry cases (1,360 invocations). The corrected two-node
trial and independent audit pass 32 cases, 96 invocations and 32 actual SIGTERMs,
with shared recovery, node-local refusal and recovery on the live origin.
Retention verifies five source archives, 125 live joins, 94 dependency files and
all 28,664 payloads in four raw archives, including the original regression and
failed address-based trial. **Scoped workspace-parent namespace identity and
cross-host recovery gate: met 2026-11-14.** Mac/Windows native proof and broader
physical/descriptor ownership remain open.

The automatic HTTP gzip (local evidence: `benchmarks/benchmark_results/v4-http-codec-20261114/README.md`)
and ENA body/model/sample (local evidence: `benchmarks/benchmark_results/v4-ena-20261114/README.md`)
checkpoints also pass their scoped full/native/retention checks.
Parsed/redirect-URL/request-header/TLS and raw CLI/clap preprocessing,
complete descriptor/physical ownership, public resource controls, final-source
performance and Git/CI remain open. Release ritual requires explicit user
approval.

The preceding
cache-recovery checkpoint (local evidence: `benchmarks/benchmark_results/v4-cache-recovery-20261114/README.md`)
adds registered same-filesystem staging, fixed checksummed receipts, retained
backups and repeatable cold restoration. Full checks pass 570 Rust tests per
hasher and 225 Python tests; all 100 hard-limit cache/HTTP component tests pass.
The corrected native matrix and independent audit pass 544 cases/1,488
invocations, including 212 actual SIGTERMs. Retained verification checks three
exact source archives, 116 current-source joins and all 16,196 raw payloads.
**Scoped registered-cache publication/restoration I/O/SIGTERM gate: met
2026-11-14.** Codec/transport workspace, complete process descriptor/physical
filesystem ownership, public resource controls, final-source performance and
Git/CI remain open. No whole-system power-loss claim follows from these trials.

The preceding
cache-data/workspace-I/O checkpoint (local evidence: `benchmarks/benchmark_results/v4-cache-data-20261114/README.md`)
passes 562 Rust tests per hasher, 225 Python tests and all development quality
checks. All 84 cache/HTTP component tests pass under actual 64-descriptor/512 MiB
limits. Independent native auditing passes 544 registered-workspace I/O cases
(1,360 invocations) under actual 64-descriptor/3 GiB limits. Retained verification
checks the exact 113-file source archive and all 13,342 raw payloads.
**Scoped compressed-cache extent and registered-workspace I/O gates: met
2026-11-14.**

The preceding
registered count/replay checkpoint (local evidence: `benchmarks/benchmark_results/v4-registered-workspace-20261113/README.md`)
passes 551 Rust tests per hasher, 225 Python tests and all development quality
checks. Independently audited native trials pass 240 cases/760 invocations,
including 240 real SIGTERMs and 128 strict quota refusals, preserving accepted
sequences and removing registered leftovers. Retained verification checks four
source archives, 8,112 raw payloads and 50 quality/provenance files. **Scoped
registered-workspace SIGTERM/quota gate: met 2026-11-14.** A later receipt guard
passes 81 output/recovery tests and 4 startup tests plus Clippy per hasher. The
cache-metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-cache-metadata-20261114/README.md`)
includes that guard in full checks: 560 Rust tests per hasher, 225 Python tests
and all development quality checks pass. All 80 cache/HTTP component tests
also pass under real 64-descriptor/512 MiB limits and independent auditing.
**Scoped cache-metadata gate: met 2026-11-14.**

The staging preparation checkpoint (local evidence: `benchmarks/benchmark_results/v4-preparation-20261113/README.md`)
passes 541 Rust tests per hasher, 225 Python tests and all development quality
checks. All 400 native cases/1,232 invocations and independent audits pass,
including 528 actual SIGTERMs and the corrected earlier-restore-copy regression.
Retained verification checks five exact source archives, 28,331 raw payloads and
70 sealed quality/provenance files. **Scoped registered-staging SIGTERM cleanup
gate: met 2026-11-13.** Whole-system power loss is outside this signal proof.

The temporary-metadata checkpoint (local evidence: `benchmarks/benchmark_results/v4-metadata-disk-20261113/README.md`)
charges new manifest/journal bytes and retains separate journal ownership
through cleanup. Count-disk peaks now describe counting/finalization rather
than completed startup. All native compatibility, I/O/SIGTERM, quota,
completed-replacement and new metadata-error/measurement audits pass with
the original ceilings. Combined retained verification passes 187 evidence
files in 11 groups, six source archives, 110 current-source joins and all raw
payload checks. **Scoped temporary-metadata and count-disk measurement gate:
met 2026-11-13.**

The corrected recovery-storage checkpoint (local evidence: `benchmarks/benchmark_results/v4-recovery-completion-20261113/README.md`)
now charges interrupted private data and planned restoration copies together,
with in-process credit and Unix hard-link deduplication. Independent native
compatibility, publication/staged/restoration I/O-fault and actual SIGTERM
checks pass. The initial large rollback-quota trial restored accepted bytes but
exposed a private-directory cleanup failure, preserved in its failed artifacts. A
correction preplans every private removal path and also accounts for restored
old public files until manifest acceptance. Its full local, native compatibility
and publication/staged/restoration I/O/SIGTERM checks pass. Large quota,
restored-public and completed-replacement checks also pass their independent
audits with the original ceilings and no private orphans. Combined retained
verification passes 127 evidence files in nine groups, three source archives,
109 current-source joins and every payload in eight new raw archives plus the
reused historical archive. **Scoped receipted-recovery storage and planned
cleanup gate: met 2026-11-13.** The schema-4 annotation migration is
implemented and locally verified under #159.

The frozen normal-default snapshot (local evidence: `benchmarks/benchmark_results/v4-production-defaults-20261113/README.md`)
passes the original 100-run speed/RSS comparison and 80-case thread/hasher
matrix. It predates later ownership changes; final-source performance and
integrated quantitative promotion remain **not met**.

Continue with final-source application-resource and correctness checks,
publication/recovery, public resource controls and performance verification.
Exhaustive pool/resolver/socket/TLS, CLI-library, allocator and filesystem
internal accounting are outside this release's scope. Use measured headroom
and existing platform safeguards; preserve the application limits and owned
disk/descriptor controls. See [HTTP scope](HTTP_WORKSPACE.md).
Panel collection and ENA resolution use the prepared shared ledger; preserve
standalone validation, dry-run/early exits, URL timeouts and input semantics
when extending ownership. Git/CI integration remains open.
Release ritual requires explicit user approval.

This plan tracks unfinished v4.0 and later work. The
[v3.2 release plan and gate history](https://github.com/caseywdunn/sharkmer/blob/v3.2.0/dev_docs/PLAN.md),
[v3.2 review](REVIEW_v3.2.md), and
[final comparison](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v3.2-final-check-20261108/RESULTS.md)
preserve the completed release's evidence and decisions. [The v3.0 plan/session log](https://github.com/caseywdunn/sharkmer/blob/v3.2.1/dev_docs/PLAN_v3.md)
retains the earlier v3 plan; [ROADMAP.md](../ROADMAP.md) describes the release sequence.
Tracking issue: [#152](https://github.com/caseywdunn/sharkmer/issues/152).

The v3.2 release retains all 86 historical v3.1 sequence identities in the
scoped comparison, including 22 now carrying warnings. Its inclusive output
preserves bounded hypotheses without asserting biological correctness. Follow
[the sequence-annotation conventions](SEQUENCE_ANNOTATIONS.md) for the v4.0
interface cleanup under [#159](https://github.com/caseywdunn/sharkmer/issues/159).

Independent held-out registration and the full nine-panel biological campaign
were deferred from v3.2 by the user, not completed. They remain required under
[#129](https://github.com/caseywdunn/sharkmer/issues/129) before further
assembly-policy tuning or broad accuracy claims. Additional emitted records
do not establish more validated loci, haplotypes, or correct repeats.

Check items off only when code, validation, documentation, and acceptance
criteria are complete. Issue bodies carry implementation scope and dependencies.
The user and coding agent assess release evidence; no external collaborator
review is required.

## Maintenance and distribution follow-up

- [x] Bioconda 3.2.1 shipped through [PR #71402](https://github.com/bioconda/bioconda-recipes/pull/71402).
  The [channel](https://anaconda.org/bioconda/sharkmer/files) lists Linux x86-64,
  Linux ARM64 and macOS ARM64 builds. PR #71407 is a duplicate update.

## Release strategy

Ship independently useful releases as the counting, read-evidence, and
metagenomic inference work matures. Version numbers are planning targets,
not promised dates. Opt-in additions can ship as v4.x minor releases; reassess
the major version if an implementation requires incompatible behavior or formats.

| Release | Primary result | Boundary |
| --- | --- | --- |
| v4.0 | Faster exact counting and bounded application buffers for ingestion/lookup | Existing standard sPCR remains the biological baseline; library/runtime headroom is not a whole-process RAM guarantee |
| v4.1 | Read-supported low-coverage nuclear recovery | Explicit partial/ambiguous outcomes; optional local multi-k remains experimental until validated |
| v4.2 | Recovery of supported metagenomic template diversity | Rare-template retention, phase uncertainty, and defensible abundance reporting |
| v5.0+ | Additional targeting modes | Single oligos, restriction sites, UCE probes, reference/profile seeds; separate future scope |

Increasing input capacity in v4.0 should already help some nuclear targets.
The new low-coverage policy and linkage-dependent sensitivity belong in v4.1.

## Incremental k-mer counting decision

**Retain the feature as an isolated legacy analysis path through v4.x.
Remove its machinery from the normal sPCR counting path.**

- Preserve explicit `--chunks > 0` use, existing histogram formats, and the
  existing viewer workflow. Maintenance and correctness fixes continue;
  new rarefaction features are outside this program.
- `--chunks 0` must create one logical exact count store, without chunk
  distribution, incremental histogram snapshots, or copying a sole table.
- Share byte parsing/encoding where useful, but do not require every counter
  backend to support incremental snapshots or legacy chunk merging.
- Preserve combined legacy histogram/PCR invocation through an adapter to the
  finalized counts. Document unsupported new-backend/legacy combinations
  explicitly rather than silently changing their semantics or resource limits.
- A final count histogram is distinct from incremental rarefaction and may be
  supported by the new counter without reintroducing incremental machinery.
- Removing incremental counting entirely is not authorized by this plan.
  Any later removal needs an explicit decision, migration/deprecation plan,
  and treatment of `sharkmer_viewer`.

## Architecture and resource contracts

Agree these contracts in [#139](https://github.com/caseywdunn/sharkmer/issues/139) before independent implementations.

1. **Exact counts and membership.** Preserve canonical identity, all observed
   accepted k-mers, and abundance through ingestion/finalization. Specify the
   maximum exact count and observable overflow behavior. Do not default to
   saturating u8 counts or probabilistic membership. Packed u16 with exact
   overflow storage is a benchmark candidate, not a predetermined winner.
2. **Canonical keys versus oriented identities.** Count keys and directional
   primer/graph identities are separate concepts and should have explicit
   APIs/types. The current graph uses (k-1)-mer nodes and k-mer edges.
3. **One pass over original input.** Counting must support non-replayable
   sources. Optional local spooling/replay for evidence and external counting
   does not violate this contract. Never require a second remote download.
4. **Bounded application resources.** Budget count-table growth/resizing, queues,
   concurrent target graphs, evidence, and lookup caches together. An external
   counter followed by a full in-memory final table is not a bounded-memory
   solution. Library/runtime/network internals use measured process headroom
   and existing platform safeguards, outside per-allocation admission.
5. **Hardware baseline.** Use a 16 GB laptop as the initial target, reserving
   operating-system headroom; validate larger jobs with disk-backed storage.
   A 32–64 GB machine is an additional benchmark class, not a prerequisite for
   the architecture. Record SSD/temp-space requirements and peak RSS.
6. **Replay and evidence.** Preserve the exact selected read subset, mate
   identity, read orientation, base positions, qualities, and invalid-base
   discontinuities. Stream against bounded batches of target graphs and keep
   compact evidence rather than all reads in RAM.
7. **Target-local graph work.** Benchmark unitigs and neighbor caching on
   target graphs first. Avoid constructing global background topology merely
   to answer a small panel.
8. **Explicit uncertainty.** No-data, contradictory evidence, unresolved
   phase/repeats, partial sequence, and exhausted search are different
   outcomes. Observed k-mers do not prove an observed full-length haplotype.
9. **Reproducible storage and execution.** Persist k, count precision,
   canonicalization, preprocessing, input/subset identity, and schema version.
   Fingerprint tested binaries and record actual parameters.

## v4.0 — Scalable exact counting and replay

- [x] [#171](https://github.com/caseywdunn/sharkmer/issues/171) — Clarify current developer guidance and preserve historical entry points. **Documentation gate: met 2026-11-09.** [Index](README.md), source-checked [overview](OVERVIEW.md), and relative-link/archive-target validation are complete.

- [ ] [#159](https://github.com/caseywdunn/sharkmer/issues/159) — Complete the
  [sequence-annotation interface cleanup](SEQUENCE_ANNOTATIONS.md#v40-interface-and-migration).
  Development writers/readers use `target`, `sequence_index`, `sequence_length`
  and structured `assembly_rank_score`; stable IDs and `kmer_count_*` remain.
  **Local compatibility gate: met 2026-11-12.** The
  integrated verification (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
  passes 316 Rust tests per hasher, 223 Python tests and 160 independently
  audited native invocations, preserving measurements and provenance and
  testing historical/new consumers. Git integration and CI remain open before
  closing the issue. The v3.2 warning revision did not implement or waive this work; default read-support
  assessment and new biological claims remain separate inference scope.
- [ ] [#139](https://github.com/caseywdunn/sharkmer/issues/139) — Define exact count-store, replay-source, and assembly evidence interfaces. [Contracts and implementation status](COUNT_STORE_CONTRACTS.md). Typed identity/extraction gate: **met 2026-11-09**, [correctness evidence](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-identity-contracts-20261109/README.md). Exact store/shared-ledger conformance gate: **met 2026-11-10**, [evidence](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-exact-store-contracts-20261110/README.md); Assembly-access gate: **met 2026-11-10**, [exact-store PCR, u64 counts and propagated-error evidence](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-assembly-access-20261110/README.md). Production counter dispatch, replay and combined-budget integration remain open. [Frozen baseline protocol and promotion criteria](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-performance-baseline-20261109/README.md) precede backend selection.
  Typed-exhaustion conformance (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
  and real temporary-file ownership checks (local evidence: `benchmarks/benchmark_results/v4-owned-disk-file-20261112/README.md`)
  now pass under both hashers. They establish configured-limit error handling
  and logical file-extent ownership; end-to-end storage integrity, replay and
  whole-run resource enforcement remain open.
  Primary graph-buffer/admission checks (local evidence: `benchmarks/benchmark_results/v4-graph-buffers-20261113/README.md`)
  pass 401 Rust tests per hasher and 224 Python tests, including eight exact
  candidates at threads 1/2/4/8 and real growth refusal. These hidden opt-in
  controls cover primary graph buffers, not auxiliary/evidence or whole-run
  memory. Rerun preservation before publication (local evidence: `benchmarks/benchmark_results/v4-rerun-preservation-20261113/README.md`)
  passes 404 Rust tests per hasher and 224 Python tests, preserving all accepted
  output bytes after a real graph-budget failure. The
  durable publication recovery implementation (local evidence: `benchmarks/benchmark_results/v4-publication-recovery-20261113/README.md`)
  passes full quality checks (409 Rust tests per hasher, 224 Python tests).
  **Scoped publication failure/restart gate: met 2026-11-13.** All 160 native
  EIO/ENOSPC/SIGTERM cases (480 invocations) and their independent audit pass
  under both hashers and threads 1/2/4/8, preserving exact accepted output bytes.
  Pre-journal orphan reclamation and whole-run budgets remain open.
  Node-lookup ownership (local evidence: `benchmarks/benchmark_results/v4-graph-lookups-20261113/README.md`)
  now shares the target ceiling, charging old/replacement buffers before growth
  and preserving stable IDs. Full local checks pass 412 Rust tests per hasher
  and 224 Python tests. The subsequent
  extension/pruning/repeat workspace migration (local evidence: `benchmarks/benchmark_results/v4-graph-workspaces-20261113/README.md`)
  passes full quality checks (421 Rust tests per hasher and 224 Python tests).
  **Scoped workspace conformance gate: met 2026-11-13.** All 128 native
  candidate cases plus two frozen previous-build references and independent
  auditing pass, preserving FASTA/PCR and exact selected-read/count identities.
  Primer-buffer ownership and bounded candidate selection (local evidence: `benchmarks/benchmark_results/v4-primer-buffers-20261113/README.md`)
  now pass 428 Rust tests per hasher and 224 Python tests, plus 512 native
  candidate cases/eight prior-build references and independent auditing.
  **Scoped primer buffer/selection gate: met 2026-11-13.** Exact/IUPAC primers
  and caps retain the same FASTA/PCR evidence, mismatch priority and counts.
  Path-search buffer ownership (local evidence: `benchmarks/benchmark_results/v4-path-buffers-20261113/README.md`)
  passes full local quality checks (434 Rust tests per hasher and 224 Python
  tests), including the unchanged 8 MiB long-target gate. DFS state, retained
  candidate paths and metric scratch now share the target ceiling. All 512
  native candidate cases, eight frozen prior-build references and independent
  auditing pass, preserving complete FASTA/PCR and count/read identities.
  **Scoped path-buffer conformance gate: met 2026-11-13.**
  Read-evidence/bubble ownership (local evidence: `benchmarks/benchmark_results/v4-evidence-buffers-20261113/README.md`)
  passes full local quality checks (440 Rust tests per hasher and 224 Python
  tests). Threading and bubble buffers share the target ceiling; support
  counters reject overflow. All 768 native candidate cases, 12 frozen prior-build
  references and independent auditing pass, including paired reverse mates,
  N gaps and an unmatched tail. **Scoped read-evidence/bubble conformance gate:
  met 2026-11-13.** Ordinary sequence-result ownership (local evidence: `benchmarks/benchmark_results/v4-result-records-20261113/README.md`)
  passes full local quality checks (447 Rust tests per hasher and 224 Python
  tests). A separate run-level result class
  preserves the combined ceiling; collected ordinary records release target
  admission. All 768 native candidate cases, 12 frozen prior-build references
  and independent auditing pass with identical FASTA/PCR evidence and read/count
  identities. **Scoped ordinary-sequence ownership gate: met 2026-11-13.**
  Diagnostic/warning payload ownership (local evidence: `benchmarks/benchmark_results/v4-warning-payloads-20261113/README.md`)
  passes full quality checks (451 Rust tests per hasher and 224 Python tests).
  Parent-run guards cover retained
  sequences, hashes, markers, source/product/key buffers and hidden copies;
  fixed labels are borrowed. All 1,024 native candidate cases, 16 frozen prior-build
  references and independent auditing pass, including hidden-diagnostic modes.
  **Scoped diagnostic/warning ownership gate: met 2026-11-13.**
  Stats metadata and parallel collection ownership (local evidence: `benchmarks/benchmark_results/v4-stats-metadata-20261113/README.md`)
  passes full quality checks (457 Rust tests per hasher and 224 Python tests).
  Target/run/input metadata and
  aggregate vectors retain parent guards; indexed parallel output reserves
  before callbacks. All 1,024 native candidate cases, 16 frozen prior-build
  references and independent auditing pass, preserving full FASTA/PCR and
  run/input/command metadata. **Scoped stats-metadata ownership gate: met
  2026-11-13.** Pinned YAML workspace accounting (local evidence: `benchmarks/benchmark_results/v4-yaml-workspace-20261113/README.md`)
  passes full quality checks (464 Rust tests per hasher and 224 Python tests),
  including allocation-free preflight and 80 isolated allocator observations.
  The reservation precedes stats staging and outlives emitter cleanup; fields
  and YAML bytes remain unchanged. All 1,024 native candidate cases, 16 frozen
  references and independent count/input-frame auditing pass. **Scoped YAML
  workspace ownership gate: met 2026-11-13.**
  Pinned deduplication scratch (local evidence: `benchmarks/benchmark_results/v4-dedup-workspace-20261113/README.md`)
  passes full quality checks (470 Rust tests per hasher and 224 Python tests).
  Both arrays are reserved before the distance call and released through its
  normal return/unwinding, with 48 exact allocator observations and independent
  full-matrix/previous-adapter distance agreement. Scratch refusal cleans up
  kept and unconsumed owned records. All 1,024 native candidates, 16 frozen
  reference runs and independent count/input-frame auditing pass. **Scoped
  deduplication workspace ownership gate: met 2026-11-13.**
  Shared startup ledger and manifest workspace (local evidence: `benchmarks/benchmark_results/v4-startup-ledger-20261113/README.md`)
  pass full quality checks (479 Rust tests per hasher and 224 Python tests).
  One source-checked numeric plan preserves all resource formulas, prepares the
  parent ledger before output resources, and defers counter allocation until
  ingestion. Count-phase peaks exclude completed startup work without resetting
  lifetime peaks. Explicit planned output-lease allowance preserves the prior
  descriptor pool; manifest workspace can refuse before temporary creation and
  preserve an existing completed generation. Native hard/soft-64-descriptor
  comparisons pass all 1,024 candidates and 16 frozen references; independent
  auditing confirms the unchanged pool and metadata transition. **Scoped
  shared-ledger/manifest-workspace gate: met 2026-11-13.**
  New-generation publication metadata (local evidence: `benchmarks/benchmark_results/v4-publication-metadata-20261113/README.md`)
  passes full quality checks (485 Rust tests per hasher and 224 Python tests)
  and 160 independently audited publication-fault cases (480 processes).
  Transaction paths, labels, staged receipts and hashes share the parent;
  manifests and recovery journals serialize borrowed metadata under reserved
  emitter workspace. All 1,024 native candidates and 16 frozen references pass
  complete FASTA/PCR, count/input and stable metadata checks under actual
  hard/soft-64-descriptor limits. Independent and archive verification pass.
  **Scoped new-generation publication metadata ownership gate: met
  2026-11-13.**
  Prior-output verification and recovery workspace (local evidence: `benchmarks/benchmark_results/v4-recovery-workspaces-20261113/README.md`)
  passes full quality checks (489 Rust tests per hasher, 224 Python tests).
  Cleanup paths and duplicate tables share the parent ledger; recovery uses
  sorted borrowed receipt indexes and admits its complete restore plan before
  public mutation. Typed refusal preserves all public files and backups, and
  recovery retries after competing memory releases. All 1,024 native candidates
  and 16 frozen references pass complete FASTA/PCR, count/input identity and
  stable metadata checks under actual hard/soft-64-descriptor limits. All 160
  fault/restart cases (480 processes), independent audits and archive
  verification pass. **Scoped recovery workspace ownership gate: met
  2026-11-13.**
  Historical input and typed-model ownership (local evidence: `benchmarks/benchmark_results/v4-parsed-metadata-20261113/README.md`)
  passes full quality (493 Rust tests per hasher, 224 Python tests). Serde seeds
  admit original input, scalar copies and receipt-vector growth, including
  aliases, and retain aggregate guards through model destruction. Historical
  parser forms and typed refusal preserve completed outputs. All 1,024 native
  candidates and 16 frozen references agree on complete FASTA/PCR, identities
  and stable metadata under actual hard/soft-64-descriptor limits. All 160
  fault cases (480 processes), independent audits and archive checks pass.
  **Scoped input/typed-metadata ownership gate: met 2026-11-13.**
  Historical loader workspace (local evidence: `benchmarks/benchmark_results/v4-loader-workspace-20261113/README.md`)
  now reserves those buffers before decoding under the same parent ledger,
  including tag-prefix reuse, ignored values and rejected second documents.
  Full checks pass 495 Rust tests per hasher and 224 Python tests; 56 isolated
  loader allocation observations cover the audited bound. All 1,024 native
  candidates and 16 frozen references agree; 160 fault cases (480 processes),
  independent audits and archive verification pass. **Scoped loader workspace
  ownership gate: met 2026-11-13.** The
  bound can conservatively overreserve scalar text; no original protocol ceiling
  is raised. Returned errors and other process/filesystem owners remain open.
  The combined-budget million-read trial (local evidence: `benchmarks/benchmark_results/v4-combined-budget-20261113/README.md`)
  and independent audit also pass with graph/result owners enabled: unchanged
  256 MiB address-space, table/sort/disk/timeout ceilings plus actual hard/soft
  64-descriptor limits. All seven FASTA outputs and complete PCR evidence agree;
  each count/replay file exceeds the address-space cap. Native RSS is 31,760
  KiB and every target fits the original 8 MiB allowance. This is scoped
  integration evidence, not a whole-run resource or production-default gate.
  Publication copy workspace (local evidence: `benchmarks/benchmark_results/v4-publication-copy-20261113/README.md`)
  now preadmits one 64 KiB backing allocation before backup creation or public
  recovery mutation. Targeted checks pass 38 output tests and Clippy per hasher,
  including refusal/retry and exact isolated heap accounting. Full quality
  passes 498 Rust tests per hasher and 224 Python tests. All 1,040 native
  comparisons, 160 fault cases (480 processes), independent audits and retained
  archive checks pass. **Scoped copy workspace ownership gate: met 2026-11-13.**
  Original protocol ceilings remain unchanged; broader publication/process
  ownership is still open.
  Optional Linux address-space safeguard (local evidence: `benchmarks/benchmark_results/v4-process-memory-20261113/README.md`)
  passes 501 Rust tests per hasher, 224 Python tests, 1,040 native comparisons,
  160 fault cases and 78 cap/startup probes. The program also completes the
  million-read reference trial under its own 256 MiB cap. All 92 retained files,
  four source archives and source/binary/audit joins pass verification.
  **Scoped optional process-cap gate: met 2026-11-13.** This optional control
  does not choose production defaults or waive the remaining owner contracts.
  Other library/process workspace and publication accounting remain under the
  [remaining owner contracts](COUNT_STORE_CONTRACTS.md).
  Remote cache-miss record buffers (local evidence: `benchmarks/benchmark_results/v4-cache-records-20261113/README.md`)
  now borrow the exact counter's bounded reader and shared input ledger. All
  32 cache tests and Clippy pass under both hashers, including real HTTP
  refusal/preservation. Full checks pass 507 Rust tests per hasher and 225 Python
  tests, including cache miss/hit/count/replay integration. Native auditing and
  retention pass. **Scoped cache-record buffer gate: met 2026-11-13.**
  Cache-metadata ownership (local evidence: `benchmarks/benchmark_results/v4-cache-metadata-20261114/README.md`)
  now bounds sidecar input/output, parser/emitter workspace and surviving strings;
  its closed temporary remains charged through directory sync. Full checks
  pass 560 Rust tests per hasher and 225 Python tests. Four independently audited
  native component invocations pass 80 tests under actual 64-descriptor/512 MiB
  limits. **Scoped cache-metadata gate: met 2026-11-14.**
  Compressed cache extents (local evidence: `benchmarks/benchmark_results/v4-cache-data-20261114/README.md`)
  now use the common owned writer through gzip/buffer finalization, file close
  and cache-directory sync. Typed quota refusal preserves prior cache pairs
  and foreign files. Full checks pass 562 Rust tests per hasher and 225 Python
  tests; all 84 cache/HTTP component tests pass under actual hard limits.
  The same snapshot passes 544 independently audited registered-workspace
  I/O/retry cases (1,360 native invocations), covering EIO, ENOSPC, EINTR and
  short writes at 11 boundaries. Source/archive/raw-payload retention passes.
  **Scoped compressed-cache extent and registered-workspace I/O gates: met
  2026-11-14.**
  Registered cache publication/restoration (local evidence: `benchmarks/benchmark_results/v4-cache-recovery-20261114/README.md`)
  now binds both generations with fixed checksummed receipts, keeps old backups
  through failed restoration, and admits the common private inode inventory
  before recovery mutation. Full checks pass 570 Rust tests per hasher and
  225 Python tests; 100 hard-limit cache/HTTP component tests pass. The native
  matrix and independent audit pass 544 cases/1,488 invocations and 212 actual
  SIGTERMs. Retained verification checks three source archives and all 16,196
  raw payloads, including failed source/harness versions. **Scoped registered-cache
  publication/restoration I/O/SIGTERM gate: met 2026-11-14.** Codec/transport
  workspace and broader process/filesystem ownership remain open.
  Application FASTQ-buffer/gzip backing (local evidence: `benchmarks/benchmark_results/v4-codec-20261114/README.md`)
  now preadmits pinned backend state, nested compression allocations, optional
  member-header growth overlap and source/flush buffers. Full checks pass 576
  Rust tests per hasher and 225 Python tests. All 112 hard-limit component tests,
  544 native cache fault cases and independent audits pass. Isolated allocator
  observations stay inside the envelope and release all allocations; source,
  dependency and all 15,447 raw-payload archive checks pass. **Scoped application
  FASTQ-buffer/gzip backing gate: met 2026-11-14.** HTTP/TLS/content-encoding
  and preprocessing, broader descriptor/physical ownership and final-source
  promotion remain open.
  ENA response/model/sample ownership (local evidence: `benchmarks/benchmark_results/v4-ena-20261114/README.md`)
  now shares the original count ledger before metadata resolution, preserves
  the existing 10 MiB body limit, lossy UTF-8 and first-row selection, and keeps
  URL owners through read-plan transfer. Full checks pass 583 Rust tests per
  hasher and 225 Python tests; all 126 hard-limit component tests and 544 native
  workspace I/O/retry cases (1,360 invocations) pass independent auditing.
  Retention verifies both source archives, 76 dependency files and all 13,342
  raw payloads. The first full check caught unused overflowing count settings
  breaking dry runs; the preserved correction passes the unchanged early-exit
  test under both hashers. **Scoped ENA body/model/sample ownership gate: met
  2026-11-14.** HTTP/TLS/content-encoding, panel/CLI preprocessing, broader
  descriptor/physical ownership and final-source promotion remain open.
  Automatic HTTP gzip ownership (local evidence: `benchmarks/benchmark_results/v4-http-codec-20261114/README.md`)
  now preadmits ureq's eager content decoder and keeps it through body-reader
  and retained HTTP-error disposal, separately from application FASTQ gzip.
  Full checks pass 589 Rust tests per hasher and 225 Python tests; all 138
  hard-limit component tests and 544 cache fault cases (1,488 invocations,
  212 actual SIGTERMs) pass independent auditing. Source/dependency retention
  and all 15,447 raw payload checks pass. Isolated allocation peaks stay within
  the 305,671-byte envelope and release. **Scoped automatic HTTP gzip gate: met
  2026-11-14.** Request/header/TLS, earlier panel/CLI preprocessing and the
  broader process/filesystem and final-source gates remain open.
  Workspace-parent namespace identity (local evidence: `benchmarks/benchmark_results/v4-workspace-identity-20261114/README.md`)
  now checks kernel objects before accepting a missing child or cleaning private
  data. Full checks pass 595 Rust tests per hasher and 225 Python tests; all 168
  hard-limit component tests and 544 workspace I/O/retry cases pass. The two-node
  trial passes 32 cases/96 invocations and 32 actual SIGTERMs, independently
  auditing shared recovery, local refusal and recovery on the live origin.
  Five source archives, 94 dependency files and all 28,664 raw payloads pass
  retention, including the original missing-parent and endpoint-identity failures.
  **Scoped workspace-parent namespace identity and cross-host recovery gate:
  met 2026-11-14.**
  Markerless development-journal recovery (local evidence: `benchmarks/benchmark_results/v4-legacy-journal-20261114/README.md`)
  now persists ownership and accepted-manifest anchors before journal removal.
  Full checks pass 597 Rust tests per hasher and 225 Python tests; independent
  auditing passes 346 hard-limit component tests and 1,088 native fault cases
  (2,184 invocations, 276 actual SIGTERMs). Three source archives, 126 live joins,
  94 dependency files and all 30,344 raw payloads pass retention. **Scoped
  markerless development-journal interruption gate: met 2026-11-14.** Non-Linux
  native proof, broader physical/descriptor ownership, library/preprocessing
  resources and final-source promotion remain open.
- [ ] [#140](https://github.com/caseywdunn/sharkmer/issues/140) — Remove redundant
  table consolidation and per-read allocations from sPCR counting.
  Direct-ingestion correctness gate: **met 2026-11-12**,
  [verification](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-direct-refinement-20261112/README.md).
  Default selection and quantitative promotion gates remain **not met**.
  [Custom-table trials](https://github.com/caseywdunn/sharkmer-benchmarks/blob/main/benchmarks/benchmark_results/v4-direct-comparison-20261112/README.md)
  and Hashbrown trials through V7 (local evidence: `benchmarks/benchmark_results/v4-hashbrown-comparison-20261112/README.md`)
  retain complete paired audits and their failed speed decisions. V7 agrees on
  all 100 invocations but misses the geometric-mean speed gate and the Human
  wall allowance. Its conditional matrix is explicitly skipped.
  V8's reusable key buffer (local evidence: `benchmarks/benchmark_results/v4-reusable-buffer-checks-20261112/README.md`)
  passes conformance but is slower in its one-input attribution pilot; no
  promotion follows. Integrated annotation/counter checks (local evidence: `benchmarks/benchmark_results/v4-integrated-checks-20261112/README.md`)
  pass both backends/hashers at four thread counts. The comparison selector
  stays hidden, explicit incremental counting stays isolated, and the original
  speed/RSS thresholds remain unchanged.
  Complete-record parsing (local evidence: `benchmarks/benchmark_results/v4-buffered-records-20261112/README.md`)
  and batch-validated counting (local evidence: `benchmarks/benchmark_results/v4-unit-batches-20261112/README.md`)
  pass correctness but show no Human speed gain; both production optimizations
  were removed. Isolated attribution (local evidence: `benchmarks/benchmark_results/v4-identity-cost-20261112/README.md`)
  locates about 0.27 s in selected-record hashing on CPUs without SHA
  instructions; a newer library does not improve it. The hidden
  ordered identity/count pipeline (local evidence: `benchmarks/benchmark_results/v4-identity-pipeline-20261112/README.md`)
  now passes 374 Rust tests per hasher, 224 Python tests and 1,008 independently
  audited native invocations, including 223 actual spills. Its Human pilot
  improves on the previous parallel producer but still trails v3.2.1.
  Expanded native failure checks pass 375 Rust tests per hasher and 224 Python
  tests, covering post-spill disk/invalid-base errors and corrupt concatenated
  gzip input at threads 1/2/4/8, with failed publication and owned-file cleanup.
  The full paired audit preserves all 287 records/counts/source identities and
  passes RSS allowances, but the 1.0139565284 wall ratio and Human allowance
  fail promotion; its conditional matrix is skipped. The pipeline's
  hard-limit replay/PCR trial failed on allocation during assembly despite low
  RSS; its failure is retained. A
  coordinator-participation correction (local evidence: `benchmarks/benchmark_results/v4-current-thread-20261112/README.md`)
  passes 376 Rust tests per hasher and 224 Python tests, including actual OS
  thread limits across ingestion/replay/assembly at 1/2/4/8. Its unchanged-cap
  million-read retry and independent audit pass with seven byte-identical
  FASTA outputs, complete reference PCR evidence and 12,828 KiB peak RSS.
  Its Human pilot meets the absolute per-input allowance but does not establish
  the complete frozen performance gate. The
  hybrid owner-transfer experiment (local evidence: `benchmarks/benchmark_results/v4-hybrid-owners-20261112/README.md`)
  now passes 388 Rust tests per hasher, 224 Python tests, 1,512 independently
  audited native cases (56 transfers, 202 spills) and phase-start/I/O failures.
  Its full paired comparison (local evidence: `benchmarks/benchmark_results/v4-hybrid-owners-20261112/comparison-v1/summary.json`)
  passes all 100 correctness checks and per-input wall/RSS allowances, but the
  0.9656945923 geometric-mean wall ratio misses the 0.95 criterion. Its conditional
  matrix is explicitly skipped because that prerequisite failed.
  The default-handoff hard-limit trial (local evidence: `benchmarks/benchmark_results/v4-hybrid-owners-20261112/bounded/summary.json`)
  passes counting and sPCR under a 256 MiB address-space cap, with exact full
  input/count/PCR/FASTA agreement and at most two sampled OS threads.
  Isolated SHA provider profiling (local evidence: `benchmarks/benchmark_results/v4-sha-backends-20261112/README.md`)
  finds lower hash cost with an already locked TLS dependency. A
  native provider trial (local evidence: `benchmarks/benchmark_results/v4-native-sha-20261113/README.md`)
  passes 390 Rust tests per hasher and 49 independently audited Human runs,
  with identical full identity/counts/FASTA/PCR evidence and 6–7% lower wall
  time than the current exact candidate. The
  integrated revision (local evidence: `benchmarks/benchmark_results/v4-selected-sha-20261113/README.md`)
  passes combined checks (390 Rust tests per hasher, 224 Python tests) and
  independent native owner-start/I/O-failure auditing. Its 1,512-case native
  input matrix and actual-default-handoff hard-limit replay/PCR trial also pass,
  including complete reference PCR/FASTA agreement under a 256 MiB cap.
  The full paired comparison (local evidence: `benchmarks/benchmark_results/v4-selected-sha-20261113/comparison-v1/summary.json`)
  passes all 100 correctness checks and RSS allowances but misses aggregate
  speed (0.9632352084 versus 0.95) and Human wall time (0.67 versus 0.66 seconds
  allowed). Its conditional matrix is explicitly skipped; default promotion
  and full quantitative acceptance remain open.
  A 27-run ownership diagnostic (local evidence: `benchmarks/benchmark_results/v4-owner-attribution-20261113/README.md`)
  confirms exact agreement but finds fixed counting ownership slower than
  hybrid on Human after the provider change; it does not justify replacement.
  Native CPU profiling (local evidence: `benchmarks/benchmark_results/v4-native-profile-20261113/README.md`)
  and its independent audit locate measurable adaptive-wrapper overhead. An
  isolated inlining candidate (local evidence: `benchmarks/benchmark_results/v4-adaptive-inline-20261113/README.md`)
  passes all 495 Rust tests per hasher and a 20-run independently audited Human
  pilot (0.52 versus 0.57 s for the frozen build, unchanged exact results).
  Its full ten-input comparison passes all 100 invocations and original speed,
  per-input wall and RSS criteria (geometric-mean ratio 0.6799958680). The
  independent 80-case thread/hasher audit also passes. The attribute is now
  applied to the working source. A
  resource-envelope pilot (local evidence: `benchmarks/benchmark_results/v4-default-resources-20261113/README.md`)
  independently verifies all ten inputs and 287 records under the proposed
  16 GB hardware configuration. Normal commands now select the exact hybrid
  counter, shared graph accounting and Linux address-space safeguard; the
  integrated default checkpoint (local evidence: `benchmarks/benchmark_results/v4-production-defaults-20261113/README.md`)
  passes its full checks and unchanged 100-run performance/80-case matrix
  criteria. **Scoped normal-default snapshot gate: met 2026-11-13.** Its
  geometric-mean elapsed ratio is 0.7760642119 with every original per-input
  wall/RSS gate met. Final-source default promotion remains **not met** pending
  remaining ownership work, integrated revalidation and Git/CI.
- [ ] [#141](https://github.com/caseywdunn/sharkmer/issues/141) — Parallelize ingestion with bounded queues and measured shard ownership.
  A hidden coordinator/owner candidate retains disjoint exact tables without
  consolidation and reuses three fixed key buffers per owner. Local checks pass
  340 Rust tests per hasher and 223 Python tests, including one-shot clean/warning
  PCR at threads 1/2/4/8, owner failure propagation and joined cleanup. The full
  native input matrix (local evidence: `benchmarks/benchmark_results/v4-parallel-counts-20261112/README.md`)
  passes all 1,512 invocations and independent raw-artifact auditing. The frozen
  paired comparison passes independent auditing of all 100 invocations. Its
  0.8587603682 geometric-mean wall ratio and all RSS allowances pass, but Human
  takes 1.05 s against a 0.97 s allowance: promotion remains **not met**. A compact
  eight-byte worker-key queue revision passes 352 Rust tests per hasher and 224
  Python tests; its Human pilot shows no ingestion improvement. A canonical-only
  extraction revision passes 354 Rust tests per hasher and 224 Python tests,
  but also trails the released baseline in its Human diagnostic. Neither pilot
  replaces the frozen promotion protocol. Production selection and whole-process
  integration remain open.
- [ ] [#142](https://github.com/caseywdunn/sharkmer/issues/142) — Benchmark packed k-mer tables while retaining exact abundance.
  A hidden byte-column candidate retains complete canonical keys and exact u64
  counts, with charged old/new resize allocations and move finalization.
  Local integration passes 352 Rust tests per hasher and 224 Python tests.
  The native input and synthetic attribution audits (local evidence: `benchmarks/benchmark_results/v4-packed-counts-20261112/README.md`)
  pass 1,512 invocations and 96 measurements respectively. Packing uses about
  24% less owned table storage at k=19 but has slower synthetic insertion/lookup.
  **Local evaluation gate: met 2026-11-12; default promotion: not met.** The
  unchanged real-data protocol and independent audit pass all 100 invocations,
  exact counts and 287 records. Every RSS allowance passes, but every input
  wall-time allowance fails; the geometric-mean wall ratio is 1.3477940354.
  This candidate is not selected. Git integration and CI remain open.
- [ ] [#143](https://github.com/caseywdunn/sharkmer/issues/143) — Count non-replayable datasets exactly within a memory budget.
  The [external development path](COUNT_STORE_CONTRACTS.md) now counts one-shot
  stdin with a fixed sort buffer and bounded streaming merges. Local checks
  pass 336 Rust tests per hasher and 223 Python tests; both candidate input-form
  matrices pass 1,512 independently checked invocations. **Scoped beyond-budget
  count/PCR gate: met 2026-11-12.** The frozen million-read stdin trial completes
  counting and PCR under an observed hard 256 MiB address-space limit with
  13,480 KiB peak RSS and a 1.19 GB count file. The
  independent audit and retained evidence (local evidence: `benchmarks/benchmark_results/v4-external-counts-20261112/README.md`)
  verify exact counts, seven sequences and full reference PCR evidence.
  A hidden serial adaptive candidate now transfers accepted in-memory counts
  into the external producer on typed count-growth pressure, without rereading
  the source. Local spill conformance (local evidence: `benchmarks/benchmark_results/v4-adaptive-spill-20261112/README.md`)
  passes 361 Rust tests per hasher and 224 Python tests, including complete
  one-shot PCR/replay agreement, failed-transfer cleanup and combined-memory
  pressure while another stage retains ownership. The independently audited
  native matrix passes all 1,008 invocations, including 224 actual spills.
  Independently spillable parallel owners (local evidence: `benchmarks/benchmark_results/v4-parallel-spill-20261112/README.md`)
  pass 365 Rust tests per hasher, 224 Python tests and all 1,008 native matrix
  invocations; independent auditing verifies counts, selection and storage
  composition, including 214 actual spill invocations. The hard 256 MiB
  parallel trial **failed** on allocation despite 47,132 KiB peak RSS. Deferred
  assembly-pool startup (local evidence: `benchmarks/benchmark_results/v4-deferred-pool-20261112/README.md`)
  passes 366 Rust tests per hasher and 224 Python tests, including an actual
  thread-lifetime probe. **Scoped parallel beyond-budget gate: met 2026-11-12.**
  The unchanged-cap retry and independent audit complete one million stdin reads
  with three disk shards, original-pass replay and complete reference PCR
  agreement under the observed hard 256 MiB address-space limit (55,332 KiB
  peak RSS). Source-bound archives retain the failure and successful retry.
  Default dispatch, production
  resource integration and full quantitative acceptance remain open.
  Native storage-fault auditing (local evidence: `benchmarks/benchmark_results/v4-storage-faults-20261112/README.md`)
  passes 48 controlled I/O-error cases and 24 SIGTERM/fresh-restart cases
  across two hashers and threads 1/2/4/8. Partial writes followed by ENOSPC and
  sync EIO cannot publish outputs. Interrupted count/replay files stay
  identifiable and unchanged; automatic reclamation/resume remains unimplemented.
  A shared count/replay descriptor cap (local evidence: `benchmarks/benchmark_results/v4-file-descriptors-20261112/README.md`)
  reserves before creating each owned file and releases after descriptor
  closure. Full checks pass 383 Rust tests per hasher and 224 Python tests;
  native I/O-fault/restart auditing also passes. Its own unchanged hard-limit
  trial and independent audit now pass with seven byte-identical outputs,
  full reference PCR evidence and 18,876 KiB peak RSS; native samples observe
  at most 21 descriptors and four OS threads. A
  process-aware revision (local evidence: `benchmarks/benchmark_results/v4-process-file-limits-20261113/README.md`)
  now clamps the pool to the OS soft limit, Linux startup census and additional
  headroom. It passes 392 Rust tests per hasher and 224 Python tests, including
  real low-limit/inherited-descriptor children, completed spill/replay/PCR and
  pre-creation refusal. All 1,512 independently audited input-matrix cases
  also pass under actual hard/soft 64-descriptor limits, with 56 handoffs and
  224 spills. Individual input/output/library ownership and complete process
  enforcement remain open.
- [ ] [#144](https://github.com/caseywdunn/sharkmer/issues/144) — Provide exact bounded-memory random lookup over external k-mer counts.
  Full-key checksummed disk lookup now feeds ordinary PCR directly. Local
  conformance covers corruption, foreign pages, sorted membership, reentrant
  lookup and charged concurrent scans. The scoped beyond-budget count/PCR gate
  passes as recorded under #143; broader storage/performance validation remains
  open, and the backend stays hidden.
- [ ] [#145](https://github.com/caseywdunn/sharkmer/issues/145) — Add bounded replay/spooling with consistent paired input semantics.
  Exact development candidates now spool selected original FASTQ fields during
  ingestion and stream them into opt-in read threading, including stdin. Focused
  checks pass both hashers for format/integrity, paired selection, one-shot PCR,
  disk-budget failure and cleanup. **Local integration gate: met 2026-11-12**,
  348 Rust tests per hasher and 224 Python tests (local evidence: `benchmarks/benchmark_results/v4-replay-20261112/README.md`).
  **Scoped beyond-budget replay/PCR gate: met 2026-11-12.** Independent auditing
  verifies the million-read stdin trial with a 1.19 GB count file and 429 MB
  replay file under a hard 256 MiB address-space limit, matching the complete
  threaded reference with 13,412 KiB peak RSS. The expanded native replay matrix
  passes all 2,688 valid invocations and independent auditing; complete archives
  preserve the bounded trial and all unsuccessful matrix setups. Legacy source reopening/read retention,
  production dispatch and whole-run/evidence budgets remain open.

Sequence: contracts → simple counter → bounded parallelism → packed-table
evaluation. External counting and replay can progress independently after their
declared prerequisites. Final disk lookup depends on external counts.
The packed-table experiment is not a dependency of the first external backend.

Compare existing ahash/Fx configurations, owner-managed shards, bounded local
aggregation, and a specialized shared table on representative workloads.
Do not allocate a complete per-thread map without measuring duplication.
Jellyfish is a relevant exact-counting comparator; KMC is a relevant external
counter. Their existence does not prescribe Sharkmer's internal design.

Start external storage with the simplest measured design. Minimizer/signature
partitioning and super-k-mers are candidates for reducing temporary I/O;
partition skew, exact global counts, and final random access remain required.

### Release gate

- Exact count and product agreement across supported threads and backends.
- Standard `--chunks 0` has no incremental counting overhead; legacy
  histogram/combined-PCR compatibility remains tested separately.
- Record stage throughput, table bytes per distinct k-mer, resize/finalization
  peaks, RSS, disk use, and lookup latency. Freeze quantitative performance
  gates against the repaired baseline before choosing production designs.
- Demonstrate a one-shot dataset whose global count store exceeds the chosen
  RAM budget completing counting **and sPCR**, using disk-backed lookup.
- Verify limits, paired selection, malformed gzip, interruption, and disk-full
  behavior without losing previously accepted counts or exposing partial files.
- Keep remote-input failures within the same accepted-output boundary. Measure
  RSS/headroom under existing process limits; do not require client-internal
  allocation proofs or promise graceful recovery from every allocation failure.
- Standard sPCR accuracy does not regress while counting implementation changes.

## v4.1 — Read evidence and low-coverage nuclear recovery

- [ ] [#129](https://github.com/caseywdunn/sharkmer/issues/129) — Complete the
  held-out registration and full nine-panel biological evaluation deferred
  from v3.2, before further assembly-policy tuning or broad accuracy claims.
  Preserve independent truth, callable targets, negative controls, and explicit
  unavailable/unassessed outcomes; regression success is not a substitute.
  Obtain an independently annotated Agalma ITS source without filling the gap
  with assembled products or misannotated sources. Review gene annotation
  independently of public sequence origin, and add independent same-sample
  truth and read-backed assessment of divergent/repeat regions. Cross-individual
  differences alone do not establish incorrect products.

- [ ] [#146](https://github.com/caseywdunn/sharkmer/issues/146) — Batch primer discovery without combinatorial variant enumeration.
- [ ] [#147](https://github.com/caseywdunn/sharkmer/issues/147) — Build length-bounded target graphs and compact nonbranching paths.
- [ ] [#148](https://github.com/caseywdunn/sharkmer/issues/148) — Use read-spanning branch histories to resolve repeats and constrain paths.
  Integrate bounded quality-preserving replay and positioned histories under
  #145/#116 behind an experimental opt-in path, with feature-off compatibility
  against its frozen baseline. Specify ambiguity-preserving variant-aware
  mapping before reconstructing alternatives. Marginal coverage does not
  justify repairing missing edges or relaxing repeat gates. Admission changes
  require cause-specific positive evidence, repeat/mixture controls, and fresh
  released-baseline high-copy comparisons that retain historical v3.1 identity
  checks.
- [ ] [#149](https://github.com/caseywdunn/sharkmer/issues/149) — Add evidence-based low-coverage nuclear recovery and partial outcomes.
- [ ] [#116](https://github.com/caseywdunn/sharkmer/issues/116) — Stream reads through bounded batches of target graphs; replace the
  old per-gene full-input scan proposal and all-read RAM vector.
- [ ] [#101](https://github.com/caseywdunn/sharkmer/issues/101) — Apply mate orientation and insert-size constraints to paths.
- [ ] [#99](https://github.com/caseywdunn/sharkmer/issues/99) — Use conservative evidence-aware pruning; lack of evidence is not
  automatically evidence of an erroneous edge.
- [ ] [#106](https://github.com/caseywdunn/sharkmer/issues/106) — Add explicit primer alternatives with current CLI/schema semantics
  and selection that does not erase minority alternatives.
- [ ] [#127](https://github.com/caseywdunn/sharkmer/issues/127) — Expose suggested-k guidance without silently switching global k.
- [ ] [#114](https://github.com/caseywdunn/sharkmer/issues/114) — Re-profile the proposed scoring optimization; implement only if
  it remains material after graph changes.

Dependency order: seeding/target graphs and replay → [#116](https://github.com/caseywdunn/sharkmer/issues/116) and [#101](https://github.com/caseywdunn/sharkmer/issues/101) →
linked-path evidence → low-coverage policy and [#99](https://github.com/caseywdunn/sharkmer/issues/99). [#106](https://github.com/caseywdunn/sharkmer/issues/106) and [#127](https://github.com/caseywdunn/sharkmer/issues/127) support
primer usability; [#114](https://github.com/caseywdunn/sharkmer/issues/114) is an optional measured optimization.

Maintain an explicit standard policy and add an opt-in low-coverage policy.
Use per-target/per-haplotype coverage, realistic intron spans, and
paralog/orthology evaluation. Admit singletons only under a specified evidence
model, not by globally relaxing all pruning. Return partial assemblies where
useful, with an explicit distinction from complete amplicons.

Evaluate local multi-k assembly on recruited/replayed reads after single-k
behavior is validated. Smaller k can restore overlap and larger k can resolve
some repeats; neither establishes unobserved bases or unsupported phase.
Do not make this experiment a requirement for shipping proven v4.1 gains.

### Release gate

- Interior reads and both strands contribute equivalent valid evidence.
- Evidence buffering is independent of total input length, with graph batch
  size and retained links bounded or spilled explicitly.
- Repeat lengths require read/fragment evidence; insufficient linkage yields
  ambiguity rather than a fabricated complete product.
- Publish held-out callable recall and sequence precision across nuclear depth,
  heterozygosity, intron length, and paralog similarity.
- Demonstrate that singleton rescue improves recall without an unreviewed
  increase in false products; physical gaps remain partial/unknown.
- Standard high-copy behavior remains a protected benchmark class.

## v4.2 — Metagenomic diversity

- [ ] [#150](https://github.com/caseywdunn/sharkmer/issues/150) — Preserve rare templates through sPCR seed, threshold, and path search.
- [ ] [#151](https://github.com/caseywdunn/sharkmer/issues/151) — Report supported haplotypes, uncertainty, and defensible target abundance.
- [ ] [#68](https://github.com/caseywdunn/sharkmer/issues/68) — Document downstream classifier/workflow integration using current
  flags and explicit limitations on marker resolution and organism abundance.

Add an explicit diversity policy. Explore plausible components and coverage
levels after a dominant product is found. Preserve exact distinct sequences by
default; optional similarity clustering is separate from assembly.

Use branch histories and mate constraints to preserve supported linkage.
If two distant variants have no distinguishing read/fragment evidence, report
phase blocks or ambiguity instead of enumerating possible combinations as
observed haplotypes. Estimate abundance from discriminating sequence and
compatible reads; do not equate shared marker coverage with organism counts.

### Release gate

- Recover both supported products in the 100:4 mixture with common primer ends.
- Evaluate 1 bp differences, indels, multiple similar templates, >40 seed
  candidates, and >20 genuine products under configured limits.
- Publish precision, rare-template recall, chimera rate, phasing accuracy, and
  abundance error stratified by target coverage and divergence.
- Every product carries evidence/completeness information; truncation is visible.
- No confident organism count or full haplotype is claimed when the marker or
  linkage is uninformative.

## Benchmark matrix and acceptance policy

Keep calibration and held-out evaluation separate. Freeze fixture identities,
truth, and target-specific tolerances before parameter sweeps.

- Existing rRNA/organelle successes, repetitive/AT-rich targets, and known
  failures; compare whole products rather than FASTA counts.
- Nuclear per-target/per-haplotype depths such as 1, 2, 3, 5, 10, 20, and 40x,
  including heterozygosity, introns, paralogs, and absent targets.
- Known metagenomic mixtures spanning abundance ratios and sequence divergence;
  vary absolute minority depth separately from abundance ratio.
- Error-free controls, realistic quality/error profiles, primer mismatches,
  physical gaps, and deliberately unresolvable repeat/phase cases.
- Plain/gzip/multi-member gzip, paired/interleaved/single-end, local/replayed/
  one-shot inputs; multiple k values and 1/2/4/8 threads where supported.
- Small exact-oracle tests, the included fixture, approximately 1M-read common
  workloads, and deep datasets that exceed the in-memory budget.

Record decompression/parsing, emission, table updates, merge/finalization,
seed discovery, graph building, threading, search, total wall time, peak RSS,
allocator heap separately, disk usage, and count-store lookup performance.
Benchmark without graph dumps for ordinary runtime; measure diagnostics
separately. Do not promise a speedup before the relevant measurements exist.

Required correctness gates are absolute. Biological improvements need reviewed
precision/recall thresholds on callable truth sets; sequencing depth alone
cannot guarantee complete recovery. Performance changes require reproducible
end-to-end improvement, not just a faster isolated hash benchmark.

## Implementation handoff

- Begin with the v4.0 contracts under #139 and annotation interface decisions
  under #159. Freeze regression fixtures and quantitative baseline gates
  before choosing production designs.
- An issue is complete only when its acceptance criteria, focused tests,
  affected documentation, and release benchmark evidence are satisfied.
- Follow [CONTRIBUTING.md](../CONTRIBUTING.md) for implementation branches,
  commits, quality gates, and release process.
- Update this plan and the associated issue in the same implementation change.
  Do not infer completion from a similarly named historical closed issue.
- Existing [#121](https://github.com/caseywdunn/sharkmer/issues/121)–[#127](https://github.com/caseywdunn/sharkmer/issues/127) and [#76](https://github.com/caseywdunn/sharkmer/issues/76) include stale or overlapping status; audit them
  against code before closing or duplicating panel/cleanup work.
- This program does not remove incremental counting or implement v5 targeting.

## Design references

These are comparison/design references, not dependencies on their software.

- [Jellyfish: exact parallel counting with a specialized shared table](https://pmc.ncbi.nlm.nih.gov/articles/PMC3051319/).
- [KMC 3: counting and manipulating k-mer statistics](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399).
- [Linked de Bruijn graphs: retaining read-scale connectivity](https://pubmed.ncbi.nlm.nih.gov/29554215/).
- [Bifrost: compacted graph construction and indexing](https://pmc.ncbi.nlm.nih.gov/articles/PMC7499882/).
- [metaSPAdes: local coverage and metagenomic assembly challenges](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

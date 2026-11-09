# v3.2.1 security maintenance release

Updated 2026-11-09. **Release authorization, engineering checks, and HTTPS
compatibility are met; package verification and publication are pending.** The user explicitly requested that the maintenance release
proceed now. No further release confirmation is required. Track the work in
[#167](https://github.com/caseywdunn/sharkmer/issues/167).

## Scope

Update locked Rustls 0.23.37 → 0.23.46 and its required rustls-webpki dependency
0.103.13 → 0.103.14, together with v3.2.1 release metadata. The inspected lockfile
diff changes only those two dependency versions and checksums. This patch makes
no assembly-policy, sequence-annotation, CLI, or output-format change. Keep the
published v3.2.0 tag immutable and carry the maintenance fix into development.

The update addresses four version-applicable Rustls advisories:

| Advisory | Applicability to v3.2.0 |
| --- | --- |
| [GHSA-2mjx-qc3c-rqvc](https://github.com/rustls/rustls/security/advisories/GHSA-2mjx-qc3c-rqvc), alert #11 | Applies to the enabled TLS 1.3 download path; some peer handshake messages are accepted at the wrong encryption level. |
| GHSA-h2v8-394m-p7vv | Triggering Brotli certificate decompression is disabled. |
| GHSA-w44r-93j7-qv28 | The selected ring provider has no hybrid key-exchange groups. |
| GHSA-c8m5-c952-6cm4 | Sharkmer uses no custom session store and does not log TLS session values. |

Alert #11 surfaced after v3.2.0 publication and prompted this maintenance
release. The maintainer states that transcript authentication still prevents
an on-path attacker from altering or completing the handshake; this is not a
demonstrated authentication bypass. The absent preconditions for the other
three advisories do not exempt alert #11. All four are addressed by the same
narrow dependency update.

## Release gates

| Gate | Status | Required evidence |
| --- | --- | --- |
| User authorization | **Met, 2026-11-09** | The user's explicit instruction to proceed with v3.2.1 now. |
| Final source and dependency review | **Met for the tested candidate, 2026-11-09** | [Source receipts](../benchmarks/benchmark_results/v3.2.1-release-20261109/README.md) confirm exactly two dependency updates plus root version metadata; no production-source, panel, schema, test, or check-script changes from v3.2.0. Preserve this scope through the release commit. |
| Full engineering release check | **Met, 2026-11-09** | [Job 32671325](../benchmarks/benchmark_results/v3.2.1-release-20261109/README.md) passes all eight groups: 283 Rust tests per backend, 213 Python tests without skips, formatting/Clippy, 158 references, ten schemas, nine exports, 19 loader dry runs, and fixture sequence oracles. |
| HTTPS download compatibility | **Met, 2026-11-09** | [Job 32674065](../benchmarks/benchmark_results/v3.2.1-release-20261109/https-v2/README.md) passes with the engineering-tested binary: pinned HTTPS panel over TLS 1.3, matching three-read local control, and untrusted-certificate rejection before any HTTP request. The earlier harness version assertion failed before network activity and is retained separately. |
| Local package build and smoke tests | Pending | Build the v3.2.1 recipe and test the installed `--help`, `--version`, `--list-panels`, and `--cite` commands. Bind the tested source archive, recipe, environment, package, and logs. |
| Publication and upstream delivery | Pending | Record the release commit, merge, tag and tag CI, GitHub release, maintenance/development integration, final source-archive checksum, and upstream Bioconda recipe status. |

The [exact-version advisory audit](../benchmarks/benchmark_results/v3.2.1-release-20261109/advisory-audit/advisory-audit.json)
finds no matching reviewed or unreviewed GitHub advisory for either selected
dependency version at the recorded query time. Its raw responses are retained;
this is a dated database check, not a guarantee against unknown vulnerabilities.

The full engineering check and bounded HTTPS check are the verification scope
for this dependency patch. The user-accepted v3.2 biological scope and its
limitations carry forward. The
[32-invocation comparison](../benchmarks/benchmark_results/v3.2-final-check-20261108/RESULTS.md)
and [v3.2 review](REVIEW_v3.2.md) retain the prior evidence; this patch does not
repeat that comparison or add a broad biological accuracy or performance claim.
Independent held-out registration and the full nine-panel biological campaign
remain deferred, not passed, under
[#129](https://github.com/caseywdunn/sharkmer/issues/129). They remain required
before further assembly-policy tuning or broad accuracy claims. This focused
maintenance checklist supplements the release procedure without reopening the
completed v3.2 implementation plan.

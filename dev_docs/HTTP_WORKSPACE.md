# HTTP scope and failure handling

**User-approved v4 scope, 2026-11-16:** delegate networking to the existing
client. V4 bounds application data buffers; it does not account for every
allocation inside URL parsing, response headers, connection pools, DNS, sockets,
TLS or automatic HTTP content decoding. These use process headroom and existing
platform limits. See the [resource boundary](COUNT_STORE_CONTRACTS.md#v4-resource-boundary).

## Current implementation

`src/http.rs` only delegates GET requests and response readers to `ureq`.
The primer-panel path retains its existing 10-second connect and 30-second
read timeouts. Ordinary FASTQ/ENA requests retain the client's defaults,
redirects, status errors and automatic gzip behavior. Sharkmer does not add
another transport, custom client leases or allocation callbacks.

Sharkmer-owned FASTQ record/source buffers, application gzip decoding,
ENA response/model storage, panel input/model storage, count/replay/lookup and
assembly buffers retain their limits. Panel and ENA body limits are unchanged.
The application ledger reports its own reservations, not client-internal heap
or process RSS. Existing process headroom and Linux address-space safeguards
remain separate from application budgets.

Network, timeout, HTTP-status and decoding errors fail the operation through
the caller's normal error path. Infallible allocation exhaustion may terminate
the process; graceful memory-error reporting and immediate cleanup are not
promised. Previously accepted data must survive and interrupted work must not
be recognized as complete. Cache/publication recovery remains tested.

## Retained dependency patch

The pinned client remains the licensed source at `vendor/ureq-2.12.1`, from
upstream commit `ba4555d54fd72e16fc737f7e8fcc812efb131260`. Its chunk-size parser
previously collected arbitrarily long framing into a Vec, exhausting memory
for a tiny body. The minimal streaming fix preserves accepted framing and
avoids that demonstrated failure. It is covered by
[`tests/http_chunk_framing.rs`](../tests/http_chunk_framing.rs), using a licensed
unmodified reference. A redundant middleware type parenthesis is removed to
keep warnings clean. Other client files are restored to upstream.

Do not extend this patch into a networking ownership framework. Prefer an
upstream fix or dependency update for future concrete client bugs, with the
smallest compatible local patch only when needed. Supported source builds
must preserve the existing patch and licenses. Formatting checks target only
the application, patched decoder and middleware; untouched upstream files
are not reformatted.

## Verification and historical evidence

The simplification checks (local evidence: `benchmarks/benchmark_results/v4-http-scope-20261116/README.md`)
pass 625 Rust tests per hasher, 225 Python tests and all quality checks. Four
transport regressions pass separately per hasher at actual 64-descriptor/512-MiB
limits. The unchanged cache fault/recovery matrix and independent audit pass
544 cases, 1,488 invocations and 212 actual SIGTERMs at those same ceilings.
Earlier test counts and sealed snapshots remain historical evidence.

The chunk-framing checkpoint (local evidence: `benchmarks/benchmark_results/v4-http-chunk-framing-20261114/README.md`)
retains the original failure, streaming fix, native low-memory replay and
source/dependency evidence. The former detailed accounting experiments are
preserved in the codec (local evidence: `benchmarks/benchmark_results/v4-http-codec-20261114/README.md`),
response (local evidence: `benchmarks/benchmark_results/v4-http-metadata-20261114/README.md`),
raw request (local evidence: `benchmarks/benchmark_results/v4-http-request-text-20261115/README.md`),
URL (local evidence: `benchmarks/benchmark_results/v4-http-url-workspace-20261115/README.md`),
derived request (local evidence: `benchmarks/benchmark_results/v4-http-request-buffers-20261116/README.md`)
and Agent (local evidence: `benchmarks/benchmark_results/v4-http-agent-backing-20261116/README.md`)
checkpoints. Their per-allocation formulas and hooks are retired. Preserve
these frozen records; do not revive their pending owner lists as release gates.

# Sharkmer's pinned client patch

This directory contains the licensed ureq 2.12.1 crate source. Its original
MIT and Apache licenses, including the chunk decoder's notices and licenses,
are retained. `.cargo_vcs_info.json` records upstream provenance. Cargo uses
this path dependency for the documented Git checkout/source builds.

The framing patch is in `src/chunked/decoder.rs`. A separate redundant
parenthesis is removed from `src/middleware.rs` to keep compiler warnings
clean on the project toolchain; its type is unchanged. Chunk-size
framing now uses a fixed-size UTF-8 scratch array, grammar state and checked
integer accumulation. It preserves the former UTF-8, Unicode trim, optional
plus, leading-zero, overflow, extension and error-consumption behavior while
avoiding an unbounded Vec for the size token. A 328-byte returned panel could
previously abort allocation when its chunk-size line had a long zero prefix.

`src/request.rs` and `src/unit.rs` also carry an optional URL-admission callback,
exported from `src/lib.rs`. It runs immediately before initial URL parsing and
each redirect join, and follows Unit reconstruction and cloning. Without a
callback the original parse/redirect/error path is unchanged. Sharkmer uses it
to reserve source-derived URL/IDNA workspace and retained clones on its shared
ledger, holding the owner through response/error destruction. Client limits,
configuration, versions and features are unchanged. The callback does not parse
or normalize URLs itself, and no new unsafe is introduced.

The callback's optional request-construction method also receives borrowed
parsed URL/configuration and actual header capacities before Unit creation,
both initially and after redirects. `header.rs` supplies a nonallocating backing
capacity accessor. Sharkmer admits method/header/authentication/prelude buffers
on a separate reservation before construction, with a fixed allowance for the
automatic gzip header that precedes URL parsing. The consuming call releases
that reservation before returning a body or retaining an HTTP error.

Status/header, URL, raw request text and gzip owners remain separate. The
remaining Agent/pool/resolver/socket/TLS and process overhead
requirements are outside these patches. See [HTTP workspace](../../dev_docs/HTTP_WORKSPACE.md)
for the audit and current verification status. Do not remove the path pin or
replace this directory without rerunning its parser, native HTTP and resource
checks. Registry/crate publication is not a documented distribution path for
this repository; any future such packaging must explicitly retain this patch.

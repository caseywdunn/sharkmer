# v3.2.1 maintenance release checks — 2026-11-09

The stable v3.2.1 candidate passes all eight engineering-check groups in Slurm
job **32671325**. The separately frozen HTTPS compatibility check passes in
job **32674065**. Package validation and publication remain pending. The focused checklist is
[RELEASE_v3.2.1.md](../../../dev_docs/RELEASE_v3.2.1.md).

The only dependency updates from released v3.2.0 are Rustls 0.23.37 → 0.23.46
and rustls-webpki 0.103.13 → 0.103.14. The root package version changes to 3.2.1.
[dependency-metadata.diff](dependency-metadata.diff) records the exact changes
against `e914ea728fc5ffc253273560e3de2cfc9395893b`; production Rust sources,
panels, schemas, tests, and check scripts are unchanged. Release documentation
and recipe metadata are separate from this production comparison.

The check passes 283 Rust tests per backend, 213 Python tests without skips,
formatting and Clippy on both backends, 158 reference checks, ten schemas,
nine exports, 19 loader dry runs, and fixture sequence oracles. It ran on
Bouchet's `devel` partition with two CPUs and a 16 GiB limit, completing in
1 m 49 s with 859,676 KiB peak batch RSS. These are engineering-job resource
receipts, not an application performance comparison. Tool versions were Rust
1.98.1 and Python 3.14.3; this does not test the declared Rust 1.88 minimum.

The default-backend binary SHA-256 is
`0d57521d5c6060df3e954d366b72fdccb0b3e40f9532162b400d62b14a59109c`.
The [raw log](logs/slurm-32671325.out) is copied byte-for-byte, including its
original line endings. It binds Cargo metadata and the check runner before
and after execution. [source-files.sha256](source-files.sha256) records source,
fixture, panel, schema, and tool hashes after the job. This is a working-tree
receipt, not a release-commit identifier. [receipt.json](receipt.json) contains
the machine-readable outcome; [sacct.psv](sacct.psv) records scheduler status
and resource use in America/New_York time (UTC−05:00).

[checks.sh](checks.sh) is the exact batch launcher. The copied
[release-check](release-check) and [release_check.py](release_check.py) retain
the runner; their original execution context is the repository root. Compiled
binaries and dependency caches are not copied here. Verify the retained bundle
with `sha256sum -c SHA256SUMS` from this directory; verify source receipts by
running `sha256sum -c` on `source-files.sha256` from the repository root.

The [HTTPS rerun](https-v2/README.md) verifies a pinned live HTTPS panel over
TLS 1.3, agreement with its local-panel control for three reads, and rejection
of an untrusted certificate before any HTTP request. It uses the same binary
as the engineering suite. The [first attempt](https-attempt-1/README.md) failed
only the harness's version-output assertion before any network request; its
raw evidence and the unavailable original-script limitation are preserved.

The [exact-version advisory audit](advisory-audit/advisory-audit.json) found no
reviewed or unreviewed GitHub advisory matching Rustls 0.23.46 or webpki
0.103.14 at its recorded query time. Raw package and version queries are
retained with hashes; this is not a guarantee against unknown vulnerabilities.

This dependency patch does not rerun the
[32-invocation v3.2 comparison](../v3.2-final-check-20261108/RESULTS.md) or expand
its biological/performance claims. The user-accepted held-out and full-nine-panel
biological deferrals remain unfinished work under #129.

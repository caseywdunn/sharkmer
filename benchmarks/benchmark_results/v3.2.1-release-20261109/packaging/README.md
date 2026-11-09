# Local v3.2.1 Bioconda validation

Job **32677573 passed** on 2026-11-09: 5 m 53 s, two CPUs, and 8,092,596 KiB
peak batch RSS (about 7.72 GiB). Conda-build 26.9.1 with Rust 1.99.1 built the
published source archive for linux-64. All four installed-command smoke tests
passed: `--help`, `--version`, `--list-panels`, and `--cite`.

An independent check extracted the resulting package and executed its binary.
[package-version-verification.json](package-version-verification.json) preserves
the exact full banner, `sharkmer 3.2.1 (https://github.com/caseywdunn/sharkmer)`,
exit status, and package/binary hashes. The package SHA-256 is
`321586cb9823c3eefdf4017a3b38cfdd92009e1bf07f799fe787087b46b1ceeb`.
Its hash was verified before retaining this bundle; the 2,143,313-byte package
is not copied here.

The build receipt binds the final source archive, recipe, original submitted
batch/environment files, environment specification, and raw job log. The
receipt's environment path records the reused v3.2.0 build environment; its
byte-identical specification is copied here. Files prefixed `submitted-` are
the archived job inputs, rather than later edits of the working scripts.

The job disabled uploads. Local validation does not imply upstream builds,
other-platform support, or channel availability. It does not test Rust 1.88
or repeat biological validation. See [publication status](../PUBLICATION.md).
This subfolder has independent checksums; verify with `sha256sum -c SHA256SUMS`.

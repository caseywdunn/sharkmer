# Postrelease 4.0.0-dev integration check

Job **32677833 passed all eight release-check groups** on 2026-11-09 after
merging the v3.2.1 maintenance release with the development branch and restoring
version 4.0.0-dev. Rustls 0.23.46 and rustls-webpki 0.103.14 remain locked.
No production-source, panel, schema, test, or check-runner change occurs from
released v3.2.1; the tested checkout was an uncommitted integration merge.

The job completed in 1 m 54 s with two CPUs and 1,473,496 KiB peak batch RSS.
The batch limit was 16 GiB. The raw log binds Cargo metadata and the check
runner before and after execution and records the resulting binary hash.
These are integration-check receipts, not a new biological or performance
comparison. The existing held-out/all-panel deferrals remain unchanged.

This subfolder has independent checksums and leaves the frozen release bundle
untouched. Verify from this directory with `sha256sum -c SHA256SUMS`.

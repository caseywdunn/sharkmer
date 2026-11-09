# v3.2.0 publication and distribution follow-up

Updated 2026-11-09. This continuing receipt is separate from the frozen
[engineering-check bundle](README.md) and its `SHA256SUMS`. It records later
publication events without changing the original logs or source receipts.

[v3.2.0 is published](https://github.com/caseywdunn/sharkmer/releases/tag/v3.2.0) at commit
`e914ea728fc5ffc253273560e3de2cfc9395893b`. Release preparation was merged in
[PR #165](https://github.com/caseywdunn/sharkmer/pull/165), followed by the
[release merge, PR #166](https://github.com/caseywdunn/sharkmer/pull/166).
[Tag CI passed](https://github.com/caseywdunn/sharkmer/actions/runs/42625940326). The `v3` maintenance branch was fast-forwarded to the
release commit.

The final GitHub source archive used by the local Bioconda recipe has SHA-256
`858dac22ea1d05a98a4917706de95e1012385e3d10078f42f7f1f98b0d883e0b`.
The recipe remains at version 3.2.0. Active Cargo development has advanced to
4.0.0-dev, with an empty Unreleased changelog section and the unfinished v4+
plan retained in [PLAN.md](../../../dev_docs/PLAN.md).

| Follow-up | Status |
| --- | --- |
| Local Linux Bioconda package build and installed-command tests | **Passed**, job 32662859: 12 m 47 s, two CPUs, about 7.64 GiB peak RSS; Rust 1.99.1 and all four smoke tests. [Packaging receipts](packaging/README.md). |
| Upstream Bioconda recipe pull request | [#71407](https://github.com/bioconda/bioconda-recipes/pull/71407) submitted; upstream CI/review and channel publication pending |
| Rustls advisory [#167](https://github.com/caseywdunn/sharkmer/issues/167) | Prompt v3.2.1 maintenance follow-up: newly surfaced alert #11, GHSA-2mjx-qc3c-rqvc, applies to the enabled TLS 1.3 path. |
| Independent held-out biological registration and full nine-panel campaign, [#129](https://github.com/caseywdunn/sharkmer/issues/129) | Explicitly deferred from v3.2, not completed; required before further assembly-policy tuning or broad accuracy claims |

Local package validation and upstream submission are complete. Upstream
CI/review and channel availability remain pending. The release retains the scoped biological limitations in
[REVIEW_v3.2.md](../../../dev_docs/REVIEW_v3.2.md).

The new Rustls alert #11 (GHSA-2mjx-qc3c-rqvc) surfaced after publication and
applies to Sharkmer's enabled TLS 1.3 path. It concerns malformed-handshake
handling; transcript authentication remains intact. Issue #167 tracks prompt
v3.2.1 maintenance. The other three reviewed advisories require configurations
absent from the current program; that assessment does not exempt alert #11.

The returned `4.0.0-dev` checkout also passes all eight release-check groups
in Slurm job 32665891. Its [log](postrelease/slurm-32665891.out) and
[accounting](postrelease/sacct.psv) are retained separately. Final PR/ancestry
verification is recorded on [release issue #164](https://github.com/caseywdunn/sharkmer/issues/164).

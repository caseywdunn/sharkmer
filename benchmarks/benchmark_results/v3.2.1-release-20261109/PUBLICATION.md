# v3.2.1 publication and distribution follow-up

Updated 2026-11-09. This continuing receipt is separate from the frozen
[engineering/HTTPS/advisory bundle](README.md) and its `SHA256SUMS`.

[v3.2.1 is published](https://github.com/caseywdunn/sharkmer/releases/tag/v3.2.1) at commit
`ca7fc4c9d35ae57c767186c475b029358af8c165` through
[PR #169](https://github.com/caseywdunn/sharkmer/pull/169).
[Tag CI passed](https://github.com/caseywdunn/sharkmer/actions/runs/42652146001); the `v3` maintenance branch has advanced to the release.
The v3.2.0 tag remains unchanged.

The final source archive for the v3.2.1 recipe has SHA-256
`e8031a7b346303f05e553886d473781b01a1a9c0dda0752d5b4573af6ad6061c`.
The 4.0.0-dev integration checkout retains Rustls 0.23.46 and webpki 0.103.14;
its [postmerge release checks](postrelease/README.md) pass all eight groups in
job 32677833. The completed security patch under
[#167](https://github.com/caseywdunn/sharkmer/issues/167) is removed from active
implementation work in [PLAN.md](../../../dev_docs/PLAN.md).

| Follow-up | Status |
| --- | --- |
| Local v3.2.1 Bioconda build and four installed-command tests | **Passed**, job 32677573; [packaging and independent version receipts](packaging/README.md). |
| Update [Bioconda PR #71407](https://github.com/bioconda/bioconda-recipes/pull/71407) to v3.2.1 | Submitted at recipe commit `2fe50affda66c4504bbdd36741889970cd7ceff3`; upstream CI/review and channel availability remain pending |
| 4.0.0-dev postmerge release checks | **Passed**, job 32677833; [receipts](postrelease/README.md). Final integration commit/PR and ancestry verification are recorded on [#167](https://github.com/caseywdunn/sharkmer/issues/167). |
| Upstream Bioconda builds/review and channel availability | Pending |

The source release does not imply channel availability. The biological
registration and all-nine-panel campaign remain deferred under #129, not passed;
this dependency maintenance does not expand the biological or performance
claims in the v3.2 review.

The [post-publication security receipt](publication-security/security-status.json)
reports zero open Dependabot alerts at 2026-11-09 16:43:18 UTC; alerts #8–#11
are explicitly fixed at 16:10:58 UTC. This records the observed repository
state without changing the frozen advisory audit.

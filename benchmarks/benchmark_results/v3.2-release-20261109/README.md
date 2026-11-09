# v3.2.0 release checks — 2026-11-09

The unchanged development candidate and stable-version candidate both pass all
eight release-check groups. The stable rerun includes one equivalent nested-`if`
rewrite requested by Clippy after correcting the declared minimum Rust version.
Packaging, merge/tag CI, release assets, and publication are still pending.

| Job | Candidate | Outcome | Elapsed | Batch peak RSS |
| --- | --- | --- | --- | --- |
| 32654564 | `3.2.0-dev`, unchanged production sources | All eight groups pass | 1 m 37 s | 1,534,548 KiB |
| 32657614 | `3.2.0`, declared Rust minimum 1.88 | Stops at `clippy::collapsible_if`; later groups not run | 8 s | 870,120 KiB |
| 32659082 | `3.2.0`, equivalent `if` rewrite | All eight groups pass | 1 m 33 s | 897,796 KiB |

Each job used the Bouchet `devel` partition, two CPUs, and a 16 GiB limit.
These are engineering-check resource receipts, not application performance
comparisons. Slurm times in [sacct.psv](sacct.psv) are America/New_York (UTC−05:00).
The full logs preserve the failed attempt as well as both successful runs.

## Checks and source relationship

The eight groups are metadata, formatting, Clippy on both backends, Rust tests
on both backends, Python tests, reference audit, panels, and fixtures. Each
successful run passes 283 Rust tests per backend (253 unit + 30 integration),
213 Python tests without skips, 158 public reference receipts, ten panel/schema
checks, nine embedded exports, and 19 loader dry runs. The fixture gate verifies
complete ordinary 18S and 28S_2 sequences; Rust tests also cover the counting
oracle. The Python log includes expected error messages from negative fixtures;
the suite concludes `OK`.

The base commit is `95f5390535690f0a9add083f48a5077590096e40`. The accepted
`3.2.0-dev` binary is reproduced exactly with SHA-256
`7ef94fedafec49300146c0acf52086915fa78a82a8786dc11413e399a339dd9a`.
For the stable candidate, Cargo package metadata changes to `3.2.0` and the
declared Rust minimum changes from 1.85 to 1.88. That declared minimum enables
Clippy's suggestion to combine `if let Some(count)` with the same `count >=
min_count` guard. The lookup, comparison, and admitted-candidate writes keep
their order and conditions. [production.diff](production.diff) records these
changes; no panel, test, schema, dependency version, or other production source
changes occur relative to the accepted candidate.

The stable default-backend binary SHA-256 is
`cf4d8f4edf28fddd06e064d6c8d77a8450189a2710c159e584d0a34ebe59652d`.
Both successful logs bind the Cargo manifest, lockfile, and check runner before
and after execution. [source-files.sha256](source-files.sha256) captures the
production, test, panel, schema, tooling, workflow, and recipe files after the
final run; it is a working-tree receipt, not a release-commit identifier.
The logs use Rust 1.98.1 and Python 3.14.3; they do not test Rust 1.88 itself.

## Scope and retained artifacts

These finite engineering gates do not rerun the 32-invocation biological and
performance comparison. That accepted evidence remains in the
[focused final-check bundle](../v3.2-final-check-20261108/RESULTS.md). The version
metadata and equivalent guard rewrite do not introduce an assembly-policy
change. The held-out registration and full nine-panel biological campaign
remain explicitly deferred, not passed; no broader biological accuracy or
performance claim follows from this bundle.

[receipt.json](receipt.json) contains machine-readable outcomes and binary
identities. [checks.sh](checks.sh) is the exact batch launcher;
[release-check](release-check) and [release_check.py](release_check.py) retain
the runner. [toolchain-activate.sh](toolchain-activate.sh) records the local Rust
activation paths. These copied scripts are provenance; run the repository's
`scripts/release-check --offline` from its root with the documented dependencies.
The binary and shared dependency caches are not copied into this small bundle.

Verify bundle files from this directory with `sha256sum -c SHA256SUMS`. To check
the bound source tree, run `sha256sum -c` on `source-files.sha256` from the
repository root. The scheduler receipt was retrieved using:

```sh
sacct -j 32654564,32657614,32659082 --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,Start,End --parsable2
```

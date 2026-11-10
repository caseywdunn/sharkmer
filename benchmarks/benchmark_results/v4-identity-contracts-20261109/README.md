# v4 identity/extraction correctness — 2026-11-09

Slurm job 32722125 passed all seven selected development-check groups on
frozen sources: formatting, Clippy and Rust tests on both backends, Python
validation, reference audit, panels, and the pinned sequence fixture.
There are 286 Rust tests per backend (256 unit + 30 integration), 213 Python
tests, 158 audited references, ten schemas, nine exports, and 19 loader checks.
The new scalar oracle independently checks positioned rolling extraction.

The two-CPU/16-GiB job completed in 1m03s; batch peak RSS was 860508 KiB.
These are check-job resources, not an application performance comparison.
Rust 1.98.1 and Python 3.14.3 were used. The raw log binds changed extraction
sources, Cargo files, and the checker before and after execution. The default
binary SHA is in receipt.json. Source receipts include all tracked production,
test, panel, schema, and script files plus the new identity module.

The retained launcher uses repository-local toolchain/cache paths; it is a
receipt of this run. Use the normal project checks for a new environment.
Run `sha256sum -c SHA256SUMS` here to verify bundle files, and check
`source-files.sha256` from the repository root to verify bound sources.

This is the first increment of #139. The exact-store, replay, resource-budget,
and backend integration gates remain open in the active plan. The existing
biological-validation deferrals remain in force. No release ritual was started.

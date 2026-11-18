# Direct-ingestion candidate verification — 2026-11-11

Slurm job 32945384 passed the normal development checks: 310 Rust tests per
backend (278 unit and 32 integration), 213 Python tests without skips, both
Clippy configurations, formatting, 158 reference audits, ten schemas, nine
embedded exports, nineteen loader dry runs and the pinned full-sequence fixture.
The before/after source observations match the retained source manifest.

The new CLI comparisons require complete fixture FASTA agreement between legacy
and exact counting at 1/2/4/8 threads and through one-shot stdin. A separate
all-key/count fingerprint agrees between the two backend implementations.
Scalar-oracle tests cover direct ingestion across N boundaries and more than
one legacy batch boundary without retaining sequence batches. Additional tests
cover record IDs, qualities, sources, mate designation/order, caps, CRLF, uneven
paired EOF selection, multi-member gzip corruption, buffer reuse, truncation,
and count/record limits failing without publishing partial FASTA.

Exact counting updates valid-base and occurrence totals while streaming,
then moves its table into assembly without consolidation. Shared reservations
cover count buffers and record arenas, including simultaneous resize buffers.
These are owned-buffer measurements; decoder/transport buffers, process overhead,
graphs and replay are not yet covered by a complete run limit.

The validated default-build executable SHA is
`2e050adbc3018058711d5e992090c78faa880401f2392df975332525b558f255`.
The comparison path remains behind hidden development controls. Default
selection still requires the frozen quantitative performance gates; job
32947879 is evaluating fresh paired comparisons separately. No speedup or
biological-accuracy claim is established by these correctness checks.

Verify `SHA256SUMS` here and `source-files.sha256` from the repository root at
this increment. The Slurm batch peak includes compilation and helper processes;
it is not per-application RSS evidence. The launcher retains local toolchain
and dependency paths as execution provenance. The benchmark harness separately
passes Ruff 0.17.4 formatting and lint checks.

# Roadmap

This document describes the planned development of sharkmer across upcoming
major releases.

## v2.0 — Cleanup and polish

The goal of v2.0 is to improve code quality, usability, and developer
infrastructure without changing core algorithms or results. This release
may change CLI flags, output format, and error behavior, so it should ship
quickly after the v1.0 manuscript to avoid workflows being built against
patterns that will soon change.

### Developer infrastructure

- Fix existing test compilation errors
- Add CI via GitHub Actions (fmt, clippy, test across platforms and hash backends)
- Add CLAUDE.md for AI-assisted development context
- Add a workspace Cargo.toml at the repo root so `cargo test` works from the
  top level
- Add test fixtures derived from synthetic data used in the manuscript
- Expand integration test coverage

### Error handling

- Replace panics and `.unwrap()` in library code with `Result` types
- Consolidate all exit-with-error logic in `main()`
- Use `anyhow` or `thiserror` for error context
- Return meaningful exit codes

### Logging and output

- Adopt `log` + `env_logger`, replacing ad-hoc `println!` and the
  `--verbosity` integer
- Separate data output (stdout) from logs/progress (stderr)
- Support `-v`/`-vv`/`-vvv` and `--quiet` flags
- Produce structured (JSON or clean TSV) stats output

### CLI improvements

- Native gzip input support via `flate2` (replaces `zcat | sharkmer` pattern)
- Rename `-n` to `--chunks`
- Replace `--pcr panels` with `--list-panels`
- Replace `--verbosity N` with `-v` stacking
- Add `--pcr-file` flag to load primer panels from a TSV file, enabling users
  to share and reuse primer definitions without code changes

### Performance (low-hanging fruit)

- Use a `HashMap<u64, NodeIndex>` for node lookup during graph construction
  instead of linear scan
- Avoid storing reverse complements in the PCR kmer hash table; look up both
  orientations at query time
- Pre-size hash maps based on expected data size
- Evaluate `ahash` as a drop-in replacement for `FxHash`

### Code cleanup

- Deduplicate `is_valid_nucleotide` (defined in both `main.rs` and
  `pcr/mod.rs`)
- Deduplicate `KmerCounts` impl blocks across hash backends using a trait
  or macro
- Remove dead code (`ParameterValue` enum, unused functions)
- Fix `bp_length_to_kmer_length` comparison of `usize <= 0`

## v3.0 — Graph traversal

The goal of v3.0 is to replace the current ad-hoc graph extension and pruning
heuristics with principled algorithms from the assembler literature. This will
change results but should have minimal impact on CLI usage or output format.

Target improvements:

- Better recovery of single-copy nuclear genes at lower coverage
- Replace heuristic ballooning detection with coverage-aware graph
  simplification (tip clipping, bubble popping)
- Replace backward-degree-based termination with principled traversal
- Improved handling of repeats (current cycle avoidance is too aggressive)
- Coverage-weighted best-path algorithm to replace `all_simple_paths`
  enumeration
- Better path scoring beyond `kmer_min_count` ordering
- Paired-end read awareness
- Preparation for v4.0 metagenomics support

## v4.0 — Metagenomics

The goal of v4.0 is to reliably produce multiple distinct products from mixed
samples.

Target improvements:

- Deconvolution of mixed-organism samples based on coverage profiles
- Support for reporting multiple distinct amplicons per gene
- Per-product coverage profiles to distinguish true variants from assembly
  artifacts
- Configurable limits for number of amplicons per gene
- Validation against known metabarcoding datasets

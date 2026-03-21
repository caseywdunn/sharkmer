# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - v2.0.0-alpha

Development in progress on the `dev` branch. See [ROADMAP.md](ROADMAP.md) and
the [issue tracker](https://github.com/caseywdunn/sharkmer/issues?q=label%3Av2.0)
for the full list of planned changes.

### Added

- ROADMAP.md with plans for v2.0, v3.0, and v4.0 releases
- CHANGELOG.md

## [1.0.1] - 2025-01-23

### Fixed

- Added missing MIT LICENSE file
- Fixed build.sh path for bioconda build environment

## [1.0.0] - 2025-01-22

### Added

- Initial stable release
- Kmer counting with incremental histogram generation
- In silico PCR (sPCR) with seeded de Bruijn graph assembly
- Preconfigured primer panels: cnidaria, human, teleostei, angiospermae,
  insecta, bacteria, metazoa
- Manual primer pair specification via `--pcr` flag
- Support for reading FASTQ from files or stdin
- Configurable kmer length, chunk count, thread count
- Multiple hash map backends via feature flags (fxhashmap, intmap, nohashmap)
- Bioconda recipe for binary distribution
- sharkmer_viewer Python tool for histogram visualization

[Unreleased]: https://github.com/caseywdunn/sharkmer/compare/v1.0.1...HEAD
[1.0.1]: https://github.com/caseywdunn/sharkmer/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/caseywdunn/sharkmer/releases/tag/v1.0.0

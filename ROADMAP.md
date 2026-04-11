# Roadmap

This document describes the planned development of sharkmer across upcoming
major releases. See the [issue tracker](https://github.com/caseywdunn/sharkmer/issues)
for detailed specifications.

## v2.0 — Cleanup and polish (released)

Released. See [CHANGELOG.md](CHANGELOG.md) for details.

## v3.0 — Graph traversal (released)

Released. See [CHANGELOG.md](CHANGELOG.md) for details.

Key improvements: bidirectional graph extension, frontier-queue extension
eliminating O(n²) node scanning, reachability pruning replacing ad-hoc
heuristics, coverage-aware tip clipping, read threading with bubble
resolution, composite path scoring, mismatch-aware primer kmer cap,
dynamic node budget, read caching, and panel versioning and validation
infrastructure. New `c_elegans` panel.

Paired-end phasing (#101) is deferred to v4.0: the infrastructure
(`PairedEndLink`, `thread_reads_paired()`) is built but downstream
consumption is not yet wired.

## v4.0 — Metagenomics

The goal of v4.0 is to reliably produce multiple distinct products from mixed
samples.

Target improvements:

- Paired-end phasing (#101): wire `PairedEndLink` data into bubble
  resolution and path scoring (infrastructure built in v3.0)
- Read-supported graph pruning (#99): prune edges/nodes lacking read
  support to reduce chimeric joins
- Deconvolution of mixed-organism samples based on coverage profiles
- Support for reporting multiple distinct amplicons per gene
- Per-product coverage profiles to distinguish true variants from assembly
  artifacts
- Configurable limits for number of amplicons per gene
- Validation against known metabarcoding datasets
- Document QIIME2/SILVA workflow and consider convenience flags (#68)

## v5.0 — New targeting paradigms

The goal of v5.0 is to generalize beyond the PCR metaphor. The paired-primer
sPCR workflow remains a primary feature, but v5.0 adds additional targeting
modes built on the underlying seeded assembly engine.

### Oligo-based extensions

- Single-primer (forward-only) seeded assembly — relax the requirement for a
  reverse primer, extending outward from a single seed until coverage drops or
  a length limit is reached
- Synthetic ddRADseq — seed with restriction enzyme recognition sites to
  produce in silico RAD-like fragments
- Synthetic targeted enrichment — seed with UCE (ultraconserved element) oligo
  probe sequences to recover UCE loci and their flanking regions from genome
  skimming data

### Exotic targeting

Novel approaches to isolating homologous genome regions across samples that
have no parallel in benchtop methods but are natural extensions of seeded
assembly:

- Seed generation from reference sequences, alignments, or gene models
  (may require helper scripts or a companion tool)
- Cross-sample kmer comparison to identify conserved seeds without a reference
- Iterative seeding where products from one round seed the next

### Architecture considerations

- The targeting mode (paired-primer, single-oligo, restriction site, UCE probe)
  should be a parameter of the panel/primer YAML format introduced in v2.0,
  not a separate code path — the seeded assembly engine is shared
- v3.0 graph traversal improvements and v4.0 multi-product support are
  prerequisites — single-primer extension and UCE enrichment both produce
  variable-length products that need robust graph handling

# Roadmap

This document describes the planned development of sharkmer across upcoming
major releases. See the [issue tracker](https://github.com/caseywdunn/sharkmer/issues)
for detailed specifications.

## v2.0 — Cleanup and polish

The goal of v2.0 is to improve code quality, usability, and developer
infrastructure without changing core algorithms or results. This release
changes CLI flags, output format, and error behavior, so it should ship
quickly after the v1.0 manuscript to avoid workflows being built against
patterns that will soon change.

### Developer infrastructure

- Fix existing test compilation errors (#17)
- Add CI via GitHub Actions (#19)
- Add CLAUDE.md (#18)
- Add workspace Cargo.toml at repo root (#20)
- Add test fixtures from ERR571460 (#30)
- Expand integration test coverage (#31)
- Regression benchmarking suite with real-world SRA datasets (#50)
- CHANGELOG.md

### Bug fixes

- Fix unsigned integer underflow risks (#55)
- Fix edge cases: empty inputs, degenerate primers, division by zero (#60)
- Ensure deterministic output across runs (#38)
- FASTQ parser does not validate format (#57)
- Deduplication may silently discard real variants (#56)
- Ensure FASTA output uses standard 80-char line wrapping (#42)

### Error handling

- Replace panics and unwraps with proper Result types (#21)
- Remove hardcoded singleton filtering before sPCR (#52)
- Rename misleading `get_canonical` method (#58)

### Logging and output

- Adopt `log` + `env_logger`, separate stdout/stderr (#22)
- Structured YAML stats output with PCR results (#32)
- Improved FASTA header format with key=value metadata (#34)
- Print concise summary line at completion (#43)
- Progress indicators for long-running steps (#37)
- Polished CLI output formatting (#49)
- Add `--color auto|always|never` and respect `NO_COLOR` (#44)
- Report peak memory usage (#47)

### CLI improvements

- Split `--pcr` into `--pcr-panel`, `--pcr-panel-file`, `--pcr-primers` (#24)
- Native gzip input support with auto-detection (#23)
- Default to no incremental kmer counting (`--chunks 0`) (#35)
- Add `--pcr-panel-file` for sideloading YAML primer panels (#25)
- Add `--list-panels`, `--export-panel`, `--help-pcr` (#24)
- Make `--sample` required, `--max-reads` optional (#24)
- Organize `--help` with grouped headings (#24)
- Add `--cite` flag (#40)
- Add `--dry-run` mode (#46)
- Generate shell completions (#39)
- Error when stdin is a terminal with no input files (#45)
- Warn when overwriting existing output files (#48)
- User-facing warnings for common mistakes (#59)
- Validate input early before ingestion (#41)

### Primer panels

- Refactor panels from Rust code to YAML files via `include_str!()` (#33)
- YAML format with panel-level metadata and nested primers (#33)
- Shared parser for built-in and user-sideloaded panels (#25, #33)
- Include a `mode` field per primer (default `paired-primer`) to prepare for
  v5.0 targeting modes. Only `paired-primer` accepted in v2.0; unknown modes
  produce a clear error. (#33)

### Performance

- HashMap for node lookup in graph construction (#27)
- Avoid storing reverse complements in PCR kmer table (#28)
- Eliminate double hash lookups in graph extension (#54)
- Free chunk hash tables after merging (#53)
- Pre-size hash maps, evaluate ahash (#29)
- Parallelize sPCR across genes with rayon (#36)

### Code cleanup

- Deduplicate code: `is_valid_nucleotide`, `KmerCounts` impls, dead code (#26)
- Update sharkmer_viewer for new histogram/stats format (#61)

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
- Preparation for v4.0 metagenomics support

### Read-based graph refinement

Use reads (not just kmers) to improve the assembly graph. Currently reads
are discarded after kmer extraction — only the kmer count table is used for
graph construction. Threading reads back through the graph enables:

- **Graph pruning**: Edges and nodes not supported by any full read path are
  likely artifacts. Read support can distinguish real low-frequency paths
  from chimeric joins, improving results at shallow coverage.
- **Bubble resolution**: When two paths diverge and reconverge (heterozygous
  sites, sequencing errors), reads that span the bubble can resolve which
  path(s) are real.
- **Coverage estimation per path**: Read threading gives per-edge coverage
  that reflects actual read support, not just kmer frequency. This is more
  informative for path scoring than raw kmer counts.
- **Paired-end constraints**: When mate pairs are available, the insert size
  distribution constrains which paths are physically plausible. A read pair
  where R1 maps to one branch and R2 maps to another implies those branches
  are linked within the insert distance.

Implementation approach: two-pass over input data. Pass 1 counts kmers and
builds the graph (as now). Pass 2 re-reads the FASTQ and threads each read
through the graph, annotating edges with read support. For `--sra` input,
the second pass re-streams from ENA. For file input, seeks back to the
start. For stdin, either buffer reads in memory or require file input when
read threading is enabled.

Paired-end awareness requires knowing which files are R1 vs R2. The v2.0
convention (positional files in R1, R2 order; `--sra` downloads in order)
is sufficient — v3.0 can infer pairing from file order, with a `--paired`
flag to explicitly enable pair-aware refinement.

### Graph construction efficiency

- Build a single graph seeded with all forward primer kmers simultaneously,
  instead of building a separate graph per forward primer kmer. Currently
  10 forward primer kmers means 10 nearly-identical graphs built from
  scratch — the largest performance waste in the tool.
- Extend graphs incrementally across coverage threshold steps instead of
  rebuilding from scratch at each threshold. The lower-threshold graph is a
  strict superset of the higher-threshold one.
- Incremental histogram updates after each chunk merge, instead of
  reiterating the full consolidated kmer table each time.

### Kmer pipeline optimizations (deferred from v2.0)

- Extract kmers directly from ASCII sequence in a single pass, eliminating
  the Read struct encoding/decoding round-trip. The 2-bit `Read` encoding was
  designed for storing reads in memory, but reads are not stored — kmers are
  extracted and the read is discarded. A direct ASCII-to-kmer sliding window
  avoids the intermediate encoding step entirely.
- Use `u32` for kmer counts instead of `u64` — counts rarely exceed 4 billion,
  and this saves ~33% memory in the hash table. Important for low-coverage
  work where every byte of hash table capacity matters.

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

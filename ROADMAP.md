# Roadmap

This document describes the planned development of sharkmer across upcoming
releases. See [dev_docs/PLAN.md](dev_docs/PLAN.md) for execution order, linked
issues, architecture decisions, and acceptance gates. Version numbers below
are planning targets, not promised dates.

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

Paired-end phasing (#101) is deferred to the read-evidence release: the infrastructure
(`PairedEndLink`, `thread_reads_paired()`) is built but downstream
consumption is not yet wired.

## v3.2 — Correctness and trustworthy validation

Fix gzip truncation, premature threshold stopping, collapsed-repeat output,
threading orientation/gap/selection errors, stale outputs, primer allocation
bounds, and cache ownership/concurrency. Repair all-product validation and
benchmark provenance before using them to select new algorithms.

## v4.0 — Scalable exact counting

Build a one-pass-over-original-input, allocation-light, parallel exact counter.
Benchmark compact tables while preserving abundance. Add external counting,
bounded-memory final random lookup, and streaming replay/spooling so deeper
datasets can run on a 16 GB laptop with appropriate temporary storage.

Retain incremental k-mer counting as an isolated legacy analysis path through
v4.x, including explicit `--chunks` use and histogram compatibility. Remove its
chunk/snapshot/merge overhead from normal sPCR. Full removal would require a
separate decision and migration plan.

## v4.1 — Low-coverage nuclear recovery

Improve primer discovery and bounded target graphs; compact nonbranching paths
where measurements justify it. Stream evidence against batches of target
graphs (#116), consume paired-end constraints (#101), and use conservative
read-supported pruning (#99). Add evidence-based singleton rescue and explicit
partial/ambiguous outcomes. Evaluate local multi-k reconstruction as an
experiment after the single-k baseline is sound.

## v4.2 — Metagenomic diversity

Preserve rare templates through seed, threshold, and path search; use exact
deduplication by default in an explicit diversity mode. Report supported
haplotypes, phase uncertainty, search truncation, and defensible abundance.
Validate every product against known mixtures and negative controls, measuring
rare-template recall and chimeras as well as sequence precision. Document
downstream classification workflows (#68).

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
- v3.0 graph traversal improvements and v4.2 supported multi-product recovery are
  prerequisites — single-primer extension and UCE enrichment both produce
  variable-length products that need robust graph handling

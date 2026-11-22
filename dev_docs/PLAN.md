# Development plan: scalable, reliable sPCR

**Active version: 4.1.0-dev. Released baseline: [v4.0.0](https://github.com/caseywdunn/sharkmer/releases/tag/v4.0.0).**

Start v4.1 with independent held-out registration and the nine-panel
biological evaluation under [#129](https://github.com/caseywdunn/sharkmer/issues/129).
Use the [v4.0 release plan](https://github.com/caseywdunn/sharkmer/blob/v4.0.0/dev_docs/PLAN.md)
and [final comparison](https://github.com/caseywdunn/sharkmer-benchmarks/tree/main/benchmarks/benchmark_results/v4-final-integration-20261117)
for the shipped counting and resource gates, historical checkpoints, and
comparison evidence. The [v3.2 release assessment](REVIEW_v3.2.md) records
the earlier biological evidence and deferrals.

Keep explicit incremental counting and its histogram workflow supported
through v4.x. Preserve the standard sPCR baseline as new read-evidence
policies are evaluated. Use the [resource contracts](COUNT_STORE_CONTRACTS.md)
and [sequence annotation conventions](SEQUENCE_ANNOTATIONS.md) for those changes.

The v4.0.0 source asset provides a checksum-pinned URL for Bioconda.
Protect existing package source pins before the history cutover below.

## Repository size and benchmark storage

The [public companion repository](https://github.com/caseywdunn/sharkmer-benchmarks)
holds browseable reports and protocols, including the final v4 verification.
Its [historical archive](https://github.com/caseywdunn/sharkmer-benchmarks/releases/tag/legacy-2026-11)
preserves all 1,538 result files at `7bd1ac3` and a benchmark-only Git bundle
covering both historical results directories. The downloaded archives matched
their checksums. The original Git bundle and unpublished result files are
backed up separately in the development workspace.

The cleanup keeps the benchmark tools in production and gives their unit tests
a small synthetic fixture. Removing the committed results reduces the tracked
checkout from 252 MiB to about 21 MiB. The 32-ref rehearsal removed
`benchmarks/benchmark_results/` and its earlier name, `benchmarks/results/`.
Every other file matched at each branch and tag tip. The Git pack fell from
187 MiB to 32.8 MiB. This measures branch and tag history in the rehearsal.

The companion holds [the original source tarballs](https://github.com/caseywdunn/sharkmer-benchmarks/releases/tag/original-sources-2026-11)
for all seven pre-v4 release tags. Before applying the rewrite on GitHub,
protect package source URLs. Bioconda
3.2.1 shipped via [PR #71402](https://github.com/bioconda/bioconda-recipes/pull/71402).
It pins the v3.2.1 generated tarball at SHA-256
`e8031a7b346303f05e553886d473781b01a1a9c0dda0752d5b4573af6ad6061c`.
The current recipe needs to point to a stable copy of those exact bytes
before the v3.2.1 tag can move. Open PR #71407 is a separate recipe update.
Review affected release tags, downstream pins and the ref map together before
cutover. Developers will need to re-clone or transplant their work onto the
rewritten history. The original Git bundle provides a rollback point.

**History cutover gate: not met.** The rehearsal and archive checks are
complete. The published tags and live package source pin still need the
compatibility step before a coordinated remote rewrite.

## v4.1 — Read evidence and low-coverage nuclear recovery

- [ ] [#129](https://github.com/caseywdunn/sharkmer/issues/129) — Complete the
  held-out registration and full nine-panel biological evaluation deferred
  from v3.2, before further assembly-policy tuning or broad accuracy claims.
  Preserve independent truth, callable targets, negative controls, and explicit
  unavailable/unassessed outcomes; regression success is not a substitute.
  Obtain an independently annotated Agalma ITS source without filling the gap
  with assembled products or misannotated sources. Review gene annotation
  independently of public sequence origin, and add independent same-sample
  truth and read-backed assessment of divergent/repeat regions. Cross-individual
  differences alone do not establish incorrect products.

- [ ] [#146](https://github.com/caseywdunn/sharkmer/issues/146) — Batch primer discovery without combinatorial variant enumeration.
- [ ] [#147](https://github.com/caseywdunn/sharkmer/issues/147) — Build length-bounded target graphs and compact nonbranching paths.
- [ ] [#148](https://github.com/caseywdunn/sharkmer/issues/148) — Use read-spanning branch histories to resolve repeats and constrain paths.
  Integrate bounded quality-preserving replay and positioned histories under
  #145/#116 behind an experimental opt-in path, with feature-off compatibility
  against its frozen baseline. Specify ambiguity-preserving variant-aware
  mapping before reconstructing alternatives. Marginal coverage does not
  justify repairing missing edges or relaxing repeat gates. Admission changes
  require cause-specific positive evidence, repeat/mixture controls, and fresh
  released-baseline high-copy comparisons that retain historical v3.1 identity
  checks.
- [ ] [#149](https://github.com/caseywdunn/sharkmer/issues/149) — Add evidence-based low-coverage nuclear recovery and partial outcomes.
- [ ] [#116](https://github.com/caseywdunn/sharkmer/issues/116) — Stream reads through bounded batches of target graphs; replace the
  old per-gene full-input scan proposal and all-read RAM vector.
- [ ] [#101](https://github.com/caseywdunn/sharkmer/issues/101) — Apply mate orientation and insert-size constraints to paths.
- [ ] [#99](https://github.com/caseywdunn/sharkmer/issues/99) — Use conservative evidence-aware pruning; lack of evidence is not
  automatically evidence of an erroneous edge.
- [ ] [#106](https://github.com/caseywdunn/sharkmer/issues/106) — Add explicit primer alternatives with current CLI/schema semantics
  and selection that does not erase minority alternatives.
- [ ] [#127](https://github.com/caseywdunn/sharkmer/issues/127) — Expose suggested-k guidance without silently switching global k.
- [ ] [#114](https://github.com/caseywdunn/sharkmer/issues/114) — Re-profile the proposed scoring optimization; implement only if
  it remains material after graph changes.

Dependency order: seeding/target graphs and replay → [#116](https://github.com/caseywdunn/sharkmer/issues/116) and [#101](https://github.com/caseywdunn/sharkmer/issues/101) →
linked-path evidence → low-coverage policy and [#99](https://github.com/caseywdunn/sharkmer/issues/99). [#106](https://github.com/caseywdunn/sharkmer/issues/106) and [#127](https://github.com/caseywdunn/sharkmer/issues/127) support
primer usability; [#114](https://github.com/caseywdunn/sharkmer/issues/114) is an optional measured optimization.

Maintain an explicit standard policy and add an opt-in low-coverage policy.
Use per-target/per-haplotype coverage, realistic intron spans, and
paralog/orthology evaluation. Admit singletons only under a specified evidence
model, not by globally relaxing all pruning. Return partial assemblies where
useful, with an explicit distinction from complete amplicons.

Evaluate local multi-k assembly on recruited/replayed reads after single-k
behavior is validated. Smaller k can restore overlap and larger k can resolve
some repeats; neither establishes unobserved bases or unsupported phase.
Do not make this experiment a requirement for shipping proven v4.1 gains.

### Release gate

- Interior reads and both strands contribute equivalent valid evidence.
- Evidence buffering is independent of total input length, with graph batch
  size and retained links bounded or spilled explicitly.
- Repeat lengths require read/fragment evidence; insufficient linkage yields
  ambiguity rather than a fabricated complete product.
- Publish held-out callable recall and sequence precision across nuclear depth,
  heterozygosity, intron length, and paralog similarity.
- Demonstrate that singleton rescue improves recall without an unreviewed
  increase in false products; physical gaps remain partial/unknown.
- Standard high-copy behavior remains a protected benchmark class.

## v4.2 — Metagenomic diversity

- [ ] [#150](https://github.com/caseywdunn/sharkmer/issues/150) — Preserve rare templates through sPCR seed, threshold, and path search.
- [ ] [#151](https://github.com/caseywdunn/sharkmer/issues/151) — Report supported haplotypes, uncertainty, and defensible target abundance.
- [ ] [#68](https://github.com/caseywdunn/sharkmer/issues/68) — Document downstream classifier/workflow integration using current
  flags and explicit limitations on marker resolution and organism abundance.

Add an explicit diversity policy. Explore plausible components and coverage
levels after a dominant product is found. Preserve exact distinct sequences by
default; optional similarity clustering is separate from assembly.

Use branch histories and mate constraints to preserve supported linkage.
If two distant variants have no distinguishing read/fragment evidence, report
phase blocks or ambiguity instead of enumerating possible combinations as
observed haplotypes. Estimate abundance from discriminating sequence and
compatible reads; do not equate shared marker coverage with organism counts.

### Release gate

- Recover both supported products in the 100:4 mixture with common primer ends.
- Evaluate 1 bp differences, indels, multiple similar templates, >40 seed
  candidates, and >20 genuine products under configured limits.
- Publish precision, rare-template recall, chimera rate, phasing accuracy, and
  abundance error stratified by target coverage and divergence.
- Every product carries evidence/completeness information; truncation is visible.
- No confident organism count or full haplotype is claimed when the marker or
  linkage is uninformative.

## Benchmark matrix and acceptance policy

Keep calibration and held-out evaluation separate. Freeze fixture identities,
truth, and target-specific tolerances before parameter sweeps.

- Existing rRNA/organelle successes, repetitive/AT-rich targets, and known
  failures; compare whole products rather than FASTA counts.
- Nuclear per-target/per-haplotype depths such as 1, 2, 3, 5, 10, 20, and 40x,
  including heterozygosity, introns, paralogs, and absent targets.
- Known metagenomic mixtures spanning abundance ratios and sequence divergence;
  vary absolute minority depth separately from abundance ratio.
- Error-free controls, realistic quality/error profiles, primer mismatches,
  physical gaps, and deliberately unresolvable repeat/phase cases.
- Plain/gzip/multi-member gzip, paired/interleaved/single-end, local/replayed/
  one-shot inputs; multiple k values and 1/2/4/8 threads where supported.
- Small exact-oracle tests, the included fixture, approximately 1M-read common
  workloads, and deep datasets that exceed the in-memory budget.

Record decompression/parsing, emission, table updates, merge/finalization,
seed discovery, graph building, threading, search, total wall time, peak RSS,
allocator heap separately, disk usage, and count-store lookup performance.
Benchmark without graph dumps for ordinary runtime; measure diagnostics
separately. Do not promise a speedup before the relevant measurements exist.

Required correctness gates are absolute. Biological improvements need reviewed
precision/recall thresholds on callable truth sets; sequencing depth alone
cannot guarantee complete recovery. Performance changes require reproducible
end-to-end improvement, not just a faster isolated hash benchmark.

## Design references

These are comparison/design references, not dependencies on their software.

- [Jellyfish: exact parallel counting with a specialized shared table](https://pmc.ncbi.nlm.nih.gov/articles/PMC3051319/).
- [KMC 3: counting and manipulating k-mer statistics](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399).
- [Linked de Bruijn graphs: retaining read-scale connectivity](https://pubmed.ncbi.nlm.nih.gov/29554215/).
- [Bifrost: compacted graph construction and indexing](https://pmc.ncbi.nlm.nih.gov/articles/PMC7499882/).
- [metaSPAdes: local coverage and metagenomic assembly challenges](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

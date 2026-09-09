# Development plan: scalable, reliable sPCR

Status: v3.2 implementation in progress following the 2026-09-08 review of `5a66468`
(Sharkmer 3.1.0). Tracking issue: [#152](https://github.com/caseywdunn/sharkmer/issues/152).

This is the active execution plan. [ROADMAP.md](../ROADMAP.md) describes the
release sequence; [PLAN_v3.md](PLAN_v3.md) preserves the historical v3 plan.
Check items off only when code, validation, documentation, and acceptance
criteria are complete. Issue bodies carry implementation scope and dependencies.

## Release strategy

Ship several independently useful releases rather than hold correctness fixes
until metagenomic inference is finished. Version numbers are planning targets,
not promised dates. Opt-in additions can ship as v4.x minor releases; reassess
the major version if an implementation requires incompatible behavior or formats.

| Release | Primary result | Boundary |
| --- | --- | --- |
| v3.2 | Correct existing sPCR behavior and establish trustworthy validation | No new assembly accuracy claims based only on more FASTA products |
| v4.0 | Faster exact counting and complete bounded-memory ingestion/lookup | Existing standard sPCR remains the biological baseline |
| v4.1 | Read-supported low-coverage nuclear recovery | Explicit partial/ambiguous outcomes; optional local multi-k remains experimental until validated |
| v4.2 | Recovery of supported metagenomic template diversity | Rare-template retention, phase uncertainty, and defensible abundance reporting |
| v5.0+ | Additional targeting modes | Single oligos, restriction sites, UCE probes, reference/profile seeds; separate future scope |

Increasing input capacity in v4.0 should already help some nuclear targets.
The new low-coverage policy and linkage-dependent sensitivity belong in v4.1.
Do not delay v3.2 waiting for a custom counter, unitigs, or a new inference model.

## Incremental k-mer counting decision

**Retain the feature as an isolated legacy analysis path through v4.x.
Remove its machinery from the normal sPCR counting path.**

- Preserve explicit `--chunks > 0` use, existing histogram formats, and the
  existing viewer workflow. Maintenance and correctness fixes continue;
  new rarefaction features are outside this program.
- `--chunks 0` must create one logical exact count store, without chunk
  distribution, incremental histogram snapshots, or copying a sole table.
- Share byte parsing/encoding where useful, but do not require every counter
  backend to support incremental snapshots or legacy chunk merging.
- Preserve combined legacy histogram/PCR invocation through an adapter to the
  finalized counts. Document unsupported new-backend/legacy combinations
  explicitly rather than silently changing their semantics or resource limits.
- A final count histogram is distinct from incremental rarefaction and may be
  supported by the new counter without reintroducing incremental machinery.
- Removing incremental counting entirely is not authorized by this plan.
  Any later removal needs an explicit decision, migration/deprecation plan,
  and treatment of `sharkmer_viewer`.

## Architecture and resource contracts

Agree these contracts in [#139](https://github.com/caseywdunn/sharkmer/issues/139) before independent implementations.

1. **Exact counts and membership.** Preserve canonical identity, all observed
   accepted k-mers, and abundance through ingestion/finalization. Specify the
   maximum exact count and observable overflow behavior. Do not default to
   saturating u8 counts or probabilistic membership. Packed u16 with exact
   overflow storage is a benchmark candidate, not a predetermined winner.
2. **Canonical keys versus oriented identities.** Count keys and directional
   primer/graph identities are separate concepts and should have explicit
   APIs/types. The current graph uses (k-1)-mer nodes and k-mer edges.
3. **One pass over original input.** Counting must support non-replayable
   sources. Optional local spooling/replay for evidence and external counting
   does not violate this contract. Never require a second remote download.
4. **Bounded end-to-end resources.** Budget count-table growth/resizing, queues,
   concurrent target graphs, evidence, and lookup caches together. An external
   counter followed by a full in-memory final table is not a bounded-memory
   solution.
5. **Hardware baseline.** Use a 16 GB laptop as the initial target, reserving
   operating-system headroom; validate larger jobs with disk-backed storage.
   A 32–64 GB machine is an additional benchmark class, not a prerequisite for
   the architecture. Record SSD/temp-space requirements and peak RSS.
6. **Replay and evidence.** Preserve the exact selected read subset, mate
   identity, read orientation, base positions, qualities, and invalid-base
   discontinuities. Stream against bounded batches of target graphs and keep
   compact evidence rather than all reads in RAM.
7. **Target-local graph work.** Benchmark unitigs and neighbor caching on
   target graphs first. Avoid constructing global background topology merely
   to answer a small panel.
8. **Explicit uncertainty.** No-data, contradictory evidence, unresolved
   phase/repeats, partial sequence, and exhausted search are different
   outcomes. Observed k-mers do not prove an observed full-length haplotype.
9. **Reproducible storage and execution.** Persist k, count precision,
   canonicalization, preprocessing, input/subset identity, and schema version.
   Fingerprint tested binaries and record actual parameters.

## v3.2 — Correctness and reliable evaluation

- [ ] [#129](https://github.com/caseywdunn/sharkmer/issues/129) — Make sPCR validation and benchmark provenance trustworthy.
  - [x] Implement per-product validation, executable/input provenance, current-run manifests, benchmark metrics, and offline regression checks; independent Astra review approved on 2026-09-08.
  - [x] Verify 145 Rust unit tests, 22 integration tests, 20 Python regressions, exact fixture counts/sequences, real BLAST ambiguity/split cases, and bounded ENA cache checks.
  - [ ] Register independent biological held-out inputs and truth before assembly-policy tuning. Existing datasets remain calibration/regression data; see [the dataset policy](../benchmarks/DATASETS.md). This acceptance gate keeps the issue open.
- [x] [#130](https://github.com/caseywdunn/sharkmer/issues/130) — Read every member of concatenated gzip FASTQ inputs. Terra implementation, Sol review; eight regressions cover ingestion, paired reads, limits, replay, and corrupt-member cache publication.
- [x] [#131](https://github.com/caseywdunn/sharkmer/issues/131) — Continue coverage thresholds until a valid amplicon is recovered. Sol implementation, Terra review; full per-threshold evaluation, authoritative length checks, explicit search limits, and an absolute primer-count floor preserve shared primer backbones without overriding the configured coverage ratio. Independent CLI mixture probe recovers the exact 180 bp target; existing high-copy oracle is unchanged.
- [x] [#132](https://github.com/caseywdunn/sharkmer/issues/132) — Prevent confident amplicon output with collapsed homopolymer lengths. Sol implementation, Terra review; path-local repeat markers survive pruning, iterative SCC detection is stack-safe, and ambiguous products are withheld with explicit reasons. Independent homopolymer/tandem CLI probes pass; short-repeat and existing high-copy controls remain unchanged. Intentional removal of shortened products is documented in README.
- [ ] [#133](https://github.com/caseywdunn/sharkmer/issues/133) — Correct read selection, strand handling, and gap continuity in threading.
- [ ] [#134](https://github.com/caseywdunn/sharkmer/issues/134) — Publish current-run amplicons and stats without stale FASTA results.
- [ ] [#135](https://github.com/caseywdunn/sharkmer/issues/135) — Bound primer ambiguity and mismatch expansion before allocation.
- [x] [#136](https://github.com/caseywdunn/sharkmer/issues/136) — Fix integer median rounding for even-sized k-mer count sets. Overflow-safe floor average, boundary and exhaustive small-pair regressions, caller/documentation audit; Terra review approved.
- [ ] [#137](https://github.com/caseywdunn/sharkmer/issues/137) — Make cache publication concurrent-safe and clearing ownership-aware.
- [ ] [#138](https://github.com/caseywdunn/sharkmer/issues/138) — Align sPCR documentation and diagnostics with current behavior.

Start with the benchmark/provenance issue and independent correctness fixes.
For homopolymers, v3.2 can conservatively withhold a falsely complete product;
full read-supported repeat reconstruction is a v4.1 task.

### Review evidence to preserve

The existing suite passed: 145 unit tests and 22 integration tests. Additional
temporary probes found gaps that must become durable regression fixtures:

| Case | Observed at review | Required behavior |
| --- | --- | --- |
| Two concatenated gzip members, one record each | One read counted, exit success | Both records counted |
| Shared primer ends; 400 bp at count 100 plus 180 bp at count 4; allowed 170–210 bp | No valid product; rare target succeeds alone | Recover the supported in-range target |
| 160 bp error-free target containing a long A run, k=19 | 137 bp product reported | Supported length or explicit uncertainty |
| Oriented primer seed TTTGA, k=5 | Its matching read rejected by canonical filter | Orientation-independent membership |
| AACGATTCCG versus reverse complement CGGAATCGTT | Different phasing evidence | Equivalent oriented graph evidence |
| AACGANACGAT on adjacent AACGA/ACGAT graph edges | False link across N | No continuity across the gap |
| Successful run followed by failed rerun at same output prefix | Old FASTA remains | Current-run manifest cannot report stale success |
| Counts [3,5] | Integer median 3 | Median 4 |

The fixture-only counting measurement was 100,000 reads, 10,698,337 distinct
k-mers, 1.4 s ingestion, 1.1 s redundant consolidation, 2.52 s total, and
562,688 KiB peak RSS on the review environment. Treat this as motivation,
not a portable performance target or a measured improvement.

### Release gate

- All reproduced correctness failures have regression coverage.
- Validate all products, including negative controls and whole-product
  alignment coverage; wrong-gene/short-fragment matches cannot validate a
  product.
- Existing high-copy recovery remains correct; intentional removal of
  previously incorrect products is documented rather than scored as a loss
  to undo.
- CLI examples, stats, failure reasons, cache ownership, and rerun semantics
  agree with the executable.

## v4.0 — Scalable exact counting and replay

- [ ] [#139](https://github.com/caseywdunn/sharkmer/issues/139) — Define exact count-store, replay-source, and assembly evidence interfaces.
- [ ] [#140](https://github.com/caseywdunn/sharkmer/issues/140) — Remove redundant table consolidation and per-read allocations from sPCR counting.
- [ ] [#141](https://github.com/caseywdunn/sharkmer/issues/141) — Parallelize ingestion with bounded queues and measured shard ownership.
- [ ] [#142](https://github.com/caseywdunn/sharkmer/issues/142) — Benchmark packed k-mer tables while retaining exact abundance.
- [ ] [#143](https://github.com/caseywdunn/sharkmer/issues/143) — Count non-replayable datasets exactly within a memory budget.
- [ ] [#144](https://github.com/caseywdunn/sharkmer/issues/144) — Provide exact bounded-memory random lookup over external k-mer counts.
- [ ] [#145](https://github.com/caseywdunn/sharkmer/issues/145) — Add bounded replay/spooling with consistent paired input semantics.

Sequence: contracts → simple counter → bounded parallelism → packed-table
evaluation. External counting and replay can progress independently after their
declared prerequisites. Final disk lookup depends on external counts.
The packed-table experiment is not a dependency of the first external backend.

Compare existing ahash/Fx configurations, owner-managed shards, bounded local
aggregation, and a specialized shared table on representative workloads.
Do not allocate a complete per-thread map without measuring duplication.
Jellyfish is a relevant exact-counting comparator; KMC is a relevant external
counter. Their existence does not prescribe Sharkmer's internal design.

Start external storage with the simplest measured design. Minimizer/signature
partitioning and super-k-mers are candidates for reducing temporary I/O;
partition skew, exact global counts, and final random access remain required.

### Release gate

- Exact count and product agreement across supported threads and backends.
- Standard `--chunks 0` has no incremental counting overhead; legacy
  histogram/combined-PCR compatibility remains tested separately.
- Record stage throughput, table bytes per distinct k-mer, resize/finalization
  peaks, RSS, disk use, and lookup latency. Freeze quantitative performance
  gates against the repaired baseline before choosing production designs.
- Demonstrate a one-shot dataset whose global count store exceeds the chosen
  RAM budget completing counting **and sPCR**, using disk-backed lookup.
- Verify limits, paired selection, malformed gzip, interruption, and disk-full
  behavior without losing previously accepted counts or exposing partial files.
- Standard sPCR accuracy does not regress while counting implementation changes.

## v4.1 — Read evidence and low-coverage nuclear recovery

- [ ] [#146](https://github.com/caseywdunn/sharkmer/issues/146) — Batch primer discovery without combinatorial variant enumeration.
- [ ] [#147](https://github.com/caseywdunn/sharkmer/issues/147) — Build length-bounded target graphs and compact nonbranching paths.
- [ ] [#148](https://github.com/caseywdunn/sharkmer/issues/148) — Use read-spanning branch histories to resolve repeats and constrain paths.
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

## Implementation handoff

- Begin with v3.2 evaluation/correctness issues. Do not launch a whole rewrite
  before regression fixtures and interface decisions exist.
- An issue is complete only when its acceptance criteria, focused tests,
  affected documentation, and release benchmark evidence are satisfied.
- Follow [CONTRIBUTING.md](../CONTRIBUTING.md) for implementation branches,
  commits, quality gates, and release process.
- Update this plan and the associated issue in the same implementation change.
  Do not infer completion from a similarly named historical closed issue.
- Existing [#121](https://github.com/caseywdunn/sharkmer/issues/121)–[#127](https://github.com/caseywdunn/sharkmer/issues/127) and [#76](https://github.com/caseywdunn/sharkmer/issues/76) include stale or overlapping status; audit them
  against code before closing or duplicating panel/cleanup work.
- This program does not remove incremental counting or implement v5 targeting.

## Design references

These are comparison/design references, not dependencies on their software.

- [Jellyfish: exact parallel counting with a specialized shared table](https://pmc.ncbi.nlm.nih.gov/articles/PMC3051319/).
- [KMC 3: counting and manipulating k-mer statistics](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399).
- [Linked de Bruijn graphs: retaining read-scale connectivity](https://pubmed.ncbi.nlm.nih.gov/29554215/).
- [Bifrost: compacted graph construction and indexing](https://pmc.ncbi.nlm.nih.gov/articles/PMC7499882/).
- [metaSPAdes: local coverage and metagenomic assembly challenges](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/).

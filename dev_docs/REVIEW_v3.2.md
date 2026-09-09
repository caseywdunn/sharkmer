# v3.2 pre-release review

Prepared 2026-09-09. **Development work only: awaiting user review, not release approval.**

## Scope and status

The v3.2 correctness implementation covers #129–138. The reviewed code point
is `03c0fc62ee0521f85a527b86f78252a60be1500b`, package `3.2.0-dev`.
Later review/evidence commits contain documentation and benchmark artifacts,
not a different counting or assembly implementation.

Incremental `--chunks` counting and histogram compatibility remain supported.
The fast exact counter, bounded-memory replay/storage, low-coverage nuclear
policy, and metagenomic diversity policy remain v4.0/v4.1/v4.2 work in
[PLAN.md](PLAN.md). This release does not claim to solve those harder inference
problems or to improve speed/RAM through a new counter.

Implementation and independent reviews:

| Issue | Implementer | Reviewer | Dev implementation commit |
| --- | --- | --- | --- |
| #129 validation/provenance | Sol | Astra | `6991d73` |
| #130 concatenated gzip | Terra | Sol | `1903d87` |
| #136 integer median | Root | Terra | `d47c3c1` |
| #131 threshold retry | Sol | Terra | `9284916` |
| #132 repeat uncertainty | Sol | Terra | `b728f14` |
| #135 primer budgets | Terra | Sol | `b41f6f8` |
| #133 oriented read evidence | Sol | Terra | `fb86ccd` |
| #134 current-run outputs | Sol; root validator integration | Terra | `2183d5e` |
| #137 cache ownership/coordination | Terra | Sol | `973edbb` |
| #138 documentation/diagnostics | Terra; root changelog/version | Sol | `03c0fc6` |

## Correctness evidence

- Release Rust tests pass with both default `ahashmap` and
  `--no-default-features --features fxhashmap`: **213 unit + 22 integration
  tests for each backend**. Formatting and both Clippy configurations pass.
- **39 Python regressions** pass. Both implementer and reviewer validated all
  nine built-in panels plus the reference example against the JSON Schema.
  Built-in loader and reference CLI validation pass. Six README sPCR commands
  pass dry-run parsing with the included fixture substituted for example local
  filenames and the original cnidaria panel substituted for its exported copy.
- Gzip regressions cover multiple members, corruption in a later member,
  record limits, paired reads, replay, and failed cache publication.
- An independent CLI mixture with shared primer ends, 400 bp at count 100 and
  180 bp at count 4, previously failed the 170–210 bp search window; the repaired
  implementation returns exactly the 180 bp target.
- At k=19, an A18 repeat control remains exactly 138 bp. A19, A40, and AC40
  targets that previously collapsed to 138 bp are now withheld with explicit
  repeat uncertainty rather than emitted at the wrong length. This is an
  intentional output reduction, not repeat-copy reconstruction.
- Threading tests cover interior reads, reverse complements, invalid-base
  continuity breaks, inverted motifs, strand ties, and repeated/overlapping
  evidence. The pinned high-copy fixture is unchanged with threading enabled.
- An independent real-SIGKILL probe confirms incomplete output cannot validate
  as current success, restart recovers the exact target, and an unrelated FASTA
  survives. Rust tests additionally cover partial publication and collisions.
- An independent CLI cache-clear probe waits for an existing lease through a
  directory alias, removes a verified legacy pair, and preserves modified,
  orphaned, symlink, and unrelated entries plus the cache directory/lock.
- Pathological ambiguous primers fail before ENA resolution under a 128 MiB
  address-space limit and five-second timeout. Integer median boundary and
  exhaustive small-pair regressions protect overflow-safe floor rounding.

These synthetic/fixture regressions protect specified behaviors. They are not
independent biological sensitivity or haplotype-accuracy measurements.

## Final calibration

All **13 historical samples across six panels** completed successfully at k=19,
two threads, and a requested cap of one million FASTQ records. All read/kmer
occurrence totals and all **71 product lengths/SHA-256 values** match the warm
repaired baseline. Selected input checksums/subsets/source plans and panel and
reference identities also match.

The baseline is the preserved executable historically labeled `b728f14`, used
after the threshold and repeat fixes, **not untouched 3.1.0**. The label is not
verified build provenance. It was selected explicitly;
its recorded workspace observations do not establish its build source. The
binary SHA is `70acbf1be0334f29f5d157e288dc2d509a5864cc17cf61a9e7e9ab6c7b8a31f1`.
The final executable was built by the harness immediately before fingerprinting
from clean `03c0fc6`, using release/default `ahashmap`; SHA
`387af7374e567828ed334ac488a032d79933f4f34a1da3823edc6c3ea7275b82`.

Final per-product reference-classifier outcomes are: 39 confirmed products,
2 confirmed gene/other taxon, 8 no reference, 6 insufficient alignment,
3 split/chimeric alignment, 5 ambiguous gene, 3 wrong gene, and 5 no significant
hit. **71 emitted products does not mean 71 validated biological products.**
The only classification changes are five human products corrected from
other-taxon to confirmed after registering the previously missing
`Homo sapiens` sample metadata. That is a metadata correction, not improved
sequence recovery. ENA scientific-name metadata was checked for SRR17535371;
the explicit taxon and provenance note are in `benchmarks/benchmark.yaml`.

| Panel | Accession | Actual records | Baseline wall s | Final wall s | Final peak RSS GiB |
| --- | --- | ---: | ---: | ---: | ---: |
| angiospermae | ERR14009273 | 1,000,000 | 31.95 | 32.99 | 2.132 |
| angiospermae | SRR25378184 | 1,000,000 | 62.21 | 63.07 | 4.257 |
| bacteria | ERR2596344 | 59,704 | 1.37 | 1.26 | 0.074 |
| bacteria | SRR19418213 | 525,982 | 14.69 | 14.84 | 1.070 |
| bacteria | SRR24806237 | 167,686 | 2.31 | 2.23 | 0.140 |
| cnidaria | SRR25099394 | 1,000,000 | 48.46 | 49.33 | 4.257 |
| cnidaria | SRR8617500 | 1,000,000 | 56.67 | 54.24 | 4.257 |
| cnidaria | SRR9278435 | 1,000,000 | 54.16 | 53.94 | 4.257 |
| human | SRR17535371 | 108,518 | 1.05 | 0.99 | 0.024 |
| insecta | SRR1057608 | 1,000,000 | 45.03 | 45.23 | 2.132 |
| insecta | SRR27962769 | 1,000,000 | 99.41 | 99.55 | 4.257 |
| insecta | SRR31887760 | 1,000,000 | 60.66 | 61.23 | 2.132 |
| teleostei | SRR22396603 | 1,000,000 | 46.52 | 48.05 | 4.257 |

These are single-run diagnostics, not evidence of a speedup. Both comparisons
use a warm application cache; OS page-cache state is unknown. The host has
12 CPUs and 62.6 GiB reported RAM, running Linux, Rust 1.98.1, and BLAST 2.17.0.
The final run had a 12 GiB address-space limit, which is not an RSS measurement
or a real 16 GB laptop validation. Four inputs exhaust below the requested
cap. Timings cover the Sharkmer invocation, not the subsequent BLAST analysis.
Stage/throughput, allocator, RSS, and final-disk metrics remain separate in the
raw results; peak temporary disk is not measured.

The pinned ERR571460 100k fixture passes both oracles: k=19 has exactly
13,197,385 kmer occurrences; k=31 retains exactly the established 18S (1,783 bp)
and 28S_2 (430 bp) products and has 11,996,833 occurrences. Threading preserves
those exact product hashes. The real-BLAST classifications for this fixture
remain other-taxon for 18S and no-significant-hit for 28S_2, distinct from the
exact-sequence regression oracle.

See the [machine-readable comparison](../benchmarks/benchmark_results/v3.2-review-20260909/comparison.json)
and [artifact index](../benchmarks/benchmark_results/v3.2-review-20260909/README.md)
for raw baseline/final results, input/binary provenance, and reproduction
commands. This is calibration/regression evidence only, not held-out validation.

## Behavior changes to review

- Standard sPCR retries lower thresholds after invalid candidates, but stops at
  the first threshold producing a valid product. It is not exhaustive
  rare-template recovery. Unresolved repeat-touched products are withheld.
- Optional read threading supplies local support, not whole-haplotype proof or
  molecule deduplication. Distinct FASTQ records can still represent duplicate
  molecules. Paired long-range links are not consumed by path selection.
- FASTA/stats outputs now require a complete, matching, checksummed current-run
  manifest. Reruns invalidate only verified prior outputs owned by that sample.
  Legacy/foreign/modified/symlink collisions are preserved and refused; use a
  fresh output directory when upgrading or retaining earlier results.
- Output publication is a manifest-mediated transaction, not an atomic rename
  of an entire multi-file directory. Readers must honor the manifest. Histogram
  and DOT outputs remain outside this transaction. Argument-validation failures
  occur before a transaction and leave previous results untouched.
- Cache leases serialize cooperating processes sharing a directory for the
  entire run, including replay. This prioritizes correctness over concurrent
  throughput. Older clients ignoring the lease are not coordinated.
- Normal failed cache replacement preserves the previous verified generation.
  Forced interruption can leave ambiguous/orphan entries; lookup refuses them
  and clear preserves them with warnings instead of guessing ownership. A fresh
  directory or `--no-cache` is the conservative recovery path.
- PCR requires k >= 2 and positive effective primer trim; count-only k=1 remains
  supported. Primer expansion has explicit checked work/allocation-count limits.
- `n_subreads_ingested` retains its existing legacy record-count semantics; it
  does not count N-split segments. `n_kmers` counts occurrences, not distinct
  keys. Allocator peak, process RSS, and final disk size are different metrics.

## Remaining release gates

1. **User review and release authorization.** Do not merge to `master`, tag,
   create a GitHub release, or publish the Bioconda recipe in this task.
2. **Independent held-out registration (#129).** No held-out biological inputs
   or truth have been supplied. Freeze accessions, checksums, subsets, callable
   targets, truth, and tolerances before further assembly-policy tuning. Existing
   historical data remain calibration/regression data; #129 stays open.
3. **Full release-depth and panel validation.** The bounded review uses a
   one-million-record cap, not the configured 2/4/8M depth sweeps. The historical
   matrix covers six panels, not biological validation of all nine. Run the
   broader release procedure in CONTRIBUTING.md and inspect every product and
   negative control before publishing; primer/schema loading alone is not
   biological validation.
4. **Dependency alert triage.** Four existing default-branch Dependabot alerts
   were observed during this review: rustls-webpki malformed CRL BIT STRING
   panic/DoS (high, #7), rand custom-logger unsoundness (low, #6), and
   rustls-webpki wildcard/URI name-constraint handling (low, #5/#4). These were
   not introduced or remediated by this work. Check applicability and patch or
   explicitly disposition them before publication; this was not a security audit.

   Read-only lockfile/advisory triage confirms the affected packages remain in
   dev. `ureq 2.12.1 -> rustls 0.23.37 -> rustls-webpki 0.103.10` is the runtime
   path; GitHub lists `0.103.13` as the smallest single patched target covering
   all three webpki alerts. Direct runtime `rand 0.8.5` is also affected; GitHub
   lists `0.8.6` as the first patched version for its release line. These are
   documented update targets, not tested dependency changes or an assessment of
   exploitability in Sharkmer. Primary advisories:
   [CRL panic](https://github.com/advisories/GHSA-82j2-j2ch-gfr8),
   [rand logger](https://github.com/advisories/GHSA-cq8v-f236-94qc),
   [wildcard constraints](https://github.com/advisories/GHSA-xgp8-3hg3-c2mh), and
   [URI constraints](https://github.com/advisories/GHSA-965h-392x-2mh5).

## Historical issues audited, not silently closed

Current-code audit comments were recorded on #121–126; all remain open:

| Issue | Remaining scope |
| --- | --- |
| #121 | `requires_sharkmer` absent; missing `panel_version` warns rather than enforcing the complete proposed schema contract. |
| #122 | Validator `--output-dir` controls Markdown only; machine results and run artifacts use separate locations. |
| #123 | Positive indices and complete indexing within mixed indexed/unindexed groups are not enforced; `gene_name` compatibility remains. |
| #124 | Optional annotation fields exist, but gene/copy metadata are incomplete in human, angiospermae, and bacteria panels. Do not infer missing biology. |
| #125 | Removing parsed-but-ignored `expected` is a future schema v3 migration, not this compatibility-preserving release. |
| #126 | Reference-sequence `notes` is absent in the Rust model and v2 schema. |

## Next implementation cycle

After this review, begin v4.0 with the count-store/replay/evidence contracts in
#139. Retain exact counts and one-pass ingestion; measure bounded end-to-end
resources before selecting packed counters or external storage. Keep legacy
incremental counting isolated rather than deleting it. Do not use a higher
FASTA product count as evidence of better nuclear/metagenomic inference.

# Independent minimum protocol review

2026-09-29. Design requirements only; this is not execution clearance for an
unreviewed new harness. No historical source, receipt, or measurement was edited.

## Comparator and frozen scope

Use pristine released-v3.1.0 commit
`5a664680c91ad59f32b8c2a847b8fef37f34a0ae` and the integrated inclusive-output
checkpoint `ed04bde` (resolve and pin its full commit). Build both from clean
exports with the same concrete Rust/Cargo 1.98.1 toolchain, locked dependencies,
release profile and default ahashmap features; pin source, lockfile, compiler,
configuration and resulting executable hashes. A new release-source build is
not the historical executable merely because its source commit matches.

Use the September 11 final ten-sample/five-panel/60-invocation scope, not the
broader September 9 thirteen-sample/deeper-input protocol. Preserve three
alternating release/candidate pairs per sample, k=19, two physical cores,
chunks=0, unpaired local uncompressed original ordered prefixes, threading off.
These are usually first-million R1 prefixes; Liriodendron contains 656,431
complete R1 records followed by 343,569 R2 records. The human input consumes
three complete ordered streams and exhausts at 108,518 records; all other
fixed prefixes contain 1,000,000. Pin the exact original prefix bytes, bases and record counts; a
restored download is eligible only when those identities match.

Freeze the actual five historical full-panel files used by both roles and the
original target-scope map. This preserves the comparison's assembly definition
and the distinction between metadata-selected high-copy targets and deferred
single-copy targets. Historical baseline expectations are 71 high-copy unique
literal target/sequence identities and 86 all-gene identities, including 15
deferred Yp2 identities. These are reproducibility checks, not correctness
labels or guaranteed new outcomes.

## Identity and recovery accounting

Read the complete candidate schema-2 transaction and mixed FASTAs, validate
class/order/count/length/provenance consistency, and independently verify that
the normal consumer selects exactly standard records. Validate the pristine
release's legacy output without imposing a schema or header field it never
produced. Use class-aware records instead of applying the old standard-only
analyzer to uncertain headers: uncertain products do not carry the ordinary
minimum/median k-mer support fields and none should be invented.

Report emitted records separately from unique `(target, literal sequence)`
identities. Record reverse-complement equivalence separately from literal
identity; do not silently canonicalize the original 71-to-64 metric. For every
baseline identity retain exact and RC evidence and summarize its disposition
as standard-exact, standard-RC-only, uncertain-exact, uncertain-RC-only, or
unobserved. Report candidate gains by class and separately for high-copy and
deferred targets. Product indices and FASTA order are not biological identities.
Do not multiply unique product totals by three timing repetitions.

Require within-role repeat stability of sequences, classifications, inclusive
counts and relevant search/retention assessments; require exact aggregate
read/base/ingested-k-mer count parity across roles. Report gene recovery in two
explicit forms: standard recovery and any emitted sequence availability.
Uncertain availability does not restore standard admission or validate repeat
length. An empty bounded collection does not establish biological absence.

## Reference evidence

Full BLAST reclassification of 60 duplicate repetitions is not necessary for
the authorized sequence-availability/resource comparison. Identity/resource
claims can be completed independently of new biological classification.
Existing per-product evidence may be reused only with matching sequence,
target, taxon, reference/catalog, classifier/tool and parameter provenance.
Otherwise mark it unassessed; missing evaluation is not a no-hit result.

If new reference annotation is performed, classify unique sample/target/
sequence combinations once, consistently across both roles, outside timing.
Pin the current verified reference catalog and classifier as a separate
posthoc layer; do not use obsolete historical bootstrap references as
independent truth. Historical assembly panels and current annotation catalogs
may deliberately differ, with the distinction recorded. Existing verified
reference provenance need not be reacquired or re-audited without a change.

## Measurement and integrity

Freeze the exact 60-command schedule before observations. Run serially on the
same allocated host and two pinned, distinct physical cores; record actual
CPU/socket/core topology and observed child affinity, not only a requested
thread count. Prewarm the same consumed prefix before each invocation outside
timing. Input verification, compilation, classification and result parsing
also stay outside the measured command. Time the full ordinary CLI, including
its inclusive publication overhead; this is not an isolated algorithm study.

Use GNU time's whole-command elapsed time and peak process RSS for both roles,
retain raw logs and process exit status, and separately retain observer time.
Apply the same 1,800-second timeout and 40-GiB address-space ceiling with
process-group cleanup. The latter is not a 40-GiB RSS cap or laptop-memory
validation. Preserve all attempts; do not silently replace a failed/slower
replicate. Failures, timeouts, count mismatches or incomplete provenance
invalidate the affected comparison rather than disappearing from summaries.

Report each sample's three measurements, median/range, paired ratios and
absolute deltas. Report every increase; >5% per-cell median runtime and >2%
median RSS remain investigation priorities, not acceptable-regression
allowances or automatic pass criteria. Tiny human timings require their
absolute differences alongside percentages. These repeats support descriptive
same-host results, not statistical equivalence or universal performance claims.

Pin and bookend sources, binaries, panels, inputs, tools, runtime, every output
and receipt. Use new no-clobber artifact paths. Archive the exact checked
source payloads, all outputs, raw timing logs and attempts; preserve external
input/compiler identities without nesting old archives. Existing full tests
on the exact integrated source need not be rerun solely because of a fast-forward
integration; new harness controls and concrete build/source binding are needed.

## Separate Drosophila diagnostic

The proposed scratch-only observer should hash every completed repeat-withheld
path before any retention early return, using the exact emitted-oriented
spelling plus its reverse complement, threshold and completion ordinal. A
bounded ledger must flag overflow, and its observed count per threshold must
equal the existing withheld-completion counter before absence among observed
paths can be inferred. Bind retained diagnostic sequences back to the ledger.

Require the 475-bp 12S and 590-bp 16S positive identities, synthetic paths
excluded by each retention cap, and unchanged ordinary outputs, counts,
threshold/search counters and collection budgets. Keep observer timing out of
the release performance matrix. Even a complete ledger with no 487/499 match
means those identities were not observed in this bounded enumeration; it does
not establish graph or biological absence.

## Local sources consulted

- `dev_docs/BENCHMARK_high_copy_followups.md`
- `benchmarks/benchmark_results/v3.2-high-copy-final-20260911/provenance.json`
- `benchmarks/benchmark_results/v3.2-high-copy-final-20260911/README.md`
- `benchmarks/benchmark_results/v3.2-high-copy-final-20260911/tools/analyze.py`
- `benchmarks/benchmark_results/v3.2-high-copy-final-20260911/tools/reclassify_audit.py`
- `benchmarks/benchmark_results/v3.2-release-comparison-20260909/protocol.json`
- `benchmarks/benchmark_results/v3.2-release-comparison-20260909/tools/driver.py`
- `dev_docs/REFERENCE_PROVENANCE.md`
- `dev_docs/PUBLIC_REFERENCE_REPLACEMENTS.md`

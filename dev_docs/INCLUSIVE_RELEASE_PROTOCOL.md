# Inclusive v3.2 release comparison

Authorized 2026-09-29: integrate the completed inclusive output into `dev`,
repeat the original high-copy comparison against released v3.1.0, and examine
the two remaining Drosophila 12S identities. `dev` was fast-forwarded to
`ed04bde`; `baseline/pre-inclusive-dev-20260929` preserves `956e771`.
This authorizes development integration and assessment, not publication.

## Frozen comparison scope

Repeat the September 11 final high-copy matrix: ten samples across cnidaria,
human, teleostei, angiospermae and insecta, three alternating baseline/candidate
pairs per sample (60 invocations). Preserve the original sample and pair order,
full panel bytes and high-copy/deferred target map from that final evidence.
The 15 historical Gryllus Yp2 identities remain outside the high-copy gate.

Use pristine git-exported released v3.1.0 (`5a664680c91ad59f32b8c2a847b8fef37f34a0ae`)
and inclusive `ed04bde`, compiled with the same project-local Rust/Cargo 1.98.1,
locked offline release builds and default ahashmap features. Record exact full
commit IDs, compiler/configuration identities, source payloads and binary hashes.
The candidate's already completed checks bind this exact implementation; the
new comparison harness requires its own controls and independent review.

Use the exact previously consumed uncompressed prefixes: one million records
for nine samples; human has 108,518 records from the archived ordered
concatenation of three EOF-verified files. Restore absent prefixes from their
recorded public sources and require exact archived byte hashes, lengths, record
counts and base counts. A source URL or accession alone is insufficient.
Use identical historical assembly panels for both versions. Their old reference
annotations are historical metadata, not independent biological truth.

Run k=19, two distinct physical cores, chunks=0, max-reads=1,000,000, no read
threading, one invocation at a time. Prewarm and verify inputs outside timing.
Measure the entire CLI with GNU time wall seconds and process peak RSS;
distinguish allocator peaks. Enforce 40 GiB address-space and 1,800-second
per-invocation limits and record actual child affinity/limits. Preserve failures
and partial attempts without replacing them or silently selecting repetitions.

## Outcomes

Check exact aggregate read/base/k-mer occurrence parity, repeated sequence
stability, complete per-target stats/file agreement, candidate schema-2 manifests
and status annotations. Validate the complete candidate output before checking
the standard-only consumer selection. Keep legacy v3.1 header parsing separate,
including its half-integer medians. No invented support fields for uncertain
records are permitted.

Report emitted records and unique target/sequence identities separately,
partitioned into standard and uncertain classes and high-copy/deferred scope.
For every released identity, record standard exact, standard reverse-complement
only, uncertain exact, uncertain reverse-complement only, or unobserved, with
the underlying literal/RC evidence. Report candidate gains separately. Five
previously recovered historical identities are expectations to check, not
biological correctness claims; keep the 487/499 bp Drosophila losses explicit.

Report all per-sample replicate timings/RSS, medians and paired changes. Flag
every loss or resource increase; >5% median runtime or >2% median RSS increases
are priority review triggers, not permissible-loss allowances. Three same-host
replicates are descriptive, not statistical equivalence or all-platform claims.

This study measures identity availability and resources. No new BLAST
classification is required; any absent assessment remains unassessed. Do not
promote obsolete bootstrap-reference labels into independent truth. Wider
panels, held-out biological validation and dependency review remain separate
release work.

## Focused Drosophila diagnostic

Use a separately built scratch copy of `ed04bde` with a bounded pre-retention
hash observer. Keep the production search, admission and retention limits
unchanged. Observe actual oriented complete withheld sequences before a
retention return, with threshold, ordinal, literal/RC hashes and disposition.
Apply the observer uniformly to the full insect panel; do not embed historical
target hashes in the executable. Compare those identities only after execution.

Require controls for retained and dropped paths, orientation, duplicates,
empty ledgers, overflow and safe publication. Bind each ledger to the actual
evaluated threshold and require complete Drosophila 12S accounting before
inferring absence among enumerated paths. Check ordinary products, inclusive
outputs, aggregate counts and all existing search/retention diagnostics against
the uninstrumented candidate. Any incomplete ledger limits interpretation;
absence from a complete ledger means not enumerated under these search bounds,
not absence from the graph or the biological sample. The diagnostic cannot
silently change the production candidate or the timed comparison.

# Public-reference audit and variation-aware reclassification

Date: 2026-09-13. This is an evaluation/provenance update under #129, not a
release or a new throughput benchmark. See
[the findings and contract](../../../dev_docs/REFERENCE_PROVENANCE.md).

## Findings

- All 159 old built-in/example entries are accounted for: 100 exact public
  source matches and 59 source mismatches; two of the 100 are additionally
  quarantined for a known 16S/mutS annotation conflict.
- Final panels have 99 independent public regions: 98 retained old identities
  rebuilt from source plus one new Hydra documentation example. Primers and
  validation samples are unchanged. Historical sequences remain regression-only.
- The 109 unique public records contain 2,123,261 bases. Raw transfer has 111
  FASTA records / 2,143,980 bases because two identical records occur twice.
- Reclassification verifies 60 previous result/FASTA receipts, then assesses
  150 representative products after exact three-replicate fingerprint checks.
  Both versions retain the same 30 gene-supported product identities. The
  seven known high-copy losses remain unresolved by these limited references;
  absent support is not proof of biological incorrectness.
- Four fresh offline assembly smoke runs (before/after, k19/k31) preserve all
  five ingestion-count fields and every product sequence hash. No counting,
  traversal, repeat-gate, or default-k policy changes are made.

## Contents

`evidence.tar.gz` contains:

- `source-archive/`: raw public NCBI metadata and FASTA acquisition.
- `catalog_nonunique_final.json`: source-only audit before annotation decisions;
  its 100/59 counts are not the final number of active references.
- `legacy_references.json`: all 159 original entries, explicitly regression-only,
  with final 98 retained / 59 source / two annotation dispositions.
- `migration_proposal_nonunique_final.json` and `annotation_decisions.json`:
  exact regions, versions, alternative repeat placements, source/sample taxon
  differences, and the additional known annotation exclusions.
- `frozen_before/panels/`: original panels from commit
  `6fe7658151fdda0ce758b6db328619f3bac1de45`.
- `source_tree/`: final evaluated panels/catalog, validator source, relevant
  schema/tests, and the Rust metadata-only change, fingerprinted per file.
- `reclassification/`: frozen script, 60 source-result receipts, pinned BLAST
  executable and reference database metadata receipts, old/new per-product
  classifications, supported-identity comparison, and migration invariance.
- `assembly-smoke/`, `assembly_smoke.py`, and test logs: fresh offline
  before/after run receipts, FASTAs/manifests, count/sequence comparisons, and
  223 Rust unit +22 integration +77 Python test results.

The final catalog SHA-256 is
`f6abbfa3ffab56cc0be24a0a9b560b33033c499de7c949ae6fff988f8f91f6c1`.
`annotation_conflicts` is reviewed local metadata, not a raw NCBI field.
Checksums bind the trusted snapshot; they are not live public-record authentication.

## Integrity and reuse

From this directory:

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` records every member's size and SHA-256. Extract into
a fresh scratch directory. The source-only audit can run offline against
`frozen_before/panels` and `source-archive`; apply `annotation_decisions.json`
afterwards to reproduce the final eligibility policy. The normal repository
offline check is `python scripts/audit_panel_references.py --verify-existing`.

Reclassification reuses the immutable raw outputs in the
[September 11 archive](../v3.2-high-copy-final-20260911/README.md); those outputs
are not regenerated or overwritten here. The original receipt paths assume
restoring that archive under `/tmp/sharkmer-high-copy-20260911`. The archived
`reclassify_saved_outputs.py` accepts `--repository-root`, `--results`,
`--catalog`, and `--output`; use the frozen `source_tree` for reproducible code
and panel identity, and put the fingerprinted BLAST+ 2.17.0 tools on `PATH`.
Database indexes are rebuilt and verified during each run; timestamp-bearing
binary-index hashes are not required to match across builds. Original binaries
and large reads are not embedded; their hashes identify them.

Temporary absolute paths and original timestamps document the measured
environment. Do not edit historical receipts to make a new run appear old.
Neither these public regions nor the reused samples are newly held-out truth.
Read support, variant phasing, gene-annotation completion, and release approval
remain separate work.

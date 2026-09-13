# Annotation-reviewed public-reference replacements

Date: 2026-09-13. Follow-up to #129 and #153; **not a release or a new
throughput benchmark**. See the [findings and accession-level source
map](../../../dev_docs/PUBLIC_REFERENCE_REPLACEMENTS.md).

## Findings

- The 61 quarantined entries are not 61 proven Sharkmer products: 59 lacked
  exact cited-source extraction and two had a reviewed 16S/mutS annotation
  conflict. Historical strings stay regression-only.
- 58 have public alternatives: 57 newly covered legacy IDs share 56 new
  reference identities, and the documentation example already uses public
  Hydra 16S. The two H3 IDs share one independently annotated genomic locus.
- Two plant trnV-atpE references and Agalma ITS remain unresolved. Related
  taxa, partial loci, predicted transcripts, and actual source organisms are
  explicit. Primer/target review for the plant conflicts is tracked in #156.
- An additional Heliconius ND1 reference and expansion of Gryllus CO1 bring
  active entries from 99 to 156, all exact extractions from a 121-record pinned
  catalog. Six panel patch versions change; primers and sample settings do not.
- Reclassification verifies 60 prior invocations / 516 raw-file receipts and
  20 stable triplicates, then assesses 150 representative products. Both
  versions have the same 51 gene-supported identities (38 expected taxon,
  13 other), up from the same 30 in the earlier limited-reference evaluation.
  This is better evaluation coverage, not new assembly recovery.
- Seven historical high-copy losses remain unresolved, not proven wrong.
  Gryllus 373 bp CO1 and Heliconius 263 bp ND1 have insufficient full-query
  coverage; five other losses have no significant hit. Retained Gryllus
  352 bp CO1 now has 100% identity / 99.72% query coverage against public CO1.
- 97 Python tests including real BLAST, 223 Rust unit tests, 22 integration
  tests, all 10 panel schemas, source audit, formatting, and Clippy pass.

## Contents

`evidence.tar.gz` includes:

- `cnidaria_other/`, `insecta/`, `plants_nematode/`: raw public GenBank
  records/feature tables, bounded searches, requests/receipts, extraction
  proposals, and scratch validators. The root-level request and receipt files
  accompany cnidarian and other source acquisitions.
- `integration_receipt.json`: every legacy/additional decision, source scope,
  coordinates/strand, hashes, shared identities, exclusions, and proposal
  bindings. The documentation example explicitly keeps prior Hydra evidence
  rather than the proposed Xenia substitute.
- `insecta/nadh_primer_mapping.json`: corrected gene-selected assay mapping;
  the initial array-index-selected CO1_2 analysis is superseded, not endorsed.
- `frozen_before/`, `proposed_after/`, `integration.patch`, and
  `prepare_integration.py`: reference-only changes from source commit
  `66704f95d8a36ce9fcdee1cf1371eae67d8b3e6a`. Unchanged panels need not occur in
  these delta directories; `source_tree/panels/` is the full evaluated set.
- `source_tree/`: evaluated panel/catalog, validator, tests, schemas,
  documentation, and CI snapshots. No Rust assembly implementation changed.
- `reclassification/`: frozen wrapper, command, summary, and receipt-bound
  per-product results. The wrapper's current target-mapping argument is the
  only compatibility update to the earlier saved-output assessment script.
- `active_reference_audit.json`, `panel_invariants.json`, and test logs:
  exact-source verification, schema/unchanged-parameter checks, and test runs.

The catalog SHA-256 is
`bc277321831258c00b6794a51b0289445c837f99b2399e54aa43fe5e60b0f49a`.
It preserves all old source sequences/taxa and reviewed annotation exclusions.
Its original `entries`, `summary`, `panel_receipts`, and `source_archive` remain
initial-audit metadata, not the current 156-entry inventory; use the new active
audit and integration decisions for current counts.

## Integrity and reuse

From this directory:

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

`ARCHIVE_CONTENTS.json` binds every member by size and SHA-256. Source-origin
hashes are consistency checks, not live NCBI authentication or absolute
biological truth. Extract into a fresh scratch directory; do not overwrite
the older evidence archive or original source results.

The offline repository audit is
`python scripts/audit_panel_references.py --verify-existing` (PyYAML required).
For reclassification, restore the immutable raw outputs from the
[September 11 evidence](../v3.2-high-copy-final-20260911/README.md) at their
receipt paths under `/tmp/sharkmer-high-copy-20260911`, put the fingerprinted
BLAST+ 2.17.0 tools on `PATH`, and run the archived wrapper with
`--repository-root` pointing to extracted `source_tree`, `--catalog` to its
`panels/reference_sources.json.gz`, `--results` to those original results,
and `--output` to a new JSON file. Database binary hashes can differ across
rebuilds; query-time hashes are checked, and the original initial executable
receipts are retained. No after-query executable attestation is claimed.

Acquisition/preparation scripts retain original scratch paths and rely on the
prior regression-only legacy manifest. Requests are preserved for review;
do not silently refresh source records to recreate a historical snapshot.
The old [reference-audit archive](../v3.2-reference-provenance-20260913/README.md)
is unchanged. These public references and reused samples are not newly held-out
truth. Population variation, sample haplotypes, read support, and release
approval remain separate questions.

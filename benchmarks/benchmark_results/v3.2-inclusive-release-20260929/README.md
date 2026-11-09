# Inclusive release comparison, 2026-09-30

The original ten-sample/five-panel matrix was repeated against released
v3.1.0 (`5a664680c91ad59f32b8c2a847b8fef37f34a0ae`) and inclusive output
(`ed04bdec5178d055b47cc8085135ee5ae2600998`): all **60 invocations pass**.
The candidate is integrated locally into `dev`; the previous baseline is
preserved at `baseline/pre-inclusive-dev-20260929` (`956e771`).

Default output contains **126 records: 64 standard and 62 uncertain**.
It retains 69/71 released high-copy identities (64 standard, five uncertain).
Two high-copy identities and 12/15 deferred Gryllus Yp2 identities remain
absent. All 54 new identities are uncertain high-copy hypotheses. This is
identity availability, not validated repeat structure or release approval.

All ten median wall times increase; Human, Liriodendron, Rhopilema and
Heliconius exceed the 5% investigation threshold. No median peak-RSS increase
exceeds 2%. The sum of median wall times is 434.01 → 450.07 seconds (+3.70%).
The comparison includes all intervening release changes and cannot assign
this overhead specifically to inclusive emission.

The separate scratch observer finds exact Drosophila 12S 499/487 bp paths
at completion ordinals 13/19, discarded by the fixed run-marker share.
All 26 completed 12S paths are accounted for. Normal production output is
unchanged; this diagnosis does not re-emit or biologically validate them.

## Reports and evidence

- [Frozen protocol](../../../dev_docs/INCLUSIVE_RELEASE_PROTOCOL.md),
  [full results](../../../dev_docs/INCLUSIVE_RELEASE_RESULTS.md), and
  [Drosophila retention diagnosis](../../../dev_docs/DROS12S_RETENTION_RESULTS.md).
- [Input inventory](inputs/inventory.json),
  [original historical sequence identities](inputs/historical_identities.json),
  and [input provenance notes](inputs/INVENTORY.md).
- [Posthoc report](results/report-v2/report.json), with tables for
  [counts](results/report-v2/sample_counts.tsv),
  [all released identities](results/report-v2/released_identity_dispositions.tsv),
  [new identities](results/report-v2/literal_gains.tsv),
  [the seven prior losses](results/report-v2/historical_seven.tsv),
  [resources](results/report-v2/resources.tsv) and
  [paired resources](results/report-v2/paired_resources.tsv).
- Independent [raw outcome review](review/raw_outcome_review.json),
  [Drosophila outcome review](review/dros12s_outcome_review.json),
  and final report/document review in `review/`.
- `evidence.tar.gz`, `archive_receipt.json`, and `SHA256SUMS` preserve
  the expanded source snapshots, all measured outputs and logs, exact
  source/receipt identities, controls, failed-attempt lineage and reviews.
  The archive embeds its member manifest and checksums. External raw inputs,
  toolchains, binaries, build caches and nested evidence archives are omitted;
  their identities remain in receipts. Small checked-source FASTQ fixtures
  are retained as source.

Local execution roots are `tmp/inclusive-release-20260929/build-v1/`,
`comparison-v1/`, `restore-v1/`, `restore-v2/` and `dros12s/`.
The production comparison is Slurm job `27936421`; the scratch replay is
`27930708`. The archive preserves portable payloads after those scratch
locations are removed.

## Attempt lineage and checks

Both clean builds pass with Rust/Cargo 1.98.1. The exact candidate's
[preceding compiled/source validation](../v3.2-inclusive-output-20260929/README.md)
remains applicable: both hash backends pass 247 unit and 29 integration
tests, formatting and Clippy; Python has 197 tests with one existing skip;
all 156 references verify. The new main harness has 15 author/independent
controls. The scratch observer passes 250 unit and 29 integration tests,
formatting, Clippy and ten Python controls.

The initial restoration attempt verified eight inputs and rejected two
because the recorded local gzip-cache SHA was compared with upstream gzip
bytes. Versioned retry code checks recorded upstream MD5/size and EOF while
still requiring the exact original raw prefix SHA, size, records and bases.
All ten inputs then pass; both attempts and their code remain preserved.
The corrected restoration and merge have 24 controls.

Posthoc report attempt `report-v1` preserves a numeric-check failure caused
by Python 3.9/3.12 summation differing by 5.68e-14 seconds. Separately named
`report_comparison_v2.py` uses `math.fsum` and a 1e-9-second absolute tolerance,
with zero relative tolerance. It rejects material changes; the measured runs
are unchanged. The original reporter and failed report remain preserved.

No production retention caps, search behavior or admission rules were
changed for this comparison. Remaining decisions concern bounded retention,
runtime regressions, biological validation and release review.

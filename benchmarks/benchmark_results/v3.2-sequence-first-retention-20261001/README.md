# Sequence-first uncertain retention assessment

This bundle evaluates a bounded retention change on
`experiment/sequence-first-retention`, separate from `dev`. The code reserves
the old detailed-path allowance and keeps additional completed repeat-uncertain
sequences under a separate bounded allowance. It marks omitted detail in the FASTA header and schema-2
uncertainty assessment, and leaves standard admission unchanged.

All 60 corrected old-inclusive/candidate invocations pass. The corrected
build retains all 126 old records and adds 200 uncertain identities, for
64 standard and 262 uncertain records (326 total). The sum of ten median wall
times rises 1.8%; no sample crosses the 5% wall or 2% RSS review thresholds.
The corrected Drosophila replay passes: 26 12S paths now appear in the normal
target FASTA versus three before; exact historical 499 bp and 487 bp
identities are present. The 23 newly emitted paths are labelled uncertain
with summary-only marker detail. Nine standard sequences are unchanged.
This is sequence availability, not validated repeat structure.

The [protocol](../../../dev_docs/SEQUENCE_RETENTION_PROTOCOL.md) and
[results](../../../dev_docs/SEQUENCE_RETENTION_RESULTS.md) describe the
change, bounds, tests, and interpretation. The independent
[full raw-output review](full-review-v3.json),
[released-identity review](released-identity-review-v3.json),
[compact report](report-v3/report.json), and
[count](report-v3/sample_counts.tsv), [resource](report-v3/resources.tsv), and
[gained-identity](report-v3/uncertain_gains.tsv) tables support the result.
The focused [Drosophila review](dros-review-v3.json) and
[candidate 12S FASTA](evidence/dros12s_candidate_v3.fasta)
make the focused result inspectable. The
[old inclusive 12S FASTA](evidence/dros12s_old_inclusive.fasta) and
[structured Drosophila summary](evidence/dros-summary-v3.json) are controls.
The corrected [plant review](plant-review-v3.json) confirms preservation of
all 27 old uncertain identities while adding 64. The passing compiled-check
logs and source/binary SHA receipts are in `evidence/checks-v3/`. The first
replay harness attempt failed before invocation
on an import-path error; its Slurm log is preserved under
`tmp/sequence-retention-20261001/slurm/dros-v1-28034088.log`. The corrected
Slurm replay is `28034278`.

The first full-matrix attempt (`28034535`) was stopped at the Liriodendron
pair because the shared path cap displaced one old uncertain identity. Its
raw outputs and receipts are preserved in
`tmp/sequence-retention-20261001/full-v1/`. The corrected candidate reserves
old detailed paths separately. Its checks and focused plant/Drosophila replays
pass. The corrected full comparison is complete as Slurm `28037442` under
`tmp/sequence-retention-20261001/full-v3/`.

The separate diagnostic-only opt-out control (`28045941`) passes with nine
identical old/new standard records and 50 unchanged schema-1 hidden threshold
payloads. Its raw receipts are under
`tmp/sequence-retention-20261001/optout-v1/`.

`evidence.tar.gz`, `archive_receipt.json`, and `SHA256SUMS` bind both attempts,
the corrected raw outputs, scripts, checks, and review files. The archive
omits the large FASTQ inputs; their exact SHA-256, sizes, record and base
counts remain in the run receipts and prior input inventory. All archive
members are rehashed after packaging.

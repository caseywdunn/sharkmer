# Independent preregistration review

Reviewed 2026-09-27: `dev_docs/UNCERTAIN_PRODUCTS.md`, 7,694 bytes, SHA-256
`762c4cd0e12dc8070e12dd98a77cec17afea730d4fbeeb8bf5b634f8c65a6092`.

The protocol matches the inspected collector, uncertainty finalizer,
publication, and validator. No correction is required at this checkpoint.

- The three cap dimensions match source constants. The frozen insect panel
  has 21 targets: the first four receive 13 run-share path records and the
  other 17 receive 12. Gene and threshold caps remain additional limits.
- Exact oriented deduplication, retained reverse-complement alternatives,
  final-standard exact/RC exclusions, zero-candidate assessments, and separate
  status/count semantics match the implementation.
- Gene-wide search and retention summaries are distinguished from emitted
  sequence length and from biological repeat length. Metadata validation is
  correctly described as source-local checking, without graph replay or
  biological support claims.
- The nine planned invocations hold inputs, full panel, search parameters,
  binary identity for off/on, and resource envelope fixed. Ordered standard
  FASTA/header and existing count/threshold parity are acceptance conditions.
  Historical oriented matches, RC matches, and bounded non-observation are
  distinct; approximate similarity is not an exact restoration.
- The assessment is explicitly descriptive and uses calibration samples.
  It does not establish biological correctness, full repeat enumeration,
  independent read support, or release readiness.

The code reviewer raised a potential empty-marker case at k=1. The
implementation author identified existing `effective_trim` validation in
`src/pcr/primers.rs` and its k=0/k=1 tests: PCR rejects zero effective trim
before graph construction. Counting-only k=1 does not invoke PCR. No new
emission restriction was necessary.

This is protocol alignment review, not compiled-code acceptance or biological
execution clearance. Final harness, executable, input, and source freeze
review remains outstanding.

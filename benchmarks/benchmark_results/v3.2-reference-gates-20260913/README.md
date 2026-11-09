# Accessioned-reference repeat-gate assessment

See [the report](../../../dev_docs/BENCHMARK_reference_gates.md) and open
[release gate #153](https://github.com/caseywdunn/sharkmer/issues/153).
These are correctness experiments, not comparative performance measurements.
No production code, defaults, release tags, or package publication change.

## Evidence

`evidence.tar.gz` preserves raw FASTAs, stats, manifests and logs; frozen
synthetic fixtures and protocol; reference inventories and exact-match results;
independent audits; source/binary build receipts and compiler logs; helper
scripts; and the panel files used as reference oracles. `ARCHIVE_CONTENTS.json`
records each member's size and SHA-256. `SHA256SUMS` covers the top-level files.

The real-data contrast uses pre-#132 commit
`928491610822ffdb236d0326b32c8dd19fecf66a`, immediate post-#132
`b728f14219945567b3d78f3602590864a256c809`, and current production code
`0c38d6a051082cc2a1eb961b4206837f82948fda` on identical Gryllus input at k=19.
The raw target-specific evidence proves that #132 removes the exact 978 bp
panel-embedded ITS_2 sequence. All five ingestion/count fields agree. A public
record check discovers that AK281180.1 is only 256 bp: the 978 bp panel entry is
a bootstrap amplicon, not an independently deposited full-length reference.
Therefore this contrast alone does not establish a biological false negative.

The native-endpoint hydrozoa diagnostic changes primer sequences and length
bounds while leaving accessioned sequence content unchanged. It must not be
interpreted as recovery with the production hydrozoa panel. Ambiguous and
partial reference eligibility, exact full-reference products, internal matches,
and noncontiguous products are kept distinct.

All 150 synthetic runs complete with version-matched counts. All 23 hydrozoa
templates reconstruct exactly at k19/k31 in every version; 22 match independently
deposited regions, whereas EU293971.1 is a one-base-deleted panel derivative.
The Gryllus bootstrap template also reconstructs from clean reads. A18 stays
correct; the known A40 collapsed product is removed by the gate. There is no
demonstrated correct-template suppression in these controlled experiments.

## Verification and replay

Verify top-level checksums with `sha256sum -c SHA256SUMS`. Extract the archive
into a new directory and check every member against `ARCHIVE_CONTENTS.json`
before analysis. Raw emitted products are retained; reference classification
is posthoc and never supplies sequence to Sharkmer's assembly algorithm.

Original absolute paths are retained in provenance rather than rewritten to
pretend a different execution location. The scripts describe the original
`/tmp/sharkmer-reference-gates-20260913` workspace and its dependencies;
they are evidence/replay recipes, not a relocatable installed test suite.
The real input FASTQ and compiled binaries are intentionally not duplicated
in Git. Restore or regenerate the recorded checksum-matching input and build
the pinned commits with the archived locked/offline clean-export builder.
Source archives are reproducible from those Git commits; their hashes and
clean source-tree/build receipts are preserved, not duplicate source tarballs
or compiled target directories.

Canonical synthetic evidence is `protocol-v2.json`, `inventory-v3.json`,
`fixtures-v2/fixtures.json`, and `results-v2/summary.json`. The real-data
contrast and independent receipt audit are under `real-counterfactual/`.
`reference_gate_evidence_v3.json` is the initial offline embedded-sequence
comparison, not proof that every labeled sequence came from its public accession.
The authoritative public checks are `public-records/external_corroboration_final.json`,
`public-records/hydrozoa_public_provenance_final.json`, and
`public-records/gry_co1_public_alignment_final.json`. Superseded draft inventories
and public comparison summaries are not included. `astra_review.json` is a
recorded independent review attestation, not a substitute for raw evidence.

The original staged input and old-product classification evidence are in the
[earlier high-copy archive](../v3.2-high-copy-final-20260911). The previous
[k sweep](../v3.2-k-sweep-20260911) remains unchanged. These panel references
and biological inputs are calibration data, not independent held-out truth.

# Longer-k high-copy comparison, 2026-09-11

This is calibration evidence for open issue #153, not release approval.
The [findings document](../../../dev_docs/BENCHMARK_high_copy_followups.md#longer-k-results)
contains the interpreted sequence changes, resource tables, and diagnostics.

## Result

Seventy-two timed invocations compare released v3.1.0 with current production
code at k=19/23/27/31 on Drosophila SRR31887760, Heliconius SRR1057608, and
Gryllus SRR27962769: identical first 1M records and unchanged whole insect
panel, three alternating paired repeats. All run, count-parity, and within-version
sequence/classification stability checks pass.

| k | Released/dev high-copy products | Exact targeted losses restored by dev | Dev median-time sum (s) |
| --- | --- | ---: | ---: |
| 19 | 33 / 26 | 0 | 198.87 |
| 23 | 34 / 27 | 0 | 188.82 |
| 27 | 28 / 28 | 0 | 179.21 |
| 31 | 37 / 30 | 0 | 173.96 |

No longer k restores any of the seven exact k19 release-baseline losses.
Released v3.1.0 also loses all seven at longer k: k27 parity is not rescue.
No setting qualifies under the frozen selection rule, so the conditional
42 wider confirmation invocations were not run. The confirmation template
is preserved as planning input, not an executed experiment.

Larger k does add some reference-confirmed products (Drosophila CO2_1 and
Gryllus 28S), retains all 12 k19 dev `confirmed_product` sequences, and reduces
the three-sample sum of dev median runtimes by up to 12.53%. It also loses
two k19 dev sequences, including other-taxon-supported Drosophila ND4; do not
claim all reference-supported products survive. RAM is effectively unchanged.
None of this changes the default k, assembly policy, or release status.

## Identities and protocol

- Released source: `5a664680c91ad59f32b8c2a847b8fef37f34a0ae` (v3.1.0).
- Candidate source: `0c38d6a051082cc2a1eb961b4206837f82948fda`; production
  files match documentation-only dev `9105f29`. Existing clean-build
  receipts are reused and binaries reverified, not silently rebuilt.
- Raw discovery protocol SHA-256:
  `464f70848477c5580482d67963774b9b0182a82df0ab45d09e996f61ff7646af`.
- Registration precedes measurement in [issue #153](https://github.com/caseywdunn/sharkmer/issues/153#issuecomment-5637615559).
- Fixed k19 reference evidence is in the prior
  [high-copy archive](../v3.2-high-copy-final-20260911).
- Two physical cores, unpaired inputs, chunks=0, threading off, prewarmed
  uncompressed prefixes, 1,800 s timeout, 40 GiB virtual-address ceiling.
  K order rotates across pair/sample blocks; paired version runs are serial.
  All discovery timing finishes before any reference classification.
- BLAST/makeblastdb 2.17.0+ executable, database, source-result, binary, and
  source identities have before/after classification receipts. Actual-k
  command/stats/current-manifest checks pass at all four k values.

## Contents and integrity

- `overview.json`: compact machine-readable totals, selection outcome, and
  source-file hashes. The raw protocol hash and analyzer canonical-JSON hash
  are intentionally different representations, identified in their receipts.
- `evidence.tar.gz`: 1,352 original files, including all measurements and raw
  outputs, classified results and databases, complete/detailed analysis,
  protocol/build/helper receipts, and eight actual-k adapter smokes.
- `ARCHIVE_CONTENTS.json`: every archived file's size and SHA-256.
- `tools/`: reviewed driver, analyzer, tests, smoke, summary, reference-probe,
  preparation, and packaging scripts. Historical helper dependencies are
  frozen in the execution receipts inside the archive.
- `validation-tests.log`: ten driver and four analyzer tests pass.
- `reference-repeat-probe.json`: reference-only oriented-substring counts,
  not an actual read graph or read-supported repeat-copy validation.
- `development-race-note.md`: a reported transient unused-helper edit during
  approval handoff. Approved launch/checkpoint and all five frozen driver
  hashes agree; this does not assert continuous working-file immutability
  or an independently fingerprinted running Python bytecode image.

```bash
sha256sum -c SHA256SUMS
tar -tzf evidence.tar.gz
```

Large input reads, source exports, and executable binaries are not embedded;
their identities are receipted. Absolute paths record the original host and
are dependencies for direct script replay, not portable installation paths.
Use a new frozen protocol/output directory for another experiment. This sweep
does not establish held-out sensitivity, all-panel parity, deeper-input or
gzip/network performance, or a 16 GB laptop hardware guarantee.

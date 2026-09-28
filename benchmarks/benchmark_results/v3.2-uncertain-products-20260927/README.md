# Uncertain-product emission: implementation and comparison evidence

Completed 2026-09-28 for [issue #157](https://github.com/caseywdunn/sharkmer/issues/157)
on `experiment/uncertain-products`. All nine frozen comparisons pass:
**five of seven historical lost identities are emitted as explicitly uncertain
FASTA; ordinary sequences, counts, gene outcomes, and threshold diagnostics
remain unchanged.** The feature remains opt-in and awaits an integration decision.

Read the [results report](../../../dev_docs/UNCERTAIN_PRODUCTS_RESULTS.md) and
the pre-execution [output contract/protocol](../../../dev_docs/UNCERTAIN_PRODUCTS.md).
Re-emission makes historical graph-path hypotheses available for inspection;
it does not validate repeat lengths or add read evidence. Production
read-history reconstruction remains separate from this output feature.

## Results

| Sample | Ordinary products in baseline/off/on | Additional uncertain hypotheses |
| --- | --- | ---: |
| Gryllus SRR27962769 | 11 / 11 / 11 | 17 |
| Drosophila SRR31887760 | 9 / 9 / 9 | 4 |
| Heliconius SRR1057608 | 6 / 6 / 6 | 13 |

Exact historical matches are Gryllus CO1_1 373 bp and ITS_2 978 bp,
Heliconius ND1 263 bp, and Drosophila 12S 475 bp and 16S_2 590 bp.
The Drosophila 12S 487/499 bp identities are unobserved under bounded
search/retention. No additional historical identity matches only by reverse
complement, and none appears in ordinary output. Drosophila 12S retains three
of 26 observed withheld paths; the other 23 exceed the marker budget.

Machine-readable findings and their derivation:

- [summary.json](results/summary.json): combined assessment and provenance.
- [sample_runs.tsv](results/sample_runs.tsv): nine runtime/RSS observations.
- [uncertain_genes.tsv](results/uncertain_genes.tsv): all 63 sample/target rows,
  including empty assessments, lengths, and retention accounting.
- [historical_losses.tsv](results/historical_losses.tsv): all seven exact
  identity classifications, hashes, and bounds.
- [derive_summary.py](results/derive_summary.py): derivation from final receipts.

Single rotated trios are descriptive resource observations, not statistical
performance equivalence. These three full-panel samples do not replace the
earlier scoped 71-to-64 high-copy release benchmark.

## Accepted execution and checks

The comparison used two CPUs, k=19, chunks=0, the full frozen 21-target insect
panel, and the exact first-million R1 prefixes without read threading.
The wrapper records actual affinity, a 40 GiB address-space limit, and a
1,800-second timeout for each sequential invocation. Job `27727704` completed.
Pre-execution freeze SHA-256:
`ab139c02c2a14874fb66c41322268bdf5c690b837120c113f95f2ef2fe4cbb79`.
Source, runtime, input, and executable identities match their final bookends.

Both hash backends pass formatting, Clippy with warnings denied, and 246 unit
plus 26 integration tests each. Python ran 184 tests with one existing skip
and no failures; the reference audit verifies 156 references across ten panels.
The selected receipt is `checks-v3/checks.json`; earlier failing checks remain
preserved as history. Their two failures were CLI test-fixture argument errors,
corrected before the successful complete checks and candidate build.

[compiled_summary.json](validation/compiled_summary.json) identifies the actual
check/build receipts. [final_execution_review.json](review/final_execution_review.json)
binds the compiled prerequisites. The independent
[final results review](results_review/final_results_review.json) binds outputs
and final documentation. Independent parser, harness, packaging, and outcome
controls remain in their respective review directories; copies of independent
tests in CI are not additional independent evidence.

## Inputs and build identities

`prepare_inputs.py` extracted the frozen panel, seven historical identities,
and supporting provenance from the preceding evidence archive.
[inputs/provenance.json](inputs/provenance.json) binds that derivation.
[current_prefix_availability.json](inputs/current_prefix_availability.json)
records exact current input hashes and byte counts. Restoration receipts
preserve the successful downloads and failed first Drosophila attempt.
Partial ENA downloads verify consumed prefix bytes against archived SHA-256;
they do not verify the complete compressed upstream file's MD5.

| Accession | Current prefix path relative to repository |
| --- | --- |
| SRR27962769 | `data/its-recruitment-20260925/SRR27962769.first1M_R1.retry1.fastq` |
| SRR31887760 | `tmp/uncertain-products-20260927/reads/SRR31887760.first1M_R1.retry1.fastq` |
| SRR1057608 | `tmp/uncertain-products-20260927/reads/SRR1057608.first1M_R1.fastq` |

Pristine baseline commit:
`956e771437a2410a7bcbb1504b033c61dd1e11b6`.
Baseline source archive SHA-256:
`6566e3d92f7f14c3aa8f6966d2ab1340d6f5ece1453b2c3829edd075ec488937`.
Baseline executable SHA-256:
`bbf418d8e3b83f7eb7bf94f294f906a666fcf92babede2c46772a5d31e3caec6`.
Candidate executable SHA-256:
`ef6074b4e817dcafce209b9ccce40ab33c70a6ebb1c47dd3b6240acc3f6c3ca2`.
The exact candidate compiled-source snapshot is archived; its export receipt
is `attempts/builds/candidate-v1/snapshot.json`.

Both builds used project-local Rust/Cargo 1.98.1, release/default `ahashmap`,
locked offline dependencies, no Cargo configuration, and no custom Rust flags.
The concrete compiler identities and full command/environment are recorded
in `attempts/builds/{baseline-v3,candidate-v1}/build.json`.
Baseline-v1/v2 are superseded exports, not accepted compiled builds.

Activate the current local installation from the repository root:

```sh
source tmp/uncertain-products-20260927/toolchain/activate.sh
```

The activation sets isolated Cargo/Rustup homes and PATH; global shell setup
was not changed. [toolchain_install.json](validation/toolchain_install.json)
records installer integrity and versions. The system module's Rust 1.86
could format the code but could not compile existing baseline let chains.
Those probe results are retained separately from successful modern builds.

## Archive and reconstruction

The durable bundle uses `evidence.tar.gz`, `ARCHIVE_PLAN.json`,
`ARCHIVE_CONTENTS.json`, `ARCHIVE_RECEIPT.json`, and `SHA256SUMS`.
Verify external checksums, then verify the inventory's exact member set,
uncompressed sizes, and member hashes before using extracted evidence.
The archive has these principal groups:

- `comparison/`: freeze, nine complete receipts and outputs, final assessment,
  constraints, timings, and orchestration receipts.
- `candidate_source/`: exact compiled candidate inputs, including the small
  committed integration FASTQ fixture.
- `checked_source/`: additional checked source payloads, including the Python
  validator, bootstrap, audit, and test dependencies. The plan's
  `checked_source_members` maps all 85 checked original paths to archived
  members; files already present under `bundle/` reuse those members.
- `attempts/`: successful/failed checks, probes, build/export receipts, and
  Slurm logs. Presence of a failed attempt does not count as a pass.
- `bundle/`: inputs, harnesses, controls, derived results, and reviews.
- `documentation/` and `packaging/`: final documentation and packaging source.

External million-read FASTQ prefixes, all compiled executables, toolchain and
dependency caches, Cargo targets, generated panel reports, and nested older
benchmark archives are intentionally excluded. The baseline source tree/archive
is also excluded; reconstruct it from the recorded repository commit:

```sh
git archive --format=tar 956e771437a2410a7bcbb1504b033c61dd1e11b6 > baseline-source.tar
sha256sum baseline-source.tar
```

Check the archive hash above before extracting to a fresh directory. Rebuild
the selected sources with their Cargo.lock and recorded compiler/profile:

```sh
cargo build --locked --offline --release --message-format=json-render-diagnostics -j 2
```

Populate the locked dependency cache first if needed. Preserve recorded build
environments and verify rebuilt executable identities; a build on another
host/path must not be assumed byte-identical. Restored raw inputs also require
hash and record/base checks. Original absolute paths are provenance, not
portable execution paths or trusted identities.

Restore the checked-source mapping to its repository-relative layout when
reconstructing the Python checks and harness. Packaging attempt v1 passed byte
integrity but omitted these Python payloads; its receipt and omission review
are retained as history. Use the current archive receipt and independent
`ARCHIVE_REVIEW.json` for the accepted package.

For a new comparison, use fresh snapshot/check/attempt directories and the
reviewed `build_snapshot.py`, `run_compiled_checks.py`, and
`run_allocated_comparison.py` workflows. Load Python 3.12 and PyYAML
(`Python/3.12.3-GCCcore-13.3.0`, `PyYAML/6.0.2-GCCcore-13.3.0` here).
The original direct Slurm command is preserved in
[comparison_submission-v1.json](validation/comparison_submission-v1.json).
Do not overwrite original outcomes or reuse their success receipts for new
source/inputs. New runtime/source identities require fresh review and freeze.

Issue #157 remains open for integration. The read-supported reconstruction,
biological-loss, held-out validation, and release decisions remain separate.

# Independent offline harness review

The revised harness passes this source and synthetic-control checkpoint.
This is **not biological execution clearance** and does not establish that
the Rust implementation compiles or passes its integration tests.

`offline_harness_review-v2.json` binds the current reviewed source hashes
before and after the recorded controls, superseding the earlier receipt for
source binding. The independent harness suite passed 20 test methods with 33
table-driven subcases; the new independent orchestration suite passed eight
methods with 21 subcases. The author's separate 24-method harness and
eight-method orchestration suites also passed when rerun by this reviewer.
There were no skips. Repeated controls are not additional independent
evidence. No biological input was scanned and no assembler was invoked.

## Findings resolved before this checkpoint

- The builder now selects and fingerprints concrete compiler/toolchain
  executables, verifies the fixed release command/profile and exporter/archive
  lineage, captures Cargo configuration, rejects uncontrolled compiler/profile
  overrides, compares baseline and candidate configurations, and preserves
  probe/build failure receipts. Timeout handling kills the build process group.
  Source/tool/bookend verification failures also leave a failed receipt.
- Hashing rejects nonregular final paths and detects identity changes while
  opening, reading, and closing files. Generated Python cache files are excluded
  from candidate source exports; integration tests are included. Generated
  `panels/validation_runs`, `panels/validation_results`, and
  `panels/validation_reports` are excluded from exports/check inventories;
  actual panel YAML, reference sources, examples, and test fixtures remain.
- All additional reviewed files enter the frozen source map. Execution checks
  the review gate, mandatory source coverage, sample/build roles, sample prefix
  metadata, typed distinct two-CPU allocation, exact command order/flags,
  limits, and the externally supplied freeze hash.
- A frozen launcher measures actual address-space limits and CPU affinity
  before replacing itself with the assembler. Synthetic children confirmed
  both constraints survive `exec` with the same PID; wrong expectations prevent
  execution. These mechanism tests use the available test allocation, while
  the real comparison gate separately requires two CPUs.
- Capture checks each gene's declared count and lengths against its FASTA,
  rejects failed-gene standard publication and unexpected threading/hidden
  diagnostics, and validates uncertain sidecars and the output transaction.
  Optional singleton/multicopy k-mer counts participate in presence/value
  parity. Record order, complete headers, and pre-existing diagnostics remain
  part of standard equality.
- The compiled-check orchestration records the exact seven-command plan,
  Rust tools, Python/PyYAML runtime, source inventory, and shared Cargo
  configuration. The comparison wrapper rejects altered commands, feature
  flags, runtime/configuration, logs, or tools and requires passing compiled
  receipts from the current source. Timeout/spawn-failure calls retain
  immutable receipts. A real synthetic descendant was killed with its process
  group on timeout. Check log size is accurately described as a polled
  termination threshold, which can overshoot, rather than a strict file cap.

## Scope of the controls

The independent capture and frozen-plan unit controls mock earlier
manifest/build-validation stages to isolate the newly reviewed contracts.
They do not simulate successful Rust builds. The separate independent
validator suite exercises actual metadata, sequence, source-local geometry,
FASTA, and manifest checks; it is recorded in
`independent_validator_review.json`. Neither control set reconstructs a graph
or validates biological repeat length.

The harness retains separate exact-oriented, RC-only, and bounded
non-observation categories for the seven fixed historical losses. It does
not promote approximate matches or uncertain records into standard recovery.
The prior focused diagnostic allocation differed from the unchanged full
21-target allocation, so non-observation cannot be treated as biological or
graph-path absence.

## Remaining execution gate

Require successful compiled Rust checks, immutable baseline/candidate build
receipts with the same effective toolchain/configuration, final review coverage
of the concrete source/test snapshots, and all three exact consumed R1 prefix
identities. Freeze and execute within the intended two-CPU allocation, with
no concurrent build/test workload during timing. Review the resulting nine
run receipts, output ledgers, standard parity, and seven-loss classifications
before interpreting the results.

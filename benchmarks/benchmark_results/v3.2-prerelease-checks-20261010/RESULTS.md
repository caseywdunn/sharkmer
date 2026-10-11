# Final local validation

Job **29147703** completed successfully in **10 minutes 13 seconds** on
`c1102u05n01`. All eleven declared gates passed. No gate was skipped or rerun.
The final receipt confirms all **75 source/config/fixture files** and **178
tool/helper bindings** match their frozen pre-run hashes.

## Actual results

| Gate | Result |
|---|---|
| Default release build | Passed; executable copied before alternate-backend builds |
| Plant comparison | All eight invocations and four before/after comparisons passed |
| Default Rust backend | 251 unit + 23 sPCR integration + 6 uncertain-output integration tests passed |
| Fxhash Rust backend | Same 280 tests passed; zero failed or ignored |
| Clippy, both backends | All targets passed with warnings denied |
| Formatting | Passed |
| Python | 200 passed, zero skipped; real BLAST database/query test executed |
| Public references | 158/158 verified across nine built-ins and the example; no exclusions |
| Explicit JSON Schema | All ten panels passed with jsonschema 4.26.0 |
| Final-binary panels | Nine embedded exports match source panels; 19 loader/dry-run checks passed |
| Fixture oracles | k19 count and k31 exact standard-product identities passed |

The Rust suite also executes the existing independently established Jellyfish
histogram count regression on both backends. The additional fixture probes bind
the exact included FASTQ to `benchmarks/known_truth.yaml`: k19 reports 13,197,385
kmers, and k31 emits exactly the pinned standard 18S (1,783 bp) and 28S_2 (430 bp)
sequences. All fixture records and full SHA-256 identities are preserved. The
k31 count is reported in the raw evidence rather than asserted from an unbound
historical count.

The eight plant runs use the reviewed before/after panel snapshots and the same
frozen final executable. They cover two synthetic public-source inputs and the
two saved real input prefixes. Their biological/reference interpretation is
reported separately in the plant assessment; passing execution is not evidence
of full sample phasing.

## Frozen executable and limits

The default-feature binary is preserved at `results-v1/bin/sharkmer` (7,040,024
bytes), SHA-256:

```text
3a6505da2750982c333b077a065ad9e0cd0e6607c6963fa7140a0b44da71a899
```

It was copied immediately after the locked, offline default release build, then
rehash-checked through the remaining gates. Alternate-backend test builds cannot
replace this saved binary. The source freeze includes the corrected plant panel
and Cargo files with rand 0.8.6 and rustls-webpki 0.103.13.

The Slurm request was two physical CPUs, 16 GiB and two hours. Actual selected
CPUs were `[25,40]`, on distinct cores. Each outer gate ran under a 12 GiB address
space limit; the plant runner further restricted its synthetic invocations to
2 GiB and its real invocations to 12 GiB, with one CPU and a 600-second limit per
plant invocation. Gates ran sequentially, and Rust tests used one test at a time.

The largest GNU-time gate RSS was 4,357.13 MiB, during the plant comparison.
Slurm's batch memory accounting reported 5,483,320 KiB. These are different
measurements; neither is the address-space limit or the memory request. This
validation is not a final performance benchmark or a 16 GiB laptop certification.

## Evidence and reproduction

- `plan-v2.json` and `independent_plan_v2_review.json`: executable plan, source
  hashes, tool bindings and pre-run approval. `plan-v1.json` is an unused draft
  preserved with its overly broad historical-output inventory.
- `results-v1/result.json`: overall status, frozen binary, source checks and all
  gate receipts; `sources.before.json` and `sources.after.json` are identical.
- `results-v1/gates/`: complete commands, enforced constraints, GNU time metrics,
  stdout/stderr and receipts for each gate.
- `results-v1/plant/`: all eight plant inputs' commands, output files, metrics,
  identities, bounds and comparison results.
- `results-v1/{panels,references,fixtures}/`: explicit checks and outputs.
- `summary-v1/summary.json` and `gates.tsv`: recomputed counts and compact timing
  table, produced by `summarize_checks.py` without another test execution.
- `slurm-29147703-accounting.json`: terminal Slurm state and resource accounting.

`run_final_checks.sbatch` uses the project-local Rust toolchain and the
`sharkmer-bench` Python 3.14/BLAST environment, with the workspace jsonschema
installation on `PYTHONPATH`. Reproduction requires a fresh output directory.
No release, push or collaborator performance acceptance is implied by these
local validation results.

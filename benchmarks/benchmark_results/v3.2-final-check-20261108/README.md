# Focused v3.2 final-candidate check

Authorized 2026-11-08 after the user approved deferring the new held-out
biological collection and full nine-panel biological campaign from v3.2.
Those studies remain unfinished follow-up before further assembly tuning or
general accuracy claims. This check measures sequence compatibility and
resource use on existing calibration inputs, not biological correctness.

## Frozen scope

Use the exact previously tested `3.2.0-dev` candidate whose SHA-256 is
`7ef94fedafec49300146c0acf52086915fa78a82a8786dc11413e399a339dd9a`.
The source is unchanged from the November annotation verification and is
integrated through `d81f863`; later documentation edits do not alter it.
Bind the executable to that validation receipt, the current source files,
panels, schema, and the scripts used for this check. The baseline is pristine
v3.1.0 (`5a664680c91ad59f32b8c2a847b8fef37f34a0ae`), SHA-256
`e751e5eeca9150495e8d74d95bf27e4109466cad4ccf862e470fe271351d5899`,
with its preserved build and source receipts. Both builds use the pinned
Rust/Cargo 1.98.1 toolchain and default hash backend.

The fixed schedule contains 32 invocations:

- Three alternating baseline/candidate pairs for Human, Rhopilema,
  Liriodendron, and Heliconius: 24 timed invocations, always using the
  historical panel bytes for both versions.
- One candidate replay for each of the six other historical samples,
  completing the ten-sample sequence inventory.
- Separate candidate replays for Liriodendron and Acer with the corrected
  plant panel. These are functional checks against saved corrected outputs;
  their one-run timings do not establish a performance change.

The input inventory independently checks saved receipt outputs against the
tracked 86-record released inventory, 326-record historical-panel candidate
inventory, and 287-record corrected-panel inventory. Expected sequences,
classes, ordering, ordinary fields, target outcomes, and aggregate counts
are fixed before execution. Old `uncertain` and current `warning` labels
are normalized only for class comparison. The current format must separately
pass schema-3 header, structured-evidence, and publication-receipt validation.

## Execution and bounds

Use the original exact uncompressed input prefixes: 1,000,000 records for
nine samples and the EOF-verified 108,518-record Human input. The Liriodendron
prefix retains its original ordered R1/R2 concatenation; no paired-fragment
inference is added. Full byte hashes, record/base counts, and input sizes must
match before use. No input downloads or new reference collection are planned.

Run sequentially in one Slurm allocation with two distinct physical cores,
k=19, two threads, chunks=0, the one-million-record cap, default inclusive
output, and read threading off. Output uses fresh directories on the same
NFS filesystem. Record node/topology, mount details, actual child affinity,
address-space limit, commands, executable identities, and environment.

The batch requests 16 GiB RAM and 2h15m; each child has a 40 GiB virtual
address-space ceiling and a 900-second timeout. The runner stops starting
commands after 6,900 seconds. The Slurm allocation is not a physical-laptop
performance claim. Preserve partial outputs and failure receipts; never
replace an earlier attempt with a successful retry.

Read and verify inputs outside timing to provide the same warm-input intent.
Do not drop caches or claim that operating-system or NFS cache state is
controlled. No assembly warm-up or tracing is part of the timed measurement.
GNU time measures the full CLI wall/user/system time and process peak RSS;
keep allocator peaks distinct. Validation and hashing are outside the timer.

## Review criteria

Require all invocations to complete within bounds and all frozen sequence,
class, ordering, target, aggregate-count, and publication checks to pass.
Require repeated outputs to agree and source/protocol/input bindings to
remain unchanged. The current-panel inventory must retain all 86 historical
identities, with the expected 22 carrying warnings. Report additions and
warning fractions without treating them as validated biological products.

Report every paired measurement, absolute and percentage changes, medians,
and ranges. Preserve the earlier protocol's >5% median wall-time and >2%
median RSS investigation flags. These are review triggers, not automatic
acceptance limits or statistical significance tests. Inspect flagged costs
and record their disposition before recommending release. Three pairs on
one host do not establish equivalence across hardware or storage systems.

An independent reviewer checks raw output, receipts, source bindings, and
resource summaries. No sequence change is silently excused by the biological
study deferral. Outcomes and the performance disposition will be recorded
separately from this frozen protocol after execution.

## Reproduction

`prepare_inputs.py` reconstructs the small oracle from saved artifacts;
`prepare_protocol.py` binds it, the source, and the fixed schedule. The
preserved temporary input and executable paths must exist and match their
fingerprints; absent data must not be silently substituted.

From the repository root, after preparation and independent protocol review:

```sh
mkdir -p tmp/final-check-20261108/logs
sbatch benchmarks/benchmark_results/v3.2-final-check-20261108/batch.sh \
  --output "$PWD/tmp/final-check-20261108/run-v1"
```

Use a fresh output directory for every attempt. Retain the protocol's hash
and source commit with the Slurm submission and the resulting receipts.

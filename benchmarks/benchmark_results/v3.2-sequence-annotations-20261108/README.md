# Sequence annotation compatibility check

This bounded synthetic check compares the new schema-3 sequence annotations
with the frozen October 11 candidate. It does not evaluate biological
correctness or establish performance equivalence.

## Verified result (2026-11-08)

The [final receipt](validation/result.json) records Slurm job `32580418`,
two CPUs and 16 GiB requested, with unchanged source bindings and all eight
gates passing: 283 Rust tests per backend (253 unit + 30 integration), 213
Python tests with no skips, both Clippy checks, formatting, all 158 reference
regions, and the independent parity harness. Individual logs are retained
under `validation/`.

All [four paired modes](validation/parity/result.json) pass. Inclusive modes
retain one clean record and one warning-bearing record; opt-out modes retain
the clean record only. The sequence bytes, identifiers, ordering, counts,
warning membership, ordinary header fields, and full source evidence match
the reviewed baseline. Source/input/executable hashes and all eight invocation
receipts and outputs are retained. The candidate executable itself is not
duplicated in this bundle; its SHA-256 is in the receipts.

The [superseded first attempt](superseded-attempt.json) also passed every
individual gate, but its final source-binding check detected the two Python
files edited to fix legacy FASTA IDs containing `=` during the run. It is not
the accepted frozen-source validation. A [second check](prior-check-before-quoted-key-fix.json)
passed on frozen sources before final review found quoted annotation keys in
the identifier position could bypass the guard. The final run above covers
both identifier fixes, including quoted keys and unmatched quotes.
Slurm accounting reports the accepted job completed in 46 seconds, with
698,920 KiB peak batch-step RSS. These describe verification resource use,
not assembly performance equivalence.

To repeat the full checks, use the existing project toolchain/conda environment
paths in `batch.sh` and submit from the repository root:

```sh
mkdir -p tmp/annotation-conventions-20261108/logs
sbatch --partition=devel benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/batch.sh /absolute/fresh/results-directory
```

The runner
records commands, source bindings, log files, and results. Do not reuse an
existing results directory. The existing cached Cargo target is used serially.

## Parity protocol

The stored input follows the existing Rust CLI fixture in
`tests/uncertain_products.rs`: deterministic seeds 157 and 158 flank a
40-base A run; seed 159 produces an independent 160-base clean control.
Four identical full-length reads are supplied for each target. No reads
are downloaded, and no public reference is used as sample truth.

Run in the project's compute allocation with Python and PyYAML:

```sh
python benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/parity.py \
  --candidate /absolute/path/to/new/sharkmer \
  --output /absolute/path/to/new/run-directory
```

The default baseline is the preserved executable at
`benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/results-v1/bin/sharkmer`.
Its SHA-256 must be
`3a6505da2750982c333b077a065ad9e0cd0e6607c6963fa7140a0b44da71a899`.
Use a new output directory for each attempt; earlier attempts are preserved.

Four paired modes exercise inclusive output and the opt-out, each with and
without hidden diagnostics. Every invocation has a 60-second timeout.
Executable/input/harness hashes, commands, exit status, logs, output
receipts, and descriptive wall times are retained. Each run's normal stats
also retain its stage timings.

The check requires exact sequence bytes, order, identifiers, class counts,
warning membership, and all released ordinary-header tokens except the
removed `assembly_status=standard`. New warning headers must contain exactly
`sample`, `gene`, `product`, `length`, `assembly_status=warning`, and the
sorted unique `assembly_status_causes` list derived independently from the
source marker counts. The same causes must be recorded in structured evidence.

Full structured statistics are compared after removing only runtime
provenance (`command`, `run_id`, allocator peak, timings) and intentional
format changes (top-level output version/layout, uncertainty-assessment
version, product descriptions, and the new product cause list). Source
markers, positional evidence, hidden diagnostics, all limits, counters,
sequences, hashes, and other fields must match exactly. Both sets of
publication receipts must verify. The result is saved as `result.json`.

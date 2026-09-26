# Restored sequence hypotheses and controls

The six local classes and five full candidate sequences in `manifest.json`
come unchanged from the archived SCC audit. The archive and all 205 regular
members pass the original `SHA256SUMS` and `ARCHIVE_CONTENTS.json` checks.
Seven byte-identical source files are preserved under `originals/`, including
graph class memberships, the unchanged CO1 control, and prior input identity.
The full restore is in `tmp/its-recruitment-20260925/prior` outside durable Git
evidence. These checks inspect archived hypotheses only; they do not scan new
biological reads or establish biological correctness.

`provenance.json` records each source member, source candidate/event,
zero-based half-open interval, sequence identity, and self-context. In
particular, `region_1_class_1` is the reverse complement of historical
ITS `[72,133)`; `region_2_class_1` is the forward ITS `[255,317)` interval.
The extra-cycle classes have no exact full-ITS self-context. The CO1 local
windows are retained `[289,349)` and lost `[289,370)`, both forward.

The finite decoys are the archived 978-base ITS, 534-base 16S, 423-base 12S,
352-base retained CO1, and 373-base lost CO1 candidate sequences. They are
hypotheses, not independent reference sequences or genome-wide competitors.

`verification_receipt.json` records an independent verification using
`verify_inputs.py`, which does not import the derivation helper. It verifies
every archive and restored member; original copies; manifest/provenance/helper
hashes; all six literal classes; historical self-contexts; the five decoys;
and all eight archived pre/post-prune graph route memberships. Each graph
route is independently spelled from its overlapping 18-base node identities
and compared in its recorded canonical orientation.

`derivation_reproduction_receipt.json` records a fresh execution of
`restore_derive_inputs.py` in `/tmp`: the derived manifest and seven original
copies are byte-identical to this directory. Paths in regenerated provenance
will reflect the chosen destinations.

From the repository root, use fresh output directories to reproduce:

```bash
python3 benchmarks/benchmark_results/v3.2-its-recruitment-20260925/inputs/restore_derive_inputs.py \
  --archive-dir benchmarks/benchmark_results/v3.2-its-scc-audit-20260913 \
  --restore tmp/its-input-reproduction/prior \
  --output tmp/its-input-reproduction/inputs
```

To independently recheck the published input bundle and existing full restore,
choose a receipt path that does not yet exist:

```bash
python3 benchmarks/benchmark_results/v3.2-its-recruitment-20260925/inputs/verify_inputs.py \
  --archive-dir benchmarks/benchmark_results/v3.2-its-scc-audit-20260913 \
  --restore tmp/its-recruitment-20260925/prior \
  --inputs benchmarks/benchmark_results/v3.2-its-recruitment-20260925/inputs \
  --receipt /tmp/its-input-verification.json
```

Both helpers refuse to overwrite their outputs. Restoring the first million
R1 FASTQ records and verifying their exact byte identity is a separate step;
neither this bundle nor its receipts attest that the FASTQ is staged.

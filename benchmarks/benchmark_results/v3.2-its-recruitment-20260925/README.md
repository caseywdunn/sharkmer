# Gryllus ITS variation-aware read recruitment

Started 2026-09-25; frozen and executed 2026-09-26 on
`experiment/its-read-threading`, baseline `956e771`.
See [findings](../../../dev_docs/ITS_READ_RECRUITMENT_RESULTS.md) and
[protocol](../../../dev_docs/ITS_READ_RECRUITMENT_PROTOCOL.md).
Production assembly and repeat gates remain unchanged; #148/#153/#129 stay open.

The exact original first million Gryllus R1 records recruit 27 ITS1 and 28
ITS2 records. Neither historical short class has a Q20-qualified alignment.
ITS1 has one low-quality longer-class alignment. ITS2 has two records aligning
to the longer class, but their complete read sequences and unfiltered
placements are identical. One qualified span excludes a Q11 endpoint through
an additional terminal deletion. These are one observed sequence pattern,
not two independent footprints or proof of a full product. The unchanged CO1
control reproduces 27 exact / 28 sensitivity retained records and zero lost.

## Evidence layout

- `freeze.json` binds protocol, sources, input identities, independent review,
  and exact argv before new biological scanning. Its SHA-256 is
  `cb68a200f410a89980380cf7d4c290c3d825c7b501f43acf16736ec954abe600`.
- `inputs/` preserves source manifests, six derived local classes, all five
  finite-context competitors, verified self-context coordinates/orientations,
  reproduction helpers, and verification receipts for all 205 prior members.
- `review/` preserves pre-outcome synthetic and independent algorithm/CLI/
  execution-wrapper checks. The documented 256-to-4,096 pair-cap amendment
  was driven by synthetic repeated reads before biological outcomes.
- `restoration/` preserves both the sandbox DNS failure and the successful
  byte-identical prefix download. `restore_prefix.py` reproduces that prefix.
- `results_review/` contains the independent post-run checks and interpretation.
- `evidence.tar.gz` includes all of the above under `experiment/`, complete
  raw scan reports and recruited-read ledgers under `run-v1/`, the frozen
  source snapshots, command receipts/logs, and repository documentation
  snapshots. Raw FASTQ, binaries, and Python caches are excluded.

The directory keeps its start date. Original absolute paths in receipts are
provenance, not input identities. The 365,356,898-byte prefix hash is
`6aef48e6ed467e7f1ff25a3cd2a21d51dc2378985bdd855e796b539f0fc0c794`;
this prefix-only restoration is distinct from the older unused 8M-record
R1/R2 stage. Nothing here adds R2 or new coverage to the scientific assay.

## Verification and rerunning

Run `sha256sum -c SHA256SUMS` in this directory. `ARCHIVE_CONTENTS.json`
binds every archived member's size and SHA-256. Extract into a fresh location
and verify every member before use; never overwrite original outcomes.

For reproduction, restore the prefix with `restore_prefix.py --output` and
a fresh destination, verify the recorded prefix hash, and restore the source
versions from the archive. `run_experiment.py --freeze <new.json> --fastq
<restored.fastq> --run-dir <fresh-directory>` makes a new source/input/argv
freeze; execute it with `--freeze <new.json> --execute` on a compute node.
Use one CPU, 40 GiB memory, and at least two hours of job wall time to cover
the three separately bounded commands. A moved checkout produces a new
attestation with new paths; it does not rewrite the archived one. Retain
failures and do not tune biological thresholds against these outcomes.

All 29 diagnostic tests, independent distance/seed/context checks, and 13
wrapper tests pass. The three scans completed with matching source/input
bookends. Diagnostic timings do not measure assembly performance or release
parity. Integration requires a separately reviewed production proposal and
fresh high-copy regression comparison; no merge or release is approved here.

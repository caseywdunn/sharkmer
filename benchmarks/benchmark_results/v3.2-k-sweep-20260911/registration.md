## Preregistered longer-k comparison — 2026-09-11

User authorized testing k=19/23/27/31 before revising repeat policy. No release,
production code, primer-panel, or default-k change is authorized.

- Discovery: 72 invocations = three affected insect samples x four k values x
  two versions x three paired repeats; identical frozen 1M-record prefixes and
  unchanged whole panels. Same prior two-core, unpaired, threading-off, chunks=0,
  warmed-local-input, 1,800 s/40 GiB address-space protocol except k.
- Reuse verified pristine released-v3.1.0 (`5a66468`) and reviewed candidate
  (`0c38d6a`) binaries; candidate production files match current dev `9105f29`.
- Rotate k order across pair/sample blocks; keep same-k version pairs contiguous
  and alternate their order. All discovery timings precede classification.
- Require actual command/stats/manifest k, same-k aggregate count parity,
  immutable input/build identities and pinned classification tool/DB receipts.
- Evaluate absolute exact sequence retention against BOTH released-k19 and
  current-k19 outputs, not merely same-k parity or total product counts. Track
  all seven lost high-copy products, especially the exact 978 bp Gryllus ITS_2.
  Full-sequence changes remain visible; any demonstrated boundary equivalence
  is auxiliary and does not count as exact rescue.
- Choose at most one k>19 only if valid stable results restore at least one of
  seven exact missing sequences in every replicate and preserve every existing
  k19 candidate reference-confirmed sequence in the affected samples. Rank by
  exact ITS_2 recovery, number of seven recovered, released-k19 high-copy
  retention, then lower k. No timing-driven selection or same-loss tie with k19.
- If eligible, perform 42 confirmation invocations (three paired repeats for
  seven remaining non-insect samples). If none is eligible, stop at discovery.
  Every remaining high-copy loss and resource increase still requires review.

Machine-readable protocol SHA-256:
`464f70848477c5580482d67963774b9b0182a82df0ab45d09e996f61ff7646af`.
Findings will append to `dev_docs/BENCHMARK_high_copy_followups.md` and the
evidence archive; `dev_docs/PLAN.md` tracks completion. Sol prepares measurement,
Terra analysis, Astra independent review. Historical data remain calibration,
not held-out truth. #153 remains open until user review dispositions its gates.

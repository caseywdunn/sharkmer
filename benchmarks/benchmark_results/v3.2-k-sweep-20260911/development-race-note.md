# Measurement driver checkpoint provenance

The driver approved for this sweep has SHA-256
`6c49f74ca76187ce1c83a1c40b45803219bc80c86451fc57b191e199be481046`.
Root's launch command printed this digest before starting Python. Root and
Astra subsequently verified that the master frozen driver, all four child
frozen drivers, and the current working driver also have this exact digest.

After launch, Sol reported that an unused `primer_boundary_audit` helper had
briefly been added to the working driver around the approval handoff and then
removed, restoring the approved bytes. Sol reports that no invocation logic
or result file was changed. The transient version was not retained, so this
statement is an implementation report, not an independently reconstructed diff.

The hashes establish checkpoint and frozen-receipt agreement, not continuous
working-file immutability or an independent fingerprint of the running Python
bytecode. Executable/input/result checks remain in the normal run receipts.
Astra's bounded provenance review finds no basis to discard or replace these
measurements; they are preserved without rerunning. This note narrows the
driver attestation and records the reported race rather than silently omitting it.

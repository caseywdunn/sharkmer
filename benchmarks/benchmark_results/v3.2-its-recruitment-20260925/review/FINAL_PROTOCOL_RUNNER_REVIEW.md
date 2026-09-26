# Final protocol and execution-wrapper review

Date: 2026-09-26. This supplement resolves the pre-execution protocol conditions
in `INDEPENDENT_REVIEW.md`; the earlier report and binding remain unchanged.
No biological FASTQ records or new biological outcomes were examined by this
reviewer or the independent wrapper-test agent.

The final protocol has SHA-256
`21e713458e802d2a60cdc2c42873ab9de82706aa1bb25b4b4e65548c154c856d`.
It now specifies 4,096 ordered pairs, the independent synthetic 315-pair rationale,
per-region local uniqueness, and all implemented resource/evidence limits.
The scanner remains the previously reviewed source:
`5c6268e840a626eb980713da8b0a1992030c50bad0e258241dcce703b3f0a75d`.
These changes satisfy the prior protocol conditions without changing the
biological sensitivity, quality thresholds, or interpretation scope.

The execution wrapper has SHA-256
`aa5151784d0b08b9242b9e9bd36f319c16602692d95c473c128fdfcbf4280743`.
Review fixed gaps that could lose failure receipts after verification or launch
errors, and prevented `/usr/bin/time` from overwriting existing timing artifacts.
JSON publication now uses an fsynced temporary file and an atomic no-clobber
hard link. Fingerprinting rejects files changed during hashing. Each command
records its actual arguments, exit/timeout state, elapsed time, input/source
checks, and output hashes. The aggregate execution receipt also retains command
outcomes if a per-command receipt cannot be published. Input/source/freeze
verification failures stop subsequent commands. Timeouts kill the child's process
group and preserve a failed receipt rather than treating partial output as success.

The wrapper applies the frozen address-space cap to subprocesses and retains
the 1,800-second limit per command. The allocation must cover all three commands
plus hashing and publication; two hours is sufficient for the declared maxima.
Uncatchable termination such as SIGKILL still requires preservation of the
external scheduler/job logs and whatever artifacts already exist.

Independent synthetic wrapper controls are recorded in
`runner_synthetic_results.json`; they use only temporary generated inputs and
commands. The final execution review binding names their exact script and
results hashes. Compilation and whitespace checks also pass. The unchanged
29 scanner tests and earlier independent DP/seed/CLI evidence remain applicable.

The reviewed protocol, scanner, wrapper, input provenance, test evidence, and
actual argument vectors must be frozen together before the first biological
scan. Execution clearance remains limited to this diagnostic; local evidence
cannot establish repeat copy number, full-product phase, genome-wide uniqueness,
or eligibility for production admission. Unchanged CO1 compatibility results
must still be checked before interpreting the new biological output.

# First HTTPS harness attempt

Job 32671996 failed before any HTTPS request. The harness expected bare version
text, while the unchanged application correctly returned
`sharkmer 3.2.1 (https://github.com/caseywdunn/sharkmer)`. This is a harness
assertion error, not an observed application or TLS regression.

The raw job log, execution receipt, and version-command outputs are retained
byte-for-byte. The receipt identifies the executed script by SHA-256
`f14064e82e3d41aa45f4acf9272aea87fce830d4df4e5a9dade04aca5b9e2ba1`.
That mutable script was edited before archival, so its exact executed source
is not retained here. A later script copy is not presented as the failed
attempt's source. The corrected rerun uses a separately frozen script.

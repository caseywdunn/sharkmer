# HTTPS compatibility — successful frozen rerun

Job **32674065 passed** on 2026-11-09 in two seconds, using the same executable
SHA-256 as the complete engineering suite. The frozen script SHA-256 is
`831f498eaf0fb14a8952e13588d1929ba73c180de070ce4ea544080b49d56c6d`.
Its source, batch launcher, scheduler receipt, raw log, and all small outputs
are retained here. The job used two CPUs and a 16 GiB limit; batch peak RSS was
27,160 KiB. These resource receipts are not a performance benchmark.

Sharkmer loaded a checksum-pinned panel from an immutable GitHub HTTPS URL
using TLS 1.3 and normal certificate verification. Its bounded local-read run
matches the local-panel control: three reads, 300 bases, and 246 accepted k-mer
occurrences. Both runs completed with valid manifests. A local server with an
untrusted certificate was rejected with `UnknownIssuer` before any HTTP request.
The public test certificate is retained; its temporary private key was deleted
by the harness and is not archived.

This exercises HTTPS panel loading through the updated TLS path, bounded local
ingestion, and certificate rejection. It does not test ENA FASTQ download/cache
behavior, the malformed-handshake advisory itself, or biological accuracy.
The [first harness attempt](../https-attempt-1/README.md) failed a version-output
assertion before any network request and is retained separately.

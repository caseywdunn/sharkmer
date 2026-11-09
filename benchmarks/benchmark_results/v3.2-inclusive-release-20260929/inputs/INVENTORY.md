# Historical benchmark input inventory

`inventory.json` pins the original final 60-invocation scope, consumed FASTQ
identities, archived source order and five assembly panels. The panels are exact
`ba64f573048b6c19a528278028810af5bd475b81` blobs. All 52 primer mappings also
equal the current panels; full YAML bytes differ because reference metadata has
changed. Both fresh binaries must receive the pinned historical assembly panels.

Three insect prefixes are already staged at the paths in `current_file`. This
inventory inspected their sizes and prior successful receipts without rescanning
them. Seven prefixes need restoration. The compute-node helper verifies raw
structure, record/base/byte counts and SHA256 for all ten, including reused files.

Source order is mechanically derived from archived `source_consumption`:

* Human (`SRR17535371`): complete unsuffixed, R1 and R2 streams in that order;
  30,024 + 39,247 + 39,247 = 108,518 records, with EOF for each.
* Liriodendron (`SRR25378184`): complete R1 (656,431 records), then the first
  343,569 R2 records, totaling one million.
* Remaining samples: the first one million records of the first source.

Full compressed SHA256/MD5/size is checked only for sources consumed through
required EOF. Partial final sources are bound by the exact uncompressed prefix,
without claiming complete upstream archive verification. The rejected initial
Liriodendron plan is preserved in `inventory.prereview-v1.json`; use only
`inventory.json` (SHA256
`0e6e6edbebc7430c70519c37fd2bcdd8c2312ddc1b977dfb8b83200e5c1e9ddc`).

`restore_inputs.sbatch` requests two CPUs, 48 GiB and four hours. The helper runs
outside benchmark timing, allows 1,800 seconds per sample and at most 512 MiB of
compressed input per source, and refuses existing attempt paths. Every sample
gets a receipt; failed partials remain available for inspection. Successful
`availability.json` rows identify the actual verified FASTQ paths. The author
controls in `restoration_controls.json` passed 16 synthetic methods without
network access or biological inputs.

`historical_identities.json` authenticates archived FASTAs and establishes exact
ordered per-target sequence parity across all three historical replicates. It
contains 86 baseline records (71 high-copy and 15 deferred Gryllus Yp2 records)
and 64 previous candidate records (all high-copy). Literal and strand-canonical
SHA256 identities are provided; these are identity expectations, not independent
biological validity claims. `target_scope.json` preserves the frozen grouping.

Released v3.1.0 peels to
`5a664680c91ad59f32b8c2a847b8fef37f34a0ae`; its historical build and the new
local toolchain both report Rust/Cargo 1.98.1. All 238 locked registry package
filenames were present in the local cache during inventory. Fresh source exports
and builds still require their own receipts. Historical source exports, binaries
and original `/tmp` staging directories are unavailable.

No BLAST run is needed for this identity/resource comparison. If later requested,
the available module is `BLAST+/2.16.0-gompi-2024a`, unlike the historical 2.17.0
executable. Current verified-public-reference classification must pin its own
catalog, database and binaries and must not rebrand old bootstrap labels as truth.

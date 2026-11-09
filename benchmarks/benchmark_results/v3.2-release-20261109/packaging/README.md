# Local Bioconda build of v3.2.0

Slurm job **32662859 completed successfully** on 2026-11-09: 12 m 47 s,
two CPUs, and 8,013,864 KiB peak batch RSS (about 7.64 GiB). The build used
conda-build 26.9.1 and conda Rust 1.99.1 on linux-64. All four installed-command
smoke tests passed: `sharkmer --help`, `--version`, `--list-panels`, and `--cite`.
This does not test the declared Rust 1.88 minimum or other platforms.

The exact published v3.2.0 source archive, recipe, environment, build scripts,
job accounting, and logs are bound by [build-receipt.json](build-receipt.json).
The source SHA-256 is
`858dac22ea1d05a98a4917706de95e1012385e3d10078f42f7f1f98b0d883e0b`.
The resulting `sharkmer-3.2.0-h9bf148f_0.conda` package is 2,136,614 bytes with
SHA-256 `66c1ec7e4621dca04d0fc855e9c18c880cd8653ebb4a1d6f174d664faa3f6d0c`.
The package hash and size were checked before retaining these receipts; the
binary package and dependency caches are not copied into this small bundle.

The batch script disabled uploads. Successful local testing does not imply
that upstream Bioconda builds have passed or the package is available through
the channel. [Publication follow-up](../PUBLICATION.md) tracks that status.
This package check does not repeat the biological comparison or change its
held-out/all-panel deferrals.

Verify the retained artifacts from this directory with
`sha256sum -c SHA256SUMS`. This directory has its own checksums and does not
alter the frozen engineering-check bundle one level above. Original absolute
paths in the scripts and logs record the execution environment; they are not
portable commands without adapting the workspace location.

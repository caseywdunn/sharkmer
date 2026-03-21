# Contributing

## Development

After cloning the repo, set up the git hooks:

    git config core.hooksPath .githooks

This enables:
- **Pre-commit**: runs `cargo fmt --check` (instant)
- **Pre-push**: runs `cargo clippy -- -D warnings` and `cargo test`

Some common rust tasks in development:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release

### Test data

This repository includes a test [dataset from Thermus thermophilus](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5324768&display=metadata) for simple tests and to develop against.

After cloning the repo, gunzip the `data` in the data dir:

    cd sharkmer/data/ # Note that this is the sharkmer folder within the sharkmer repository
    gunzip -c SRR5324768_pass_1.fastq.gz > SRR5324768_pass_1.fastq

## Branching model

### Branches

- `master` — tagged releases only. Never commit directly to master.
- `dev` — active development for the next unreleased version. All feature
  work targets this branch.
- `vN` (e.g. `v2`, `v3`) — release maintenance branches, created from the
  release tag when a major version ships. Used only for patches to released
  versions.
- Issue branches — created from `dev` (or from `vN` for patches), named after
  the issue they address.

### Normal workflow

1. Create an issue branch from `dev` (`git checkout -b issue-NNN dev`)
2. Do the work, ensure it passes quality gates
3. Merge to `dev` and push (`git checkout dev && git merge issue-NNN && git push`)
4. Delete the issue branch
5. When ready to release: merge `dev` to `master`, tag the release, create
   a `vN` branch from the tag

Always merge and push completed issue branches to `dev` before starting the
next issue. This keeps `dev` up to date and avoids dependency tangles when
later issue branches need earlier work.

### Patching a released version

If a bug is found in v2.0 while v3.0 is in development on `dev`:

1. Create an issue branch from `v2`
2. Fix the bug
3. Merge to `v2`, tag `v2.0.1`
4. Merge `v2` to `master`
5. Cherry-pick the fix into `dev` so the next version gets it too

### Quality gates

A branch is considered passing when:

- `cargo test` passes
- `cargo clippy` reports no warnings
- `cargo fmt --check` reports no changes needed

## Bioconda recipe

Bioconda release follows their [contribution workflow](https://bioconda.github.io/contributor/index.html).

The bioconda recipe files are in `sharkmer/` (`meta.yaml` and `build.sh`). These are the source of truth and are copied to the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) repository for releases.

Before submitting a pull request, test the recipe locally with:

    # first time only
    conda create -n bioconda-test -c conda-forge -c bioconda bioconda-utils 

    # each test
    conda activate bioconda-test
    cd bioconda-recipes
    conda build recipes/sharkmer -c conda-forge -c bioconda
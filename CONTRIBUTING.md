# Contributing

### Installation from a specific branch on the repository

This allows you to test pre-release features.

Install sharkmer from a specific branch on the repo:

    cargo install --git https://github.com/caseywdunn/sharkmer.git --branch dev

Re-run the install command to pull and build the latest changes:

    cargo install --git https://github.com/caseywdunn/sharkmer.git --branch dev --force

To remove sharkmer:

    cargo uninstall sharkmer

## Development

After cloning the repo, install the pre-commit hook:

    ln -sf ../../scripts/pre-commit .git/hooks/pre-commit

This runs the same quality gates as CI before every commit: `cargo fmt --all --check`,
`cargo clippy -- -D warnings`, and `cargo test --release`.

Some common rust tasks in development:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release

### Test data

Integration tests use a 100k-read fixture from ERR571460 (Porites lutea),
stored gzipped at `tests/fixtures/ERR571460_100k_R1.fastq.gz`. No manual
setup is needed — `cargo test` handles it automatically.

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

- `cargo test --release` passes
- `cargo clippy -- -D warnings` reports no warnings
- `cargo fmt --all --check` reports no changes needed

## Releasing

These steps assume all work is complete on `dev`, quality gates pass, and
regression benchmarks have been run.

### 1. Run and commit benchmarks

Before preparing the release, run the full regression benchmark suite on `dev`
and commit the results:

    conda activate sharkmer-bench
    python benchmarks/run_benchmark.py

Review the results in `benchmarks/benchmark_results/` and compare against prior versions
with `benchmarks/compare.py`. If there are regressions, fix them before
proceeding. Commit the benchmark results to `dev`:

    git add benchmarks/benchmark_results/
    git commit -m "Add benchmark results for vX.Y.Z release"

### 2. Prepare the release

- Verify `version` in `Cargo.toml` matches the intended release (e.g. `3.0.0`).
- Verify `CHANGELOG.md` has an entry for this version with the correct date.
- Verify `meta.yaml` version and sha256 are updated (sha256 can only be
  finalized after the GitHub release tarball exists — use `PLACEHOLDER` until
  then).
- Update `sharkmer_version` in all `panels/*.yaml` changelog entries to match
  the release version.

### 3. Merge to master and push

    git checkout master
    git merge dev
    git push origin master

### 4. Tag the release

    git tag -a vX.Y.Z -m "vX.Y.Z"
    git push origin vX.Y.Z

### 5. Create the GitHub release

Go to <https://github.com/caseywdunn/sharkmer/releases/new>, select the tag
you just pushed, and create a release. Use the CHANGELOG entry as the release
notes.

Alternatively, use the CLI (replace `X.Y.Z` with the version number):

    VERSION=X.Y.Z
    gh release create "v${VERSION}" --title "v${VERSION}" --notes-file - <<< "$(sed -n "/^## \[${VERSION}\]/,/^## \[/{ /^## \[${VERSION}\]/d; /^## \[/d; p; }" CHANGELOG.md)"

### 6. Create the release maintenance branch

    git checkout -b vN vX.Y.Z   # e.g. git checkout -b v3 v3.0.0
    git push origin vN

This branch is used for future patch releases (vX.Y.1, etc.) without
pulling in unreleased work from `dev`.

### 7. Update the bioconda recipe

After the GitHub release is created, get the sha256 of the source tarball:

    curl -sL https://github.com/caseywdunn/sharkmer/archive/refs/tags/vX.Y.Z.tar.gz | sha256sum

Update `meta.yaml` with the real sha256, then follow the bioconda submission
steps in the Bioconda section below.

### 8. Resume development

    git checkout dev

Update `Cargo.toml` version to the next development version (e.g. `3.1.0-dev`)
and add an `[Unreleased]` section to `CHANGELOG.md`.

## Regression benchmarks

A benchmark suite in `benchmarks/` runs sharkmer against 14 real-world SRA
datasets across multiple primer panels. Run after each development phase or
before a release to check for regressions.

Set up the benchmark conda environment (first time only):

    conda init zsh      # or bash — writes conda hooks to your shell rc file
    # Restart your shell (or open a new terminal) for the hooks to take effect.
    # In a fresh shell that has not yet been restarted you can source the hooks
    # directly instead:
    #   source /opt/conda/etc/profile.d/conda.sh
    conda env create -f benchmarks/environment.yaml

Then activate it before running benchmarks:

    conda activate sharkmer-bench

    # Run the full benchmark (builds sharkmer, downloads data if needed, runs all samples)
    python benchmarks/run_benchmark.py

    # Run specific samples only
    python benchmarks/run_benchmark.py --samples Porites_lutea Agalma_elegans

    # Pre-download all sample data without running benchmarks
    python benchmarks/run_benchmark.py --download-only

Results are written as YAML to `benchmarks/benchmark_results/` with the date, version,
and git commit in the filename. Use `benchmarks/compare.py` to diff results
across versions.

Sample data (~1M reads each) is cached in `benchmarks/data/` and downloaded
from ENA on first run. The download streams and truncates early, so it does
not fetch full runs. Configuration (samples, panels, read counts) is in
`benchmarks/config.yaml`.

## Bioconda recipe

Bioconda release follows their [contribution workflow](https://bioconda.github.io/contributor/index.html).

The bioconda recipe files are at the repo root (`meta.yaml` and `build.sh`). These are the source of truth and are copied to the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) repository for releases.

Before submitting a pull request, test the recipe locally with:

    # first time only
    conda create -n bioconda-test -c conda-forge -c bioconda bioconda-utils 

    # each test
    conda activate bioconda-test
    cd bioconda-recipes
    conda build recipes/sharkmer -c conda-forge -c bioconda
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

This runs the baseline Rust checks before every commit: `cargo fmt --all --check`,
`cargo clippy -- -D warnings`, and `cargo test --release`. The complete
release check below also covers both backends, Python, panels, and metadata.

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
3. Push the issue branch and open a pull request into `dev`; merge after CI passes
4. Delete the issue branch
5. When ready to release, follow the release ritual below, including updating
   the existing `vN` maintenance branch or creating it for a new major version

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

**v3.2 scope exception approved by the user, 2026-11-08:** follow the retained
[release assessment](dev_docs/REVIEW_v3.2.md#decision-ownership-2026-11-08) and
[final-build results](benchmarks/benchmark_results/v3.2-final-check-20261108/RESULTS.md)
in place of the full benchmark and biological sample campaigns in steps 1–3.
Use the retained release-baseline evidence and focused check of saved sequences,
warning classifications, and the four flagged performance comparisons; keep
the intentional plant-assay change separate. The new held-out dataset and
full nine-panel biological evaluation are deferred, not passed. Quality gates,
reference-provenance checks, and schema/export/loader checks for all built-in
panels and the reference example remain required. Retain the deferred biological
work before further assembly tuning or broad accuracy claims. The general
procedure below continues to apply to other releases.

### 1. Run and commit benchmarks

Before preparing the release, run the full regression benchmark suite on `dev`
and commit the results:

    conda activate sharkmer-bench
    python benchmarks/run_benchmark.py

Review the results in `benchmarks/benchmark_results/` and compare the structured
per-product YAML records, classifications, and performance fields with prior
runs. There is no separate comparison script. Investigate regressions before
proceeding. Commit the benchmark results to `dev`:

    git add benchmarks/benchmark_results/
    git commit -m "Add benchmark results for vX.Y.Z release"

### 2. Validate all built-in panels

Run the full panel validation suite to confirm all built-in panels load,
parse, and produce amplicons against their declared validation samples.
On a local machine:

    conda activate sharkmer-bench
    bash scripts/validate_all_panels.sh

On the Bouchet HPC cluster (submits one SLURM job per panel):

    ml Rust GCC
    cargo build --release
    bash scripts/validate_all_panels_slurm.sh

If any panel fails, fix the discrepancy before proceeding. Pass `--no-blast`
to skip BLAST validation if a local BLAST database is unavailable.

### 3. Validate the reference panel

`panels/examples/reference.yaml` is a documentation panel that exercises
every schema v2 feature — multiple primer indices, deprecated entries,
compartment/gene_type/copy_number metadata, and so on. It is not embedded
in the binary but must parse and validate cleanly before every release:

    conda activate sharkmer-bench
    python scripts/validate_panel.py panels/examples/reference.yaml

If the reference panel fails to load or validate, the schema implementation
or documentation is out of sync. Fix the discrepancy before proceeding.

### 4. Check readiness and prepare the release

Confirm readiness gates in `dev_docs/PLAN.md` are met or explicitly waived;
complete the final packaging/publication gate through the steps below.
Activate the benchmark environment (which supplies Python, PyYAML, jsonschema,
and BLAST+) and use a Rust toolchain with rustfmt and Clippy. Run the same
release check used by CI on pull requests into `master`:

    git fetch origin
    git switch dev
    git merge --ff-only origin/dev
    conda activate sharkmer-bench
    scripts/release-check

The check runs formatting, Clippy and Rust tests with both backends, the full
Python suite, the offline reference audit, all nine built-in panels and the
reference example (schemas, exports, and loader dry runs), and the pinned
known-answer CLI fixture. It checks Cargo/lockfile/changelog agreement and,
for release versions, the recipe version. It does not run the biological
campaigns in steps 1–3. Use `--offline` when Cargo dependencies are cached;
`--only` selects individual gates for diagnosis. Run the full check in a Slurm
allocation on the cluster, not in the coding-agent allocation.

Prepare a release branch from the checked `dev` (the examples use v3.2.0):

    VERSION=3.2.0
    git switch -c "release-${VERSION}" origin/dev

Set the package version in `Cargo.toml` and `Cargo.lock` to `${VERSION}` without
`-dev`. Rename the changelog's `[Unreleased]` heading to
`[3.2.0] - YYYY-MM-DD` with the actual release date. Set `meta.yaml` to the same
version and use `PLACEHOLDER` for its checksum until the tagged archive exists.
Run `scripts/release-check` again, then commit these changes. Save the release
changelog section for the pull requests and GitHub release:

    awk -v heading="## [${VERSION}] -" 'index($0, heading)==1 {printing=1; next} /^## / && printing {exit} printing {print}' CHANGELOG.md > /tmp/sharkmer-release-notes.md
    git push -u origin "release-${VERSION}"
    gh pr create --base dev --head "release-${VERSION}" --title "Prepare ${VERSION}" --body-file /tmp/sharkmer-release-notes.md
    gh pr checks "release-${VERSION}" --watch
    gh pr merge "release-${VERSION}" --merge

### 5. Merge the release into master

Open the release pull request from `dev`. After CI passes and release is
authorized, merge with a merge commit so the development history is preserved:

    gh pr create --base master --head dev --title "Release ${VERSION}" --body-file /tmp/sharkmer-release-notes.md
    gh pr checks dev --watch
    gh pr merge dev --merge

The user already authorized the v3.2 release ritual; that decision is recorded
in the retained release assessment linked above.

### 6. Tag the exact released commit

    git fetch origin
    git tag -a "v${VERSION}" origin/master -m "Release ${VERSION}"
    git push origin "v${VERSION}"

### 7. Create the GitHub release

    gh release create "v${VERSION}" --verify-tag --title "${VERSION}" --notes-file /tmp/sharkmer-release-notes.md

### 8. Advance the release maintenance branch

For v3.2, fast-forward the existing `v3` branch to the release tag:

    git switch v3
    git merge --ff-only origin/v3
    git merge --ff-only "v${VERSION}"
    git push origin v3

If there is no local `v3`, first use `git switch --track origin/v3`. Create a
new `vN` branch from the tag only when releasing a new major version that has
no maintenance branch. A failed fast-forward requires investigating the branch
history; do not recreate or force-push the existing branch.

### 9. Finalize the bioconda recipe

Download the exact URL used by `meta.yaml`. Hash the saved archive only after
the HTTP request succeeds:

    curl --fail --location --output "/tmp/sharkmer-${VERSION}.tar.gz" "https://github.com/caseywdunn/sharkmer/archive/refs/tags/v${VERSION}.tar.gz" &&
        sha256sum "/tmp/sharkmer-${VERSION}.tar.gz"

Replace `PLACEHOLDER` in `meta.yaml` with this checksum on the post-release
branch below. Test and submit the recipe through the Bioconda workflow below;
record the package pull request and build outcome with the release evidence.

### 10. Resume development and verify ancestry

Start from `origin/master` so the next development branch carries the actual
release merge commit:

    git fetch origin
    git switch -c "post-release-${VERSION}" origin/master

For v3.2, set the package version in `Cargo.toml` and `Cargo.lock` to
`4.0.0-dev`, add an empty `[Unreleased]` changelog section, and prune the released
section from `dev_docs/PLAN.md` while preserving its evidence and deferred work.
Keep the recipe at `3.2.0` and insert the checksum from step 9. Run
`scripts/release-check`, commit, and return these changes to `dev` through a PR:

    git push -u origin "post-release-${VERSION}"
    gh pr create --base dev --head "post-release-${VERSION}" --title "Start development after ${VERSION}" --body "Preserve release ancestry, start 4.0.0-dev, and finalize the released source checksum."
    gh pr checks "post-release-${VERSION}" --watch
    gh pr merge "post-release-${VERSION}" --merge
    git fetch origin
    git merge-base --is-ancestor "v${VERSION}" origin/dev
    gh run list --workflow tag-check.yml --branch "v${VERSION}"

The ancestry check must succeed. Verify that the tag workflow passed; it
requires the tag, Cargo versions, recipe version, and dated changelog to agree.

## Regression benchmarks

A benchmark suite in `benchmarks/` runs sharkmer against 13 real-world SRA
datasets across 6 of the 9 built-in primer panels. Run after each development
phase or before a release to check for regressions.

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

    # Run the included offline counting oracle
    python benchmarks/run_benchmark.py --config benchmarks/known_truth.yaml -k 19 --scope counting-only --no-blast

Per-panel YAML and Markdown results and a commit-labelled combined summary are
written to `benchmarks/benchmark_results/`. Compare the structured per-product
records and performance fields between runs; the repository does not currently
provide an automated cross-run comparator.

Dataset roles, offline known-truth commands, and the registration gate for
future held-out data are defined in `benchmarks/DATASETS.md`. Results record
the executable fingerprint and workspace source observation, effective
command, input and reference checksums, cache protocol, stage timings,
throughput, table capacity, allocator peak, peak RSS when available, and final
run-directory size. Unavailable metrics remain null with an availability
statement.

Sample data is cached in `benchmarks/data/` and downloaded from ENA on first
run. The 1M setting is a requested FASTQ-record cap, not a promise that every
input reaches it; four calibration inputs exhaust earlier. The download streams
and truncates early, so it does not fetch full runs. Configuration (samples,
panels, read counts) is in `benchmarks/benchmark.yaml`.

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

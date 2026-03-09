# Contributing

## Development

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

All work is done on issue branches that correspond to entries in the
[issue tracker](https://github.com/caseywdunn/sharkmer/issues). Issue branches
are merged into `dev` only when passing. `dev` is merged into `master` only for
releases.

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
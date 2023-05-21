# sharkmer

A kmer analysis tool written in rust

## Development

After cloning the repo, gunzip the `data` in the data dir:

    gunzip -c SRR5324768_pass_1.fastq.gz > SRR5324768_pass_1.fastq

Some common tasks in development:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release

### Test data

This repository includes a test [dataset from Thermus thermophilus](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5324768&display=metadata) to develop against.
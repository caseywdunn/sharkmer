# sharkmer

`sharkmer` is a kmer counter designed from the ground up for rarefaction analyses. It counts kmers
on subsets of the data, and then builds incremental histograms from these counts. This allows you
to see how the results change as data are added. This builds more insight from the same data, and 
helps with practical questions such as figuring out how much coverage you need to get good 
genome size estimates in your organism and deciding whether to collect more data.

There are two components to sharkmer:
- The `sharkmer` executable, written in rust, that inputs reads, counts kmers, and outputs histograms.
- The `sharkmer_viewer.py` python script for viewing and analyzing the histogram

Here is an overview of how `sharkmer` works:
1. fastq data are ingested one read at a time and recoded as 8 bit integers, with 2 bits per base. Reads 
   are broken into subreads at any instances of `N`, since 2 bit encoding only covers the 4 unambiguous 
   bases and kmers can't span them anyway. This encoding can only store bases in multiples of 4, so the 
   last 0-3 bases in each subread are discarded. This encoding strategy greatly improves memory efficiency 
   and speed at the cost of discarding a small fraction of bases (about 1.3% for 150bp reads). The discarded
   bases are at the ends of reads, which tend to be lower quality anyway.
2. The order of the subreads is shuffled.
2. A bloom filter containing all kmers that occur more than once is 
   constructed.
3. The subreads are broken into `n` chunks of subreads. Within each chunk, kmers in the bloom filter 
   (ie, kmers that were observed more than once) are counted in a hashmap. This excludes singleton reads,
   which are abundant and almost all sequencing errors, from counting and greatly reduces the size of these 
   hashmaps. This improves 
   RAM and computational performance. The counting of kmers in parallelized across chunks, allowing multiple
   threads to be used.
4. The hashmaps for the `n` chunks are summed one by one, and a histogram is generated after each chunk of 
   counts is added in. This produces `n` histograms, each summarizing more reads than the last.


A few notes:
- The read data must be uncompressed before analysis. No `.fastq.gz` files, just `.fastq`.
- All the read data are stored in RAM in a compressed integer format. This is a tradeoff that improves speed at the cost of requiring more memory. Every 4 gigabases of sequence reads will need about 1 GB of RAM to store. This means that a 2 gigabase genome with 60x coverage (120 gigabases of reads) will need 30GB of RAM just to store the reads. Additional memory is needed for the hashmaps.

## Installation

Download and install [Rust and Cargo](https://www.rust-lang.org/tools/install).

Then, in this repo build the binary:

    cargo build --release

The executable is then at `target/release/sharkmer`. Move it somewhere into your path.

## Usage

To get full usage information, run

    sharkmer --help

An example analysis would look like this:

    sharkmer -k 21 -n 10 --histo-max 10000 -o Agalma-elegans agalma_1_R1.fastq agalma_1_R2.fastq agalma_2_R1.fastq agalma_2_R2.fastq

The incremental histogram files in this case would be:

    Agalma-elegans.histo # All the incremental histograms, each in their own column. Suitable for analysis with `sharkmer_viewer.py`.

    Agalma-elegans_final.histo # Just the final histogram with all data. Suitable for analysis with genomescope and other tools.


Then to explore the results:
    python sharkmer_viewer.py Agalma-elegans.histo

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
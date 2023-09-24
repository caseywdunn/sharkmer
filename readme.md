# sharkmer - a kmer counter and analysis tool

Functionalities of sharkmer include:

- Incremental kmer counting. This allows you to run kmers on incrementally larger subsets of your data. Applications include assessing the robustness of genome size estimates to sequencing depth.
- in silico PCR. This allows you to supply primer pairs and a fastq file, and get a fasta file of the amplicons that would be produced by PCR. This is useful for assembling and isolating particular genes from raw genome skimming data.  

There are two components to sharkmer:

- The `sharkmer` executable, written in rust, that inputs reads, counts kmers, and outputs histograms.
- The optional `sharkmer_viewer.py` python script for viewing and analyzing histograms

`sharkmer` counts kmers
on subsets of the data, and then builds incremental histograms from these counts. This is more efficient than counting kmers on all the data at once, and allows you
to see how the results change as data are added. This builds more insight from the same data, and 
helps with practical questions such as figuring out how much coverage you need to get good 
genome size estimates in your organism and deciding whether to collect more data.

Here is an overview of how kmer counting works in `sharkmer`:

1. fastq data are ingested one read at a time and recoded as 8 bit integers, with 2 bits per base. Reads 
   are broken into subreads at any instances of `N`, since 2 bit encoding only covers the 4 unambiguous 
   bases and kmers can't span them anyway. This encoding can only store bases in multiples of 4, so the 
   last 0-3 bases in each subread are discarded. This encoding strategy greatly improves memory efficiency 
   and speed at the cost of discarding a small fraction of bases (about 1.3% for 150bp reads). The discarded
   bases are at the ends of reads, which tend to be lower quality anyway.
2. The order of the subreads is shuffled.
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
- All the read data are stored in RAM in a compressed integer format. This is a tradeoff that improves speed at the cost of requiring more memory. Every 4 gigabases of sequence reads will need about 1 GB of RAM to store. This means that a 2 gigabase genome with 60x coverage (120 gigabases of reads) will need 30GB of RAM just to store the reads. Additional memory is needed for the hashmaps. For some analyses, such as **in silico PCR**, you can use a small subset of reads and easily run analyses on a laptop.

## Installation

This repository includes sharkmer, which is written in rust, and some helper programs written in python. The python components are only needed for some followup analyses.

### Rust components

#### From source

First, [install the rust build tools](https://www.rust-lang.org/tools/install).

Then clone the respository and build sharkmer (note that the sharkmer rust code is in the `sharkmer/sharkmer` folder, not the top level `sharkmer` folder):

    git clone https://github.com/caseywdunn/sharkmer.git
    cd sharkmer/sharkmer
    cargo build --release

The executable will be in `sharkmer/target/release/sharkmer`. Move it to a location in your path.

### Python components

You can create a conda environment with all needed python components as follows:

    conda create -n shark -c conda-forge -c bioconda python==3.10 ffmpeg genomescope2
    conda activate shark
    cd sharkmer_viewer/
    pip install .

## Test data

### **Thermus thermophilus**

This repository includes a test [dataset from Thermus thermophilus](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5324768&display=metadata) for simple tests and to develop against.

After cloning the repo, gunzip the `data` in the data dir:

    cd sharkmer/data/ # Note that this is the sharkmer folder within the sharkmer repository
    gunzip -c SRR5324768_pass_1.fastq.gz > SRR5324768_pass_1.fastq

### **Cordagalma ordinatum**

To get a sense of the tool it is best to grab larger datasetsets. The examples below will use a dataset from the siphonophore **Cordagalma ordinatum**. This dataset is available from the NCBI SRA as [SRR23143278](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR23143278), and is from the manuscript:

> Ahuja, N., Cao, X., Schultz, D. T., Picciani, N., Lord, A., Shao, S., Burdick, D. R., Haddock, S. H. D., Li, Y., & Dunn, C. W. (2023). Giants among Cnidaria: large nuclear genomes and rearranged mitochondrial genomes in siphonophores. bioRxiv. https://doi.org/10.1101/2023.05.12.540511

To retrieve these data, first download and install the [sra toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit). Details on downloading data are provided in the examples below.



Then run:

    fastq-dump SRR23143278

If you want just a million reads to try out sPCR, then you can get them without first downloading the whole dataset:

    fastq-dump -X 1000000 SRR23143278

## Usage

To get full usage information, run

    sharkmer --help

### Incremental kmer counting

Incremental kmer counting for genome size estimation takes a lot of data (about 50x coverage of the genome), and a large amount of RAM. So this example isn't practical on most laptops, given their disk and RAM limitations, and will require a workstation or cluster.

First download the data as follows (this is the full dataset, so it will be quite large):

    fastq-dump SRR23143278

An example analysis would look like this:

    sharkmer -k 21 -n 10 --histo-max 10000 -o Cordagalma-ordinatum SRR23143278.fastq

The incremental histogram files in this case would be:

    Cordagalma-ordinatum.histo # All the incremental histograms, each in their own column. Suitable for analysis with `sharkmer_viewer.py`.

    Cordagalma-ordinatum.final.histo # Just the final histogram with all data. Suitable for analysis with genomescope and other tools.

Then to explore the results:

    sharkmer_viewer Cordagalma-ordinatum.histo

The final histogram on all the data is also written to its own file, and you can view that with, for example, [GenomeScope2](https://github.com/tbenavi1/genomescope2.0):

    genomescope2 -i Cordagalma-ordinatum.final.histo -o Cordagalma-ordinatum -k 21

The included `genomemovie.sh` script will generate a movie of the incremental GenomeScope2 histograms. For example, to create a movie of the `Cordagalma` test dataset:

    conda activate shark
    bash genomescopemovie.sh Cordagalma-ordinatum.histo Cordagalma-ordinatum.output

### **in silico** PCR (sPCR)

It is often very useful to pull small genome regions out of genome skimming data, for example to blast a commonly sequenced gene to verify that the sample you sequenced is the species you expected. This common task is surprisingly challenging in practice, though. You can map reads to known sequences and then collapse them into a sequence prediction, but this does not always work well across species and can miss variable regions. You can assemble all the reads and then pull out the region of interest, but this is computationally expensive and often the region of interest is not assembled well given how shallow skimming data often are.

**in silico** PCR (sPCR) is a new alternative approach. You specify file with raw reads and one or more primer pairs, and sharkmer outputs a fasta file with the sequence of the region that would be amplified by PCR on the genome the reads are derived from. There are multiple advantages to this approach:

- sPCR directly leverages the decades of work that have been done to optimize PCR primers that work well across species and span informative gene regions. These primers tend to bind conserved regions that flank variable informative regions and have minimal off-target binding.
- Because it is primer based, you can use it to obtain the exact same gene regions (co1, 16s, 18s, 28s, etc...) that have been PCR amplified for decades and still remain the most broadly sampled across species in public databases.
- sPCR doesn't take much data. You can use small datasets, or analyze small (eg one million read) subsets of your data.
- sPCR is fast and has minimum computational requirements. It can be run on a laptop in a couple minutes on a million reads.
- sPCR requires a single tool (sharkmer), not complex workflows with multiple tools.

sPCR is useful when you want specific genes from skimming datasets you have collected for other purposes. But in some cases it may be more cost effective and easier to skim and apply sPCR than to use traditional PCR. With sPCR you sequence once and then pull out as many gene regions as you want, as opposed to PCR where you amplify and sequence each region separately. There is little additional computational cost for each added primer pair, since most of the work is counting kmers and this is done once for all primer pairs. So the cost of sPCR is fixed and does not depend on the number of genes considered.

#### **in silico** PCR example

sPCR of nuclear ribosomal RNA genes (eg animal 16s, 18s) and mitochondrial genes does not take much coverage, given the relatively high copy number of these genes in skimming data. For Illumina raw reads, 0.25x coverage of the genomes is sufficient. The **Cordagalma ordinatum** is 700Mb, so 1 million 150bp reads, a total of 150Mb of data, is sufficient.

Download one million reads of the **Cordagalma ordinatum** dataset:

    fastq-dump -X 1000000 SRR23143278

Then run sPCR on the downloaded reads by specifying primer pairs with the `--pcr` argument:

    sharkmer \
      -k 31 -n 100 -t 4 -o Cordagalma_CWD6 --coverage 3 \
      --pcr "GACTGTTTACCAAAAACATA_AATTCAACATCGAGG_1000_16s" \
      --pcr "TCATAAAGATATTGG_ATGCCCGAAAAACCA_2000_co1" \
      --pcr "AACCTGGTTGATCCTGCCAGT_TGATCCTTCTGCAGGTTCACCTAC_2500_18s" \
      --pcr "CCYYAGTAACGGCGAGT_SWACAGATGGTAGCTTCG_4000_28s"  \
      SRR23143278.fastq

The `--pcr` argument passes a string with the format `forward_reverse_max-length_name`.

This analysis will generate one fasta file for each primer pair. These fasta files are named with the argument passed to `--pcr`. If no product was found, the fasta file is not present. The fasta file can contain more than one sequence.

### Reading compressed data

`sharkmer` does not read compressed data directly, but it can read uncompressed data from `stdin`.
So, for example, you could `gunzip` files and pipe them to `sharkmer`:

    zcat agalma_*.fastq.gz | sharkmer -k 21 -n 10 --histo-max 10000 -o Agalma-elegans

Decompressing files takes quite a bit of compute. Handling decompression outside of `sharkmer` allows you to 
use whichever approach you prefer on your system, for example parallel tools such as `pigz`. 

## Development

### Rust 

Some common tasks in development:

    cargo test
    cargo clippy
    cargo clippy --fix
    cargo fmt
    cargo build # Debug
    cargo build --release

### Docker

Example:
    cd docker
    docker build -t shark .
    docker run -itv /myhome/repos/sharkmer:/sharkmer shark /bin/bash

# sharkmer - a kmer analysis tool

Functionalities of sharkmer include:

- Incremental kmer counting. This allows you to run kmers on incrementally larger subsets of your data. Applications include assessing the robustness of genome size estimates to sequencing depth.
- *in silico* PCR (sPCR). You supply a fastq file of whole genome shotgun reads and sequences for primer pairs, and get a fasta file of the amplicons that would be produced by PCR on the genome. This is useful for assembling and isolating particular genes from raw genome skimming data.  

There are two components to sharkmer:

- The `sharkmer` executable, written in rust, that inputs reads, counts kmers, and outputs histograms.
- The optional `sharkmer_viewer.py` python script for viewing and analyzing histograms from incremental kmer counting

`sharkmer` counts kmers
on subsets of the data, and then builds incremental histograms from these counts. This is more efficient than counting kmers on all the data at once, and allows you
to see how the results change as data are added. This builds more insight from the same data. It also 
helps with practical questions, such as figuring out how much coverage you need to get good 
genome size estimates in your organism and deciding whether to collect more data.

Here is an overview of how kmer counting works in `sharkmer`:

1. fastq data are ingested one read at a time and recoded as 8 bit integers, with 2 bits per base. Reads 
   are broken into subreads at any instances of `N`, since 2 bit encoding only covers the 4 unambiguous 
   bases and kmers can't span N anyway. The subreads are distributed across `n` chunks of subreads as thew are read. 
2. Within each chunk, kmers are counted in a hashmap. 
   The counting of kmers in parallelized across chunks, allowing multiple threads to be used.
3. The hashmaps for the `n` chunks are summed one by one, and a histogram is generated after each chunk of 
   counts is added in. This produces `n` histograms, each summarizing more reads than the last.

A few notes:

- The read data must be uncompressed before analysis. No `.fastq.gz` files, just `.fastq`.

## Installation

This repository includes sharkmer, which is written in rust, and some helper programs written in python. The python components are optional tools for some followup analyses.

### Rust components

#### From source

First, [install the rust build tools](https://www.rust-lang.org/tools/install).

Then clone this repository and build sharkmer (note that the sharkmer rust code is in the `sharkmer/sharkmer` folder, not the top level `sharkmer` folder):

    git clone https://github.com/caseywdunn/sharkmer.git
    cd sharkmer/sharkmer
    cargo build --release

The executable will be in `sharkmer/target/release/sharkmer`. Move it to a location in your path.

### Python components

You can create a conda environment with all optional python components as follows:

    conda create -n shark -c conda-forge -c bioconda python==3.10 ffmpeg genomescope2
    conda activate shark
    cd sharkmer_viewer/
    pip install .

## Test data

### *Thermus thermophilus*

This repository includes a test [dataset from Thermus thermophilus](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5324768&display=metadata) for simple tests and to develop against.

After cloning the repo, gunzip the `data` in the data dir:

    cd sharkmer/data/ # Note that this is the sharkmer folder within the sharkmer repository
    gunzip -c SRR5324768_pass_1.fastq.gz > SRR5324768_pass_1.fastq

### Additional datasets

Additional datasets are available from the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra). Install the [sra toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to download these datasets as described below.


## Usage

To get full usage information, run

    sharkmer --help

### Incremental kmer counting

Incremental kmer counting for genome size estimation takes a lot of data (about 50x coverage of the genome), and a large amount of RAM. So this example isn't practical on most laptops, given their disk and RAM limitations, and will require a workstation or cluster.

We will need a relatively large dataset. Download this *Cordagalma ordinatum* dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/SRX10340700) using the [sra toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit):

    fastq-dump --split-files SRR23143278

Now perform incremental kmer counting on the downloaded reads:

    sharkmer -o output/ -s Cordagalma-ordinatum data/SRR23143278_1.fastq data/SRR23143278_2.fastq

Notice that you can specify multiple fastq files, in this case the R1 and R2 reads. The `-o` argument specifies the output directory. The `-s` argument specifies a prefix for the output files.

The incremental histogram files in this case would be:

    output/Cordagalma-ordinatum.histo # All the incremental histograms, each in their own column. Suitable for analysis with `sharkmer_viewer.py`.

    output/Cordagalma-ordinatum.final.histo # Just the final histogram with all data. Suitable for analysis with genomescope and other tools.

Then to explore the results with the optional python tools:

    conda activate shark
    sharkmer_viewer Cordagalma-ordinatum.histo

The final histogram on all the data is in a standard format and can be examined with a variety of tools, for example, [GenomeScope2](https://github.com/tbenavi1/genomescope2.0):

    genomescope2 -i Cordagalma-ordinatum.final.histo -s Cordagalma-ordinatum -k 21

The included `genomemovie.sh` script will generate a movie of the incremental GenomeScope2 histograms. For example, to create a movie of the `Cordagalma` test dataset:

    conda activate shark
    bash genomescopemovie.sh Cordagalma-ordinatum.histo Cordagalma-ordinatum.output

### *in silico* PCR (sPCR)

Investigators often want to pull small genome regions out of genome skimming data, for example to blast a commonly sequenced gene to verify that the sample is the expected species. This task is surprisingly challenging in practice, though. You can map reads to known sequences and then collapse them into a sequence prediction, but this does not always work well across species and can miss variable regions. You can assemble all the reads and then pull out the region of interest, but this is computationally expensive and often the region of interest is not assembled well given how shallow skimming data often are.

*in silico* PCR (sPCR) is a new alternative approach. You specify a file with raw reads and one or more primer pairs, and sharkmer outputs a fasta file with the sequence of the region that would be amplified by PCR on the genome the reads are derived from. There are multiple advantages to this approach:

- sPCR directly leverages the decades of work that have been done to optimize PCR primers that work well across species and span informative gene regions. These primers tend to bind conserved regions that flank variable informative regions and have minimal off-target binding.
- Because sPCR is primer based, you can use it to obtain the exact same gene regions (co1, 16s, 18s, 28s, etc...) that have been PCR amplified for decades and still remain the most broadly sampled across species in public databases.
- sPCR doesn't take much data. You can use small datasets, or analyze small (eg one million read) subsets of your data.
- sPCR is fast and has minimum computational requirements. It can be run on a laptop in a couple minutes on a million reads.
- sPCR requires a single tool (sharkmer), not complex workflows with multiple tools.

sPCR is useful when you want specific genes from skimming datasets you have collected for other purposes. But in some cases it may be more cost effective and easier to skim and apply sPCR than to use traditional PCR. With sPCR you sequence once and then pull out as many gene regions as you want, as opposed to PCR where you amplify and sequence each region separately. There is little additional computational cost for each added primer pair, since most of the work is counting kmers and this is done once for all primer pairs. So the cost of sPCR is fixed and does not depend on the number of genes considered.

#### *in silico* PCR example

sPCR of nuclear ribosomal RNA genes (eg animal 16s, 18s) and mitochondrial genes does not take much coverage, given the relatively high copy number of these genes in skimming data. For Illumina raw reads, 0.25x coverage of the genomes is sufficient.

We will download a smaller dataset from the coral *Stenogorgia casta*:

    fasterq-dump SRR26955578

Run sPCR on the downloaded reads by specifying that we want to run a panel of cnidarian primers:

    sharkmer --max-reads 1000000 -s Stenogorgia_casta -o output/ --pcr cnidaria SRR26955578_1.fastq SRR26955578_2.fastq


This is equivalent to specifying the primer pairs manually, also with the `--pcr` argument:

    sharkmer \
      --max-reads 1000000 \
      -s Cordagalma_CWD6 -o output/ \
      --pcr "GRCTGTTTACCAAAAACATA,AATTCAACATMGAGG,700,16s,min_length=500" \
      --pcr "TCATAARGATATHGG,RTGNCCAAAAAACCA,800,co1,min_length=600" \
      --pcr "AACCTGGTTGATCCTGCCAGT,TGATCCTTCTGCAGGTTCACCTAC,2000,18s,min_length=1600" \
      --pcr "CCYYAGTAACGGCGAGT,SWACAGATGGTAGCTTCG,3500,28s,min_length=2900"  \
      --pcr "TACACACCGCCCGTCGCTACTA,ACTCGCCGTTACTRRGG,1000,ITSfull,min_length=600" \
      SRR26955578_1.fastq SRR26955578_2.fastq



The `--pcr` argument passes a string with the format `forward,reverse,max-length,gene-name`. Note that commas delimit fields. `max-length` should be greater than the expect PCR product size. It indicates the furthest distance from the forward primer that sharkmer should search for a reverse primer.

The `--max-reads 1000000` arguments indicates that the first million reads should be used. THis is plenty for nuclear rRNA sequences 18s, 28s, and ITS, since it occurs in many copies in the genome, and mitochondrial sequences 16s and co1. Single copy nuclear genes would require more data.

This analysis will generate one fasta file for each primer pair. These fasta files are named with the argument passed to `--pcr`. If no product was found, the fasta file is not present. The fasta file can contain more than one sequence.

There are a limited set of preconfigured primer panels available at this time. You can see where they are hard coded [here](https://github.com/caseywdunn/sharkmer/blob/dev/sharkmer/src/pcr/preconfigured.rs). If you have other primers that you would like to have added to the tool, or suggestings for optimizing the primers that are already there, please open an issue in the [issue tracker](https://github.com/caseywdunn/sharkmer/issues).

### Reading compressed data

`sharkmer` does not read compressed data directly, but it can read uncompressed data from `stdin`.
So, for example, you could `gunzip` files and pipe them to `sharkmer`:

    zcat agalma_*.fastq.gz | sharkmer -k 21 -n 10 --histo-max 10000 -s Agalma-elegans

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

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
3. The hashmaps for the `n` chunks are summed one by one, and a histogram is generated after each chunk of
   counts is added in. This produces `n` histograms, each summarizing more reads than the last.

A few notes:

- The read data must be uncompressed before analysis. No `.fastq.gz` files, just `.fastq`. But you can uncompress data with an external tool and pipe it to `sharkmer`, see the [Reading compressed data](#reading-compressed-data) section below.

## Citing

We are working on a paper that describes sharkmer. In the meantime, please cite the two papers below. I wrote the tool for these projects, though we don't describe it in detail there.

For incremental kmer counting:

> N Ahuja, X Cao, DT Schultz, N Picciani, A Lord, S Shao, K Jia, DR Burdick, SH D Haddock, Y Li, CW Dunn (2024) Giants among Cnidaria: Large Nuclear Genomes and Rearranged Mitochondrial Genomes in Siphonophores. Genome Biology and Evolution, 16(3). [doi:10.1093/gbe/evae048](https://doi.org/10.1093/gbe/evae048).

For *in silico* PCR:

> Church et al. (2024) Global genomics of the man-o’-war (Physalia) reveals biodiversity at the ocean surface. [doi:10.1101/2024.07.10.602499](https://doi.org/10.1101/2024.07.10.602499)

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

These are optional. But if you do want to install them you can create a conda environment with all optional python components as follows:

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

Additional datasets are available from the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra). Install the [sra toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to download these datasets with `fasterq-dump` as described below.

## Usage

To get full usage information, run

    sharkmer --help

### Incremental kmer counting

Incremental kmer counting for genome size estimation takes a lot of data (about 50x coverage of the genome), and a large amount of RAM. So this example isn't practical on most laptops, given their disk and RAM limitations, and will require a workstation or cluster.

We will need a relatively large dataset. Download this *Cordagalma ordinatum* dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/SRX10340700) using the [sra toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit):

    cd data
    prefetch SRR23143278
    fasterq-dump --split-files SRR23143278
    cd ..

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

*in silico* PCR (sPCR) is a new alternative approach. You specify a file with raw reads, and specify one or more primer pairs. sharkmer outputs a fasta file with the sequence of the region that would be amplified by PCR on the genome the reads are derived from. There are multiple advantages to this approach:

- sPCR directly leverages the decades of work that have been done to optimize PCR primers that work well across species and span informative gene regions. These primers tend to bind conserved regions that flank variable informative regions and have minimal off-target binding.
- Because sPCR is primer based, you can use it to obtain the exact same gene regions (co1, 16s, 18s, 28s, etc...) that have been PCR amplified for decades and still remain the most broadly sampled across species in public databases.
- sPCR often doesn't take much data. For small genomes or high copy sequences (such as rRNA), you can analyze small (eg one million read) subsets of your data.
- sPCR is fast and has minimum computational requirements. It can be run on a laptop in a couple minutes on a million reads.
- sPCR requires a single tool (sharkmer), not complex workflows with multiple tools.

sPCR is useful when you want specific genes from skimming datasets you have collected for other purposes. But in some cases it may be more cost effective and easier to skim and apply sPCR than to use traditional PCR. With sPCR you sequence once and then pull out as many gene regions as you want, as opposed to PCR where you amplify and sequence each region separately. There is little additional computational cost for each added primer pair, since most of the work is counting kmers and this is done once for all primer pairs. So the cost of sPCR is fixed and does not depend on the number of genes considered.

#### *in silico* PCR example

sPCR of nuclear ribosomal RNA genes (eg animal 16s, 18s) and mitochondrial genes does not take much coverage, given the relatively high copy number of these genes in skimming data. For Illumina raw reads, 0.25x coverage of the genomes is sufficient.

We will download a smaller dataset from the coral *Stenogorgia casta*:

    cd data
    prefetch SRR26955578
    fasterq-dump --split-files SRR26955578
    cd ..

Run sPCR on the downloaded reads by specifying that we want to run a panel of cnidarian primers:

    sharkmer --max-reads 1000000 -s Stenogorgia_casta -o output/ --pcr cnidaria data/SRR26955578_1.fastq data/SRR26955578_2.fastq

This is equivalent to specifying the primer pairs manually, also with the `--pcr` argument:

    sharkmer \
      --max_reads 1000000 \
      -s Cordagalma_CWD6 -o output/ \
      --pcr "forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16s,min-length=500" \
      --pcr "forward=TCATAARGATATHGG,reverse=RTGNCCAAAAAACCA,max-length=800,name=co1,min-length=600" \
      --pcr "forward=AACCTGGTTGATCCTGCCAGT,reverse=TGATCCTTCTGCAGGTTCACCTAC,max-length=2000,name=18s,min-length=1600" \
      --pcr "forward=CCYYAGTAACGGCGAGT,reverse=SWACAGATGGTAGCTTCG,max-length=3500,name=28s,min-length=2900"  \
      --pcr "forward=TACACACCGCCCGTCGCTACTA,reverse=ACTCGCCGTTACTRRGG,max-length=1000,name=ITSfull,min-length=600" \
      data/SRR26955578_1.fastq data/SRR26955578_2.fastq

The `--pcr` argument passes a string with the format `key1=value1,key2=value2,key3=value3,...`, where the required keys are `forward`, `reverse`, `name`, and `max-length`. Note that commas delimit fields. `max_-_length` should be greater than the expect PCR product size. It indicates the furthest distance from the forward primer that sharkmer should search for a reverse primer.

The `--max_reads 1000000` arguments indicates that the first million reads should be used. This is plenty for nuclear rRNA sequences 18s, 28s, and ITS, since it occurs in many copies in the genome, and mitochondrial sequences 16s and co1. Single copy nuclear genes would require more data.

This analysis will generate one fasta file for each primer pair. These fasta files are named with the argument passed to `--pcr`. If no product was found, the fasta file is not present. The fasta file can contain more than one sequence.

There are a limited set of preconfigured primer panels available at this time. You can see where they are hard coded [here](https://github.com/caseywdunn/sharkmer/blob/dev/sharkmer/src/pcr/preconfigured.rs). If you have other primers that you would like to have added to the tool, or suggestings for optimizing the primers that are already there, please open an issue in the [issue tracker](https://github.com/caseywdunn/sharkmer/issues).

#### Optimizing *in silico* PCR (sPCR)

There are a few different strategies to take if you are not getting a sPCR product, or it is working inconsistently. If you get things working, please let me know in in the [issue tracker](https://github.com/caseywdunn/sharkmer/issues) so I can improve the default primer sets and help other users. If you hit a wall, please also met me know in the issue tracker.

The things you should try first are:

- Increase the `--pcr` parameter `max-length`. It may be that the amplified region is longer than expected. And if there are introns, the amlified product can be much longer than the coding sequence.

- Optimize your primer sequences. Make some multiple sequence alignments of the desired sequence region from several closely related species, and refine the primer sequences to be more specific to the target region. You can use [degenerate nucleotide symbols](https://en.wikipedia.org/wiki/Nucleic_acid_notation), such as R for A or G, a variable sites in the site where you would like the sequence to bind. Remember that, just as in real PCR, the reverse primer should be reverse complemented.

- Adjust the number of reads. If you are not getting a product, try increasing the number of reads. You can do this with the `--max_reads` argument. Likewise, if you are getting products and want to speed things up, or you are getting many products for a gene, reduce the number of reads.

- Adjust the `--pcr` parameter `trim`. The default is 15. This is the max number of bases to keep at the 3' end of each primer. Primers used for real PCR tend to be longer than what is required for successful specific sPCR. This is becuase they are lengthened to adjust melting temperature. If you are getting no product, try reducing this value. If you are getting too many spurious products, try increasing this value.

If these do not work, then you can try adjusting other parameters.

- Specify a reasonable `--pcr` parameter `min-length`. This value defaults to 0, but raising it can get rid of small spurious products.

- Adjust the `--pcr` parameter `coverage`. This is the depth of kmer coverage required to extend a sPCR product. It defaults to 3, and must be aat least 2 to avoid once off errors. You can try raising it to 4 or 5 to get only the best supported products, but this requires more data. For example, `--pcr "forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16s,min-length=500,coverage=5"`

- Adjust the `--pcr` parameter `mismatches`. This defaults to 2. You can try raising it to 3 or 4 if you aren't getting the desired product, but this may increase the number of spurious paths that need to be traversed and bog down the run.

Keep in mind that there is no way to assemble a sPCR product without coverage along its full length that meets or exceeds the `coverage` parameter. The tool cannot output sequences that are not in the input. If you are trying to amplify a single copy nuclear gene, that means your sequencing depth (average coverage) of the genome will need to be higher than the `coverage` parameter, since there will be fluctuations in coverage along the length of the region. If covereage at each site is independently distributed, then to have a 95% chance of coverage $\geq 3$ at each site in a region of nucleotide length $n$, you would need a sequencing depth of 15x for a 1000bp region, and 14x for a 500bp region. That is on the order of 30 million 150 bp reads for a 300Mb genome. This may place single copy nuclear genes out of reach for some organisms with larger genomes, especially if computer RAM limits the number of reads that can be processed.

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

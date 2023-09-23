# Sharkmer

A collection of kmer counting and analysis tools.

Functionalities include:

- Incremental kmer counting. This allows you to run kmers on incrementally larger subsets of your data. Applications include assessessing the robustness of genome size estimates to sequencing depth.
- in silico PCR. This allows you to supply primer pairs and a fastq file, and get a fasta file of the amplicons that would be produced by PCR. This is useful for assembling and isolating particular genes from raw genome skimming data.  

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

### Optional python components
#### Conda

To build and install the conda package locally:

    conda build .
    conda install miniconda3/envs/devel/conda-bld/osx-arm64/sharkmer_rust-0.1.0-0.tar.bz2 # exact path wil come from above command

## Usage

### Incremental kmer counting

### **in silico** PCR (sPCR)

It is often very useful to pull small genome regions out of genome skimming data, for example to blast a commonly sequenced gene to verify that the sample you sequenced is the species you expected. This common task is surprisingly challenging in practice, though. You can map reads to known sequences and then collapse them into a sequence prediction, but this does not always work well across species and can miss variable regions. You can assemble all the reads and then pull out the region of interest, but this is computationally expensive and often the region of interest is not assembled well given how shallow skimming data often are.

**in silico** PCR (sPCR) is a new alternative approach. You specify file with raw reads and one or more primer pairs, and sharkmer outputs a fasta file with the sequence of the region that would be amplified by PCR on the genome the reads are derived from. There are multiple advantages to this approach:

- sPCR directly leverages the decades of work that have been done to optimize PCR primers that work well across species and span informative gene regions. These primers tend to bind conserved regions that flank variable informative regions and have minimal off-target binding.
- Because it is primer based, you can use it to obtain the exact same gene regions (co1, 16s, 18s, 28s, etc...) that have been PCR amplified for decades and still remain the most broadly sampled across species in public databases.
- sPCR doesn't take much data. You can use small datasets, or analyze small (eg one million read) subsets of your data.
- sPCR is fast and has minimum computational requirements. It can be run on a laptop in a couple minutes on a million reads.
- sPCR requires a single tool (sharkmer), not complex workflows with multiple tools.

sPCR is useful when you want specific genes from skimming datasets you have collected for other purposes. But in some cases it may be more cost effective and easier to skim and apply sPCR than to use traditional PCR. With sPCR you sequence once and then pull out as many gene regions as you want, as opposed to PCR where you amplify and sequence each region separately. There is little additional computational cost for each added primer pair, since most of the work is counting kmers and this is done once for all primer pairs. So the cost of sPCR is fixed and does not depend on the number of genes considered.

## Development
### Debugging

To debug in vs code:

- Go to File > add folder to workspace
- Select `sharkmer/sharkmer` folder

This opens just the rust code, and debugging works fine. 
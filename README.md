# sharkmer - a kmer analysis tool

[![CI](https://github.com/caseywdunn/sharkmer/actions/workflows/ci.yml/badge.svg)](https://github.com/caseywdunn/sharkmer/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19020708.svg)](https://doi.org/10.5281/zenodo.19020708)
[![Bioconda version](https://img.shields.io/conda/vn/bioconda/sharkmer.svg?label=Bioconda)](https://anaconda.org/bioconda/sharkmer)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/sharkmer.svg?label=downloads)](https://anaconda.org/bioconda/sharkmer)

`sharkmer` is a kmer counter and seeded *de Bruijn* graph assembler. `sharkmer` features include:

- *in silico* PCR (sPCR). You supply a fastq file of whole genome shotgun reads and sequences for primer pairs, and get a fasta file of the amplicons that would be produced by PCR on the genome. This is useful for assembling and isolating particular genes from raw genome skimming data.  
- Incremental kmer counting. This allows you to run kmers on incrementally larger subsets of your data. Applications include assessing the robustness of genome size estimates to sequencing depth.

## Citing

If you use sharkmer in published work, please cite:

For *in silico* PCR:

> Dunn CW, Church SH (2026) Sharkmer: repurposing PCR primers for targeted genome assembly using in silico PCR. Bioinformatics, btag163. [doi:10.1093/bioinformatics/btag163](https://doi.org/10.1093/bioinformatics/btag163)

For incremental kmer counting:

> N Ahuja, X Cao, DT Schultz, N Picciani, A Lord, S Shao, K Jia, DR Burdick, SH D Haddock, Y Li, CW Dunn (2024) Giants among Cnidaria: Large Nuclear Genomes and Rearranged Mitochondrial Genomes in Siphonophores. Genome Biology and Evolution, 16(3). [doi:10.1093/gbe/evae048](https://doi.org/10.1093/gbe/evae048).

## Installation

### With conda

This is the simplest way to install sharkmer.

If you don't have it already, [install miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install).

To install sharkmer into an existing active conda environment:

    conda install -c bioconda -c conda-forge sharkmer

To create a new environment (here called `sharkmer_env`, but you can name it anything you like) and install sharkmer in it:

    conda create -n sharkmer_env -c bioconda -c conda-forge sharkmer

### From the repository

If you don't have them already, [install the rust build tools](https://www.rust-lang.org/tools/install).

Then clone this repository and build sharkmer:

    git clone https://github.com/caseywdunn/sharkmer.git   # Or download and expand the zip file
    cd sharkmer
    cargo build --release

You can then install the built binary into your Cargo bin directory (usually `~/.cargo/bin/`), and later uninstall it:

    cargo install --path .
    cargo uninstall sharkmer

Alternatively, you can copy the compiled executable from `target/release/sharkmer` to a directory already in your `PATH`, or add `target/release/` to your `PATH`.

## Usage

To get full usage information, run

    sharkmer --help

`sharkmer` reads FASTQ input from files (`.fastq` or `.fastq.gz`) or from stdin. Gzip-compressed files are detected and decompressed automatically.

    sharkmer --max-reads 1000000 -s Agalma-elegans -o output/ --pcr-panel cnidaria agalma_*.fastq.gz

Raw reads, including [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) reads, can be streamed directly from [ENA](https://www.ebi.ac.uk/ena/browser/home) by accession using `--ena`, without requiring any additional tools:

    sharkmer --max-reads 1000000 -s Agalma-elegans -o output/ --pcr-panel cnidaria --ena SRR25099394

Uncompressed FASTQ data can also be piped via stdin. For example, to trim adapters and low quality regions with [fastp](https://github.com/OpenGene/fastp) before counting kmers:

    fastp -i agalma_1.fastq.gz -I agalma_2.fastq.gz --stdout | sharkmer --max-reads 1000000 -s Agalma-elegans -o output/ --pcr-panel cnidaria

## *in silico* PCR (sPCR)

Investigators often want to extract small genome regions from genome skimming data, for example to blast a commonly sequenced gene to verify that the sample is the expected species. This task is surprisingly challenging in practice, though. You can map reads to known sequences and then collapse them into a sequence prediction, but this does not always work well across species and can miss variable regions. You can assemble all the reads and then pull out the region of interest, but this is computationally expensive and often the region of interest is not assembled well.

*in silico* PCR (sPCR) is a new alternative approach. You specify a file with raw reads, and provide one or more primer pairs. `sharkmer` outputs a fasta file with the sequence of the region that would be amplified by PCR on the genome the reads are derived from. There are multiple advantages to this approach:

- sPCR directly leverages the decades of work to optimize PCR primers that work well across species and span informative gene regions. These primers tend to bind conserved regions that flank variable informative regions and have minimal off-target binding.
- Because sPCR is primer based, you can use it to obtain the exact same gene regions (co1, 16s, 18s, 28s, etc...) that have been PCR amplified for decades and still remain the most broadly sampled across species in public databases.
- sPCR often doesn't take much data. For small genomes or high copy regions (such as rRNA or mitochondrial genes), you can analyze small (eg one million read) subsets of your data.
- sPCR is fast and has minimum computational requirements. It can be run on a laptop in a couple minutes on a million reads.
- sPCR requires a single tool (sharkmer), not complex workflows with multiple tools.

sPCR is useful when you want specific genes from skimming datasets you have collected for other purposes. But in some cases it may be more cost effective and easier to skim and apply sPCR than to use traditional PCR. With sPCR you sequence once and then pull out as many gene regions as you want, as opposed to PCR where you amplify and sequence each region separately. There is little additional computational cost for each added primer pair, since most of the work is counting kmers and this is done once for all primer pairs. So the cost of sPCR is fixed and does not depend on the number of genes considered.

### *in silico* PCR example

sPCR of nuclear ribosomal RNA genes (eg animal 28s, 18s, ITS) and mitochondrial genes does not take much sequence data, given the relatively high copy number of these genes. For Illumina raw reads, 0.25x average sequencing depth of the genome is often sufficient.

We will use the coral *Stenogorgia casta* as an example. With `--ena`, sharkmer downloads the SRA reads directly from ENA — no extra tools needed:

    sharkmer --max-reads 1000000 -s Stenogorgia_casta -o output/ --pcr-panel cnidaria --ena SRR26955578

Alternatively, if you have local FASTQ files (gzipped or uncompressed):

    sharkmer --max-reads 1000000 -s Stenogorgia_casta -o output/ --pcr-panel cnidaria data/SRR26955578_1.fastq.gz data/SRR26955578_2.fastq.gz

This is equivalent to specifying the primer pairs manually with `--pcr-primers`:

    sharkmer \
      --max-reads 1000000 \
      -s Stenogorgia_casta -o output/ \
      --pcr-primers "forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16S,min-length=500" \
      --pcr-primers "forward=WAAYCATAAAGATAT,reverse=GGRTGMCCAAAAAACCARA,max-length=800,name=CO1,min-length=600" \
      --pcr-primers "forward=AACCTGGTTGATCCTGCCAGT,reverse=TGATCCTTCTGCAGGTTCACCTAC,max-length=2000,name=18S,min-length=1600" \
      --pcr-primers "forward=CCYYAGTAACGGCGAGT,reverse=SWACAGATGGTAGCTTCG,max-length=3500,name=28S,min-length=2900" \
      --pcr-primers "forward=TACACACCGCCCGTCGCTACTA,reverse=ACTCGCCGTTACTRRGG,max-length=1000,name=ITS,min-length=600" \
      --pcr-primers "forward=GTAGGTGAACCTGCAGAAGGATCA,reverse=ACTCGCCGTTACTRRGG,max-length=1000,name=ITS-v2,min-length=600" \
      --pcr-primers "forward=CGTGAAACCGYTRRAAGGG,reverse=TTGGTCCGTGTTTCAAGACG,max-length=700,name=28S-v2,min-length=300" \
      --ena SRR26955578

The `--pcr-primers` argument takes a string with the format `key1=value1,key2=value2,...`, where the required keys are `forward`, `reverse`, and `name`. Run `sharkmer --help-pcr` for details on all available keys.

The `--max-reads 1000000` argument indicates that the first million reads (individual reads, not read pairs) should be used. This is plenty for nuclear rRNA sequences 18s, 28s, and ITS, since they occur in many copies in the genome, and mitochondrial sequences 16s and co1. Single copy nuclear genes require more data.

This analysis will generate one fasta file for each primer pair, named `{sample}_{gene}.fasta` (e.g., `Stenogorgia_casta_18S.fasta`). If no product was found, the fasta file is not generated. The fasta file can contain more than one sequence when multiple products are found. A YAML stats file (`{sample}.stats.yaml`) is also produced with run statistics and per-gene PCR results.

Reads downloaded via `--ena` are cached (SHA-256 verified) under a shared cache directory so repeated runs on the same accession do not re-download. Use `--cache-dir` to override the location, `--no-cache` to stream directly without touching the cache, or `--clear-cache` to delete cached reads.

You can see all the available built-in PCR panels with:

    sharkmer --list-panels

To export a built-in panel as YAML (for customization or sideloading with `--pcr-panel-file`):

    sharkmer --export-panel cnidaria > octocorallia.yaml

You can then edit the exported YAML to add genes, remove genes, or optimize primer sequences for your study. To run sharkmer with your custom panel, use `--pcr-panel-file`:

    sharkmer --max-reads 1000000 -s Stenogorgia_casta -o output/ --pcr-panel-file octocorallia.yaml --ena SRR26955578

See [PCR.md](PCR.md) for the full panel file schema (versioning, maintainers, changelog, validation block) and the panel development cycle, including how to validate a panel against declared ENA/SRA samples with `scripts/validate_panel.py`. Contributions of new panels and primer additions are welcome — PCR.md describes the expected PR contents.

### Optimizing *in silico* PCR (sPCR)

There are a few different strategies to take if you are not getting a sPCR product, or it is working inconsistently.

The things you should try first are:

- Specify a reasonable value for the `--pcr-primers` parameter `max-length`, which is the maximum expected length of the amplified product, based on what is known about the gene. It does not need to be exact. It does need to be longer than the amplified product, but if it is much too long it may result in amplification of spurious products or runs that take too long without finding anything. Keep in mind that if there are introns, the amplified product can be much longer than the coding sequence.

- Optimize your primer sequences. See [PCR.md](PCR.md) for more information about refining primer sequences and developing primer panels. Make some multiple sequence alignments of the desired sequence region from several closely related species, and refine the primer sequences to be more specific to the target region. You can use [degenerate nucleotide symbols](https://en.wikipedia.org/wiki/Nucleic_acid_notation), such as R for A or G, a variable sites in the site where you would like the sequence to bind. Remember that, just as in real PCR, the forward and reverse primers bind opposite strands with their 3' ends pointing toward each other. The reverse primer sequence should be reverse complemented relative to the reference strand. Adding more degeneracy broadens the taxonomic range a primer can hit, but at a cost: each ambiguous base expands the set of primer variants sharkmer must track, which slows runs down and increases the chance of off-target binding that can cause a run to fail on any particular sample. If your work is focused on a single clade, you will usually get better and faster results by building a panel tailored to that clade that uses only the minimum degeneracy needed to span it. 

- Pick new primers that shorten the region you are trying to amplify, i.e. primers that are closer together in the genome sequence. Shorter amplification fragments tend to require fewer reads to assemble.

- Adjust the `--pcr-primers` parameter `trim`. The default is 15. This is the max number of bases to keep at the 3' end of each primer. Primers used for real PCR tend to be longer than what is required for them to be unique within the genome. This is because they are lengthened to increase melting temperature. If you don't get a product, try reducing `trim` to reduce the specificity of the primer. If you are getting too many spurious products, try increasing this value. Modifying `trim` rather than adjusting the primer sequence makes subsequent adjustments easier (since you don't have to look up the primer sequence again) and also makes the provenance of primer sequences clearer.

- Adjust the number of reads. If you are not getting a product, try increasing the number of reads. You can do this with the `--max-reads` argument. Likewise, if you are getting products and want to speed things up, or you are getting many products for a gene, reduce the number of reads.

### Node budget

The *node budget* (`--node-budget-global`) controls how many nodes the de Bruijn graph is allowed to grow before sharkmer stops extending and moves on. It is the main knob that trades sensitivity against runtime:

- **Too low:** sharkmer gives up before it can bridge the full amplicon, so genes that are present in the data are missed.
- **Too high:** when a gene is absent or primers bind off-target, sharkmer explores a large, fruitless graph before stopping, which wastes time.

The optimal value depends on several factors. Primer specificity is the biggest: longer primers with less degeneracy produce fewer off-target seeds and need less budget. Genome size and complexity also matter — larger, more repetitive genomes generate more graph branches per seed.

By default the budget is set dynamically based on the number of bases ingested, scaling linearly from 100,000 nodes (at ~1 M reads) to 500,000 nodes (at ~5 M reads or more). The auto-selected value is logged at the start of each run. You can override it:

```bash
# Raise the budget for a large, complex genome with degenerate primers
sharkmer --node-budget-global 750000 ...

# Lower it for a quick scan with specific primers on a small genome
sharkmer --node-budget-global 50000 ...
```

If you find that runs are slow but not producing products, lowering the budget will make them fail faster. If runs succeed for some genes but miss others that you expect to be present, raising the budget (and possibly increasing `--max-reads`) is worth trying.

### Other parameters

If these do not work, then you can try adjusting other parameters.

- Specify a reasonable `--pcr-primers` parameter `min-length`. This value defaults to 0, but raising it can get rid of small spurious products.

- Adjust `--min-kmer-count`. This globally filters the kmer table before sPCR, removing all kmers with counts below this value. It defaults to 2, and should generally be at least 2 to avoid one-off sequencing errors. Raising it to 3 or 4 reduces noise but requires more sequencing data.

- Adjust the per-primer `min-count` parameter in `--pcr-primers`. This sets the minimum kmer count required to extend the de Bruijn graph during sPCR for a specific primer pair. It must be at least as high as `--min-kmer-count` (since lower-count kmers have already been filtered). Raising it for a particular gene can help when that gene tends to produce spurious products. For example, `--pcr-primers "forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16S,min-length=500,min-count=4"`

- Adjust the `--pcr-primers` parameter `mismatches`. This defaults to 2. You can try raising it to 3 or 4 if you aren't getting the desired product. This reduces specificity, but this may increase the number of spurious paths that need to be traversed and bog down the run.

By default, sharkmer stops searching as soon as the first connected primer-binding component yields a product (`--pcr-stopping-criteria first-product`). If a gene is producing the wrong product on a difficult sample, try `--pcr-stopping-criteria connected-only` (explores all promising seed components) or `all-components` (exhaustive, slowest but highest sensitivity) — this can surface the correct product when an off-target binding site wins the race at the default setting.

Keep in mind that there is no way to assemble a sPCR product without kmer counts along its full length that meet or exceed the `min-count` parameter. The tool cannot output assembled sequences in the fasta file that are not in the input raw reads from the fastq file. If you are trying to amplify a single copy nuclear gene, that means your sequencing depth (average coverage) of the genome will need to be quite a bit higher than the `min-count` parameter, since there will be fluctuations in coverage along the length of the target region. If coverage at each site is independently distributed, then to have a 95% chance of coverage $\geq 2$ at each site in a region of length $n$, you would need a sequencing depth of 13x for a 1000bp region. That is on the order of 26 million 150 bp reads for a 300Mb genome. This may place single copy nuclear genes out of reach for some organisms with larger genomes, especially if computer RAM limits the number of reads that can be processed.

### Working with complex samples

Some samples contain multiple similar templates that you want to recover as distinct products rather than collapse into a single consensus: metagenomic samples with several related taxa, heterozygous individuals where allelic variants matter, multi-copy gene families where paralogs differ by a handful of bases, or pooled samples. The defaults in sharkmer are tuned for the common "one sample, one product per gene" case and will actively work against you here — they stop at the first product and aggressively collapse near-identical paths.

To pull multiple similar products out of a complex sample:

- Change `--pcr-stopping-criteria` from the `first-product` default to `connected-only` or `all-components` so sharkmer does not bail out after the first successful primer-binding component. `all-components` is exhaustive and slowest but will explore every seed.
- Lower the per-primer `dedup-edit-threshold`. Sharkmer's final deduplication collapses any two output records within this many edits (Levenshtein distance) of each other; the default is 10, which is appropriate when you expect at most one true product per gene but discards closely related variants. Set it lower (e.g. `dedup-edit-threshold=2` or `dedup-edit-threshold=0`) to retain products that differ by only a few bases. For example:

  ```
  --pcr-primers "forward=...,reverse=...,name=CO1,max-length=800,dedup-edit-threshold=2"
  ```

These two settings work together: the stopping criterion decides how many seed regions are explored, and the dedup threshold decides how aggressively the final product set is collapsed. Changing one alone is often not enough for complex samples.

## Incremental kmer counting

Kmer analyses have become an essential component of many routine genomic analyses ([Manekar and Sathe, 2018](https://doi.org/10.1093/gigascience/giy125)). There are multiple excellent highly optimized stand-alone kmer counters, including [Jellyfish](https://github.com/gmarcais/Jellyfish) ([Marçais and Kingsford, 2011](https://doi.org/10.1093/bioinformatics/btr011)) and [KMC](https://github.com/refresh-bio/KMC) ([Kokot et al., 2017](https://doi.org/10.1093/bioinformatics/btx304)). These tools ingest sequence data, such as raw reads, and generate a variety of intermediate products, including count tables of all observed kmers and histograms of observed kmer counts. A variety of downstream tools are available for specific biological analyses. These include genome size estimation with GenomeScope2.0 ([Ranallo-Benavidez et al., 2020](https://doi.org/10.1038/s41467-020-14998-3)).

Kmer spectra analyses typically focus on a single snapshot of data - the complete set of kmers at the time the analysis is performed. Many questions that motivate kmer spectrum analyses are about what happens as data are added. Rather than performing a single analysis on all the data, one can rarefy the data and look at progressively larger nested subsets of reads. This provides more insight into the data in hand, builds better intuition for what changes as data are added, and allows the investigator to better understand what would happen if more data were added. But reanalyzing nested subsets is computationally expensive, since all the data shared across nested subsets are reanalyzed.

Incremental kmer counting, as implemented in `sharkmer`, allows investigators to efficiently investigate the effects of adding data without needing to re-analyze nested subsets. Instead, the data are broken into exclusive subsets. kmer spectra are calculated once for each, and then incrementally combined. 

### Incremental kmer counting example - genome size estimation

Genome size estimation takes a lot of data (about 50x coverage of the genome), and a large amount of RAM. So this example isn't practical on most laptops, given their disk and RAM limitations, and will require a workstation or cluster.

We will use a *Cordagalma ordinatum* dataset from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/SRX10340700). With `--ena`, sharkmer downloads the reads directly:

    sharkmer --chunks 10 -o output/ -s Cordagalma-ordinatum --ena SRR23143278

The `--chunks 10` argument specifies breaking the reads into 10 incremental subsets (histograms are only produced when `--chunks` is greater than 0). The `-o` argument specifies the output directory. The `-s` argument specifies a prefix for the output files.

Alternatively, with local files (gzipped or uncompressed):

    sharkmer --chunks 10 -o output/ -s Cordagalma-ordinatum data/SRR23143278_1.fastq.gz data/SRR23143278_2.fastq.gz

The incremental histogram files in this case will be:

    output/Cordagalma-ordinatum.histo # All the incremental histograms, each in their own column. Suitable for analysis with `sharkmer_viewer.py`.

    output/Cordagalma-ordinatum.final.histo # Just the final histogram with all data. Suitable for analysis with genomescope and other tools.

### Optional Python visualization tools

Optional Python tools for visualizing incremental kmer counting results are available in the [`sharkmer_viewer`](sharkmer_viewer/) folder. See the [sharkmer_viewer readme](sharkmer_viewer/readme.md) for installation and usage instructions.

### Implementation of incremental kmer counting

Here is an overview of how kmer counting works in `sharkmer`:

1. Fastq reads are ingested, split at any `N` bases, and distributed across `n` chunks.
2. Within each chunk, kmers are counted in a hashmap.
3. The hashmaps for the `n` chunks are summed one by one, and a histogram is generated after each chunk of
   counts is added in. This produces `n` histograms, each summarizing more reads than the last.

Incremental kmer counting is used when ingesting reads for other features, including *in silico* PCR.

## Tab completion

`sharkmer` can generate tab-completion scripts for your shell so that flag
names, panel names, and other arguments are completed when you press Tab.

If you installed via **bioconda** (`conda install sharkmer`), completions are
installed automatically — no extra steps needed.

For other installation methods (cargo install, manual build), you need to
generate and install the completion script yourself. First, check which shell
you are using:

```bash
echo $SHELL
```

Then follow the instructions for your shell:

**zsh** (default on macOS, common on Linux):

```bash
# Create the completions directory if it doesn't exist
mkdir -p ~/.zfunc

# Generate the completion script
sharkmer --completions zsh > ~/.zfunc/_sharkmer

# Add to your ~/.zshrc (only needed once):
#   fpath=(~/.zfunc $fpath)
#   autoload -Uz compinit && compinit
```

**bash:**

```bash
# System-wide (may require sudo):
sharkmer --completions bash > /etc/bash_completion.d/sharkmer

# Or per-user:
mkdir -p ~/.local/share/bash-completion/completions
sharkmer --completions bash > ~/.local/share/bash-completion/completions/sharkmer
```

**fish:**

```bash
mkdir -p ~/.config/fish/completions
sharkmer --completions fish > ~/.config/fish/completions/sharkmer.fish
```

After installing, open a new terminal session (or run `source ~/.zshrc` etc.)
for completions to take effect. Completion scripts are a static snapshot of the
available flags, so if you upgrade sharkmer you should regenerate them by
re-running the same command above. Bioconda users don't need to worry about
this — completions are updated automatically with each package upgrade.

## PCR Primer Panels

Primer panels are YAML files that bundle a set of primer pairs for a target
clade (e.g. `cnidaria`, `insecta`, `teleostei`) so they can be run together
with a single `--pcr-panel` flag. Built-in panels ship with sharkmer and can
be listed with `sharkmer --list-panels`; custom panels can be loaded from a
file with `--pcr-panel-file`. See [PCR.md](PCR.md) for the full panel
file schema, validation workflow, and guidance on contributing new panels.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for software development and testing information.

See [PCR.md](PCR.md) for primer, and primer panel, development, testing, and troubleshooting information.
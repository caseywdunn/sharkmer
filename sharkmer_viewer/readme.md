# sharkmer_viewer

Python tools for visualizing genome size estimates from incremental
kmer counting. These are not required to run other features of
`sharkmer`.

Outputs include animations of GenomeScope2 model fits as data are
added:

[![Example GenomeScope2 movie](https://img.youtube.com/vi/SeX5j_xMv5k/0.jpg)](https://www.youtube.com/shorts/SeX5j_xMv5k)

## Manual installation

Create a conda environment with the required dependencies:

    conda create -n shark -c conda-forge -c bioconda python==3.10 ffmpeg genomescope2
    conda activate shark
    cd sharkmer_viewer/
    pip install .

Activate the environment before running any commands below:

    conda activate shark

## Docker installation

The Dockerfile builds an image with all dependencies (GenomeScope2, R,
ffmpeg, and the `sharkmer_viewer` Python package) pre-installed and
ready to use. No conda activation is needed inside the container.

    cd sharkmer_viewer/
    docker build -t sharkmer_viewer .

To run commands inside the container, mount the directory with your
sharkmer output files as `/data`:

    docker run -v /path/to/sharkmer_output:/data sharkmer_viewer <command>

For example:

    docker run -v $(pwd):/data sharkmer_viewer sharkmer_viewer -d Nanomia-bijuga.histo -s Nanomia-bijuga.stats.yaml

Or open an interactive shell:

    docker run -it -v $(pwd):/data sharkmer_viewer bash

All of the usage examples below show the command by itself. If you are
using Docker, prepend `docker run -v $(pwd):/data sharkmer_viewer` to
each command.

## Use

These tools visualize the output of `sharkmer` incremental kmer
counting. To generate that input, run `sharkmer` with the `--chunks`
flag:

    sharkmer --chunks 10 -s Cordagalma-ordinatum -o output/ data/*.fastq.gz

This produces two histogram files in the output directory:

- `Cordagalma-ordinatum.histo` — the incremental histogram, a
  tab-separated file where column 1 is the kmer coverage and columns
  2 through N+1 are the kmer frequency counts for each of the N
  cumulative chunks. Each successive column includes more reads than
  the last.
- `Cordagalma-ordinatum.final.histo` — the final histogram (all reads
  combined), in standard two-column format compatible with GenomeScope2
  and other tools.

It also produces a stats file (`Cordagalma-ordinatum.stats.yaml`) with
run metadata including total bases read.

### Interactive HTML visualization

`sharkmer_viewer` reads the incremental histogram file and the stats
file and produces interactive HTML plots that you can open in a browser.

    sharkmer_viewer -d Cordagalma-ordinatum.histo -s Cordagalma-ordinatum.stats.yaml

Options:

- `-d` — histogram file (required)
- `-s` — stats file, `.stats.yaml` (new format) or `.stats` TSV (legacy)
- `-n` — display name for plots
- `-o` — output file base name
- `-g` — reference genome size in Mb (for scaling)

This produces two HTML files:

- `{name}.html` — animated histogram with play/pause controls. Shows
  the kmer spectrum at each incremental chunk, with peak and valley
  markers. Hover over data points for exact values.
- `{name}_genome_size.html` — line plots of heterozygous and
  homozygous genome size estimates as a function of the number of reads
  included.

### GenomeScope2 movie

The `genomescopemovie.sh` script takes the same incremental histogram
file and runs GenomeScope2 on each chunk to produce an MP4 movie of
the model fit evolving as data are added.

    bash genomescopemovie.sh -i Cordagalma-ordinatum.histo -t 8 -o gs_output -k 21

Options:

- `-i` — input histogram file (required)
- `-t` — number of parallel GenomeScope2 jobs (default: 1)
- `-o` — output directory (default: input filename without extension)
- `-k` — kmer size (default: 21)

The script runs four stages:

**1. Split histograms.** The incremental histogram file is split into
one two-column file per chunk (`sample_0001.histo`,
`sample_0002.histo`, ...). Each contains the kmer coverage and the
cumulative frequency counts for that chunk.

**2. Run GenomeScope2.** GenomeScope2 is run independently on each
per-chunk histogram. It fits the kmer spectrum model and estimates
genome size, heterozygosity, and repeat content. Runs up to `-t` jobs
in parallel. Each run produces:

- `sample_XXXX_linear_plot.png` — linear-scale kmer spectrum plot
- `sample_XXXX_log_plot.png` — log-scale kmer spectrum plot
- `sample_XXXX_summary.txt` — text summary with estimated parameters

**3. Consolidate statistics.** The per-chunk summary files are parsed
into a single TSV (`{name}_genomescope_stats.tsv`) with columns for
each GenomeScope2 parameter (homozygous/heterozygous coverage, haploid
genome length, repeat length, unique length, model fit, read error
rate), each with min and max values.

**4. Generate movies.** The per-chunk PNG plots are combined with
ffmpeg into two MP4 movies:

- `{name}_linear_plot.mp4` — animation of linear-scale plots
- `{name}_log_plot.mp4` — animation of log-scale plots

Each frame corresponds to one incremental chunk, so the movie shows the
GenomeScope2 fit stabilizing as more reads are added.

### Standalone GenomeScope2

The final histogram can also be analyzed directly with GenomeScope2
without the movie pipeline:

    genomescope2 -i Cordagalma-ordinatum.final.histo -s Cordagalma-ordinatum -k 21

## Development

Follows model at https://github.com/shchurch/countland/blob/main/development.md

See also https://github.com/dunnlab/phylopytho

To debug in VS Code, open the `sharkmer_viewer` subdirectory rather
than the whole repo.

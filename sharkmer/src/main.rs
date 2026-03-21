use anyhow::{bail, ensure, Context, Result};
use bio::io::fasta;
use clap::Parser;
use log::{debug, info, warn};
use pcr::preconfigured;
use std::io::BufRead;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

use crate::kmer::Chunk;
use crate::kmer::KmerCounts;

mod kmer;
mod pcr;

const N_READS_PER_BATCH: u64 = 1000;

#[allow(dead_code)]
fn is_valid_nucleotide(c: char) -> bool {
    match c {
        'A' => true,
        'C' => true,
        'G' => true,
        'T' => true,
        'R' => true, // A or G
        'Y' => true, // C or T
        'S' => true, // C or G
        'W' => true, // A or T
        'K' => true, // G or T
        'M' => true, // A or C
        'B' => true, // C or G or T
        'D' => true, // A or G or T
        'H' => true, // A or C or T
        'V' => true, // A or C or G
        'N' => true, // A or C or G or T
        _ => false,
    }
}

pub fn parse_pcr_string(pcr_string: &str) -> Result<Vec<pcr::PCRParams>> {
    if pcr_string.is_empty() {
        bail!("Invalid empty pcr string");
    }

    // Split the string on commas
    let split: Vec<&str> = pcr_string.split(',').collect();

    // If the string is a single element, check if it is a preconfigured panel
    if split.len() == 1 {
        // If OK, validate and return the panel
        // If not OK, return an error

        match preconfigured::get_panel(pcr_string) {
            Ok(params) => return Ok(params),
            Err(_) => bail!("Invalid preconfigured PCR panel: {}", pcr_string),
        }
    }

    let mut forward_seq = "".to_string();
    let mut reverse_seq = "".to_string();
    let mut gene_name = "".to_string();
    let mut max_length = 10000;
    let mut min_length = 0;
    let mut min_coverage = 2;
    let mut mismatches = 2;
    let mut trim = 15;
    let mut citation = "".to_string();
    let mut notes = "".to_string();
    let mut dedup_edit_threshold = pcr::DEFAULT_DEDUP_EDIT_THRESHOLD;

    // Loop over additional parameters, which are of the form key=value and are separated by commas
    for item in split.iter() {
        let key_value: Vec<&str> = item.split('=').collect();
        if key_value.len() != 2 {
            bail!("Invalid parameter, should be in format key=value: {}", item);
        }

        let key = key_value[0].to_lowercase();
        let key = key.as_str();
        let value = key_value[1];

        match key {
            "name" => {
                gene_name = value.to_string();
            }
            "forward" => {
                forward_seq = value.to_string().to_uppercase();
            }
            "reverse" => {
                reverse_seq = value.to_string().to_uppercase();
            }
            "max-length" => {
                max_length = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "min-length" => {
                min_length = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "min-coverage" => {
                min_coverage = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "mismatches" => {
                mismatches = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "trim" => {
                trim = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "citation" => {
                citation = value.to_string();
            }
            "notes" => {
                notes = value.to_string();
            }
            "dedup-edit-threshold" => {
                dedup_edit_threshold = value
                    .parse()
                    .with_context(|| format!("Invalid value for {}: {}", key, value))?;
            }
            "panel" => {
                if split.len() > 1 {
                    bail!("Invalid --pcr arguments, if a panel is specified it should the the only argument: {}", pcr_string);
                }
                match preconfigured::get_panel(value) {
                    Ok(params) => return Ok(params),
                    Err(_) => bail!("Invalid preconfigured PCR panel: {}", value),
                }
            }
            _ => {
                bail!("Unexpected parameter: {}", key);
            }
        }
    }

    let pcr_params = pcr::PCRParams {
        forward_seq,
        reverse_seq,
        min_length,
        max_length,
        gene_name,
        min_coverage,
        mismatches,
        trim,
        citation,
        notes,
        dedup_edit_threshold,
    };

    pcr::validate_pcr_params(&pcr_params)?;

    Ok(vec![pcr_params])
}

/// A collection of kmer counting and analysis tools
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// k-mer length
    #[arg(short, default_value_t = 21)]
    k: usize,

    /// Maximum value for histogram
    #[arg(long, default_value_t = 10000)]
    histo_max: u64,

    /// Number of chunks to divide the data into
    #[arg(short, default_value_t = 100)]
    n: usize,

    /// Maximum number of reads to process
    #[arg(short = 'm', long, default_value_t = 0)]
    max_reads: u64,

    /// Number of threads to use
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Name of the sample, could be species and sample ID eg Nanomia-bijuga-YPMIZ035039.
    /// Will be used as the prefix for output files
    #[arg(short, long, default_value_t = String::from("sample") )]
    sample: String,

    /// Directory for output files, defaults to current directory
    #[arg(short, long, default_value_t = String::from("./") )]
    outdir: String,

    /// Input files, fastq. If no files are specified, data will be read
    /// from stdin. This can be used to uncompress a gz file and send them
    /// to sharkmer.
    #[arg()]
    input: Option<Vec<String>>,

    /// Optional primer pairs and parameters for in silico PCR (sPCR).
    ///
    /// To see a list of available preconfigured panels and exit, use:
    ///
    ///    --pcr panels
    ///
    /// To use a specific preconfigured panel, specify it by name with panel= . For example:
    ///
    ///    --pcr panel=cnidaria
    ///
    /// To manually specify a single primer pair, use the format:
    ///
    ///    --pcr "key1=value1,key2=value2,key3=value3,..."
    ///
    /// For example:
    ///
    ///    --pcr "forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16s,min-length=500"
    ///
    /// Where required keys are:
    ///
    ///   forward is the forward primer sequence in 5' to 3' orientation
    ///
    ///   reverse is the reverse primer sequence in 5' to 3' orientation
    ///     along the opposite strand as the forward primer, so that the
    ///     primers are in the same orientation that you would use in an
    ///     actual in vitro PCR reaction (3' ends facing each other).
    ///
    ///   name is a unique name for the primer pair or amplified gene
    ///     region. This will be used to specify amplified regions in
    ///     the output fasta file.
    ///
    /// The following keys are optional:
    ///
    ///    min-length is the minimum length of the PCR product, including
    ///     the primers. Default is 0.
    ///
    ///    max-length is the maximum length of the PCR product, including
    ///     the primers. Default is 10000.
    ///
    ///    min-coverage: minimum coverage for a kmer to be included in the
    ///      amplified region. Default is 2.
    ///
    ///    mismatches: maximum number of mismatches allowed between the
    ///      primer and the kmer. Default is 2.
    ///
    ///    trim: number of bases to keep at the 3' end of each primer.
    ///      Default is 15.
    ///
    /// More than one primer pair can be specified, for example:
    ///
    ///    --pcr "..." --pcr "..."
    #[arg(short = 'p', long)]
    pcr: Vec<String>,

    /// Load primer panel(s) from YAML file(s). Can be specified multiple times.
    #[arg(long)]
    pcr_file: Vec<String>,

    /// Validate FASTQ format every N records (0 = first record only)
    #[arg(long, default_value_t = 0)]
    validate_every: u64,

    /// Increase verbosity (-v info, -vv debug, -vvv trace)
    #[arg(short = 'v', long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Print citation information and exit
    #[arg(long)]
    cite: bool,
}

const FASTA_LINE_WIDTH: usize = 80;

/// Write a FASTA record with sequence lines wrapped at FASTA_LINE_WIDTH characters.
fn write_fasta_record(writer: &mut impl Write, record: &fasta::Record) -> std::io::Result<()> {
    write!(writer, ">{}", record.id())?;
    if let Some(desc) = record.desc() {
        write!(writer, " {}", desc)?;
    }
    writeln!(writer)?;
    for chunk in record.seq().chunks(FASTA_LINE_WIDTH) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Validate a FASTQ record (4 lines). Returns an error if the record is malformed.
fn validate_fastq_record(
    header: &str,
    _sequence: &str,
    separator: &str,
    quality: &str,
    sequence_len: usize,
    record_num: u64,
) -> Result<()> {
    ensure!(
        header.starts_with('@'),
        "FASTQ record {} has invalid header (expected '@', got '{}'): {}",
        record_num + 1,
        header.chars().next().unwrap_or(' '),
        header,
    );
    ensure!(
        separator.starts_with('+'),
        "FASTQ record {} has invalid separator line (expected '+', got '{}'): {}",
        record_num + 1,
        separator.chars().next().unwrap_or(' '),
        separator,
    );
    ensure!(
        quality.len() == sequence_len,
        "FASTQ record {} has mismatched sequence ({}) and quality ({}) lengths",
        record_num + 1,
        sequence_len,
        quality.len(),
    );
    Ok(())
}

/// Mutable state for FASTQ reading, shared across multiple input sources.
struct FastqReadState {
    chunks: Vec<kmer::Chunk>,
    chunk_index: usize,
    reads: Vec<kmer::Read>,
    n_reads_read: u64,
    n_bases_read: u64,
}

/// Read FASTQ records from a buffered reader, ingesting sequences into chunks.
///
/// Reads 4 lines at a time (header, sequence, separator, quality) and validates
/// the format. Returns true if max_reads was reached.
fn read_fastq<R: BufRead>(
    reader: R,
    state: &mut FastqReadState,
    max_reads: u64,
    validate_every: u64,
    source_name: &str,
) -> Result<bool> {
    let mut lines = reader.lines();
    let n_chunks = state.chunks.len();

    while let Some(header_result) = lines.next() {
        // Read 4 lines for one FASTQ record
        let header =
            header_result.with_context(|| format!("Failed to read header from {}", source_name))?;

        let sequence = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read sequence from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing sequence line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let separator = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read separator from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing separator line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        let quality = match lines.next() {
            Some(line) => {
                line.with_context(|| format!("Failed to read quality from {}", source_name))?
            }
            None => bail!(
                "Truncated FASTQ record at record {} in {}: missing quality line",
                state.n_reads_read + 1,
                source_name,
            ),
        };

        // Validate the record
        let should_validate = state.n_reads_read == 0
            || (validate_every > 0 && state.n_reads_read.is_multiple_of(validate_every));
        if should_validate {
            validate_fastq_record(
                &header,
                &sequence,
                &separator,
                &quality,
                sequence.len(),
                state.n_reads_read,
            )?;
        }

        // Process the sequence
        state.n_bases_read += sequence.len() as u64;
        let new_reads = kmer::seq_to_reads(&sequence)?;
        state.reads.extend(new_reads);
        state.n_reads_read += 1;

        // If we have read enough reads, ingest them into current chunk
        if state.n_reads_read.is_multiple_of(N_READS_PER_BATCH) {
            state.chunks[state.chunk_index].ingest_reads(&state.reads)?;
            state.chunk_index += 1;
            if state.chunk_index == n_chunks {
                state.chunk_index = 0;
            }
            state.reads.clear();
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

fn main() -> Result<()> {
    let start_run = std::time::Instant::now();

    // Ingest command line arguments
    let args = Args::parse();

    // Initialize logging based on verbosity level
    let log_level = match args.verbose {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    env_logger::Builder::new()
        .filter_level(log_level)
        .format_timestamp(None)
        .init();

    // Print citation information and exit if --cite is specified
    if args.cite {
        println!("{} {}\n", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
        println!("If you use sharkmer, please cite:\n");
        println!("  Dunn et al. (2026) sharkmer: a kmer counter and seeded de Bruijn graph");
        println!("  assembler for in silico PCR and incremental kmer counting.");
        println!("  doi: 10.xxxx/xxxxx\n");
        println!("BibTeX:");
        println!("  @article{{dunn2026sharkmer,");
        println!("    title={{sharkmer: a kmer counter and seeded de Bruijn graph assembler");
        println!("           for in silico PCR and incremental kmer counting}},");
        println!("    author={{Dunn, Casey W.}},");
        println!("    year={{2026}},");
        println!("    doi={{10.xxxx/xxxxx}}");
        println!("  }}");
        std::process::exit(0);
    }

    info!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    debug!("{:?}", args);

    // Parse the outdir path and sample, create directories if necessary
    let path = PathBuf::from(&args.outdir);

    // Create the output directory if it does not exist
    let directory = format!(
        "{}/",
        path.to_str()
            .context("Output directory path contains invalid UTF-8")?
    );
    std::fs::create_dir_all(&directory)
        .with_context(|| format!("Failed to create output directory: {}", directory))?;

    let k = args.k;

    // Check that the arguments are valid
    ensure!(
        k < 32,
        "k must be less than 32 due to use of 64 bit integers to encode kmers"
    );
    ensure!(k > 0, "k must be greater than 0");
    ensure!(k % 2 == 1, "k must be odd");
    ensure!(args.histo_max > 0, "histo_max must be greater than 0");
    ensure!(
        args.histo_max <= 1_000_000,
        "histo_max must not exceed 1000000, got {}",
        args.histo_max
    );
    ensure!(args.n > 0, "n must be greater than 0");

    // Create an empty data frame for pcr runs
    let mut pcr_runs: Vec<pcr::PCRParams> = Vec::new();

    // Loop over the pcr strings, check that they are valid, and add each to the pcr_runs vector
    for pcr_string in args.pcr.iter() {
        if pcr_string == "panels" {
            println!("Available preconfigured PCR panels:");
            preconfigured::print_pcr_panels();
            std::process::exit(0);
        }
        let pcr_params = parse_pcr_string(pcr_string)
            .with_context(|| format!("Error parsing pcr string: {}", pcr_string))?;
        pcr_runs.extend(pcr_params);
    }

    // Load primer panels from YAML files
    for pcr_file in args.pcr_file.iter() {
        let pcr_params = preconfigured::load_panel_file(pcr_file)
            .with_context(|| format!("Error loading panel file: {}", pcr_file))?;
        pcr_runs.extend(pcr_params);
    }

    // Check that there are no duplicate gene names
    let mut gene_names: Vec<String> = Vec::new();
    for pcr_params in pcr_runs.iter() {
        ensure!(
            !gene_names.contains(&pcr_params.gene_name),
            "Duplicate gene name: {}",
            pcr_params.gene_name
        );
        gene_names.push(pcr_params.gene_name.clone());
    }

    // Validate input files exist before starting processing
    if let Some(input_files) = &args.input {
        for file_name in input_files.iter() {
            let file_path = Path::new(file_name);
            ensure!(
                file_path.exists(),
                "Input file does not exist: {}",
                file_name
            );
            ensure!(
                file_path.is_file(),
                "Input path is not a file: {}",
                file_name
            );
        }
    }

    // Validate pcr-file paths exist
    for pcr_file in args.pcr_file.iter() {
        let file_path = Path::new(pcr_file);
        ensure!(
            file_path.exists(),
            "PCR panel file does not exist: {}",
            pcr_file
        );
    }

    // Set the number of threads for Rayon to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Ingest the fastq data
    let start = std::time::Instant::now();
    info!("Ingesting reads...");

    let mut chunks: Vec<kmer::Chunk> = Vec::new();
    for _ in 0..args.n {
        chunks.push(Chunk::new(&k));
    }

    let mut state = FastqReadState {
        chunks,
        chunk_index: 0,
        reads: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
    };

    match &args.input {
        Some(input_files) => {
            // read from one or more files
            for file_name in input_files.iter() {
                let file_path = Path::new(&file_name);
                let file = std::fs::File::open(file_path)
                    .with_context(|| format!("Failed to open file: {}", file_name))?;

                let reader: Box<dyn BufRead> = if file_name.ends_with(".gz") {
                    Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
                } else {
                    Box::new(std::io::BufReader::new(file))
                };
                let reached_max = read_fastq(
                    reader,
                    &mut state,
                    args.max_reads,
                    args.validate_every,
                    file_name,
                )?;
                if reached_max {
                    break;
                }
            }
        }
        None => {
            // read from stdin
            let stdin = std::io::stdin();
            let handle = stdin.lock();

            read_fastq(
                handle,
                &mut state,
                args.max_reads,
                args.validate_every,
                "stdin",
            )?;
        }
    }

    // Ingest any remaining reads
    state.chunks[state.chunk_index].ingest_reads(&state.reads)?;
    state.reads.clear();

    let mut n_reads_ingested: u64 = 0;
    let mut n_bases_ingested: u64 = 0;
    let mut n_kmers_ingested: u64 = 0;
    for chunk in state.chunks.iter() {
        n_reads_ingested += chunk.get_n_reads();
        n_bases_ingested += chunk.get_n_bases();
        n_kmers_ingested += chunk.get_n_kmers();
    }

    info!(
        "Read {} reads, {} bases",
        state.n_reads_read, state.n_bases_read
    );
    info!(
        "Ingested {} subreads, {} bases, {} kmers",
        n_reads_ingested, n_bases_ingested, n_kmers_ingested
    );
    info!("Time to ingest reads: {:?}", start.elapsed());

    if n_reads_ingested == 0 {
        warn!("No reads were ingested. All output will be empty.");
    }

    // Create the histograms
    info!("Consolidating chunks and creating histograms...");
    let start = std::time::Instant::now();

    let mut kmer_counts: KmerCounts = KmerCounts::new(&k);
    let mut histos: Vec<kmer::Histogram> = Vec::with_capacity(args.n);

    // Drain the chunks so each chunk's hash table is freed after merging
    for chunk in state.chunks.drain(..) {
        kmer_counts.extend(chunk.get_kmer_counts())?;
        drop(chunk);

        let histo = kmer::Histogram::from_kmer_counts(&kmer_counts, &args.histo_max);

        histos.push(histo);
    }
    info!("Chunks consolidated, time: {:?}", start.elapsed());

    let n_hashed_kmers: u64 = kmer_counts.get_n_kmers();
    info!(
        "{} unique kmers with a total count of {} were found",
        kmer_counts.get_n_unique_kmers(),
        n_hashed_kmers
    );

    ensure!(
        n_hashed_kmers == n_kmers_ingested,
        "The total count of hashed kmers ({}) does not equal the number of ingested kmers ({})",
        n_hashed_kmers,
        n_kmers_ingested,
    );

    // Write the histograms to a tab delimited file, with the first column being the count
    // Skip the first row, which is the count of 0. Do not include a header
    info!("Writing histograms to file...");
    let mut file = std::fs::File::create(format!("{}{}.histo", directory, args.sample))
        .context("Failed to create histogram file")?;
    for i in 1..args.histo_max as usize + 2 {
        let mut line = format!("{}", i);
        for histo in histos.iter() {
            let histo_vec = kmer::Histogram::get_vector(histo)?;
            line = format!("{}\t{}", line, histo_vec[i]);
        }
        line = format!("{}\n", line);
        file.write_all(line.as_bytes())
            .context("Failed to write histogram data")?;
    }

    // Write the final histogram to a file, ready for GenomeScope2 etc...
    info!("Writing final histogram to file...");

    let last_histo = &histos[histos.len() - 1];
    let last_histo_vec = kmer::Histogram::get_vector(last_histo)?;

    let mut file = std::fs::File::create(format!("{}{}.final.histo", directory, args.sample))
        .context("Failed to create final histogram file")?;
    for (i, value) in last_histo_vec
        .iter()
        .enumerate()
        .skip(1)
        .take(args.histo_max as usize + 1)
    {
        let line = format!("{}\t{}\n", i, value);
        file.write_all(line.as_bytes())
            .context("Failed to write final histogram data")?;
    }

    let n_singleton_kmers = last_histo_vec[1];
    let n_unique_kmers_histo: u64 = last_histo.get_n_unique_kmers();
    let n_kmers_histo: u64 = last_histo.get_n_kmers();

    debug!("{} unique kmers in histogram", n_unique_kmers_histo);
    debug!("{} kmers in histogram", n_kmers_histo);

    ensure!(
        n_kmers_histo == n_kmers_ingested,
        "The total count of kmers in the histogram ({}) does not equal the total expected count of kmers ({})",
        n_kmers_histo,
        n_kmers_ingested,
    );

    ensure!(
        n_unique_kmers_histo == kmer_counts.get_n_unique_kmers(),
        "The total count of unique kmers in the histogram ({}) does not equal the total count of hashed kmers ({})",
        n_unique_kmers_histo,
        kmer_counts.get_n_unique_kmers(),
    );

    info!("Writing stats to file...");
    let mut file_stats = std::fs::File::create(format!("{}{}.stats", directory, args.sample))
        .context("Failed to create stats file")?;
    let mut line = format!("arguments\t{:?}\n", args);
    line = format!("{}kmer_length\t{}\n", line, args.k);
    line = format!("{}n_reads_read\t{}\n", line, state.n_reads_read);
    line = format!("{}n_bases_read\t{}\n", line, state.n_bases_read);
    line = format!("{}n_subreads_ingested\t{}\n", line, n_reads_ingested);
    line = format!("{}n_bases_ingested\t{}\n", line, n_bases_ingested);
    line = format!("{}n_kmers\t{}\n", line, n_kmers_ingested);
    line = format!(
        "{}n_multi_kmers\t{}\n",
        line,
        n_kmers_ingested.saturating_sub(n_singleton_kmers)
    );
    line = format!("{}n_singleton_kmers\t{}\n", line, n_singleton_kmers);

    file_stats
        .write_all(line.as_bytes())
        .context("Failed to write stats file")?;
    if !pcr_runs.is_empty() {
        info!("Running in silico PCR...");

        // Prep kmer counts for in silico PCR. Remove singleton kmers (to reduce size) and
        // add reverse complements
        let min_count: u64 = 2;
        let kmer_counts_pcr = kmer_counts.get_pcr_kmers(&min_count);

        for pcr_params in pcr_runs.iter() {
            let fasta = pcr::do_pcr(&kmer_counts_pcr, &args.sample, pcr_params)?;

            if !fasta.is_empty() {
                let fasta_path = format!(
                    "{}{}_{}.fasta",
                    directory, args.sample, pcr_params.gene_name
                );
                let mut file = std::fs::File::create(&fasta_path)
                    .with_context(|| format!("Failed to create FASTA file: {}", fasta_path))?;
                for record in fasta {
                    write_fasta_record(&mut file, &record)
                        .context("Failed to write FASTA record")?;
                }
            }
        }

        info!("Done running in silico PCR");
    }

    info!("Total run time: {:?}", start_run.elapsed());
    Ok(())
}

#[cfg(test)]
mod tests {}

use anyhow::{bail, ensure, Context, Result};
use bio::io::fasta;
use clap::{CommandFactory, Parser};
use clap_complete::Shell;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, info, warn};
use pcr::preconfigured;
use peak_alloc::PeakAlloc;
use rayon::prelude::*;
use serde::Serialize;
use std::io::BufRead;
use std::io::IsTerminal;
use std::io::Write;
use std::path::PathBuf;

use crate::kmer::Chunk;
use crate::kmer::KmerCounts;

mod kmer;
mod pcr;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

const N_READS_PER_BATCH: u64 = 1000;

/// Stats for a single PCR gene result, included in the YAML stats output.
#[derive(Serialize)]
struct PcrGeneResult {
    gene_name: String,
    status: String,
    n_products: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    product_lengths: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    output_file: Option<String>,
}

/// Structured run statistics, serialized as YAML.
#[derive(Serialize)]
struct RunStats {
    sharkmer_version: String,
    command: String,
    sample: String,
    kmer_length: usize,
    chunks: usize,
    n_reads_read: u64,
    n_bases_read: u64,
    n_subreads_ingested: u64,
    n_bases_ingested: u64,
    n_kmers: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    n_multi_kmers: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    n_singleton_kmers: Option<u64>,
    peak_memory_bytes: u64,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pcr_results: Vec<PcrGeneResult>,
}

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

/// Parse an inline primer specification string (key=value,key=value,...) into PCRParams.
pub fn parse_pcr_primers_string(pcr_string: &str) -> Result<pcr::PCRParams> {
    ensure!(!pcr_string.is_empty(), "Invalid empty primer specification");

    let split: Vec<&str> = pcr_string.split(',').collect();

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

    Ok(pcr_params)
}

/// A tool for kmer counting and in silico PCR (sPCR)
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(after_help = "\
Output files:\n  \
  {outdir}/{sample}_{panel}_{gene}.fasta  sPCR products per gene\n  \
  {outdir}/{sample}.histo                 Incremental histograms (--chunks > 0)\n  \
  {outdir}/{sample}.final.histo           Final histogram (--chunks > 0)\n  \
  {outdir}/{sample}.stats.yaml             Run statistics (YAML)")]
struct Args {
    /// FASTQ input files (.fastq or .fastq.gz). If omitted, reads from stdin.
    #[arg(help_heading = "Input")]
    input: Option<Vec<PathBuf>>,

    /// Stream reads directly from ENA by SRA/ENA accession (e.g. SRR5324768)
    #[arg(long, help_heading = "Input")]
    sra: Option<String>,

    /// Sample name (used as output file prefix, required for processing)
    #[arg(short, long, help_heading = "Output")]
    sample: Option<String>,

    /// Output directory
    #[arg(short, long, default_value = "./", help_heading = "Output")]
    outdir: PathBuf,

    /// Use a preconfigured primer panel (repeatable)
    #[arg(long, help_heading = "PCR")]
    pcr_panel: Vec<String>,

    /// Load primers from a YAML file (repeatable)
    #[arg(long, help_heading = "PCR")]
    pcr_file: Vec<PathBuf>,

    /// Specify a primer pair inline (repeatable, see --help-pcr)
    #[arg(long, help_heading = "PCR")]
    pcr_primers: Vec<String>,

    /// List available primer panels and exit
    #[arg(long, help_heading = "PCR")]
    list_panels: bool,

    /// Export a built-in panel as YAML to stdout and exit
    #[arg(long, help_heading = "PCR")]
    export_panel: Option<String>,

    /// Show detailed help for --pcr-primers format
    #[arg(long, help_heading = "PCR")]
    help_pcr: bool,

    /// Kmer length
    #[arg(short, default_value_t = 21, help_heading = "Kmer counting")]
    k: usize,

    /// Number of incremental chunks (0 = skip histograms)
    #[arg(long, default_value_t = 0, help_heading = "Kmer counting")]
    chunks: usize,

    /// Maximum histogram count value
    #[arg(long, default_value_t = 10000, help_heading = "Kmer counting")]
    histo_max: u64,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 1, help_heading = "General")]
    threads: usize,

    /// Maximum number of reads to process (default: all)
    #[arg(short = 'm', long, help_heading = "General")]
    max_reads: Option<u64>,

    /// Minimum kmer count for sPCR (filters low-count kmers before PCR)
    #[arg(long, default_value_t = 2, help_heading = "General")]
    min_kmer_count: u64,

    /// Validate FASTQ format every N records (0 = first record only)
    #[arg(long, default_value_t = 0, help_heading = "General")]
    validate_every: u64,

    /// Increase verbosity (-v info, -vv debug, -vvv trace)
    #[arg(short = 'v', long, action = clap::ArgAction::Count, help_heading = "General")]
    verbose: u8,

    /// Control color output: auto (default), always, never
    #[arg(long, default_value = "auto", help_heading = "General")]
    color: ColorMode,

    /// Print citation information and exit
    #[arg(long, help_heading = "General")]
    cite: bool,

    /// Generate shell completions and exit (bash, zsh, fish, elvish, powershell)
    #[arg(long, help_heading = "General", value_name = "SHELL")]
    generate_completions: Option<Shell>,

    /// Validate inputs and print what would happen, then exit
    #[arg(long, help_heading = "General")]
    dry_run: bool,
}

/// Color output mode for log messages.
#[derive(Debug, Clone, clap::ValueEnum)]
enum ColorMode {
    Auto,
    Always,
    Never,
}

const FASTA_LINE_WIDTH: usize = 80;

/// Query ENA filereport API to get FASTQ URLs for an SRA/ENA accession.
/// Returns URLs ordered by mate number (R1 first, then R2 if paired-end).
fn get_ena_fastq_urls(accession: &str) -> Result<Vec<String>> {
    let url = format!(
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=run_accession,fastq_ftp",
        accession
    );

    info!("Querying ENA for accession {}...", accession);
    let response = ureq::get(&url)
        .call()
        .with_context(|| format!("Failed to query ENA API for accession {}", accession))?;

    let body = response
        .into_string()
        .context("Failed to read ENA API response")?;

    // Response is TSV: header row, then data rows
    // Fields: run_accession\tfastq_ftp
    let lines: Vec<&str> = body.lines().collect();
    ensure!(
        lines.len() >= 2,
        "ENA returned no results for accession '{}'. Check that the accession is valid.",
        accession
    );

    // Parse the fastq_ftp field (semicolon-separated URLs)
    let data_line = lines[1];
    let fields: Vec<&str> = data_line.split('\t').collect();
    ensure!(
        fields.len() >= 2 && !fields[1].is_empty(),
        "ENA returned no FASTQ URLs for accession '{}'. The run may not have public FASTQ files.",
        accession
    );

    let urls: Vec<String> = fields[1]
        .split(';')
        .map(|u| {
            if u.starts_with("ftp://") || u.starts_with("http") {
                u.to_string()
            } else {
                format!("http://{}", u)
            }
        })
        .collect();

    info!(
        "Found {} FASTQ file(s) for {}: {}",
        urls.len(),
        accession,
        urls.join(", ")
    );

    Ok(urls)
}

/// Format a count with comma-separated thousands (e.g. 1,234,567).
fn format_count(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Format bytes as a human-readable string (e.g. "1.2 GB").
fn format_bytes(bytes: u64) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    const GB: f64 = MB * 1024.0;
    let b = bytes as f64;
    if b < KB {
        format!("{} B", bytes)
    } else if b < MB {
        format!("{:.1} KB", b / KB)
    } else if b < GB {
        format!("{:.1} MB", b / MB)
    } else {
        format!("{:.1} GB", b / GB)
    }
}

/// Warn if an output file already exists (it will be overwritten).
fn warn_if_exists(path: &str) {
    if std::path::Path::new(path).exists() {
        warn!("Overwriting existing file {}", path);
    }
}

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
    if header.starts_with('>') {
        bail!(
            "Input appears to be FASTA format, not FASTQ (record {} starts with '>'). \
             sharkmer requires FASTQ input with quality scores.",
            record_num + 1,
        );
    }
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
    progress: &ProgressBar,
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
            progress.set_position(state.n_reads_read);
        }

        if max_reads > 0 && state.n_reads_read >= max_reads {
            progress.set_position(state.n_reads_read);
            return Ok(true); // reached max reads
        }
    }

    Ok(false) // did not reach max reads (EOF)
}

fn main() -> Result<()> {
    let start_run = std::time::Instant::now();

    // Ingest command line arguments
    let args = Args::parse();

    // Initialize logging based on verbosity level and color mode
    let log_level = match args.verbose {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    let write_style = match args.color {
        ColorMode::Auto => env_logger::WriteStyle::Auto,
        ColorMode::Always => env_logger::WriteStyle::Always,
        ColorMode::Never => env_logger::WriteStyle::Never,
    };
    env_logger::Builder::new()
        .filter_level(log_level)
        .write_style(write_style)
        .format_timestamp(None)
        .init();

    // Handle early-exit flags before any processing

    // Generate shell completions and exit
    if let Some(shell) = args.generate_completions {
        let mut cmd = Args::command();
        clap_complete::generate(shell, &mut cmd, "sharkmer", &mut std::io::stdout());
        std::process::exit(0);
    }

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

    // List available panels and exit
    if args.list_panels {
        println!("Available preconfigured PCR panels:");
        preconfigured::print_pcr_panels();
        std::process::exit(0);
    }

    // Export a built-in panel as YAML and exit
    if let Some(panel_name) = &args.export_panel {
        let yaml = preconfigured::export_panel_yaml(panel_name)
            .with_context(|| format!("Failed to export panel: {}", panel_name))?;
        print!("{}", yaml);
        std::process::exit(0);
    }

    // Show detailed help for --pcr-primers format
    if args.help_pcr {
        println!("Inline primer specification format for --pcr-primers:\n");
        println!("  --pcr-primers \"key1=value1,key2=value2,...\"\n");
        println!("Example:");
        println!("  --pcr-primers \"forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16s,min-length=500\"\n");
        println!("Required keys:");
        println!("  forward       Forward primer sequence (5' to 3')");
        println!("  reverse       Reverse primer sequence (5' to 3' on opposite strand)");
        println!("  name          Unique name for the primer pair or gene region\n");
        println!("Optional keys:");
        println!("  min-length              Minimum product length including primers [0]");
        println!("  max-length              Maximum product length including primers [10000]");
        println!("  min-coverage            Minimum kmer coverage [2]");
        println!("  mismatches              Maximum primer-kmer mismatches [2]");
        println!("  trim                    Bases to keep at 3' end of each primer [15]");
        println!("  dedup-edit-threshold    Levenshtein distance for deduplication [10]\n");
        println!("Multiple primer pairs can be specified by repeating the flag:");
        println!("  --pcr-primers \"...\" --pcr-primers \"...\"");
        std::process::exit(0);
    }

    // Validate that --sample is provided for processing (not needed for early-exit flags above)
    let sample = args
        .sample
        .as_ref()
        .context("--sample is required. Provide a sample name as output file prefix.")?
        .clone();

    // Validate sample name for safe use in filenames
    ensure!(
        sample
            .chars()
            .all(|c| c.is_alphanumeric() || c == '_' || c == '-' || c == '.'),
        "Sample name '{}' contains characters that are unsafe for filenames. \
         Use only alphanumeric characters, hyphens, underscores, and periods.",
        sample
    );

    info!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    debug!("{:?}", args);

    // Create the output directory if it does not exist
    let directory = format!(
        "{}/",
        args.outdir
            .to_str()
            .context("Output directory path contains invalid UTF-8")?
    );
    std::fs::create_dir_all(&directory)
        .with_context(|| format!("Failed to create output directory: {}", directory))?;

    let k = args.k;

    // Validate arguments
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
    ensure!(
        args.min_kmer_count >= 1,
        "min-kmer-count must be at least 1"
    );

    // Collect PCR primer specifications from all sources
    let mut pcr_runs: Vec<pcr::PCRParams> = Vec::new();

    // Load preconfigured panels by name
    for panel_name in args.pcr_panel.iter() {
        let pcr_params = preconfigured::get_panel(panel_name)
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Error loading panel: {}", panel_name))?;
        pcr_runs.extend(pcr_params);
    }

    // Load primer panels from YAML files
    for pcr_file in args.pcr_file.iter() {
        let path_str = pcr_file
            .to_str()
            .context("PCR panel file path contains invalid UTF-8")?;
        let pcr_params = preconfigured::load_panel_file(path_str)
            .with_context(|| format!("Error loading panel file: {}", path_str))?;
        pcr_runs.extend(pcr_params);
    }

    // Parse inline primer specifications
    for pcr_string in args.pcr_primers.iter() {
        let pcr_params = parse_pcr_primers_string(pcr_string)
            .with_context(|| format!("Error parsing primer specification: {}", pcr_string))?;
        pcr_runs.push(pcr_params);
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

    // Warn if no output will be produced
    if args.chunks == 0 && pcr_runs.is_empty() {
        warn!("No --pcr-panel/--pcr-file/--pcr-primers and --chunks is 0: only a stats file will be produced");
    }

    // Validate that --sra is not combined with input files
    if args.sra.is_some() && args.input.is_some() {
        bail!("--sra cannot be combined with input files. Use one or the other.");
    }

    // Validate input files exist before starting processing
    if let Some(input_files) = &args.input {
        for file_path in input_files.iter() {
            ensure!(
                file_path.exists(),
                "Input file does not exist: {}",
                file_path.display()
            );
            ensure!(
                file_path.is_file(),
                "Input path is not a file: {}",
                file_path.display()
            );
        }

        // Check for duplicate input files (after canonicalization)
        let mut canonical_paths: Vec<PathBuf> = Vec::new();
        for file_path in input_files.iter() {
            let canonical = file_path
                .canonicalize()
                .with_context(|| format!("Failed to resolve path: {}", file_path.display()))?;
            if canonical_paths.contains(&canonical) {
                warn!(
                    "Duplicate input file: {} (same as previous entry after path resolution)",
                    file_path.display()
                );
            }
            canonical_paths.push(canonical);
        }
    }

    // Validate pcr-file paths exist
    for pcr_file in args.pcr_file.iter() {
        ensure!(
            pcr_file.exists(),
            "PCR panel file does not exist: {}",
            pcr_file.display()
        );
    }

    // Dry-run: print what would happen and exit
    if args.dry_run {
        eprintln!("sharkmer {} (dry run)", env!("CARGO_PKG_VERSION"));
        eprintln!();

        // Input sources
        eprintln!("Input:");
        if let Some(accession) = &args.sra {
            eprintln!("  SRA accession: {}", accession);
        } else if let Some(input_files) = &args.input {
            for f in input_files {
                eprintln!("  {}", f.display());
            }
        } else {
            eprintln!("  stdin");
        }

        eprintln!();
        eprintln!("Configuration:");
        eprintln!("  Sample:         {}", sample);
        eprintln!("  Output dir:     {}", directory);
        eprintln!("  Kmer length:    {}", k);
        eprintln!("  Chunks:         {}", args.chunks);
        eprintln!("  Threads:        {}", args.threads);
        eprintln!("  Min kmer count: {}", args.min_kmer_count);
        if let Some(max) = args.max_reads {
            eprintln!("  Max reads:      {}", max);
        }

        eprintln!();
        eprintln!("Output files:");
        eprintln!("  {}{}.stats.yaml", directory, sample);
        if args.chunks > 0 {
            eprintln!("  {}{}.histo", directory, sample);
            eprintln!("  {}{}.final.histo", directory, sample);
        }
        for pcr_params in &pcr_runs {
            eprintln!("  {}{}_{}.fasta", directory, sample, pcr_params.gene_name);
        }

        if !pcr_runs.is_empty() {
            eprintln!();
            eprintln!(
                "PCR primers ({} gene{}):",
                pcr_runs.len(),
                if pcr_runs.len() == 1 { "" } else { "s" }
            );
            for pcr_params in &pcr_runs {
                eprintln!(
                    "  {} (fwd: {}...{}, rev: {}...{}, len: {}-{})",
                    pcr_params.gene_name,
                    &pcr_params.forward_seq[..pcr_params.forward_seq.len().min(8)],
                    &pcr_params.forward_seq[pcr_params.forward_seq.len().saturating_sub(4)..],
                    &pcr_params.reverse_seq[..pcr_params.reverse_seq.len().min(8)],
                    &pcr_params.reverse_seq[pcr_params.reverse_seq.len().saturating_sub(4)..],
                    pcr_params.min_length,
                    pcr_params.max_length,
                );
            }
        }

        std::process::exit(0);
    }

    // Set the number of threads for Rayon to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // Ingest the fastq data
    let start = std::time::Instant::now();
    info!("Ingesting reads...");

    // When chunks == 0, allocate 1 internal chunk for kmer storage (needed for sPCR)
    let n_chunks = if args.chunks == 0 { 1 } else { args.chunks };
    let mut chunks: Vec<kmer::Chunk> = Vec::new();
    for _ in 0..n_chunks {
        chunks.push(Chunk::new(&k));
    }

    let mut state = FastqReadState {
        chunks,
        chunk_index: 0,
        reads: Vec::new(),
        n_reads_read: 0,
        n_bases_read: 0,
    };

    let max_reads = args.max_reads.unwrap_or(0);

    // Create progress indicator: bar with ETA if max_reads is known, spinner otherwise.
    // Only show at info verbosity or higher.
    let show_progress = args.verbose >= 1 && std::io::stderr().is_terminal();
    let progress = if show_progress && max_reads > 0 {
        let pb = ProgressBar::new(max_reads);
        pb.set_style(
            ProgressStyle::with_template(
                "Ingesting reads {bar:30} {human_pos}/{human_len} [{per_sec}] [ETA {eta}]",
            )
            .expect("valid progress template"),
        );
        pb
    } else if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::with_template("Ingesting reads {spinner} {human_pos} [{per_sec}]")
                .expect("valid progress template"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    if let Some(accession) = &args.sra {
        // Stream reads from ENA
        let urls = get_ena_fastq_urls(accession)?;
        for url in &urls {
            info!("Streaming from {}...", url);
            let response = ureq::get(url)
                .call()
                .with_context(|| format!("Failed to download {}", url))?;
            let reader =
                std::io::BufReader::new(flate2::read::GzDecoder::new(response.into_reader()));
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                url,
                &progress,
            )?;
            if reached_max {
                break;
            }
        }
    } else if let Some(input_files) = &args.input {
        // Read from one or more files
        for file_path in input_files.iter() {
            let file = std::fs::File::open(file_path)
                .with_context(|| format!("Failed to open file: {}", file_path.display()))?;

            let file_name = file_path.to_string_lossy();
            let has_gz_ext = file_name.ends_with(".gz") || file_name.ends_with(".gzip");

            // Detect gzip magic bytes for files without .gz extension
            let use_gzip = if has_gz_ext {
                true
            } else {
                let mut buf_reader = std::io::BufReader::new(file);
                let magic = buf_reader.fill_buf().context("Failed to peek at file")?;
                let is_gzip = magic.len() >= 2 && magic[0] == 0x1f && magic[1] == 0x8b;
                if is_gzip {
                    warn!(
                        "File '{}' appears to be gzipped (magic bytes detected) but lacks a .gz extension. Reading as gzipped.",
                        file_path.display()
                    );
                }
                // We need to re-open the file since BufReader consumed ownership
                drop(buf_reader);
                is_gzip
            };

            let file = std::fs::File::open(file_path)
                .with_context(|| format!("Failed to open file: {}", file_path.display()))?;
            let reader: Box<dyn BufRead> = if use_gzip {
                Box::new(std::io::BufReader::new(flate2::read::GzDecoder::new(file)))
            } else {
                Box::new(std::io::BufReader::new(file))
            };
            let reached_max = read_fastq(
                reader,
                &mut state,
                max_reads,
                args.validate_every,
                &file_name,
                &progress,
            )?;
            if reached_max {
                break;
            }
        }
    } else {
        // Read from stdin
        let stdin = std::io::stdin();
        ensure!(
            !stdin.is_terminal(),
            "No input files specified and stdin is a terminal.\n\
             Provide FASTQ files as arguments, use --sra, or pipe data via stdin.\n\
             Example: sharkmer -s sample -k 21 reads.fastq\n\
             Example: sharkmer -s sample --sra SRR5324768\n\
             Example: zcat reads.fastq.gz | sharkmer -s sample -k 21"
        );
        let handle = stdin.lock();

        read_fastq(
            handle,
            &mut state,
            max_reads,
            args.validate_every,
            "stdin",
            &progress,
        )?;
    }

    progress.finish_and_clear();

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

    // Warn if very few reads for sPCR
    if !pcr_runs.is_empty() && state.n_reads_read < 10_000 {
        warn!(
            "Only {} reads ingested. sPCR typically needs many more reads to produce results.",
            state.n_reads_read
        );
    }
    info!("Time to ingest reads: {:?}", start.elapsed());

    if n_reads_ingested == 0 {
        warn!("No reads were ingested. All output will be empty.");
    }

    // Consolidate chunks and create histograms
    info!("Consolidating chunks...");
    let start = std::time::Instant::now();

    let mut kmer_counts: KmerCounts = KmerCounts::new(&k);
    let mut n_singleton_kmers: Option<u64> = None;

    let histo_comment = format!(
        "# sharkmer {} k={} chunks={}",
        env!("CARGO_PKG_VERSION"),
        args.k,
        args.chunks
    );

    if args.chunks > 0 {
        // Incremental histogram mode: build histograms as chunks are merged
        let mut histos: Vec<kmer::Histogram> = Vec::with_capacity(args.chunks);

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

        // Write incremental histograms with header
        info!("Writing histograms to file...");
        let histo_path = format!("{}{}.histo", directory, sample);
        warn_if_exists(&histo_path);
        let mut file =
            std::fs::File::create(&histo_path).context("Failed to create histogram file")?;

        // Comment line with version and parameters
        writeln!(file, "{}", histo_comment).context("Failed to write histogram comment")?;

        // Header row
        let mut header = "count".to_string();
        for i in 1..=histos.len() {
            header = format!("{}\tchunk_{}", header, i);
        }
        writeln!(file, "{}", header).context("Failed to write histogram header")?;

        // Data rows
        for i in 1..args.histo_max as usize + 2 {
            let mut line = format!("{}", i);
            for histo in histos.iter() {
                let histo_vec = kmer::Histogram::get_vector(histo)?;
                line = format!("{}\t{}", line, histo_vec[i]);
            }
            writeln!(file, "{}", line).context("Failed to write histogram data")?;
        }

        // Write the final histogram with header
        info!("Writing final histogram to file...");
        let last_histo = &histos[histos.len() - 1];
        let last_histo_vec = kmer::Histogram::get_vector(last_histo)?;

        let final_histo_path = format!("{}{}.final.histo", directory, sample);
        warn_if_exists(&final_histo_path);
        let mut file = std::fs::File::create(&final_histo_path)
            .context("Failed to create final histogram file")?;

        writeln!(file, "{}", histo_comment).context("Failed to write final histogram comment")?;
        writeln!(file, "count\tfrequency").context("Failed to write final histogram header")?;

        for (i, value) in last_histo_vec
            .iter()
            .enumerate()
            .skip(1)
            .take(args.histo_max as usize + 1)
        {
            writeln!(file, "{}\t{}", i, value).context("Failed to write final histogram data")?;
        }

        n_singleton_kmers = Some(last_histo_vec[1]);

        // Warn if singleton rate is very high (>95% of unique kmers)
        let n_unique = last_histo.get_n_unique_kmers();
        if n_unique > 0 {
            let singleton_rate = last_histo_vec[1] as f64 / n_unique as f64;
            if singleton_rate > 0.95 {
                warn!(
                    "Very high singleton rate ({:.1}%). This may indicate very low coverage \
                     or contamination. sPCR results may be unreliable.",
                    singleton_rate * 100.0
                );
            }
        }

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
    } else {
        // No histogram mode: merge all chunks into a single kmer count table
        for chunk in state.chunks.drain(..) {
            kmer_counts.extend(chunk.get_kmer_counts())?;
            drop(chunk);
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
    }

    // Run sPCR and collect results (parallelized across genes)
    let mut pcr_results: Vec<PcrGeneResult> = Vec::new();

    if !pcr_runs.is_empty() {
        info!("Running in silico PCR...");

        // Prep kmer counts for in silico PCR. Filter low-count kmers and add reverse complements
        info!(
            "Filtering kmers with count < {} before PCR",
            args.min_kmer_count
        );
        let kmer_counts_pcr = kmer_counts.get_pcr_kmers(&args.min_kmer_count);

        // Run PCR for each gene in parallel; kmer_counts_pcr is read-only and shared
        let pcr_fasta_results: Vec<_> = pcr_runs
            .par_iter()
            .map(|pcr_params| {
                let fasta = pcr::do_pcr(&kmer_counts_pcr, &sample, pcr_params);
                (pcr_params, fasta)
            })
            .collect();

        // Write output files sequentially to maintain deterministic order
        for (pcr_params, fasta_result) in pcr_fasta_results {
            let fasta = fasta_result?;

            if !fasta.is_empty() {
                let fasta_path = format!("{}{}_{}.fasta", directory, sample, pcr_params.gene_name);
                warn_if_exists(&fasta_path);
                let mut file = std::fs::File::create(&fasta_path)
                    .with_context(|| format!("Failed to create FASTA file: {}", fasta_path))?;

                let product_lengths: Vec<usize> = fasta.iter().map(|r| r.seq().len()).collect();

                for record in fasta.iter() {
                    write_fasta_record(&mut file, record)
                        .context("Failed to write FASTA record")?;
                }

                pcr_results.push(PcrGeneResult {
                    gene_name: pcr_params.gene_name.clone(),
                    status: "success".to_string(),
                    n_products: product_lengths.len(),
                    product_lengths,
                    output_file: Some(fasta_path),
                });
            } else {
                pcr_results.push(PcrGeneResult {
                    gene_name: pcr_params.gene_name.clone(),
                    status: "fail".to_string(),
                    n_products: 0,
                    product_lengths: Vec::new(),
                    output_file: None,
                });
            }
        }

        info!("Done running in silico PCR");
    }

    // Reconstruct command line string
    let command = std::env::args().collect::<Vec<String>>().join(" ");

    // Write stats as YAML
    info!("Writing stats to file...");
    let run_stats = RunStats {
        sharkmer_version: env!("CARGO_PKG_VERSION").to_string(),
        command,
        sample: sample.clone(),
        kmer_length: args.k,
        chunks: args.chunks,
        n_reads_read: state.n_reads_read,
        n_bases_read: state.n_bases_read,
        n_subreads_ingested: n_reads_ingested,
        n_bases_ingested,
        n_kmers: n_kmers_ingested,
        n_multi_kmers: n_singleton_kmers.map(|s| n_kmers_ingested.saturating_sub(s)),
        n_singleton_kmers,
        peak_memory_bytes: PEAK_ALLOC.peak_usage() as u64,
        pcr_results,
    };

    let stats_path = format!("{}{}.stats.yaml", directory, sample);
    warn_if_exists(&stats_path);
    let file_stats = std::fs::File::create(&stats_path).context("Failed to create stats file")?;
    serde_yaml::to_writer(file_stats, &run_stats).context("Failed to write stats YAML")?;

    // Print summary line (warn level = always visible unless --quiet)
    let elapsed = start_run.elapsed();
    let elapsed_secs = elapsed.as_secs_f64();
    let elapsed_str = if elapsed_secs < 60.0 {
        format!("{:.1}s", elapsed_secs)
    } else if elapsed_secs < 3600.0 {
        let mins = (elapsed_secs / 60.0).floor() as u64;
        let secs = elapsed_secs - (mins as f64 * 60.0);
        format!("{}m {:.0}s", mins, secs)
    } else {
        let hours = (elapsed_secs / 3600.0).floor() as u64;
        let remaining = elapsed_secs - (hours as f64 * 3600.0);
        let mins = (remaining / 60.0).floor() as u64;
        format!("{}h {}m", hours, mins)
    };

    let reads_str = format_count(run_stats.n_reads_read);

    if !run_stats.pcr_results.is_empty() {
        let n_success = run_stats
            .pcr_results
            .iter()
            .filter(|r| r.status == "success")
            .count();
        let n_total = run_stats.pcr_results.len();
        let gene_names: Vec<&str> = run_stats
            .pcr_results
            .iter()
            .filter(|r| r.status == "success")
            .map(|r| r.gene_name.as_str())
            .collect();
        let genes_detail = if gene_names.is_empty() {
            String::new()
        } else {
            format!(" ({})", gene_names.join(", "))
        };
        warn!(
            "sharkmer complete: {} reads, {}/{} genes amplified{}, peak mem {}, {}",
            reads_str,
            n_success,
            n_total,
            genes_detail,
            format_bytes(run_stats.peak_memory_bytes),
            elapsed_str
        );
    } else {
        let kmers_str = format_count(run_stats.n_kmers);
        warn!(
            "sharkmer complete: {} reads, {} kmers, {} chunks, peak mem {}, {}",
            reads_str,
            kmers_str,
            run_stats.chunks,
            format_bytes(run_stats.peak_memory_bytes),
            elapsed_str
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {}

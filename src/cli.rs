use anyhow::{Context, Result, bail, ensure};
use clap::{CommandFactory, Parser};
use clap_complete::Shell;
use log::warn;
use std::path::PathBuf;

use crate::io::EnaResult;
use crate::pcr;
use crate::pcr::preconfigured;

/// Parse an inline primer specification string (key=value,key=value,...) into PCRParams.
pub fn parse_pcr_primers_string(pcr_string: &str) -> Result<pcr::PCRParams> {
    ensure!(!pcr_string.is_empty(), "Invalid empty primer specification");

    let split: Vec<&str> = pcr_string.split(',').collect();

    let mut forward_seq = String::new();
    let mut reverse_seq = String::new();
    let mut gene_name = String::new();
    let mut max_length = 10000;
    let mut min_length = 0;
    let mut min_count = 2;
    let mut mismatches = 2;
    let mut trim = 15;
    let mut citation = String::new();
    let mut notes = String::new();
    let mut dedup_edit_threshold = pcr::DEFAULT_DEDUP_EDIT_THRESHOLD;

    for item in split.iter() {
        let key_value: Vec<&str> = item.split('=').collect();
        if key_value.len() != 2 {
            bail!(
                "Invalid parameter (should be key=value): '{}'\n\
                 Commas are not allowed in field values. \
                 Use --pcr-panel-file with a YAML panel for complex metadata.",
                item
            );
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
            "min-count" => {
                min_count = value
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
        min_count,
        mismatches,
        trim,
        citation,
        notes,
        dedup_edit_threshold,
        source: format!("--pcr-primers \"{}\"", pcr_string),
    };

    Ok(pcr_params)
}

/// A tool for kmer counting and in silico PCR (sPCR)
#[derive(Parser, Debug)]
#[command(
    author,
    version = concat!(env!("CARGO_PKG_VERSION"), " (", env!("CARGO_PKG_REPOSITORY"), ")"),
    about,
    long_about = None,
    arg_required_else_help = true
)]
#[command(after_help = "\
Example:\n  \
  Extract cnidarian genes from ENA reads (downloads automatically):\n  \
  sharkmer --ena SRR23143286 --pcr-panel cnidaria -m 1000000 -o output\n\
\n\
Output files:\n  \
  {outdir}/{sample}.stats.yaml             Run statistics (always produced)\n\
\n  \
  PCR:\n  \
  {outdir}/{sample}_{panel}_{gene}.fasta   sPCR products per gene\n\
\n  \
  Incremental counting (--chunks > 0):\n  \
  {outdir}/{sample}.histo                  All incremental histograms\n  \
  {outdir}/{sample}.final.histo            Final histogram")]
pub(crate) struct Args {
    /// FASTQ input files (.fastq or .fastq.gz). Reads from stdin if omitted
    #[arg(help_heading = "Input", display_order = 0)]
    pub(crate) input: Option<Vec<PathBuf>>,

    /// Stream reads directly from ENA by accession (e.g. SRR5324768)
    #[arg(long, help_heading = "Input", display_order = 1)]
    pub(crate) ena: Option<String>,

    /// Sample name (output file prefix; required unless --ena derives it)
    #[arg(short, long, help_heading = "Output")]
    pub(crate) sample: Option<String>,

    /// Output directory
    #[arg(short, long, default_value = "./", help_heading = "Output")]
    pub(crate) outdir: PathBuf,

    /// Use a preconfigured primer panel (repeatable)
    #[arg(long, help_heading = "PCR")]
    pub(crate) pcr_panel: Vec<String>,

    /// Load a primer panel from a YAML file or URL (repeatable)
    #[arg(long, help_heading = "PCR")]
    pub(crate) pcr_panel_file: Vec<String>,

    /// Specify a primer pair inline (repeatable, see --help-pcr)
    #[arg(long, help_heading = "PCR")]
    pub(crate) pcr_primers: Vec<String>,

    /// List available primer panels and exit
    #[arg(long, help_heading = "PCR info")]
    pub(crate) list_panels: bool,

    /// Export a built-in panel as YAML to stdout and exit
    #[arg(long, help_heading = "PCR info")]
    pub(crate) export_panel: Option<String>,

    /// Show detailed help for --pcr-primers format
    #[arg(long, help_heading = "PCR info")]
    pub(crate) help_pcr: bool,

    /// Kmer length
    #[arg(short, default_value_t = 21, help_heading = "Kmer")]
    pub(crate) k: usize,

    /// Number of incremental chunks (0 = skip histograms)
    #[arg(long, default_value_t = 0, help_heading = "Incremental counting")]
    pub(crate) chunks: usize,

    /// Maximum histogram count value
    #[arg(long, default_value_t = 10000, help_heading = "Incremental counting")]
    pub(crate) histo_max: u64,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 1, help_heading = "General")]
    pub(crate) threads: usize,

    /// Maximum number of reads to process (default: all)
    #[arg(short = 'm', long, help_heading = "General")]
    pub(crate) max_reads: Option<u64>,

    /// Minimum kmer count for sPCR (filters low-count kmers before PCR)
    #[arg(long, default_value_t = 2, help_heading = "General")]
    pub(crate) min_kmer_count: u32,

    /// Validate FASTQ format every N records (0 = first record only)
    #[arg(long, default_value_t = 0, help_heading = "General")]
    pub(crate) validate_every: u64,

    /// Increase verbosity (-v info, -vv debug, -vvv trace)
    #[arg(short = 'v', long, action = clap::ArgAction::Count, help_heading = "General")]
    pub(crate) verbose: u8,

    /// Suppress all output except errors
    #[arg(short = 'q', long, help_heading = "General")]
    pub(crate) quiet: bool,

    /// Color output
    #[arg(long, default_value = "auto", help_heading = "General")]
    pub(crate) color: ColorMode,

    /// Print citation information and exit
    #[arg(long, help_heading = "General")]
    pub(crate) cite: bool,

    /// Print shell tab-completion script and exit
    #[arg(long, help_heading = "General", value_name = "SHELL")]
    pub(crate) completions: Option<Shell>,

    /// Write assembly graphs as annotated DOT (Graphviz) files
    #[arg(long, help_heading = "PCR info")]
    pub(crate) dump_graph: bool,

    /// Validate primer panels/primers and exit
    #[arg(long, help_heading = "PCR info")]
    pub(crate) validate_panels: bool,

    /// Validate inputs and print what would happen, then exit
    #[arg(long, help_heading = "General")]
    pub(crate) dry_run: bool,

    /// Override cache directory for remote reads
    #[arg(long, help_heading = "Cache")]
    pub(crate) cache_dir: Option<PathBuf>,

    /// Disable read caching (stream directly, do not read from or write to cache)
    #[arg(long, help_heading = "Cache")]
    pub(crate) no_cache: bool,

    /// Delete the read cache directory and exit
    #[arg(long, help_heading = "Cache")]
    pub(crate) clear_cache: bool,
}

/// Color output mode for log messages.
#[derive(Debug, Clone, clap::ValueEnum)]
pub(crate) enum ColorMode {
    Auto,
    Always,
    Never,
}

/// Initialize env_logger with the given verbosity level, quiet flag, and color mode.
pub(crate) fn init_logging(verbose: u8, quiet: bool, color: &ColorMode) {
    let log_level = if quiet {
        log::LevelFilter::Error
    } else {
        match verbose {
            0 => log::LevelFilter::Warn,
            1 => log::LevelFilter::Info,
            2 => log::LevelFilter::Debug,
            _ => log::LevelFilter::Trace,
        }
    };
    let write_style = match color {
        ColorMode::Auto => env_logger::WriteStyle::Auto,
        ColorMode::Always => env_logger::WriteStyle::Always,
        ColorMode::Never => env_logger::WriteStyle::Never,
    };
    env_logger::Builder::new()
        .filter_level(log_level)
        .write_style(write_style)
        .format_timestamp(None)
        .format(|buf, record| {
            use env_logger::fmt::style::{AnsiColor, Style};
            use std::io::Write as _;
            let level = record.level();
            let style = match level {
                log::Level::Error => Style::new().fg_color(Some(AnsiColor::Red.into())).bold(),
                log::Level::Warn => Style::new().fg_color(Some(AnsiColor::Green.into())).bold(),
                log::Level::Info => Style::new().fg_color(Some(AnsiColor::Cyan.into())).bold(),
                log::Level::Debug => Style::new().dimmed(),
                log::Level::Trace => Style::new().dimmed(),
            };
            let label = match level {
                log::Level::Error => "error",
                log::Level::Warn => "",
                log::Level::Info => "info",
                log::Level::Debug => "debug",
                log::Level::Trace => "trace",
            };
            if label.is_empty() {
                // Warn-level: no prefix, just the message (for status output)
                writeln!(buf, "{}", record.args())
            } else {
                writeln!(
                    buf,
                    "{}{}{} {}",
                    style.render(),
                    label,
                    style.render_reset(),
                    record.args()
                )
            }
        })
        .init();
}

/// Handle early-exit flags (--completions, --cite, --list-panels, --export-panel, --help-pcr).
/// Each prints output and calls process::exit(0).
pub(crate) fn handle_early_exits(args: &Args) -> Result<()> {
    // Clear read cache and exit
    if args.clear_cache {
        crate::cache::CacheConfig::clear(args.cache_dir.as_deref())?;
        println!("Cache cleared.");
        std::process::exit(0);
    }

    // Generate shell completions and exit
    if let Some(shell) = args.completions {
        let mut cmd = Args::command();
        clap_complete::generate(shell, &mut cmd, "sharkmer", &mut std::io::stdout());
        std::process::exit(0);
    }

    // Print citation information and exit if --cite is specified
    if args.cite {
        println!("{} {}\n", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
        println!("A paper describing sharkmer is in preparation. In the meantime,");
        println!("please cite the following papers:\n");
        println!("For in silico PCR:\n");
        println!("  Church et al. (2025) Population genomics of a sailing siphonophore");
        println!("  reveals genetic structure in the open ocean.");
        println!("  Current Biology, 35(15), 3556-3569.e6.");
        println!("  doi: 10.1016/j.cub.2025.05.066\n");
        println!("  @article{{church2025physalia,");
        println!("    title={{Population genomics of a sailing siphonophore reveals");
        println!("           genetic structure in the open ocean}},");
        println!("    author={{Church, Samuel H. and Abedon, River B. and Ahuja, Namrata and");
        println!("            Anthony, Colin J. and Destanovi\\'{{c}}, Dalila and others}},");
        println!("    journal={{Current Biology}},");
        println!("    volume={{35}},");
        println!("    number={{15}},");
        println!("    pages={{3556--3569.e6}},");
        println!("    year={{2025}},");
        println!("    doi={{10.1016/j.cub.2025.05.066}}");
        println!("  }}\n");
        println!("For incremental kmer counting:\n");
        println!("  Ahuja et al. (2024) Giants among Cnidaria: Large Nuclear Genomes and");
        println!("  Rearranged Mitochondrial Genomes in Siphonophores.");
        println!("  Genome Biology and Evolution, 16(3).");
        println!("  doi: 10.1093/gbe/evae048\n");
        println!("  @article{{ahuja2024siphonophore,");
        println!("    title={{Giants among Cnidaria: Large Nuclear Genomes and Rearranged");
        println!("           Mitochondrial Genomes in Siphonophores}},");
        println!("    author={{Ahuja, N. and Cao, X. and Schultz, D. T. and Picciani, N. and");
        println!("            Lord, A. and Shao, S. and Jia, K. and Burdick, D. R. and");
        println!("            Haddock, S. H. D. and Li, Y. and Dunn, C. W.}},");
        println!("    journal={{Genome Biology and Evolution}},");
        println!("    volume={{16}},");
        println!("    number={{3}},");
        println!("    year={{2024}},");
        println!("    doi={{10.1093/gbe/evae048}}");
        println!("  }}");
        std::process::exit(0);
    }

    // List available panels and exit
    if args.list_panels {
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
        println!(
            "  --pcr-primers \"forward=GRCTGTTTACCAAAAACATA,reverse=AATTCAACATMGAGG,max-length=700,name=16s,min-length=500\"\n"
        );
        println!("Required keys:");
        println!("  forward       Forward primer sequence (5' to 3')");
        println!("  reverse       Reverse primer sequence (5' to 3' on opposite strand)");
        println!("  name          Unique name for the primer pair or gene region\n");
        println!("Optional keys:");
        println!("  min-length              Minimum product length including primers [0]");
        println!("  max-length              Maximum product length including primers [10000]");
        println!("  min-count               Minimum kmer count for graph extension [2]");
        println!("  mismatches              Maximum primer-kmer mismatches [2]");
        println!("  trim                    Bases to keep at 3' end of each primer [15]");
        println!("  dedup-edit-threshold    Levenshtein distance for deduplication [10]\n");
        println!("Primer sequences support IUPAC ambiguity codes:");
        println!("  R (A/G)  Y (C/T)  S (G/C)  W (A/T)  K (G/T)  M (A/C)");
        println!("  B (C/G/T)  D (A/G/T)  H (A/C/T)  V (A/C/G)  N (A/C/G/T)\n");
        println!("Multiple primer pairs can be specified by repeating the flag:");
        println!("  --pcr-primers \"...\" --pcr-primers \"...\"\n");
        println!("Note: when using --pcr-panel or --pcr-panel-file, gene names in output");
        println!("files are prefixed with the panel name (e.g., cnidaria_18S).");
        println!("Inline --pcr-primers gene names are used as-is.");
        std::process::exit(0);
    }

    Ok(())
}

/// Collect PCR primer specifications from all sources (--pcr-panel, --pcr-panel-file,
/// --pcr-primers), validate them, check for duplicates, and clamp min-count.
pub(crate) fn collect_pcr_params(args: &Args) -> Result<Vec<pcr::PCRParams>> {
    let mut pcr_runs: Vec<pcr::PCRParams> = Vec::new();

    // Load preconfigured panels by name
    for panel_name in args.pcr_panel.iter() {
        let mut pcr_params = preconfigured::get_panel(panel_name)
            .map_err(|e| anyhow::anyhow!(e))
            .with_context(|| format!("Failed to load panel: {}", panel_name))?;
        for p in pcr_params.iter_mut() {
            p.source = format!("built-in panel '{}'", panel_name);
        }
        pcr_runs.extend(pcr_params);
    }

    // Load primer panels from YAML files or URLs
    for panel_source in args.pcr_panel_file.iter() {
        let mut pcr_params = preconfigured::load_panel_source(panel_source)
            .with_context(|| format!("Failed to load panel: {}", panel_source))?;
        let source_label = if preconfigured::is_url(panel_source) {
            format!("panel URL '{}'", panel_source)
        } else {
            format!("panel file '{}'", panel_source)
        };
        for p in pcr_params.iter_mut() {
            p.source = source_label.clone();
        }
        pcr_runs.extend(pcr_params);
    }

    // Parse inline primer specifications
    for pcr_string in args.pcr_primers.iter() {
        let pcr_params = parse_pcr_primers_string(pcr_string)
            .with_context(|| format!("Could not parse primer specification: \"{}\"", pcr_string))?;
        pcr_runs.push(pcr_params);
    }

    // Validate all primer parameters, collecting all errors
    let mut total_errors = 0usize;
    let mut error_report = String::new();
    for pcr_params in pcr_runs.iter() {
        let errors = pcr::validate_pcr_params(pcr_params);
        if !errors.is_empty() {
            total_errors += errors.len();
            error_report.push_str(&format!(
                "\n  {} ({}):\n",
                pcr_params.gene_name, pcr_params.source
            ));
            for (err, suggestion) in &errors {
                error_report.push_str(&format!(
                    "    - {}\n      Suggestion: {}\n",
                    err, suggestion
                ));
            }
        }
    }
    if total_errors > 0 {
        bail!(
            "Primer validation failed ({} error{}):\n{}",
            total_errors,
            if total_errors == 1 { "" } else { "s" },
            error_report
        );
    }

    // Warn and clamp if any primer's min-count < --min-kmer-count
    for pcr_params in pcr_runs.iter_mut() {
        if pcr_params.min_count < args.min_kmer_count {
            warn!(
                "{}: min-count ({}) is less than --min-kmer-count ({}). \
                 Kmers below {} have already been filtered. Using {} as effective min-count.",
                pcr_params.gene_name,
                pcr_params.min_count,
                args.min_kmer_count,
                args.min_kmer_count,
                args.min_kmer_count
            );
            pcr_params.min_count = args.min_kmer_count;
        }
    }

    // Check that there are no duplicate gene names
    let mut gene_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    for pcr_params in pcr_runs.iter() {
        ensure!(
            gene_names.insert(pcr_params.gene_name.clone()),
            "Duplicate gene name '{}' (from {})",
            pcr_params.gene_name,
            pcr_params.source
        );
    }

    Ok(pcr_runs)
}

/// Print validation output for all collected primer pairs and exit.
pub(crate) fn handle_validate_panels(pcr_runs: &[pcr::PCRParams]) -> Result<()> {
    if pcr_runs.is_empty() {
        bail!(
            "--validate-panels requires at least one of --pcr-panel, --pcr-panel-file, or --pcr-primers"
        );
    }
    println!("Validated {} primer pairs:\n", pcr_runs.len());
    for p in pcr_runs {
        println!("  {}", p.gene_name);
        println!(
            "    forward:  {} ({} bp)",
            p.forward_seq,
            p.forward_seq.len()
        );
        println!(
            "    reverse:  {} ({} bp)",
            p.reverse_seq,
            p.reverse_seq.len()
        );
        println!("    length:   {}-{} bp", p.min_length, p.max_length);
        println!("    min-count: >= {}", p.min_count);
        println!(
            "    mismatches: {}, trim: {}, dedup-edit-threshold: {}",
            p.mismatches, p.trim, p.dedup_edit_threshold
        );
    }
    println!("\nAll primers valid.");
    std::process::exit(0);
}

/// Derive the sample name from --sample or --ena ENA metadata.
/// Returns (sample_name, optional_cached_ena_result).
pub(crate) fn resolve_sample_name(args: &Args) -> Result<(String, Option<EnaResult>)> {
    let mut cached_ena_result: Option<EnaResult> = None;
    let sample = if let Some(ref s) = args.sample {
        s.clone()
    } else if let Some(ref accession) = args.ena {
        let ena_result = crate::io::get_ena_fastq_urls(accession)?;
        let derived = if let Some(ref sci_name) = ena_result.scientific_name {
            let genus_species = sci_name.replace(' ', "_");
            format!("{}_{}", genus_species, accession)
        } else {
            accession.clone()
        };
        warn!(
            "No --sample provided, using '{}' derived from ENA metadata",
            derived
        );
        cached_ena_result = Some(ena_result);
        derived
    } else {
        anyhow::bail!(
            "--sample is required. Provide a sample name as output file prefix.\n\
             When using --ena, the sample name can be derived automatically from ENA metadata."
        );
    };

    // Validate sample name for safe use in filenames
    ensure!(
        sample
            .chars()
            .all(|c| c.is_alphanumeric() || c == '_' || c == '-' || c == '.'),
        "Sample name '{}' contains characters that are unsafe for filenames. \
         Use only alphanumeric characters, hyphens, underscores, and periods.",
        sample
    );

    Ok((sample, cached_ena_result))
}

/// Validate CLI arguments: k, histo_max, min_kmer_count, input file existence,
/// --ena vs input conflicts, and panel file paths.
pub(crate) fn validate_args(args: &Args, pcr_runs: &[pcr::PCRParams]) -> Result<()> {
    let k = args.k;

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

    // Validate that --ena is not combined with input files
    if args.ena.is_some() && args.input.is_some() {
        bail!("--ena cannot be combined with input files. Use one or the other.");
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
        let mut canonical_paths: std::collections::HashSet<PathBuf> =
            std::collections::HashSet::new();
        for file_path in input_files.iter() {
            let canonical = file_path
                .canonicalize()
                .with_context(|| format!("Failed to resolve path: {}", file_path.display()))?;
            if !canonical_paths.insert(canonical) {
                warn!(
                    "Duplicate input file: {} (same as previous entry after path resolution)",
                    file_path.display()
                );
            }
        }
    }

    // Validate pcr-panel-file paths exist (skip URL sources — validated at download time)
    for panel_source in args.pcr_panel_file.iter() {
        if !preconfigured::is_url(panel_source) {
            let path = std::path::Path::new(panel_source);
            ensure!(
                path.exists(),
                "PCR panel file does not exist: {}",
                panel_source
            );
        }
    }

    // Warn if no output will be produced (after all validation passes)
    if args.chunks == 0 && pcr_runs.is_empty() {
        warn!(
            "No --pcr-panel/--pcr-panel-file/--pcr-primers and --chunks is 0: only a stats file will be produced"
        );
    }

    Ok(())
}

/// Print dry-run information (what would happen) and exit.
pub(crate) fn handle_dry_run(
    args: &Args,
    sample: &str,
    directory: &str,
    pcr_runs: &[pcr::PCRParams],
) {
    let k = args.k;

    eprintln!("sharkmer {} (dry run)", env!("CARGO_PKG_VERSION"));
    eprintln!();

    // Input sources
    eprintln!("Input:");
    if let Some(accession) = &args.ena {
        eprintln!("  ENA accession: {}", accession);
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
    for pcr_params in pcr_runs {
        eprintln!("  {}{}_{}.fasta", directory, sample, pcr_params.gene_name);
    }

    if !pcr_runs.is_empty() {
        eprintln!();
        eprintln!(
            "PCR primers ({} gene{}):",
            pcr_runs.len(),
            if pcr_runs.len() == 1 { "" } else { "s" }
        );
        for pcr_params in pcr_runs {
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

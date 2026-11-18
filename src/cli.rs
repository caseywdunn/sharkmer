use anyhow::{Context, Result, bail, ensure};
use clap::{CommandFactory, Parser};
use clap_complete::Shell;
use log::warn;
use std::path::PathBuf;

use crate::configuration::{Owned, Owner};
use crate::io::EnaResult;
use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use crate::pcr;
use crate::pcr::preconfigured;
use std::sync::Arc;

/// Parse an inline primer specification string (key=value,key=value,...) into PCRParams.
fn parse_pcr_primers_owned(
    pcr_string: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Owned<pcr::PCRParams>> {
    let mut owner = Owner::new(budget)?;
    ensure!(!pcr_string.is_empty(), "Invalid empty primer specification");

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
    let mut seen_keys = BudgetMap::in_class(budget, MemoryClass::RetainedResults);

    for item in pcr_string.split(',') {
        // split_once(|c| == '=') keeps everything after the first '=' as the
        // value. Values may legitimately contain '=' — citations include
        // URLs, notes may include equations — so we do not reject them.
        // Commas are still forbidden because they are the field separator
        // at the outer level.
        let (raw_key, value) = match item.split_once('=') {
            Some(kv) => kv,
            None => bail!(
                "Invalid parameter (should be key=value): '{}'\n\
                 Commas are not allowed in field values. \
                 Use --pcr-panel-file with a YAML panel for complex metadata.",
                item
            ),
        };

        let key_owned = owner.case(raw_key, false)?;
        if seen_keys.get(key_owned.as_str()).is_some() {
            bail!(
                "Duplicate parameter '{}' in primer specification '{}'. Each key may appear at most once.",
                key_owned,
                pcr_string
            );
        }
        let key = key_owned.as_str();

        match key {
            "name" => {
                gene_name = owner.text(format_args!("{value}"))?;
            }
            "forward" => {
                forward_seq = owner.case(value, true)?;
            }
            "reverse" => {
                reverse_seq = owner.case(value, true)?;
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
                citation = owner.text(format_args!("{value}"))?;
            }
            "notes" => {
                notes = owner.text(format_args!("{value}"))?;
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
        seen_keys.insert(key_owned, ())?;
    }

    let pcr_params = pcr::PCRParams {
        forward_seq,
        reverse_seq,
        min_length,
        max_length,
        gene_name,
        gene: None,
        region: None,
        index: None,
        compartment: None,
        gene_type: None,
        copy_number: None,
        deprecated: false,
        deprecated_by: None,
        deprecated_reason: None,
        min_count,
        mismatches,
        trim,
        expected_length: None,
        citation,
        notes,
        dedup_edit_threshold,
        source: owner.text(format_args!("--pcr-primers \"{}\"", pcr_string))?,
        max_dfs_states: pcr::DEFAULT_MAX_DFS_STATES,
        max_paths_per_pair: pcr::DEFAULT_MAX_PATHS_PER_PAIR,
        max_node_visits: pcr::DEFAULT_MAX_NODE_VISITS,
        max_primer_kmers: pcr::DEFAULT_MAX_NUM_PRIMER_KMERS,
        high_coverage_ratio: pcr::DEFAULT_HIGH_COVERAGE_RATIO,
        tip_coverage_fraction: pcr::DEFAULT_TIP_COVERAGE_FRACTION,
        diagnose_withheld_paths: false,
        emit_uncertain: false,
    };

    drop(seen_keys);
    owner.retain_bytes(pcr_params.configuration_bytes()?)?;
    Ok(owner.finish(pcr_params))
}

#[cfg(test)]
pub(crate) fn parse_pcr_primers_string(pcr_string: &str) -> Result<pcr::PCRParams> {
    Ok(parse_pcr_primers_owned(pcr_string, None)?.into_parts().0)
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
  {outdir}/{sample}_{gene}.fasta           sPCR products per gene\n\
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
    #[arg(short, default_value_t = 19, help_heading = "Kmer")]
    pub(crate) k: usize,

    /// Number of incremental chunks (0 = skip histograms)
    #[arg(long, default_value_t = 0, help_heading = "Incremental counting")]
    pub(crate) chunks: usize,

    /// Development override; normal sPCR uses the measured exact counter.
    #[arg(long, hide = true, value_parser = ["legacy", "memory", "external", "parallel", "packed", "adaptive", "adaptive-parallel", "pipeline", "hybrid"])]
    pub(crate) count_backend: Option<String>,

    /// Verify every count pair in a separate, opt-in scan; excluded from timing.
    #[arg(long, hide = true)]
    pub(crate) count_fingerprint: bool,

    /// Experimental two-thread worker transfer at this current table size.
    #[arg(long, hide = true, default_value_t = 16)]
    pub(crate) count_owner_threshold_mib: u64,

    /// Count-table memory allowance in MiB; full tables spill to temporary storage
    #[arg(long, help_heading = "Resources", default_value_t = 8192)]
    pub(crate) count_table_budget_mib: u64,

    /// Memory allowance in MiB for sorting spilled counts
    #[arg(long, help_heading = "Resources", default_value_t = 64)]
    pub(crate) count_sort_buffer_mib: u64,

    /// Shared PCR working-memory allowance in MiB (normal default: 1024; 0 disables)
    #[arg(long, help_heading = "Resources")]
    pub(crate) graph_buffer_budget_mib: Option<u64>,

    /// PCR working-memory ceiling per target in MiB
    #[arg(long, help_heading = "Resources", default_value_t = 256)]
    pub(crate) graph_buffer_target_mib: u64,

    /// Temporary-file allowance in MiB for counts, replay, staged output and recovery
    #[arg(long, help_heading = "Resources", default_value_t = 65536)]
    pub(crate) count_disk_budget_mib: u64,

    /// Directory for temporary count/replay files (default: OS temporary directory)
    #[arg(long, help_heading = "Resources")]
    pub(crate) count_temp_dir: Option<PathBuf>,

    /// Maximum four-line FASTQ record size, including line delimiters.
    #[arg(long, help_heading = "Resources", default_value_t = 16777216)]
    pub(crate) fastq_record_bytes: usize,

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

    /// Global node budget: maximum total nodes in the assembly graph.
    /// If not set, scales automatically with data volume (100K–500K).
    #[arg(long, help_heading = "PCR", hide = true)]
    pub(crate) node_budget_global: Option<usize>,

    /// Maximum DFS states to explore per start node during path finding
    #[arg(long, default_value_t = crate::pcr::DEFAULT_MAX_DFS_STATES, help_heading = "PCR", hide = true)]
    pub(crate) max_dfs_states: usize,

    /// Maximum paths to find per primer pair
    #[arg(long, default_value_t = crate::pcr::DEFAULT_MAX_PATHS_PER_PAIR, help_heading = "PCR", hide = true)]
    pub(crate) max_paths_per_pair: usize,

    /// Maximum times a node can appear in a single path (cycle tolerance)
    #[arg(long, default_value_t = crate::pcr::DEFAULT_MAX_NODE_VISITS, help_heading = "PCR", hide = true)]
    pub(crate) max_node_visits: usize,

    /// Maximum number of primer-matching kmers to retain per primer
    #[arg(long, default_value_t = crate::pcr::DEFAULT_MAX_NUM_PRIMER_KMERS, help_heading = "PCR", hide = true)]
    pub(crate) max_primer_kmers: usize,

    /// Skip edges during extension whose count exceeds this multiple of the median
    #[arg(long, default_value_t = crate::pcr::DEFAULT_HIGH_COVERAGE_RATIO, help_heading = "PCR", hide = true)]
    pub(crate) high_coverage_ratio: f64,

    /// Tip coverage fraction below which a dead-end tip is pruned
    #[arg(long, default_value_t = crate::pcr::DEFAULT_TIP_COVERAGE_FRACTION, help_heading = "PCR", hide = true)]
    pub(crate) tip_coverage_fraction: f64,

    /// Enable read threading (Pass 2): re-read FASTQ to annotate graph edges with read support.
    /// Hidden: works but loads all reads into RAM before sPCR; needs to be
    /// moved into do_pcr() for per-gene streaming before public exposure.
    #[arg(long, help_heading = "PCR", hide = true)]
    pub(crate) read_threading: bool,

    #[arg(
        long,
        help = "Retain bounded diagnostic records for repeat-withheld paths",
        help_heading = "PCR",
        hide = true
    )]
    pub(crate) diagnose_withheld_paths: bool,

    /// Include bounded sequences with assembly warnings (already enabled by default)
    #[arg(
        long,
        help_heading = "PCR",
        hide = true,
        conflicts_with = "no_emit_uncertain"
    )]
    pub(crate) emit_uncertain: bool,

    /// Exclude sequences with assembly warnings from FASTA output
    #[arg(long, help_heading = "PCR", conflicts_with = "emit_uncertain")]
    pub(crate) no_emit_uncertain: bool,

    /// Treat input as paired-end reads (exactly 2 files required: R1, R2).
    /// Hidden: ingestion works, but the paired-end phasing pipeline it is
    /// meant to unlock is not yet wired into branch ranking. The
    /// `PairedEndLink` data structure is populated during threading but
    /// no downstream consumer reads it. See caseywdunn/sharkmer#101.
    #[arg(long, help_heading = "Input", hide = true)]
    pub(crate) paired: bool,

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

impl Args {
    /// Select exact counting for ordinary runs and legacy counting for histograms.
    pub(crate) fn count_backend(&self) -> &str {
        self.count_backend
            .as_deref()
            .unwrap_or(if self.chunks == 0 { "hybrid" } else { "legacy" })
    }

    /// Explicit development backend selectors retain their independent controls.
    pub(crate) fn uses_production_defaults(&self) -> bool {
        self.count_backend.is_none() && self.chunks == 0
    }

    /// Explicit zero disables graph accounting; omission enables normal defaults.
    pub(crate) fn graph_buffer_budget_mib(&self) -> u64 {
        self.graph_buffer_budget_mib
            .unwrap_or(if self.uses_production_defaults() {
                1024
            } else {
                0
            })
    }
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
        println!("If you use sharkmer in published work, please cite:\n");
        println!("For in silico PCR:\n");
        println!("  Dunn and Church (2026) Sharkmer: repurposing PCR primers for");
        println!("  targeted genome assembly using in silico PCR.");
        println!("  Bioinformatics, btag163.");
        println!("  doi: 10.1093/bioinformatics/btag163\n");
        println!("  @article{{Dunn2026sharkmer,");
        println!("    author = {{Dunn, Casey W and Church, Samuel H}},");
        println!("    title = {{Sharkmer: repurposing PCR primers for targeted genome");
        println!("             assembly using in silico PCR}},");
        println!("    journal = {{Bioinformatics}},");
        println!("    pages = {{btag163}},");
        println!("    year = {{2026}},");
        println!("    month = {{04}},");
        println!("    issn = {{1367-4811}},");
        println!("    doi = {{10.1093/bioinformatics/btag163}}");
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
/// Flat primers and their separate metadata owners, all on the original ledger.
/// The parameter values drop before owners for parsed strings and arrays.
pub(crate) struct CollectedPrimers {
    values: BudgetVec<pcr::PCRParams>,
    owners: BudgetVec<Owner>,
}
impl std::ops::Deref for CollectedPrimers {
    type Target = [pcr::PCRParams];
    fn deref(&self) -> &Self::Target {
        &self.values
    }
}
impl std::fmt::Debug for CollectedPrimers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.values.fmt(f)
    }
}
impl CollectedPrimers {
    /// Runtime overrides change scalar tuning values, retaining metadata owners.
    pub(crate) fn iter_mut(&mut self) -> std::slice::IterMut<'_, pcr::PCRParams> {
        self.values.iter_mut()
    }
    fn append(&mut self, group: Owned<Vec<pcr::PCRParams>>) -> Result<()> {
        self.owners.reserve_one()?;
        let retained = group.iter().try_fold(0_u64, |total, value| {
            total
                .checked_add(value.configuration_bytes()?)
                .context("Primer metadata size overflow")
        })?;
        let (values, owner) = group.into_parts();
        self.owners.push(owner)?;
        for value in values {
            self.values.push(value)?;
        }
        self.owners
            .last_mut()
            .context("Missing primer metadata owner")?
            .retain_bytes(retained)?;
        Ok(())
    }
    fn push(&mut self, value: Owned<pcr::PCRParams>) -> Result<()> {
        self.owners.reserve_one()?;
        let (value, owner) = value.into_parts();
        self.owners.push(owner)?;
        self.values.push(value)
    }
}

pub(crate) fn collect_pcr_params(
    args: &Args,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<CollectedPrimers> {
    let groups = args
        .pcr_panel
        .len()
        .checked_add(args.pcr_panel_file.len())
        .and_then(|n| n.checked_add(args.pcr_primers.len()))
        .context("Primer owner capacity overflow")?;
    let mut pcr_runs = CollectedPrimers {
        values: BudgetVec::in_class(budget, MemoryClass::RetainedResults),
        owners: BudgetVec::with_capacity_in_class(groups, budget, MemoryClass::RetainedResults)?,
    };
    for panel_name in &args.pcr_panel {
        pcr_runs.append(
            preconfigured::bounded::builtin(panel_name, budget)
                .with_context(|| format!("Failed to load panel: {panel_name}"))?,
        )?;
    }
    for source in &args.pcr_panel_file {
        pcr_runs.append(
            preconfigured::bounded::source(source, budget)
                .with_context(|| format!("Failed to load panel: {source}"))?,
        )?;
    }
    for specification in &args.pcr_primers {
        pcr_runs.push(
            parse_pcr_primers_owned(specification, budget).with_context(|| {
                format!("Could not parse primer specification: \"{specification}\"")
            })?,
        )?;
    }

    // Validate all primer parameters, collecting all errors
    let mut total_errors = 0usize;
    let mut error_report = BudgetString::new(budget, MemoryClass::RetainedResults);
    for pcr_params in pcr_runs.iter() {
        let errors = pcr::validate_pcr_params(pcr_params, budget)?;
        if !errors.is_empty() {
            total_errors += errors.len();
            error_report.push_str(&BudgetString::format(
                format_args!("\n  {} ({}):\n", pcr_params.gene_name, pcr_params.source),
                budget,
                MemoryClass::RetainedResults,
            )?)?;
            for (err, suggestion) in errors.iter() {
                error_report.push_str(&BudgetString::format(
                    format_args!("    - {}\n      Suggestion: {}\n", err, suggestion),
                    budget,
                    MemoryClass::RetainedResults,
                )?)?;
            }
        }
    }
    if total_errors > 0 {
        return Err(anyhow::Error::msg(BudgetString::format(
            format_args!(
                "Primer validation failed ({} error{}):\n{}",
                total_errors,
                if total_errors == 1 { "" } else { "s" },
                error_report
            ),
            budget,
            MemoryClass::RetainedResults,
        )?));
    }

    for pcr_params in pcr_runs.iter() {
        pcr::validate_primer_expansion(pcr_params, args.k).map_err(|error| {
            anyhow::anyhow!(
                "Primer expansion validation failed for {} ({}): {error:#}",
                pcr_params.gene_name,
                pcr_params.source
            )
        })?;
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
    let mut gene_names = BudgetMap::in_class(budget, MemoryClass::RetainedResults);
    for pcr_params in pcr_runs.iter() {
        ensure!(
            gene_names.get(pcr_params.gene_name.as_str()).is_none(),
            "Duplicate gene name '{}' (from {})",
            pcr_params.gene_name,
            pcr_params.source
        );
        gene_names.insert(pcr_params.gene_name.as_str(), ())?;
    }
    drop(gene_names);

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
pub(crate) fn resolve_sample_name(
    args: &Args,
    budget: Option<&std::sync::Arc<crate::kmer::budget::ResourceBudget>>,
) -> Result<(crate::kmer::budget::text::BudgetString, Option<EnaResult>)> {
    crate::io::resolve_input_sample(args, budget)
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
        args.fastq_record_bytes > 0,
        "FASTQ record byte limit must be positive"
    );
    ensure!(
        args.count_table_budget_mib > 0,
        "Count table budget must be positive"
    );
    ensure!(
        args.count_owner_threshold_mib > 0,
        "Count owner transfer threshold must be positive"
    );
    if args.graph_buffer_budget_mib() > 0 {
        ensure!(
            args.count_backend() != "legacy" && args.chunks == 0,
            "Graph-buffer budgeting requires an exact counter and --chunks 0"
        );
        ensure!(
            args.graph_buffer_target_mib > 0
                && args.graph_buffer_target_mib <= args.graph_buffer_budget_mib(),
            "Per-target graph-buffer limit must be positive and fit the shared graph-buffer budget"
        );
        args.graph_buffer_budget_mib()
            .checked_mul(1 << 20)
            .context("Graph-buffer budget overflow")?;
        args.graph_buffer_target_mib
            .checked_mul(1 << 20)
            .context("Per-target graph-buffer budget overflow")?;
    }
    ensure!(
        args.count_backend() == "legacy" || args.chunks == 0,
        "An exact count backend requires --chunks 0; incremental histograms use the legacy counter"
    );
    ensure!(
        args.count_sort_buffer_mib > 0,
        "External sort budget must be positive"
    );
    ensure!(
        args.count_disk_budget_mib > 0,
        "External disk budget must be positive"
    );
    if matches!(
        args.count_backend(),
        "external" | "adaptive" | "adaptive-parallel" | "pipeline" | "hybrid"
    ) {
        let directory = args
            .count_temp_dir
            .clone()
            .unwrap_or_else(std::env::temp_dir);
        ensure!(
            directory.is_dir(),
            "External count temporary directory does not exist: {}",
            directory.display()
        );
    }
    ensure!(
        !args.paired || args.max_reads != Some(u64::MAX),
        "Paired --max-reads cannot be rounded above u64::MAX"
    );
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

    // Validate --paired: requires exactly 2 input files
    if args.paired {
        if args.ena.is_some() {
            bail!("--paired cannot be combined with --ena. Provide exactly 2 local input files.");
        }
        match &args.input {
            Some(files) if files.len() == 2 => {} // ok
            Some(files) => bail!(
                "--paired requires exactly 2 input files (R1, R2), got {}",
                files.len()
            ),
            None => bail!("--paired requires exactly 2 input files (R1, R2). Cannot use stdin."),
        }
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
    eprintln!(
        "  PCR output:     {}",
        if args.no_emit_uncertain {
            "records without assembly warnings only"
        } else {
            "records with and without assembly warnings"
        }
    );
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn normal_defaults_preserve_incremental_and_explicit_development_paths() {
        let normal = Args::parse_from(["sharkmer", "-s", "sample"]);
        assert_eq!(normal.count_backend(), "hybrid");
        assert!(normal.uses_production_defaults());
        assert_eq!(normal.graph_buffer_budget_mib(), 1024);
        assert_eq!(normal.graph_buffer_target_mib, 256);
        validate_args(&normal, &[]).unwrap();

        let incremental = Args::parse_from(["sharkmer", "-s", "sample", "--chunks", "2"]);
        assert_eq!(incremental.count_backend(), "legacy");
        assert!(!incremental.uses_production_defaults());
        assert_eq!(incremental.graph_buffer_budget_mib(), 0);
        validate_args(&incremental, &[]).unwrap();

        for backend in ["legacy", "memory", "hybrid"] {
            let explicit =
                Args::parse_from(["sharkmer", "-s", "sample", "--count-backend", backend]);
            assert_eq!(explicit.count_backend(), backend);
            assert!(!explicit.uses_production_defaults());
            assert_eq!(explicit.graph_buffer_budget_mib(), 0);
            validate_args(&explicit, &[]).unwrap();
        }
        let zero = Args::parse_from(["sharkmer", "-s", "sample", "--graph-buffer-budget-mib", "0"]);
        assert_eq!(zero.count_backend(), "hybrid");
        assert_eq!(zero.graph_buffer_budget_mib(), 0);
    }

    #[test]
    fn test_parse_pcr_primers_duplicate_key_rejected() {
        let spec = "forward=AAAA,reverse=TTTT,forward=CCCC,name=X";
        let err = parse_pcr_primers_string(spec)
            .expect_err("duplicate forward key must produce an error");
        let msg = format!("{}", err);
        assert!(
            msg.contains("Duplicate parameter"),
            "error should mention duplicate parameter: {}",
            msg
        );
    }

    #[test]
    fn test_parse_pcr_primers_unique_keys_ok() {
        let spec = "forward=AAAA,reverse=TTTT,name=X,max-length=500";
        let params = parse_pcr_primers_string(spec).expect("unique keys should parse");
        assert_eq!(params.gene_name, "X");
        assert_eq!(params.forward_seq, "AAAA");
        assert_eq!(params.reverse_seq, "TTTT");
        assert_eq!(params.max_length, 500);
    }

    /// Values may legitimately contain '=' (URLs in citations, equations in
    /// notes). Only the first '=' separates key from value.
    #[test]
    fn test_parse_pcr_primers_value_with_equals() {
        let spec = "forward=AAAA,reverse=TTTT,name=X,citation=https://doi.org/10.1/abc?x=1&y=2";
        let params = parse_pcr_primers_string(spec)
            .expect("values containing '=' should parse via splitn(2)");
        assert_eq!(params.gene_name, "X");
        assert_eq!(params.citation, "https://doi.org/10.1/abc?x=1&y=2");
    }

    #[test]
    fn test_primer_budget_fails_during_collection_before_input_resolution() {
        let primer = format!("{}{}{}", "N".repeat(5), "B".repeat(2), "A".repeat(23));
        let args = Args::try_parse_from([
            "sharkmer",
            "--ena",
            "SRR000000",
            "-k",
            "31",
            "--pcr-primers",
            &format!("name=budget,forward={primer},reverse=ACGTACGT,trim=30,mismatches=1"),
        ])
        .unwrap();

        let error = collect_pcr_params(&args, None).unwrap_err();
        assert!(format!("{error:#}").contains("generated mismatch candidates"));
    }

    #[test]
    fn test_builtin_panel_passes_primer_budget_validation() {
        let args = Args::try_parse_from(["sharkmer", "--pcr-panel", "cnidaria"]).unwrap();
        let params = collect_pcr_params(&args, None).unwrap();

        assert!(!params.is_empty());
    }

    #[test]
    fn direct_counter_controls_are_hidden_and_reject_incompatible_legacy_mode() {
        let help = Args::try_parse_from(["sharkmer", "--help"])
            .unwrap_err()
            .to_string();
        assert!(!help.contains("--count-backend"));
        let args = Args::try_parse_from(["sharkmer", "--count-backend", "memory", "--chunks", "2"])
            .unwrap();
        assert!(
            validate_args(&args, &[])
                .unwrap_err()
                .to_string()
                .contains("requires --chunks 0")
        );
        let args = Args::try_parse_from([
            "sharkmer",
            "--paired",
            "--max-reads",
            "18446744073709551615",
        ])
        .unwrap();
        assert!(
            validate_args(&args, &[])
                .unwrap_err()
                .to_string()
                .contains("rounded above u64::MAX")
        );
    }

    #[test]
    fn test_counting_only_kmer_length_one_remains_valid() {
        let args = Args::try_parse_from(["sharkmer", "-k", "1"]).unwrap();
        let params = collect_pcr_params(&args, None).unwrap();

        validate_args(&args, &params).unwrap();
    }

    #[test]
    fn test_withheld_path_diagnostics_flag_is_hidden_but_parseable() {
        let args =
            Args::try_parse_from(["sharkmer", "--diagnose-withheld-paths", "reads.fastq"]).unwrap();

        assert!(args.diagnose_withheld_paths);
        let help = Args::try_parse_from(["sharkmer", "--help"])
            .unwrap_err()
            .to_string();
        assert!(!help.contains("--diagnose-withheld-paths"));
    }

    #[test]
    fn uncertain_delivery_defaults_on_with_explicit_opt_out_and_legacy_opt_in() {
        let default_args = Args::try_parse_from(["sharkmer", "reads.fastq"]).unwrap();
        assert!(!default_args.no_emit_uncertain);
        let enabled =
            Args::try_parse_from(["sharkmer", "--emit-uncertain", "reads.fastq"]).unwrap();
        assert!(enabled.emit_uncertain);
        assert!(!enabled.no_emit_uncertain);
        assert!(!enabled.diagnose_withheld_paths);
        let disabled =
            Args::try_parse_from(["sharkmer", "--no-emit-uncertain", "reads.fastq"]).unwrap();
        assert!(disabled.no_emit_uncertain);
        assert!(
            Args::try_parse_from(["sharkmer", "--emit-uncertain", "--no-emit-uncertain"]).is_err()
        );
        let help = Args::try_parse_from(["sharkmer", "--help"])
            .unwrap_err()
            .to_string();
        assert!(help.contains("--no-emit-uncertain"));
    }

    #[test]
    fn combined_output_allows_legacy_gene_names_ending_in_uncertain() {
        let specifications = [
            "name=A,forward=ACGTACGTACGTACGTACGT,reverse=TGCATGCATGCATGCATGCA",
            "name=A.uncertain,forward=ACGTACGTACGTACGTACGT,reverse=TGCATGCATGCATGCATGCA",
        ];
        let args = Args::try_parse_from([
            "sharkmer",
            "--ena",
            "SRR1",
            "--pcr-primers",
            specifications[0],
            "--pcr-primers",
            specifications[1],
        ])
        .unwrap();
        assert_eq!(collect_pcr_params(&args, None).unwrap().len(), 2);
    }
}

#[cfg(test)]
mod ownership_tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};
    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }
    #[test]
    fn collected_inline_configuration_keeps_exact_live_owners_and_typed_refusal() {
        let specification = "name=one,forward=TAACATATGCTTGTCTCAAAG,reverse=CCTGTATTGTTATTTTTCGTCAC,trim=15,mismatches=0";
        let args =
            Args::try_parse_from(["sharkmer", "-s", "sample", "--pcr-primers", specification])
                .unwrap();
        let ledger = budget(1 << 20);
        let params = collect_pcr_params(&args, Some(&ledger)).unwrap();
        assert_eq!(params.len(), 1);
        let actual = params[0].configuration_bytes().unwrap()
            + (params.values.capacity() * size_of::<pcr::PCRParams>()) as u64
            + (params.owners.capacity() * size_of::<Owner>()) as u64;
        assert_eq!(ledger.snapshot().reserved_memory_bytes, actual);
        assert_eq!(
            params[0].source,
            format!("--pcr-primers \"{specification}\"")
        );
        drop(params);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        let small = budget(16);
        let error = collect_pcr_params(&args, Some(&small)).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(small.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn invalid_primer_diagnostics_keep_order_and_reservations_until_error_drop() {
        let args = Args::try_parse_from([
            "sharkmer",
            "-s",
            "sample",
            "--pcr-primers",
            "name=bad,forward=XZ,reverse=XZ,min-count=1",
        ])
        .unwrap();
        let ledger = budget(1 << 20);
        let error = collect_pcr_params(&args, Some(&ledger)).unwrap_err();
        let text = error.to_string();
        assert!(
            text.contains("Primer validation failed (4 errors)"),
            "{text}"
        );
        let forward = text
            .find("Invalid nucleotide(s) X, Z in forward primer XZ")
            .unwrap();
        let reverse = text
            .find("Invalid nucleotide(s) X, Z in reverse primer XZ")
            .unwrap();
        assert!(forward < reverse);
        assert!(text.contains("min-count is 1, must be at least 2"));
        assert!(text.contains("Forward and reverse primers are identical: XZ"));
        assert!(ledger.snapshot().reserved_memory_bytes > 0);
        drop(error);
        assert_eq!(ledger.snapshot().reserved_memory_bytes, 0);
        assert!(ledger.snapshot().valid);
    }
}

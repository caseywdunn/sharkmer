use anyhow::{Context, Result};
use clap::Parser;
use log::{debug, info};
use peak_alloc::PeakAlloc;
use std::io::IsTerminal;

mod cli;
mod format;
mod io;
mod kmer;
mod pcr;
mod stats;

use cli::Args;
use stats::RunStats;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> Result<()> {
    let start_run = std::time::Instant::now();

    // Parse CLI arguments
    let args = Args::parse();

    // Initialize logging
    cli::init_logging(args.verbose, &args.color);

    // Handle early-exit flags (--completions, --cite, --list-panels, --export-panel, --help-pcr)
    cli::handle_early_exits(&args)?;

    // Collect and validate PCR primer specifications from all sources
    let pcr_runs = cli::collect_pcr_params(&args)?;

    // Handle --validate-panels (prints and exits)
    if args.validate_panels {
        cli::handle_validate_panels(&pcr_runs)?;
    }

    // Derive sample name from --sample or --sra ENA metadata
    let (sample, cached_ena_result) = cli::resolve_sample_name(&args)?;

    info!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
    debug!("{:?}", args);

    // Resolve the output directory path (creation deferred until after dry-run check)
    let outdir_str = args
        .outdir
        .to_str()
        .context("Output directory path contains invalid UTF-8")?;
    let directory = if outdir_str.ends_with('/') {
        outdir_str.to_string()
    } else {
        format!("{}/", outdir_str)
    };

    let k = args.k;

    // Validate arguments
    cli::validate_args(&args, &pcr_runs)?;

    // Handle --dry-run (prints and exits)
    if args.dry_run {
        cli::handle_dry_run(&args, &sample, &directory, &pcr_runs);
    }

    // Create the output directory now that we know this is not a dry run
    std::fs::create_dir_all(&directory)
        .with_context(|| format!("Failed to create output directory: {}", directory))?;

    // Set the number of threads for Rayon to use
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    let show_progress = std::io::stderr().is_terminal();

    // Ingest FASTQ reads from all input sources
    let (state, n_reads_ingested, n_bases_ingested, n_kmers_ingested) =
        io::ingest_reads(&args, k, cached_ena_result, show_progress)?;

    // Consolidate chunks and optionally write histograms
    let mut state = state;
    let (kmer_counts, n_singleton_kmers) = io::consolidate_and_histogram(
        &mut state,
        &args,
        k,
        &sample,
        &directory,
        n_kmers_ingested,
        show_progress,
    )?;

    // Run in silico PCR
    let pcr_results = stats::run_pcr(
        &kmer_counts,
        &pcr_runs,
        &sample,
        &directory,
        args.min_kmer_count,
        show_progress,
    )?;

    // Build and write run statistics
    let command = std::env::args().collect::<Vec<String>>().join(" ");
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

    stats::write_stats(&run_stats, &directory, &sample)?;

    // Print final summary
    stats::print_summary(&run_stats, start_run.elapsed());

    Ok(())
}

#[cfg(test)]
mod tests {}

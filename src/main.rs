#![allow(clippy::manual_is_multiple_of)]

use anyhow::{Context, Result};
use clap::Parser;
use log::{debug, info};
use peak_alloc::PeakAlloc;
use std::io::IsTerminal;

mod cache;
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
    cli::init_logging(args.verbose, args.quiet, &args.color);

    // Handle early-exit flags (--completions, --cite, --list-panels, --export-panel, --help-pcr)
    cli::handle_early_exits(&args)?;

    // Collect and validate PCR primer specifications from all sources
    let mut pcr_runs = cli::collect_pcr_params(&args)?;

    // Apply CLI tuning overrides to all PCR params
    for p in &mut pcr_runs {
        p.max_dfs_states = args.max_dfs_states;
        p.max_paths_per_pair = args.max_paths_per_pair;
        p.max_node_visits = args.max_node_visits;
        p.max_primer_kmers = args.max_primer_kmers;
        p.max_seed_nodes = args.max_seed_nodes;
        p.high_coverage_ratio = args.high_coverage_ratio;
        p.tip_coverage_fraction = args.tip_coverage_fraction;
    }

    // Handle --validate-panels (prints and exits)
    if args.validate_panels {
        cli::handle_validate_panels(&pcr_runs)?;
    }

    // Derive sample name from --sample or --ena ENA metadata
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

    // Pre-encode primer Oligos for Pass 1 read retention (opt-in via --read-eval)
    let oligo_filter = if args.read_eval && !pcr_runs.is_empty() {
        let sets = pcr::preprocess_primer_oligos(&pcr_runs, k)?;
        info!(
            "Read-backed seed evaluation enabled: pre-encoded primer Oligos for {} gene(s)",
            sets.len()
        );
        Some(io::OligoFilter::new(&sets, k))
    } else {
        None
    };

    // Set up read cache for remote downloads
    let cache_config = if !args.no_cache && args.ena.is_some() {
        let cc = cache::CacheConfig::new(args.cache_dir.as_deref())?;
        info!("Read cache: {}", cc.cache_dir.display());
        Some(cc)
    } else {
        None
    };

    // Ingest FASTQ reads from all input sources (Pass 1: kmer counting + read retention)
    let (state, n_reads_ingested, n_bases_ingested, n_kmers_ingested, read_plan) =
        io::ingest_reads(
            &args,
            k,
            cached_ena_result,
            cache_config.as_ref(),
            show_progress,
            oligo_filter.as_ref(),
        )?;

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

    // Pass 2: re-read sequences for read threading (opt-in via --read-threading)
    let threading_reads = if args.read_threading && !pcr_runs.is_empty() {
        match &read_plan.source {
            io::ReadSourcePlan::Unavailable => {
                info!("Read threading unavailable (stdin input); using kmer-only scoring");
                None
            }
            _ => Some(io::reread_sequences(&read_plan, show_progress)?),
        }
    } else {
        None
    };

    // Run in silico PCR
    let pcr_results = stats::run_pcr(
        &kmer_counts,
        &pcr_runs,
        &sample,
        &directory,
        args.min_kmer_count,
        args.dump_graph,
        show_progress,
        threading_reads.as_deref(),
        &state.retained_reads,
        args.max_nodes,
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

#![allow(clippy::manual_is_multiple_of)]

use anyhow::{Context, Result};
use clap::Parser;
use log::{debug, info};
use peak_alloc::PeakAlloc;
use std::io::IsTerminal;

mod cache;
mod cli;
mod codec;
mod configuration;
mod fastq;
mod format;
mod http;
mod input_identity;
mod io;
mod kmer;
mod output;
mod pcr;
mod process_memory;
mod replay;
mod stats;

use cli::Args;
use stats::{InputSourceStats, RunStats, StageTimings};

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

fn main() -> Result<()> {
    let start_run = std::time::Instant::now();

    let process_address_space = process_memory::from_environment()?;

    // Parse CLI arguments
    let args = Args::parse();
    let process_address_space =
        process_memory::for_normal_run(process_address_space, args.uses_production_defaults())?;

    // Initialize logging
    cli::init_logging(args.verbose, args.quiet, &args.color);

    // Handle early-exit flags (--completions, --cite, --list-panels, --export-panel, --help-pcr)
    cli::handle_early_exits(&args)?;

    // Numeric preparation creates no tables, workers or private paths. Reuse
    // this same ledger for metadata/sample ownership and later exact ingestion.
    let ingestion = if args.dry_run || args.validate_panels {
        None
    } else {
        Some(io::prepare_ingestion(&args)?)
    };
    let preview_budget = if ingestion
        .as_ref()
        .is_none_or(|plan| plan.input_budget().is_none())
    {
        Some(io::metadata_preview_budget(
            process_memory::metadata_ceiling(process_address_space.as_ref()),
        )?)
    } else {
        None
    };

    // Reject an impossible normal process plan before any panel body is read.
    if args.uses_production_defaults()
        && let Some(plan) = &ingestion
    {
        plan.check_memory_envelope(
            process_address_space
                .as_ref()
                .map_or(process_memory::NORMAL_ADDRESS_SPACE_BYTES, |cap| {
                    cap.effective_limit_bytes()
                }),
            process_memory::NORMAL_PROCESS_HEADROOM_BYTES,
        )?;
    }

    // Configuration is admitted on the original run ledger before decoding.
    let mut pcr_runs = cli::collect_pcr_params(
        &args,
        ingestion
            .as_ref()
            .and_then(|plan| plan.input_budget())
            .or(preview_budget.as_ref()),
    )?;

    // Apply CLI tuning overrides to all PCR params.
    //
    // Global CLI values always win: unset flags supply their hardcoded defaults
    // and overwrite these panel compatibility fields. This is intentional for
    // hidden engine diagnostics, not per-gene tuning. dedup_edit_threshold
    // remains the supported per-primer runtime tuning field.
    for p in pcr_runs.iter_mut() {
        p.max_dfs_states = args.max_dfs_states;
        p.max_paths_per_pair = args.max_paths_per_pair;
        p.max_node_visits = args.max_node_visits;
        p.max_primer_kmers = args.max_primer_kmers;
        p.high_coverage_ratio = args.high_coverage_ratio;
        p.tip_coverage_fraction = args.tip_coverage_fraction;
        p.diagnose_withheld_paths = args.diagnose_withheld_paths;
        p.emit_uncertain = args.emit_uncertain || !args.no_emit_uncertain;
    }

    // Handle --validate-panels early: pure primer check, no sample/output needed
    if args.validate_panels {
        cli::handle_validate_panels(&pcr_runs)?;
    }

    // Derive sample name from --sample or --ena ENA metadata
    let (sample, cached_ena_result) = cli::resolve_sample_name(
        &args,
        ingestion
            .as_ref()
            .and_then(|plan| plan.input_budget())
            .or(preview_budget.as_ref()),
    )?;

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
    let ingestion = ingestion.context("Prepared ingestion is unavailable")?;

    // Create the output directory now that we know this is not a dry run
    std::fs::create_dir_all(&directory)
        .with_context(|| format!("Failed to create output directory: {}", directory))?;
    let mut output_transaction = output::OutputTransaction::begin_with_budget(
        &args.outdir,
        &sample,
        ingestion.result_budget(),
    )?;

    let ingestion = if ingestion.needs_workspace() {
        ingestion.with_workspace(output_transaction.workspace_root(args.count_temp_dir.as_deref())?)
    } else {
        ingestion
    };

    let show_progress = std::io::stderr().is_terminal();

    // Set up read cache for remote downloads
    let cache_config = if !args.no_cache && args.ena.is_some() {
        let cc = cache::CacheConfig::new(args.cache_dir.as_deref())?;
        info!("Read cache: {}", cc.cache_dir.display());
        Some(cc)
    } else {
        None
    };

    // Ingest FASTQ reads from all input sources (Pass 1: kmer counting)
    let ingest_start = std::time::Instant::now();
    let (state, n_reads_ingested, n_bases_ingested, n_kmers_ingested, mut read_plan) =
        io::ingest_reads(
            ingestion,
            cached_ena_result,
            cache_config.as_ref(),
            show_progress,
        )?;
    let read_ingest_s = ingest_start.elapsed().as_secs_f64();
    // Consolidate chunks and optionally write histograms
    let mut state = state;
    let finalize_start = std::time::Instant::now();
    let (kmer_counts, n_singleton_kmers) = io::consolidate_and_histogram(
        &mut state,
        &args,
        k,
        &sample,
        &directory,
        n_kmers_ingested,
        show_progress,
    )?;
    let count_finalize_s = finalize_start.elapsed().as_secs_f64();
    let count_table_capacity = kmer_counts.capacity();

    let result_budget = if args.graph_buffer_budget_mib() > 0 {
        kmer_counts.resource_budget()
    } else {
        None
    };
    let mut input_source_inputs = crate::kmer::budget::collections::BudgetVec::in_class(
        result_budget,
        crate::kmer::budget::MemoryClass::RetainedResults,
    );
    let input_source_kind = match &read_plan.source {
        io::ReadSourcePlan::LocalFiles(paths) | io::ReadSourcePlan::CachedRemote(paths) => {
            for path in paths {
                input_source_inputs.push(crate::kmer::budget::text::BudgetString::format(
                    format_args!("{}", path.display()),
                    result_budget,
                    crate::kmer::budget::MemoryClass::RetainedResults,
                )?)?;
            }
            if matches!(&read_plan.source, io::ReadSourcePlan::LocalFiles(_)) {
                "local_files"
            } else {
                "cached_remote"
            }
        }
        io::ReadSourcePlan::UncachedRemote(urls) => {
            for url in urls {
                input_source_inputs.push(crate::kmer::budget::text::BudgetString::from_str(
                    url,
                    result_budget,
                    crate::kmer::budget::MemoryClass::RetainedResults,
                )?)?;
            }
            "uncached_remote"
        }
        io::ReadSourcePlan::Unavailable => "unavailable",
    };
    let input_source = InputSourceStats {
        kind: input_source_kind,
        inputs: input_source_inputs,
        paired: read_plan.paired,
        max_reads: read_plan.max_reads,
    };

    // Pass 2: re-read sequences for read threading (opt-in via --read-threading)
    let threading_input_start = std::time::Instant::now();
    let threading_reads = if args.read_threading && !pcr_runs.is_empty() {
        if let Some(replay) = read_plan.replay.take() {
            Some(io::ThreadingInput::Replay(Box::new(replay)))
        } else {
            match &read_plan.source {
                io::ReadSourcePlan::Unavailable => {
                    info!("Read threading unavailable (stdin input); using kmer-only scoring");
                    None
                }
                _ => Some(io::ThreadingInput::Memory(io::reread_sequences(
                    &read_plan,
                    show_progress,
                )?)),
            }
        }
    } else {
        None
    };
    let read_threading_input_s = if args.read_threading && !pcr_runs.is_empty() {
        Some(threading_input_start.elapsed().as_secs_f64())
    } else {
        None
    };

    // Resolve the per-gene/threshold node budget from user input or data volume
    let node_budget_global = match args.node_budget_global {
        Some(budget) => {
            info!(
                "Per-gene/threshold node budget: {} (pinned via --node-budget-global)",
                budget
            );
            budget
        }
        None => {
            let budget = pcr::compute_node_budget(n_bases_ingested);
            info!(
                "Per-gene/threshold node budget: {} (auto, based on {} bp ingested; ~{}M reads at 150bp). Pin with --node-budget-global.",
                budget,
                n_bases_ingested,
                n_bases_ingested / 150_000_000
            );
            budget
        }
    };

    // Counting owners have joined before assembly workers start. Keeping an
    // idle assembly pool during ingestion wastes stacks and allocator arenas.
    // The coordinator participates in assembly rather than waiting outside an
    // additional full-sized pool. This keeps OS threads within --threads and
    // avoids an extra stack/allocator arena at the counting/assembly boundary.
    // Counting-only runs never need an assembly pool.
    if !pcr_runs.is_empty() {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .use_current_thread()
            .thread_name(|index| format!("assembly-{index}"))
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    let graph_buffers = if args.graph_buffer_budget_mib() > 0 && !pcr_runs.is_empty() {
        let budget = kmer_counts
            .resource_budget()
            .context("Graph-buffer budgeting requires exact counts")?;
        Some(pcr::memory::admission::GraphBufferPool::new(
            std::sync::Arc::clone(budget),
            args.graph_buffer_target_mib
                .checked_mul(1 << 20)
                .context("Per-target graph-buffer budget overflow")?,
        )?)
    } else {
        None
    };

    // Run in silico PCR
    let pcr_start = std::time::Instant::now();
    let pcr_results = stats::run_pcr(
        &kmer_counts.filtered_view(0),
        &pcr_runs,
        &sample,
        &directory,
        args.min_kmer_count,
        args.dump_graph,
        show_progress,
        threading_reads.as_ref().map(io::ThreadingInput::source),
        node_budget_global,
        &mut output_transaction,
        graph_buffers.as_ref(),
    )?;
    let pcr_s = pcr_start.elapsed().as_secs_f64();

    // Build and write run statistics
    let mut command = crate::kmer::budget::text::BudgetString::new(
        result_budget,
        crate::kmer::budget::MemoryClass::RetainedResults,
    );
    for (index, argument) in std::env::args().enumerate() {
        if index > 0 {
            command.push_str(" ")?;
        }
        command.push_str(&argument)?;
    }
    let run_stats = RunStats {
        process_address_space,
        assembly_graph_buffers: graph_buffers
            .as_ref()
            .map(pcr::memory::admission::GraphBufferPool::summary)
            .transpose()?,
        read_replay: threading_reads
            .as_ref()
            .map(|input| input.summary(result_budget))
            .transpose()?
            .flatten(),
        sharkmer_version: env!("CARGO_PKG_VERSION"),
        command,
        sample: crate::kmer::budget::text::BudgetString::from_str(
            &sample,
            result_budget,
            crate::kmer::budget::MemoryClass::RetainedResults,
        )?,
        run_id: crate::kmer::budget::text::BudgetString::from_str(
            output_transaction.run_id(),
            result_budget,
            crate::kmer::budget::MemoryClass::RetainedResults,
        )?,
        output_manifest: crate::kmer::budget::text::BudgetString::from_str(
            output_transaction.manifest_basename(),
            result_budget,
            crate::kmer::budget::MemoryClass::RetainedResults,
        )?,
        run_status: "complete",
        product_output_schema_version: 4,
        product_output_layout: "combined_warnings",
        input_source,
        kmer_length: args.k,
        chunks: args.chunks,
        n_reads_read: state.n_reads_read,
        n_bases_read: state.n_bases_read,
        n_subreads_ingested: n_reads_ingested,
        n_bases_ingested,
        n_kmers: n_kmers_ingested,
        n_multi_kmers: n_singleton_kmers.map(|s| n_kmers_ingested.saturating_sub(s)),
        n_singleton_kmers,
        count_table_capacity,
        count_store: kmer_counts.summary(result_budget)?,
        count_fingerprint: if args.count_fingerprint {
            Some(kmer_counts.fingerprint(result_budget)?)
        } else {
            None
        },
        peak_memory_bytes: PEAK_ALLOC.peak_usage() as u64,
        stage_timings: StageTimings {
            read_ingest_s,
            count_finalize_s,
            read_threading_input_s,
            pcr_s,
            pipeline_total_s: start_run.elapsed().as_secs_f64(),
        },
        pcr_results,
    };

    stats::write_stats(&run_stats, &mut output_transaction)?;
    if let Some(reads) = threading_reads {
        reads.close()?;
    }
    if let Some(unused_replay) = read_plan.replay.take() {
        unused_replay.close()?;
    }
    kmer_counts.close()?;
    output_transaction.commit()?;

    // Print final summary
    stats::print_summary(&run_stats, start_run.elapsed());

    Ok(())
}

#[cfg(test)]
mod tests {}

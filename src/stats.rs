use crate::format::{format_bytes, format_count, format_duration};
use crate::io::{warn_if_exists, write_fasta_record};
use crate::kmer::KmerCounts;
use crate::pcr;
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use rayon::prelude::*;
use serde::Serialize;

/// Stats for a single PCR gene result, included in the YAML stats output.
#[derive(Serialize)]
pub(crate) struct PcrGeneResult {
    pub(crate) gene_name: String,
    pub(crate) status: String,
    pub(crate) n_products: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) product_lengths: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_file: Option<String>,
}

/// Structured run statistics, serialized as YAML.
#[derive(Serialize)]
pub(crate) struct RunStats {
    pub(crate) sharkmer_version: String,
    pub(crate) command: String,
    pub(crate) sample: String,
    pub(crate) kmer_length: usize,
    pub(crate) chunks: usize,
    pub(crate) n_reads_read: u64,
    pub(crate) n_bases_read: u64,
    pub(crate) n_subreads_ingested: u64,
    pub(crate) n_bases_ingested: u64,
    pub(crate) n_kmers: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) n_multi_kmers: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) n_singleton_kmers: Option<u64>,
    pub(crate) peak_memory_bytes: u64,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) pcr_results: Vec<PcrGeneResult>,
}

/// Run in silico PCR for all primer pairs, write FASTA output files, and print results.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_pcr(
    kmer_counts: &KmerCounts,
    pcr_runs: &[pcr::PCRParams],
    sample: &str,
    directory: &str,
    min_kmer_count: u32,
    dump_graph: bool,
    show_progress: bool,
    reads: Option<&[crate::io::ReadRecord]>,
    retained_reads: &crate::io::RetainedReads,
    max_nodes: usize,
) -> Result<Vec<PcrGeneResult>> {
    let mut pcr_results: Vec<PcrGeneResult> = Vec::new();

    if pcr_runs.is_empty() {
        return Ok(pcr_results);
    }

    info!("Running in silico PCR...");
    let spinner_style =
        ProgressStyle::with_template("{spinner:.cyan} {msg}").expect("valid spinner template");
    let pcr_spinner = if show_progress {
        let sp = ProgressBar::new_spinner();
        sp.set_style(spinner_style);
        sp.set_message("Running in silico PCR...");
        sp.enable_steady_tick(std::time::Duration::from_millis(80));
        sp
    } else {
        ProgressBar::hidden()
    };

    // Create a filtered view of kmer counts for PCR (no data copied)
    info!("Filtering kmers with count < {} before PCR", min_kmer_count);
    let kmer_counts_pcr = kmer_counts.filtered_view(min_kmer_count);

    // Collect all retained read sequences (shared across genes — each gene's
    // divergence check filters to reads matching its own seed sub_kmers)
    let all_retained: Vec<&str> = retained_reads
        .reads
        .iter()
        .map(|r| r.sequence.as_str())
        .collect();

    // Run PCR for each gene in parallel; kmer_counts_pcr and reads are read-only and shared
    let pcr_fasta_results: Vec<_> = pcr_runs
        .par_iter()
        .map(|pcr_params| {
            let fasta = pcr::do_pcr(
                &kmer_counts_pcr,
                sample,
                pcr_params,
                dump_graph,
                directory,
                reads,
                &all_retained,
                max_nodes,
            );
            (pcr_params, fasta)
        })
        .collect();
    pcr_spinner.finish_and_clear();

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
                write_fasta_record(&mut file, record).context("Failed to write FASTA record")?;
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

    // Print gene result table
    let (sym_pass, sym_fail) = if show_progress {
        ("\u{2714}", "\u{2718}") // checkmark, x-mark
    } else {
        ("+", "-")
    };
    for result in &pcr_results {
        if result.status == "success" {
            let lengths: Vec<String> = result
                .product_lengths
                .iter()
                .map(|l| l.to_string())
                .collect();
            warn!(
                "  {} {} ({} product{}, {} bp)",
                sym_pass,
                result.gene_name,
                result.n_products,
                if result.n_products == 1 { "" } else { "s" },
                lengths.join(", ")
            );
        } else {
            warn!("  {} {} (no products)", sym_fail, result.gene_name);
        }
    }

    info!("Done running in silico PCR");

    Ok(pcr_results)
}

/// Write run statistics as a YAML file.
pub(crate) fn write_stats(run_stats: &RunStats, directory: &str, sample: &str) -> Result<()> {
    info!("Writing stats to file...");
    let stats_path = format!("{}{}.stats.yaml", directory, sample);
    warn_if_exists(&stats_path);
    let file_stats = std::fs::File::create(&stats_path).context("Failed to create stats file")?;
    serde_yml::to_writer(file_stats, run_stats).context("Failed to write stats YAML")?;
    Ok(())
}

/// Print the final summary line to stderr.
pub(crate) fn print_summary(run_stats: &RunStats, elapsed: std::time::Duration) {
    let elapsed_str = format_duration(elapsed);
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
        if run_stats.chunks > 0 {
            warn!(
                "sharkmer complete: {} reads, {} kmers, {} chunks, peak mem {}, {}",
                reads_str,
                kmers_str,
                run_stats.chunks,
                format_bytes(run_stats.peak_memory_bytes),
                elapsed_str
            );
        } else {
            warn!(
                "sharkmer complete: {} reads, {} kmers, peak mem {}, {}",
                reads_str,
                kmers_str,
                format_bytes(run_stats.peak_memory_bytes),
                elapsed_str
            );
        }
    }
}

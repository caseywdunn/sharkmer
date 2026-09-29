use crate::format::{format_bytes, format_count, format_duration};
use crate::io::write_fasta_record;
use crate::kmer::KmerCounts;
use crate::output::OutputTransaction;
use crate::pcr;
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use rayon::prelude::*;
use serde::Serialize;

/// Provenance for the uncertain subset of a target's combined FASTA.
#[derive(Serialize)]
pub(crate) struct UncertainProductOutput {
    pub(crate) n_products: usize,
    pub(crate) product_lengths: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_file: Option<String>,
    #[serde(flatten)]
    pub(crate) assessment: pcr::UncertainProducts,
}

/// Stats for a single PCR gene result, included in the YAML stats output.
#[derive(Serialize)]
pub(crate) struct PcrGeneResult {
    pub(crate) gene_name: String,
    pub(crate) status: String,
    pub(crate) n_products: usize,
    pub(crate) n_standard_products: usize,
    pub(crate) n_uncertain_products: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) product_lengths: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) failure_reason: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub(crate) threshold_diagnostics: Vec<pcr::PcrThresholdDiagnostic>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) uncertain_products: Option<UncertainProductOutput>,
}

#[derive(Serialize)]
pub(crate) struct StageTimings {
    pub(crate) read_ingest_s: f64,
    pub(crate) count_finalize_s: f64,
    pub(crate) read_threading_input_s: Option<f64>,
    pub(crate) pcr_s: f64,
    pub(crate) pipeline_total_s: f64,
}

#[derive(Serialize)]
pub(crate) struct InputSourceStats {
    pub(crate) kind: String,
    pub(crate) inputs: Vec<String>,
    pub(crate) paired: bool,
    pub(crate) max_reads: u64,
}

/// Structured run statistics, serialized as YAML.
#[derive(Serialize)]
pub(crate) struct RunStats {
    pub(crate) sharkmer_version: String,
    pub(crate) command: String,
    pub(crate) sample: String,
    pub(crate) run_id: String,
    pub(crate) output_manifest: String,
    pub(crate) run_status: String,
    pub(crate) product_output_schema_version: u32,
    pub(crate) product_output_layout: String,
    pub(crate) input_source: InputSourceStats,
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
    pub(crate) count_table_capacity: usize,
    pub(crate) peak_memory_bytes: u64,
    pub(crate) stage_timings: StageTimings,
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
    max_nodes: usize,
    output_transaction: &mut OutputTransaction,
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

    // Run PCR for each gene in parallel; kmer_counts_pcr and reads are read-only and shared
    let pcr_fasta_results: Vec<_> = pcr_runs
        .par_iter()
        .enumerate()
        .map(|(gene_index, pcr_params)| {
            let diagnostic_limits =
                pcr::withheld_diagnostic_limits_for_gene(gene_index, pcr_runs.len());
            let fasta = pcr::do_pcr(
                &kmer_counts_pcr,
                sample,
                pcr_params,
                dump_graph,
                directory,
                reads,
                max_nodes,
                diagnostic_limits,
            );
            (pcr_params, fasta)
        })
        .collect();
    pcr_spinner.finish_and_clear();

    // Write output files sequentially to maintain deterministic order
    for (pcr_params, fasta_result) in pcr_fasta_results {
        let outcome = fasta_result?;
        let pcr::PcrOutcome {
            records,
            failure_reason,
            threshold_diagnostics,
            uncertain_products,
        } = outcome;

        let n_standard_products = records.len();
        let fasta_basename = format!("{}_{}.fasta", sample, pcr_params.gene_name);
        let uncertain_products = uncertain_products
            .map(|assessment| -> Result<UncertainProductOutput> {
                let n_products = assessment.products.len();
                anyhow::ensure!(
                    n_products == assessment.unique_products,
                    "Uncertain product count differs from its retained records"
                );
                let product_lengths = assessment
                    .products
                    .iter()
                    .map(|product| product.oriented_sequence.len())
                    .collect();
                let output_file = (n_products > 0).then(|| fasta_basename.clone());
                Ok(UncertainProductOutput {
                    n_products,
                    product_lengths,
                    output_file,
                    assessment,
                })
            })
            .transpose()?;
        let n_uncertain_products = uncertain_products
            .as_ref()
            .map_or(0, |uncertain| uncertain.n_products);
        let n_products = n_standard_products + n_uncertain_products;
        let mut product_lengths: Vec<usize> =
            records.iter().map(|record| record.seq().len()).collect();
        if let Some(uncertain) = &uncertain_products {
            product_lengths.extend(&uncertain.product_lengths);
        }
        let output_file = if n_products > 0 {
            output_transaction.stage_file(&fasta_basename, |file| {
                for record in &records {
                    let description = match record.desc() {
                        Some(description) => format!("{description} assembly_status=standard"),
                        None => "assembly_status=standard".to_string(),
                    };
                    let classified = bio::io::fasta::Record::with_attrs(
                        record.id(),
                        Some(&description),
                        record.seq(),
                    );
                    write_fasta_record(file, &classified)
                        .context("Failed to write standard FASTA record")?;
                }
                if let Some(uncertain) = &uncertain_products {
                    for product in &uncertain.assessment.products {
                        let record = bio::io::fasta::Record::with_attrs(
                            &product.id,
                            Some(&product.description),
                            product.oriented_sequence.as_bytes(),
                        );
                        write_fasta_record(file, &record)
                            .context("Failed to write uncertain FASTA record")?;
                    }
                }
                Ok(())
            })?;
            Some(fasta_basename)
        } else {
            None
        };
        let status = if n_standard_products > 0 {
            "success"
        } else if n_uncertain_products > 0 {
            "uncertain"
        } else {
            "fail"
        };
        // An uncertain-only output retains the reason standard admission failed.
        let failure_reason = if n_standard_products > 0 {
            None
        } else {
            failure_reason
                .or_else(|| Some("unknown (no reason reported by PCR pipeline)".to_string()))
        };
        pcr_results.push(PcrGeneResult {
            gene_name: pcr_params.gene_name.clone(),
            status: status.to_string(),
            n_products,
            n_standard_products,
            n_uncertain_products,
            product_lengths,
            output_file,
            failure_reason,
            threshold_diagnostics,
            uncertain_products,
        });
    }

    // Print gene result table
    let (sym_pass, sym_fail) = if show_progress {
        ("\u{2714}", "\u{2718}") // checkmark, x-mark
    } else {
        ("+", "-")
    };
    for result in &pcr_results {
        if result.n_products > 0 {
            let lengths: Vec<String> = result
                .product_lengths
                .iter()
                .map(|l| l.to_string())
                .collect();
            warn!(
                "  {} {} ({} product{}: {} standard, {} uncertain; {} bp)",
                if result.n_standard_products > 0 {
                    sym_pass
                } else {
                    "?"
                },
                result.gene_name,
                result.n_products,
                if result.n_products == 1 { "" } else { "s" },
                result.n_standard_products,
                result.n_uncertain_products,
                lengths.join(", ")
            );
            if result.n_uncertain_products > 0 {
                warn!(
                    "    Uncertain records have unresolved repeat structure and unevaluated read support."
                );
            }
            if let Some(reason) = &result.failure_reason {
                warn!("    Standard admission: {}", reason);
            }
        } else {
            if let Some(ref reason) = result.failure_reason {
                warn!(
                    "  {} {} (no products, {})",
                    sym_fail, result.gene_name, reason
                );
            } else {
                warn!("  {} {} (no products)", sym_fail, result.gene_name);
            }
        }
    }

    info!("Done running in silico PCR");

    Ok(pcr_results)
}

/// Write run statistics as a YAML file.
pub(crate) fn write_stats(
    run_stats: &RunStats,
    output_transaction: &mut OutputTransaction,
) -> Result<()> {
    info!("Writing stats to file...");
    let stats_basename = format!("{}.stats.yaml", run_stats.sample);
    output_transaction.stage_file(&stats_basename, |file| {
        serde_yaml_ng::to_writer(file, run_stats).context("Failed to write stats YAML")
    })
}

/// Print the final summary line to stderr.
pub(crate) fn print_summary(run_stats: &RunStats, elapsed: std::time::Duration) {
    let elapsed_str = format_duration(elapsed);
    let reads_str = format_count(run_stats.n_reads_read);

    if !run_stats.pcr_results.is_empty() {
        let n_with_output = run_stats
            .pcr_results
            .iter()
            .filter(|r| r.n_products > 0)
            .count();
        let n_total = run_stats.pcr_results.len();
        let gene_names: Vec<&str> = run_stats
            .pcr_results
            .iter()
            .filter(|r| r.n_products > 0)
            .map(|r| r.gene_name.as_str())
            .collect();
        let genes_detail = if gene_names.is_empty() {
            String::new()
        } else if gene_names.len() <= 10 {
            format!(" ({})", gene_names.join(", "))
        } else {
            format!(
                " ({}, +{} more)",
                gene_names[..10].join(", "),
                gene_names.len() - 10
            )
        };
        let n_standard: usize = run_stats
            .pcr_results
            .iter()
            .map(|result| result.n_standard_products)
            .sum();
        let n_uncertain: usize = run_stats
            .pcr_results
            .iter()
            .map(|result| result.n_uncertain_products)
            .sum();
        warn!(
            "sharkmer complete: {} reads, {}/{} targets with output{}, {} records ({} standard, {} uncertain), allocator peak {}, {}",
            reads_str,
            n_with_output,
            n_total,
            genes_detail,
            n_standard + n_uncertain,
            n_standard,
            n_uncertain,
            format_bytes(run_stats.peak_memory_bytes),
            elapsed_str
        );
    } else {
        let kmers_str = format_count(run_stats.n_kmers);
        if run_stats.chunks > 0 {
            warn!(
                "sharkmer complete: {} reads, {} kmers, {} chunks, allocator peak {}, {}",
                reads_str,
                kmers_str,
                run_stats.chunks,
                format_bytes(run_stats.peak_memory_bytes),
                elapsed_str
            );
        } else {
            warn!(
                "sharkmer complete: {} reads, {} kmers, allocator peak {}, {}",
                reads_str,
                kmers_str,
                format_bytes(run_stats.peak_memory_bytes),
                elapsed_str
            );
        }
    }
}

use crate::format::{format_bytes, format_count, format_duration};
use crate::io::write_fasta_fields;
use crate::kmer::FilteredKmerCounts;
use crate::kmer::budget::MemoryClass;
use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::text::BudgetString;
use crate::output::OutputTransaction;
use crate::pcr;
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, warn};
use rayon::prelude::*;
use serde::Serialize;

pub(crate) mod yaml;

/// Provenance for the uncertain subset of a target's combined FASTA.
#[derive(Serialize)]
pub(crate) struct UncertainProductOutput {
    pub(crate) n_products: usize,
    pub(crate) product_lengths: BudgetVec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_file: Option<BudgetString>,
    #[serde(flatten)]
    pub(crate) assessment: pcr::UncertainProducts,
}

/// Stats for a single PCR gene result, included in the YAML stats output.
#[derive(Serialize)]
pub(crate) struct PcrGeneResult {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) graph_buffers: Option<pcr::memory::admission::TargetGraphBuffers>,
    pub(crate) gene_name: BudgetString,
    pub(crate) status: &'static str,
    pub(crate) n_products: usize,
    pub(crate) n_standard_products: usize,
    pub(crate) n_uncertain_products: usize,
    #[serde(skip_serializing_if = "BudgetVec::is_empty")]
    pub(crate) product_lengths: BudgetVec<usize>,
    #[serde(skip_serializing_if = "BudgetVec::is_empty")]
    pub(crate) sequence_annotations: BudgetVec<pcr::SequenceAnnotation>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) output_file: Option<BudgetString>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) failure_reason: Option<BudgetString>,
    #[serde(skip_serializing_if = "BudgetVec::is_empty")]
    pub(crate) threshold_diagnostics: BudgetVec<pcr::PcrThresholdDiagnostic>,
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
    pub(crate) kind: &'static str,
    pub(crate) inputs: BudgetVec<BudgetString>,
    pub(crate) paired: bool,
    pub(crate) max_reads: u64,
}

/// Count identity and owned buffer measurements; these are not process RSS.
#[derive(Serialize)]
pub(crate) struct CountStoreSummary {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) process_file_limits: Option<crate::kmer::budget::ProcessFileLimits>,
    pub(crate) format_version: u32,
    pub(crate) backend: &'static str,
    pub(crate) hash_algorithm: &'static str,
    pub(crate) count_width_bits: u8,
    pub(crate) overflow_policy: &'static str,
    pub(crate) preprocessing: &'static str,
    pub(crate) canonicalization: &'static str,
    pub(crate) selected_input_sha256: BudgetString,
    pub(crate) distinct_keys: u64,
    pub(crate) occurrences: u64,
    pub(crate) owned_table_bytes: usize,
    pub(crate) peak_table_bytes: usize,
    pub(crate) peak_record_buffer_bytes: u64,
    pub(crate) peak_kmer_buffer_bytes: u64,
    pub(crate) peak_count_input_reserved_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) external: Option<ExternalCountSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) input_identity_worker: Option<crate::input_identity::IdentityWorkerSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) count_owner_transfer: Option<crate::io::CountOwnerTransferSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) adaptive_parallel: Option<AdaptiveParallelSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) parallel: Option<crate::kmer::parallel::ParallelMetrics>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) packed: Option<crate::kmer::packed::PackedMetrics>,
}

/// Actual storage of independently spillable counting owners.
#[derive(Serialize)]
pub(crate) struct AdaptiveParallelSummary {
    pub shards: crate::kmer::parallel::AdaptiveShardMetrics,
    pub peak_sort_buffer_bytes: u64,
    pub peak_temporary_disk_bytes: u64,
}

/// Measurements emitted only when exact external counting was performed.
#[derive(Serialize)]
pub(crate) struct ExternalCountSummary {
    pub(crate) count_file_bytes: u64,
    pub(crate) lookup_cache_bytes: usize,
    pub(crate) peak_count_page_bytes: u64,
    pub(crate) peak_sort_buffer_bytes: u64,
    // Shared temporary allowance during counting/finalization, including replay
    // when performed. Completed startup and later publication are excluded.
    pub(crate) peak_temporary_disk_bytes: u64,
    pub(crate) sorted_run_count: u64,
    pub(crate) merge_count: u64,
}

/// Optional all-entry equivalence evidence; absent unless explicitly requested.
#[derive(Serialize)]
pub(crate) struct CountFingerprint {
    pub(crate) algorithm: &'static str,
    pub(crate) distinct_keys: u64,
    pub(crate) occurrences: u64,
    pub(crate) pair_digest: BudgetString,
}

/// Structured run statistics, serialized as YAML.
#[derive(Serialize)]
pub(crate) struct RunStats {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) process_address_space: Option<crate::process_memory::ProcessAddressSpace>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) assembly_graph_buffers: Option<pcr::memory::admission::GraphBufferSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) read_replay: Option<crate::replay::ReplaySummary>,
    pub(crate) sharkmer_version: &'static str,
    pub(crate) command: BudgetString,
    pub(crate) sample: BudgetString,
    pub(crate) run_id: BudgetString,
    pub(crate) output_manifest: BudgetString,
    pub(crate) run_status: &'static str,
    pub(crate) product_output_schema_version: u32,
    pub(crate) product_output_layout: &'static str,
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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) count_store: Option<CountStoreSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) count_fingerprint: Option<CountFingerprint>,
    pub(crate) peak_memory_bytes: u64,
    pub(crate) stage_timings: StageTimings,
    #[serde(skip_serializing_if = "BudgetVec::is_empty")]
    pub(crate) pcr_results: BudgetVec<PcrGeneResult>,
}

/// Run in silico PCR for all primer pairs, write FASTA output files, and print results.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_pcr(
    kmer_counts: &FilteredKmerCounts,
    pcr_runs: &[pcr::PCRParams],
    sample: &str,
    directory: &str,
    min_kmer_count: u32,
    dump_graph: bool,
    show_progress: bool,
    reads: Option<crate::io::ThreadingSource<'_>>,
    max_nodes: usize,
    output_transaction: &mut OutputTransaction,
    graph_buffers: Option<&pcr::memory::admission::GraphBufferPool>,
) -> Result<BudgetVec<PcrGeneResult>> {
    let result_budget = graph_buffers.map(pcr::memory::admission::GraphBufferPool::result_budget);
    let mut pcr_results = BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);

    if pcr_runs.is_empty() {
        return Ok(pcr_results);
    }

    info!("Running in silico PCR...");
    let spinner_style = ProgressStyle::with_template("{spinner:.cyan} {msg}")
        .context("Invalid PCR progress style")?;
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
    let kmer_counts_pcr = kmer_counts.with_minimum(u64::from(min_kmer_count));

    // Run PCR for each gene in parallel; kmer_counts_pcr and reads are read-only and shared
    let pcr_fasta_results = BudgetVec::collect_indexed(
        pcr_runs
            .par_iter()
            .enumerate()
            .map(|(gene_index, pcr_params)| {
                let diagnostic_limits =
                    pcr::withheld_diagnostic_limits_for_gene(gene_index, pcr_runs.len());
                let fasta = graph_buffers
                    .map(|pool| pool.acquire())
                    .transpose()
                    .and_then(|lease| {
                        let result = if let Some(lease) = &lease {
                            pcr::do_pcr_with_graph_budget(
                                &kmer_counts_pcr,
                                sample,
                                pcr_params,
                                dump_graph,
                                directory,
                                reads,
                                max_nodes,
                                diagnostic_limits,
                                Some(lease.budget()),
                                result_budget,
                            )
                        } else {
                            pcr::do_pcr(
                                &kmer_counts_pcr,
                                sample,
                                pcr_params,
                                dump_graph,
                                directory,
                                reads,
                                max_nodes,
                                diagnostic_limits,
                            )
                        };
                        let summary = lease
                            .as_ref()
                            .map(pcr::memory::admission::TargetGraphLease::summary);
                        result.map(|outcome| (outcome, summary))
                    });
                (pcr_params, fasta)
            }),
        result_budget,
        MemoryClass::RetainedResults,
    )?;
    pcr_spinner.finish_and_clear();

    // Write output files sequentially to maintain deterministic order
    for (pcr_params, fasta_result) in pcr_fasta_results {
        let (outcome, graph_buffers) = fasta_result?;
        let pcr::PcrOutcome {
            records,
            mut sequence_annotations,
            failure_reason,
            threshold_diagnostics,
            uncertain_products,
        } = outcome;

        let n_standard_products = records.len();
        let fasta_basename = BudgetString::format(
            format_args!("{}_{}.fasta", sample, pcr_params.gene_name),
            result_budget,
            MemoryClass::RetainedResults,
        )?;
        let uncertain_products = uncertain_products
            .map(|assessment| -> Result<UncertainProductOutput> {
                let n_products = assessment.products.len();
                anyhow::ensure!(
                    n_products == assessment.unique_products,
                    "Uncertain product count differs from its retained records"
                );
                let product_lengths = BudgetVec::collect_in_class(
                    assessment
                        .products
                        .iter()
                        .map(|product| product.oriented_sequence.len()),
                    result_budget,
                    MemoryClass::RetainedResults,
                )?;
                let output_file = (n_products > 0)
                    .then(|| fasta_basename.try_clone())
                    .transpose()?;
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
        if let Some(uncertain) = &uncertain_products {
            for (index, product) in uncertain.assessment.products.iter().enumerate() {
                sequence_annotations.push(pcr::SequenceAnnotation {
                    sequence_id: BudgetString::from_str(
                        &product.id,
                        result_budget,
                        MemoryClass::RetainedResults,
                    )?,
                    sequence_index: index,
                    sequence_length: product.oriented_sequence.len(),
                    assembly_rank_score: None,
                })?;
            }
        }
        let mut product_lengths = BudgetVec::collect_in_class(
            records.iter().map(|record| record.seq().len()),
            result_budget,
            MemoryClass::RetainedResults,
        )?;
        if let Some(uncertain) = &uncertain_products {
            for &length in &uncertain.product_lengths {
                product_lengths.push(length)?;
            }
        }
        let output_file = if n_products > 0 {
            output_transaction.stage_file(&fasta_basename, |file| {
                for record in &records {
                    write_fasta_fields(file, record.id(), record.desc(), record.seq())
                        .context("Failed to write standard FASTA record")?;
                }
                if let Some(uncertain) = &uncertain_products {
                    for product in &uncertain.assessment.products {
                        write_fasta_fields(
                            file,
                            &product.id,
                            Some(&product.description),
                            product.oriented_sequence.as_bytes(),
                        )
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
            match failure_reason {
                Some(reason) => Some(reason),
                None => Some(BudgetString::from_str(
                    "unknown (no reason reported by PCR pipeline)",
                    result_budget,
                    MemoryClass::RetainedResults,
                )?),
            }
        };
        pcr_results.push(PcrGeneResult {
            graph_buffers,
            gene_name: BudgetString::from_str(
                &pcr_params.gene_name,
                result_budget,
                MemoryClass::RetainedResults,
            )?,
            status,
            n_products,
            n_standard_products,
            n_uncertain_products,
            product_lengths,
            sequence_annotations,
            output_file,
            failure_reason,
            threshold_diagnostics,
            uncertain_products,
        })?;
    }

    // Print gene result table
    let (sym_pass, sym_fail) = if show_progress {
        ("\u{2714}", "\u{2718}") // checkmark, x-mark
    } else {
        ("+", "-")
    };
    for result in &pcr_results {
        if result.n_products > 0 {
            let mut lengths = BudgetString::new(result_budget, MemoryClass::RetainedResults);
            for (index, length) in result.product_lengths.iter().enumerate() {
                if index > 0 {
                    lengths.push_str(", ")?;
                }
                let value = BudgetString::format(
                    format_args!("{length}"),
                    result_budget,
                    MemoryClass::RetainedResults,
                )?;
                lengths.push_str(&value)?;
            }
            warn!(
                "  {} {} ({} record{}: {} without assembly warnings, {} with assembly warnings; {} bp)",
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
                lengths
            );
            if result.n_uncertain_products > 0 {
                warn!(
                    "    Assembly warnings flag unresolved or suspected repeats; see assembly_status_causes in FASTA headers."
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
        if let Some(uncertain) = &result.uncertain_products {
            if uncertain.assessment.search_limits_reached {
                warn!(
                    "    Target search limits reached; additional sequences may remain unexplored."
                );
            }
            if uncertain.assessment.retention_truncated {
                warn!(
                    "    Target retention limits reached; additional sequence hypotheses may have been omitted."
                );
            }
            if uncertain.assessment.marker_detail_truncated {
                warn!(
                    "    Target diagnostic detail limits reached; some records retain summary evidence only."
                );
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
    let stats_basename = BudgetString::format(
        format_args!("{}.stats.yaml", run_stats.sample),
        run_stats.sample.budget(),
        MemoryClass::RetainedResults,
    )?;
    let yaml = yaml::prepare(run_stats)?;
    output_transaction.stage_file(&stats_basename, |file| {
        yaml.write(file).context("Failed to write stats YAML")
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
            "sharkmer complete: {} reads, {}/{} targets with output{}, {} records ({} without assembly warnings, {} with assembly warnings), allocator peak {}, {}",
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceBudget, ResourceLimits};
    use std::sync::Arc;
    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::TargetGraphs as usize] = 8 << 20;
        classes[MemoryClass::RetainedResults as usize] = bytes;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: 32 << 20,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }
    fn dna(mut state: u64, length: usize) -> String {
        (0..length)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[((state >> 32) & 3) as usize] as char
            })
            .collect()
    }
    fn fixture() -> (crate::kmer::KmerCounts, Vec<pcr::PCRParams>) {
        let mut counts = crate::kmer::KmerCounts::new(&19);
        let mut targets = Vec::new();
        for (name, sequence) in [
            ("repeat", dna(157, 60) + &"A".repeat(40) + &dna(158, 60)),
            ("clean", dna(159, 160)),
        ] {
            for _ in 0..4 {
                counts.ingest_seq(&sequence).unwrap();
            }
            let reverse = String::from_utf8(bio::alphabets::dna::revcomp(
                &sequence.as_bytes()[sequence.len() - 15..],
            ))
            .unwrap();
            let mut params=crate::cli::parse_pcr_primers_string(&format!("name={name},forward={},reverse={reverse},trim=15,mismatches=0,min-length=100,max-length=220,dedup-edit-threshold=0",&sequence[..15])).unwrap();
            params.emit_uncertain = true;
            params.diagnose_withheld_paths = true;
            targets.push(params);
        }
        let mut missing=crate::cli::parse_pcr_primers_string("name=missing,forward=CCCCCCCCCCCCCCC,reverse=GGGGGGGGGGGGGGG,trim=15,mismatches=0,min-length=100,max-length=220").unwrap();
        missing.emit_uncertain = true;
        targets.push(missing);
        (counts, targets)
    }
    #[test]
    fn target_stats_keep_owned_metadata_through_serialization_after_workers_finish() {
        let (counts, targets) = fixture();
        let budget = budget(1 << 20);
        let pool =
            pcr::memory::admission::GraphBufferPool::new(Arc::clone(&budget), 8 << 20).unwrap();
        let directory = tempfile::tempdir().unwrap();
        let mut transaction = OutputTransaction::begin(directory.path(), "stats").unwrap();
        let workers = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let results = workers
            .install(|| {
                run_pcr(
                    &counts.filtered_view(1),
                    &targets,
                    "stats",
                    directory.path().to_str().unwrap(),
                    1,
                    false,
                    false,
                    None,
                    100_000,
                    &mut transaction,
                    Some(&pool),
                )
            })
            .unwrap();
        assert_eq!(results.len(), 3);
        assert!(results[0].n_uncertain_products > 0);
        assert_eq!(results[1].n_standard_products, 1);
        assert_eq!(results[2].status, "fail");
        assert_eq!(pool.summary().unwrap().peak_active_targets, 1);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        assert!(budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize] > 0);
        struct CheckedWriter<'a> {
            data: Vec<u8>,
            budget: &'a ResourceBudget,
        }
        impl std::io::Write for CheckedWriter<'_> {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                assert!(
                    self.budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize]
                        > 0
                );
                self.data.extend_from_slice(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let mut writer = CheckedWriter {
            data: Vec::new(),
            budget: &budget,
        };
        yaml::Prepared::new(&results, Some(&budget))
            .unwrap()
            .write(&mut writer)
            .unwrap();
        let value: serde_yaml_ng::Value = serde_yaml_ng::from_slice(&writer.data).unwrap();
        assert_eq!(value[0]["gene_name"].as_str(), Some("repeat"));
        assert_eq!(value[1]["product_lengths"][0].as_u64(), Some(160));
        assert_eq!(
            value[2]["failure_reason"].as_str(),
            Some("forward and reverse primers not found")
        );
        drop(results);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
    #[test]
    fn stats_allocation_refusal_releases_collected_and_unconsumed_target_results() {
        let (counts, targets) = fixture();
        let workers = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let (mut refused, mut succeeded) = (0, 0);
        for bytes in [
            1,
            64,
            256,
            512,
            1024,
            2048,
            4096,
            8192,
            16384,
            32768,
            65536,
            1 << 20,
        ] {
            let budget = budget(bytes);
            let pool =
                pcr::memory::admission::GraphBufferPool::new(Arc::clone(&budget), 8 << 20).unwrap();
            let directory = tempfile::tempdir().unwrap();
            let mut transaction = OutputTransaction::begin(directory.path(), "stats").unwrap();
            match workers.install(|| {
                run_pcr(
                    &counts.filtered_view(1),
                    &targets,
                    "stats",
                    directory.path().to_str().unwrap(),
                    1,
                    false,
                    false,
                    None,
                    100_000,
                    &mut transaction,
                    Some(&pool),
                )
            }) {
                Ok(results) => {
                    assert_eq!(results.len(), 3);
                    drop(results);
                    succeeded += 1;
                }
                Err(error) => {
                    assert!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some(),
                        "{error:#}"
                    );
                    refused += 1;
                }
            }
            let summary = pool.summary().unwrap();
            if bytes == 1 {
                assert_eq!(summary.peak_active_targets, 0);
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0, "budget {bytes}");
            drop(transaction);
            assert!(!directory.path().join("stats_repeat.fasta").exists());
            assert!(!directory.path().join("stats_clean.fasta").exists());
        }
        assert!(refused > 4 && succeeded > 0);
    }
}

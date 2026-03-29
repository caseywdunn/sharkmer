// pcr/paths.rs — path finding, sequence extraction, dedup

use anyhow::{Context, Result};
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use log::debug;
use petgraph::algo::all_simple_paths;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;

use crate::kmer::FilteredKmerCounts;

use super::graph::{compute_mean, compute_median, get_end_nodes, get_start_nodes};
use super::{AssemblyRecord, DBEdge, DBNode, PCRParams};

/// In cases where a primer pair has many possible paths, limit the number of paths to consider
const MAX_NUM_PATHS_PER_PAIR: usize = 20;

/// The maximum number of fasta records to return
const MAX_NUM_AMPLICONS: usize = 20;

// Given the length of a sequence in bp, return the number of kmers needed to cover it
// If the length is 0, return 0
// If the length is greater than 0 and less than or equal to k, return 1
fn bp_length_to_kmer_length(bp_length: usize, k: usize) -> usize {
    if bp_length == 0 {
        0
    } else if bp_length <= k {
        1
    } else {
        (bp_length - k) + 1
    }
}

pub fn get_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
) -> Vec<Vec<NodeIndex>> {
    let mut all_paths = Vec::new();

    for start in get_start_nodes(graph) {
        for end in get_end_nodes(graph) {
            let paths_for_this_pair =
                all_simple_paths::<Vec<NodeIndex>, &StableDiGraph<DBNode, DBEdge>>(
                    graph,
                    start,
                    end,
                    bp_length_to_kmer_length(params.min_length, kmer_counts.get_k()),
                    Some(bp_length_to_kmer_length(
                        params.max_length,
                        kmer_counts.get_k(),
                    )),
                );

            // Limit the number of paths to consider to MAX_NUM_PATHS_PER_PAIR
            all_paths.extend(
                paths_for_this_pair
                    .take(MAX_NUM_PATHS_PER_PAIR)
                    .collect::<Vec<Vec<NodeIndex>>>(),
            );
        }
    }

    all_paths
}

/// Extract sequences from graph paths, producing FASTA assembly records.
/// Returns the generated records and the updated amplicon index.
pub(super) fn generate_sequences_from_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    all_paths: Vec<Vec<NodeIndex>>,
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    mut amplicon_index: usize,
) -> Result<(Vec<AssemblyRecord>, usize)> {
    let mut assembly_records: Vec<AssemblyRecord> = Vec::new();

    for path in all_paths.into_iter() {
        let mut sequence = String::new();
        let mut edge_counts: Vec<u64> = Vec::new(); // u64 for compute_mean/median compatibility
        let mut parent_node: NodeIndex = NodeIndex::new(0);
        for node in path.iter() {
            let node_data = graph
                .node_weight(*node)
                .context("Node not found in graph during sequence generation")?;
            let subread = crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(kmer_counts.get_k() - 1));
            if sequence.is_empty() {
                sequence = subread;
                parent_node = *node;
            } else {
                let last_char = subread
                    .chars()
                    .last()
                    .context("Empty subread during sequence generation")?;
                sequence.push(last_char);

                let edge = graph
                    .find_edge(parent_node, *node)
                    .context("Edge not found between path nodes")?;
                let edge_data = graph.edge_weight(edge).context("Edge weight not found")?;
                edge_counts.push(edge_data.count as u64);
                parent_node = *node;
            }
        }

        if sequence.len() < params.min_length {
            debug!(
                "  Path length is {} bp, shorter than min-length {}. Skipping.",
                sequence.len(),
                params.min_length
            );
            continue;
        }

        let count_mean = compute_mean(&edge_counts);
        let count_median = compute_median(&edge_counts);
        let count_min = edge_counts
            .iter()
            .min()
            .context("No edge counts found for path")?;
        let count_max = edge_counts
            .iter()
            .max()
            .context("No edge counts found for path")?;

        let id = format!("{}_{}_{}", sample_name, params.gene_name, amplicon_index);
        let desc = format!(
            "sample={} gene={} product={} length={} kmer_count_mean={:.2} kmer_count_median={} kmer_count_min={} kmer_count_max={}",
            sample_name,
            params.gene_name,
            amplicon_index,
            sequence.len(),
            count_mean,
            count_median,
            count_min,
            count_max
        );

        amplicon_index += 1;

        debug!(">{} {}", id, desc);
        let record = fasta::Record::with_attrs(&id, Some(&desc), sequence.as_bytes());
        assembly_records.push(AssemblyRecord {
            fasta_record: record,
            kmer_min_count: *count_min as u32,
        });
    }

    Ok((assembly_records, amplicon_index))
}

pub(super) fn pairwise_sequence_distances(
    records: &[fasta::Record],
    threshold: u32,
) -> Vec<Vec<Option<u32>>> {
    // https://docs.rs/bio/latest/bio/alignment/distance/simd/fn.bounded_levenshtein.html

    let n = records.len();
    let mut matrix: Vec<Vec<Option<u32>>> = vec![vec![Some(0); n]; n];

    for i in 0..n {
        for j in i + 1..n {
            let dist = bounded_levenshtein(records[i].seq(), records[j].seq(), threshold);
            matrix[i][j] = dist;
            matrix[j][i] = dist; // Symmetric matrix
        }
    }

    matrix
}

/// Sort assembly records by descending kmer_min_count, deduplicate near-identical
/// sequences, and truncate to the maximum number of amplicons.
pub(super) fn sort_and_deduplicate(
    assembly_records: Vec<AssemblyRecord>,
    params: &PCRParams,
) -> Vec<fasta::Record> {
    let mut sorted = assembly_records;
    sorted.sort_by(|a, b| {
        b.kmer_min_count
            .cmp(&a.kmer_min_count)
            .then_with(|| a.fasta_record.seq().cmp(b.fasta_record.seq()))
    });

    let mut records: Vec<fasta::Record> = sorted.into_iter().map(|ar| ar.fasta_record).collect();

    let num_records_all = records.len();

    // Greedy clustering: keep each record only if it is not within
    // dedup_edit_threshold edits of any previously kept record.
    let distances = pairwise_sequence_distances(&records, params.dedup_edit_threshold);
    let n = records.len();
    let mut keep = vec![true; n];
    for i in 0..n {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..n {
            if !keep[j] {
                continue;
            }
            if distances[i][j].is_some() {
                keep[j] = false;
            }
        }
    }
    let mut idx = 0;
    records.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });

    if num_records_all == records.len() {
        gene_info!(
            params.gene_name,
            "{} PCR products were generated and retained.",
            num_records_all
        );
    } else {
        gene_info!(
            params.gene_name,
            "{} PCR products were generated and {} were retained ({} removed as near-duplicates within edit distance {}).",
            num_records_all,
            records.len(),
            num_records_all - records.len(),
            params.dedup_edit_threshold
        );
    }

    if records.len() > MAX_NUM_AMPLICONS {
        gene_warn!(
            params.gene_name,
            "There are {} PCR products. This exceeds the maximum of {}. Retaining only the first {} records.",
            records.len(),
            MAX_NUM_AMPLICONS,
            MAX_NUM_AMPLICONS
        );
        records.truncate(MAX_NUM_AMPLICONS);
    }

    records
}

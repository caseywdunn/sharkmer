// pcr/paths.rs — path finding, sequence extraction, dedup

use anyhow::{Context, Result};
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use log::debug;
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;

use crate::kmer::FilteredKmerCounts;

use super::graph::{compute_mean, compute_median, get_end_nodes, get_start_nodes};
use super::{AssemblyRecord, DBEdge, DBNode, PCRParams};

/// Budget: maximum paths to find per start-end pair
const MAX_NUM_PATHS_PER_PAIR: usize = 20;

/// The maximum number of fasta records to return
const MAX_NUM_AMPLICONS: usize = 20;

/// Maximum number of times a node can appear in a single path.
/// 1 = no revisiting (classic simple path). 2 = allow one repeat
/// traversal (handles tandem duplications).
const MAX_NODE_VISITS: usize = 2;

/// Maximum number of DFS states to explore per start node before
/// giving up. Prevents combinatorial explosion in complex graphs.
const MAX_DFS_STATES: usize = 100_000;

/// Find paths from start nodes to end nodes using coverage-weighted DFS.
/// At each branch point, outgoing edges are explored in descending order
/// of edge count, so the highest-coverage paths are found first.
pub fn get_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
) -> Vec<Vec<NodeIndex>> {
    let min_path_nodes = if params.min_length <= kmer_counts.get_k() {
        1
    } else {
        (params.min_length - kmer_counts.get_k()) + 1
    };
    let max_path_nodes = if params.max_length <= kmer_counts.get_k() {
        1
    } else {
        (params.max_length - kmer_counts.get_k()) + 1
    };

    let end_nodes: std::collections::HashSet<NodeIndex> =
        get_end_nodes(graph).into_iter().collect();
    let mut all_paths = Vec::new();

    for start in get_start_nodes(graph) {
        let mut paths_from_start = 0;
        let mut states_explored: usize = 0;

        // DFS stack: (current_path, visit_counts)
        // visit_counts tracks how many times each node appears in the path
        let mut stack: Vec<(Vec<NodeIndex>, std::collections::HashMap<NodeIndex, usize>)> =
            Vec::new();
        let mut initial_visits = std::collections::HashMap::new();
        initial_visits.insert(start, 1);
        stack.push((vec![start], initial_visits));

        while let Some((path, visit_counts)) = stack.pop() {
            states_explored += 1;
            if paths_from_start >= MAX_NUM_PATHS_PER_PAIR || states_explored >= MAX_DFS_STATES {
                break;
            }

            let current = *path.last().unwrap();
            let path_len = path.len();

            // Check if we reached an end node with valid length
            if end_nodes.contains(&current) && path_len >= min_path_nodes {
                all_paths.push(path.clone());
                paths_from_start += 1;
                // Don't continue extending past end nodes
                continue;
            }

            // Don't extend past max length
            if path_len >= max_path_nodes {
                continue;
            }

            // Get outgoing edges sorted by count (descending) for coverage-weighted exploration
            let mut outgoing: Vec<(NodeIndex, u32)> = graph
                .edges_directed(current, Direction::Outgoing)
                .map(|e| (e.target(), e.weight().count))
                .collect();
            // Sort ascending so that when pushed to stack, highest-count is popped first
            outgoing.sort_by(|a, b| a.1.cmp(&b.1));

            for (neighbor, _count) in outgoing {
                let current_visits = visit_counts.get(&neighbor).copied().unwrap_or(0);
                if current_visits < MAX_NODE_VISITS {
                    let mut new_path = path.clone();
                    new_path.push(neighbor);
                    let mut new_visits = visit_counts.clone();
                    *new_visits.entry(neighbor).or_insert(0) += 1;
                    stack.push((new_path, new_visits));
                }
            }
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

        // Compute coverage consistency (coefficient of variation)
        let coverage_cv = if count_mean > 0.0 {
            let variance = edge_counts
                .iter()
                .map(|&c| {
                    let diff = c as f64 - count_mean;
                    diff * diff
                })
                .sum::<f64>()
                / edge_counts.len() as f64;
            variance.sqrt() / count_mean
        } else {
            0.0
        };

        // Compute max coverage ratio from edge annotations
        let max_coverage_ratio = {
            let mut max_ratio: f64 = 0.0;
            let mut parent = path[0];
            for &node in &path[1..] {
                if let Some(edge_idx) = graph.find_edge(parent, node) {
                    let ratio = graph
                        .edge_weight(edge_idx)
                        .map_or(0.0, |e| e.coverage_ratio);
                    if ratio > max_ratio {
                        max_ratio = ratio;
                    }
                }
                parent = node;
            }
            max_ratio
        };

        let score = super::PathScore {
            kmer_min_count: *count_min as u32,
            kmer_median_count: count_median,
            coverage_cv,
            max_coverage_ratio,
        };

        let id = format!("{}_{}_{}", sample_name, params.gene_name, amplicon_index);
        let desc = format!(
            "sample={} gene={} product={} length={} kmer_count_mean={:.2} kmer_count_median={} kmer_count_min={} kmer_count_max={} score={:.2}",
            sample_name,
            params.gene_name,
            amplicon_index,
            sequence.len(),
            count_mean,
            count_median,
            count_min,
            count_max,
            score.composite()
        );

        amplicon_index += 1;

        debug!(">{} {}", id, desc);
        let record = fasta::Record::with_attrs(&id, Some(&desc), sequence.as_bytes());
        assembly_records.push(AssemblyRecord {
            fasta_record: record,
            score,
        });
    }

    Ok((assembly_records, amplicon_index))
}

/// Sort assembly records by descending composite score, deduplicate near-identical
/// sequences, and truncate to the maximum number of amplicons.
pub(super) fn sort_and_deduplicate(
    assembly_records: Vec<AssemblyRecord>,
    params: &PCRParams,
) -> Vec<fasta::Record> {
    let mut sorted = assembly_records;
    sorted.sort_by(|a, b| {
        b.score
            .composite()
            .partial_cmp(&a.score.composite())
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.fasta_record.seq().cmp(b.fasta_record.seq()))
    });

    let records: Vec<fasta::Record> = sorted.into_iter().map(|ar| ar.fasta_record).collect();
    let num_records_all = records.len();

    // Greedy clustering: keep each record only if it is not within
    // dedup_edit_threshold edits of any previously kept record.
    // Distances are computed on-the-fly to avoid O(N^2) memory.
    let mut kept: Vec<fasta::Record> = Vec::new();
    for record in records {
        let is_duplicate = kept.iter().any(|kept_record| {
            bounded_levenshtein(record.seq(), kept_record.seq(), params.dedup_edit_threshold)
                .is_some()
        });
        if !is_duplicate {
            kept.push(record);
        }
    }
    let mut records = kept;

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

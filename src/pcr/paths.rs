// pcr/paths.rs — path finding, sequence extraction, dedup

use anyhow::{Context, Result};
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use log::debug;
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;

use crate::kmer::FilteredKmerCounts;

use super::graph::{compute_mean, compute_median, get_end_nodes, get_start_nodes};
use super::{AssemblyRecord, DBEdge, DBNode, PCRParams};

/// The maximum number of fasta records to return
const MAX_NUM_AMPLICONS: usize = 20;

/// Get outgoing edges sorted by score (ascending so pop gives highest first).
fn sorted_children(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    edge_preferences: Option<&std::collections::HashMap<EdgeIndex, f64>>,
) -> Vec<(NodeIndex, f64)> {
    let mut outgoing: Vec<(NodeIndex, f64)> = graph
        .edges_directed(node, Direction::Outgoing)
        .map(|e| {
            let base_score = e.weight().count as f64;
            let pref = edge_preferences
                .and_then(|prefs| prefs.get(&e.id()))
                .copied()
                .unwrap_or(1.0);
            (e.target(), base_score * pref)
        })
        .collect();
    outgoing.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    outgoing
}

/// Find paths from start nodes to end nodes using coverage-weighted DFS.
/// At each branch point, outgoing edges are explored in descending order
/// of edge count, so the highest-coverage paths are found first.
/// When `edge_preferences` are provided (from bubble resolution), they
/// boost the ordering of read-supported edges at branch points.
pub fn get_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
    edge_preferences: Option<&std::collections::HashMap<EdgeIndex, f64>>,
) -> Vec<Vec<NodeIndex>> {
    // A path of N nodes produces a sequence of (k-1) + (N-1) = N+k-2 bases
    // (first node contributes k-1 bases via its sub_kmer, each subsequent
    // node extends by one base). Inverting: N = seq_length - k + 2.
    //
    // min_path_nodes is used only as a pre-filter during DFS; the exact
    // lower bound is re-checked in `generate_sequences_from_paths` against
    // the actual sequence length.
    //
    // max_path_nodes is the authoritative upper bound — the DFS caps paths
    // at this value and there is no secondary check downstream, so an
    // off-by-one here directly rejects valid amplicons at the user's
    // declared max-length.
    let k = kmer_counts.get_k();
    let min_path_nodes = if params.min_length <= k {
        1
    } else {
        params.min_length - k + 2
    };
    let max_path_nodes = if params.max_length <= k {
        1
    } else {
        params.max_length - k + 2
    };

    let end_nodes: std::collections::HashSet<NodeIndex> =
        get_end_nodes(graph).into_iter().collect();
    let mut all_paths = Vec::new();

    for start in get_start_nodes(graph) {
        let mut paths_from_start = 0;
        let mut states_explored: usize = 0;

        // Stack-based DFS with push/pop backtracking instead of cloning.
        // `path` and `visit_counts` are maintained incrementally.
        // `child_stack` holds remaining children to explore at each depth.
        let mut path: Vec<NodeIndex> = vec![start];
        let mut visit_counts: std::collections::HashMap<NodeIndex, usize> =
            std::collections::HashMap::new();
        visit_counts.insert(start, 1);

        // Compute sorted children for the start node
        let children = sorted_children(graph, start, edge_preferences);
        let mut child_stack: Vec<Vec<(NodeIndex, f64)>> = vec![children];

        loop {
            if paths_from_start >= params.max_paths_per_pair
                || states_explored >= params.max_dfs_states
            {
                break;
            }

            let depth = child_stack.len() - 1;

            if let Some((neighbor, _score)) = child_stack[depth].pop() {
                states_explored += 1;

                let current_visits = visit_counts.get(&neighbor).copied().unwrap_or(0);
                if current_visits >= params.max_node_visits {
                    continue;
                }

                // Push neighbor onto path
                path.push(neighbor);
                *visit_counts.entry(neighbor).or_insert(0) += 1;

                let path_len = path.len();

                // Check if we reached an end node with valid length
                if end_nodes.contains(&neighbor) && path_len >= min_path_nodes {
                    all_paths.push(path.clone());
                    paths_from_start += 1;
                    // Backtrack: undo push
                    *visit_counts.get_mut(&neighbor).unwrap() -= 1;
                    path.pop();
                    continue;
                }

                // Don't extend past max length
                if path_len >= max_path_nodes {
                    *visit_counts.get_mut(&neighbor).unwrap() -= 1;
                    path.pop();
                    continue;
                }

                // Explore deeper: push children frame
                let children = sorted_children(graph, neighbor, edge_preferences);
                child_stack.push(children);
            } else {
                // No more children at this depth, backtrack
                child_stack.pop();
                if child_stack.is_empty() {
                    break;
                }
                let backtrack_node = path.pop().expect("BUG: path empty during DFS backtrack");
                *visit_counts
                    .get_mut(&backtrack_node)
                    .expect("BUG: backtrack node missing from visit_counts") -= 1;
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
    threading: Option<&super::threading::ThreadingAnnotations>,
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
            if sequence.is_empty() {
                sequence =
                    crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(kmer_counts.get_k() - 1));
                parent_node = *node;
            } else {
                sequence.push(crate::kmer::kmer_last_base(&node_data.sub_kmer));

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

        if edge_counts.is_empty() {
            debug!("  Path has no edges (single node). Skipping.");
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

        // Compute read-support metrics from threading annotations
        let (zero_support_edges, median_unambiguous_support, edge_support_fraction) =
            if let Some(ann) = threading {
                let mut total_edges = 0u32;
                let mut supported_edges = 0u32;
                let mut zero_count = 0u32;
                let mut unambiguous_counts: Vec<u64> = Vec::new();

                let mut parent = path[0];
                for &node in &path[1..] {
                    if let Some(edge_idx) = graph.find_edge(parent, node) {
                        total_edges += 1;
                        match ann.edge_support.get(&edge_idx) {
                            Some(s) if s.read_support_total > 0 => {
                                supported_edges += 1;
                                unambiguous_counts.push(s.read_support_unambiguous as u64);
                            }
                            _ => {
                                zero_count += 1;
                                unambiguous_counts.push(0);
                            }
                        }
                    }
                    parent = node;
                }

                let frac = if total_edges > 0 {
                    supported_edges as f64 / total_edges as f64
                } else {
                    0.0
                };
                let median_unamb = if unambiguous_counts.is_empty() {
                    0.0
                } else {
                    compute_median(&unambiguous_counts)
                };

                (Some(zero_count), Some(median_unamb), Some(frac))
            } else {
                (None, None, None)
            };

        let score = super::PathScore {
            kmer_min_count: *count_min as u32,
            kmer_median_count: count_median,
            coverage_cv,
            max_coverage_ratio,
            zero_support_edges,
            median_unambiguous_support,
            edge_support_fraction,
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
            // Byte-level sequence comparison as a deterministic tiebreaker
            // when two records tie on composite score. Not biologically
            // meaningful — its only purpose is to make dedup order stable
            // across runs (otherwise iteration order would depend on
            // HashMap seeds deeper in path enumeration).
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

#[cfg(test)]
mod tests {
    use super::super::{
        DBEdge, DBNode, DEFAULT_DEDUP_EDIT_THRESHOLD, DEFAULT_HIGH_COVERAGE_RATIO,
        DEFAULT_MAX_DFS_STATES, DEFAULT_MAX_NODE_VISITS, DEFAULT_MAX_NUM_PRIMER_KMERS,
        DEFAULT_MAX_PATHS_PER_PAIR, DEFAULT_MIN_COMPONENT_BUDGET, DEFAULT_TIP_COVERAGE_FRACTION,
        StoppingCriteria,
    };
    use super::*;

    fn mk_node(sub_kmer: u64, is_start: bool, is_end: bool) -> DBNode {
        DBNode {
            sub_kmer,
            is_start,
            is_end,
            is_terminal: false,
            visited: false,
        }
    }

    fn mk_edge(count: u32) -> DBEdge {
        DBEdge {
            count,
            coverage_ratio: 1.0,
        }
    }

    fn test_params(min_length: usize, max_length: usize) -> PCRParams {
        PCRParams {
            forward_seq: "ACGT".to_string(),
            reverse_seq: "TGCA".to_string(),
            min_length,
            max_length,
            gene_name: "test".to_string(),
            min_count: 2,
            mismatches: 0,
            trim: 0,
            expected_length: None,
            citation: String::new(),
            notes: String::new(),
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            source: "test".to_string(),
            max_dfs_states: DEFAULT_MAX_DFS_STATES,
            max_paths_per_pair: DEFAULT_MAX_PATHS_PER_PAIR,
            max_node_visits: DEFAULT_MAX_NODE_VISITS,
            max_primer_kmers: DEFAULT_MAX_NUM_PRIMER_KMERS,
            high_coverage_ratio: DEFAULT_HIGH_COVERAGE_RATIO,
            tip_coverage_fraction: DEFAULT_TIP_COVERAGE_FRACTION,
            stopping_criteria: StoppingCriteria::AllComponents,
            min_component_budget: DEFAULT_MIN_COMPONENT_BUDGET,
        }
    }

    /// Linear graph: start -> a -> b -> end. Should find exactly one path.
    #[test]
    fn test_linear_path() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let e = graph.add_node(mk_node(3, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, e, mk_edge(10));

        // k=3, path of 4 nodes = 4 + 3 - 2 = 5 bases
        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0], vec![s, a, b, e]);
    }

    /// Diamond graph: start -> {a, b} -> end. Should find two paths.
    #[test]
    fn test_diamond_finds_both_paths() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let e = graph.add_node(mk_node(3, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(s, b, mk_edge(5));
        graph.add_edge(a, e, mk_edge(10));
        graph.add_edge(b, e, mk_edge(5));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert_eq!(paths.len(), 2);
    }

    /// Empty graph should produce no paths.
    #[test]
    fn test_no_start_nodes_gives_empty() {
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert!(paths.is_empty());
    }

    /// max_path_nodes caps path length correctly.
    #[test]
    fn test_max_length_caps_paths() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let c = graph.add_node(mk_node(3, false, false));
        let e = graph.add_node(mk_node(4, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, c, mk_edge(10));
        graph.add_edge(c, e, mk_edge(10));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        // 5 nodes = 5 + 3 - 2 = 6 bases. Set max_length to 5 (too short).
        // max_path_nodes = 5 - 3 + 2 = 4. Path needs 5 nodes, so no paths found.
        let params = test_params(0, 5);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert!(paths.is_empty());
    }

    /// DFS state budget limits exploration.
    #[test]
    fn test_dfs_budget_limits_exploration() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let e = graph.add_node(mk_node(2, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, e, mk_edge(10));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let mut params = test_params(0, 100);
        params.max_dfs_states = 0; // no exploration allowed
        let paths = get_assembly_paths(&graph, &fkc, &params, None);
        assert!(paths.is_empty());
    }

    /// sorted_children returns edges in ascending score order (so pop gives highest).
    #[test]
    fn test_sorted_children_order() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let lo = graph.add_node(mk_node(1, false, false));
        let hi = graph.add_node(mk_node(2, false, false));
        graph.add_edge(s, lo, mk_edge(1));
        graph.add_edge(s, hi, mk_edge(100));

        let children = sorted_children(&graph, s, None);
        assert_eq!(children.len(), 2);
        // Ascending: low first, high second (so pop gives high)
        assert_eq!(children[0].0, lo);
        assert_eq!(children[1].0, hi);
    }
}

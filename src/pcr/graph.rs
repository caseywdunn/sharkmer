// pcr/graph.rs — graph data structures and construction

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use log::debug;
// Used only by `summarize_extension` for a single debug-level diagnostic
// about whether the extended graph contains cycles. Not load-bearing for
// correctness — if this import ever becomes unused again, the summarize
// call has been removed and both can go.
use petgraph::algo::is_cyclic_directed;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use std::collections::VecDeque;

use crate::kmer::FilteredKmerCounts;
use crate::kmer::KmerCounts;

use super::{DBEdge, DBNode, PCRParams};

/// How frequently to check extension graph for ballooning growth
const EXTENSION_EVALUATION_FREQUENCY: usize = 1_000;

/// The graph depth over which to evaluate ballooning growth
pub(super) const EXTENSION_EVALUATION_DEPTH: usize = 4;

/// Maximum global node budget (used at high data volumes)
pub const DEFAULT_MAX_NUM_NODES: usize = 500_000;

/// Minimum global node budget (used at low data volumes)
const MIN_NODE_BUDGET: usize = 100_000;

/// Bases ingested below which the minimum budget is used (~1M × 150bp reads)
const BUDGET_LERP_LOW_BP: u64 = 150_000_000;

/// Bases ingested above which the maximum budget is used (~5M × 150bp reads)
const BUDGET_LERP_HIGH_BP: u64 = 750_000_000;

/// Compute the global node budget based on data volume.
///
/// Lerps from 100K nodes (≤150M bp, ~1M reads at 150bp) to 500K nodes
/// (≥750M bp, ~5M reads at 150bp). At low coverage, a smaller budget
/// avoids wasting time on off-target components that can't produce
/// products anyway. At high coverage, more seeds survive and need
/// more budget to be tried.
pub fn compute_node_budget(n_bases_ingested: u64) -> usize {
    if n_bases_ingested <= BUDGET_LERP_LOW_BP {
        MIN_NODE_BUDGET
    } else if n_bases_ingested >= BUDGET_LERP_HIGH_BP {
        DEFAULT_MAX_NUM_NODES
    } else {
        let fraction = (n_bases_ingested - BUDGET_LERP_LOW_BP) as f64
            / (BUDGET_LERP_HIGH_BP - BUDGET_LERP_LOW_BP) as f64;
        let budget =
            MIN_NODE_BUDGET as f64 + fraction * (DEFAULT_MAX_NUM_NODES - MIN_NODE_BUDGET) as f64;
        budget as usize
    }
}

// HIGH_COVERAGE_RATIO_THRESHOLD is now read from params.high_coverage_ratio

// A mask that can be used to isolate the last k-1 nucleotides of a kmer
pub(super) fn get_suffix_mask(k: &usize) -> u64 {
    debug_assert!(*k > 0 && *k <= 32, "k must be in 1..=32, got {}", k);
    (1 << (2 * (*k - 1))) - 1
}

/// Compute the median of `counts` as an f64, using `select_nth_unstable`.
///
/// This is expected-linear time, vs `sort()`'s O(n log n) — a measurable
/// difference because `compute_median_edge_count` runs every
/// `EXTENSION_EVALUATION_FREQUENCY` nodes during graph extension, where
/// `counts.len()` grows with the edge count of the working graph (up to
/// hundreds of thousands of edges on dense panels).
///
/// Returns `None` if `counts` is empty. The even-length case avoids the
/// extra full sort: after partitioning at `mid`, the smaller half is
/// `counts[..mid]` (any order, but bounded above by `counts[mid]`), so a
/// single linear scan over that half finds the (mid-1)-th order statistic.
pub(super) fn median_u32_f64(mut counts: Vec<u32>) -> Option<f64> {
    if counts.is_empty() {
        return None;
    }
    let mid = counts.len() / 2;
    if counts.len() % 2 == 0 {
        counts.select_nth_unstable(mid);
        let upper_min = counts[mid] as f64;
        let lower_max = *counts[..mid]
            .iter()
            .max()
            .expect("non-empty lower half when len >= 2") as f64;
        Some((lower_max + upper_min) / 2.0)
    } else {
        Some(*counts.select_nth_unstable(mid).1 as f64)
    }
}

fn compute_median_edge_count(graph: &StableDiGraph<DBNode, DBEdge>, default: f64) -> f64 {
    let counts: Vec<u32> = graph.edge_indices().map(|e| graph[e].count).collect();
    median_u32_f64(counts).unwrap_or(default)
}

pub fn get_dbedge(kmer: &u64, kmer_counts: &FilteredKmerCounts) -> DBEdge {
    DBEdge {
        count: kmer_counts.get_canonical_count(kmer),
        coverage_ratio: 0.0, // set during annotation pass
    }
}

/// Reconstruct the full kmer for an edge from its source and target node sub_kmers.
/// The edge kmer is: (source.sub_kmer << 2) | (target.sub_kmer & 3).
pub(super) fn reconstruct_edge_kmer(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edge_idx: petgraph::graph::EdgeIndex,
) -> u64 {
    // Safe: callers pass edge indices from graph.edge_indices() / edge_references()
    let (src, tgt) = graph.edge_endpoints(edge_idx).unwrap();
    (graph[src].sub_kmer << 2) | (graph[tgt].sub_kmer & 3)
}

pub fn get_start_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    graph
        .node_indices()
        .filter(|&node| graph[node].is_start)
        .collect()
}

pub fn get_end_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    graph
        .node_indices()
        .filter(|&node| graph[node].is_end)
        .collect()
}

// Find the descendants of a node, each descendent no more than `depth` edges away
pub(super) fn descendants(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    depth: usize,
) -> AHashSet<NodeIndex> {
    let mut visited: AHashSet<NodeIndex> = AHashSet::new();
    let mut queue = VecDeque::new();

    queue.push_back((node, 0));
    visited.insert(node);

    let mut descendants: AHashSet<NodeIndex> = AHashSet::new();

    while let Some((current_node, current_depth)) = queue.pop_front() {
        if current_depth >= depth {
            continue;
        }

        for neighbor in graph.neighbors(current_node) {
            if visited.insert(neighbor) {
                // insert returns false if the item was already in the set
                descendants.insert(neighbor);
                queue.push_back((neighbor, current_depth + 1));
            }
        }
    }

    descendants
}

pub fn compute_mean(numbers: &[u64]) -> f64 {
    if numbers.is_empty() {
        return 0.0;
    }
    let sum: u64 = numbers.iter().sum();
    sum as f64 / numbers.len() as f64
}

pub fn compute_median(numbers: &[u64]) -> f64 {
    if numbers.is_empty() {
        return 0.0;
    }
    let mut sorted = numbers.to_vec();
    sorted.sort();

    let mid = sorted.len() / 2;

    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) as f64 / 2.0
    } else {
        sorted[mid] as f64
    }
}

pub(super) fn create_seed_graph(
    forward_primer_kmers: &KmerCounts,
    reverse_primer_kmers: &KmerCounts,
    kmer_counts: &FilteredKmerCounts,
) -> (StableDiGraph<DBNode, DBEdge>, AHashMap<u64, NodeIndex>) {
    let mut seed_graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
    let mut node_lookup: AHashMap<u64, NodeIndex> = AHashMap::new();

    let k = kmer_counts.get_k();
    let suffix_mask = get_suffix_mask(&k);

    // Seed the graph with primer kmer nodes.
    //
    // Both forward and reverse primer kmers arrive from find_oligos_in_kmers()
    // oriented with the primer at the 5' (START) of the kmer. For forward
    // primers this is the sense strand; for reverse primers this is the
    // antisense strand.
    //
    // All graph nodes must be on the same strand for extension and path
    // finding to work (nodes are (k-1)-mers; two nodes representing the
    // same position on opposite strands have different sub_kmer values and
    // would never converge). We use the sense strand throughout.
    //
    // Forward seed sub_kmer derivation:
    //   Kmer is already sense-strand:  5' PPPPPPPPPPPPPPP|XXXXXX 3'
    //   sub_kmer = prefix (kmer >> 2): 5' PPPPPPPPPPPPPPPXXXXX  3'
    //   Forward extension appends bases at the 3' end (rightward).
    //
    // Reverse seed sub_kmer derivation (strand normalization):
    //   Kmer is antisense-strand:      5' RRRRRRRRRRRRRRR|XXXXXX 3'
    //   Revcomp to sense strand:       5' X'X'X'X'X'X'|R'R'R'R'R'R'R'R'R'R'R'R'R'R'R' 3'
    //   sub_kmer = suffix (& mask):    5' X'X'X'X'X'R'R'R'R'R'R'R'R'R'R'R'R'R'R'R'    3'
    //   This is the sense-strand (k-1)-mer at the 3' end of the amplicon,
    //   where forward extension arrives. Reverse extension prepends bases
    //   at the 5' end (leftward on sense strand).

    // Forward primer seeds
    let mut sorted_forward_kmers = forward_primer_kmers.kmers();
    sorted_forward_kmers.sort();
    for kmer in sorted_forward_kmers {
        let sub_kmer = kmer >> 2;

        if let Some(&existing) = node_lookup.get(&sub_kmer) {
            seed_graph[existing].is_start = true;
        } else {
            let node = seed_graph.add_node(DBNode {
                sub_kmer,
                is_start: true,
                is_end: false,
            });
            node_lookup.insert(sub_kmer, node);
        }
    }

    // Reverse primer seeds (strand-normalized to sense strand)
    let mut sorted_reverse_kmers = reverse_primer_kmers.kmers();
    sorted_reverse_kmers.sort();
    for kmer in sorted_reverse_kmers {
        let rc_kmer = crate::kmer::revcomp_kmer(&kmer, &k);
        let sub_kmer = rc_kmer & suffix_mask;

        if let Some(&existing) = node_lookup.get(&sub_kmer) {
            seed_graph[existing].is_end = true;
        } else {
            let node = seed_graph.add_node(DBNode {
                sub_kmer,
                is_start: false,
                is_end: true,
            });
            node_lookup.insert(sub_kmer, node);
        }
    }

    // Loop over the nodes and see if any of the nodes are start and end nodes
    // If so, print a warning
    for node in seed_graph.node_indices() {
        if seed_graph[node].is_start && seed_graph[node].is_end {
            debug!("node {} is both a start and end node", node.index());
        }
    }

    (seed_graph, node_lookup)
}

// Extend graph by adding edges and new nodes to non-terminal nodes.
// Terminal nodes have any of the following properties:
// - is_end = true
// - kmers that extend the node have been searched for but not found
// - The node has a path length longer than max-length - k + 1 from a start node

// Extension continues until all nodes have been visited

// Extension proceeds by:
// - For each non-terminal node, find the kmers that extend the node
// - Get the suffix of each kmer that extends the node
// - If a node with sub_kmer == suffix already exists, add an edge to the existing node
// - Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node

// where:
// - The prefix of the kmer of the node is the sub_kmer of the parent node in the graph
// - The suffix of the kmer of the node is the sub_kmer of the new node in the graph
// - If a node with the sub_kmer already exists, add a new edge to the existing node

/// Direction tag for entries in the unified bidirectional frontier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExtDir {
    Forward,
    Reverse,
}

/// Unified bidirectional graph extension. Processes forward and reverse seeds
/// in a single pass with an interleaved frontier. Each frontier entry is
/// tagged with its extension direction; new nodes inherit the direction of
/// their parent. Forward and reverse extension share the same global node
/// budget and operate on the same shared graph.
///
/// Returns `(graph, node_lookup, found_path)`. `found_path` is true when
/// forward and reverse extensions meet — i.e., a forward extension reaches a
/// node previously added by reverse extension (including `is_end` seed nodes),
/// or vice versa. When this happens, the graph contains a potential
/// start-to-end path.
///
/// Symmetric: swapping the "forward" and "reverse" labels on a primer pair
/// does not change behavior — both directions extend in lockstep.
#[allow(clippy::type_complexity)]
pub(super) fn extend_graph(
    mut graph: StableDiGraph<DBNode, DBEdge>,
    mut node_lookup: AHashMap<u64, NodeIndex>,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
    params: &PCRParams,
    max_num_nodes: usize,
) -> Result<(
    StableDiGraph<DBNode, DBEdge>,
    AHashMap<u64, NodeIndex>,
    bool,
)> {
    let suffix_mask: u64 = get_suffix_mask(&kmer_counts.get_k());
    let k = kmer_counts.get_k();
    let prefix_shift = 2 * (k - 1);
    let mut found_path = false;

    let mut last_check: usize = 0;
    let mut median_edge_count = compute_median_edge_count(&graph, *min_count as f64);
    let mut last_median_check: usize = 0;

    // Build the initial frontier from all seeds. Forward seeds get a Forward
    // entry; reverse seeds get a Reverse entry. Dual-role nodes (overlapping
    // primers, both is_start and is_end) get both entries so they can be
    // extended in both directions.
    let mut frontier: VecDeque<(NodeIndex, ExtDir)> = VecDeque::new();
    for node in graph.node_indices().collect::<Vec<_>>() {
        if graph[node].is_start {
            frontier.push_back((node, ExtDir::Forward));
        }
        if graph[node].is_end {
            frontier.push_back((node, ExtDir::Reverse));
        }
    }

    // Track per-direction processing. The DBNode visited flag is single-bit
    // and can't distinguish direction, so we use side tables here.
    let mut processed_fwd: AHashSet<NodeIndex> = AHashSet::new();
    let mut processed_rev: AHashSet<NodeIndex> = AHashSet::new();

    // Track which direction added each node. Seeds are tagged by their type
    // (is_start = Forward, is_end = Reverse, dual = both). New nodes inherit
    // their parent's direction. When forward extension reaches a node in the
    // Reverse set (or vice versa), the extensions have met in the middle and
    // a path is potentially complete.
    let mut added_by_fwd: AHashSet<NodeIndex> = AHashSet::new();
    let mut added_by_rev: AHashSet<NodeIndex> = AHashSet::new();
    for node in graph.node_indices() {
        if graph[node].is_start {
            added_by_fwd.insert(node);
        }
        if graph[node].is_end {
            added_by_rev.insert(node);
        }
    }

    while let Some((node, dir)) = frontier.pop_front() {
        // Skip if already processed in this direction
        let already_processed = match dir {
            ExtDir::Forward => !processed_fwd.insert(node),
            ExtDir::Reverse => !processed_rev.insert(node),
        };
        if already_processed {
            continue;
        }

        let n_nodes = graph.node_count();

        if n_nodes > max_num_nodes {
            gene_info!(
                params.gene_name,
                "There are {} nodes in the graph. This exceeds the maximum of {}, abandoning search.",
                n_nodes,
                max_num_nodes
            );
            break;
        }

        // Recompute cached median edge count periodically
        if n_nodes > last_median_check
            && (n_nodes - last_median_check) > EXTENSION_EVALUATION_FREQUENCY
        {
            median_edge_count = compute_median_edge_count(&graph, *min_count as f64);
            last_median_check = n_nodes - (n_nodes % EXTENSION_EVALUATION_FREQUENCY);
        }

        // Periodically log extension progress
        if (n_nodes > last_check) && ((n_nodes - last_check) > EXTENSION_EVALUATION_FREQUENCY) {
            last_check = n_nodes - (n_nodes % EXTENSION_EVALUATION_FREQUENCY);
            gene_info!(params.gene_name, "Evaluating bidirectional extension:");
            summarize_extension(&graph, "    ");
        }

        let sub_kmer = graph[node].sub_kmer;

        // Find candidate kmers based on direction
        let mut candidate_kmers: [Option<u64>; 4] = [None; 4];
        let mut n_candidates = 0;
        for base in 0..4u64 {
            let kmer = match dir {
                ExtDir::Forward => (sub_kmer << 2) | base,
                ExtDir::Reverse => (base << prefix_shift) | sub_kmer,
            };
            if let Some(count) = kmer_counts.get_canonical(&kmer) {
                if count >= *min_count {
                    candidate_kmers[n_candidates] = Some(kmer);
                    n_candidates += 1;
                }
            }
        }

        if n_candidates == 0 {
            // Dead end in this direction; the node may still extend in the
            // other direction (if it's a dual-role seed) so we don't mark
            // is_terminal here. The visited side-table prevents reprocessing
            // in this direction.
            continue;
        }

        for kmer in candidate_kmers.iter().flatten() {
            let new_sub_kmer = match dir {
                ExtDir::Forward => kmer & suffix_mask,
                ExtDir::Reverse => kmer >> 2,
            };

            // Self-loop: node would extend to itself
            if new_sub_kmer == sub_kmer {
                continue;
            }

            if let Some(&existing_node) = node_lookup.get(&new_sub_kmer) {
                // Connect to existing node. Edge direction depends on extension dir.
                let edge_check = match dir {
                    ExtDir::Forward => graph.find_edge(node, existing_node).is_none(),
                    ExtDir::Reverse => graph.find_edge(existing_node, node).is_none(),
                };
                if edge_check {
                    let edge = get_dbedge(kmer, kmer_counts);
                    match dir {
                        ExtDir::Forward => {
                            graph.add_edge(node, existing_node, edge);
                            // Path found if forward extension reaches any node
                            // that was added by reverse extension (or is_end seed)
                            if added_by_rev.contains(&existing_node) {
                                if !found_path {
                                    gene_info!(
                                        params.gene_name,
                                        "Forward and reverse extensions met (forward reached a reverse-added node)."
                                    );
                                }
                                found_path = true;
                            }
                        }
                        ExtDir::Reverse => {
                            graph.add_edge(existing_node, node, edge);
                            // Path found if reverse extension reaches any node
                            // that was added by forward extension (or is_start seed)
                            if added_by_fwd.contains(&existing_node) {
                                if !found_path {
                                    gene_info!(
                                        params.gene_name,
                                        "Forward and reverse extensions met (reverse reached a forward-added node)."
                                    );
                                }
                                found_path = true;
                            }
                        }
                    }
                }
            } else {
                let edge = get_dbedge(kmer, kmer_counts);
                let edge_count = edge.count;

                // Skip high-coverage edges (likely repetitive)
                if (edge_count as f64) > (median_edge_count * params.high_coverage_ratio) {
                    continue;
                }

                let new_node = graph.add_node(DBNode {
                    sub_kmer: new_sub_kmer,
                    is_start: false,
                    is_end: false,
                });
                node_lookup.insert(new_sub_kmer, new_node);

                // Tag the new node with the direction that added it
                match dir {
                    ExtDir::Forward => {
                        added_by_fwd.insert(new_node);
                        graph.add_edge(node, new_node, edge);
                    }
                    ExtDir::Reverse => {
                        added_by_rev.insert(new_node);
                        graph.add_edge(new_node, node, edge);
                    }
                }

                // Always add to frontier — max-length enforcement happens at
                // path finding time, not during extension. This avoids a
                // per-node BFS that was the dominant cost in bidirectional
                // extension.
                frontier.push_back((new_node, dir));
            }
        }
    }

    Ok((graph, node_lookup, found_path))
}

/// Annotate each edge with its coverage ratio: count / global median.
/// High ratios (>> 1.0) indicate potentially repetitive edges.
/// Low ratios (<< 1.0) indicate potential sequencing errors.
pub(super) fn annotate_coverage_ratios(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    let counts: Vec<u32> = graph.edge_indices().map(|e| graph[e].count).collect();
    let Some(median) = median_u32_f64(counts) else {
        return;
    };
    if median <= 0.0 {
        return;
    }

    for edge_idx in graph.edge_indices().collect::<Vec<_>>() {
        let count = graph[edge_idx].count as f64;
        graph[edge_idx].coverage_ratio = count / median;
    }
}

pub fn summarize_extension(graph: &StableDiGraph<DBNode, DBEdge>, pad: &str) {
    if !log::log_enabled!(log::Level::Debug) {
        return;
    }

    // Print the number of nodes and edges in the graph
    debug!("{}There are {} nodes in the graph", pad, graph.node_count());
    debug!("{}There are {} edges in the graph", pad, graph.edge_count());

    let has_cycles = is_cyclic_directed(graph);
    if has_cycles {
        debug!("{}The graph has cycles", pad);
    } else {
        debug!("{}The graph does not have cycles", pad);
    }

    // Print the mean, median, and max degree of all nodes
    let degrees: Vec<usize> = graph
        .node_indices()
        .map(|node| graph.neighbors(node).count())
        .collect();

    if degrees.is_empty() {
        debug!(
            "{}There are no nodes in the graph, terminating summary.",
            pad
        );
        return;
    }

    let max_degree = degrees
        .iter()
        .max()
        .expect("degrees is non-empty (checked above)");
    let degrees_u64: Vec<u64> = degrees.iter().map(|&x| x as u64).collect();
    let mean_degree = compute_mean(&degrees_u64);
    let median_degree = compute_median(&degrees_u64);

    debug!(
        "{}Max node degree {}, mean {:.2} median {:.1}",
        pad, max_degree, mean_degree, median_degree
    );

    // Print the mean, median, and max count of edges
    let counts: Vec<u64> = graph
        .edge_indices()
        .map(|e| graph[e].count as u64)
        .collect();

    if counts.is_empty() {
        debug!(
            "{}There are no edges in the graph, terminating summary.",
            pad
        );
        return;
    }

    let max_count = counts
        .iter()
        .max()
        .expect("counts is non-empty (checked above)");
    let mean_count = compute_mean(&counts);
    let median_count = compute_median(&counts);

    debug!(
        "{}Max edge count {}, mean {:.2} median {:.1}",
        pad, max_count, mean_count, median_count
    );

    // Create a vector of n_descendants for each node given a depth of EXTENSION_EVALUATION_DEPTH
    let n_descendants_vec: Vec<u64> = graph
        .node_indices()
        .map(|node| descendants(graph, node, EXTENSION_EVALUATION_DEPTH).len() as u64)
        .collect();

    let max_n_descendants = n_descendants_vec
        .iter()
        .max()
        .expect("n_descendants_vec is non-empty (checked above)");
    let mean_n_descendants = compute_mean(&n_descendants_vec);
    let median_n_descendants = compute_median(&n_descendants_vec);

    debug!(
        "{}Number of descendants to a depth of {}, max {} mean {:.2} median {:.1}",
        pad,
        EXTENSION_EVALUATION_DEPTH,
        max_n_descendants,
        mean_n_descendants,
        median_n_descendants
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_node(sub_kmer: u64, is_start: bool, is_end: bool) -> DBNode {
        DBNode {
            sub_kmer,
            is_start,
            is_end,
        }
    }

    fn mk_edge(count: u32) -> DBEdge {
        DBEdge {
            count,
            coverage_ratio: 1.0,
        }
    }

    #[test]
    fn test_compute_node_budget_low() {
        // Below low threshold: returns minimum budget
        assert_eq!(compute_node_budget(0), MIN_NODE_BUDGET);
        assert_eq!(compute_node_budget(BUDGET_LERP_LOW_BP), MIN_NODE_BUDGET);
    }

    #[test]
    fn test_compute_node_budget_high() {
        // Above high threshold: returns maximum budget
        assert_eq!(
            compute_node_budget(BUDGET_LERP_HIGH_BP),
            DEFAULT_MAX_NUM_NODES
        );
        assert_eq!(compute_node_budget(u64::MAX), DEFAULT_MAX_NUM_NODES);
    }

    #[test]
    fn test_compute_node_budget_interpolates() {
        let mid = (BUDGET_LERP_LOW_BP + BUDGET_LERP_HIGH_BP) / 2;
        let budget = compute_node_budget(mid);
        assert!(budget > MIN_NODE_BUDGET);
        assert!(budget < DEFAULT_MAX_NUM_NODES);
    }

    #[test]
    fn test_get_suffix_mask_values() {
        // k=3: mask for 2 nucleotides = 0b1111
        assert_eq!(get_suffix_mask(&3), 0b1111);
        // k=2: mask for 1 nucleotide = 0b11
        assert_eq!(get_suffix_mask(&2), 0b11);
    }

    #[test]
    fn test_compute_median_edge_count_empty() {
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        assert_eq!(compute_median_edge_count(&graph, 42.0), 42.0);
    }

    #[test]
    fn test_compute_median_edge_count_values() {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let a = graph.add_node(mk_node(0, true, false));
        let b = graph.add_node(mk_node(1, false, false));
        let c = graph.add_node(mk_node(2, false, true));
        graph.add_edge(a, b, mk_edge(5));
        graph.add_edge(b, c, mk_edge(15));
        graph.add_edge(a, c, mk_edge(10));
        // Sorted: [5, 10, 15], median = 10
        assert_eq!(compute_median_edge_count(&graph, 0.0), 10.0);
    }

    #[test]
    fn test_median_u32_f64_empty() {
        assert_eq!(median_u32_f64(vec![]), None);
    }

    #[test]
    fn test_median_u32_f64_single() {
        assert_eq!(median_u32_f64(vec![42]), Some(42.0));
    }

    #[test]
    fn test_median_u32_f64_odd() {
        // Sorted: [1, 5, 9] → median 5
        assert_eq!(median_u32_f64(vec![9, 1, 5]), Some(5.0));
    }

    #[test]
    fn test_median_u32_f64_even() {
        // Sorted: [1, 5, 9, 11] → median (5 + 9) / 2 = 7
        assert_eq!(median_u32_f64(vec![11, 1, 9, 5]), Some(7.0));
    }

    #[test]
    fn test_median_u32_f64_even_two() {
        // Smallest even case: [3, 7] → median 5
        assert_eq!(median_u32_f64(vec![7, 3]), Some(5.0));
    }

    #[test]
    fn test_median_u32_f64_large_values_no_overflow() {
        // u32::MAX-adjacent values whose pairwise sum would overflow u32
        // — verify we never wrap (the helper casts to f64 before adding).
        let large = u32::MAX - 1;
        let huge = u32::MAX;
        let m = median_u32_f64(vec![large, huge]).unwrap();
        // (large + huge) / 2 in f64 — within fp precision of (large+huge)/2.
        let expected = (large as f64 + huge as f64) / 2.0;
        assert!((m - expected).abs() < 1.0);
    }
}

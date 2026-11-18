// pcr/graph.rs — graph data structures and construction

use super::memory::buffers::{BudgetNodeSet, BudgetQueue, BudgetSet, BudgetVec};
use anyhow::Result;
use log::debug;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;

use crate::kmer::FilteredKmerCounts;
#[cfg(test)]
use crate::kmer::KmerCounts;
use crate::kmer::PrimerKmers;

use super::memory::AssemblyGraph;
use super::memory::lookup::NodeLookup;
use super::{DBEdge, DBNode, PCRParams};
use crate::kmer::budget::ResourceBudget;
use std::sync::Arc;

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

/// Generic median core: compute the median of `counts` as an f64, using
/// `select_nth_unstable`.
///
/// Expected-linear time, vs `sort()`'s O(n log n) — a measurable difference
/// because `compute_median_edge_count` runs every
/// `EXTENSION_EVALUATION_FREQUENCY` nodes during graph extension, where
/// `counts.len()` grows with the edge count of the working graph (up to
/// hundreds of thousands of edges on dense panels).
///
/// `to_f64` lets the caller convert each element to f64 with their own cast
/// (`x as f64`). This is necessary because `f64::From<u64>` does not exist
/// — `u64` is not losslessly convertible to f64 — so we can't just bound
/// `T: Into<f64>`. The cast is the same one the dedicated wrappers used to
/// inline; pulling it through a callback lets `median_u32_f64` and
/// `median_u64_f64` share the partition logic.
///
/// Returns `None` if `counts` is empty. The even-length case avoids the
/// extra full sort: after partitioning at `mid`, the lower half is
/// `counts[..mid]` (any order, but bounded above by `counts[mid]`), so a
/// single linear scan over that half finds the (mid-1)-th order statistic.
fn median_via_select<T>(counts: &mut [T], to_f64: impl Fn(T) -> f64) -> Option<f64>
where
    T: Copy + Ord,
{
    if counts.is_empty() {
        return None;
    }
    let mid = counts.len() / 2;
    if counts.len() % 2 == 0 {
        counts.select_nth_unstable(mid);
        let upper_min = to_f64(counts[mid]);
        let lower_max = to_f64(*counts[..mid].iter().max()?);
        Some((lower_max + upper_min) / 2.0)
    } else {
        Some(to_f64(*counts.select_nth_unstable(mid).1))
    }
}

#[cfg(test)]
pub(super) fn median_u32_f64(mut counts: Vec<u32>) -> Option<f64> {
    median_via_select(&mut counts, |x| x as f64)
}

#[cfg(test)]
pub(super) fn median_u64_f64(mut counts: Vec<u64>) -> Option<f64> {
    median_via_select(&mut counts, |x| x as f64)
}

pub(super) fn median_u64_slice(counts: &mut [u64]) -> Option<f64> {
    median_via_select(counts, |x| x as f64)
}

#[cfg(test)]
fn compute_median_edge_count(graph: &StableDiGraph<DBNode, DBEdge>, default: f64) -> f64 {
    median_edge_count_budgeted(graph, default, None).unwrap()
}

fn median_edge_count_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    default: f64,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<f64> {
    let mut counts = BudgetVec::collect(graph.edge_indices().map(|e| graph[e].count), budget)?;
    Ok(median_u64_slice(&mut counts).unwrap_or(default))
}

fn high_coverage_limit(
    median_edge_count: f64,
    high_coverage_ratio: f64,
    coverage_count_floor: u64,
) -> f64 {
    (median_edge_count * high_coverage_ratio).max(coverage_count_floor as f64)
}

pub fn get_dbedge(kmer: &u64, kmer_counts: &FilteredKmerCounts) -> Result<DBEdge> {
    Ok(DBEdge {
        count: kmer_counts.get_canonical_count(kmer)?,
        coverage_ratio: 0.0, // set during annotation pass
    })
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

#[cfg(test)]
pub fn get_start_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    graph
        .node_indices()
        .filter(|&node| graph[node].is_start)
        .collect()
}

#[cfg(test)]
pub fn get_end_nodes(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<NodeIndex> {
    graph
        .node_indices()
        .filter(|&node| graph[node].is_end)
        .collect()
}

// Find the descendants of a node, each descendent no more than `depth` edges away
#[cfg(test)]
pub(super) fn descendants(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    depth: usize,
) -> BudgetNodeSet {
    descendants_budgeted(graph, node, depth, None).unwrap()
}

fn descendants_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    depth: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetNodeSet> {
    let mut visited = BudgetNodeSet::new(budget);
    let mut queue = BudgetQueue::new(budget);
    queue.push_back((node, 0))?;
    visited.insert(node)?;
    let mut descendants = BudgetNodeSet::new(budget);
    while let Some((current_node, current_depth)) = queue.pop_front() {
        if current_depth >= depth {
            continue;
        }
        for neighbor in graph.neighbors(current_node) {
            if visited.insert(neighbor)? {
                descendants.insert(neighbor)?;
                queue.push_back((neighbor, current_depth + 1))?;
            }
        }
    }
    Ok(descendants)
}

pub fn compute_mean(numbers: &[u64]) -> f64 {
    if numbers.is_empty() {
        return 0.0;
    }
    let sum: u128 = numbers.iter().map(|&number| u128::from(number)).sum();
    sum as f64 / numbers.len() as f64
}

#[cfg(test)]
pub fn compute_median(numbers: &[u64]) -> f64 {
    // Empty input returns 0.0 (preserved from the previous behavior — callers
    // in path scoring and `summarize_extension` rely on a sentinel rather than
    // an Option).
    median_u64_f64(numbers.to_vec()).unwrap_or(0.0)
}

pub(super) fn create_seed_graph(
    forward_primer_kmers: &PrimerKmers,
    reverse_primer_kmers: &PrimerKmers,
    kmer_counts: &FilteredKmerCounts,
    graph_budget: Option<&Arc<ResourceBudget>>,
) -> Result<(AssemblyGraph, NodeLookup)> {
    let mut seed_graph = AssemblyGraph::new(graph_budget)?;
    let mut node_lookup = NodeLookup::new(graph_budget);

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
    let mut sorted_forward_kmers = BudgetVec::collect(
        forward_primer_kmers.iter().map(|(key, _)| *key),
        graph_budget,
    )?;
    sorted_forward_kmers.sort_unstable();
    for &kmer in sorted_forward_kmers.iter() {
        let sub_kmer = kmer >> 2;

        if let Some(&existing) = node_lookup.get(&sub_kmer) {
            seed_graph[existing].is_start = true;
        } else {
            let node = seed_graph.add_node(DBNode {
                sub_kmer,
                is_start: true,
                is_end: false,
            })?;
            node_lookup.insert(sub_kmer, node)?;
        }
    }

    drop(sorted_forward_kmers);

    // Reverse primer seeds (strand-normalized to sense strand)
    let mut sorted_reverse_kmers = BudgetVec::collect(
        reverse_primer_kmers.iter().map(|(key, _)| *key),
        graph_budget,
    )?;
    sorted_reverse_kmers.sort_unstable();
    for &kmer in sorted_reverse_kmers.iter() {
        let rc_kmer = crate::kmer::revcomp_kmer(&kmer, &k);
        let sub_kmer = rc_kmer & suffix_mask;

        if let Some(&existing) = node_lookup.get(&sub_kmer) {
            seed_graph[existing].is_end = true;
        } else {
            let node = seed_graph.add_node(DBNode {
                sub_kmer,
                is_start: false,
                is_end: true,
            })?;
            node_lookup.insert(sub_kmer, node)?;
        }
    }

    // Loop over the nodes and see if any of the nodes are start and end nodes
    // If so, print a warning
    for node in seed_graph.node_indices() {
        if seed_graph[node].is_start && seed_graph[node].is_end {
            debug!("node {} is both a start and end node", node.index());
        }
    }

    Ok((seed_graph, node_lookup))
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

#[derive(Debug, Default)]
pub(super) struct ExtensionRepeatMarkers {
    pub omitted_self_loop_sub_kmers: BudgetSet<u64>,
    pub retained_collision_edges: BudgetSet<EdgeIndex>,
}

#[derive(Debug, Default)]
pub(super) struct RepeatMarkers {
    pub omitted_self_loop_sub_kmers: BudgetSet<u64>,
    pub cyclic_scc_sub_kmers: BudgetSet<u64>,
    pub retained_collision_edges: BudgetSet<EdgeIndex>,
}

/// Unified bidirectional graph extension. Processes forward and reverse seeds
/// in a single pass with an interleaved frontier. Each frontier entry is
/// tagged with its extension direction; new nodes inherit the direction of
/// their parent. Forward and reverse extension share the same global node
/// budget and operate on the same shared graph.
///
/// Returns the graph, node lookup, connectivity status, repeat markers, and
/// node-budget status. Connectivity is true when
/// forward and reverse extensions meet — i.e., a forward extension reaches a
/// node previously added by reverse extension (including `is_end` seed nodes),
/// or vice versa. When this happens, the graph contains a potential
/// start-to-end path.
///
/// Symmetric: swapping the "forward" and "reverse" labels on a primer pair
/// does not change behavior — both directions extend in lockstep.
#[allow(clippy::type_complexity)]
pub(super) fn extend_graph(
    graph: impl Into<AssemblyGraph>,
    mut node_lookup: NodeLookup,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u64,
    params: &PCRParams,
    coverage_count_floor: u64,
    max_num_nodes: usize,
) -> Result<(
    AssemblyGraph,
    NodeLookup,
    bool,
    ExtensionRepeatMarkers,
    bool,
)> {
    let mut graph = graph.into();
    let suffix_mask: u64 = get_suffix_mask(&kmer_counts.get_k());
    let k = kmer_counts.get_k();
    let prefix_shift = 2 * (k - 1);
    let mut found_path = false;
    let mut node_budget_reached = false;
    let mut repeat_markers = ExtensionRepeatMarkers {
        omitted_self_loop_sub_kmers: BudgetSet::new(graph.budget()),
        retained_collision_edges: BudgetSet::new(graph.budget()),
    };

    let mut last_check: usize = 0;
    let mut median_edge_count =
        median_edge_count_budgeted(&graph, *min_count as f64, graph.budget())?;
    let mut last_median_check: usize = 0;

    // Build the initial frontier from all seeds. Forward seeds get a Forward
    // entry; reverse seeds get a Reverse entry. Dual-role nodes (overlapping
    // primers, both is_start and is_end) get both entries so they can be
    // extended in both directions.
    let mut frontier = BudgetQueue::new(graph.budget());
    for node in graph.node_indices() {
        if graph[node].is_start {
            frontier.push_back((node, ExtDir::Forward))?;
        }
        if graph[node].is_end {
            frontier.push_back((node, ExtDir::Reverse))?;
        }
    }

    // Track per-direction processing. The DBNode visited flag is single-bit
    // and can't distinguish direction, so we use side tables here.
    let mut processed_fwd = BudgetNodeSet::new(graph.budget());
    let mut processed_rev = BudgetNodeSet::new(graph.budget());

    // Track which direction added each node. Seeds are tagged by their type
    // (is_start = Forward, is_end = Reverse, dual = both). New nodes inherit
    // their parent's direction. When forward extension reaches a node in the
    // Reverse set (or vice versa), the extensions have met in the middle and
    // a path is potentially complete.
    let mut added_by_fwd = BudgetNodeSet::new(graph.budget());
    let mut added_by_rev = BudgetNodeSet::new(graph.budget());
    for node in graph.node_indices() {
        if graph[node].is_start {
            added_by_fwd.insert(node)?;
        }
        if graph[node].is_end {
            added_by_rev.insert(node)?;
        }
    }

    while let Some((node, dir)) = frontier.pop_front() {
        // Skip if already processed in this direction
        let already_processed = match dir {
            ExtDir::Forward => !processed_fwd.insert(node)?,
            ExtDir::Reverse => !processed_rev.insert(node)?,
        };
        if already_processed {
            continue;
        }

        let n_nodes = graph.node_count();

        if n_nodes > max_num_nodes {
            node_budget_reached = true;
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
            median_edge_count =
                median_edge_count_budgeted(&graph, *min_count as f64, graph.budget())?;
            last_median_check = n_nodes - (n_nodes % EXTENSION_EVALUATION_FREQUENCY);
        }

        // Periodically log extension progress
        if (n_nodes > last_check) && ((n_nodes - last_check) > EXTENSION_EVALUATION_FREQUENCY) {
            last_check = n_nodes - (n_nodes % EXTENSION_EVALUATION_FREQUENCY);
            gene_info!(params.gene_name, "Evaluating bidirectional extension:");
            summarize_extension(&graph, "    ", graph.budget())?;
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
            if let Some(count) = kmer_counts.get_canonical(&kmer)?
                && count >= *min_count
            {
                candidate_kmers[n_candidates] = Some(kmer);
                n_candidates += 1;
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
                repeat_markers
                    .omitted_self_loop_sub_kmers
                    .insert(sub_kmer)?;
                continue;
            }

            let edge = get_dbedge(kmer, kmer_counts)?;
            let exceeds_coverage_limit = (edge.count as f64)
                > high_coverage_limit(
                    median_edge_count,
                    params.high_coverage_ratio,
                    coverage_count_floor,
                );

            if let Some(&existing_node) = node_lookup.get(&new_sub_kmer) {
                // Connect to existing node. Edge direction depends on extension dir.
                let existing_edge = match dir {
                    ExtDir::Forward => graph.find_edge(node, existing_node),
                    ExtDir::Reverse => graph.find_edge(existing_node, node),
                };
                let collision_edge = if let Some(edge_id) = existing_edge {
                    edge_id
                } else {
                    match dir {
                        ExtDir::Forward => {
                            let edge_id = graph.add_edge(node, existing_node, edge)?;
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
                            edge_id
                        }
                        ExtDir::Reverse => {
                            let edge_id = graph.add_edge(existing_node, node, edge)?;
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
                            edge_id
                        }
                    }
                };
                if exceeds_coverage_limit {
                    repeat_markers
                        .retained_collision_edges
                        .insert(collision_edge)?;
                }
            } else {
                // Skip high-coverage edges (likely repetitive)
                if exceeds_coverage_limit {
                    continue;
                }

                let new_node = graph.add_node(DBNode {
                    sub_kmer: new_sub_kmer,
                    is_start: false,
                    is_end: false,
                })?;
                node_lookup.insert(new_sub_kmer, new_node)?;

                // Tag the new node with the direction that added it
                match dir {
                    ExtDir::Forward => {
                        added_by_fwd.insert(new_node)?;
                        graph.add_edge(node, new_node, edge)?;
                    }
                    ExtDir::Reverse => {
                        added_by_rev.insert(new_node)?;
                        graph.add_edge(new_node, node, edge)?;
                    }
                }

                // Always add to frontier — max-length enforcement happens at
                // path finding time, not during extension. This avoids a
                // per-node BFS that was the dominant cost in bidirectional
                // extension.
                frontier.push_back((new_node, dir))?;
            }
        }
    }

    Ok((
        graph,
        node_lookup,
        found_path,
        repeat_markers,
        node_budget_reached,
    ))
}

#[cfg(test)]
pub(super) fn repeat_markers_before_pruning(
    graph: &StableDiGraph<DBNode, DBEdge>,
    extension_markers: ExtensionRepeatMarkers,
) -> RepeatMarkers {
    repeat_markers_before_pruning_budgeted(graph, extension_markers, None).unwrap()
}

pub(super) fn repeat_markers_before_pruning_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    extension_markers: ExtensionRepeatMarkers,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<RepeatMarkers> {
    debug_assert!(
        extension_markers
            .retained_collision_edges
            .iter()
            .all(|edge| graph.edge_weight(*edge).is_some())
    );
    Ok(RepeatMarkers {
        omitted_self_loop_sub_kmers: extension_markers.omitted_self_loop_sub_kmers,
        cyclic_scc_sub_kmers: cyclic_scc_sub_kmers_budgeted(graph, budget)?,
        retained_collision_edges: extension_markers.retained_collision_edges,
    })
}

/// Iterative Kosaraju scan with owned buffers; no uncharged component arrays.
fn cyclic_scc_sub_kmers_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetSet<u64>> {
    let mut visited = BudgetNodeSet::new(budget);
    let mut stack = BudgetVec::new(budget);
    let mut order = BudgetVec::new(budget);
    for root in graph.node_indices() {
        if visited.contains(&root) {
            continue;
        }
        stack.push((root, false))?;
        while let Some((node, exiting)) = stack.pop() {
            if exiting {
                order.push(node)?;
            } else if visited.insert(node)? {
                stack.push((node, true))?;
                for neighbor in graph.neighbors_directed(node, petgraph::Direction::Outgoing) {
                    if !visited.contains(&neighbor) {
                        stack.push((neighbor, false))?;
                    }
                }
            }
        }
    }
    visited.clear();
    let mut component = BudgetVec::new(budget);
    let mut cyclic = BudgetSet::new(budget);
    while let Some(root) = order.pop() {
        if visited.contains(&root) {
            continue;
        }
        component.clear();
        stack.push((root, false))?;
        while let Some((node, _)) = stack.pop() {
            if visited.insert(node)? {
                component.push(node)?;
                for neighbor in graph.neighbors_directed(node, petgraph::Direction::Incoming) {
                    if !visited.contains(&neighbor) {
                        stack.push((neighbor, false))?;
                    }
                }
            }
        }
        if component.len() > 1 || graph.find_edge(root, root).is_some() {
            for &node in component.iter() {
                cyclic.insert(graph[node].sub_kmer)?;
            }
        }
    }
    Ok(cyclic)
}

/// Annotate each edge with its coverage ratio: count / global median.
/// High ratios (>> 1.0) indicate potentially repetitive edges.
/// Low ratios (<< 1.0) indicate potential sequencing errors.
pub(super) fn annotate_coverage_ratios(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<()> {
    let median = median_edge_count_budgeted(graph, 0.0, budget)?;
    if median <= 0.0 {
        return Ok(());
    }
    // Weight iteration changes scalars only and needs no copied edge IDs.
    for edge in graph.edge_weights_mut() {
        edge.coverage_ratio = edge.count as f64 / median;
    }
    Ok(())
}

pub fn summarize_extension(
    graph: &StableDiGraph<DBNode, DBEdge>,
    pad: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<()> {
    if !log::log_enabled!(log::Level::Debug) {
        return Ok(());
    }
    debug!("{}There are {} nodes in the graph", pad, graph.node_count());
    debug!("{}There are {} edges in the graph", pad, graph.edge_count());
    let mut degrees = BudgetVec::collect(
        graph
            .node_indices()
            .map(|node| graph.neighbors(node).count() as u64),
        budget,
    )?;
    let Some(max_degree) = degrees.iter().max().copied() else {
        debug!(
            "{}There are no nodes in the graph, terminating summary.",
            pad
        );
        return Ok(());
    };
    let mean_degree = compute_mean(&degrees);
    let median_degree = median_u64_slice(&mut degrees).unwrap_or(0.0);
    debug!(
        "{}Max node degree {}, mean {:.2} median {:.1}",
        pad, max_degree, mean_degree, median_degree
    );
    drop(degrees);
    let mut counts =
        BudgetVec::collect(graph.edge_indices().map(|edge| graph[edge].count), budget)?;
    let Some(max_count) = counts.iter().max().copied() else {
        debug!(
            "{}There are no edges in the graph, terminating summary.",
            pad
        );
        return Ok(());
    };
    let mean_count = compute_mean(&counts);
    let median_count = median_u64_slice(&mut counts).unwrap_or(0.0);
    debug!(
        "{}Max edge count {}, mean {:.2} median {:.1}",
        pad, max_count, mean_count, median_count
    );
    drop(counts);
    let mut descendants = BudgetVec::new(budget);
    for node in graph.node_indices() {
        descendants.push(
            descendants_budgeted(graph, node, EXTENSION_EVALUATION_DEPTH, budget)?.len() as u64,
        )?;
    }
    let max_descendants = descendants.iter().max().copied().unwrap_or(0);
    let mean_descendants = compute_mean(&descendants);
    let median_descendants = median_u64_slice(&mut descendants).unwrap_or(0.0);
    debug!(
        "{}Number of descendants to a depth of {}, max {} mean {:.2} median {:.1}",
        pad, EXTENSION_EVALUATION_DEPTH, max_descendants, mean_descendants, median_descendants
    );
    Ok(())
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

    fn mk_edge(count: u64) -> DBEdge {
        DBEdge {
            count,
            coverage_ratio: 1.0,
        }
    }

    fn workspace_budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; crate::kmer::budget::MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    #[test]
    fn owned_scc_scan_matches_petgraph_on_directed_graphs_with_holes_and_parallel_edges() {
        let budget = workspace_budget(1 << 20);
        let mut random = 157u64;
        for round in 0..96 {
            let mut graph = StableDiGraph::new();
            for index in 0..33 {
                graph.add_node(mk_node(round * 100 + index, index == 0, index == 32));
            }
            graph.remove_node(NodeIndex::new(5));
            graph.remove_node(NodeIndex::new(17));
            let nodes: Vec<_> = graph.node_indices().collect();
            for _ in 0..80 {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let from = nodes[(random >> 32) as usize % nodes.len()];
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let to = nodes[(random >> 32) as usize % nodes.len()];
                graph.add_edge(from, to, mk_edge(10));
            }
            let expected: std::collections::BTreeSet<_> = petgraph::algo::kosaraju_scc(&graph)
                .into_iter()
                .filter(|group| group.len() > 1 || graph.find_edge(group[0], group[0]).is_some())
                .flatten()
                .map(|node| graph[node].sub_kmer)
                .collect();
            let result = cyclic_scc_sub_kmers_budgeted(&graph, Some(&budget)).unwrap();
            assert_eq!(
                result
                    .iter()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                expected
            );
            drop(result);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn scc_and_median_workspace_exhaustion_preserve_graph_and_release_all_buffers() {
        let mut graph = StableDiGraph::new();
        let a = graph.add_node(mk_node(1, true, false));
        let b = graph.add_node(mk_node(2, false, true));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, a, mk_edge(20));
        let tiny = workspace_budget(1);
        let error = cyclic_scc_sub_kmers_budgeted(&graph, Some(&tiny)).unwrap_err();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        assert!(
            median_edge_count_budgeted(&graph, 0.0, Some(&tiny))
                .unwrap_err()
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        assert_eq!((graph.node_count(), graph.edge_count()), (2, 2));
        let enough = workspace_budget(4096);
        let result = cyclic_scc_sub_kmers_budgeted(&graph, Some(&enough)).unwrap();
        assert_eq!(result.len(), 2);
        drop(result);
        assert_eq!(
            median_edge_count_budgeted(&graph, 0.0, Some(&enough)).unwrap(),
            15.0
        );
        assert_eq!(enough.snapshot().reserved_memory_bytes, 0);
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
    fn test_high_coverage_limit_uses_absolute_floor() {
        assert_eq!(high_coverage_limit(2.0, 10.0, 104), 104.0);
        assert_eq!(high_coverage_limit(100.0, 10.0, 104), 1000.0);
        assert_eq!(high_coverage_limit(100.0, 2.0, 104), 200.0);
    }

    #[test]
    fn retained_high_coverage_collision_marks_new_edge_not_endpoints() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(
            crate::kmer::seq_to_kmer("AA").unwrap(),
            true,
            false,
        ));
        let end = graph.add_node(mk_node(
            crate::kmer::seq_to_kmer("AC").unwrap(),
            false,
            true,
        ));
        let lookup =
            NodeLookup::from_iter([(graph[start].sub_kmer, start), (graph[end].sub_kmer, end)]);
        let kmer_length = 3;
        let mut counts = KmerCounts::new(&kmer_length);
        for _copy in 0..100 {
            counts.ingest_seq("AAC").unwrap();
        }
        let params = crate::cli::parse_pcr_primers_string(
            "name=collision,forward=AA,reverse=GT,min-length=0,max-length=100,mismatches=0",
        )
        .unwrap();

        let (extended, _, found, markers, _) =
            extend_graph(graph, lookup, &counts.filtered_view(1), &1, &params, 0, 100).unwrap();
        let collision_edge = extended.find_edge(start, end).unwrap();

        assert!(found);
        assert!(markers.omitted_self_loop_sub_kmers.is_empty());
        assert_eq!(markers.retained_collision_edges.len(), 1);
        assert!(markers.retained_collision_edges.contains(&collision_edge));
    }

    #[test]
    fn retained_high_coverage_collision_marks_already_present_edge() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(
            crate::kmer::seq_to_kmer("AA").unwrap(),
            true,
            false,
        ));
        let end = graph.add_node(mk_node(
            crate::kmer::seq_to_kmer("AC").unwrap(),
            false,
            true,
        ));
        let existing_edge = graph.add_edge(start, end, mk_edge(1));
        let lookup =
            NodeLookup::from_iter([(graph[start].sub_kmer, start), (graph[end].sub_kmer, end)]);
        let kmer_length = 3;
        let mut counts = KmerCounts::new(&kmer_length);
        for _copy in 0..100 {
            counts.ingest_seq("AAC").unwrap();
        }
        let params = crate::cli::parse_pcr_primers_string(
            "name=collision,forward=AA,reverse=GT,min-length=0,max-length=100,mismatches=0",
        )
        .unwrap();

        let (_, _, _, markers, _) =
            extend_graph(graph, lookup, &counts.filtered_view(1), &1, &params, 0, 100).unwrap();

        assert!(markers.omitted_self_loop_sub_kmers.is_empty());
        assert_eq!(markers.retained_collision_edges.len(), 1);
        assert!(markers.retained_collision_edges.contains(&existing_edge));
    }

    #[test]
    fn test_repeat_scan_handles_deep_graph_on_small_stack() {
        let mut graph = StableDiGraph::new();
        let nodes = (0..20_000)
            .map(|sub_kmer| graph.add_node(mk_node(sub_kmer, false, false)))
            .collect::<Vec<_>>();
        for pair in nodes.windows(2) {
            graph.add_edge(pair[0], pair[1], mk_edge(1));
        }

        let unresolved_count = std::thread::Builder::new()
            .stack_size(64 * 1024)
            .spawn(move || cyclic_scc_sub_kmers_budgeted(&graph, None).unwrap().len())
            .unwrap()
            .join()
            .unwrap();

        assert_eq!(unresolved_count, 0);
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

    #[test]
    fn test_median_u64_f64_smoke() {
        assert_eq!(median_u64_f64(vec![]), None);
        assert_eq!(median_u64_f64(vec![42]), Some(42.0));
        // Odd: sorted [1, 5, 9] → 5
        assert_eq!(median_u64_f64(vec![9, 1, 5]), Some(5.0));
        // Even: sorted [1, 5, 9, 11] → (5+9)/2 = 7
        assert_eq!(median_u64_f64(vec![11, 1, 9, 5]), Some(7.0));
    }

    #[test]
    fn test_compute_median_matches_helper() {
        // compute_median is the legacy entry point used from path scoring
        // and summarize_extension; verify it still produces the same value
        // as the helper after the rewire.
        let v: Vec<u64> = vec![3, 1, 4, 1, 5, 9, 2, 6];
        // Sorted: [1, 1, 2, 3, 4, 5, 6, 9] → median (3+4)/2 = 3.5
        assert_eq!(compute_median(&v), 3.5);
        // Empty stays 0.0 (sentinel).
        assert_eq!(compute_median(&[]), 0.0);
    }
}

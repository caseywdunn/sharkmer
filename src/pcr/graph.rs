// pcr/graph.rs — graph data structures and construction

use anyhow::{Context, Result};
use log::{debug, trace};
use petgraph::Direction;
use petgraph::algo::is_cyclic_directed;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use crate::kmer::FilteredKmerCounts;
use crate::kmer::KmerCounts;

use super::{DBEdge, DBNode, PCRParams};

/// How frequently to check extension graph for ballooning growth
const EXTENSION_EVALUATION_FREQUENCY: usize = 1_000;

/// The graph depth over which to evaluate ballooning growth
pub(super) const EXTENSION_EVALUATION_DEPTH: usize = 4;

/// Default node budget — give up if the graph gets too large
pub const DEFAULT_MAX_NUM_NODES: usize = 200_000;

// HIGH_COVERAGE_RATIO_THRESHOLD is now read from params.high_coverage_ratio

// A mask that can be used to isolate the last k-1 nucleotides of a kmer
pub(super) fn get_suffix_mask(k: &usize) -> u64 {
    (1 << (2 * (*k - 1))) - 1
}

pub fn n_nonterminal_nodes_in_graph(graph: &StableDiGraph<DBNode, DBEdge>) -> usize {
    graph
        .node_indices()
        .filter(|&node| !graph[node].is_terminal)
        .count()
}

pub fn n_unvisited_nodes_in_graph(graph: &StableDiGraph<DBNode, DBEdge>) -> usize {
    graph
        .node_indices()
        .filter(|&node| !graph[node].visited)
        .count()
}

fn compute_median_edge_count(graph: &StableDiGraph<DBNode, DBEdge>, default: f64) -> f64 {
    let mut counts: Vec<u32> = graph.edge_indices().map(|e| graph[e].count).collect();
    if counts.is_empty() {
        default
    } else {
        counts.sort();
        let mid = counts.len() / 2;
        if counts.len() % 2 == 0 {
            (counts[mid - 1] as f64 + counts[mid] as f64) / 2.0
        } else {
            counts[mid] as f64
        }
    }
}

pub fn get_path_length(
    graph: &StableDiGraph<DBNode, DBEdge>,
    new_node: NodeIndex,
) -> Result<Option<usize>> {
    // Get the length of the path from the start node to the new node.
    // Use a bounded counter instead of a HashSet to detect cycles:
    // the graph has at most graph.node_count() nodes, so any path
    // longer than that must contain a cycle.
    let max_steps = graph.node_count();
    let mut path_length = 0;
    let mut current_node = new_node;

    loop {
        // If we've taken more steps than there are nodes, we're in a cycle
        if path_length > max_steps {
            return Ok(None);
        }

        let current_node_data = graph
            .node_weight(current_node)
            .context("Node not found in graph during path length calculation")?;
        if current_node_data.is_start {
            break;
        }

        path_length += 1;
        if let Some(next_node) = graph
            .neighbors_directed(current_node, Direction::Incoming)
            .next()
        {
            current_node = next_node;
        } else {
            // Node has no incoming edge — may be a reverse-extended node
            // that hasn't connected to a start node yet
            return Ok(None);
        }
    }
    Ok(Some(path_length))
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
) -> HashSet<NodeIndex> {
    let mut visited: HashSet<NodeIndex> = HashSet::new();
    let mut queue = VecDeque::new();

    queue.push_back((node, 0));
    visited.insert(node);

    let mut descendants: HashSet<NodeIndex> = HashSet::new();

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

// Get a vector of edge counts by traversing the graph backwards from the focal node
#[allow(dead_code)]
pub(super) fn get_backward_edge_counts(
    graph: &StableDiGraph<DBNode, DBEdge>,
    focal_node: NodeIndex,
    depth: usize,
) -> Vec<u32> {
    let mut edge_counts = Vec::new();
    let mut current_node = focal_node;
    let mut current_depth = 0;

    while current_depth < depth {
        if let Some(edge) = graph
            .edges_directed(current_node, Direction::Incoming)
            .next()
        {
            edge_counts.push(edge.weight().count);
            current_node = edge.source();
        } else {
            break;
        }

        current_depth += 1;
    }

    edge_counts
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
) -> (StableDiGraph<DBNode, DBEdge>, HashMap<u64, NodeIndex>) {
    let mut seed_graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
    let mut node_lookup: HashMap<u64, NodeIndex> = HashMap::new();

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
                is_terminal: false,
                visited: false,
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
                is_terminal: false,
                visited: false,
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
/// Reset terminal and visited flags on nodes that may have new successors or
/// predecessors at a lower count threshold. Nodes that are both start AND end
/// (overlapping primers) stay terminal. Nodes at max-length from their
/// respective direction's seed stay terminal. All other terminal nodes are
/// un-terminaled and un-visited so extension can retry them.
pub(super) fn prepare_for_lower_threshold(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
) {
    // Collect nodes to reset (cannot mutate graph while iterating)
    let nodes_to_reset: Vec<NodeIndex> = graph
        .node_indices()
        .filter(|&node| {
            let data = &graph[node];
            if !data.is_terminal {
                return false;
            }
            // Nodes that are both start and end (overlapping primers) stay terminal
            if data.is_start && data.is_end {
                return false;
            }
            // Keep terminal if at max-length from start (forward direction)
            if !data.is_end {
                if let Ok(Some(path_length)) = get_path_length(graph, node) {
                    if path_length + kmer_counts.get_k() > params.max_length {
                        return false;
                    }
                }
            }
            // Keep terminal if at max-length from end (reverse direction)
            if !data.is_start {
                if let Ok(Some(path_length)) = get_path_length_from_end(graph, node) {
                    if path_length + kmer_counts.get_k() > params.max_length {
                        return false;
                    }
                }
            }
            true
        })
        .collect();

    for node in nodes_to_reset {
        graph[node].is_terminal = false;
        graph[node].visited = false;
    }
}

/// Mask all seed nodes NOT in the given component set by marking them
/// visited and terminal. Returns saved states for later restoration.
pub(super) fn mask_other_components(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    component_seeds: &std::collections::HashSet<NodeIndex>,
) -> Vec<(NodeIndex, bool, bool)> {
    let mut saved = Vec::new();
    for node in graph.node_indices().collect::<Vec<_>>() {
        let data = &graph[node];
        if (data.is_start || data.is_end) && !component_seeds.contains(&node) {
            saved.push((node, data.visited, data.is_terminal));
            graph[node].visited = true;
            graph[node].is_terminal = true;
        }
    }
    saved
}

/// Restore seed node states that were masked by `mask_other_components`.
pub(super) fn unmask_nodes(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    saved: &[(NodeIndex, bool, bool)],
) {
    for &(node, visited, is_terminal) in saved {
        if graph.contains_node(node) {
            graph[node].visited = visited;
            graph[node].is_terminal = is_terminal;
        }
    }
}

#[allow(clippy::type_complexity)]
pub(super) fn extend_graph(
    mut graph: StableDiGraph<DBNode, DBEdge>,
    mut node_lookup: HashMap<u64, NodeIndex>,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
    params: &PCRParams,
    max_num_nodes: usize,
    component_budget: Option<usize>,
) -> Result<(StableDiGraph<DBNode, DBEdge>, HashMap<u64, NodeIndex>, bool)> {
    let suffix_mask: u64 = get_suffix_mask(&kmer_counts.get_k());
    let mut found_end_node = false;
    let nodes_at_start = graph.node_count();

    // Mark end-only nodes as visited so forward extension skips them.
    // They will be un-visited by reverse extension.
    for node in graph.node_indices().collect::<Vec<_>>() {
        if graph[node].is_end && !graph[node].is_start {
            graph[node].visited = true;
        }
    }

    let mut last_check: usize = 0;
    let mut median_edge_count = compute_median_edge_count(&graph, *min_count as f64);
    let mut last_median_check: usize = 0;

    while n_unvisited_nodes_in_graph(&graph) > 0 {
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

        if let Some(budget) = component_budget {
            if n_nodes.saturating_sub(nodes_at_start) >= budget {
                gene_info!(
                    params.gene_name,
                    "Component budget exhausted ({} nodes added, budget {}), moving to next component.",
                    n_nodes - nodes_at_start,
                    budget
                );
                break;
            }
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

            gene_info!(params.gene_name, "Evaluating extension:");
            summarize_extension(&graph, "    ");
        }

        // Iterate over the nodes
        let node_indices: Vec<_> = graph.node_indices().collect();
        for node in node_indices {
            if !(graph[node].visited) {
                // Get the suffix of the kmer of the node
                let sub_kmer = graph[node].sub_kmer;

                trace!(
                    "  {} sub_kmer being extended for node {}.",
                    crate::kmer::kmer_to_seq(&sub_kmer, &(kmer_counts.get_k() - 1)),
                    node.index()
                );

                // Get the kmers that could extend the node (at most 4, one per base)
                let mut candidate_kmers: [Option<u64>; 4] = [None; 4];
                let mut n_candidates = 0;

                for base in 0..4u64 {
                    let kmer = (sub_kmer << 2) | base;

                    // If the kmer is in the kmer_counts hash (either orientation)
                    // and has count >= min_count, add it to the candidate kmers
                    if let Some(count) = kmer_counts.get_canonical(&kmer)
                        && count >= *min_count
                    {
                        candidate_kmers[n_candidates] = Some(kmer);
                        n_candidates += 1;
                    }
                }

                trace!("There are {} candidate kmers for extension.", n_candidates);

                // If there are no candidate kmers, the node is terminal
                if n_candidates == 0 {
                    trace!(
                        "Marking node as terminal because there are no candidates for extension."
                    );
                    graph[node].is_terminal = true;
                    graph[node].visited = true;
                    continue;
                }

                // Add new nodes if needed, and new edges
                for kmer in candidate_kmers.iter().flatten() {
                    let suffix = kmer & suffix_mask;

                    // Check if the node extends by itself and mark it as terminal if it does
                    if suffix == sub_kmer {
                        graph[node].is_terminal = true;
                        graph[node].visited = true;
                        trace!("Node {} extends itself. Marking as terminal.", node.index());
                        break;
                    }

                    // If the node with sub_kmer == suffix already exists, add an edge to the existing node
                    // Otherwise, create a new node with sub_kmer == suffix, and add an edge to the new node

                    if let Some(&existing_node) = node_lookup.get(&suffix) {
                        // Allow edges that would form cycles — cycles are
                        // handled during path finding via bounded revisitation.
                        // Only add the edge if it doesn't already exist.
                        if graph.find_edge(node, existing_node).is_none() {
                            let edge = get_dbedge(kmer, kmer_counts);
                            graph.add_edge(node, existing_node, edge);
                            if graph[existing_node].is_end {
                                found_end_node = true;
                                gene_info!(
                                    params.gene_name,
                                    "End node incorporated into graph, complete PCR product found."
                                );
                            }
                        }
                    } else {
                        let edge = get_dbedge(kmer, kmer_counts);
                        let edge_count = edge.count;

                        // Skip edges with count far above the median — they
                        // likely lead into repetitive regions
                        if (edge_count as f64) > (median_edge_count * params.high_coverage_ratio) {
                            trace!(
                                "Skipping high-coverage edge (count {} > {:.0} × {:.0} median). kmer {}",
                                edge_count,
                                params.high_coverage_ratio,
                                median_edge_count,
                                crate::kmer::kmer_to_seq(kmer, &kmer_counts.get_k()),
                            );
                            continue;
                        }

                        let new_node = graph.add_node(DBNode {
                            sub_kmer: suffix,
                            is_start: false,
                            is_end: false,
                            is_terminal: false,
                            visited: false,
                        });
                        node_lookup.insert(suffix, new_node);

                        graph.add_edge(node, new_node, edge);

                        trace!(
                            "Added sub_kmer {} for new node {} with edge kmer count {}.",
                            crate::kmer::kmer_to_seq(&suffix, &(kmer_counts.get_k() - 1)),
                            new_node.index(),
                            edge_count
                        );

                        // Check if the new node is max-length - k + 1 from a start node
                        // If so, mark the new node as terminal
                        let path_length = get_path_length(&graph, new_node)?;

                        // If the path length is None, the node is part of a cycle and is marked terminal.
                        // If the path length is Some, is marked terminal if the path length is >= max-length - k + 1
                        if let Some(path_length) = path_length {
                            trace!("Path length is {}.", path_length);

                            if path_length + kmer_counts.get_k() > params.max_length {
                                graph[new_node].is_terminal = true;
                                graph[new_node].visited = true;
                                trace!(
                                    "Marking new node {} as terminal because it exceeds max-length from start.",
                                    new_node.index()
                                );
                            }
                        } else {
                            graph[new_node].is_terminal = true;
                            trace!(
                                "Marking new node {} as terminal because it is part of a cycle.",
                                new_node.index()
                            );
                        }
                    }
                }
                graph[node].visited = true;

                trace!(
                    "There are now {} unvisited and {} non-terminal nodes in the graph.",
                    n_unvisited_nodes_in_graph(&graph),
                    n_nonterminal_nodes_in_graph(&graph)
                );
            }
        }
    }

    Ok((graph, node_lookup, found_end_node))
}

/// Trace forward (outgoing edges) from a node to find the shortest path to an
/// end node. Returns `None` if a cycle is encountered before reaching an end.
pub fn get_path_length_from_end(
    graph: &StableDiGraph<DBNode, DBEdge>,
    new_node: NodeIndex,
) -> Result<Option<usize>> {
    let max_steps = graph.node_count();
    let mut path_length = 0;
    let mut current_node = new_node;

    loop {
        if path_length > max_steps {
            return Ok(None);
        }

        let current_node_data = graph
            .node_weight(current_node)
            .context("Node not found in graph during path length calculation")?;
        if current_node_data.is_end {
            break;
        }

        path_length += 1;
        if let Some(next_node) = graph
            .neighbors_directed(current_node, Direction::Outgoing)
            .next()
        {
            current_node = next_node;
        } else {
            // Node has no outgoing edge and is not an end node — during
            // reverse extension, new nodes may not yet have outgoing paths.
            return Ok(None);
        }
    }
    Ok(Some(path_length))
}

/// Extend graph in reverse from end nodes by adding predecessor edges/nodes.
/// Mirrors `extend_graph` but extends leftward: for each unvisited end node,
/// finds candidate predecessor kmers and adds them to the graph.
pub(super) fn extend_graph_reverse(
    mut graph: StableDiGraph<DBNode, DBEdge>,
    mut node_lookup: HashMap<u64, NodeIndex>,
    kmer_counts: &FilteredKmerCounts,
    min_count: &u32,
    params: &PCRParams,
    max_num_nodes: usize,
    component_budget: Option<usize>,
) -> Result<(StableDiGraph<DBNode, DBEdge>, HashMap<u64, NodeIndex>)> {
    let k = kmer_counts.get_k();
    let prefix_shift = 2 * (k - 1);
    let nodes_at_start = graph.node_count();

    let mut last_check: usize = 0;
    let mut median_edge_count = compute_median_edge_count(&graph, *min_count as f64);
    let mut last_median_check: usize = 0;

    // Un-visit end-only nodes so reverse extension processes them.
    // Forward extension marked them visited to skip them.
    for node in graph.node_indices().collect::<Vec<_>>() {
        if graph[node].is_end && !graph[node].is_start && !graph[node].is_terminal {
            graph[node].visited = false;
        }
    }
    loop {
        // Find unvisited nodes that should be reverse-extended.
        // These are nodes that are either end nodes or were created by reverse extension.
        let unvisited: Vec<NodeIndex> = graph
            .node_indices()
            .filter(|&node| !graph[node].visited)
            .collect();

        if unvisited.is_empty() {
            break;
        }

        let n_nodes = graph.node_count();

        if n_nodes > max_num_nodes {
            gene_info!(
                params.gene_name,
                "Reverse extension: {} nodes exceeds maximum of {}, stopping.",
                n_nodes,
                max_num_nodes
            );
            break;
        }

        if let Some(budget) = component_budget {
            if n_nodes.saturating_sub(nodes_at_start) >= budget {
                gene_info!(
                    params.gene_name,
                    "Reverse extension: component budget exhausted ({} nodes added, budget {}), moving to next component.",
                    n_nodes - nodes_at_start,
                    budget
                );
                break;
            }
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
            gene_info!(params.gene_name, "Evaluating reverse extension:");
            summarize_extension(&graph, "    ");
        }

        for node in unvisited {
            // Skip if already visited (may have been visited during this iteration)
            if graph[node].visited {
                continue;
            }

            let sub_kmer = graph[node].sub_kmer;

            trace!(
                "  Reverse extending node {} with sub_kmer {}.",
                node.index(),
                crate::kmer::kmer_to_seq(&sub_kmer, &(k - 1)),
            );

            // Find candidate predecessor kmers: (base << prefix_shift) | sub_kmer
            let mut candidate_kmers: [Option<u64>; 4] = [None; 4];
            let mut n_candidates = 0;

            for base in 0u64..4 {
                let kmer = (base << prefix_shift) | sub_kmer;

                if let Some(count) = kmer_counts.get_canonical(&kmer)
                    && count >= *min_count
                {
                    candidate_kmers[n_candidates] = Some(kmer);
                    n_candidates += 1;
                }
            }

            trace!("There are {} candidate predecessor kmers.", n_candidates);

            // If there are no candidates, the node is terminal
            if n_candidates == 0 {
                trace!(
                    "Marking node {} as terminal (no reverse candidates).",
                    node.index()
                );
                graph[node].is_terminal = true;
                graph[node].visited = true;
                continue;
            }

            for kmer in candidate_kmers.iter().flatten() {
                let prefix = kmer >> 2;

                // Self-loop detection
                if prefix == sub_kmer {
                    graph[node].is_terminal = true;
                    graph[node].visited = true;
                    trace!(
                        "Node {} reverse-extends to itself. Marking as terminal.",
                        node.index()
                    );
                    break;
                }

                if let Some(&existing_node) = node_lookup.get(&prefix) {
                    // Add edge from existing predecessor TO current node
                    if graph.find_edge(existing_node, node).is_none() {
                        let edge = get_dbedge(kmer, kmer_counts);
                        graph.add_edge(existing_node, node, edge);
                        if graph[existing_node].is_start {
                            gene_info!(
                                params.gene_name,
                                "Start node reached from end, reverse extension convergence."
                            );
                        }
                    }
                } else {
                    let edge = get_dbedge(kmer, kmer_counts);
                    let edge_count = edge.count;

                    // Skip edges with count far above the median
                    if (edge_count as f64) > (median_edge_count * params.high_coverage_ratio) {
                        trace!(
                            "Skipping high-coverage reverse edge (count {} > {:.0} x {:.0} median). kmer {}",
                            edge_count,
                            params.high_coverage_ratio,
                            median_edge_count,
                            crate::kmer::kmer_to_seq(kmer, &k),
                        );
                        continue;
                    }

                    let new_node = graph.add_node(DBNode {
                        sub_kmer: prefix,
                        is_start: false,
                        is_end: false,
                        is_terminal: false,
                        visited: false,
                    });
                    node_lookup.insert(prefix, new_node);

                    // Edge goes from new predecessor TO current node
                    graph.add_edge(new_node, node, edge);

                    trace!(
                        "Added reverse sub_kmer {} for new node {} with edge kmer count {}.",
                        crate::kmer::kmer_to_seq(&prefix, &(k - 1)),
                        new_node.index(),
                        edge_count
                    );

                    // Check max-length from end node
                    let path_length = get_path_length_from_end(&graph, new_node);

                    match path_length {
                        Ok(Some(path_length)) => {
                            trace!("Reverse path length to end is {}.", path_length);
                            if path_length + k > params.max_length {
                                graph[new_node].is_terminal = true;
                                graph[new_node].visited = true;
                                trace!(
                                    "Marking new reverse node {} as terminal (exceeds max-length from end).",
                                    new_node.index()
                                );
                            }
                        }
                        Ok(None) => {
                            graph[new_node].is_terminal = true;
                            trace!(
                                "Marking new reverse node {} as terminal (cycle detected).",
                                new_node.index()
                            );
                        }
                        Err(_) => {
                            // Node has no outgoing path to end — this can happen when
                            // reverse extension creates a node whose only outgoing
                            // neighbor was already visited and is not an end node.
                            // Leave it non-terminal so it continues extending in reverse.
                        }
                    }
                }
            }
            graph[node].visited = true;

            trace!(
                "Reverse extension: {} unvisited nodes remaining.",
                n_unvisited_nodes_in_graph(&graph),
            );
        }
    }

    Ok((graph, node_lookup))
}

/// Annotate each edge with its coverage ratio: count / global median.
/// High ratios (>> 1.0) indicate potentially repetitive edges.
/// Low ratios (<< 1.0) indicate potential sequencing errors.
pub(super) fn annotate_coverage_ratios(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    let mut counts: Vec<u32> = graph.edge_indices().map(|e| graph[e].count).collect();
    if counts.is_empty() {
        return;
    }
    counts.sort();
    let mid = counts.len() / 2;
    let median = if counts.len() % 2 == 0 {
        (counts[mid - 1] as f64 + counts[mid] as f64) / 2.0
    } else {
        counts[mid] as f64
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

/// Log debug diagnostics about start/end nodes in the extended graph.
pub(super) fn log_extended_graph_diagnostics(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
) {
    let mut start_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();
    let mut end_nodes_map: HashMap<NodeIndex, usize> = HashMap::new();

    for node in graph.node_indices() {
        if graph[node].is_start {
            start_nodes_map.insert(
                node,
                graph.neighbors_directed(node, Direction::Outgoing).count(),
            );
        }
        if graph[node].is_end {
            end_nodes_map.insert(
                node,
                graph.neighbors_directed(node, Direction::Incoming).count(),
            );
        }
    }

    for (node, edge_count) in &start_nodes_map {
        debug!(
            "  Start node {} with subkmer {} has {} edges",
            node.index(),
            crate::kmer::kmer_to_seq(&graph[*node].sub_kmer, &(kmer_counts.get_k() - 1)),
            edge_count
        );
    }

    for (node, edge_count) in &end_nodes_map {
        debug!(
            "  End node {} with subkmer {} has {} edges",
            node.index(),
            crate::kmer::kmer_to_seq(&graph[*node].sub_kmer, &(kmer_counts.get_k() - 1)),
            edge_count
        );
    }
}

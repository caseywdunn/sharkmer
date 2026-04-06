// pcr/pruning.rs — graph cleanup

use log::{debug, trace};
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;

use super::{DBEdge, DBNode};

/// Remove dead-end tips that are short (< k nodes) AND have low coverage
/// relative to the local median. These are almost certainly sequencing errors.
/// Tips with meaningful coverage are preserved — they may represent real
/// variants whose graph extension was truncated.
///
/// With reverse extension, dead-ends can occur in both directions:
/// - Forward tips: nodes with no outgoing edges (not end nodes)
/// - Backward tips: nodes with no incoming edges (not start nodes)
pub fn remove_low_coverage_tips(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    k: &usize,
    tip_coverage_fraction: f64,
) {
    // Compute the coverage reference once, from the full graph, before any
    // tips are removed. Using the pre-pruning median is both faster (no
    // O(n log n) resort per iteration) and more stable: as low-coverage
    // tips fall away, the median would drift upward, making the threshold
    // more aggressive on later iterations — an arbitrary dependence on
    // pruning order rather than on the amplicon's actual coverage profile.
    let median_count = global_median_edge_count(graph).unwrap_or(1.0);
    let min_tip_count = (median_count * tip_coverage_fraction).max(1.0);

    let mut removed = 1;
    while removed > 0 {
        removed = 0;

        let nodes_to_remove: Vec<NodeIndex> = graph
            .node_indices()
            .filter(|&node| {
                // Don't remove end nodes or start nodes
                if graph[node].is_end || graph[node].is_start {
                    return false;
                }

                let no_outgoing = graph.neighbors_directed(node, Direction::Outgoing).count() == 0;
                let no_incoming = graph.neighbors_directed(node, Direction::Incoming).count() == 0;

                // Must be a dead end in at least one direction
                if !no_outgoing && !no_incoming {
                    return false;
                }

                if no_outgoing {
                    // Forward tip: trace back to branch point
                    let tip_len = tip_length_backward(graph, node);
                    if tip_len >= *k {
                        return false;
                    }
                    let max_incoming_count = graph
                        .edges_directed(node, Direction::Incoming)
                        .map(|e| e.weight().count)
                        .max()
                        .unwrap_or(0);
                    if (max_incoming_count as f64) >= min_tip_count {
                        return false;
                    }
                }

                if no_incoming {
                    // Backward tip: trace forward to branch point
                    let tip_len = tip_length_forward(graph, node);
                    if tip_len >= *k {
                        return false;
                    }
                    let max_outgoing_count = graph
                        .edges_directed(node, Direction::Outgoing)
                        .map(|e| e.weight().count)
                        .max()
                        .unwrap_or(0);
                    if (max_outgoing_count as f64) >= min_tip_count {
                        return false;
                    }
                }

                true
            })
            .collect();

        for node in nodes_to_remove {
            trace!("Removing low-coverage tip node {}", node.index());
            graph.remove_node(node);
            removed += 1;
        }
    }
}

/// Trace back from a dead-end node (no outgoing edges) to count how many
/// nodes until a branch point (node with out-degree > 1) or a start node.
fn tip_length_backward(graph: &StableDiGraph<DBNode, DBEdge>, node: NodeIndex) -> usize {
    let mut length = 0;
    let mut current = node;
    loop {
        length += 1;
        let incoming: Vec<NodeIndex> = graph
            .neighbors_directed(current, Direction::Incoming)
            .collect();
        if incoming.len() != 1 {
            break;
        }
        let parent = incoming[0];
        if graph
            .neighbors_directed(parent, Direction::Outgoing)
            .count()
            > 1
        {
            break;
        }
        if graph[parent].is_start {
            break;
        }
        current = parent;
    }
    length
}

/// Trace forward from a dead-end node (no incoming edges) to count how many
/// nodes until a branch point (node with in-degree > 1) or an end node.
fn tip_length_forward(graph: &StableDiGraph<DBNode, DBEdge>, node: NodeIndex) -> usize {
    let mut length = 0;
    let mut current = node;
    loop {
        length += 1;
        let outgoing: Vec<NodeIndex> = graph
            .neighbors_directed(current, Direction::Outgoing)
            .collect();
        if outgoing.len() != 1 {
            break;
        }
        let child = outgoing[0];
        if graph.neighbors_directed(child, Direction::Incoming).count() > 1 {
            break;
        }
        if graph[child].is_end {
            break;
        }
        current = child;
    }
    length
}

fn global_median_edge_count(graph: &StableDiGraph<DBNode, DBEdge>) -> Option<f64> {
    let mut counts: Vec<u32> = graph.edge_indices().map(|e| graph[e].count).collect();
    if counts.is_empty() {
        return None;
    }
    counts.sort();
    let mid = counts.len() / 2;
    if counts.len() % 2 == 0 {
        Some((counts[mid - 1] as f64 + counts[mid] as f64) / 2.0)
    } else {
        Some(counts[mid] as f64)
    }
}

/// Remove nodes that cannot be part of any start-to-end path.
///
/// A node is reachable from start if a forward BFS from any start node visits
/// it. A node can reach an end if a backward BFS from any end node visits it.
/// Nodes that fail either criterion are removed. This also removes
/// disconnected components and subsumes `remove_orphan_nodes` behavior.
pub fn reachability_pruning(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    use std::collections::HashSet;

    // Forward reachability from all start nodes
    let mut forward_reachable: HashSet<NodeIndex> = HashSet::new();
    for node in graph.node_indices() {
        if graph[node].is_start {
            let mut stack = vec![node];
            while let Some(n) = stack.pop() {
                if forward_reachable.insert(n) {
                    for neighbor in graph.neighbors_directed(n, Direction::Outgoing) {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }

    // Backward reachability from all end nodes
    let mut backward_reachable: HashSet<NodeIndex> = HashSet::new();
    for node in graph.node_indices() {
        if graph[node].is_end {
            let mut stack = vec![node];
            while let Some(n) = stack.pop() {
                if backward_reachable.insert(n) {
                    for neighbor in graph.neighbors_directed(n, Direction::Incoming) {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }

    // Remove nodes not in the intersection
    let nodes_to_remove: Vec<NodeIndex> = graph
        .node_indices()
        .filter(|n| !forward_reachable.contains(n) || !backward_reachable.contains(n))
        .collect();

    if !nodes_to_remove.is_empty() {
        debug!(
            "Reachability pruning: removing {} of {} nodes not on any start-to-end path",
            nodes_to_remove.len(),
            graph.node_count()
        );
    }

    for node in nodes_to_remove {
        graph.remove_node(node);
    }
}

#[cfg(test)]
mod tests {
    use super::super::{DBEdge, DBNode};
    use super::*;

    fn mk_node(is_start: bool, is_end: bool) -> DBNode {
        DBNode {
            sub_kmer: 0,
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

    #[test]
    fn test_global_median_edge_count_empty() {
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        assert!(global_median_edge_count(&graph).is_none());
    }

    #[test]
    fn test_global_median_edge_count_odd() {
        let mut graph = StableDiGraph::new();
        let a = graph.add_node(mk_node(true, false));
        let b = graph.add_node(mk_node(false, false));
        let c = graph.add_node(mk_node(false, true));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, c, mk_edge(20));
        graph.add_edge(a, c, mk_edge(30));
        assert_eq!(global_median_edge_count(&graph).unwrap(), 20.0);
    }

    #[test]
    fn test_remove_low_coverage_tip_forward() {
        // start -> a -> b -> end
        //               \-> tip (low coverage, dead end)
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(true, false));
        let a = graph.add_node(mk_node(false, false));
        let b = graph.add_node(mk_node(false, false));
        let end = graph.add_node(mk_node(false, true));
        let tip = graph.add_node(mk_node(false, false)); // dead end, not end node

        graph.add_edge(start, a, mk_edge(100));
        graph.add_edge(a, b, mk_edge(100));
        graph.add_edge(b, end, mk_edge(100));
        graph.add_edge(b, tip, mk_edge(1)); // low coverage tip

        let k = 3;
        remove_low_coverage_tips(&mut graph, &k, 0.1);

        // tip should be removed (short, low coverage), 4 nodes remain
        assert_eq!(graph.node_count(), 4);
    }

    #[test]
    fn test_preserve_high_coverage_tip() {
        // Same structure but tip has high coverage — should be preserved
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(true, false));
        let a = graph.add_node(mk_node(false, false));
        let end = graph.add_node(mk_node(false, true));
        let tip = graph.add_node(mk_node(false, false));

        graph.add_edge(start, a, mk_edge(10));
        graph.add_edge(a, end, mk_edge(10));
        graph.add_edge(a, tip, mk_edge(10)); // same coverage as main path

        let k = 3;
        remove_low_coverage_tips(&mut graph, &k, 0.1);

        // tip should be preserved (high coverage)
        assert_eq!(graph.node_count(), 4);
    }

    #[test]
    fn test_reachability_pruning_removes_orphan() {
        // start -> a -> end, plus orphan disconnected node
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(true, false));
        let a = graph.add_node(mk_node(false, false));
        let end = graph.add_node(mk_node(false, true));
        let orphan = graph.add_node(mk_node(false, false));

        graph.add_edge(start, a, mk_edge(10));
        graph.add_edge(a, end, mk_edge(10));
        // orphan has no edges

        reachability_pruning(&mut graph);
        assert_eq!(graph.node_count(), 3);
        assert!(graph.contains_node(start));
        assert!(!graph.contains_node(orphan));
    }

    #[test]
    fn test_reachability_pruning_removes_dead_branch() {
        // start -> a -> end
        // start -> b (dead end, no path to end)
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(true, false));
        let a = graph.add_node(mk_node(false, false));
        let end = graph.add_node(mk_node(false, true));
        let b = graph.add_node(mk_node(false, false));

        graph.add_edge(start, a, mk_edge(10));
        graph.add_edge(a, end, mk_edge(10));
        graph.add_edge(start, b, mk_edge(10));
        // b has no path to end

        reachability_pruning(&mut graph);
        assert_eq!(graph.node_count(), 3);
        assert!(!graph.contains_node(b));
    }

    #[test]
    fn test_reachability_pruning_empty_graph() {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        reachability_pruning(&mut graph);
        assert_eq!(graph.node_count(), 0);
    }
}

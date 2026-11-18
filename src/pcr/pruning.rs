// pcr/pruning.rs — graph cleanup

use super::memory::buffers::{BudgetNodeSet, BudgetVec};
use crate::kmer::budget::ResourceBudget;
use anyhow::Result;
use log::{debug, trace};
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;
use std::sync::Arc;

use super::{DBEdge, DBNode};

/// Remove dead-end tips that are short (< k nodes) AND have low coverage
/// relative to the local median. These are almost certainly sequencing errors.
/// Tips with meaningful coverage are preserved — they may represent real
/// variants whose graph extension was truncated.
///
/// With reverse extension, dead-ends can occur in both directions:
/// - Forward tips: nodes with no outgoing edges (not end nodes)
/// - Backward tips: nodes with no incoming edges (not start nodes)
#[cfg(test)]
pub fn remove_low_coverage_tips(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    k: &usize,
    fraction: f64,
) {
    remove_low_coverage_tips_budgeted(graph, k, fraction, None).unwrap();
}

pub(super) fn remove_low_coverage_tips_budgeted(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    k: &usize,
    tip_coverage_fraction: f64,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<()> {
    // Compute the coverage reference once, from the full graph, before any
    // tips are removed. Using the pre-pruning median is both faster (no
    // O(n log n) resort per iteration) and more stable: as low-coverage
    // tips fall away, the median would drift upward, making the threshold
    // more aggressive on later iterations — an arbitrary dependence on
    // pruning order rather than on the amplicon's actual coverage profile.
    let median_count = global_median_edge_count(graph, budget)?.unwrap_or(1.0);
    let min_tip_count = (median_count * tip_coverage_fraction).max(1.0);

    let mut removed = 1;
    while removed > 0 {
        removed = 0;

        let nodes_to_remove = BudgetVec::collect(
            graph.node_indices().filter(|&node| {
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
            }),
            budget,
        )?;

        for &node in nodes_to_remove.iter() {
            trace!("Removing low-coverage tip node {}", node.index());
            graph.remove_node(node);
            removed += 1;
        }
    }
    Ok(())
}

/// Trace back from a dead-end node (no outgoing edges) to count how many
/// nodes until a branch point (node with out-degree > 1) or a start node.
fn tip_length_backward(graph: &StableDiGraph<DBNode, DBEdge>, node: NodeIndex) -> usize {
    let mut length = 0;
    let mut current = node;
    loop {
        length += 1;
        let mut incoming = graph.neighbors_directed(current, Direction::Incoming);
        let Some(parent) = incoming.next() else {
            break;
        };
        if incoming.next().is_some() {
            break;
        }
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
        let mut outgoing = graph.neighbors_directed(current, Direction::Outgoing);
        let Some(child) = outgoing.next() else {
            break;
        };
        if outgoing.next().is_some() {
            break;
        }
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

fn global_median_edge_count(
    graph: &StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Option<f64>> {
    let mut counts = BudgetVec::collect(graph.edge_indices().map(|e| graph[e].count), budget)?;
    Ok(super::graph::median_u64_slice(&mut counts))
}

/// Remove nodes that cannot be part of any start-to-end path.
///
/// A node is reachable from start if a forward BFS from any start node visits
/// it. A node can reach an end if a backward BFS from any end node visits it.
/// Nodes that fail either criterion are removed. This also removes
/// disconnected subgraphs and orphan branches in one pass.
///
/// Both BFS passes are run as a single union BFS — all start nodes are
/// pushed onto one shared stack before the loop begins, instead of running
/// a separate per-start BFS that restarted the stack each time. This is
/// equivalent to the previous behavior (each start's reachable set is
/// still merged into `forward_reachable`) but does one full sweep instead
/// of N starts × restart overhead, and reads more directly. Same shape
/// for the backward BFS over end nodes.
#[cfg(test)]
pub fn reachability_pruning(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    reachability_pruning_budgeted(graph, None).unwrap();
}

pub(super) fn reachability_pruning_budgeted(
    graph: &mut StableDiGraph<DBNode, DBEdge>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<()> {
    // Forward reachability: union BFS from all start nodes at once.
    let mut forward_reachable = BudgetNodeSet::new(budget);
    let mut stack =
        BudgetVec::collect(graph.node_indices().filter(|&n| graph[n].is_start), budget)?;
    while let Some(n) = stack.pop() {
        if forward_reachable.insert(n)? {
            for neighbor in graph.neighbors_directed(n, Direction::Outgoing) {
                stack.push(neighbor)?;
            }
        }
    }

    // Backward reachability: union BFS from all end nodes at once.
    let mut backward_reachable = BudgetNodeSet::new(budget);
    stack.clear();
    for node in graph.node_indices().filter(|&n| graph[n].is_end) {
        stack.push(node)?;
    }
    while let Some(n) = stack.pop() {
        if backward_reachable.insert(n)? {
            for neighbor in graph.neighbors_directed(n, Direction::Incoming) {
                stack.push(neighbor)?;
            }
        }
    }

    // Remove nodes not in the intersection
    let nodes_to_remove = BudgetVec::collect(
        graph
            .node_indices()
            .filter(|n| !forward_reachable.contains(n) || !backward_reachable.contains(n)),
        budget,
    )?;

    if !nodes_to_remove.is_empty() {
        debug!(
            "Reachability pruning: removing {} of {} nodes not on any start-to-end path",
            nodes_to_remove.len(),
            graph.node_count()
        );
    }

    for &node in nodes_to_remove.iter() {
        graph.remove_node(node);
    }
    Ok(())
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
    fn charged_reachability_matches_independent_path_oracle_with_stable_holes() {
        let budget = workspace_budget(1 << 20);
        let mut random = 159u64;
        for _ in 0..64 {
            let mut graph = StableDiGraph::new();
            for index in 0..11 {
                graph.add_node(mk_node(index == 0 || index == 4, index == 10));
            }
            graph.remove_node(NodeIndex::new(3));
            graph.remove_node(NodeIndex::new(7));
            let nodes: Vec<_> = graph.node_indices().collect();
            for _ in 0..18 {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let from = nodes[(random >> 32) as usize % nodes.len()];
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let to = nodes[(random >> 32) as usize % nodes.len()];
                graph.add_edge(from, to, mk_edge(10));
            }
            let expected: Vec<_> = graph
                .node_indices()
                .filter(|&node| {
                    nodes.iter().any(|&start| {
                        graph[start].is_start
                            && petgraph::algo::has_path_connecting(&graph, start, node, None)
                    }) && nodes.iter().any(|&end| {
                        graph[end].is_end
                            && petgraph::algo::has_path_connecting(&graph, node, end, None)
                    })
                })
                .collect();
            reachability_pruning_budgeted(&mut graph, Some(&budget)).unwrap();
            assert_eq!(graph.node_indices().collect::<Vec<_>>(), expected);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn pruning_workspace_refusal_is_typed_before_removing_nodes() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(true, false));
        let end = graph.add_node(mk_node(false, true));
        let tip = graph.add_node(mk_node(false, false));
        graph.add_edge(start, end, mk_edge(100));
        graph.add_edge(start, tip, mk_edge(1));
        let tiny = workspace_budget(1);
        for error in [
            remove_low_coverage_tips_budgeted(&mut graph, &3, 0.1, Some(&tiny)).unwrap_err(),
            reachability_pruning_budgeted(&mut graph, Some(&tiny)).unwrap_err(),
        ] {
            assert!(
                error
                    .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                    .is_some()
            );
            assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        }
        assert_eq!((graph.node_count(), graph.edge_count()), (3, 2));
        assert!(graph.node_weight(tip).is_some());
        let enough = workspace_budget(4096);
        remove_low_coverage_tips_budgeted(&mut graph, &3, 0.1, Some(&enough)).unwrap();
        assert!(graph.node_weight(tip).is_none());
        assert_eq!(enough.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn test_global_median_edge_count_empty() {
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        assert!(global_median_edge_count(&graph, None).unwrap().is_none());
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
        assert_eq!(global_median_edge_count(&graph, None).unwrap(), Some(20.0));
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

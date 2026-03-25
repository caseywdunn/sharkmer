// pcr/pruning.rs — graph pruning

use log::{debug, info, trace, warn};
use petgraph::Direction;
use petgraph::graph::NodeIndex;
use petgraph::stable_graph::StableDiGraph;

use super::graph::{EXTENSION_EVALUATION_DEPTH, descendants, get_descendants};
use super::{DBEdge, DBNode};

const EXTENSION_EVALUATION_DIFF: usize = 1;

pub(super) fn pop_balloons(graph: &mut StableDiGraph<DBNode, DBEdge>, k: &usize) {
    let mut to_clip: Vec<NodeIndex> = Vec::new();
    for node in graph.node_indices() {
        let d = descendants(graph, node, EXTENSION_EVALUATION_DEPTH);
        let n = d.len();
        // The maximum number of descendants would be 4^EXTENSION_EVALUATION_DEPTH
        if n > 4_usize.pow((EXTENSION_EVALUATION_DEPTH) as u32) {
            let seq = crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k - 1));
            warn!(
                "Node {} with sequence {} has {} descendants at a depth of {}. This exceeds the maximum of 4^{}={} that is expected",
                node.index(),
                seq,
                n,
                EXTENSION_EVALUATION_DEPTH,
                EXTENSION_EVALUATION_DEPTH,
                4_usize.pow((EXTENSION_EVALUATION_DEPTH) as u32)
            );

            // Get a vector of sequences of the descendants
            let mut seqs: Vec<String> = Vec::new();
            for descendant in d {
                seqs.push(crate::kmer::kmer_to_seq(
                    &graph[descendant].sub_kmer,
                    &(*k - 1),
                ));
            }

            debug!("  Sequences: {:?}", seqs);
        }

        if n > 4_usize.pow((EXTENSION_EVALUATION_DEPTH - EXTENSION_EVALUATION_DIFF) as u32) {
            to_clip.push(node);
            trace!(
                "  Node {} with sequence {} has {} descendants at a depth of {}, descendants will be clipped",
                node.index(),
                crate::kmer::kmer_to_seq(&graph[node].sub_kmer, &(*k - 1)),
                n,
                EXTENSION_EVALUATION_DEPTH
            );
        }
    }

    // Mark the clipped nodes as terminal
    for node in &to_clip {
        graph[*node].is_terminal = true;
    }

    // Vector to hold nodes to be pruned
    let mut to_prune: Vec<NodeIndex> = Vec::new();
    for node in &to_clip {
        // Node may have been removed already, so check if it is in graph
        if graph.node_weight(*node).is_some() {
            // prune away all the descendants of the node, but keep the node
            let mut descendants = get_descendants(graph, *node);
            to_prune.append(&mut descendants);
        }
    }

    if !to_prune.is_empty() {
        info!(
            "    Removing {} nodes descended from {} nodes with ballooning graph extension",
            to_prune.len(),
            to_clip.len()
        );
    }

    // StableDiGraph preserves indices on removal, so order doesn't matter,
    // but sorting for deterministic behavior
    to_prune.sort_by(|a, b| b.cmp(a));

    // Remove the nodes in to_prune
    for node in to_prune {
        graph.remove_node(node);
    }
}

// Iteratively remove nodes that do not have outgoing edges and are not end nodes
// These are terminal side branches.
pub fn remove_side_branches(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    let mut removed_nodes = 1;
    while removed_nodes > 0 {
        removed_nodes = 0;

        let nodes_to_remove: Vec<_> = graph
            .node_indices()
            .filter(|&node| {
                if graph[node].is_end {
                    false
                } else {
                    graph.neighbors_directed(node, Direction::Outgoing).count() == 0
                }
            })
            .collect();

        for node in nodes_to_remove {
            graph.remove_node(node);
            removed_nodes += 1;
        }
    }
}

// remove end nodes that have no incoming edges
pub fn remove_orphan_nodes(graph: &mut StableDiGraph<DBNode, DBEdge>) {
    let nodes_to_remove: Vec<NodeIndex> = graph
        .node_indices()
        .filter(|&node| {
            if graph[node].is_end {
                graph.neighbors_directed(node, Direction::Incoming).count() == 0
            } else {
                false
            }
        })
        .collect();

    for node in nodes_to_remove {
        graph.remove_node(node);
    }
}

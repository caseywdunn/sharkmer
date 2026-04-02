// pcr/bubble.rs — Read-aware bubble detection and resolution
//
// Detects simple bubbles in the de Bruijn graph (diverge at one node,
// converge at another) and ranks alternative branches using read-support
// and branch-point phasing data from threading annotations.

use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;
use std::collections::{HashMap, HashSet};

use super::threading::{BranchLink, ThreadingAnnotations};
use super::{DBEdge, DBNode};

/// Maximum depth to search for bubble convergence from a divergence point.
const MAX_BUBBLE_DEPTH: usize = 50;

/// Result of resolving one bubble: branches ranked by read support.
#[derive(Debug)]
#[allow(dead_code)]
pub struct BubbleResolution {
    /// The node where branches diverge
    pub source: NodeIndex,
    /// The node where branches converge
    #[allow(dead_code)]
    pub sink: NodeIndex,
    /// Branches ranked by total read support (best first)
    pub ranked_branches: Vec<BranchRanking>,
}

/// A single branch through a bubble with its support evidence.
#[derive(Debug)]
pub struct BranchRanking {
    /// Edges along this branch
    pub edges: Vec<EdgeIndex>,
    /// Total read support summed across all edges in this branch
    pub total_read_support: u32,
    /// Branch-point phasing support (reads linking this branch's entry to its exit)
    pub phasing_support: u32,
}

/// Detect simple bubbles and resolve them using threading annotations.
///
/// A "simple bubble" is a pattern where:
/// - A source node has out-degree >= 2
/// - Multiple branches emerge and converge at the same sink node
/// - Each branch is a short linear path (no internal branching)
///
/// Returns edge preferences: a map from EdgeIndex to a priority score.
/// Higher scores indicate more read-supported edges at branch points.
pub fn resolve_bubbles(
    graph: &StableDiGraph<DBNode, DBEdge>,
    annotations: &ThreadingAnnotations,
) -> HashMap<EdgeIndex, f64> {
    let mut edge_preferences: HashMap<EdgeIndex, f64> = HashMap::new();

    let bubbles = detect_simple_bubbles(graph);

    for bubble in &bubbles {
        let rankings = rank_branches(graph, bubble, annotations);

        // Assign preference scores: best branch gets 1.0, worst gets lowest
        if rankings.len() >= 2 {
            let max_support = rankings
                .iter()
                .map(|r| r.total_read_support + r.phasing_support)
                .max()
                .unwrap_or(0);

            for ranking in &rankings {
                let support = ranking.total_read_support + ranking.phasing_support;
                let preference = if max_support > 0 {
                    support as f64 / max_support as f64
                } else {
                    1.0 // no read data: all branches equal
                };

                for &edge in &ranking.edges {
                    edge_preferences.insert(edge, preference);
                }
            }
        }
    }

    edge_preferences
}

/// A raw detected bubble before ranking.
struct Bubble {
    source: NodeIndex,
    #[allow(dead_code)]
    sink: NodeIndex,
    branches: Vec<Vec<EdgeIndex>>,
}

/// Detect simple bubbles in the graph.
///
/// For each node with out-degree >= 2, trace each outgoing branch
/// to see if they converge at a common sink within MAX_BUBBLE_DEPTH.
fn detect_simple_bubbles(graph: &StableDiGraph<DBNode, DBEdge>) -> Vec<Bubble> {
    let mut bubbles: Vec<Bubble> = Vec::new();

    for source in graph.node_indices() {
        let outgoing: Vec<_> = graph.edges_directed(source, Direction::Outgoing).collect();

        if outgoing.len() < 2 {
            continue;
        }

        // Trace each branch to find convergence
        let mut branch_endpoints: HashMap<NodeIndex, Vec<Vec<EdgeIndex>>> = HashMap::new();

        for edge_ref in &outgoing {
            let first_edge = edge_ref.id();
            let mut current = edge_ref.target();
            let mut branch_edges = vec![first_edge];
            let mut visited = HashSet::new();
            visited.insert(source);
            let mut depth = 0;

            // Follow linear path (single outgoing edge) until branch or convergence
            loop {
                if depth >= MAX_BUBBLE_DEPTH {
                    break;
                }
                depth += 1;

                if visited.contains(&current) {
                    break; // cycle
                }
                visited.insert(current);

                let next_edges: Vec<_> =
                    graph.edges_directed(current, Direction::Outgoing).collect();

                if next_edges.len() == 1 {
                    // Linear path: keep following
                    branch_edges.push(next_edges[0].id());
                    current = next_edges[0].target();
                } else {
                    // Branch point or dead end: record endpoint
                    break;
                }
            }

            // Record where this branch ended up
            branch_endpoints
                .entry(current)
                .or_default()
                .push(branch_edges);
        }

        // A bubble exists if 2+ branches converge at the same sink
        for (sink, branches) in branch_endpoints {
            if branches.len() >= 2 && sink != source {
                bubbles.push(Bubble {
                    source,
                    sink,
                    branches,
                });
            }
        }
    }

    bubbles
}

/// Rank branches of a bubble by read support and phasing evidence.
fn rank_branches(
    graph: &StableDiGraph<DBNode, DBEdge>,
    bubble: &Bubble,
    annotations: &ThreadingAnnotations,
) -> Vec<BranchRanking> {
    let mut rankings: Vec<BranchRanking> = Vec::new();

    for branch_edges in &bubble.branches {
        // Sum read support across branch edges
        let total_read_support: u32 = branch_edges
            .iter()
            .map(|edge_idx| {
                annotations
                    .edge_support
                    .get(edge_idx)
                    .map_or(0, |s| s.read_support_total)
            })
            .sum();

        // Check for phasing support: reads that link the entry edge to
        // the branch path through the divergence point
        let phasing_support = if let Some(&first_edge) = branch_edges.first() {
            // Look for branch links entering source → first_edge
            let incoming_to_source: Vec<EdgeIndex> = graph
                .edges_directed(bubble.source, Direction::Incoming)
                .map(|e| e.id())
                .collect();

            let mut phasing_count: u32 = 0;
            for &incoming in &incoming_to_source {
                let link = BranchLink {
                    incoming_edge: incoming,
                    outgoing_edge: first_edge,
                };
                if let Some(&count) = annotations.branch_links.get(&link) {
                    phasing_count += count;
                }
            }
            phasing_count
        } else {
            0
        };

        rankings.push(BranchRanking {
            edges: branch_edges.clone(),
            total_read_support,
            phasing_support,
        });
    }

    // Sort by total evidence (descending)
    rankings.sort_by(|a, b| {
        let a_total = a.total_read_support + a.phasing_support;
        let b_total = b.total_read_support + b.phasing_support;
        b_total.cmp(&a_total)
    });

    rankings
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bubble_graph() -> StableDiGraph<DBNode, DBEdge> {
        // Create a diamond bubble:
        //   n0 --e0--> n1 --e2--> n3
        //   n0 --e1--> n2 --e3--> n3
        let mut graph = StableDiGraph::new();
        let n0 = graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n2 = graph.add_node(DBNode {
            sub_kmer: 2,
            is_start: false,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n3 = graph.add_node(DBNode {
            sub_kmer: 3,
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });

        graph.add_edge(
            n0,
            n1,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n0,
            n2,
            DBEdge {
                count: 3,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n1,
            n3,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );
        graph.add_edge(
            n2,
            n3,
            DBEdge {
                count: 3,
                coverage_ratio: 1.0,
            },
        );

        graph
    }

    #[test]
    fn test_detect_simple_bubble() {
        let graph = make_bubble_graph();
        let bubbles = detect_simple_bubbles(&graph);
        assert_eq!(bubbles.len(), 1);
        assert_eq!(bubbles[0].branches.len(), 2);
    }

    #[test]
    fn test_resolve_with_read_support() {
        let graph = make_bubble_graph();

        let mut annotations = ThreadingAnnotations {
            edge_support: HashMap::new(),
            branch_links: HashMap::new(),
            paired_links: Vec::new(),
        };

        // Give branch 1 (through n1) more read support
        let e0 = graph.edge_indices().nth(0).unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();
        let e2 = graph.edge_indices().nth(2).unwrap();
        let e3 = graph.edge_indices().nth(3).unwrap();

        annotations.edge_support.insert(
            e0,
            super::super::threading::EdgeReadSupport {
                read_support_total: 10,
                read_support_unambiguous: 8,
            },
        );
        annotations.edge_support.insert(
            e2,
            super::super::threading::EdgeReadSupport {
                read_support_total: 10,
                read_support_unambiguous: 8,
            },
        );
        annotations.edge_support.insert(
            e1,
            super::super::threading::EdgeReadSupport {
                read_support_total: 2,
                read_support_unambiguous: 1,
            },
        );
        annotations.edge_support.insert(
            e3,
            super::super::threading::EdgeReadSupport {
                read_support_total: 2,
                read_support_unambiguous: 1,
            },
        );

        let prefs = resolve_bubbles(&graph, &annotations);

        // Branch through n1 (e0, e2) should have higher preference
        assert!(prefs[&e0] > prefs[&e1]);
        assert!(prefs[&e2] > prefs[&e3]);
    }

    #[test]
    fn test_no_bubble_linear() {
        // Linear graph has no bubbles
        let mut graph = StableDiGraph::new();
        let n0 = graph.add_node(DBNode {
            sub_kmer: 0,
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 1,
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });
        graph.add_edge(
            n0,
            n1,
            DBEdge {
                count: 5,
                coverage_ratio: 1.0,
            },
        );

        let bubbles = detect_simple_bubbles(&graph);
        assert!(bubbles.is_empty());
    }
}

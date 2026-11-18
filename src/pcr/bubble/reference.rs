//! Frozen pre-ownership bubble algorithm for differential tests only.
#![allow(dead_code)]

// pcr/bubble.rs — Read-aware bubble detection and resolution
//
// Detects simple bubbles in the de Bruijn graph (diverge at one node,
// converge at another) and ranks alternative branches using read-support
// and branch-point phasing data from threading annotations.

use ahash::{AHashMap, AHashSet};
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;

use super::{DBEdge, DBNode};
use crate::pcr::threading::reference::{BranchLink, ThreadingAnnotations};

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
) -> AHashMap<EdgeIndex, f64> {
    let mut edge_preferences: AHashMap<EdgeIndex, f64> = AHashMap::new();

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

        // Trace each branch to find convergence. Branches that terminate
        // by running out of depth budget (MAX_BUBBLE_DEPTH) are discarded
        // rather than recorded: they may or may not converge somewhere
        // past the horizon, and pretending their current position is a
        // bubble sink would fabricate a bubble from two unrelated paths
        // that happen to land on the same node at the depth cutoff.
        let mut branch_endpoints: AHashMap<NodeIndex, Vec<Vec<EdgeIndex>>> = AHashMap::new();

        for edge_ref in &outgoing {
            let first_edge = edge_ref.id();
            let mut current = edge_ref.target();
            let mut branch_edges = vec![first_edge];
            let mut visited: AHashSet<NodeIndex> = AHashSet::new();
            visited.insert(source);
            let mut depth = 0;
            let mut terminated_naturally = false;

            // Follow linear path (single outgoing edge) until branch or convergence
            loop {
                if depth >= MAX_BUBBLE_DEPTH {
                    // Hit the depth limit without reaching a branch point,
                    // dead end, or cycle. We cannot classify this branch.
                    break;
                }
                depth += 1;

                if visited.contains(&current) {
                    // Cycle detected: the branch looped back to an
                    // already-visited node. This is not a real sink,
                    // so drop the branch rather than treat the loop
                    // point as a convergence candidate.
                    break;
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
                    terminated_naturally = true;
                    break;
                }
            }

            // Only record branches that terminated naturally — at a branch
            // point, dead end, or detected cycle. Depth-limited branches
            // are dropped.
            if terminated_naturally {
                branch_endpoints
                    .entry(current)
                    .or_default()
                    .push(branch_edges);
            }
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

    // Sort by total evidence (descending). When two branches tie on
    // evidence, we need a tiebreaker that is stable across runs — otherwise
    // the winner depends on HashMap iteration order (via branch_endpoints
    // in detect_simple_bubbles), which changes across hash seeds and leaks
    // non-determinism into the output sequence. Use the canonical kmer of
    // each branch's first edge as the tiebreaker: it is a content-derived
    // value that stays the same regardless of graph insertion order.
    rankings.sort_by(|a, b| {
        let a_total = a.total_read_support + a.phasing_support;
        let b_total = b.total_read_support + b.phasing_support;
        b_total.cmp(&a_total).then_with(|| {
            let a_key = branch_sort_key(graph, &a.edges);
            let b_key = branch_sort_key(graph, &b.edges);
            a_key.cmp(&b_key)
        })
    });

    rankings
}

/// Deterministic sort key for a branch: the reconstructed directional
/// kmer of the branch's first edge. Two branches of the same bubble
/// diverge at a shared source node with distinct outgoing edges, so
/// their first edges have different target sub_kmers and therefore
/// different directional kmers — no canonicalization is needed to make
/// the key unique within a single bubble. The value is content-derived
/// (depends only on node sub_kmers, which are determined by the input
/// kmer sequences), so it is stable across runs regardless of HashMap
/// iteration order during graph construction.
fn branch_sort_key(graph: &StableDiGraph<DBNode, DBEdge>, edges: &[EdgeIndex]) -> u64 {
    let first = match edges.first() {
        Some(&e) => e,
        None => return 0,
    };
    let (src, tgt) = graph
        .edge_endpoints(first)
        .expect("edge must exist in graph");
    (graph[src].sub_kmer << 2) | (graph[tgt].sub_kmer & 3)
}

// pcr/threading.rs — Read threading through assembly graphs
//
// Maps reads to graph edges via maximal contiguous runs of adjacent
// graph kmers. Annotates edges with read support counts and records
// branch-point phasing links.

use petgraph::Direction;
use petgraph::graph::EdgeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use std::collections::HashMap;

use super::{DBEdge, DBNode};
use crate::io::{Mate, ReadRecord};
use crate::kmer::encoding::{kmers_from_ascii, revcomp_kmer};

/// Per-edge read support annotation.
#[derive(Debug, Clone, Default)]
pub struct EdgeReadSupport {
    /// Total reads whose contiguous run includes this edge
    pub read_support_total: u32,
    /// Reads mapping to a single unbranched path through this edge
    pub read_support_unambiguous: u32,
}

/// A branch-point phasing observation: a read's contiguous run passed
/// through a branch point, linking an incoming edge to an outgoing edge.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct BranchLink {
    pub incoming_edge: EdgeIndex,
    pub outgoing_edge: EdgeIndex,
}

/// A paired-end phasing link: both mates of a pair map to the same
/// amplicon graph, providing long-range phasing information.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PairedEndLink {
    /// Edges covered by R1 mate
    pub r1_edges: Vec<EdgeIndex>,
    /// Edges covered by R2 mate
    pub r2_edges: Vec<EdgeIndex>,
}

/// Complete read threading annotations for a graph.
#[derive(Debug, Clone)]
pub struct ThreadingAnnotations {
    /// Per-edge read support counts
    pub edge_support: HashMap<EdgeIndex, EdgeReadSupport>,
    /// Branch-point phasing: (incoming, outgoing) -> count
    pub branch_links: HashMap<BranchLink, u32>,
    /// Paired-end links connecting distant regions of the graph
    pub paired_links: Vec<PairedEndLink>,
}

impl ThreadingAnnotations {
    fn new() -> Self {
        ThreadingAnnotations {
            edge_support: HashMap::new(),
            branch_links: HashMap::new(),
            paired_links: Vec::new(),
        }
    }
}

/// A contiguous run of graph edges that a single read maps to.
struct ReadRun {
    edges: Vec<EdgeIndex>,
}

/// Thread reads through a graph, producing annotations.
/// This function does NOT modify the graph.
///
/// For each read:
/// 1. Extract kmers from the read sequence
/// 2. Look up each kmer as a graph edge
/// 3. Find maximal contiguous runs of adjacent edges
/// 4. For each run, annotate edges and detect branch-point crossings
pub fn thread_reads(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[&ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    let mut annotations = ThreadingAnnotations::new();

    // Build edge lookup: canonical kmer -> EdgeIndex
    let edge_lookup = build_edge_lookup(graph, k);

    for read in reads {
        let kmers = match kmers_from_ascii(&read.sequence, k) {
            Ok(k) => k,
            Err(_) => continue, // skip reads with invalid characters
        };

        // Map kmers to edge indices
        let edge_hits: Vec<Option<EdgeIndex>> = kmers
            .iter()
            .map(|kmer| edge_lookup.get(kmer).copied())
            .collect();

        // Find maximal contiguous runs
        let runs = find_contiguous_runs(&edge_hits, graph);

        // Annotate edges from each run
        for run in &runs {
            let is_unambiguous = is_run_unambiguous(graph, &run.edges);

            for &edge_idx in &run.edges {
                let support = annotations.edge_support.entry(edge_idx).or_default();
                support.read_support_total += 1;
                if is_unambiguous {
                    support.read_support_unambiguous += 1;
                }
            }

            // Record branch-point phasing links
            record_branch_links(graph, &run.edges, &mut annotations.branch_links);
        }
    }

    annotations
}

/// Thread reads with paired-end phasing support.
/// First performs standard threading, then identifies read pairs
/// where both mates map to the same graph.
pub fn thread_reads_paired(
    graph: &StableDiGraph<DBNode, DBEdge>,
    reads: &[&ReadRecord],
    k: usize,
) -> ThreadingAnnotations {
    let mut annotations = ThreadingAnnotations::new();
    let edge_lookup = build_edge_lookup(graph, k);

    // Group reads by pair index for paired-end phasing
    // Pair index = read_index / 2 for alternating R1/R2 reads
    let mut pair_runs: HashMap<u64, (Vec<EdgeIndex>, Vec<EdgeIndex>)> = HashMap::new();

    for read in reads {
        let kmers = match kmers_from_ascii(&read.sequence, k) {
            Ok(k) => k,
            Err(_) => continue,
        };

        let edge_hits: Vec<Option<EdgeIndex>> = kmers
            .iter()
            .map(|kmer| edge_lookup.get(kmer).copied())
            .collect();

        let runs = find_contiguous_runs(&edge_hits, graph);

        // Collect all edges from all runs for this read
        let mut all_edges: Vec<EdgeIndex> = Vec::new();

        for run in &runs {
            let is_unambiguous = is_run_unambiguous(graph, &run.edges);

            for &edge_idx in &run.edges {
                let support = annotations.edge_support.entry(edge_idx).or_default();
                support.read_support_total += 1;
                if is_unambiguous {
                    support.read_support_unambiguous += 1;
                }
            }

            record_branch_links(graph, &run.edges, &mut annotations.branch_links);
            all_edges.extend_from_slice(&run.edges);
        }

        // Track paired-end edges
        if !all_edges.is_empty() {
            match read.mate {
                Mate::R1 => {
                    let pair_idx = read.index / 2;
                    pair_runs.entry(pair_idx).or_default().0 = all_edges;
                }
                Mate::R2 => {
                    let pair_idx = read.index / 2;
                    pair_runs.entry(pair_idx).or_default().1 = all_edges;
                }
                Mate::Unpaired => {} // no pairing info
            }
        }
    }

    // Create paired-end links where both mates map
    for (_pair_idx, (r1_edges, r2_edges)) in pair_runs {
        if !r1_edges.is_empty() && !r2_edges.is_empty() {
            annotations
                .paired_links
                .push(PairedEndLink { r1_edges, r2_edges });
        }
    }

    annotations
}

/// Build a lookup table from canonical kmer -> EdgeIndex.
fn build_edge_lookup(graph: &StableDiGraph<DBNode, DBEdge>, k: usize) -> HashMap<u64, EdgeIndex> {
    let mut lookup: HashMap<u64, EdgeIndex> = HashMap::new();

    for edge_ref in graph.edge_references() {
        let kmer = super::graph::reconstruct_edge_kmer(graph, edge_ref.id());
        // Store canonical form (what kmers_from_ascii returns)
        let rc = revcomp_kmer(&kmer, &k);
        let canonical = kmer.min(rc);
        lookup.insert(canonical, edge_ref.id());
    }

    lookup
}

/// Find maximal contiguous runs of edges from a sequence of edge lookups.
/// Two consecutive edges are "contiguous" if they share a node
/// (target of first == source of second) in the graph.
fn find_contiguous_runs(
    edge_hits: &[Option<EdgeIndex>],
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> Vec<ReadRun> {
    let mut runs: Vec<ReadRun> = Vec::new();
    let mut current_run: Vec<EdgeIndex> = Vec::new();

    for hit in edge_hits {
        match *hit {
            Some(edge_idx) => {
                if let Some(&prev_edge) = current_run.last() {
                    // Check adjacency: target of previous edge == source of current edge
                    let (_, prev_target) = graph
                        .edge_endpoints(prev_edge)
                        .expect("edge must exist in graph");
                    let (curr_source, _) = graph
                        .edge_endpoints(edge_idx)
                        .expect("edge must exist in graph");

                    if prev_target == curr_source {
                        current_run.push(edge_idx);
                    } else {
                        // Break: start a new run
                        if !current_run.is_empty() {
                            runs.push(ReadRun {
                                edges: std::mem::take(&mut current_run),
                            });
                        }
                        current_run.push(edge_idx);
                    }
                } else {
                    current_run.push(edge_idx);
                }
            }
            None => {
                // No edge match: end current run
                if !current_run.is_empty() {
                    runs.push(ReadRun {
                        edges: std::mem::take(&mut current_run),
                    });
                }
            }
        }
    }

    // Flush final run
    if !current_run.is_empty() {
        runs.push(ReadRun { edges: current_run });
    }

    runs
}

/// Check if a run is "unambiguous": every intermediate node has
/// in-degree <= 1 and out-degree <= 1 in the full graph.
fn is_run_unambiguous(graph: &StableDiGraph<DBNode, DBEdge>, edges: &[EdgeIndex]) -> bool {
    if edges.len() < 2 {
        return true; // single edge is trivially unambiguous
    }

    // Check intermediate nodes (shared between consecutive edges)
    for window in edges.windows(2) {
        let (_, node) = graph.edge_endpoints(window[0]).expect("edge must exist");
        let in_deg = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_deg = graph.neighbors_directed(node, Direction::Outgoing).count();
        if in_deg > 1 || out_deg > 1 {
            return false;
        }
    }

    true
}

/// Record branch-point phasing links for consecutive edge pairs
/// at branch points (in-degree > 1 or out-degree > 1).
fn record_branch_links(
    graph: &StableDiGraph<DBNode, DBEdge>,
    edges: &[EdgeIndex],
    branch_links: &mut HashMap<BranchLink, u32>,
) {
    for window in edges.windows(2) {
        let incoming = window[0];
        let outgoing = window[1];

        // The shared node is the target of incoming / source of outgoing
        let (_, node) = graph.edge_endpoints(incoming).expect("edge must exist");

        let in_deg = graph.neighbors_directed(node, Direction::Incoming).count();
        let out_deg = graph.neighbors_directed(node, Direction::Outgoing).count();

        if in_deg > 1 || out_deg > 1 {
            let link = BranchLink {
                incoming_edge: incoming,
                outgoing_edge: outgoing,
            };
            *branch_links.entry(link).or_insert(0) += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::stable_graph::StableDiGraph;

    fn make_test_graph() -> StableDiGraph<DBNode, DBEdge> {
        // Create a simple linear graph: n0 --e0--> n1 --e1--> n2
        let mut graph = StableDiGraph::new();
        let n0 = graph.add_node(DBNode {
            sub_kmer: 0b0000, // AA
            is_start: true,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 0b0001, // AC
            is_start: false,
            is_end: false,
            is_terminal: false,
            visited: false,
        });
        let n2 = graph.add_node(DBNode {
            sub_kmer: 0b0110, // CG
            is_start: false,
            is_end: true,
            is_terminal: false,
            visited: false,
        });

        // Edge e0: kmer AAC (0b000001), connects n0->n1
        graph.add_edge(
            n0,
            n1,
            DBEdge {
                count: 10,
                coverage_ratio: 1.0,
            },
        );
        // Edge e1: kmer ACG (0b000110), connects n1->n2
        graph.add_edge(
            n1,
            n2,
            DBEdge {
                count: 8,
                coverage_ratio: 1.0,
            },
        );

        graph
    }

    #[test]
    fn test_build_edge_lookup() {
        let graph = make_test_graph();
        let lookup = build_edge_lookup(&graph, 3);
        // Should have 2 entries (one per edge)
        assert_eq!(lookup.len(), 2);
    }

    #[test]
    fn test_contiguous_run_linear() {
        let graph = make_test_graph();
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // Simulate two adjacent edges
        let hits = vec![Some(e0), Some(e1)];
        let runs = find_contiguous_runs(&hits, &graph);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].edges.len(), 2);
    }

    #[test]
    fn test_contiguous_run_gap() {
        let graph = make_test_graph();
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // Gap between edges
        let hits = vec![Some(e0), None, Some(e1)];
        let runs = find_contiguous_runs(&hits, &graph);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].edges.len(), 1);
        assert_eq!(runs[1].edges.len(), 1);
    }

    #[test]
    fn test_unambiguous_linear() {
        let graph = make_test_graph();
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // Linear graph: all intermediate nodes have degree 1
        assert!(is_run_unambiguous(&graph, &[e0, e1]));
    }

    #[test]
    fn test_branch_point_detection() {
        let mut graph = make_test_graph();

        // Add a branch at n1: n1 -> n3
        let n3 = graph.add_node(DBNode {
            sub_kmer: 0b1010,
            is_start: false,
            is_end: false,
            is_terminal: true,
            visited: false,
        });
        let _e2 = graph.add_edge(
            petgraph::graph::NodeIndex::new(1), // n1
            n3,
            DBEdge {
                count: 3,
                coverage_ratio: 0.3,
            },
        );

        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // n1 now has out-degree 2, so the run is NOT unambiguous
        assert!(!is_run_unambiguous(&graph, &[e0, e1]));

        // Should record a branch link
        let mut branch_links = HashMap::new();
        record_branch_links(&graph, &[e0, e1], &mut branch_links);
        assert_eq!(branch_links.len(), 1);
    }
}

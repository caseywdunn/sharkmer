// pcr/threading.rs — Read threading through assembly graphs
//
// Maps reads to graph edges via maximal contiguous runs of adjacent
// graph kmers. Annotates edges with read support counts and records
// branch-point phasing links.

use petgraph::Direction;
use petgraph::graph::EdgeIndex;
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use smallvec::{SmallVec, smallvec};
use std::collections::HashMap;

use super::{DBEdge, DBNode};
use crate::io::{Mate, ReadRecord};
use crate::kmer::encoding::{kmers_from_ascii, revcomp_kmer};

/// Candidates for a canonical-kmer lookup. At most two directional edges can
/// share a canonical key (the kmer and its reverse complement), because
/// directional kmers are unique in the graph (nodes are deduplicated by
/// sub_kmer, parallel edges are forbidden in `extend_graph`). A SmallVec
/// with inline capacity 2 stores both without heap allocation.
type EdgeCandidates = SmallVec<[EdgeIndex; 2]>;

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

    // Build edge lookup: canonical kmer -> candidate EdgeIndices
    let edge_lookup = build_edge_lookup(graph, k);

    for read in reads {
        let kmers = match kmers_from_ascii(&read.sequence, k) {
            Ok(k) => k,
            Err(_) => continue, // skip reads with invalid characters
        };

        let runs = find_contiguous_runs(&kmers, &edge_lookup, graph);

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

        let runs = find_contiguous_runs(&kmers, &edge_lookup, graph);

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

/// Build a lookup table from canonical kmer -> candidate EdgeIndices.
///
/// A canonical key can legitimately map to up to two directional edges when
/// the amplicon contains an inverted repeat of length >= k: one edge carries
/// a kmer X, the other carries rc(X), and they share the canonical key
/// `min(X, rc(X))`. This is common in rRNA and mitochondrial targets where
/// secondary-structure stems are inverted repeats. Both candidates must be
/// retained so that threading can disambiguate based on graph adjacency to
/// the previous edge in the read's run.
fn build_edge_lookup(
    graph: &StableDiGraph<DBNode, DBEdge>,
    k: usize,
) -> HashMap<u64, EdgeCandidates> {
    let mut lookup: HashMap<u64, EdgeCandidates> = HashMap::new();

    for edge_ref in graph.edge_references() {
        let kmer = super::graph::reconstruct_edge_kmer(graph, edge_ref.id());
        let rc = revcomp_kmer(&kmer, &k);
        let canonical = kmer.min(rc);
        lookup
            .entry(canonical)
            .and_modify(|cands| cands.push(edge_ref.id()))
            .or_insert_with(|| smallvec![edge_ref.id()]);
    }

    lookup
}

/// Pick the best edge candidate for the current lookup given the previous
/// edge in the run. When there is only one candidate, that candidate is
/// returned. When there are two (inverted-repeat collision), prefer the one
/// whose source node equals the previous edge's target (i.e. the one that
/// extends the current run). If neither is adjacent (start of run, or the
/// previous edge was in a different part of the graph), return the first.
///
/// Returning the first arbitrarily when no disambiguation is possible is
/// self-correcting: the adjacency check in `find_contiguous_runs` will
/// break the run on the next kmer if the wrong choice was made, producing
/// at worst one spurious singleton run of length 1.
fn resolve_candidates(
    candidates: &EdgeCandidates,
    prev_edge: Option<EdgeIndex>,
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> EdgeIndex {
    debug_assert!(!candidates.is_empty(), "lookup entry must be non-empty");
    if candidates.len() == 1 {
        return candidates[0];
    }
    if let Some(prev) = prev_edge {
        let (_, prev_target) = graph
            .edge_endpoints(prev)
            .expect("previous edge must exist in graph");
        for &cand in candidates {
            let (cand_source, _) = graph
                .edge_endpoints(cand)
                .expect("candidate edge must exist in graph");
            if cand_source == prev_target {
                return cand;
            }
        }
    }
    candidates[0]
}

/// Walk a read's kmer sequence against the edge lookup, producing maximal
/// contiguous runs. Each kmer is resolved against the previous edge's
/// adjacency to handle inverted-repeat canonical collisions.
fn find_contiguous_runs(
    kmers: &[u64],
    edge_lookup: &HashMap<u64, EdgeCandidates>,
    graph: &StableDiGraph<DBNode, DBEdge>,
) -> Vec<ReadRun> {
    let mut runs: Vec<ReadRun> = Vec::new();
    let mut current_run: Vec<EdgeIndex> = Vec::new();

    for kmer in kmers {
        let candidates = match edge_lookup.get(kmer) {
            Some(c) => c,
            None => {
                // No edge match: end current run
                if !current_run.is_empty() {
                    runs.push(ReadRun {
                        edges: std::mem::take(&mut current_run),
                    });
                }
                continue;
            }
        };

        let edge_idx = resolve_candidates(candidates, current_run.last().copied(), graph);

        if let Some(&prev_edge) = current_run.last() {
            let (_, prev_target) = graph
                .edge_endpoints(prev_edge)
                .expect("edge must exist in graph");
            let (curr_source, _) = graph
                .edge_endpoints(edge_idx)
                .expect("edge must exist in graph");

            if prev_target == curr_source {
                current_run.push(edge_idx);
            } else {
                // Adjacency broke: flush and start a new run
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

    // Flush final run
    if !current_run.is_empty() {
        runs.push(ReadRun { edges: current_run });
    }

    runs
}

/// Check if a run is "unambiguous": every *intermediate* node has
/// in-degree <= 1 and out-degree <= 1 in the full graph.
/// Entry and exit nodes are intentionally excluded — branch points at
/// run boundaries are expected and handled by the caller.
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
        });
        let n1 = graph.add_node(DBNode {
            sub_kmer: 0b0001, // AC
            is_start: false,
            is_end: false,
            is_terminal: false,
        });
        let n2 = graph.add_node(DBNode {
            sub_kmer: 0b0110, // CG
            is_start: false,
            is_end: true,
            is_terminal: false,
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

    /// Helper: compute the canonical form of an edge's kmer via the graph.
    fn canonical_edge_kmer(
        graph: &StableDiGraph<DBNode, DBEdge>,
        edge: EdgeIndex,
        k: usize,
    ) -> u64 {
        let kmer = super::super::graph::reconstruct_edge_kmer(graph, edge);
        let rc = revcomp_kmer(&kmer, &k);
        kmer.min(rc)
    }

    #[test]
    fn test_build_edge_lookup() {
        let graph = make_test_graph();
        let lookup = build_edge_lookup(&graph, 3);
        // Should have 2 entries (one per edge) in the linear test graph
        assert_eq!(lookup.len(), 2);
        // Each entry has exactly one candidate
        for cands in lookup.values() {
            assert_eq!(cands.len(), 1);
        }
    }

    #[test]
    fn test_contiguous_run_linear() {
        let graph = make_test_graph();
        let lookup = build_edge_lookup(&graph, 3);
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // Two adjacent edges in the graph, expressed as their canonical kmers
        let kmers = vec![
            canonical_edge_kmer(&graph, e0, 3),
            canonical_edge_kmer(&graph, e1, 3),
        ];
        let runs = find_contiguous_runs(&kmers, &lookup, &graph);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].edges.len(), 2);
    }

    #[test]
    fn test_contiguous_run_gap() {
        let graph = make_test_graph();
        let lookup = build_edge_lookup(&graph, 3);
        let e0 = graph.edge_indices().next().unwrap();
        let e1 = graph.edge_indices().nth(1).unwrap();

        // A kmer that is not present in the lookup creates a gap
        let gap_kmer: u64 = 0xDEAD_BEEF;
        assert!(!lookup.contains_key(&gap_kmer));

        let kmers = vec![
            canonical_edge_kmer(&graph, e0, 3),
            gap_kmer,
            canonical_edge_kmer(&graph, e1, 3),
        ];
        let runs = find_contiguous_runs(&kmers, &lookup, &graph);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].edges.len(), 1);
        assert_eq!(runs[1].edges.len(), 1);
    }

    /// Construct a graph where two directional edges share the same canonical
    /// kmer (inverted-repeat collision), then verify that:
    /// (1) the lookup retains both candidates, and
    /// (2) a contiguous run picks the candidate adjacent to the previous edge.
    #[test]
    fn test_inverted_repeat_disambiguation() {
        // We need two edges whose reconstructed kmers are reverse complements
        // of each other but whose (src, tgt) node pairs are distinct. Using
        // k=3 for ease of reasoning:
        //
        //   Stem arm 1: n_a --kmer_X--> n_b
        //   Stem arm 2: n_c --kmer_rcX--> n_d
        //
        // And a bridging edge from n_b into n_c so that a "read" can walk
        // n_a -> n_b -> n_c -> n_d, traversing both arms as separate edges.
        //
        // Pick kmer X = AAC (0b000001). rc(AAC) = GTT (0b101111).
        // Source sub_kmer of X = AA (0b0000), target sub_kmer = AC (0b0001).
        // Source sub_kmer of rc(X)=GTT = GT (0b1011), target sub_kmer = TT (0b1111).
        //
        // These four sub_kmers are all distinct, so the node uniqueness
        // invariant holds.
        let mut graph = StableDiGraph::new();
        let n_a = graph.add_node(DBNode {
            sub_kmer: 0b0000, // AA
            is_start: true,
            is_end: false,
            is_terminal: false,
        });
        let n_b = graph.add_node(DBNode {
            sub_kmer: 0b0001, // AC
            is_start: false,
            is_end: false,
            is_terminal: false,
        });
        let n_c = graph.add_node(DBNode {
            sub_kmer: 0b1011, // GT
            is_start: false,
            is_end: false,
            is_terminal: false,
        });
        let n_d = graph.add_node(DBNode {
            sub_kmer: 0b1111, // TT
            is_start: false,
            is_end: true,
            is_terminal: false,
        });

        // Edge e_x: AAC, n_a -> n_b
        let e_x = graph.add_edge(
            n_a,
            n_b,
            DBEdge {
                count: 10,
                coverage_ratio: 1.0,
            },
        );
        // Bridge edge: n_b -> n_c (kmer ACGT-like, not used for lookup in this
        // test — we just need graph connectivity for adjacency walking)
        let e_bridge = graph.add_edge(
            n_b,
            n_c,
            DBEdge {
                count: 10,
                coverage_ratio: 1.0,
            },
        );
        // Edge e_rcx: GTT, n_c -> n_d
        let e_rcx = graph.add_edge(
            n_c,
            n_d,
            DBEdge {
                count: 10,
                coverage_ratio: 1.0,
            },
        );

        let lookup = build_edge_lookup(&graph, 3);

        // AAC and GTT should share a canonical key with both edges as candidates
        let canonical = canonical_edge_kmer(&graph, e_x, 3);
        assert_eq!(canonical, canonical_edge_kmer(&graph, e_rcx, 3));
        let cands = lookup.get(&canonical).expect("collision key present");
        assert_eq!(
            cands.len(),
            2,
            "both edges must be retained under the canonical key"
        );
        assert!(cands.contains(&e_x));
        assert!(cands.contains(&e_rcx));

        // A read that walks n_a -> n_b -> n_c -> n_d produces kmers
        // [canonical(e_x), canonical(e_bridge), canonical(e_rcx)]. With
        // disambiguation by previous-edge adjacency, the run should pick e_x
        // at the start (no previous edge, but then e_bridge anchors the run
        // to n_b -> n_c, and e_rcx's source n_c is adjacent to e_bridge's
        // target, so e_rcx wins over e_x at the third position).
        let kmers = vec![
            canonical_edge_kmer(&graph, e_x, 3),
            canonical_edge_kmer(&graph, e_bridge, 3),
            canonical_edge_kmer(&graph, e_rcx, 3),
        ];
        let runs = find_contiguous_runs(&kmers, &lookup, &graph);

        // Either one run of length 3, or a fresh run starting from the
        // collision point. In either case, e_rcx (not e_x) must appear at
        // the third position and must follow e_bridge.
        let picked_third = runs
            .iter()
            .flat_map(|r| r.edges.iter())
            .nth(2)
            .copied()
            .expect("three edges must be resolved");
        assert_eq!(
            picked_third, e_rcx,
            "disambiguation must pick the edge adjacent to the previous run edge"
        );
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

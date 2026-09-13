// pcr/paths.rs — path finding, sequence extraction, dedup

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use bio::alignment::distance::simd::*;
use bio::io::fasta;
use log::debug;
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;
use sha2::{Digest, Sha256};
use smallvec::SmallVec;

use crate::kmer::FilteredKmerCounts;

use super::graph::{RepeatMarkers, compute_mean, compute_median, get_end_nodes, get_start_nodes};
use super::{
    AssemblyRecord, DBEdge, DBNode, PCRParams, WithheldDiagnosticLimits, WithheldDiagnosticUsage,
    WithheldMarkerCause, WithheldMarkerCoveredRun, WithheldMarkerIdentityKind,
    WithheldMarkerOccurrence, WithheldPathDiagnostics, WithheldPathRecord,
};

/// The maximum number of fasta records to return
const MAX_NUM_AMPLICONS: usize = 20;

/// One step along a found path: the node visited, and the edge that brought
/// us into it. The starting node of a path has `edge = None`; every other
/// step records the `EdgeIndex` of the edge that connects it to the previous
/// step. Carrying the edge alongside the node lets `generate_sequences_from_paths`
/// look up edge weights directly via `graph.edge_weight()` instead of calling
/// `graph.find_edge(parent, node)`, which is O(out_degree(parent)) per step
/// and was previously called three times per edge per path.
pub type PathStep = (NodeIndex, Option<EdgeIndex>);

pub(super) struct PathSearchResult {
    pub paths: Vec<Vec<PathStep>>,
    pub completed_candidate_path_count: usize,
    pub eligible_candidate_path_count: usize,
    pub withheld_candidate_path_count: usize,
    pub withheld_by_self_loop_node_count: usize,
    pub withheld_by_scc_node_count: usize,
    pub withheld_by_collision_edge_count: usize,
    pub repeat_tainted_out_of_range_path_count: usize,
    pub dfs_limit_reached: bool,
    pub path_limit_reached: bool,
    pub end_below_min_length: bool,
    pub max_length_reached: bool,
    pub withheld_path_diagnostics: Option<WithheldPathDiagnostics>,
}

pub(super) struct WithheldDiagnosticBudget {
    limits: WithheldDiagnosticLimits,
    usage: WithheldDiagnosticUsage,
}

impl WithheldDiagnosticBudget {
    pub(super) fn new(limits: WithheldDiagnosticLimits) -> Self {
        Self {
            limits,
            usage: WithheldDiagnosticUsage::default(),
        }
    }

    pub(super) fn empty_threshold_diagnostics(
        &self,
        node_visit_limit: usize,
    ) -> WithheldPathDiagnostics {
        WithheldPathDiagnostics {
            schema_version: 1,
            coordinate_system: "zero_based_half_open".to_string(),
            purpose: "diagnostic_only_unsupported_candidate_hypotheses".to_string(),
            marker_span_scope:
                "marker_covered_positions_not_complete_ambiguity_or_bridge_ready_intervals"
                    .to_string(),
            read_support_evaluation: "not_evaluated".to_string(),
            limits: self.limits,
            gene_usage_before_threshold: self.usage,
            gene_usage_after_threshold: self.usage,
            node_visit_limit,
            ..WithheldPathDiagnostics::default()
        }
    }
}

struct ThresholdDiagnosticCollector<'a> {
    budget: &'a mut WithheldDiagnosticBudget,
    payload: WithheldPathDiagnostics,
}

#[derive(Clone, Copy)]
struct WithheldPathCauseCounts {
    omitted_self_loop: usize,
    cyclic_scc: usize,
    retained_collision: usize,
}

impl WithheldPathCauseCounts {
    fn total(self) -> usize {
        self.omitted_self_loop + self.cyclic_scc + self.retained_collision
    }
}

impl<'a> ThresholdDiagnosticCollector<'a> {
    fn new(budget: &'a mut WithheldDiagnosticBudget, node_visit_limit: usize) -> Self {
        let payload = budget.empty_threshold_diagnostics(node_visit_limit);
        Self { budget, payload }
    }

    fn record_node_visit_skip(&mut self) {
        self.payload.node_visit_skips += 1;
    }

    fn observe_withheld_path(
        &mut self,
        graph: &StableDiGraph<DBNode, DBEdge>,
        path: &[PathStep],
        kmer_length: usize,
        repeat_markers: &RepeatMarkers,
        visit_counts: &AHashMap<NodeIndex, usize>,
        cause_counts: WithheldPathCauseCounts,
    ) {
        self.payload.observed_withheld_paths += 1;
        let sequence_length = path.len().saturating_add(kmer_length).saturating_sub(2);
        let marker_occurrence_count = cause_counts.total();
        let allocated = self.budget.limits.allocated_run_share_for_gene;
        let threshold = self.budget.limits.threshold;

        let drop_reason = if sequence_length > threshold.sequence_bases
            || sequence_length > allocated.sequence_bases
        {
            Some("oversize_sequence")
        } else if marker_occurrence_count > threshold.marker_occurrences
            || marker_occurrence_count > allocated.marker_occurrences
        {
            Some("oversize_marker_set")
        } else if self.payload.retained_paths >= threshold.paths
            || self.budget.usage.paths >= allocated.paths
        {
            Some("path_cap")
        } else if self
            .payload
            .retained_sequence_bases
            .saturating_add(sequence_length)
            > threshold.sequence_bases
            || self
                .budget
                .usage
                .sequence_bases
                .saturating_add(sequence_length)
                > allocated.sequence_bases
        {
            Some("sequence_base_cap")
        } else if self
            .payload
            .retained_marker_occurrences
            .saturating_add(marker_occurrence_count)
            > threshold.marker_occurrences
            || self
                .budget
                .usage
                .marker_occurrences
                .saturating_add(marker_occurrence_count)
                > allocated.marker_occurrences
        {
            Some("marker_cap")
        } else {
            None
        };

        if let Some(drop_reason) = drop_reason {
            self.payload.observed_not_retained_paths += 1;
            match drop_reason {
                "oversize_sequence" => self.payload.dropped_oversize_sequence += 1,
                "oversize_marker_set" => self.payload.dropped_oversize_marker_set += 1,
                "path_cap" => self.payload.dropped_path_cap += 1,
                "sequence_base_cap" => self.payload.dropped_sequence_base_cap += 1,
                "marker_cap" => self.payload.dropped_marker_cap += 1,
                _ => unreachable!(),
            }
            self.payload.retention_truncated = true;
            return;
        }

        let oriented_sequence = reconstruct_path_sequence(graph, path, kmer_length);
        let marker_occurrences = collect_marker_occurrences(
            graph,
            path,
            kmer_length,
            repeat_markers,
            &oriented_sequence,
        );
        debug_assert_eq!(marker_occurrences.len(), marker_occurrence_count);
        let marker_covered_runs = merge_marker_covered_runs(&marker_occurrences);
        let max_observed_node_visits = path
            .iter()
            .filter_map(|(node, _)| visit_counts.get(node).copied())
            .max()
            .unwrap_or(0);
        let sequence_sha256 = sha256_text(&oriented_sequence);

        self.payload.retained_paths += 1;
        self.payload.retained_sequence_bases += sequence_length;
        self.payload.retained_marker_occurrences += marker_occurrence_count;
        self.budget.usage.paths += 1;
        self.budget.usage.sequence_bases += sequence_length;
        self.budget.usage.marker_occurrences += marker_occurrence_count;
        self.payload.paths.push(WithheldPathRecord {
            oriented_sequence,
            sequence_sha256,
            sequence_length,
            node_visit_limit: self.payload.node_visit_limit,
            max_observed_node_visits,
            omitted_self_loop_marker_occurrences: cause_counts.omitted_self_loop,
            cyclic_scc_marker_occurrences: cause_counts.cyclic_scc,
            retained_collision_marker_occurrences: cause_counts.retained_collision,
            marker_occurrences,
            marker_covered_runs,
        });
    }

    fn finish(mut self) -> WithheldPathDiagnostics {
        self.payload.gene_usage_after_threshold = self.budget.usage;
        debug_assert_eq!(
            self.payload.observed_withheld_paths,
            self.payload.retained_paths + self.payload.observed_not_retained_paths
        );
        self.payload
    }
}

fn reconstruct_path_sequence(
    graph: &StableDiGraph<DBNode, DBEdge>,
    path: &[PathStep],
    kmer_length: usize,
) -> String {
    let mut sequence = String::with_capacity(path.len().saturating_add(kmer_length));
    for (path_position, (node, _)) in path.iter().enumerate() {
        if path_position == 0 {
            sequence.push_str(&crate::kmer::kmer_to_seq(
                &graph[*node].sub_kmer,
                &(kmer_length - 1),
            ));
        } else {
            sequence.push(crate::kmer::kmer_last_base(&graph[*node].sub_kmer));
        }
    }
    sequence
}

fn sha256_text(value: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(value.as_bytes());
    format!("{:x}", hasher.finalize())
}

fn collect_marker_occurrences(
    graph: &StableDiGraph<DBNode, DBEdge>,
    path: &[PathStep],
    kmer_length: usize,
    repeat_markers: &RepeatMarkers,
    oriented_sequence: &str,
) -> Vec<WithheldMarkerOccurrence> {
    let mut occurrences = Vec::new();
    for (path_position, (node, incoming_edge)) in path.iter().enumerate() {
        let node_identity = graph[*node].sub_kmer;
        let node_start = path_position;
        let node_end = path_position + kmer_length - 1;
        for cause in [
            repeat_markers
                .omitted_self_loop_sub_kmers
                .contains(&node_identity)
                .then_some(WithheldMarkerCause::PrePruneOmittedSelfLoopNode),
            repeat_markers
                .cyclic_scc_sub_kmers
                .contains(&node_identity)
                .then_some(WithheldMarkerCause::PrePruneCyclicSccNode),
        ]
        .into_iter()
        .flatten()
        {
            let identity_sequence = crate::kmer::kmer_to_seq(&node_identity, &(kmer_length - 1));
            debug_assert_eq!(
                oriented_sequence.get(node_start..node_end),
                Some(identity_sequence.as_str())
            );
            occurrences.push(WithheldMarkerOccurrence {
                cause,
                identity_kind: WithheldMarkerIdentityKind::OrientedNodeSequence,
                identity_sha256: sha256_text(&identity_sequence),
                identity_sequence,
                start: node_start,
                end: node_end,
            });
        }

        if let Some(edge) = incoming_edge
            && repeat_markers.retained_collision_edges.contains(edge)
        {
            let edge_start = path_position - 1;
            let edge_end = path_position + kmer_length - 1;
            let edge_identity = super::graph::reconstruct_edge_kmer(graph, *edge);
            let identity_sequence = crate::kmer::kmer_to_seq(&edge_identity, &kmer_length);
            debug_assert_eq!(
                oriented_sequence.get(edge_start..edge_end),
                Some(identity_sequence.as_str())
            );
            occurrences.push(WithheldMarkerOccurrence {
                cause: WithheldMarkerCause::RetainedCollisionEdge,
                identity_kind: WithheldMarkerIdentityKind::OrientedEdgeSequence,
                identity_sha256: sha256_text(&identity_sequence),
                identity_sequence,
                start: edge_start,
                end: edge_end,
            });
        }
    }
    occurrences
}

fn merge_marker_covered_runs(
    marker_occurrences: &[WithheldMarkerOccurrence],
) -> Vec<WithheldMarkerCoveredRun> {
    let mut spans: Vec<(usize, usize)> = marker_occurrences
        .iter()
        .map(|occurrence| (occurrence.start, occurrence.end))
        .collect();
    spans.sort_unstable();
    let mut runs: Vec<WithheldMarkerCoveredRun> = Vec::new();
    for (start, end) in spans {
        if let Some(previous) = runs.last_mut()
            && start <= previous.end
        {
            previous.end = previous.end.max(end);
        } else {
            runs.push(WithheldMarkerCoveredRun { start, end });
        }
    }
    runs
}

/// Children-of-a-node, sorted by score, used as a DFS frame in `child_stack`.
/// In a de Bruijn graph each node's outgoing edge set is bounded by ≤4 (one
/// per possible 2-bit suffix base), so the inline capacity exactly matches
/// the worst case and the SmallVec never spills to the heap. This replaces
/// a per-call `Vec` heap allocation that fired at every visited node during
/// DFS — up to `max_dfs_states = 100K` times per gene.
type ChildFrame = SmallVec<[(NodeIndex, EdgeIndex, f64); 4]>;

/// Get outgoing edges sorted by score (ascending so pop gives highest first).
/// Each entry carries `(target, edge_id, score)` so the DFS can record the
/// edge directly into the path without re-finding it later.
fn sorted_children(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    edge_preferences: Option<&AHashMap<EdgeIndex, f64>>,
) -> ChildFrame {
    let mut outgoing: ChildFrame = graph
        .edges_directed(node, Direction::Outgoing)
        .map(|e| {
            let base_score = e.weight().count as f64;
            let pref = edge_preferences
                .and_then(|prefs| prefs.get(&e.id()))
                .copied()
                .unwrap_or(1.0);
            (e.target(), e.id(), base_score * pref)
        })
        .collect();
    // total_cmp gives a total order on f64 (NaN sorts to the end) instead of
    // partial_cmp's Option<Ordering>. Inputs here are u32 → f64 multiplied by
    // a finite preference, so NaN is impossible — total_cmp simply removes
    // the unwrap_or(Equal) papering over an "impossible" case.
    outgoing.sort_by(|a, b| a.2.total_cmp(&b.2));
    outgoing
}

/// Find paths from start nodes to end nodes using coverage-weighted DFS.
/// At each branch point, outgoing edges are explored in descending order
/// of edge count, so the highest-coverage paths are found first.
/// When `edge_preferences` are provided (from bubble resolution), they
/// boost the ordering of read-supported edges at branch points.
///
/// Each returned path is a list of `PathStep`s. The first step's
/// `edge` is `None` (the start node has no incoming edge in the path);
/// every subsequent step's `edge` is `Some(EdgeIndex)` of the edge that
/// connects the previous node to this one. Recording the edge alongside
/// the node lets downstream code skip three `find_edge` calls per edge
/// per path in `generate_sequences_from_paths`.
#[cfg(test)]
pub fn get_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
    edge_preferences: Option<&AHashMap<EdgeIndex, f64>>,
) -> Vec<Vec<PathStep>> {
    search_assembly_paths(
        graph,
        kmer_counts,
        params,
        edge_preferences,
        &RepeatMarkers::default(),
        None,
    )
    .paths
}

pub(super) fn search_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
    edge_preferences: Option<&AHashMap<EdgeIndex, f64>>,
    repeat_markers: &RepeatMarkers,
    diagnostic_budget: Option<&mut WithheldDiagnosticBudget>,
) -> PathSearchResult {
    // A path of N nodes produces a sequence of (k-1) + (N-1) = N+k-2 bases
    // (first node contributes k-1 bases via its sub_kmer, each subsequent
    // node extends by one base). Inverting: N = seq_length - k + 2.
    //
    // min_path_nodes is used only as a pre-filter during DFS; the exact
    // lower bound is re-checked in `generate_sequences_from_paths` against
    // the actual sequence length.
    //
    // max_path_nodes is the authoritative upper bound — the DFS caps paths
    // at this value and there is no secondary check downstream, so an
    // off-by-one here directly rejects valid amplicons at the user's
    // declared max-length.
    let k = kmer_counts.get_k();
    let node_length_offset = k.saturating_sub(2);
    let min_path_nodes = params.min_length.saturating_sub(node_length_offset).max(1);
    let max_path_nodes = params.max_length.saturating_sub(node_length_offset).max(1);

    let end_nodes: AHashSet<NodeIndex> = get_end_nodes(graph).into_iter().collect();
    let mut result = PathSearchResult {
        paths: Vec::new(),
        completed_candidate_path_count: 0,
        eligible_candidate_path_count: 0,
        withheld_candidate_path_count: 0,
        withheld_by_self_loop_node_count: 0,
        withheld_by_scc_node_count: 0,
        withheld_by_collision_edge_count: 0,
        repeat_tainted_out_of_range_path_count: 0,
        dfs_limit_reached: false,
        path_limit_reached: false,
        end_below_min_length: false,
        max_length_reached: false,
        withheld_path_diagnostics: None,
    };
    let mut diagnostic_collector = diagnostic_budget
        .map(|budget| ThresholdDiagnosticCollector::new(budget, params.max_node_visits));

    for start in get_start_nodes(graph) {
        let mut paths_from_start = 0;
        let mut states_explored: usize = 0;

        // Stack-based DFS with push/pop backtracking instead of cloning.
        // `path` and `visit_counts` are maintained incrementally.
        // `child_stack` holds remaining children to explore at each depth.
        // The starting step has edge = None; all subsequent steps record
        // the edge that connects the previous node to this one.
        let mut path: Vec<PathStep> = vec![(start, None)];
        let mut visit_counts: AHashMap<NodeIndex, usize> = AHashMap::new();
        visit_counts.insert(start, 1);
        let mut omitted_self_loop_nodes_in_path = usize::from(
            repeat_markers
                .omitted_self_loop_sub_kmers
                .contains(&graph[start].sub_kmer),
        );
        let mut cyclic_scc_nodes_in_path = usize::from(
            repeat_markers
                .cyclic_scc_sub_kmers
                .contains(&graph[start].sub_kmer),
        );
        let mut retained_collision_edges_in_path = 0_usize;

        // Compute sorted children for the start node
        let children = sorted_children(graph, start, edge_preferences);
        let mut child_stack: Vec<ChildFrame> = vec![children];

        loop {
            if paths_from_start >= params.max_paths_per_pair {
                result.path_limit_reached = true;
                break;
            }
            if states_explored >= params.max_dfs_states {
                result.dfs_limit_reached = true;
                break;
            }

            let depth = child_stack.len() - 1;

            if let Some((neighbor, edge_id, _score)) = child_stack[depth].pop() {
                states_explored += 1;

                let current_visits = visit_counts.get(&neighbor).copied().unwrap_or(0);
                if current_visits >= params.max_node_visits {
                    if let Some(collector) = diagnostic_collector.as_mut() {
                        collector.record_node_visit_skip();
                    }
                    continue;
                }

                // Push neighbor onto path with the edge that connects to it
                path.push((neighbor, Some(edge_id)));
                *visit_counts.entry(neighbor).or_insert(0) += 1;
                let neighbor_has_omitted_self_loop = repeat_markers
                    .omitted_self_loop_sub_kmers
                    .contains(&graph[neighbor].sub_kmer);
                let neighbor_is_in_cyclic_scc = repeat_markers
                    .cyclic_scc_sub_kmers
                    .contains(&graph[neighbor].sub_kmer);
                let edge_is_retained_collision =
                    repeat_markers.retained_collision_edges.contains(&edge_id);
                omitted_self_loop_nodes_in_path += usize::from(neighbor_has_omitted_self_loop);
                cyclic_scc_nodes_in_path += usize::from(neighbor_is_in_cyclic_scc);
                retained_collision_edges_in_path += usize::from(edge_is_retained_collision);

                let path_len = path.len();

                // Check if we reached an end node with valid length
                if end_nodes.contains(&neighbor) {
                    if path_len >= min_path_nodes && path_len <= max_path_nodes {
                        result.completed_candidate_path_count += 1;
                        let withheld_by_self_loop = omitted_self_loop_nodes_in_path > 0;
                        let withheld_by_scc = cyclic_scc_nodes_in_path > 0;
                        let withheld_by_collision = retained_collision_edges_in_path > 0;
                        if withheld_by_self_loop || withheld_by_scc || withheld_by_collision {
                            result.withheld_candidate_path_count += 1;
                            result.withheld_by_self_loop_node_count +=
                                usize::from(withheld_by_self_loop);
                            result.withheld_by_scc_node_count += usize::from(withheld_by_scc);
                            result.withheld_by_collision_edge_count +=
                                usize::from(withheld_by_collision);
                            if let Some(collector) = diagnostic_collector.as_mut() {
                                collector.observe_withheld_path(
                                    graph,
                                    &path,
                                    k,
                                    repeat_markers,
                                    &visit_counts,
                                    WithheldPathCauseCounts {
                                        omitted_self_loop: omitted_self_loop_nodes_in_path,
                                        cyclic_scc: cyclic_scc_nodes_in_path,
                                        retained_collision: retained_collision_edges_in_path,
                                    },
                                );
                            }
                        } else {
                            result.paths.push(path.clone());
                            result.eligible_candidate_path_count += 1;
                            paths_from_start += 1;
                        }
                        omitted_self_loop_nodes_in_path -=
                            usize::from(neighbor_has_omitted_self_loop);
                        cyclic_scc_nodes_in_path -= usize::from(neighbor_is_in_cyclic_scc);
                        retained_collision_edges_in_path -= usize::from(edge_is_retained_collision);
                        *visit_counts.get_mut(&neighbor).unwrap() -= 1;
                        path.pop();
                        continue;
                    }
                    if path_len < min_path_nodes {
                        result.end_below_min_length = true;
                    } else {
                        result.max_length_reached = true;
                    }
                    if omitted_self_loop_nodes_in_path > 0
                        || cyclic_scc_nodes_in_path > 0
                        || retained_collision_edges_in_path > 0
                    {
                        result.repeat_tainted_out_of_range_path_count += 1;
                    }
                }

                // Don't extend past max length
                if path_len >= max_path_nodes {
                    result.max_length_reached = true;
                    omitted_self_loop_nodes_in_path -= usize::from(neighbor_has_omitted_self_loop);
                    cyclic_scc_nodes_in_path -= usize::from(neighbor_is_in_cyclic_scc);
                    retained_collision_edges_in_path -= usize::from(edge_is_retained_collision);
                    *visit_counts.get_mut(&neighbor).unwrap() -= 1;
                    path.pop();
                    continue;
                }

                // Explore deeper: push children frame
                let children = sorted_children(graph, neighbor, edge_preferences);
                child_stack.push(children);
            } else {
                // No more children at this depth, backtrack
                child_stack.pop();
                if child_stack.is_empty() {
                    break;
                }
                let (backtrack_node, incoming_edge) =
                    path.pop().expect("BUG: path empty during DFS backtrack");
                omitted_self_loop_nodes_in_path -= usize::from(
                    repeat_markers
                        .omitted_self_loop_sub_kmers
                        .contains(&graph[backtrack_node].sub_kmer),
                );
                cyclic_scc_nodes_in_path -= usize::from(
                    repeat_markers
                        .cyclic_scc_sub_kmers
                        .contains(&graph[backtrack_node].sub_kmer),
                );
                retained_collision_edges_in_path -=
                    usize::from(incoming_edge.is_some_and(|edge_id| {
                        repeat_markers.retained_collision_edges.contains(&edge_id)
                    }));
                *visit_counts
                    .get_mut(&backtrack_node)
                    .expect("BUG: backtrack node missing from visit_counts") -= 1;
            }
        }
    }

    result.withheld_path_diagnostics =
        diagnostic_collector.map(ThresholdDiagnosticCollector::finish);

    result
}

/// Extract sequences from graph paths, producing FASTA assembly records.
/// Returns the generated records and the updated amplicon index.
pub(super) fn generate_sequences_from_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    all_paths: Vec<Vec<PathStep>>,
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    mut amplicon_index: usize,
    threading: Option<&super::threading::ThreadingAnnotations>,
) -> Result<(Vec<AssemblyRecord>, usize)> {
    let mut assembly_records: Vec<AssemblyRecord> = Vec::new();

    for path in all_paths.into_iter() {
        let mut sequence = String::new();
        let mut edge_counts: Vec<u64> = Vec::new(); // u64 for compute_mean/median compatibility
        // Edges traversed by this path, in order. Recorded once during the
        // sequence-extraction loop and reused for max_coverage_ratio and
        // threading-support metrics below — replaces three find_edge passes.
        let mut path_edges: Vec<EdgeIndex> = Vec::with_capacity(path.len().saturating_sub(1));
        for &(node, edge_opt) in path.iter() {
            let node_data = graph
                .node_weight(node)
                .context("Node not found in graph during sequence generation")?;
            if sequence.is_empty() {
                sequence =
                    crate::kmer::kmer_to_seq(&node_data.sub_kmer, &(kmer_counts.get_k() - 1));
                debug_assert!(
                    edge_opt.is_none(),
                    "first PathStep must have edge = None (this is the start node)"
                );
            } else {
                sequence.push(crate::kmer::kmer_last_base(&node_data.sub_kmer));
                let edge_idx =
                    edge_opt.context("Non-start PathStep is missing its edge index (DFS bug)")?;
                let edge_data = graph
                    .edge_weight(edge_idx)
                    .context("Edge weight not found")?;
                edge_counts.push(edge_data.count as u64);
                path_edges.push(edge_idx);
            }
        }

        if sequence.len() < params.min_length || sequence.len() > params.max_length {
            debug!(
                "  Path length is {} bp, outside requested range {}-{}. Skipping.",
                sequence.len(),
                params.min_length,
                params.max_length
            );
            continue;
        }

        if edge_counts.is_empty() {
            debug!("  Path has no edges (single node). Skipping.");
            continue;
        }
        let count_mean = compute_mean(&edge_counts);
        let count_median = compute_median(&edge_counts);
        let count_min = edge_counts
            .iter()
            .min()
            .context("No edge counts found for path")?;
        let count_max = edge_counts
            .iter()
            .max()
            .context("No edge counts found for path")?;

        // Compute coverage consistency (coefficient of variation)
        let coverage_cv = if count_mean > 0.0 {
            let variance = edge_counts
                .iter()
                .map(|&c| {
                    let diff = c as f64 - count_mean;
                    diff * diff
                })
                .sum::<f64>()
                / edge_counts.len() as f64;
            variance.sqrt() / count_mean
        } else {
            0.0
        };

        // Compute max coverage ratio from edge annotations.
        // Reuses the path_edges list captured during sequence extraction.
        let max_coverage_ratio = path_edges
            .iter()
            .map(|&edge_idx| {
                graph
                    .edge_weight(edge_idx)
                    .map_or(0.0, |e| e.coverage_ratio)
            })
            .fold(0.0_f64, f64::max);

        // Compute read-support metrics from threading annotations.
        // Reuses path_edges; threading.edge_support is keyed by EdgeIndex
        // so a direct hash lookup replaces the previous find_edge probe.
        let (zero_support_edges, median_unambiguous_support, edge_support_fraction) =
            if let Some(ann) = threading {
                let mut total_edges = 0u32;
                let mut supported_edges = 0u32;
                let mut zero_count = 0u32;
                let mut unambiguous_counts: Vec<u64> = Vec::new();

                for &edge_idx in &path_edges {
                    total_edges += 1;
                    match ann.edge_support.get(&edge_idx) {
                        Some(s) if s.read_support_total > 0 => {
                            supported_edges += 1;
                            unambiguous_counts.push(s.read_support_unambiguous as u64);
                        }
                        _ => {
                            zero_count += 1;
                            unambiguous_counts.push(0);
                        }
                    }
                }

                let frac = if total_edges > 0 {
                    supported_edges as f64 / total_edges as f64
                } else {
                    0.0
                };
                let median_unamb = if unambiguous_counts.is_empty() {
                    0.0
                } else {
                    compute_median(&unambiguous_counts)
                };

                (Some(zero_count), Some(median_unamb), Some(frac))
            } else {
                (None, None, None)
            };

        let score = super::PathScore {
            kmer_min_count: *count_min as u32,
            kmer_median_count: count_median,
            coverage_cv,
            max_coverage_ratio,
            zero_support_edges,
            median_unambiguous_support,
            edge_support_fraction,
        };

        let id = format!("{}_{}_{}", sample_name, params.gene_name, amplicon_index);
        let desc = format!(
            "sample={} gene={} product={} length={} kmer_count_mean={:.2} kmer_count_median={} kmer_count_min={} kmer_count_max={} score={:.2}",
            sample_name,
            params.gene_name,
            amplicon_index,
            sequence.len(),
            count_mean,
            count_median,
            count_min,
            count_max,
            score.composite()
        );

        amplicon_index += 1;

        debug!(">{} {}", id, desc);
        let record = fasta::Record::with_attrs(&id, Some(&desc), sequence.as_bytes());
        assembly_records.push(AssemblyRecord {
            fasta_record: record,
            score,
        });
    }

    Ok((assembly_records, amplicon_index))
}

/// Sort assembly records by descending composite score, deduplicate near-identical
/// sequences, and truncate to the maximum number of amplicons.
pub(super) fn sort_and_deduplicate(
    assembly_records: Vec<AssemblyRecord>,
    params: &PCRParams,
) -> Vec<fasta::Record> {
    let mut sorted = assembly_records;
    sorted.sort_by(|a, b| {
        // total_cmp is a total order on f64; composite() multiplies finite
        // u32/u64-derived values, so NaN is impossible here and total_cmp
        // is strictly better than partial_cmp().unwrap_or(Equal).
        b.score
            .composite()
            .total_cmp(&a.score.composite())
            // Byte-level sequence comparison as a deterministic tiebreaker
            // when two records tie on composite score. Not biologically
            // meaningful — its only purpose is to make dedup order stable
            // across runs (otherwise iteration order would depend on
            // HashMap seeds deeper in path enumeration).
            .then_with(|| a.fasta_record.seq().cmp(b.fasta_record.seq()))
    });

    let records: Vec<fasta::Record> = sorted.into_iter().map(|ar| ar.fasta_record).collect();
    let num_records_all = records.len();

    // Greedy clustering: keep each record only if it is not within
    // dedup_edit_threshold edits of any previously kept record.
    // Distances are computed on-the-fly to avoid O(N^2) memory.
    let mut kept: Vec<fasta::Record> = Vec::new();
    for record in records {
        let is_duplicate = kept.iter().any(|kept_record| {
            bounded_levenshtein(record.seq(), kept_record.seq(), params.dedup_edit_threshold)
                .is_some()
        });
        if !is_duplicate {
            kept.push(record);
        }
    }
    let mut records = kept;

    if num_records_all == records.len() {
        gene_info!(
            params.gene_name,
            "{} PCR products were generated and retained.",
            num_records_all
        );
    } else {
        gene_info!(
            params.gene_name,
            "{} PCR products were generated and {} were retained ({} removed as near-duplicates within edit distance {}).",
            num_records_all,
            records.len(),
            num_records_all - records.len(),
            params.dedup_edit_threshold
        );
    }

    if records.len() > MAX_NUM_AMPLICONS {
        gene_warn!(
            params.gene_name,
            "There are {} PCR products. This exceeds the maximum of {}. Retaining only the first {} records.",
            records.len(),
            MAX_NUM_AMPLICONS,
            MAX_NUM_AMPLICONS
        );
        records.truncate(MAX_NUM_AMPLICONS);
    }

    records
}

#[cfg(test)]
mod tests {
    use super::super::{
        DBEdge, DBNode, DEFAULT_DEDUP_EDIT_THRESHOLD, DEFAULT_HIGH_COVERAGE_RATIO,
        DEFAULT_MAX_DFS_STATES, DEFAULT_MAX_NODE_VISITS, DEFAULT_MAX_NUM_PRIMER_KMERS,
        DEFAULT_MAX_PATHS_PER_PAIR, DEFAULT_TIP_COVERAGE_FRACTION,
        WITHHELD_DIAGNOSTIC_RUN_BASE_CAP, WITHHELD_DIAGNOSTIC_RUN_MARKER_CAP,
        WITHHELD_DIAGNOSTIC_RUN_PATH_CAP, withheld_diagnostic_limits_for_gene,
    };
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

    fn encode(sequence: &str) -> u64 {
        crate::kmer::encoding::seq_to_kmer(sequence).unwrap()
    }

    fn sequence_graph(sequence: &str, kmer_length: usize) -> StableDiGraph<DBNode, DBEdge> {
        let mut graph = StableDiGraph::new();
        let mut nodes = Vec::new();
        for position in 0..=sequence.len() - (kmer_length - 1) {
            let node_sequence = &sequence[position..position + kmer_length - 1];
            let sub_kmer = crate::kmer::encoding::seq_to_kmer(node_sequence).unwrap();
            nodes.push(graph.add_node(mk_node(
                sub_kmer,
                position == 0,
                position == sequence.len() - (kmer_length - 1),
            )));
        }
        for window in nodes.windows(2) {
            graph.add_edge(window[0], window[1], mk_edge(10));
        }
        graph
    }

    fn diagnostic_limits(
        paths: usize,
        sequence_bases: usize,
        marker_occurrences: usize,
    ) -> WithheldDiagnosticLimits {
        let cap = super::super::WithheldDiagnosticCapSet {
            paths,
            sequence_bases,
            marker_occurrences,
        };
        WithheldDiagnosticLimits {
            threshold: cap,
            gene: cap,
            run: cap,
            allocated_run_share_for_gene: cap,
            unused_share_is_not_redistributed: true,
        }
    }

    fn path_visit_counts(path: &[PathStep]) -> AHashMap<NodeIndex, usize> {
        let mut counts = AHashMap::new();
        for (node, _) in path {
            *counts.entry(*node).or_insert(0) += 1;
        }
        counts
    }

    fn cause_counts(
        omitted_self_loop: usize,
        cyclic_scc: usize,
        retained_collision: usize,
    ) -> WithheldPathCauseCounts {
        WithheldPathCauseCounts {
            omitted_self_loop,
            cyclic_scc,
            retained_collision,
        }
    }

    fn test_params(min_length: usize, max_length: usize) -> PCRParams {
        PCRParams {
            forward_seq: "ACGT".to_string(),
            reverse_seq: "TGCA".to_string(),
            min_length,
            max_length,
            gene_name: "test".to_string(),
            gene: None,
            region: None,
            index: None,
            compartment: None,
            gene_type: None,
            copy_number: None,
            deprecated: false,
            deprecated_by: None,
            deprecated_reason: None,
            min_count: 2,
            mismatches: 0,
            trim: 0,
            expected_length: None,
            citation: String::new(),
            notes: String::new(),
            dedup_edit_threshold: DEFAULT_DEDUP_EDIT_THRESHOLD,
            source: "test".to_string(),
            max_dfs_states: DEFAULT_MAX_DFS_STATES,
            max_paths_per_pair: DEFAULT_MAX_PATHS_PER_PAIR,
            max_node_visits: DEFAULT_MAX_NODE_VISITS,
            max_primer_kmers: DEFAULT_MAX_NUM_PRIMER_KMERS,
            high_coverage_ratio: DEFAULT_HIGH_COVERAGE_RATIO,
            tip_coverage_fraction: DEFAULT_TIP_COVERAGE_FRACTION,
            diagnose_withheld_paths: false,
        }
    }

    /// Linear graph: start -> a -> b -> end. Should find exactly one path.
    #[test]
    fn test_linear_path() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let e = graph.add_node(mk_node(3, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, e, mk_edge(10));

        // k=3, path of 4 nodes = 4 + 3 - 2 = 5 bases
        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert_eq!(paths.len(), 1);
        let nodes: Vec<NodeIndex> = paths[0].iter().map(|&(n, _)| n).collect();
        assert_eq!(nodes, vec![s, a, b, e]);
        // First step's edge is None (start node), the rest are Some.
        assert!(paths[0][0].1.is_none());
        assert!(paths[0][1..].iter().all(|&(_, e)| e.is_some()));
    }

    #[test]
    fn withheld_diagnostics_preserve_search_and_record_stable_marker_spans() {
        let graph = sequence_graph("ACGTA", 3);
        let path_nodes: Vec<NodeIndex> = graph.node_indices().collect();
        let collision_edge = graph.find_edge(path_nodes[0], path_nodes[1]).unwrap();
        let marked_sub_kmer = graph[path_nodes[1]].sub_kmer;
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]),
            cyclic_scc_sub_kmers: AHashSet::from_iter([marked_sub_kmer]),
            retained_collision_edges: AHashSet::from_iter([collision_edge]),
        };
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("ACGTA").unwrap();
        let filtered = counts.filtered_view(1);
        let params = test_params(5, 5);
        let without_diagnostics =
            search_assembly_paths(&graph, &filtered, &params, None, &markers, None);
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(4, 100, 20));
        let with_diagnostics = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &markers,
            Some(&mut budget),
        );

        assert_eq!(without_diagnostics.paths, with_diagnostics.paths);
        assert_eq!(without_diagnostics.completed_candidate_path_count, 1);
        assert_eq!(with_diagnostics.completed_candidate_path_count, 1);
        assert_eq!(without_diagnostics.withheld_candidate_path_count, 1);
        assert_eq!(with_diagnostics.withheld_candidate_path_count, 1);
        assert!(without_diagnostics.withheld_path_diagnostics.is_none());
        let diagnostics = with_diagnostics.withheld_path_diagnostics.unwrap();
        assert_eq!(diagnostics.observed_withheld_paths, 1);
        assert_eq!(diagnostics.retained_paths, 1);
        assert_eq!(diagnostics.observed_not_retained_paths, 0);
        let record = &diagnostics.paths[0];
        assert_eq!(record.oriented_sequence, "ACGTA");
        assert_eq!(record.sequence_sha256, sha256_text("ACGTA"));
        assert_eq!(record.sequence_length, 5);
        assert_eq!(record.omitted_self_loop_marker_occurrences, 1);
        assert_eq!(record.cyclic_scc_marker_occurrences, 1);
        assert_eq!(record.retained_collision_marker_occurrences, 1);
        assert_eq!(record.marker_occurrences.len(), 3);
        assert_eq!(record.marker_occurrences[0].identity_sequence, "CG");
        assert_eq!(
            (
                record.marker_occurrences[0].start,
                record.marker_occurrences[0].end
            ),
            (1, 3)
        );
        assert_eq!(record.marker_occurrences[1].identity_sequence, "CG");
        assert_eq!(
            (
                record.marker_occurrences[1].start,
                record.marker_occurrences[1].end
            ),
            (1, 3)
        );
        assert_eq!(record.marker_occurrences[2].identity_sequence, "ACG");
        assert_eq!(
            (
                record.marker_occurrences[2].start,
                record.marker_occurrences[2].end
            ),
            (0, 3)
        );
        assert_eq!(
            record.marker_covered_runs,
            vec![WithheldMarkerCoveredRun { start: 0, end: 3 }]
        );
        assert_eq!(
            diagnostics.marker_span_scope,
            "marker_covered_positions_not_complete_ambiguity_or_bridge_ready_intervals"
        );
        assert_eq!(diagnostics.read_support_evaluation, "not_evaluated");
    }

    #[test]
    fn withheld_diagnostics_preserve_clean_path_order_and_quota() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(encode("AC"), true, false));
        let repeat_first = graph.add_node(mk_node(encode("CG"), false, false));
        let repeat_second = graph.add_node(mk_node(encode("GT"), false, false));
        let clean_first = graph.add_node(mk_node(encode("CT"), false, false));
        let clean_second = graph.add_node(mk_node(encode("TT"), false, false));
        let end = graph.add_node(mk_node(encode("TA"), false, true));
        graph.add_edge(start, repeat_first, mk_edge(20));
        graph.add_edge(repeat_first, repeat_second, mk_edge(20));
        graph.add_edge(repeat_second, end, mk_edge(20));
        graph.add_edge(start, clean_first, mk_edge(10));
        graph.add_edge(clean_first, clean_second, mk_edge(10));
        graph.add_edge(clean_second, end, mk_edge(10));
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::new(),
            cyclic_scc_sub_kmers: AHashSet::from_iter([graph[repeat_first].sub_kmer]),
            retained_collision_edges: AHashSet::new(),
        };
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("ACGTA").unwrap();
        counts.ingest_seq("ACTTA").unwrap();
        let filtered = counts.filtered_view(1);
        let mut params = test_params(5, 5);
        params.max_paths_per_pair = 1;
        let without_diagnostics =
            search_assembly_paths(&graph, &filtered, &params, None, &markers, None);
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(8, 100, 20));
        let with_diagnostics = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &markers,
            Some(&mut budget),
        );

        assert_eq!(without_diagnostics.paths, with_diagnostics.paths);
        assert_eq!(without_diagnostics.completed_candidate_path_count, 2);
        assert_eq!(with_diagnostics.completed_candidate_path_count, 2);
        assert_eq!(without_diagnostics.eligible_candidate_path_count, 1);
        assert_eq!(with_diagnostics.eligible_candidate_path_count, 1);
        assert_eq!(without_diagnostics.withheld_candidate_path_count, 1);
        assert_eq!(with_diagnostics.withheld_candidate_path_count, 1);
        assert!(without_diagnostics.path_limit_reached);
        assert!(with_diagnostics.path_limit_reached);
        let clean_sequence = reconstruct_path_sequence(&graph, &with_diagnostics.paths[0], 3);
        assert_eq!(clean_sequence, "ACTTA");
    }

    #[test]
    fn pre_prune_scc_marker_identity_survives_cycle_pruning() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(encode("CA"), true, false));
        let repeat = graph.add_node(mk_node(encode("AA"), false, false));
        let end = graph.add_node(mk_node(encode("AT"), false, true));
        graph.add_edge(start, repeat, mk_edge(10));
        let cycle = graph.add_edge(repeat, repeat, mk_edge(10));
        graph.add_edge(repeat, end, mk_edge(10));
        let markers = super::super::graph::repeat_markers_before_pruning(
            &graph,
            super::super::graph::ExtensionRepeatMarkers::default(),
        );
        graph.remove_edge(cycle);
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("CAAT").unwrap();
        let filtered = counts.filtered_view(1);
        let params = test_params(4, 4);
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(4, 100, 20));
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &markers,
            Some(&mut budget),
        );

        let record = &result.withheld_path_diagnostics.unwrap().paths[0];
        assert_eq!(record.oriented_sequence, "CAAT");
        assert_eq!(record.cyclic_scc_marker_occurrences, 1);
        assert_eq!(record.marker_occurrences[0].identity_sequence, "AA");
        assert_eq!(
            (
                record.marker_occurrences[0].start,
                record.marker_occurrences[0].end
            ),
            (1, 3)
        );
    }

    #[test]
    fn withheld_diagnostic_sequences_preserve_graph_orientation() {
        let forward_graph = sequence_graph("ACGTA", 3);
        let reverse_graph = sequence_graph("TACGT", 3);
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("ACGTA").unwrap();
        counts.ingest_seq("TACGT").unwrap();
        let filtered = counts.filtered_view(1);
        let params = test_params(5, 5);

        let capture = |graph: &StableDiGraph<DBNode, DBEdge>| {
            let marked = graph[graph.node_indices().nth(1).unwrap()].sub_kmer;
            let markers = RepeatMarkers {
                omitted_self_loop_sub_kmers: AHashSet::from_iter([marked]),
                cyclic_scc_sub_kmers: AHashSet::new(),
                retained_collision_edges: AHashSet::new(),
            };
            let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(4, 100, 20));
            search_assembly_paths(graph, &filtered, &params, None, &markers, Some(&mut budget))
                .withheld_path_diagnostics
                .unwrap()
                .paths
                .remove(0)
        };

        let forward = capture(&forward_graph);
        let reverse = capture(&reverse_graph);
        assert_eq!(forward.oriented_sequence, "ACGTA");
        assert_eq!(reverse.oriented_sequence, "TACGT");
        assert_eq!(forward.marker_occurrences[0].identity_sequence, "CG");
        assert_eq!(reverse.marker_occurrences[0].identity_sequence, "AC");
        assert_eq!(
            bio::alphabets::dna::revcomp(forward.oriented_sequence.as_bytes()),
            reverse.oriented_sequence.as_bytes()
        );
    }

    #[test]
    fn repeated_node_marker_occurrences_keep_each_oriented_coordinate() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(encode("CA"), true, false));
        let repeat = graph.add_node(mk_node(encode("AA"), false, false));
        let end = graph.add_node(mk_node(encode("AT"), false, true));
        graph.add_edge(start, repeat, mk_edge(10));
        graph.add_edge(repeat, repeat, mk_edge(20));
        graph.add_edge(repeat, end, mk_edge(10));
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::new(),
            cyclic_scc_sub_kmers: AHashSet::from_iter([graph[repeat].sub_kmer]),
            retained_collision_edges: AHashSet::new(),
        };
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("CAAAT").unwrap();
        let filtered = counts.filtered_view(1);
        let mut params = test_params(4, 5);
        params.max_node_visits = 2;
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(8, 100, 20));
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &markers,
            Some(&mut budget),
        );

        let diagnostics = result.withheld_path_diagnostics.unwrap();
        let repeated = diagnostics
            .paths
            .iter()
            .find(|record| record.oriented_sequence == "CAAAT")
            .unwrap();
        let spans: Vec<(usize, usize)> = repeated
            .marker_occurrences
            .iter()
            .map(|occurrence| (occurrence.start, occurrence.end))
            .collect();
        assert_eq!(spans, vec![(1, 3), (2, 4)]);
        assert_eq!(repeated.max_observed_node_visits, 2);
        assert_eq!(repeated.node_visit_limit, 2);
        assert!(diagnostics.node_visit_skips > 0);
    }

    #[test]
    fn diagnostic_drop_reasons_are_exclusive_and_preserve_observed_identity() {
        let graph = sequence_graph("ACGTA", 3);
        let marked_sub_kmer = graph[graph.node_indices().nth(1).unwrap()].sub_kmer;
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]),
            cyclic_scc_sub_kmers: AHashSet::new(),
            retained_collision_edges: AHashSet::new(),
        };
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("ACGTA").unwrap();
        let filtered = counts.filtered_view(1);
        let params = test_params(5, 5);
        let mut zero_path_budget = WithheldDiagnosticBudget::new(diagnostic_limits(0, 100, 20));
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &markers,
            Some(&mut zero_path_budget),
        );
        let diagnostics = result.withheld_path_diagnostics.unwrap();

        assert_eq!(diagnostics.observed_withheld_paths, 1);
        assert_eq!(diagnostics.retained_paths, 0);
        assert_eq!(diagnostics.observed_not_retained_paths, 1);
        assert_eq!(diagnostics.dropped_path_cap, 1);
        assert_eq!(diagnostics.dropped_oversize_sequence, 0);
        assert_eq!(diagnostics.dropped_oversize_marker_set, 0);
        assert_eq!(diagnostics.dropped_sequence_base_cap, 0);
        assert_eq!(diagnostics.dropped_marker_cap, 0);
        assert!(diagnostics.paths.is_empty());
        assert!(diagnostics.retention_truncated);
    }

    #[test]
    fn diagnostic_sequence_and_marker_budgets_never_emit_partial_records() {
        let graph = sequence_graph("ACGTA", 3);
        let path_nodes: Vec<NodeIndex> = graph.node_indices().collect();
        let path = vec![
            (path_nodes[0], None),
            (path_nodes[1], graph.find_edge(path_nodes[0], path_nodes[1])),
            (path_nodes[2], graph.find_edge(path_nodes[1], path_nodes[2])),
            (path_nodes[3], graph.find_edge(path_nodes[2], path_nodes[3])),
        ];
        let marked_sub_kmer = graph[path_nodes[1]].sub_kmer;
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]),
            cyclic_scc_sub_kmers: AHashSet::new(),
            retained_collision_edges: AHashSet::new(),
        };

        let collect_twice = |limits: WithheldDiagnosticLimits| {
            let mut budget = WithheldDiagnosticBudget::new(limits);
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            let visit_counts = path_visit_counts(&path);
            for _ in 0..2 {
                collector.observe_withheld_path(
                    &graph,
                    &path,
                    3,
                    &markers,
                    &visit_counts,
                    cause_counts(1, 0, 0),
                );
            }
            collector.finish()
        };

        let mut base_limits = diagnostic_limits(4, 9, 10);
        base_limits.threshold.sequence_bases = 9;
        let base_limited = collect_twice(base_limits);
        assert_eq!(base_limited.retained_paths, 1);
        assert_eq!(base_limited.dropped_sequence_base_cap, 1);
        assert_eq!(base_limited.paths[0].oriented_sequence, "ACGTA");

        let marker_limited = collect_twice(diagnostic_limits(4, 100, 1));
        assert_eq!(marker_limited.retained_paths, 1);
        assert_eq!(marker_limited.dropped_marker_cap, 1);
        assert_eq!(marker_limited.paths[0].marker_occurrences.len(), 1);

        let oversize_sequence = collect_twice(diagnostic_limits(4, 4, 10));
        assert_eq!(oversize_sequence.retained_paths, 0);
        assert_eq!(oversize_sequence.dropped_oversize_sequence, 2);
        assert!(oversize_sequence.paths.is_empty());

        let oversize_markers = collect_twice(diagnostic_limits(4, 100, 0));
        assert_eq!(oversize_markers.retained_paths, 0);
        assert_eq!(oversize_markers.dropped_oversize_marker_set, 2);
        assert!(oversize_markers.paths.is_empty());

        let mut threshold_path_limits = diagnostic_limits(4, 100, 10);
        threshold_path_limits.threshold.paths = 1;
        let path_limited = collect_twice(threshold_path_limits);
        assert_eq!(path_limited.retained_paths, 1);
        assert_eq!(path_limited.dropped_path_cap, 1);
    }

    #[test]
    fn diagnostic_gene_budget_carries_across_thresholds() {
        let graph = sequence_graph("ACGTA", 3);
        let path_nodes: Vec<NodeIndex> = graph.node_indices().collect();
        let path = vec![
            (path_nodes[0], None),
            (path_nodes[1], graph.find_edge(path_nodes[0], path_nodes[1])),
            (path_nodes[2], graph.find_edge(path_nodes[1], path_nodes[2])),
            (path_nodes[3], graph.find_edge(path_nodes[2], path_nodes[3])),
        ];
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([graph[path_nodes[1]].sub_kmer]),
            cyclic_scc_sub_kmers: AHashSet::new(),
            retained_collision_edges: AHashSet::new(),
        };
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(1, 100, 10));
        let visit_counts = path_visit_counts(&path);

        let first = {
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            collector.observe_withheld_path(
                &graph,
                &path,
                3,
                &markers,
                &visit_counts,
                cause_counts(1, 0, 0),
            );
            collector.finish()
        };
        let second = {
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            collector.observe_withheld_path(
                &graph,
                &path,
                3,
                &markers,
                &visit_counts,
                cause_counts(1, 0, 0),
            );
            collector.finish()
        };

        assert_eq!(first.gene_usage_before_threshold.paths, 0);
        assert_eq!(first.gene_usage_after_threshold.paths, 1);
        assert_eq!(second.gene_usage_before_threshold.paths, 1);
        assert_eq!(second.gene_usage_after_threshold.paths, 1);
        assert_eq!(second.retained_paths, 0);
        assert_eq!(second.dropped_path_cap, 1);
    }

    #[test]
    fn deterministic_gene_allocations_stay_within_gene_and_run_caps() {
        let gene_count = WITHHELD_DIAGNOSTIC_RUN_PATH_CAP + 10;
        let allocations: Vec<_> = (0..gene_count)
            .map(|gene_index| withheld_diagnostic_limits_for_gene(gene_index, gene_count))
            .collect();

        assert!(allocations.iter().all(|limits| {
            limits.allocated_run_share_for_gene.paths <= limits.gene.paths
                && limits.allocated_run_share_for_gene.sequence_bases <= limits.gene.sequence_bases
                && limits.allocated_run_share_for_gene.marker_occurrences
                    <= limits.gene.marker_occurrences
        }));
        assert!(
            allocations
                .iter()
                .map(|limits| limits.allocated_run_share_for_gene.paths)
                .sum::<usize>()
                <= WITHHELD_DIAGNOSTIC_RUN_PATH_CAP
        );
        assert!(
            allocations
                .iter()
                .map(|limits| limits.allocated_run_share_for_gene.sequence_bases)
                .sum::<usize>()
                <= WITHHELD_DIAGNOSTIC_RUN_BASE_CAP
        );
        assert!(
            allocations
                .iter()
                .map(|limits| limits.allocated_run_share_for_gene.marker_occurrences)
                .sum::<usize>()
                <= WITHHELD_DIAGNOSTIC_RUN_MARKER_CAP
        );
        assert_eq!(
            allocations[WITHHELD_DIAGNOSTIC_RUN_PATH_CAP]
                .allocated_run_share_for_gene
                .paths,
            0
        );
        assert_eq!(
            allocations,
            (0..gene_count)
                .map(|gene_index| withheld_diagnostic_limits_for_gene(gene_index, gene_count))
                .collect::<Vec<_>>()
        );
    }

    /// Diamond graph: start -> {a, b} -> end. Should find two paths.
    #[test]
    fn test_diamond_finds_both_paths() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let e = graph.add_node(mk_node(3, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(s, b, mk_edge(5));
        graph.add_edge(a, e, mk_edge(10));
        graph.add_edge(b, e, mk_edge(5));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert_eq!(paths.len(), 2);
    }

    /// Empty graph should produce no paths.
    #[test]
    fn test_no_start_nodes_gives_empty() {
        let graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let params = test_params(0, 100);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert!(paths.is_empty());
    }

    /// max_path_nodes caps path length correctly.
    #[test]
    fn test_max_length_caps_paths() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let b = graph.add_node(mk_node(2, false, false));
        let c = graph.add_node(mk_node(3, false, false));
        let e = graph.add_node(mk_node(4, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, b, mk_edge(10));
        graph.add_edge(b, c, mk_edge(10));
        graph.add_edge(c, e, mk_edge(10));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        // 5 nodes = 5 + 3 - 2 = 6 bases. Set max_length to 5 (too short).
        // max_path_nodes = 5 - 3 + 2 = 4. Path needs 5 nodes, so no paths found.
        let params = test_params(0, 5);
        let paths = get_assembly_paths(&graph, &fkc, &params, None);

        assert!(paths.is_empty());
    }

    /// DFS state budget limits exploration.
    #[test]
    fn test_dfs_budget_limits_exploration() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let a = graph.add_node(mk_node(1, false, false));
        let e = graph.add_node(mk_node(2, false, true));
        graph.add_edge(s, a, mk_edge(10));
        graph.add_edge(a, e, mk_edge(10));

        let mut kc = crate::kmer::KmerCounts::new(&3);
        kc.insert(&0, &10);
        let fkc = kc.filtered_view(1);

        let mut params = test_params(0, 100);
        params.max_dfs_states = 0; // no exploration allowed
        let result =
            search_assembly_paths(&graph, &fkc, &params, None, &RepeatMarkers::default(), None);
        assert!(result.paths.is_empty());
        assert!(result.dfs_limit_reached);
        assert!(!result.path_limit_reached);
    }

    #[test]
    fn test_path_budget_limits_exploration() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let end = graph.add_node(mk_node(2, false, true));
        graph.add_edge(start, end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);

        let mut params = test_params(0, 100);
        params.max_paths_per_pair = 0;
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &params,
            None,
            &RepeatMarkers::default(),
            None,
        );
        assert!(result.paths.is_empty());
        assert!(!result.dfs_limit_reached);
        assert!(result.path_limit_reached);
    }

    #[test]
    fn test_length_search_diagnostics() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let middle = graph.add_node(mk_node(1, false, false));
        let end = graph.add_node(mk_node(2, false, true));
        graph.add_edge(start, middle, mk_edge(10));
        graph.add_edge(middle, end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);

        let too_short = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(5, 10),
            None,
            &RepeatMarkers::default(),
            None,
        );
        assert!(too_short.paths.is_empty());
        assert!(too_short.end_below_min_length);

        let too_long = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(0, 3),
            None,
            &RepeatMarkers::default(),
            None,
        );
        assert!(too_long.paths.is_empty());
        assert!(too_long.max_length_reached);
    }

    #[test]
    fn test_sequence_lengths_at_k_boundary_are_valid() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let middle = graph.add_node(mk_node(1, false, false));
        let end = graph.add_node(mk_node(2, false, true));
        graph.add_edge(start, middle, mk_edge(10));
        graph.add_edge(middle, end, mk_edge(10));
        graph.add_edge(start, end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);
        let length_k = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(3, 3),
            None,
            &RepeatMarkers::default(),
            None,
        );
        let length_k_plus_one = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(4, 4),
            None,
            &RepeatMarkers::default(),
            None,
        );

        assert_eq!(length_k.paths.len(), 1);
        assert_eq!(length_k_plus_one.paths.len(), 1);
    }

    #[test]
    fn test_sequence_generation_rejects_length_above_max() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let middle = graph.add_node(mk_node(1, false, false));
        let end = graph.add_node(mk_node(2, false, true));
        let first_edge = graph.add_edge(start, middle, mk_edge(10));
        let second_edge = graph.add_edge(middle, end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);
        let path = vec![
            (start, None),
            (middle, Some(first_edge)),
            (end, Some(second_edge)),
        ];

        let (records, _) = generate_sequences_from_paths(
            &graph,
            vec![path],
            &filtered,
            "test",
            &test_params(0, 3),
            0,
            None,
        )
        .unwrap();

        assert!(records.is_empty());
    }

    #[test]
    fn test_end_above_max_is_not_returned_by_search() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let end = graph.add_node(mk_node(1, false, true));
        graph.add_edge(start, end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(0, 2),
            None,
            &RepeatMarkers::default(),
            None,
        );

        assert!(result.paths.is_empty());
        assert!(result.max_length_reached);
    }

    #[test]
    fn test_end_below_min_does_not_stop_longer_valid_path() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let short_end = graph.add_node(mk_node(1, false, true));
        let valid_end = graph.add_node(mk_node(2, false, true));
        graph.add_edge(start, short_end, mk_edge(10));
        graph.add_edge(short_end, valid_end, mk_edge(10));

        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.insert(&0, &10);
        let filtered = counts.filtered_view(1);
        let result = search_assembly_paths(
            &graph,
            &filtered,
            &test_params(4, 4),
            None,
            &RepeatMarkers::default(),
            None,
        );

        assert_eq!(result.paths.len(), 1);
        assert!(result.end_below_min_length);
        assert_eq!(result.paths[0].last().unwrap().0, valid_end);
    }

    /// sorted_children returns edges in ascending score order (so pop gives highest).
    #[test]
    fn test_sorted_children_order() {
        let mut graph = StableDiGraph::new();
        let s = graph.add_node(mk_node(0, true, false));
        let lo = graph.add_node(mk_node(1, false, false));
        let hi = graph.add_node(mk_node(2, false, false));
        let lo_edge = graph.add_edge(s, lo, mk_edge(1));
        let hi_edge = graph.add_edge(s, hi, mk_edge(100));

        let children = sorted_children(&graph, s, None);
        assert_eq!(children.len(), 2);
        // Ascending: low first, high second (so pop gives high)
        assert_eq!(children[0].0, lo);
        assert_eq!(children[0].1, lo_edge);
        assert_eq!(children[1].0, hi);
        assert_eq!(children[1].1, hi_edge);
    }
}

// pcr/paths.rs — path finding, sequence extraction, dedup

use super::memory::buffers::BudgetNodeSet;
#[cfg(test)]
use crate::kmer::budget::collections::bytes;
use crate::kmer::budget::collections::{BudgetMap, BudgetVec};
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{MemoryClass, ResourceBudget};
#[cfg(test)]
use ahash::AHashMap;
#[cfg(test)]
use ahash::AHashSet;
use anyhow::{Context, Result};
#[cfg(test)]
use bio::io::fasta;
use log::debug;
use petgraph::Direction;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::stable_graph::StableDiGraph;
use petgraph::visit::EdgeRef;
use sha2::{Digest, Sha256};
use std::sync::Arc;

use crate::kmer::FilteredKmerCounts;

use super::graph::{RepeatMarkers, compute_mean, median_u64_slice};
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
    pub paths: BudgetVec<BudgetVec<PathStep>>,
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
    sequence_first: bool,
    legacy_limits: WithheldDiagnosticLimits,
    legacy_usage: WithheldDiagnosticUsage,
    extra_limits: WithheldDiagnosticLimits,
    extra_usage: WithheldDiagnosticUsage,
    result_budget: Option<Arc<ResourceBudget>>,
}

impl WithheldDiagnosticBudget {
    pub(super) fn new(limits: WithheldDiagnosticLimits) -> Self {
        Self {
            limits,
            usage: WithheldDiagnosticUsage::default(),
            sequence_first: false,
            legacy_limits: limits,
            legacy_usage: WithheldDiagnosticUsage::default(),
            extra_limits: WithheldDiagnosticLimits::default(),
            extra_usage: WithheldDiagnosticUsage::default(),
            result_budget: None,
        }
    }

    pub(super) fn for_emission(legacy_limits: WithheldDiagnosticLimits) -> Self {
        let mut extra_limits = legacy_limits;
        extra_limits.run.paths = extra_limits
            .run
            .paths
            .saturating_mul(super::UNCERTAIN_SEQUENCE_PATH_SHARE_MULTIPLIER);
        extra_limits.allocated_run_share_for_gene.paths = extra_limits
            .allocated_run_share_for_gene
            .paths
            .saturating_mul(super::UNCERTAIN_SEQUENCE_PATH_SHARE_MULTIPLIER)
            .min(extra_limits.gene.paths);
        Self {
            limits: super::inclusive_sequence_limits(legacy_limits),
            usage: WithheldDiagnosticUsage::default(),
            sequence_first: true,
            legacy_limits,
            legacy_usage: WithheldDiagnosticUsage::default(),
            extra_limits,
            extra_usage: WithheldDiagnosticUsage::default(),
            result_budget: None,
        }
    }

    pub(super) fn result_budget_from(&mut self, budget: Option<&Arc<ResourceBudget>>) {
        self.result_budget = budget.cloned();
    }

    pub(super) fn empty_threshold_diagnostics(
        &self,
        node_visit_limit: usize,
    ) -> WithheldPathDiagnostics {
        WithheldPathDiagnostics {
            schema_version: if self.sequence_first { 2 } else { 1 },
            coordinate_system: "zero_based_half_open",
            purpose: "diagnostic_only_unsupported_candidate_hypotheses",
            marker_span_scope: if self.sequence_first {
                "marker_covered_positions_only_for_complete_detail_sources"
            } else {
                "marker_covered_positions_not_complete_ambiguity_or_bridge_ready_intervals"
            },
            read_support_evaluation: "not_evaluated",
            limits: self.limits,
            gene_usage_before_threshold: self.usage,
            gene_usage_after_threshold: self.usage,
            node_visit_limit,
            paths: BudgetVec::in_class(self.result_budget.as_ref(), MemoryClass::RetainedResults),
            ..WithheldPathDiagnostics::default()
        }
    }
}

struct ThresholdDiagnosticCollector<'a> {
    budget: &'a mut WithheldDiagnosticBudget,
    payload: WithheldPathDiagnostics,
    legacy_threshold_usage: WithheldDiagnosticUsage,
    extra_threshold_usage: WithheldDiagnosticUsage,
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

fn retention_drop_reason(
    limits: WithheldDiagnosticLimits,
    threshold_usage: WithheldDiagnosticUsage,
    gene_usage: WithheldDiagnosticUsage,
    sequence_length: usize,
    marker_count: usize,
    require_marker_detail: bool,
) -> Option<&'static str> {
    let threshold = limits.threshold;
    let allocated = limits.allocated_run_share_for_gene;
    if sequence_length > threshold.sequence_bases || sequence_length > allocated.sequence_bases {
        Some("oversize_sequence")
    } else if require_marker_detail
        && (marker_count > threshold.marker_occurrences
            || marker_count > allocated.marker_occurrences)
    {
        Some("oversize_marker_set")
    } else if threshold_usage.paths >= threshold.paths || gene_usage.paths >= allocated.paths {
        Some("path_cap")
    } else if threshold_usage
        .sequence_bases
        .saturating_add(sequence_length)
        > threshold.sequence_bases
        || gene_usage.sequence_bases.saturating_add(sequence_length) > allocated.sequence_bases
    {
        Some("sequence_base_cap")
    } else if require_marker_detail
        && (threshold_usage
            .marker_occurrences
            .saturating_add(marker_count)
            > threshold.marker_occurrences
            || gene_usage.marker_occurrences.saturating_add(marker_count)
                > allocated.marker_occurrences)
    {
        Some("marker_cap")
    } else {
        None
    }
}

impl<'a> ThresholdDiagnosticCollector<'a> {
    fn new(budget: &'a mut WithheldDiagnosticBudget, node_visit_limit: usize) -> Self {
        let payload = budget.empty_threshold_diagnostics(node_visit_limit);
        Self {
            budget,
            payload,
            legacy_threshold_usage: WithheldDiagnosticUsage::default(),
            extra_threshold_usage: WithheldDiagnosticUsage::default(),
        }
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
        visit_counts: &BudgetMap<NodeIndex, usize>,
        cause_counts: WithheldPathCauseCounts,
    ) -> Result<()> {
        if self.budget.sequence_first {
            self.observe_sequence_first(
                graph,
                path,
                kmer_length,
                repeat_markers,
                visit_counts,
                cause_counts,
            )?;
            return Ok(());
        }
        self.payload.observed_withheld_paths += 1;
        let sequence_length = path.len().saturating_add(kmer_length).saturating_sub(2);
        let marker_occurrence_count = cause_counts.total();
        let allocated = self.budget.limits.allocated_run_share_for_gene;
        let threshold = self.budget.limits.threshold;

        let drop_reason = if sequence_length > threshold.sequence_bases
            || sequence_length > allocated.sequence_bases
        {
            Some("oversize_sequence")
        } else if !self.budget.sequence_first
            && (marker_occurrence_count > threshold.marker_occurrences
                || marker_occurrence_count > allocated.marker_occurrences)
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
        } else if !self.budget.sequence_first
            && (self
                .payload
                .retained_marker_occurrences
                .saturating_add(marker_occurrence_count)
                > threshold.marker_occurrences
                || self
                    .budget
                    .usage
                    .marker_occurrences
                    .saturating_add(marker_occurrence_count)
                    > allocated.marker_occurrences)
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
                _ => anyhow::bail!("Invalid diagnostic retention reason: {drop_reason}"),
            }
            self.payload.retention_truncated = true;
            return Ok(());
        }

        let oriented_sequence = reconstruct_path_sequence_budgeted(
            graph,
            path,
            kmer_length,
            self.budget.result_budget.as_ref(),
        )?;
        let detail_fits = self
            .payload
            .retained_marker_occurrences
            .saturating_add(marker_occurrence_count)
            <= threshold.marker_occurrences
            && self
                .budget
                .usage
                .marker_occurrences
                .saturating_add(marker_occurrence_count)
                <= allocated.marker_occurrences;
        let (marker_occurrences, marker_covered_runs) = if detail_fits {
            let occurrences = collect_marker_occurrences(
                graph,
                path,
                kmer_length,
                repeat_markers,
                &oriented_sequence,
                self.budget.result_budget.as_ref(),
            )?;
            debug_assert_eq!(occurrences.len(), marker_occurrence_count);
            let runs = merge_marker_covered_runs(&occurrences, self.budget.result_budget.as_ref())?;
            (occurrences, runs)
        } else {
            self.payload.retained_summary_only_paths += 1;
            self.payload.omitted_detailed_marker_occurrences += marker_occurrence_count;
            self.payload.marker_detail_truncated = true;
            (
                BudgetVec::in_class(
                    self.budget.result_budget.as_ref(),
                    MemoryClass::RetainedResults,
                ),
                BudgetVec::in_class(
                    self.budget.result_budget.as_ref(),
                    MemoryClass::RetainedResults,
                ),
            )
        };
        let max_observed_node_visits = path
            .iter()
            .filter_map(|(node, _)| visit_counts.get(node).copied())
            .max()
            .unwrap_or(0);
        let sequence_sha256 =
            sha256_text_budgeted(&oriented_sequence, self.budget.result_budget.as_ref())?;

        self.payload.retained_paths += 1;
        self.payload.retained_sequence_bases += sequence_length;
        if detail_fits {
            self.payload.retained_marker_occurrences += marker_occurrence_count;
            self.budget.usage.marker_occurrences += marker_occurrence_count;
        }
        self.budget.usage.paths += 1;
        self.budget.usage.sequence_bases += sequence_length;
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
        })?;
        Ok(())
    }

    fn observe_sequence_first(
        &mut self,
        graph: &StableDiGraph<DBNode, DBEdge>,
        path: &[PathStep],
        kmer_length: usize,
        repeat_markers: &RepeatMarkers,
        visit_counts: &BudgetMap<NodeIndex, usize>,
        cause_counts: WithheldPathCauseCounts,
    ) -> Result<()> {
        self.payload.observed_withheld_paths += 1;
        let sequence_length = path.len().saturating_add(kmer_length).saturating_sub(2);
        let marker_count = cause_counts.total();
        debug_assert!(marker_count > 0);
        let legacy_reason = retention_drop_reason(
            self.budget.legacy_limits,
            self.legacy_threshold_usage,
            self.budget.legacy_usage,
            sequence_length,
            marker_count,
            true,
        );
        let retained_with_detail = legacy_reason.is_none();
        if !retained_with_detail {
            let extra_reason = retention_drop_reason(
                self.budget.extra_limits,
                self.extra_threshold_usage,
                self.budget.extra_usage,
                sequence_length,
                marker_count,
                false,
            );
            if let Some(reason) = extra_reason {
                self.payload.observed_not_retained_paths += 1;
                match reason {
                    "oversize_sequence" => self.payload.dropped_oversize_sequence += 1,
                    "path_cap" => self.payload.dropped_path_cap += 1,
                    "sequence_base_cap" => self.payload.dropped_sequence_base_cap += 1,
                    _ => anyhow::bail!("Invalid sequence retention reason: {reason}"),
                }
                self.payload.retention_truncated = true;
                return Ok(());
            }
        }

        let oriented_sequence = reconstruct_path_sequence_budgeted(
            graph,
            path,
            kmer_length,
            self.budget.result_budget.as_ref(),
        )?;
        let (marker_occurrences, marker_covered_runs) = if retained_with_detail {
            let occurrences = collect_marker_occurrences(
                graph,
                path,
                kmer_length,
                repeat_markers,
                &oriented_sequence,
                self.budget.result_budget.as_ref(),
            )?;
            debug_assert_eq!(occurrences.len(), marker_count);
            let runs = merge_marker_covered_runs(&occurrences, self.budget.result_budget.as_ref())?;
            self.legacy_threshold_usage.paths += 1;
            self.legacy_threshold_usage.sequence_bases += sequence_length;
            self.legacy_threshold_usage.marker_occurrences += marker_count;
            self.budget.legacy_usage.paths += 1;
            self.budget.legacy_usage.sequence_bases += sequence_length;
            self.budget.legacy_usage.marker_occurrences += marker_count;
            self.payload.retained_marker_occurrences += marker_count;
            self.budget.usage.marker_occurrences += marker_count;
            (occurrences, runs)
        } else {
            self.extra_threshold_usage.paths += 1;
            self.extra_threshold_usage.sequence_bases += sequence_length;
            self.budget.extra_usage.paths += 1;
            self.budget.extra_usage.sequence_bases += sequence_length;
            self.payload.retained_summary_only_paths += 1;
            self.payload.omitted_detailed_marker_occurrences += marker_count;
            self.payload.marker_detail_truncated = true;
            (
                BudgetVec::in_class(
                    self.budget.result_budget.as_ref(),
                    MemoryClass::RetainedResults,
                ),
                BudgetVec::in_class(
                    self.budget.result_budget.as_ref(),
                    MemoryClass::RetainedResults,
                ),
            )
        };
        let max_observed_node_visits = path
            .iter()
            .filter_map(|(node, _)| visit_counts.get(node).copied())
            .max()
            .unwrap_or(0);
        let sequence_sha256 =
            sha256_text_budgeted(&oriented_sequence, self.budget.result_budget.as_ref())?;
        self.payload.retained_paths += 1;
        self.payload.retained_sequence_bases += sequence_length;
        self.budget.usage.paths += 1;
        self.budget.usage.sequence_bases += sequence_length;
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
        })?;
        Ok(())
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

#[cfg(test)]
fn reconstruct_path_sequence(
    graph: &StableDiGraph<DBNode, DBEdge>,
    path: &[PathStep],
    k: usize,
) -> BudgetString {
    reconstruct_path_sequence_budgeted(graph, path, k, None).unwrap()
}

fn encoded_text(
    identity: u64,
    length: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetString> {
    let mut text = BudgetString::new(budget, MemoryClass::RetainedResults);
    for position in (0..length).rev() {
        text.push_char(char::from(
            b"ACGT"[((identity >> (2 * position)) & 3) as usize],
        ))?;
    }
    Ok(text)
}

fn reconstruct_path_sequence_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    path: &[PathStep],
    k: usize,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetString> {
    let mut sequence = BudgetString::new(budget, MemoryClass::RetainedResults);
    for (position, (node, _)) in path.iter().enumerate() {
        if position == 0 {
            sequence = encoded_text(graph[*node].sub_kmer, k - 1, budget)?;
        } else {
            sequence.push_char(crate::kmer::kmer_last_base(&graph[*node].sub_kmer))?;
        }
    }
    Ok(sequence)
}

#[cfg(test)]
fn sha256_text(value: &str) -> BudgetString {
    sha256_text_budgeted(value, None).unwrap()
}

fn sha256_text_budgeted(value: &str, budget: Option<&Arc<ResourceBudget>>) -> Result<BudgetString> {
    let mut hasher = Sha256::new();
    hasher.update(value.as_bytes());
    BudgetString::format(
        format_args!("{:x}", hasher.finalize()),
        budget,
        MemoryClass::RetainedResults,
    )
}

fn collect_marker_occurrences(
    graph: &StableDiGraph<DBNode, DBEdge>,
    path: &[PathStep],
    k: usize,
    markers: &RepeatMarkers,
    sequence: &str,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<WithheldMarkerOccurrence>> {
    let mut occurrences = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
    for (position, (node, incoming)) in path.iter().enumerate() {
        let identity = graph[*node].sub_kmer;
        for cause in [
            markers
                .omitted_self_loop_sub_kmers
                .contains(&identity)
                .then_some(WithheldMarkerCause::PrePruneOmittedSelfLoopNode),
            markers
                .cyclic_scc_sub_kmers
                .contains(&identity)
                .then_some(WithheldMarkerCause::PrePruneCyclicSccNode),
        ]
        .into_iter()
        .flatten()
        {
            let text = encoded_text(identity, k - 1, budget)?;
            debug_assert_eq!(
                sequence.get(position..position + k - 1),
                Some(text.as_str())
            );
            occurrences.push(WithheldMarkerOccurrence {
                cause,
                identity_kind: WithheldMarkerIdentityKind::OrientedNodeSequence,
                identity_sha256: sha256_text_budgeted(&text, budget)?,
                identity_sequence: text,
                start: position,
                end: position + k - 1,
            })?;
        }
        if let Some(edge) = incoming
            && markers.retained_collision_edges.contains(edge)
        {
            let text = encoded_text(super::graph::reconstruct_edge_kmer(graph, *edge), k, budget)?;
            debug_assert_eq!(
                sequence.get(position - 1..position + k - 1),
                Some(text.as_str())
            );
            occurrences.push(WithheldMarkerOccurrence {
                cause: WithheldMarkerCause::RetainedCollisionEdge,
                identity_kind: WithheldMarkerIdentityKind::OrientedEdgeSequence,
                identity_sha256: sha256_text_budgeted(&text, budget)?,
                identity_sequence: text,
                start: position - 1,
                end: position + k - 1,
            })?;
        }
    }
    Ok(occurrences)
}

fn merge_marker_covered_runs(
    occurrences: &[WithheldMarkerOccurrence],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<WithheldMarkerCoveredRun>> {
    let mut spans = BudgetVec::collect_in_class(
        occurrences.iter().map(|o| (o.start, o.end)),
        budget,
        MemoryClass::RetainedResults,
    )?;
    spans.sort_unstable();
    let mut runs: BudgetVec<WithheldMarkerCoveredRun> =
        BudgetVec::in_class(budget, MemoryClass::RetainedResults);
    for (start, end) in spans {
        if let Some(previous) = runs.last_mut()
            && start <= previous.end
        {
            previous.end = previous.end.max(end);
        } else {
            runs.push(WithheldMarkerCoveredRun { start, end })?;
        }
    }
    Ok(runs)
}

/// Children-of-a-node, sorted by score, used as a DFS frame in `child_stack`.
/// In a de Bruijn graph each node's outgoing edge set is bounded by ≤4 (one
/// per possible 2-bit suffix base), so ordinary frames fit inline. Overflow
/// in manually constructed graphs has separately charged storage. This replaces
/// a per-call `Vec` heap allocation that fired at every visited node during
/// DFS — up to `max_dfs_states = 100K` times per gene.
type Child = (NodeIndex, EdgeIndex, f64);

enum ChildFrame {
    Inline { children: [Child; 4], len: usize },
    Heap(BudgetVec<Child>),
}

impl ChildFrame {
    fn new() -> Self {
        Self::Inline {
            children: [(NodeIndex::new(0), EdgeIndex::new(0), 0.0); 4],
            len: 0,
        }
    }

    fn as_slice(&self) -> &[Child] {
        match self {
            Self::Inline { children, len } => &children[..*len],
            Self::Heap(values) => values,
        }
    }

    fn as_mut_slice(&mut self) -> &mut [Child] {
        match self {
            Self::Inline { children, len } => &mut children[..*len],
            Self::Heap(values) => values,
        }
    }

    fn push(&mut self, child: Child, budget: Option<&Arc<ResourceBudget>>) -> Result<()> {
        match self {
            Self::Inline { children, len } if *len < 4 => {
                children[*len] = child;
                *len += 1;
            }
            Self::Inline { children, .. } => {
                let mut values = BudgetVec::with_capacity(8, budget)?;
                for &old in children.iter() {
                    values.push(old)?;
                }
                values.push(child)?;
                *self = Self::Heap(values);
            }
            Self::Heap(values) => values.push(child)?,
        }
        Ok(())
    }

    fn pop(&mut self) -> Option<Child> {
        match self {
            Self::Inline { children, len } if *len > 0 => {
                *len -= 1;
                Some(children[*len])
            }
            Self::Inline { .. } => None,
            Self::Heap(values) => values.pop(),
        }
    }
}

impl std::ops::Deref for ChildFrame {
    type Target = [Child];
    fn deref(&self) -> &[Child] {
        self.as_slice()
    }
}

/// Keep four ordinary de Bruijn branches inline. Manually constructed graphs
/// with more branches retain a charged spill buffer and charged sort scratch.
fn sorted_children_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    edge_preferences: Option<&BudgetMap<EdgeIndex, f64>>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<ChildFrame> {
    let mut frame = ChildFrame::new();
    for edge in graph.edges_directed(node, Direction::Outgoing) {
        let preference = edge_preferences
            .and_then(|prefs| prefs.get(&edge.id()))
            .copied()
            .unwrap_or(1.0);
        frame.push(
            (
                edge.target(),
                edge.id(),
                edge.weight().count as f64 * preference,
            ),
            budget,
        )?;
    }
    if frame.len() <= 4 {
        // Stable in-place insertion over at most four entries has no heap or
        // unbounded quadratic cost. Preserve the previous score-tie order.
        for index in 1..frame.len() {
            let mut index = index;
            while index > 0 && frame[index - 1].2.total_cmp(&frame[index].2).is_gt() {
                frame.as_mut_slice().swap(index - 1, index);
                index -= 1;
            }
        }
    } else {
        let mut indexed = BudgetVec::collect(frame.iter().copied().enumerate(), budget)?;
        indexed.sort_unstable_by(|a, b| a.1.2.total_cmp(&b.1.2).then(a.0.cmp(&b.0)));
        for (destination, (_, child)) in frame.as_mut_slice().iter_mut().zip(indexed.iter()) {
            *destination = *child;
        }
    }
    Ok(frame)
}

#[cfg(test)]
fn test_preferences(
    source: Option<&AHashMap<EdgeIndex, f64>>,
) -> Option<BudgetMap<EdgeIndex, f64>> {
    source.map(|source| {
        let mut result = BudgetMap::new(None);
        for (&edge, &score) in source {
            result.insert(edge, score).unwrap();
        }
        result
    })
}

#[cfg(test)]
fn sorted_children(
    graph: &StableDiGraph<DBNode, DBEdge>,
    node: NodeIndex,
    preferences: Option<&AHashMap<EdgeIndex, f64>>,
) -> ChildFrame {
    let preferences = test_preferences(preferences);
    sorted_children_budgeted(graph, node, preferences.as_ref(), None).unwrap()
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
    .into_iter()
    .map(|path| path.into_iter().collect())
    .collect()
}

#[cfg(test)]
pub(super) fn search_assembly_paths(
    graph: &StableDiGraph<DBNode, DBEdge>,
    counts: &FilteredKmerCounts,
    params: &PCRParams,
    preferences: Option<&AHashMap<EdgeIndex, f64>>,
    markers: &RepeatMarkers,
    diagnostics: Option<&mut WithheldDiagnosticBudget>,
) -> PathSearchResult {
    let preferences = test_preferences(preferences);
    search_assembly_paths_budgeted(
        graph,
        counts,
        params,
        preferences.as_ref(),
        markers,
        diagnostics,
        None,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
pub(super) fn search_assembly_paths_budgeted(
    graph: &StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
    params: &PCRParams,
    edge_preferences: Option<&BudgetMap<EdgeIndex, f64>>,
    repeat_markers: &RepeatMarkers,
    diagnostic_budget: Option<&mut WithheldDiagnosticBudget>,
    buffer_budget: Option<&Arc<ResourceBudget>>,
) -> Result<PathSearchResult> {
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

    let mut end_nodes = BudgetNodeSet::new(buffer_budget);
    for node in graph.node_indices().filter(|&node| graph[node].is_end) {
        end_nodes.insert(node)?;
    }
    let mut result = PathSearchResult {
        paths: BudgetVec::new(buffer_budget),
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

    for start in graph.node_indices().filter(|&node| graph[node].is_start) {
        let mut paths_from_start = 0;
        let mut states_explored: usize = 0;

        // Stack-based DFS with push/pop backtracking instead of cloning.
        // `path` and `visit_counts` are maintained incrementally.
        // `child_stack` holds remaining children to explore at each depth.
        // The starting step has edge = None; all subsequent steps record
        // the edge that connects the previous node to this one.
        let mut path = BudgetVec::new(buffer_budget);
        path.push((start, None))?;
        let mut visit_counts = BudgetMap::new(buffer_budget);
        visit_counts.insert(start, 1usize)?;
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
        let children = sorted_children_budgeted(graph, start, edge_preferences, buffer_budget)?;
        let mut child_stack = BudgetVec::new(buffer_budget);
        child_stack.push(children)?;

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
                path.push((neighbor, Some(edge_id)))?;
                visit_counts.insert(
                    neighbor,
                    current_visits
                        .checked_add(1)
                        .context("DFS visit count overflow")?,
                )?;
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
                                )?;
                            }
                        } else {
                            result.paths.push(path.try_clone()?)?;
                            result.eligible_candidate_path_count += 1;
                            paths_from_start += 1;
                        }
                        omitted_self_loop_nodes_in_path -=
                            usize::from(neighbor_has_omitted_self_loop);
                        cyclic_scc_nodes_in_path -= usize::from(neighbor_is_in_cyclic_scc);
                        retained_collision_edges_in_path -= usize::from(edge_is_retained_collision);
                        *visit_counts
                            .get_mut(&neighbor)
                            .context("DFS neighbor visit count is missing")? -= 1;
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
                    *visit_counts
                        .get_mut(&neighbor)
                        .context("DFS neighbor visit count is missing")? -= 1;
                    path.pop();
                    continue;
                }

                // Explore deeper: push children frame
                let children =
                    sorted_children_budgeted(graph, neighbor, edge_preferences, buffer_budget)?;
                child_stack.push(children)?;
            } else {
                // No more children at this depth, backtrack
                child_stack.pop();
                if child_stack.is_empty() {
                    break;
                }
                let (backtrack_node, incoming_edge) =
                    path.pop().context("DFS path is empty during backtrack")?;
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
                    .context("DFS backtrack node is missing from visit counts")? -= 1;
            }
        }
    }

    result.withheld_path_diagnostics =
        diagnostic_collector.map(ThresholdDiagnosticCollector::finish);

    Ok(result)
}

/// Extract sequences from graph paths, producing FASTA assembly records.
/// Returns the generated records and the updated amplicon index.
#[cfg(test)]
pub(super) fn generate_sequences_from_paths<P: AsRef<[PathStep]>>(
    graph: &StableDiGraph<DBNode, DBEdge>,
    all_paths: impl IntoIterator<Item = P>,
    counts: &FilteredKmerCounts,
    sample: &str,
    params: &PCRParams,
    index: usize,
    threading: Option<&super::threading::ThreadingAnnotations>,
) -> Result<(Vec<AssemblyRecord>, usize)> {
    generate_sequences_from_paths_budgeted(
        graph, all_paths, counts, sample, params, index, threading, None, None,
    )
    .map(|(records, index)| (records.into_iter().collect(), index))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn generate_sequences_from_paths_budgeted<P: AsRef<[PathStep]>>(
    graph: &StableDiGraph<DBNode, DBEdge>,
    all_paths: impl IntoIterator<Item = P>,
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    mut amplicon_index: usize,
    threading: Option<&super::threading::ThreadingAnnotations>,
    budget: Option<&Arc<ResourceBudget>>,
    result_budget: Option<&Arc<ResourceBudget>>,
) -> Result<(BudgetVec<AssemblyRecord>, usize)> {
    let mut assembly_records = BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);

    for owned_path in all_paths.into_iter() {
        let path = owned_path.as_ref();
        let mut sequence = BudgetString::new(result_budget, MemoryClass::RetainedResults);
        let mut edge_counts = BudgetVec::with_capacity(path.len().saturating_sub(1), budget)?; // u64 for compute_mean/median compatibility
        // Edges traversed by this path, in order. Recorded once during the
        // sequence-extraction loop and reused for max_coverage_ratio and
        // threading-support metrics below — replaces three find_edge passes.
        let mut path_edges = BudgetVec::with_capacity(path.len().saturating_sub(1), budget)?;
        for &(node, edge_opt) in path.iter() {
            let node_data = graph
                .node_weight(node)
                .context("Node not found in graph during sequence generation")?;
            if sequence.is_empty() {
                for position in (0..kmer_counts.get_k() - 1).rev() {
                    let base = b"ACGT"[((node_data.sub_kmer >> (position * 2)) & 3) as usize];
                    sequence.push_char(char::from(base))?;
                }
                debug_assert!(
                    edge_opt.is_none(),
                    "first PathStep must have edge = None (this is the start node)"
                );
            } else {
                sequence.push_char(crate::kmer::kmer_last_base(&node_data.sub_kmer))?;
                let edge_idx =
                    edge_opt.context("Non-start PathStep is missing its edge index (DFS bug)")?;
                let edge_data = graph
                    .edge_weight(edge_idx)
                    .context("Edge weight not found")?;
                edge_counts.push(edge_data.count)?;
                path_edges.push(edge_idx)?;
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
        let count_min = edge_counts
            .iter()
            .min()
            .copied()
            .context("No edge counts found for path")?;
        let count_max = edge_counts
            .iter()
            .max()
            .copied()
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

        // CV above consumes counts in the original path order, preserving its
        // floating-point sum. Median can now partition in-place without a copy.
        let count_median = median_u64_slice(&mut edge_counts).context("Missing path median")?;

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
                let mut unambiguous_counts = BudgetVec::with_capacity(path_edges.len(), budget)?;

                for &edge_idx in path_edges.iter() {
                    total_edges += 1;
                    match ann.edge_support.get(&edge_idx) {
                        Some(s) if s.read_support_total > 0 => {
                            supported_edges += 1;
                            unambiguous_counts.push(s.read_support_unambiguous as u64)?;
                        }
                        _ => {
                            zero_count += 1;
                            unambiguous_counts.push(0)?;
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
                    median_u64_slice(&mut unambiguous_counts)
                        .context("Missing threading-support median")?
                };

                (Some(zero_count), Some(median_unamb), Some(frac))
            } else {
                (None, None, None)
            };

        let score = super::PathScore {
            kmer_min_count: count_min,
            kmer_median_count: count_median,
            coverage_cv,
            max_coverage_ratio,
            zero_support_edges,
            median_unambiguous_support,
            edge_support_fraction,
        };

        let id = BudgetString::format(
            format_args!("{}_{}_{}", sample_name, params.gene_name, amplicon_index),
            result_budget,
            MemoryClass::RetainedResults,
        )?;
        let desc = BudgetString::format(
            format_args!(
                "sample={} target={} sequence_index={} sequence_length={} kmer_count_mean={:.2} kmer_count_median={} kmer_count_min={} kmer_count_max={}",
                sample_name,
                params.gene_name,
                amplicon_index,
                sequence.len(),
                count_mean,
                count_median,
                count_min,
                count_max
            ),
            result_budget,
            MemoryClass::RetainedResults,
        )?;

        amplicon_index += 1;

        debug!(">{} {}", id, desc);
        let record = super::sequence::SequenceRecord::from_parts(id, Some(desc), sequence);
        assembly_records.push(AssemblyRecord {
            fasta_record: record,
            score,
        })?;
    }

    Ok((assembly_records, amplicon_index))
}

/// Sort assembly records by descending composite score, deduplicate near-identical
/// sequences, and truncate to the maximum number of amplicons.
#[cfg(test)]
pub(super) fn sort_and_deduplicate(
    records: Vec<AssemblyRecord>,
    params: &PCRParams,
) -> Vec<AssemblyRecord> {
    let records = BudgetVec::collect(records, None).unwrap();
    sort_and_deduplicate_budgeted(records, params)
        .unwrap()
        .into_iter()
        .collect()
}

pub(super) fn sort_and_deduplicate_budgeted(
    assembly_records: BudgetVec<AssemblyRecord>,
    params: &PCRParams,
) -> Result<BudgetVec<AssemblyRecord>> {
    let budget = assembly_records.budget().cloned();
    let class = assembly_records.class();
    let mut sorted = BudgetVec::in_class(budget.as_ref(), class);
    for item in assembly_records.into_iter().enumerate() {
        sorted.push(item)?;
    }
    sorted.sort_unstable_by(|(a_index, a), (b_index, b)| {
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
            .then(a_index.cmp(b_index))
    });

    let records = sorted;
    let num_records_all = records.len();

    // Greedy clustering: keep each record only if it is not within
    // dedup_edit_threshold edits of any previously kept record.
    // Distances are computed on-the-fly to avoid O(N^2) memory.
    let mut kept: BudgetVec<AssemblyRecord> = BudgetVec::in_class(budget.as_ref(), class);
    for (_, record) in records {
        let mut is_duplicate = false;
        for kept_record in &kept {
            if super::distance::bounded_levenshtein(
                record.fasta_record.seq(),
                kept_record.fasta_record.seq(),
                params.dedup_edit_threshold,
                budget.as_ref(),
                class,
            )?
            .is_some()
            {
                is_duplicate = true;
                break;
            }
        }
        if !is_duplicate {
            kept.push(record)?;
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

    Ok(records)
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

    fn mk_edge(count: u64) -> DBEdge {
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

    fn path_visit_counts(path: &[PathStep]) -> BudgetMap<NodeIndex, usize> {
        let mut counts = BudgetMap::new(None);
        for (node, _) in path {
            let count = counts.get(node).copied().unwrap_or(0usize);
            counts.insert(*node, count + 1).unwrap();
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
            emit_uncertain: false,
        }
    }

    fn path_budget(bytes: u64) -> Arc<ResourceBudget> {
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
    fn child_frames_preserve_stable_score_ties_and_charge_only_overflow_storage() {
        assert!(size_of::<ChildFrame>() <= 80);
        for degree in [0, 1, 2, 4, 5, 20, 100] {
            let mut graph = StableDiGraph::new();
            let root = graph.add_node(mk_node(0, true, false));
            for index in 0..degree {
                let node = graph.add_node(mk_node(index as u64, false, false));
                graph.add_edge(root, node, mk_edge((index % 3 + 1) as u64));
            }
            if degree > 4 {
                graph.remove_edge(EdgeIndex::new(2));
                graph.add_edge(root, root, mk_edge(1));
            }
            let preferences: AHashMap<_, _> = graph
                .edge_indices()
                .map(|edge| (edge, if edge.index() % 7 == 0 { f64::NAN } else { 1.0 }))
                .collect();
            for preferences in [None, Some(&preferences)] {
                let mut expected: Vec<_> = graph
                    .edges_directed(root, Direction::Outgoing)
                    .map(|edge| {
                        let preference = preferences
                            .and_then(|prefs| prefs.get(&edge.id()))
                            .copied()
                            .unwrap_or(1.0);
                        (
                            edge.target(),
                            edge.id(),
                            edge.weight().count as f64 * preference,
                        )
                    })
                    .collect();
                expected.sort_by(|a, b| a.2.total_cmp(&b.2));
                let budget = path_budget(16384);
                let owned = test_preferences(preferences);
                let frame =
                    sorted_children_budgeted(&graph, root, owned.as_ref(), Some(&budget)).unwrap();
                assert_eq!(
                    frame.iter().map(|child| child.1).collect::<Vec<_>>(),
                    expected.iter().map(|child| child.1).collect::<Vec<_>>()
                );
                assert_eq!(matches!(frame, ChildFrame::Heap(_)), degree > 4);
                assert_eq!(
                    budget.snapshot().reserved_memory_bytes,
                    if degree > 4 {
                        match &frame {
                            ChildFrame::Heap(values) => bytes::<Child>(values.capacity()).unwrap(),
                            ChildFrame::Inline { .. } => 0,
                        }
                    } else {
                        0
                    }
                );
                drop(frame);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
    }

    #[test]
    fn child_frame_spill_or_sort_refusal_releases_all_storage() {
        let mut graph = StableDiGraph::new();
        let root = graph.add_node(mk_node(0, true, false));
        for index in 0..5 {
            let node = graph.add_node(mk_node(index, false, false));
            graph.add_edge(root, node, mk_edge(1));
        }
        for bytes in [1, 128] {
            let budget = path_budget(bytes);
            let error = sorted_children_budgeted(&graph, root, None, Some(&budget))
                .err()
                .unwrap();
            assert!(
                error
                    .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                    .is_some()
            );
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(graph.edge_count(), 5);
        }
    }

    // Recursive small-graph oracle uses raw containers and stable sorting,
    // independently of the production iterative stack/backtracking owners.
    fn reference_paths(
        graph: &StableDiGraph<DBNode, DBEdge>,
        params: &PCRParams,
        k: usize,
    ) -> Vec<Vec<PathStep>> {
        struct Reference<'a> {
            graph: &'a StableDiGraph<DBNode, DBEdge>,
            params: &'a PCRParams,
            minimum: usize,
            maximum: usize,
            paths: Vec<Vec<PathStep>>,
            states: usize,
            found: usize,
        }
        impl Reference<'_> {
            fn visit(
                &mut self,
                node: NodeIndex,
                path: &mut Vec<PathStep>,
                counts: &mut AHashMap<NodeIndex, usize>,
            ) {
                let mut edges: Vec<_> = self
                    .graph
                    .edges_directed(node, Direction::Outgoing)
                    .map(|edge| (edge.target(), edge.id(), edge.weight().count))
                    .collect();
                edges.sort_by_key(|edge| edge.2);
                for (next, edge, _) in edges.into_iter().rev() {
                    if self.states >= self.params.max_dfs_states
                        || self.found >= self.params.max_paths_per_pair
                    {
                        return;
                    }
                    self.states += 1;
                    let visits = counts.get(&next).copied().unwrap_or(0);
                    if visits >= self.params.max_node_visits {
                        continue;
                    }
                    path.push((next, Some(edge)));
                    counts.insert(next, visits + 1);
                    if self.graph[next].is_end
                        && path.len() >= self.minimum
                        && path.len() <= self.maximum
                    {
                        self.paths.push(path.clone());
                        self.found += 1;
                    } else if path.len() < self.maximum {
                        self.visit(next, path, counts);
                    }
                    counts.insert(next, visits);
                    path.pop();
                }
            }
        }
        let mut reference = Reference {
            graph,
            params,
            minimum: params.min_length.saturating_sub(k.saturating_sub(2)).max(1),
            maximum: params.max_length.saturating_sub(k.saturating_sub(2)).max(1),
            paths: Vec::new(),
            states: 0,
            found: 0,
        };
        for start in graph.node_indices().filter(|&node| graph[node].is_start) {
            reference.states = 0;
            reference.found = 0;
            reference.visit(
                start,
                &mut vec![(start, None)],
                &mut AHashMap::from_iter([(start, 1)]),
            );
        }
        reference.paths
    }

    #[test]
    fn owned_dfs_matches_recursive_path_order_quotas_and_backtracking_with_holes() {
        let mut random = 157u64;
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("AAACCCGGGTTT").unwrap();
        let budget = path_budget(1 << 20);
        for _ in 0..48 {
            let mut graph = StableDiGraph::new();
            for index in 0..9 {
                graph.add_node(mk_node(
                    index,
                    index == 0 || index == 4,
                    index == 7 || index == 8,
                ));
            }
            graph.remove_node(NodeIndex::new(3));
            let nodes: Vec<_> = graph.node_indices().collect();
            for _ in 0..16 {
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let from = nodes[(random >> 32) as usize % nodes.len()];
                random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
                let to = nodes[(random >> 32) as usize % nodes.len()];
                graph.add_edge(from, to, mk_edge(random % 5 + 1));
            }
            for states in [8, 1000] {
                let mut params = test_params(3, 8);
                params.max_paths_per_pair = 3;
                params.max_dfs_states = states;
                params.max_node_visits = 2;
                let expected = reference_paths(&graph, &params, 3);
                let actual = search_assembly_paths_budgeted(
                    &graph,
                    &counts.filtered_view(1),
                    &params,
                    None,
                    &RepeatMarkers::default(),
                    None,
                    Some(&budget),
                )
                .unwrap();
                let observed: Vec<_> = actual
                    .paths
                    .iter()
                    .map(|path| path.iter().copied().collect::<Vec<_>>())
                    .collect();
                assert_eq!(observed, expected);
                drop(actual);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
    }

    #[test]
    fn dfs_or_metric_refusal_cleans_workspace_and_success_keeps_path_ownership() {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let end = graph.add_node(mk_node(1, false, true));
        let edge = graph.add_edge(start, end, mk_edge(4));
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("AAAC").unwrap();
        let params = test_params(3, 8);
        let tiny = path_budget(1);
        let error = search_assembly_paths_budgeted(
            &graph,
            &counts.filtered_view(1),
            &params,
            None,
            &RepeatMarkers::default(),
            None,
            Some(&tiny),
        )
        .err()
        .unwrap();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        let paths = vec![vec![(start, None), (end, Some(edge))]];
        assert!(
            generate_sequences_from_paths_budgeted(
                &graph,
                paths,
                &counts.filtered_view(1),
                "sample",
                &params,
                0,
                None,
                Some(&tiny),
                None,
            )
            .err()
            .unwrap()
            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
            .is_some()
        );
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
        let enough = path_budget(4096);
        let found = search_assembly_paths_budgeted(
            &graph,
            &counts.filtered_view(1),
            &params,
            None,
            &RepeatMarkers::default(),
            None,
            Some(&enough),
        )
        .unwrap();
        assert_eq!(found.paths.len(), 1);
        assert!(enough.snapshot().reserved_memory_bytes > 0);
        let (records, _) = generate_sequences_from_paths_budgeted(
            &graph,
            found.paths,
            &counts.filtered_view(1),
            "sample",
            &params,
            0,
            None,
            Some(&enough),
            None,
        )
        .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(enough.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn ranked_records_preserve_unrounded_winner_through_deduplication_and_limits() {
        let make = |sequence: &[u8], score| AssemblyRecord {
            fasta_record: fasta::Record::with_attrs("temporary", None, sequence).into(),
            score: super::super::PathScore {
                kmer_min_count: 1,
                kmer_median_count: score,
                coverage_cv: 0.0,
                max_coverage_ratio: 1.0,
                zero_support_edges: None,
                median_unambiguous_support: None,
                edge_support_fraction: None,
            },
        };
        let mut params = test_params(0, 100);
        params.dedup_edit_threshold = 1;
        let selected = sort_and_deduplicate(
            vec![
                make(b"AAAAA", 2.0),
                make(b"AAAAT", 13.1254),
                make(b"CCCCC", 8.0001),
            ],
            &params,
        );
        assert_eq!(selected.len(), 2);
        assert_eq!(selected[0].fasta_record.seq(), b"AAAAT");
        assert_eq!(selected[0].score.composite(), 13.1254);
        assert_eq!(selected[1].score.composite(), 8.0001);

        params.dedup_edit_threshold = 0;
        let inputs = (0..25u64)
            .map(|index| {
                let sequence = (0..5)
                    .map(|position| b"ACGT"[((index >> (position * 2)) & 3) as usize])
                    .collect::<Vec<_>>();
                make(&sequence, index as f64 + 0.1234)
            })
            .collect();
        let selected = sort_and_deduplicate(inputs, &params);
        assert_eq!(selected.len(), MAX_NUM_AMPLICONS);
        assert_eq!(selected[0].score.composite(), 24.1234);
        assert_eq!(selected.last().unwrap().score.composite(), 5.1234);
    }

    #[test]
    fn dedup_scratch_refusal_releases_owned_results_across_allocation_limits() {
        let mut params = test_params(0, 100);
        params.dedup_edit_threshold = 64;
        let mut refused = 0;
        let mut scratch_refused = 0;
        let mut succeeded = 0;
        for limit in (64..=16384).step_by(64) {
            let budget = path_budget(limit);
            let result = (|| -> Result<BudgetVec<AssemblyRecord>> {
                let mut input = BudgetVec::in_class(Some(&budget), MemoryClass::RetainedResults);
                for (id, sequence, score) in [
                    ("low", "A".repeat(64), 2.0),
                    ("winner", "C".repeat(64), 13.1254),
                    ("second", "G".repeat(64), 8.0001),
                    ("third", "T".repeat(64), 3.0),
                ] {
                    let text = |s: &str| {
                        BudgetString::from_str(s, Some(&budget), MemoryClass::RetainedResults)
                    };
                    input.push(AssemblyRecord {
                        fasta_record: super::super::sequence::SequenceRecord::from_parts(
                            text(id)?,
                            None,
                            text(&sequence)?,
                        ),
                        score: super::super::PathScore {
                            kmer_min_count: 1,
                            kmer_median_count: score,
                            coverage_cv: 0.0,
                            max_coverage_ratio: 1.0,
                            zero_support_edges: None,
                            median_unambiguous_support: None,
                            edge_support_fraction: None,
                        },
                    })?;
                }
                sort_and_deduplicate_budgeted(input, &params)
            })();
            match &result {
                Ok(records) => {
                    succeeded += 1;
                    assert_eq!(records.len(), 1);
                    assert_eq!(records[0].fasta_record.id(), "winner");
                    assert_eq!(records[0].fasta_record.seq(), vec![b'C'; 64]);
                    assert_eq!(records[0].score.composite(), 13.1254);
                    assert!(
                        budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize]
                            > 0
                    );
                }
                Err(err) => {
                    refused += 1;
                    assert!(
                        err.downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some()
                    );
                    if err.to_string().contains("edit-distance scratch") {
                        scratch_refused += 1;
                    }
                }
            }
            assert_eq!(
                budget.snapshot().peak_reserved_by_class[MemoryClass::TargetGraphs as usize],
                0
            );
            drop(result);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
        assert!(refused > 0);
        assert!(
            scratch_refused > 0,
            "must exercise refusal inside the greedy comparison loop"
        );
        assert!(succeeded > 0);
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
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]).into(),
            cyclic_scc_sub_kmers: AHashSet::from_iter([marked_sub_kmer]).into(),
            retained_collision_edges: AHashSet::from_iter([collision_edge]).into(),
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
            omitted_self_loop_sub_kmers: AHashSet::new().into(),
            cyclic_scc_sub_kmers: AHashSet::from_iter([graph[repeat_first].sub_kmer]).into(),
            retained_collision_edges: AHashSet::new().into(),
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
                omitted_self_loop_sub_kmers: AHashSet::from_iter([marked]).into(),
                cyclic_scc_sub_kmers: AHashSet::new().into(),
                retained_collision_edges: AHashSet::new().into(),
            };
            let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(4, 100, 20));
            search_assembly_paths(graph, &filtered, &params, None, &markers, Some(&mut budget))
                .withheld_path_diagnostics
                .unwrap()
                .paths
                .into_iter()
                .next()
                .unwrap()
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
            omitted_self_loop_sub_kmers: AHashSet::new().into(),
            cyclic_scc_sub_kmers: AHashSet::from_iter([graph[repeat].sub_kmer]).into(),
            retained_collision_edges: AHashSet::new().into(),
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
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]).into(),
            cyclic_scc_sub_kmers: AHashSet::new().into(),
            retained_collision_edges: AHashSet::new().into(),
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
            omitted_self_loop_sub_kmers: AHashSet::from_iter([marked_sub_kmer]).into(),
            cyclic_scc_sub_kmers: AHashSet::new().into(),
            retained_collision_edges: AHashSet::new().into(),
        };

        let collect_twice = |limits: WithheldDiagnosticLimits| {
            let mut budget = WithheldDiagnosticBudget::new(limits);
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            let visit_counts = path_visit_counts(&path);
            for _ in 0..2 {
                collector
                    .observe_withheld_path(
                        &graph,
                        &path,
                        3,
                        &markers,
                        &visit_counts,
                        cause_counts(1, 0, 0),
                    )
                    .unwrap();
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
    fn inclusive_retention_keeps_sequence_when_marker_budget_is_full() {
        let graph = sequence_graph("ACGTA", 3);
        let nodes: Vec<NodeIndex> = graph.node_indices().collect();
        let path = vec![
            (nodes[0], None),
            (nodes[1], graph.find_edge(nodes[0], nodes[1])),
            (nodes[2], graph.find_edge(nodes[1], nodes[2])),
            (nodes[3], graph.find_edge(nodes[2], nodes[3])),
        ];
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([graph[nodes[1]].sub_kmer]).into(),
            cyclic_scc_sub_kmers: AHashSet::new().into(),
            retained_collision_edges: AHashSet::new().into(),
        };
        let mut limits = diagnostic_limits(1, 100, 1);
        limits.threshold.paths = 3;
        let mut budget = WithheldDiagnosticBudget::for_emission(limits);
        let visit_counts = path_visit_counts(&path);
        let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
        for _ in 0..3 {
            collector
                .observe_withheld_path(
                    &graph,
                    &path,
                    3,
                    &markers,
                    &visit_counts,
                    cause_counts(1, 0, 0),
                )
                .unwrap();
        }
        let diagnostic = collector.finish();
        assert_eq!(diagnostic.schema_version, 2);
        assert_eq!(diagnostic.observed_withheld_paths, 3);
        assert_eq!(diagnostic.retained_paths, 2);
        assert_eq!(diagnostic.dropped_path_cap, 1);
        assert_eq!(diagnostic.dropped_marker_cap, 0);
        assert_eq!(diagnostic.retained_marker_occurrences, 1);
        assert_eq!(diagnostic.retained_summary_only_paths, 1);
        assert_eq!(diagnostic.omitted_detailed_marker_occurrences, 1);
        assert!(diagnostic.marker_detail_truncated);
        assert_eq!(diagnostic.paths[0].oriented_sequence, "ACGTA");
        assert_eq!(diagnostic.paths[0].marker_occurrences.len(), 1);
        assert_eq!(diagnostic.paths[1].oriented_sequence, "ACGTA");
        assert!(diagnostic.paths[1].marker_occurrences.is_empty());
        assert!(diagnostic.paths[1].marker_covered_runs.is_empty());
        assert_eq!(diagnostic.paths[1].omitted_self_loop_marker_occurrences, 1);
    }

    #[test]
    fn extra_paths_cannot_displace_a_later_legacy_detail_path() {
        let graph = sequence_graph("ACGTA", 3);
        let nodes: Vec<NodeIndex> = graph.node_indices().collect();
        let path = vec![
            (nodes[0], None),
            (nodes[1], graph.find_edge(nodes[0], nodes[1])),
            (nodes[2], graph.find_edge(nodes[1], nodes[2])),
            (nodes[3], graph.find_edge(nodes[2], nodes[3])),
        ];
        let markers = RepeatMarkers {
            omitted_self_loop_sub_kmers: AHashSet::from_iter([graph[nodes[1]].sub_kmer]).into(),
            cyclic_scc_sub_kmers: AHashSet::new().into(),
            retained_collision_edges: AHashSet::new().into(),
        };
        let mut budget = WithheldDiagnosticBudget::for_emission(diagnostic_limits(1, 100, 1));
        let visit_counts = path_visit_counts(&path);
        let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
        for count in [2, 2, 1] {
            collector
                .observe_withheld_path(
                    &graph,
                    &path,
                    3,
                    &markers,
                    &visit_counts,
                    cause_counts(count, 0, 0),
                )
                .unwrap();
        }
        let diagnostic = collector.finish();
        assert_eq!(diagnostic.observed_withheld_paths, 3);
        assert_eq!(diagnostic.retained_paths, 2);
        assert_eq!(diagnostic.retained_summary_only_paths, 1);
        assert_eq!(diagnostic.dropped_path_cap, 1);
        assert!(diagnostic.paths[0].marker_occurrences.is_empty());
        assert_eq!(diagnostic.paths[1].marker_occurrences.len(), 1);
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
            omitted_self_loop_sub_kmers: AHashSet::from_iter([graph[path_nodes[1]].sub_kmer])
                .into(),
            cyclic_scc_sub_kmers: AHashSet::new().into(),
            retained_collision_edges: AHashSet::new().into(),
        };
        let mut budget = WithheldDiagnosticBudget::new(diagnostic_limits(1, 100, 10));
        let visit_counts = path_visit_counts(&path);

        let first = {
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            collector
                .observe_withheld_path(
                    &graph,
                    &path,
                    3,
                    &markers,
                    &visit_counts,
                    cause_counts(1, 0, 0),
                )
                .unwrap();
            collector.finish()
        };
        let second = {
            let mut collector = ThresholdDiagnosticCollector::new(&mut budget, 2);
            collector
                .observe_withheld_path(
                    &graph,
                    &path,
                    3,
                    &markers,
                    &visit_counts,
                    cause_counts(1, 0, 0),
                )
                .unwrap();
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

    #[test]
    fn inclusive_sequence_shares_are_bounded_and_fit_full_panel_12s() {
        let allocations: Vec<_> = (0..21)
            .map(|index| {
                super::super::inclusive_sequence_limits(withheld_diagnostic_limits_for_gene(
                    index, 21,
                ))
            })
            .collect();
        assert_eq!(allocations[0].run.paths, 1024);
        assert_eq!(
            allocations
                .iter()
                .map(|limits| limits.allocated_run_share_for_gene.paths)
                .sum::<usize>(),
            1024
        );
        assert!(allocations.iter().all(|limits| {
            limits.allocated_run_share_for_gene.paths >= 32
                && limits.threshold.paths == 64
                && limits.allocated_run_share_for_gene.marker_occurrences <= 1561
        }));
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
    fn diagnostic_result_budget(bytes: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = bytes;
        ResourceBudget::new(crate::kmer::budget::ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }

    fn owned_diagnostic_search(
        budget: &Arc<ResourceBudget>,
        limits: WithheldDiagnosticLimits,
    ) -> Result<PathSearchResult> {
        let mut graph = StableDiGraph::new();
        let start = graph.add_node(mk_node(0, true, false));
        let marked = graph.add_node(mk_node(1, false, false));
        let end = graph.add_node(mk_node(6, false, true));
        graph.add_edge(start, marked, mk_edge(4));
        graph.add_edge(marked, end, mk_edge(4));
        let mut markers = RepeatMarkers::default();
        markers.cyclic_scc_sub_kmers.insert(1).unwrap();
        let mut counts = crate::kmer::KmerCounts::new(&3);
        counts.ingest_seq("AACG").unwrap();
        let mut diagnostics = WithheldDiagnosticBudget::for_emission(limits);
        diagnostics.result_budget_from(Some(budget));
        search_assembly_paths_budgeted(
            &graph,
            &counts.filtered_view(1),
            &test_params(4, 8),
            None,
            &markers,
            Some(&mut diagnostics),
            None,
        )
    }

    #[test]
    fn diagnostic_payload_copies_keep_independent_run_guards_and_typed_refusal() {
        use crate::kmer::budget::collections::TryClone;
        let budget = diagnostic_result_budget(1 << 20);
        let result = owned_diagnostic_search(&budget, diagnostic_limits(4, 100, 20)).unwrap();
        let payload = result.withheld_path_diagnostics.unwrap();
        assert_eq!(payload.paths[0].oriented_sequence, "AACG");
        let used = budget.snapshot().reserved_memory_bytes;
        assert!(used > 0);
        let held = budget
            .reserve(MemoryClass::RetainedResults, (1 << 20) - used)
            .unwrap();
        let error = TryClone::try_clone(&payload).unwrap_err();
        assert!(
            error
                .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                .is_some()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, 1 << 20);
        drop(held);
        let copy = TryClone::try_clone(&payload).unwrap();
        assert_eq!(copy, payload);
        assert!(budget.snapshot().reserved_memory_bytes > used);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        drop(payload);
        assert!(budget.snapshot().reserved_memory_bytes > 0);
        drop(copy);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn diagnostic_allocation_refusal_cleans_nested_payloads_without_changing_retention_policy() {
        let (mut refused, mut succeeded) = (0, 0);
        for bytes in (1..=128).map(|n| n * 64).chain([1 << 20]) {
            let budget = diagnostic_result_budget(bytes);
            match owned_diagnostic_search(&budget, diagnostic_limits(4, 100, 20)) {
                Ok(result) => {
                    let payload = result.withheld_path_diagnostics.unwrap();
                    assert_eq!(payload.retained_paths, 1);
                    assert_eq!(payload.paths[0].oriented_sequence, "AACG");
                    succeeded += 1;
                }
                Err(error) => {
                    assert_eq!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .unwrap()
                            .resource,
                        crate::kmer::budget::ResourceKind::StageMemory(
                            MemoryClass::RetainedResults
                        )
                    );
                    refused += 1;
                }
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0, "budget {bytes}");
        }
        assert!(refused > 4 && succeeded > 0);
        let budget = diagnostic_result_budget(1);
        let result = owned_diagnostic_search(&budget, diagnostic_limits(0, 100, 20)).unwrap();
        let payload = result.withheld_path_diagnostics.unwrap();
        assert_eq!(payload.retained_paths, 0);
        assert_eq!(payload.dropped_path_cap, 1);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}

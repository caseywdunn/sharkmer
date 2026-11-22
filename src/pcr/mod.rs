// pcr/mod.rs

use crate::kmer::budget::collections::{BudgetVec, TryClone};
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use anyhow::{Context, Result};
use log::{debug, trace};
use serde::Serialize;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use crate::format::format_duration;
use crate::kmer::FilteredKmerCounts;
#[cfg(test)]
use crate::kmer::KmerCounts;

/// Log at info level with a `[gene_name]` prefix for attribution in parallel runs.
macro_rules! gene_info {
    ($gene:expr_2021, $($arg:tt)*) => {
        log::info!("[{}] {}", $gene, format!($($arg)*))
    };
}

/// Log at warn level with a `[gene_name]` prefix for attribution in parallel runs.
macro_rules! gene_warn {
    ($gene:expr_2021, $($arg:tt)*) => {
        log::warn!("[{}] {}", $gene, format!($($arg)*))
    };
}

pub mod preconfigured;

mod bubble;
mod distance;
mod graph;
mod paths;
mod primers;
mod pruning;
#[cfg(test)]
mod repeat_tests;
mod sequence;
#[cfg(test)]
mod threshold_tests;
mod uncertain;

pub use graph::compute_node_budget;
pub(crate) use uncertain::UncertainProducts;
pub(crate) mod threading;

#[allow(dead_code)]
pub(crate) mod memory;

// Constants that may require tuning

/// The multiplier for establishing when a kmer is considered to have high coverage,
/// relative to the min_count threshold. It is then used to also adjust the threshold.
const COVERAGE_MULTIPLIER: u64 = 2;

/// A multiplier for adjusting the threshold as it is applied.
const COVERAGE_STEPS: u64 = 4;

// MAX_NUM_PRIMER_KMERS is now read from params.max_primer_kmers

/// Default edit distance threshold for deduplication of sPCR products
pub const DEFAULT_DEDUP_EDIT_THRESHOLD: u32 = 10;

/// Score for a path through the assembly graph. Higher is better.
/// Used for ranking and selecting the best amplicon sequences.
#[derive(Debug, Clone, PartialEq)]
struct PathScore {
    /// Minimum kmer count along the path
    kmer_min_count: u64,
    /// Median kmer count along the path
    kmer_median_count: f64,
    /// Coverage consistency: coefficient of variation of edge counts
    /// (lower is better — a perfect single-copy path has uniform coverage)
    coverage_cv: f64,
    /// Maximum coverage ratio seen on any edge in the path
    /// (high values indicate the path traverses a repeat)
    max_coverage_ratio: f64,
    /// Number of edges with zero read support (None = threading not available)
    zero_support_edges: Option<u32>,
    /// Median unambiguous read support across path edges
    median_unambiguous_support: Option<f64>,
    /// Fraction of path edges with any read support
    edge_support_fraction: Option<f64>,
}

impl PathScore {
    /// Composite score for ranking paths. Higher is better.
    fn composite(&self) -> f64 {
        // Start with median count as the base signal (robust to outliers)
        let base = self.kmer_median_count;
        // Penalize paths with high coverage variance (likely chimeric)
        let cv_penalty = if self.coverage_cv > 1.0 {
            1.0 / self.coverage_cv
        } else {
            1.0
        };
        // Penalize paths that traverse high-coverage-ratio edges (repeats)
        let repeat_penalty = if self.max_coverage_ratio > 5.0 {
            5.0 / self.max_coverage_ratio
        } else {
            1.0
        };

        // Read support factor (only when threading data available)
        let read_support_factor = match (self.zero_support_edges, self.edge_support_fraction) {
            (Some(zero_edges), Some(support_frac)) => {
                // Penalize paths with zero-support edges: halve score per zero-support edge
                let zero_penalty = if zero_edges > 0 {
                    0.5_f64.powi(zero_edges.min(10) as i32)
                } else {
                    1.0
                };
                // Reward high support fraction (floor at 0.1 to avoid zeroing out)
                let support_bonus = support_frac.max(0.1);
                zero_penalty * support_bonus
            }
            _ => 1.0, // No threading data: neutral
        };

        base * cv_penalty * repeat_penalty * read_support_factor
    }
}

struct AssemblyRecord {
    fasta_record: sequence::SequenceRecord,
    score: PathScore,
}

/// Structured annotations joined to the final emitted FASTA identifier.
#[derive(Debug, PartialEq, Serialize)]
pub(crate) struct SequenceAnnotation {
    pub(crate) sequence_id: crate::kmer::budget::text::BudgetString,
    pub(crate) sequence_index: usize,
    pub(crate) sequence_length: usize,
    /// The existing assembly ranking heuristic, not calibrated confidence.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) assembly_rank_score: Option<f64>,
}

// Create a structure to hold a kmer representing an oligo up to 32 nucleotides long in the
// length*2 least significant bits
#[cfg(test)]
pub struct Oligo {
    pub length: usize,
    pub kmer: u64,
}

#[derive(PartialEq)]
enum PrimerDirection {
    Forward,
    Reverse,
}

// De Bruijn graph node
#[derive(Debug, Clone)]
pub struct DBNode {
    pub sub_kmer: u64,  // k-1 mer that contains overlap between kmers
    pub is_start: bool, // Contains the forward primer
    pub is_end: bool,   // Contains the reverse primer
}

#[derive(Debug, Clone)]
pub struct DBEdge {
    pub count: u64,          // Number of times this kmer was observed
    pub coverage_ratio: f64, // count / local median (set during annotation)
}

#[derive(Clone, Debug, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PCRParams {
    pub forward_seq: String,
    pub reverse_seq: String,
    #[serde(default)]
    pub min_length: usize,
    #[serde(default = "default_max_length")]
    pub max_length: usize,

    // --- Schema v2 target identification (panel YAML) ---
    // gene_name is the internal derived identifier (never in YAML; set by
    // preconfigured.rs after load, or directly by parse_pcr_primers_string).
    #[serde(skip, default)]
    pub gene_name: String,
    /// Gene being targeted (e.g. "CO1", "18S"). Required in panel YAML.
    /// Must not contain '-' or '_'.
    #[serde(default)]
    pub gene: Option<String>,
    /// Optional sub-region label (e.g. "V9", "V5-V7"). Must not contain '_'.
    #[serde(default)]
    pub region: Option<String>,
    /// Optional integer index distinguishing multiple primer pairs for the
    /// same (gene, region). Appears in output naming when set.
    #[serde(default)]
    pub index: Option<u32>,

    // --- Schema v2 target metadata ---
    /// Cellular compartment. Use INSDC /organelle vocabulary; absent = nuclear.
    #[serde(default)]
    #[allow(dead_code)]
    pub compartment: Option<String>,
    /// Gene/transcript type. Use NCBI RefSeq gene_biotype values plus
    /// rRNA_SSU, rRNA_LSU, rRNA_5S, ITS extensions.
    #[serde(default)]
    #[allow(dead_code)]
    pub gene_type: Option<String>,
    /// Approximate copy number per cell. One of: single_copy, low_copy,
    /// high_copy.
    #[serde(default)]
    #[allow(dead_code)]
    pub copy_number: Option<String>,

    // --- Deprecation ---
    #[serde(default)]
    #[allow(dead_code)]
    pub deprecated: bool,
    #[serde(default)]
    #[allow(dead_code)]
    pub deprecated_by: Option<String>,
    #[serde(default)]
    #[allow(dead_code)]
    pub deprecated_reason: Option<String>,

    #[serde(default = "default_min_count")]
    pub min_count: u32,
    #[serde(default = "default_mismatches")]
    pub mismatches: usize,
    #[serde(default = "default_trim")]
    pub trim: usize,
    /// Canonical expected amplicon length. Distinct from min_length/max_length
    /// search window. Provenance metadata; consumed by validate_panel.py.
    #[serde(default)]
    #[allow(dead_code)]
    pub expected_length: Option<usize>,
    // Provenance metadata for human readers and downstream tools.
    // Never consumed by the Rust binary; kept here because deny_unknown_fields
    // requires all YAML fields to be declared.
    #[serde(default)]
    #[allow(dead_code)]
    pub citation: String,
    #[serde(default)]
    #[allow(dead_code)]
    pub notes: String,
    #[serde(default = "default_dedup_edit_threshold")]
    pub dedup_edit_threshold: u32,
    /// Where this primer was loaded from (display only, not serialized)
    #[serde(skip)]
    pub source: String,

    // --- Runtime tuning parameters (set from CLI, not from YAML panels) ---
    #[serde(default = "default_max_dfs_states")]
    pub max_dfs_states: usize,
    #[serde(default = "default_max_paths_per_pair")]
    pub max_paths_per_pair: usize,
    #[serde(default = "default_max_node_visits")]
    pub max_node_visits: usize,
    #[serde(default = "default_max_primer_kmers")]
    pub max_primer_kmers: usize,
    #[serde(default = "default_high_coverage_ratio")]
    pub high_coverage_ratio: f64,
    #[serde(default = "default_tip_coverage_fraction")]
    pub tip_coverage_fraction: f64,
    #[serde(skip, default)]
    pub diagnose_withheld_paths: bool,
    /// Runtime-only delivery of labelled repeat hypotheses; set from the CLI.
    #[serde(skip, default)]
    pub emit_uncertain: bool,
}

impl PCRParams {
    /// All nested heap allocations retained by this configuration record.
    pub(crate) fn configuration_bytes(&self) -> Result<u64> {
        [
            Some(&self.forward_seq),
            Some(&self.reverse_seq),
            Some(&self.gene_name),
            Some(&self.citation),
            Some(&self.notes),
            Some(&self.source),
            self.gene.as_ref(),
            self.region.as_ref(),
            self.compartment.as_ref(),
            self.gene_type.as_ref(),
            self.copy_number.as_ref(),
            self.deprecated_by.as_ref(),
            self.deprecated_reason.as_ref(),
        ]
        .into_iter()
        .flatten()
        .try_fold(0_u64, |total, text| {
            total
                .checked_add(u64::try_from(text.capacity())?)
                .context("Primer configuration size overflow")
        })
    }
}

fn default_max_length() -> usize {
    10000
}
fn default_min_count() -> u32 {
    2
}
fn default_mismatches() -> usize {
    2
}
fn default_trim() -> usize {
    15
}
fn default_dedup_edit_threshold() -> u32 {
    DEFAULT_DEDUP_EDIT_THRESHOLD
}
fn default_max_dfs_states() -> usize {
    DEFAULT_MAX_DFS_STATES
}
fn default_max_paths_per_pair() -> usize {
    DEFAULT_MAX_PATHS_PER_PAIR
}
fn default_max_node_visits() -> usize {
    DEFAULT_MAX_NODE_VISITS
}
fn default_max_primer_kmers() -> usize {
    DEFAULT_MAX_NUM_PRIMER_KMERS
}
fn default_high_coverage_ratio() -> f64 {
    DEFAULT_HIGH_COVERAGE_RATIO
}
fn default_tip_coverage_fraction() -> f64 {
    DEFAULT_TIP_COVERAGE_FRACTION
}

// Default values for tuning constants (exposed as hidden CLI arguments)
pub const DEFAULT_MAX_DFS_STATES: usize = 100_000;
pub const DEFAULT_MAX_PATHS_PER_PAIR: usize = 20;
pub const DEFAULT_MAX_NODE_VISITS: usize = 2;
pub const DEFAULT_MAX_NUM_PRIMER_KMERS: usize = 40;
pub const DEFAULT_HIGH_COVERAGE_RATIO: f64 = 10.0;
pub const DEFAULT_TIP_COVERAGE_FRACTION: f64 = 0.1;
pub(crate) const WITHHELD_DIAGNOSTIC_THRESHOLD_PATH_CAP: usize = 32;
pub(crate) const WITHHELD_DIAGNOSTIC_THRESHOLD_BASE_CAP: usize = 512 * 1024;
pub(crate) const WITHHELD_DIAGNOSTIC_THRESHOLD_MARKER_CAP: usize = 4096;
pub(crate) const WITHHELD_DIAGNOSTIC_GENE_PATH_CAP: usize = 64;
pub(crate) const WITHHELD_DIAGNOSTIC_GENE_BASE_CAP: usize = 1024 * 1024;
pub(crate) const WITHHELD_DIAGNOSTIC_GENE_MARKER_CAP: usize = 8192;
pub(crate) const WITHHELD_DIAGNOSTIC_RUN_PATH_CAP: usize = 256;
// Inclusive delivery has a separate bounded sequence-only allowance alongside
// the original diagnostic reserve, so extra paths cannot displace old output.
pub(crate) const UNCERTAIN_SEQUENCE_PATH_SHARE_MULTIPLIER: usize = 3;
pub(crate) const WITHHELD_DIAGNOSTIC_RUN_BASE_CAP: usize = 4 * 1024 * 1024;
pub(crate) const WITHHELD_DIAGNOSTIC_RUN_MARKER_CAP: usize = 32768;

pub(crate) fn validate_primer_expansion(params: &PCRParams, kmer_length: usize) -> Result<()> {
    primers::validate_primer_expansion(params, kmer_length)
}

/// Result of running in silico PCR for a single gene.
pub struct PcrOutcome {
    /// FASTA records of recovered products (empty if none found)
    pub records: crate::kmer::budget::collections::BudgetVec<sequence::SequenceRecord>,
    pub(crate) sequence_annotations:
        crate::kmer::budget::collections::BudgetVec<SequenceAnnotation>,
    /// Short description of why no product was found (None if successful)
    pub failure_reason: Option<BudgetString>,
    pub threshold_diagnostics: BudgetVec<PcrThresholdDiagnostic>,
    pub(crate) uncertain_products: Option<UncertainProducts>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldDiagnosticCapSet {
    pub(crate) paths: usize,
    pub(crate) sequence_bases: usize,
    pub(crate) marker_occurrences: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldDiagnosticLimits {
    pub(crate) threshold: WithheldDiagnosticCapSet,
    pub(crate) gene: WithheldDiagnosticCapSet,
    pub(crate) run: WithheldDiagnosticCapSet,
    pub(crate) allocated_run_share_for_gene: WithheldDiagnosticCapSet,
    pub(crate) unused_share_is_not_redistributed: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldDiagnosticUsage {
    pub(crate) paths: usize,
    pub(crate) sequence_bases: usize,
    pub(crate) marker_occurrences: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum WithheldMarkerCause {
    PrePruneOmittedSelfLoopNode,
    PrePruneCyclicSccNode,
    RetainedCollisionEdge,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum WithheldMarkerIdentityKind {
    OrientedNodeSequence,
    OrientedEdgeSequence,
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldMarkerOccurrence {
    pub(crate) cause: WithheldMarkerCause,
    pub(crate) identity_kind: WithheldMarkerIdentityKind,
    pub(crate) identity_sequence: BudgetString,
    pub(crate) identity_sha256: BudgetString,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldMarkerCoveredRun {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldPathRecord {
    pub(crate) oriented_sequence: BudgetString,
    pub(crate) sequence_sha256: BudgetString,
    pub(crate) sequence_length: usize,
    pub(crate) node_visit_limit: usize,
    pub(crate) max_observed_node_visits: usize,
    pub(crate) omitted_self_loop_marker_occurrences: usize,
    pub(crate) cyclic_scc_marker_occurrences: usize,
    pub(crate) retained_collision_marker_occurrences: usize,
    pub(crate) marker_occurrences: BudgetVec<WithheldMarkerOccurrence>,
    pub(crate) marker_covered_runs: BudgetVec<WithheldMarkerCoveredRun>,
}

impl TryClone for WithheldDiagnosticCapSet {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldDiagnosticLimits {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldDiagnosticUsage {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldMarkerCause {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldMarkerIdentityKind {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldMarkerCoveredRun {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for WithheldMarkerOccurrence {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            cause: TryClone::try_clone(&self.cause)?,
            identity_kind: TryClone::try_clone(&self.identity_kind)?,
            identity_sequence: TryClone::try_clone(&self.identity_sequence)?,
            identity_sha256: TryClone::try_clone(&self.identity_sha256)?,
            start: TryClone::try_clone(&self.start)?,
            end: TryClone::try_clone(&self.end)?,
        })
    }
}
impl TryClone for WithheldPathRecord {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            oriented_sequence: TryClone::try_clone(&self.oriented_sequence)?,
            sequence_sha256: TryClone::try_clone(&self.sequence_sha256)?,
            sequence_length: TryClone::try_clone(&self.sequence_length)?,
            node_visit_limit: TryClone::try_clone(&self.node_visit_limit)?,
            max_observed_node_visits: TryClone::try_clone(&self.max_observed_node_visits)?,
            omitted_self_loop_marker_occurrences: TryClone::try_clone(
                &self.omitted_self_loop_marker_occurrences,
            )?,
            cyclic_scc_marker_occurrences: TryClone::try_clone(
                &self.cyclic_scc_marker_occurrences,
            )?,
            retained_collision_marker_occurrences: TryClone::try_clone(
                &self.retained_collision_marker_occurrences,
            )?,
            marker_occurrences: TryClone::try_clone(&self.marker_occurrences)?,
            marker_covered_runs: TryClone::try_clone(&self.marker_covered_runs)?,
        })
    }
}
impl TryClone for WithheldPathDiagnostics {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            schema_version: TryClone::try_clone(&self.schema_version)?,
            coordinate_system: TryClone::try_clone(&self.coordinate_system)?,
            purpose: TryClone::try_clone(&self.purpose)?,
            marker_span_scope: TryClone::try_clone(&self.marker_span_scope)?,
            read_support_evaluation: TryClone::try_clone(&self.read_support_evaluation)?,
            limits: TryClone::try_clone(&self.limits)?,
            gene_usage_before_threshold: TryClone::try_clone(&self.gene_usage_before_threshold)?,
            gene_usage_after_threshold: TryClone::try_clone(&self.gene_usage_after_threshold)?,
            observed_withheld_paths: TryClone::try_clone(&self.observed_withheld_paths)?,
            retained_paths: TryClone::try_clone(&self.retained_paths)?,
            observed_not_retained_paths: TryClone::try_clone(&self.observed_not_retained_paths)?,
            retained_sequence_bases: TryClone::try_clone(&self.retained_sequence_bases)?,
            retained_marker_occurrences: TryClone::try_clone(&self.retained_marker_occurrences)?,
            retained_summary_only_paths: TryClone::try_clone(&self.retained_summary_only_paths)?,
            omitted_detailed_marker_occurrences: TryClone::try_clone(
                &self.omitted_detailed_marker_occurrences,
            )?,
            marker_detail_truncated: TryClone::try_clone(&self.marker_detail_truncated)?,
            dropped_oversize_sequence: TryClone::try_clone(&self.dropped_oversize_sequence)?,
            dropped_oversize_marker_set: TryClone::try_clone(&self.dropped_oversize_marker_set)?,
            dropped_path_cap: TryClone::try_clone(&self.dropped_path_cap)?,
            dropped_sequence_base_cap: TryClone::try_clone(&self.dropped_sequence_base_cap)?,
            dropped_marker_cap: TryClone::try_clone(&self.dropped_marker_cap)?,
            node_visit_limit: TryClone::try_clone(&self.node_visit_limit)?,
            node_visit_skips: TryClone::try_clone(&self.node_visit_skips)?,
            retention_truncated: TryClone::try_clone(&self.retention_truncated)?,
            paths: TryClone::try_clone(&self.paths)?,
        })
    }
}
impl TryClone for PcrThresholdDiagnostic {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            threshold: TryClone::try_clone(&self.threshold)?,
            connectivity_found: TryClone::try_clone(&self.connectivity_found)?,
            node_budget_reached: TryClone::try_clone(&self.node_budget_reached)?,
            scc_evaluated: TryClone::try_clone(&self.scc_evaluated)?,
            path_search_evaluated: TryClone::try_clone(&self.path_search_evaluated)?,
            omitted_self_loop_nodes: TryClone::try_clone(&self.omitted_self_loop_nodes)?,
            pre_prune_scc_nodes: TryClone::try_clone(&self.pre_prune_scc_nodes)?,
            retained_collision_edges: TryClone::try_clone(&self.retained_collision_edges)?,
            dfs_limit_reached: TryClone::try_clone(&self.dfs_limit_reached)?,
            path_quota_reached: TryClone::try_clone(&self.path_quota_reached)?,
            end_below_min_length: TryClone::try_clone(&self.end_below_min_length)?,
            max_length_reached: TryClone::try_clone(&self.max_length_reached)?,
            completed_candidate_paths: TryClone::try_clone(&self.completed_candidate_paths)?,
            eligible_candidate_paths: TryClone::try_clone(&self.eligible_candidate_paths)?,
            withheld_candidate_paths: TryClone::try_clone(&self.withheld_candidate_paths)?,
            withheld_by_self_loop_node_paths: TryClone::try_clone(
                &self.withheld_by_self_loop_node_paths,
            )?,
            withheld_by_scc_node_paths: TryClone::try_clone(&self.withheld_by_scc_node_paths)?,
            withheld_by_collision_edge_paths: TryClone::try_clone(
                &self.withheld_by_collision_edge_paths,
            )?,
            repeat_tainted_out_of_range_paths: TryClone::try_clone(
                &self.repeat_tainted_out_of_range_paths,
            )?,
            generated_products: TryClone::try_clone(&self.generated_products)?,
            failure_reason: TryClone::try_clone(&self.failure_reason)?,
            withheld_path_diagnostics: TryClone::try_clone(&self.withheld_path_diagnostics)?,
        })
    }
}
fn is_zero_usize(value: &usize) -> bool {
    *value == 0
}

fn is_false(value: &bool) -> bool {
    !*value
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub(crate) struct WithheldPathDiagnostics {
    pub(crate) schema_version: u32,
    pub(crate) coordinate_system: &'static str,
    pub(crate) purpose: &'static str,
    pub(crate) marker_span_scope: &'static str,
    pub(crate) read_support_evaluation: &'static str,
    pub(crate) limits: WithheldDiagnosticLimits,
    pub(crate) gene_usage_before_threshold: WithheldDiagnosticUsage,
    pub(crate) gene_usage_after_threshold: WithheldDiagnosticUsage,
    pub(crate) observed_withheld_paths: usize,
    pub(crate) retained_paths: usize,
    pub(crate) observed_not_retained_paths: usize,
    pub(crate) retained_sequence_bases: usize,
    pub(crate) retained_marker_occurrences: usize,
    #[serde(skip_serializing_if = "is_zero_usize")]
    pub(crate) retained_summary_only_paths: usize,
    #[serde(skip_serializing_if = "is_zero_usize")]
    pub(crate) omitted_detailed_marker_occurrences: usize,
    #[serde(skip_serializing_if = "is_false")]
    pub(crate) marker_detail_truncated: bool,
    pub(crate) dropped_oversize_sequence: usize,
    pub(crate) dropped_oversize_marker_set: usize,
    pub(crate) dropped_path_cap: usize,
    pub(crate) dropped_sequence_base_cap: usize,
    pub(crate) dropped_marker_cap: usize,
    pub(crate) node_visit_limit: usize,
    pub(crate) node_visit_skips: usize,
    pub(crate) retention_truncated: bool,
    pub(crate) paths: BudgetVec<WithheldPathRecord>,
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, Default, PartialEq, Eq, Serialize)]
pub struct PcrThresholdDiagnostic {
    pub threshold: u64,
    pub connectivity_found: bool,
    pub node_budget_reached: bool,
    pub scc_evaluated: bool,
    pub path_search_evaluated: bool,
    pub omitted_self_loop_nodes: usize,
    pub pre_prune_scc_nodes: Option<usize>,
    pub retained_collision_edges: usize,
    pub dfs_limit_reached: bool,
    pub path_quota_reached: bool,
    pub end_below_min_length: bool,
    pub max_length_reached: bool,
    pub completed_candidate_paths: usize,
    pub eligible_candidate_paths: usize,
    pub withheld_candidate_paths: usize,
    pub withheld_by_self_loop_node_paths: usize,
    pub withheld_by_scc_node_paths: usize,
    pub withheld_by_collision_edge_paths: usize,
    pub repeat_tainted_out_of_range_paths: usize,
    pub generated_products: usize,
    pub failure_reason: Option<BudgetString>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) withheld_path_diagnostics: Option<WithheldPathDiagnostics>,
}

fn partition_cap(total: usize, gene_index: usize, gene_count: usize) -> usize {
    if gene_count == 0 {
        return 0;
    }
    total / gene_count + usize::from(gene_index < total % gene_count)
}

pub(crate) fn withheld_diagnostic_limits_for_gene(
    gene_index: usize,
    gene_count: usize,
) -> WithheldDiagnosticLimits {
    let run = WithheldDiagnosticCapSet {
        paths: WITHHELD_DIAGNOSTIC_RUN_PATH_CAP,
        sequence_bases: WITHHELD_DIAGNOSTIC_RUN_BASE_CAP,
        marker_occurrences: WITHHELD_DIAGNOSTIC_RUN_MARKER_CAP,
    };
    let gene = WithheldDiagnosticCapSet {
        paths: WITHHELD_DIAGNOSTIC_GENE_PATH_CAP,
        sequence_bases: WITHHELD_DIAGNOSTIC_GENE_BASE_CAP,
        marker_occurrences: WITHHELD_DIAGNOSTIC_GENE_MARKER_CAP,
    };
    WithheldDiagnosticLimits {
        threshold: WithheldDiagnosticCapSet {
            paths: WITHHELD_DIAGNOSTIC_THRESHOLD_PATH_CAP,
            sequence_bases: WITHHELD_DIAGNOSTIC_THRESHOLD_BASE_CAP,
            marker_occurrences: WITHHELD_DIAGNOSTIC_THRESHOLD_MARKER_CAP,
        },
        gene,
        run,
        allocated_run_share_for_gene: WithheldDiagnosticCapSet {
            paths: partition_cap(run.paths, gene_index, gene_count).min(gene.paths),
            sequence_bases: partition_cap(run.sequence_bases, gene_index, gene_count)
                .min(gene.sequence_bases),
            marker_occurrences: partition_cap(run.marker_occurrences, gene_index, gene_count)
                .min(gene.marker_occurrences),
        },
        unused_share_is_not_redistributed: true,
    }
}

pub(crate) fn inclusive_sequence_limits(
    mut limits: WithheldDiagnosticLimits,
) -> WithheldDiagnosticLimits {
    let extra_share = limits
        .allocated_run_share_for_gene
        .paths
        .saturating_mul(UNCERTAIN_SEQUENCE_PATH_SHARE_MULTIPLIER)
        .min(limits.gene.paths);
    limits.threshold.paths = limits.threshold.paths.saturating_mul(2);
    limits.gene.paths = limits.gene.paths.saturating_mul(2);
    limits.run.paths = limits
        .run
        .paths
        .saturating_mul(1 + UNCERTAIN_SEQUENCE_PATH_SHARE_MULTIPLIER);
    limits.allocated_run_share_for_gene.paths = limits
        .allocated_run_share_for_gene
        .paths
        .saturating_add(extra_share);
    limits.threshold.sequence_bases = limits.threshold.sequence_bases.saturating_mul(2);
    limits.gene.sequence_bases = limits.gene.sequence_bases.saturating_mul(2);
    limits.run.sequence_bases = limits.run.sequence_bases.saturating_mul(2);
    limits.allocated_run_share_for_gene.sequence_bases = limits
        .allocated_run_share_for_gene
        .sequence_bases
        .saturating_mul(2);
    limits
}

/// Validate a primer pair and return a list of (error, suggestion) pairs.
/// An empty list means the primer is valid.
pub(crate) fn validate_pcr_params(
    params: &PCRParams,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<(BudgetString, BudgetString)>> {
    let mut errors = BudgetVec::in_class(budget, MemoryClass::RetainedResults);
    macro_rules! report {
        ($message:expr, $suggestion:expr) => {
            errors.push((
                BudgetString::format($message, budget, MemoryClass::RetainedResults)?,
                BudgetString::from_str($suggestion, budget, MemoryClass::RetainedResults)?,
            ))?
        };
    }
    for (direction, sequence) in [
        ("Forward", &params.forward_seq),
        ("Reverse", &params.reverse_seq),
    ] {
        if sequence.len() < 2 {
            report!(
                format_args!("{direction} primer sequence is too short: '{sequence}'"),
                "Primer sequences must be at least 2 bases"
            );
        }
    }
    struct Invalid<'a>(&'a str);
    impl std::fmt::Display for Invalid<'_> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let mut first = true;
            for character in self.0.chars().filter(|c| !primers::is_valid_nucleotide(*c)) {
                if !first {
                    f.write_str(", ")?;
                }
                write!(f, "{character}")?;
                first = false;
            }
            Ok(())
        }
    }
    for (direction, sequence) in [
        ("forward", &params.forward_seq),
        ("reverse", &params.reverse_seq),
    ] {
        if sequence.len() >= 2 && sequence.chars().any(|c| !primers::is_valid_nucleotide(c)) {
            report!(
                format_args!(
                    "Invalid nucleotide(s) {} in {direction} primer {sequence}",
                    Invalid(sequence)
                ),
                "Valid characters: A C G T R Y W S M K B D H V N"
            );
        }
    }
    if params.min_length > params.max_length {
        report!(
            format_args!(
                "min-length ({}) is greater than max-length ({})",
                params.min_length, params.max_length
            ),
            "Swap the values or adjust the range"
        );
    }
    if params.min_count < 2 {
        report!(
            format_args!("min-count is {}, must be at least 2", params.min_count),
            "Set min-count to at least 2"
        );
    }
    if params.max_length == 0 {
        report!(
            format_args!("max-length is 0"),
            "Set max-length to a positive value"
        );
    }
    if params.gene_name.is_empty() {
        report!(
            format_args!("Gene name is empty"),
            "Provide a unique name for the primer pair via the 'name' field"
        );
    }
    if params.forward_seq == params.reverse_seq && params.forward_seq.len() >= 2 {
        report!(
            format_args!(
                "Forward and reverse primers are identical: {}",
                params.forward_seq
            ),
            "Check that forward and reverse sequences are not swapped"
        );
    }
    Ok(errors)
}

/// Compute a sequence of coverage thresholds for graph extension, stepping down
/// from a high threshold (derived from observed primer coverage) to `min_count`.
fn compute_coverage_thresholds(primer_count: u64, min_count: u64) -> Vec<u64> {
    let coverage_high_threshold = primer_count / COVERAGE_MULTIPLIER;
    let mut thresholds: Vec<u64> = Vec::new();

    if coverage_high_threshold <= min_count {
        thresholds.push(min_count);
    } else {
        let step_size = (coverage_high_threshold - min_count) / (COVERAGE_STEPS - 1);

        for i in 0..COVERAGE_STEPS {
            thresholds.push(coverage_high_threshold.saturating_sub(i * step_size));
        }

        // Make sure the last element is min_count, even if there are rounding errors
        *thresholds
            .last_mut()
            .expect("coverage_thresholds is non-empty") = min_count;
    }

    // Deduplicate consecutive identical thresholds (e.g. [4, 4, 4, 2] -> [4, 2]).
    // When step_size rounds to 0 the loop above produces repeats, which would
    // make the threshold sweep redo identical work.
    thresholds.dedup();

    thresholds
}

#[cfg(test)]
fn earlier_threshold_limit_summary(diagnostics: &[PcrThresholdDiagnostic]) -> Option<BudgetString> {
    earlier_threshold_limit_summary_budgeted(diagnostics, None).unwrap()
}

fn earlier_threshold_limit_summary_budgeted(
    diagnostics: &[PcrThresholdDiagnostic],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Option<BudgetString>> {
    let Some(last) = diagnostics.last() else {
        return Ok(None);
    };
    let mut summary = BudgetString::new(budget, MemoryClass::RetainedResults);
    for diagnostic in &diagnostics[..diagnostics.len().saturating_sub(1)] {
        let repeat = diagnostic.withheld_candidate_paths > 0
            || diagnostic.repeat_tainted_out_of_range_paths > 0;
        let last_repeat =
            last.withheld_candidate_paths > 0 || last.repeat_tainted_out_of_range_paths > 0;
        let causes = [
            (
                "node budget reached",
                diagnostic.node_budget_reached && !last.node_budget_reached,
            ),
            (
                "DFS state limit reached",
                diagnostic.dfs_limit_reached && !last.dfs_limit_reached,
            ),
            (
                "path quota reached",
                diagnostic.path_quota_reached && !last.path_quota_reached,
            ),
            ("candidate-local repeat uncertainty", repeat && !last_repeat),
        ];
        if !causes.iter().any(|(_, present)| *present) {
            continue;
        }
        if !summary.is_empty() {
            summary.push_str("; ")?;
        }
        let prefix = BudgetString::format(
            format_args!("threshold {}: ", diagnostic.threshold),
            budget,
            MemoryClass::RetainedResults,
        )?;
        summary.push_str(&prefix)?;
        let mut first = true;
        for (cause, present) in causes {
            if present {
                if !first {
                    summary.push_str(", ")?;
                }
                summary.push_str(cause)?;
                first = false;
            }
        }
    }
    Ok((!summary.is_empty()).then_some(summary))
}

struct ThresholdEvaluation {
    records: BudgetVec<AssemblyRecord>,
    failure_reason: Option<BudgetString>,
    diagnostic: PcrThresholdDiagnostic,
}

#[allow(clippy::too_many_arguments)]
fn evaluate_threshold_graph(
    mut pruned_graph: memory::AssemblyGraph,
    extension_repeat_markers: graph::ExtensionRepeatMarkers,
    min_count: u64,
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    dump_graph: bool,
    output_directory: &str,
    reads: Option<crate::io::ThreadingSource<'_>>,
    diagnostic_budget: Option<&mut paths::WithheldDiagnosticBudget>,
    result_budget: Option<&Arc<ResourceBudget>>,
) -> Result<ThresholdEvaluation> {
    let repeat_markers = graph::repeat_markers_before_pruning_budgeted(
        &pruned_graph,
        extension_repeat_markers,
        pruned_graph.budget(),
    )?;
    let omitted_self_loop_nodes = repeat_markers.omitted_self_loop_sub_kmers.len();
    let pre_prune_scc_nodes = repeat_markers.cyclic_scc_sub_kmers.len();
    let retained_collision_edges = repeat_markers.retained_collision_edges.len();
    let prune_start = std::time::Instant::now();
    gene_info!(params.gene_name, "Pruning the assembly graph...");

    pruned_graph.prune_tips(&kmer_counts.get_k(), params.tip_coverage_fraction)?;
    pruned_graph.prune_reachability()?;
    pruned_graph.annotate_coverage()?;
    gene_info!(
        params.gene_name,
        "Done. Time to prune graph: {}",
        format_duration(prune_start.elapsed())
    );

    if dump_graph || log::log_enabled!(log::Level::Trace) {
        let dot_string = write_annotated_dot(&pruned_graph, kmer_counts);
        let file_name = format!(
            "{}{}_{}_{}.dot",
            output_directory, sample_name, params.gene_name, min_count
        );
        trace!("Writing dot file {}", file_name);
        let mut file = File::create(&file_name).context("Unable to create dot file")?;
        file.write_all(dot_string.as_bytes())
            .context("Unable to write dot file")?;
    }

    let threading_annotations = if let Some(reads) = reads {
        if !reads.is_empty() {
            let threading_start = std::time::Instant::now();
            gene_info!(
                params.gene_name,
                "Threading reads through assembly graph..."
            );
            let annotations = match reads {
                crate::io::ThreadingSource::Memory(records) => {
                    if records
                        .iter()
                        .any(|read| read.mate != crate::io::Mate::Unpaired)
                    {
                        threading::thread_reads_paired_budgeted(
                            &pruned_graph,
                            records,
                            kmer_counts.get_k(),
                            pruned_graph.budget(),
                        )?
                    } else {
                        threading::thread_reads_budgeted(
                            &pruned_graph,
                            records,
                            kmer_counts.get_k(),
                            pruned_graph.budget(),
                        )?
                    }
                }
                crate::io::ThreadingSource::Replay(replay) => threading::thread_replay_budgeted(
                    &pruned_graph,
                    replay,
                    kmer_counts.get_k(),
                    pruned_graph.budget(),
                )?,
            };
            let supported_edges = annotations
                .edge_support
                .values()
                .filter(|support| support.read_support_total > 0)
                .count();
            gene_info!(
                params.gene_name,
                "Threading: {}/{} edges have read support, {} branch links, {} paired links. Time: {}",
                supported_edges,
                pruned_graph.edge_count(),
                annotations.branch_links.len(),
                annotations.paired_links.len(),
                format_duration(threading_start.elapsed())
            );
            Some(annotations)
        } else {
            None
        }
    } else {
        None
    };

    let path_start = std::time::Instant::now();
    gene_info!(
        params.gene_name,
        "Traversing the assembly graph to find paths from forward to reverse primers..."
    );

    let edge_preferences = threading_annotations
        .as_ref()
        .map(|annotations| {
            bubble::resolve_bubbles_budgeted(&pruned_graph, annotations, pruned_graph.budget())
        })
        .transpose()?;
    let mut path_search = paths::search_assembly_paths_budgeted(
        &pruned_graph,
        kmer_counts,
        params,
        edge_preferences.as_ref(),
        &repeat_markers,
        diagnostic_budget,
        pruned_graph.budget(),
    )?;

    gene_info!(
        params.gene_name,
        "Found {} paths. Time: {}",
        path_search.paths.len(),
        format_duration(path_start.elapsed())
    );

    if path_search.dfs_limit_reached {
        gene_warn!(
            params.gene_name,
            "Path search reached the DFS state limit of {} at threshold {}.",
            params.max_dfs_states,
            min_count
        );
    }
    if path_search.path_limit_reached {
        gene_warn!(
            params.gene_name,
            "Path search reached the per-primer-pair path limit of {} at threshold {}.",
            params.max_paths_per_pair,
            min_count
        );
    }

    if path_search.withheld_candidate_path_count > 0 {
        gene_warn!(
            params.gene_name,
            "Withholding {} complete candidate path(s) because repeat length is unresolved at threshold {}.",
            path_search.withheld_candidate_path_count,
            min_count
        );
    }

    let mut diagnostic = PcrThresholdDiagnostic {
        threshold: min_count,
        connectivity_found: true,
        node_budget_reached: false,
        scc_evaluated: true,
        path_search_evaluated: true,
        omitted_self_loop_nodes,
        pre_prune_scc_nodes: Some(pre_prune_scc_nodes),
        retained_collision_edges,
        dfs_limit_reached: path_search.dfs_limit_reached,
        path_quota_reached: path_search.path_limit_reached,
        end_below_min_length: path_search.end_below_min_length,
        max_length_reached: path_search.max_length_reached,
        completed_candidate_paths: path_search.completed_candidate_path_count,
        eligible_candidate_paths: path_search.eligible_candidate_path_count,
        withheld_candidate_paths: path_search.withheld_candidate_path_count,
        withheld_by_self_loop_node_paths: path_search.withheld_by_self_loop_node_count,
        withheld_by_scc_node_paths: path_search.withheld_by_scc_node_count,
        withheld_by_collision_edge_paths: path_search.withheld_by_collision_edge_count,
        repeat_tainted_out_of_range_paths: path_search.repeat_tainted_out_of_range_path_count,
        generated_products: 0,
        failure_reason: None,
        withheld_path_diagnostics: path_search.withheld_path_diagnostics.take(),
    };

    if path_search.paths.is_empty() {
        let failure_reason =
            if path_search.withheld_candidate_path_count > 0 && path_search.dfs_limit_reached {
                "DFS state limit reached; completed candidate paths had unresolved repeat length"
            } else if path_search.withheld_candidate_path_count > 0 {
                "repeat length unresolved for all complete candidate paths"
            } else if path_search.repeat_tainted_out_of_range_path_count > 0 {
                "repeat length unresolved on an out-of-range start-to-end path"
            } else if path_search.dfs_limit_reached {
                "DFS state limit reached before a valid amplicon was found"
            } else if path_search.path_limit_reached {
                "path limit reached before a valid amplicon was found"
            } else if path_search.end_below_min_length || path_search.max_length_reached {
                "connectivity found but no path satisfied the requested length range"
            } else {
                "connectivity found but no start-to-end path remained after pruning"
            };
        diagnostic.failure_reason = Some(BudgetString::from_str(
            failure_reason,
            result_budget,
            MemoryClass::RetainedResults,
        )?);
        return Ok(ThresholdEvaluation {
            records: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            failure_reason: Some(BudgetString::from_str(
                failure_reason,
                result_budget,
                MemoryClass::RetainedResults,
            )?),
            diagnostic,
        });
    }

    let (records, _) = paths::generate_sequences_from_paths_budgeted(
        &pruned_graph,
        path_search.paths,
        kmer_counts,
        sample_name,
        params,
        0,
        threading_annotations.as_ref(),
        pruned_graph.budget(),
        result_budget,
    )?;

    let failure_reason = if records.is_empty() {
        Some(BudgetString::from_str(
            "connectivity found but no path produced a valid-length amplicon",
            result_budget,
            MemoryClass::RetainedResults,
        )?)
    } else {
        None
    };
    diagnostic.generated_products = records.len();
    diagnostic.failure_reason = TryClone::try_clone(&failure_reason)?;
    Ok(ThresholdEvaluation {
        records,
        failure_reason,
        diagnostic,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn do_pcr(
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    dump_graph: bool,
    output_directory: &str,
    reads: Option<crate::io::ThreadingSource<'_>>,
    max_num_nodes: usize,
    withheld_diagnostic_limits: WithheldDiagnosticLimits,
) -> Result<PcrOutcome> {
    do_pcr_with_graph_budget(
        kmer_counts,
        sample_name,
        params,
        dump_graph,
        output_directory,
        reads,
        max_num_nodes,
        withheld_diagnostic_limits,
        None,
        None,
    )
}

/// Execute the same assembly policy with optional graph-buffer ownership.
#[allow(clippy::too_many_arguments)]
pub(crate) fn do_pcr_with_graph_budget(
    kmer_counts: &FilteredKmerCounts,
    sample_name: &str,
    params: &PCRParams,
    dump_graph: bool,
    output_directory: &str,
    reads: Option<crate::io::ThreadingSource<'_>>,
    max_num_nodes: usize,
    withheld_diagnostic_limits: WithheldDiagnosticLimits,
    graph_budget: Option<&std::sync::Arc<crate::kmer::budget::ResourceBudget>>,
    result_budget: Option<&Arc<ResourceBudget>>,
) -> Result<PcrOutcome> {
    gene_info!(params.gene_name, "Running PCR");
    let emission_limits = if params.emit_uncertain {
        inclusive_sequence_limits(withheld_diagnostic_limits)
    } else {
        withheld_diagnostic_limits
    };

    if params.max_length < kmer_counts.get_k() {
        return Ok(PcrOutcome {
            records: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            sequence_annotations: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            failure_reason: Some(BudgetString::format(
                format_args!(
                    "max-length {} is shorter than k-mer length {}",
                    params.max_length,
                    kmer_counts.get_k()
                ),
                result_budget,
                MemoryClass::RetainedResults,
            )?),
            threshold_diagnostics: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            uncertain_products: params
                .emit_uncertain
                .then(|| UncertainProducts::empty(emission_limits, result_budget)),
        });
    }

    gene_info!(params.gene_name, "Preprocessing primers");
    let (forward_primer_kmers, reverse_primer_kmers) =
        primers::get_primer_kmers(params, kmer_counts, graph_budget)?;

    let fwd_missing = forward_primer_kmers.is_empty();
    let rev_missing = reverse_primer_kmers.is_empty();
    if fwd_missing || rev_missing {
        let which = match (fwd_missing, rev_missing) {
            (true, true) => "forward and reverse primers",
            (true, false) => "forward primer",
            (false, true) => "reverse primer",
            _ => unreachable!(),
        };
        gene_info!(
            params.gene_name,
            "Binding sites were not found for the {}. Abandoning PCR.",
            which
        );
        gene_info!(
            params.gene_name,
            "Suggested actions: optimize primer sequence, or increase the number of reads."
        );
        return Ok(PcrOutcome {
            records: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            sequence_annotations: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            failure_reason: Some(BudgetString::format(
                format_args!("{} not found", which),
                result_budget,
                MemoryClass::RetainedResults,
            )?),
            threshold_diagnostics: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            uncertain_products: params
                .emit_uncertain
                .then(|| UncertainProducts::empty(emission_limits, result_budget)),
        });
    }

    // Log forward primer kmers for diagnostics
    let mut sorted_forward = memory::buffers::BudgetVec::collect(
        forward_primer_kmers.iter().map(|(&k, &v)| (k, v)),
        graph_budget,
    )?;
    sorted_forward.sort_unstable();
    for (kmer, count) in sorted_forward.iter() {
        gene_info!(
            params.gene_name,
            "Forward primer kmer {} (count {})",
            crate::kmer::kmer_to_seq(kmer, &kmer_counts.get_k()),
            count
        );
    }

    drop(sorted_forward);

    // Log reverse primer kmers for diagnostics
    let mut sorted_reverse = memory::buffers::BudgetVec::collect(
        reverse_primer_kmers.iter().map(|(&k, &v)| (k, v)),
        graph_budget,
    )?;
    sorted_reverse.sort_unstable();
    for (kmer, count) in sorted_reverse.iter() {
        gene_info!(
            params.gene_name,
            "Reverse primer kmer {} (count {})",
            crate::kmer::kmer_to_seq(kmer, &kmer_counts.get_k()),
            count
        );
    }

    drop(sorted_reverse);

    // Build a single graph seeded with all forward and reverse primer kmers
    gene_info!(
        params.gene_name,
        "Creating graph, seeding with {} forward and {} reverse primer kmer nodes...",
        forward_primer_kmers.len(),
        reverse_primer_kmers.len()
    );
    let (seed_graph, node_lookup) = graph::create_seed_graph(
        &forward_primer_kmers,
        &reverse_primer_kmers,
        kmer_counts,
        graph_budget,
    )?;

    debug!(
        "There are {} start nodes",
        seed_graph
            .node_weights()
            .filter(|node| node.is_start)
            .count()
    );
    debug!(
        "There are {} end nodes",
        seed_graph.node_weights().filter(|node| node.is_end).count()
    );

    for node in seed_graph.node_indices() {
        trace!("Node {}:", node.index());
        trace!(
            "  sub_kmer: {}",
            crate::kmer::kmer_to_seq(&seed_graph[node].sub_kmer, &(kmer_counts.get_k() - 1))
        );
        trace!("  is_start: {}", seed_graph[node].is_start);
        trace!("  is_end: {}", seed_graph[node].is_end);
    }

    let max_forward_count = forward_primer_kmers.get_max_count();
    let max_reverse_count = reverse_primer_kmers.get_max_count();
    let max_primer_count = max_forward_count.min(max_reverse_count);

    let median_forward_count = forward_primer_kmers.get_median_count()?;
    let median_reverse_count = reverse_primer_kmers.get_median_count()?;
    let median_primer_count = median_forward_count.min(median_reverse_count);

    gene_info!(
        params.gene_name,
        "Observed primer coverage: median {}, max fwd {}, max rev {}. User specified min-count is {}",
        median_primer_count,
        max_forward_count,
        max_reverse_count,
        params.min_count
    );

    drop(forward_primer_kmers);
    drop(reverse_primer_kmers);

    // Graph extension thresholds: derived from max primer kmer count.
    // Starting high ensures the graph is built from confident kmers first,
    // then steps down to include lower-coverage regions.
    let coverage_thresholds =
        compute_coverage_thresholds(max_primer_count, u64::from(params.min_count));
    debug!("Minimum kmer counts to attempt: {:?}", coverage_thresholds);

    let mut assembly_records_all = BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);
    let mut failure_reason = Some(BudgetString::from_str(
        "no path found",
        result_budget,
        MemoryClass::RetainedResults,
    )?);
    let mut threshold_diagnostics =
        BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);
    let mut withheld_diagnostic_budget = (params.diagnose_withheld_paths || params.emit_uncertain)
        .then(|| {
            if params.emit_uncertain {
                paths::WithheldDiagnosticBudget::for_emission(withheld_diagnostic_limits)
            } else {
                paths::WithheldDiagnosticBudget::new(withheld_diagnostic_limits)
            }
        });
    if let Some(budget) = withheld_diagnostic_budget.as_mut() {
        budget.result_budget_from(result_budget);
    }

    // Coverage threshold sweep: at each min_count, clone the seed graph fresh
    // and extend at that threshold. Pruning + path-finding happens after each
    // step. If a product is found, stop. This matches v2's approach: each
    // threshold step is independent — no incremental state across steps.
    gene_info!(
        params.gene_name,
        "Extending graph with thresholds {:?} (global budget {}, repeat edge count floor {})",
        coverage_thresholds,
        max_num_nodes,
        max_primer_count
    );

    let extend_start = std::time::Instant::now();
    drop(node_lookup);
    'threshold_loop: for (step_idx, min_count) in coverage_thresholds.iter().enumerate() {
        gene_info!(
            params.gene_name,
            "Threshold step {}/{} (min_count={})",
            step_idx + 1,
            coverage_thresholds.len(),
            min_count
        );

        // Start fresh from the seed graph at each threshold step
        let fresh_graph = seed_graph.try_clone()?;
        let fresh_lookup = memory::lookup::NodeLookup::from_graph(&fresh_graph)?;

        // Unified bidirectional extension: processes forward and reverse
        // seeds in one pass with an interleaved frontier. Symmetric — swapping
        // forward/reverse primer labels doesn't change behavior.
        let (final_graph, _final_lookup, found, extension_repeat_markers, node_budget_reached) =
            graph::extend_graph(
                fresh_graph,
                fresh_lookup,
                kmer_counts,
                min_count,
                params,
                max_primer_count,
                max_num_nodes,
            )?;
        // Extension lookup is not used by pruning, threading or path search.
        // Release its owned table before their transient resize peaks.
        drop(_final_lookup);
        if !found {
            let threshold_failure = if node_budget_reached {
                "node budget exceeded"
            } else {
                "no start-to-end connectivity established at this threshold"
            };
            failure_reason = Some(BudgetString::from_str(
                threshold_failure,
                result_budget,
                MemoryClass::RetainedResults,
            )?);
            threshold_diagnostics.push(PcrThresholdDiagnostic {
                threshold: *min_count,
                connectivity_found: false,
                node_budget_reached,
                omitted_self_loop_nodes: extension_repeat_markers.omitted_self_loop_sub_kmers.len(),
                retained_collision_edges: extension_repeat_markers.retained_collision_edges.len(),
                failure_reason: Some(BudgetString::from_str(
                    threshold_failure,
                    result_budget,
                    MemoryClass::RetainedResults,
                )?),
                withheld_path_diagnostics: withheld_diagnostic_budget
                    .as_ref()
                    .map(|budget| budget.empty_threshold_diagnostics(params.max_node_visits)),
                ..PcrThresholdDiagnostic::default()
            })?;
            continue;
        }
        let mut evaluation = evaluate_threshold_graph(
            final_graph,
            extension_repeat_markers,
            *min_count,
            kmer_counts,
            sample_name,
            params,
            dump_graph,
            output_directory,
            reads,
            withheld_diagnostic_budget.as_mut(),
            result_budget,
        )?;
        evaluation.diagnostic.node_budget_reached = node_budget_reached;

        if !evaluation.records.is_empty() {
            gene_info!(
                params.gene_name,
                "Obtained {} PCR product(s) at threshold {}.",
                evaluation.records.len(),
                min_count
            );
            threshold_diagnostics.push(evaluation.diagnostic)?;
            for record in evaluation.records {
                assembly_records_all.push(record)?;
            }
            failure_reason = None;
            break 'threshold_loop;
        }

        let evaluation_failure = match evaluation.failure_reason {
            Some(reason) => reason,
            None => BudgetString::from_str(
                "connectivity found but no valid amplicon was produced",
                result_budget,
                MemoryClass::RetainedResults,
            )?,
        };
        failure_reason = Some(if node_budget_reached {
            BudgetString::format(
                format_args!("node budget exceeded; {}", evaluation_failure),
                result_budget,
                MemoryClass::RetainedResults,
            )?
        } else {
            evaluation_failure
        });
        evaluation.diagnostic.failure_reason = TryClone::try_clone(&failure_reason)?;
        threshold_diagnostics.push(evaluation.diagnostic)?;
    }

    gene_info!(
        params.gene_name,
        "Done. Time to extend and evaluate graphs: {}",
        format_duration(extend_start.elapsed())
    );

    if assembly_records_all.is_empty()
        && let Some(summary) =
            earlier_threshold_limit_summary_budgeted(&threshold_diagnostics, result_budget)?
    {
        failure_reason = Some(BudgetString::format(
            format_args!(
                "{}; earlier threshold outcomes: {}",
                failure_reason.as_deref().unwrap_or("no path found"),
                summary
            ),
            result_budget,
            MemoryClass::RetainedResults,
        )?);
    }

    debug!(
        "      - The maximum count of a forward kmer is {} and of a reverse kmer is {}. Large differences in value can indicate non-specific binding of one of the primers.",
        max_forward_count, max_reverse_count
    );

    let count_threshold: u64 = 5;
    if (max_forward_count < count_threshold) || (max_reverse_count < count_threshold) {
        gene_info!(
            params.gene_name,
            "Primer kmer counts are low, in this case less than {}. Consider increasing the number of reads.",
            count_threshold
        );
    }

    gene_info!(params.gene_name, "Done.");

    if assembly_records_all.is_empty() {
        let reported_failure = failure_reason.as_deref().unwrap_or("no path found");
        gene_info!(
            params.gene_name,
            "No valid PCR product was recovered: {}. Abandoning PCR.",
            reported_failure
        );
        gene_info!(params.gene_name, "Suggested actions:");
        gene_info!(
            params.gene_name,
            "  - The max-length for the PCR product of {} may be too short. Consider increasing it.",
            params.max_length
        );
        gene_info!(
            params.gene_name,
            "  - The primers may have non-specific binding and are not close enough to generate a product. Consider increasing the primer trim length from the default to create a more specific primer."
        );

        let uncertain_products = params
            .emit_uncertain
            .then(|| {
                uncertain::collect(
                    sample_name,
                    &params.gene_name,
                    params.diagnose_withheld_paths,
                    emission_limits,
                    &[],
                    &mut threshold_diagnostics,
                    result_budget,
                )
            })
            .transpose()?;
        return Ok(PcrOutcome {
            records: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            sequence_annotations: BudgetVec::in_class(result_budget, MemoryClass::RetainedResults),
            failure_reason,
            threshold_diagnostics,
            uncertain_products,
        });
    }

    let selected = paths::sort_and_deduplicate_budgeted(assembly_records_all, params)?;
    let mut sequence_annotations = BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);
    let mut records = BudgetVec::in_class(result_budget, MemoryClass::RetainedResults);
    // Header replacement moves the original sequence instead of cloning it.
    for (index, mut selected) in selected.into_iter().enumerate() {
        let id = BudgetString::format(
            format_args!("{}_{}_{}", sample_name, params.gene_name, index),
            result_budget,
            MemoryClass::RetainedResults,
        )?;
        sequence_annotations.push(SequenceAnnotation {
            sequence_id: BudgetString::from_str(&id, result_budget, MemoryClass::RetainedResults)?,
            sequence_index: index,
            sequence_length: selected.fasta_record.seq().len(),
            assembly_rank_score: Some(selected.score.composite()),
        })?;
        let desc = if let Some(description) = selected.fasta_record.desc() {
            let mut updated = BudgetString::new(result_budget, MemoryClass::RetainedResults);
            for (field_index, field) in description.split_whitespace().enumerate() {
                if field_index > 0 {
                    updated.push_str(" ")?;
                }
                if field.starts_with("sequence_index=") {
                    let index_field = BudgetString::format(
                        format_args!("sequence_index={index}"),
                        result_budget,
                        MemoryClass::RetainedResults,
                    )?;
                    updated.push_str(&index_field)?;
                } else {
                    updated.push_str(field)?;
                }
            }
            Some(updated)
        } else {
            None
        };
        selected.fasta_record.replace_header(id, desc);
        records.push(selected.fasta_record)?;
    }

    let uncertain_products = params
        .emit_uncertain
        .then(|| {
            uncertain::collect(
                sample_name,
                &params.gene_name,
                params.diagnose_withheld_paths,
                emission_limits,
                &records,
                &mut threshold_diagnostics,
                result_budget,
            )
        })
        .transpose()?;
    Ok(PcrOutcome {
        records,
        sequence_annotations,
        failure_reason: None,
        threshold_diagnostics,
        uncertain_products,
    })
}

/// Write an annotated DOT file with node and edge attributes for visualization.
///
/// Node attributes: sub_kmer sequence, is_start, is_end.
/// Edge attributes: kmer sequence, kmer count from the graph.
fn write_annotated_dot(
    graph: &petgraph::stable_graph::StableDiGraph<DBNode, DBEdge>,
    kmer_counts: &FilteredKmerCounts,
) -> String {
    use std::fmt::Write as FmtWrite;

    let k = kmer_counts.get_k();
    let sub_k = k - 1;
    let mut dot = String::new();

    // writeln! to String is infallible; unwrap is safe
    writeln!(dot, "digraph {{").unwrap();
    writeln!(dot, "  rankdir=LR;").unwrap();

    for node_idx in graph.node_indices() {
        let node = &graph[node_idx];
        let seq = crate::kmer::kmer_to_seq(&node.sub_kmer, &sub_k);

        let mut attrs = vec![format!("label=\"{}\"", seq)];

        // Shape: start nodes are double-circle, end nodes are box, both are diamond
        if node.is_start && node.is_end {
            attrs.push("shape=diamond".to_string());
        } else if node.is_start {
            attrs.push("shape=doublecircle".to_string());
        } else if node.is_end {
            attrs.push("shape=box".to_string());
        }

        writeln!(dot, "  {} [{}];", node_idx.index(), attrs.join(", ")).unwrap();
    }

    for edge_idx in graph.edge_indices() {
        let (src, tgt) = graph.edge_endpoints(edge_idx).unwrap();
        let edge = &graph[edge_idx];
        let edge_kmer = graph::reconstruct_edge_kmer(graph, edge_idx);
        let seq = crate::kmer::kmer_to_seq(&edge_kmer, &k);
        writeln!(
            dot,
            "  {} -> {} [label=\"{} ({})\"];",
            src.index(),
            tgt.index(),
            seq,
            edge.count
        )
        .unwrap();
    }

    writeln!(dot, "}}").unwrap();
    dot
}

#[cfg(test)]
mod tests {
    use crate::kmer::seq_to_kmer;

    use super::graph::*;
    use super::paths::*;
    use super::primers::*;
    use super::*;
    use petgraph::stable_graph::StableDiGraph;
    use std::collections::HashMap;
    use std::collections::HashSet;

    // These functions are used in the tests
    fn n_combinations(n: usize, r: usize) -> usize {
        // Compute n choose r iteratively to avoid factorial overflow.
        // n*(n-1)*...*(n-r+1) / r!
        if r > n {
            return 0;
        }
        let r = r.min(n - r); // take advantage of symmetry
        let mut result: usize = 1;
        for i in 0..r {
            result = result * (n - i) / (i + 1);
        }
        result
    }

    fn expected_permutations(n: usize, r: usize) -> usize {
        // Given a sequence of length n and r sites that can be permuted,
        // there are (n! / (r!(n-r)!)) combinations of r sites in the sequence and
        // 4^r permutations of the r sites.
        // So there are (n! / (r!(n-r)!)) * 4^r permutations.
        // The original sequence is included in each combination, but we only
        // want to include it once. Se we subtract the number of combinations
        // and add 1.

        n_combinations(n, r) * 4_usize.pow(r as u32) - n_combinations(n, r) + 1
    }

    // Create a test graph of known structure
    fn create_test_graph() -> (
        StableDiGraph<DBNode, DBEdge>,
        HashMap<&'static str, petgraph::graph::NodeIndex>,
    ) {
        let mut graph: StableDiGraph<DBNode, DBEdge> = StableDiGraph::new();
        let mut nodes = HashMap::new();

        // Add nodes and store their indices in the HashMap
        nodes.insert(
            "a",
            graph.add_node(DBNode {
                sub_kmer: 0,
                is_start: true,
                is_end: false,
            }),
        );
        nodes.insert(
            "b",
            graph.add_node(DBNode {
                sub_kmer: 1,
                is_start: false,
                is_end: false,
            }),
        );
        nodes.insert(
            "c",
            graph.add_node(DBNode {
                sub_kmer: 2,
                is_start: false,
                is_end: false,
            }),
        );
        nodes.insert(
            "d",
            graph.add_node(DBNode {
                sub_kmer: 3,
                is_start: false,
                is_end: true,
            }),
        );
        nodes.insert(
            "e",
            graph.add_node(DBNode {
                sub_kmer: 4,
                is_start: false,
                is_end: true,
            }),
        );

        // Add directed edges
        graph.add_edge(
            nodes["a"],
            nodes["b"],
            DBEdge {
                count: 5,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["b"],
            nodes["c"],
            DBEdge {
                count: 10,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["d"],
            DBEdge {
                count: 4,
                coverage_ratio: 0.0,
            },
        );
        graph.add_edge(
            nodes["c"],
            nodes["e"],
            DBEdge {
                count: 1,
                coverage_ratio: 0.0,
            },
        );

        (graph, nodes)
    }

    #[test]
    fn test_n_descendants() {
        let (graph, nodes) = create_test_graph();

        // Testing using the node indices from the HashMap
        assert_eq!(graph::descendants(&graph, nodes["a"], 1).len(), 1); // Direct successor
        assert_eq!(graph::descendants(&graph, nodes["a"], 2).len(), 2); // Including b's successors
        assert_eq!(graph::descendants(&graph, nodes["a"], 3).len(), 4); // All nodes reachable from a within 3 steps
        assert_eq!(graph::descendants(&graph, nodes["a"], 4).len(), 4); // All nodes reachable from a within 4 steps
        assert_eq!(graph::descendants(&graph, nodes["b"], 2).len(), 3);
    }

    #[test]
    fn test_n_combinations() {
        assert_eq!(n_combinations(3, 2), 3);
        assert_eq!(n_combinations(5, 5), 1);
    }

    #[test]
    fn test_combinations() {
        // Make sure there are the correct numbers of combinations
        assert_eq!(combinations(3, 2).len(), n_combinations(3, 2));
        assert_eq!(combinations(8, 5).len(), n_combinations(8, 5));

        // Make sure the first combination has r elements
        assert_eq!(combinations(8, 5)[0].len(), 5);

        // Make sure that the combinations are unique
        let combos = combinations(8, 5);
        let combos_set: HashSet<Vec<usize>> = HashSet::from_iter(combos.clone());
        assert_eq!(combos.len(), combos_set.len());
    }

    #[test]
    fn test_expected_permutations() {
        assert_eq!(expected_permutations(2, 1), 7);
    }

    #[test]
    fn test_resolve_primers() {
        // Check with no ambiguous nucleotides
        let seq1 = "CGTAATGCGGCGA".to_string();
        let mut expected1: HashSet<String> = HashSet::new();
        expected1.insert(seq1.clone());
        let result1 = resolve_primer(&seq1);
        assert_eq!(result1, expected1);

        // Check with one ambiguous nucleotide
        let seq2 = "CGTAATGCGGCGN".to_string();
        let mut expected2: HashSet<String> = HashSet::new();
        expected2.insert("CGTAATGCGGCGA".to_string());
        expected2.insert("CGTAATGCGGCGC".to_string());
        expected2.insert("CGTAATGCGGCGG".to_string());
        expected2.insert("CGTAATGCGGCGT".to_string());

        let result2 = resolve_primer(&seq2);
        assert_eq!(result2, expected2);

        // Check with one ambiguous nucleotide
        let seq3 = "CGTAATRCGGCGA".to_string();
        let mut expected3: HashSet<String> = HashSet::new();
        expected3.insert("CGTAATACGGCGA".to_string());
        expected3.insert("CGTAATGCGGCGA".to_string());
        let result3 = resolve_primer(&seq3);
        assert_eq!(result3, expected3);

        // Check with two ambiguous nucleotides
        let seq4 = "CGTAATRCGGCGY".to_string();
        let mut expected4: HashSet<String> = HashSet::new();
        expected4.insert("CGTAATACGGCGC".to_string());
        expected4.insert("CGTAATGCGGCGC".to_string());
        expected4.insert("CGTAATACGGCGT".to_string());
        expected4.insert("CGTAATGCGGCGT".to_string());

        let result4 = resolve_primer(&seq4);
        assert_eq!(result4, expected4);

        // Check when the first nucleotide is ambiguous
        let seq5 = "RCGTAATCGGCGA".to_string();
        let mut expected5: HashSet<String> = HashSet::new();
        expected5.insert("ACGTAATCGGCGA".to_string());
        expected5.insert("GCGTAATCGGCGA".to_string());
        let result5 = resolve_primer(&seq5);
        assert_eq!(result5, expected5);
    }

    #[test]
    fn test_permute_sequences() {
        // Check specific permutations for a tiny example
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CG".to_string());

            let mut expected: HashSet<String> = HashSet::new();
            expected.insert("CA".to_string());
            expected.insert("CC".to_string());
            expected.insert("CG".to_string());
            expected.insert("CT".to_string());
            expected.insert("AG".to_string());
            expected.insert("GG".to_string());
            expected.insert("TG".to_string());

            let result = permute_sequences(seq, &1_usize);
            println!(
                "Permutations: {}",
                result.iter().cloned().collect::<Vec<_>>().join(", ")
            );
            assert_eq!(result, expected);
        }

        // Check sequence length 3 and r=3
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());

            let result = permute_sequences(seq, &3_usize);

            // All sites permuted, so there should be 4^n sequences
            assert_eq!(result.len(), 64);
        }

        // Construct all the permutations procedurally
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("CGT".to_string());
            let r: usize = 2;
            let n: usize = 3;

            let bases = ['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("C{}{}", bases[i], bases[j]));
                    expected.insert(format!("{}G{}", bases[i], bases[j]));
                    expected.insert(format!("{}{}T", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }

        // Procedural test where n=4 and r=2
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("ACGT".to_string());
            let r: usize = 2;
            let n: usize = 4;

            let bases = ['A', 'C', 'G', 'T'];
            let mut expected: HashSet<String> = HashSet::new();
            for i in 0..4 {
                for j in 0..4 {
                    expected.insert(format!("{}{}GT", bases[i], bases[j]));
                    expected.insert(format!("{}C{}T", bases[i], bases[j]));
                    expected.insert(format!("{}CG{}", bases[i], bases[j]));
                    expected.insert(format!("A{}{}T", bases[i], bases[j]));
                    expected.insert(format!("A{}G{}", bases[i], bases[j]));
                    expected.insert(format!("AC{}{}", bases[i], bases[j]));
                }
            }

            println!(
                "There are {} combinations for n={} r={}",
                expected.len(),
                n,
                r
            );
            let result = permute_sequences(seq, &r);
            assert_eq!(expected.len(), result.len());
        }

        // Check for inclusion of particular strings as a spot check
        {
            let mut seq: HashSet<String> = HashSet::new();
            seq.insert("TGCAGGTTCACCTAC".to_string());
            let r: usize = 2;
            let result = permute_sequences(seq, &r);

            // Check for inclusion of original sequence
            assert!(result.contains("TGCAGGTTCACCTAC"));

            // Check for inclusion of "GGCAGGTTCACCTAC"
            assert!(result.contains("GGCAGGTTCACCTAC"));
        }
    }

    #[test]
    fn test_string_to_oligo() {
        let oligo = string_to_oligo("GCGA").unwrap();
        assert_eq!(oligo.kmer, 0b1001_1000);
        assert_eq!(oligo.length, 4);
    }

    #[test]
    fn test_get_end_nodes() {
        let (graph, nodes) = create_test_graph();

        let end_nodes = get_end_nodes(&graph);
        assert_eq!(end_nodes.len(), 2);
        assert!(end_nodes.contains(&nodes["d"]));
        assert!(end_nodes.contains(&nodes["e"]));
    }

    #[test]
    fn test_get_start_nodes() {
        let (graph, nodes) = create_test_graph();

        let start_nodes = get_start_nodes(&graph);
        assert_eq!(start_nodes.len(), 1);
        assert!(start_nodes.contains(&nodes["a"]));
    }

    #[test]
    fn test_get_suffix_mask() {
        // If k is 21, then suffix should be 20 nucleotides long, and there are two bits per nucleotide
        let k: usize = 21;
        assert_eq!(
            graph::get_suffix_mask(&k),
            0b11111111_11111111_11111111_11111111_11111111
        );

        let k: usize = 3;
        assert_eq!(graph::get_suffix_mask(&k), 0b1111);
    }

    fn build_test_case() -> (String, usize, usize, KmerCounts, PCRParams) {
        // This is the 18s that was assembled from https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR26955578, per the readme tutorial
        // Padded with C's at the ends
        let read_string = "CCCCCCCCCCCCGTTGATCCTGCCAGTATCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCACTTGTACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATTGTACTCTCTTACTACTTGGATAACCGTAGTAATTCTAGAGCTAATACATGCGAAAAGTCCCGACTCTCGTGGAAGGGATGTATTTATTAGATTAAAAACCAATGCGGCTTAACGGCCGCTTACAAACTTGGTGATTCATAGTAACTGTTCGAATCGCATGGCCTTCTCTGTTCGTGCCGGCGATGTTTCATTCAAATTTCTGCCCTATCAACTGTCGATGGTAAGATAGTGGCTTACCATGGTCGCAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACGTGGGGAGGTAGTGACAAAAAATAACAATACAGGGCTTTTTTGTAGTCTTGTAATTGGAATGAGTACAATTTAAATCTCTTAACGAGGACCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATTGTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGGCTCGACGGCGACGGTCAGCCGCAAGGTATGTCACTGTCGACGTTGGCCTTCTTCGCGCAGACTTCGCGTGCTCTTAACTGAGTGTGCGTTGGATACGCGACGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCTTGTGCTTGGATACATAAGCATGGAATAATGGAATAGGACTTTGGTTCTATTTTCCGTTGGTTTCTGGAACCGAAGTAATGATTAATAGGGACAGTTGGGGGCATTCGTATTTCGTTGTCAGAGGTGAAATTCTTGGATTTACGAAAGACGAACTAATGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGATCGAAGACGATCAGATACCGTCCTAGTTCTAACCATAAACGATGCCGACTAGGGATCAGCGAGTGTTATTTGATGACCTCGTTGGCACCTTATGGGAAACCAAAGTTTTTGGGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTCACCAGGTCCAGACATAGTAAGGATTGACAGATTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTGACCTGCTAAATAGTCAGACGATTCTCGAATCGCTCTCGACTTCTTAGAGGGACTGTTGCGTGTTTAACCAAAGTCAGGAAGGCAATAACAGGTCTGTGATGCCCTTAGATGTCCTGGGCCGCACGCGCGCTACACTGACGATGGCAACGAGTCGCTCCTTCACCGAAAGGTGTGGGTAATCTTGTGAATCATCGTCGTGCTGGGGATAGATCATTGTAATTCTTGATCTTGAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGGTTTAGTGAGGCCTCCGGATTGGCACTGTCAGATGGGCTTCGGTCCATCCGACGGACGTCAAAAAGTTGGTCAAACTTGATCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGCCCCCCCCCCCC".to_string();

        let k: usize = 21;
        let replicates: usize = 10;
        let mut kmer_counts = KmerCounts::new(&k);
        for _i in 0..replicates {
            let reads = crate::kmer::seq_to_reads(&read_string).unwrap();
            kmer_counts.ingest_reads(&reads).unwrap();
        }

        let params = PCRParams {
            forward_seq: "AACCTGGTTGATCCTGCCAGT".to_string(),

            // AGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCG 3' end of 18S above
            // CGCAGGTTCACCTACGGAAACCTTGTTACGACTTTTACT revcomp of above
            reverse_seq: "TGATCCTTCTGCAGGTTCACCTAC".to_string(),
            min_length: 0,
            max_length: 2500,
            gene_name: "18s".to_string(),
            gene: None,
            region: None,
            index: None,
            compartment: None,
            gene_type: None,
            copy_number: None,
            deprecated: false,
            deprecated_by: None,
            deprecated_reason: None,
            min_count: 3,
            mismatches: 2,
            trim: 15,
            expected_length: None,
            citation: "".to_string(),
            notes: "".to_string(),
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
        };

        (read_string, k, replicates, kmer_counts, params)
    }

    /// Exercise the exact store through actual primer discovery and graph PCR,
    /// rather than converting it back to the legacy global table.
    #[test]
    fn exact_store_pcr_matches_legacy_and_preserves_large_counts() {
        use crate::kmer::identity::{CanonicalKey, KmerSize};
        use crate::kmer::store::MemoryCountBuilder;
        use std::num::NonZeroU64;

        let (_, k, _, counts, mut params) = build_test_case();
        params.emit_uncertain = true;
        let run = |view: &FilteredKmerCounts| {
            do_pcr(
                view,
                "count_store",
                &params,
                false,
                "",
                None,
                100_000,
                withheld_diagnostic_limits_for_gene(0, 1),
            )
            .unwrap()
        };
        let baseline = run(&counts.filtered_view(1));
        assert!(!baseline.records.is_empty());
        let size = KmerSize::new(k).unwrap();
        let mut builder = MemoryCountBuilder::new(size, 1 << 20, 16).unwrap();
        for (&bits, &count) in counts.iter() {
            builder
                .add(
                    CanonicalKey::new(bits, size).unwrap(),
                    NonZeroU64::new(u64::from(count)).unwrap(),
                )
                .unwrap();
        }
        let store = builder.finalize([0; 32]);
        let exact = run(&FilteredKmerCounts::from_store(&store, 1));
        let inventory = |outcome: &PcrOutcome| {
            outcome
                .records
                .iter()
                .map(|r| {
                    (
                        r.id().to_owned(),
                        r.desc().map(str::to_owned),
                        r.seq().to_vec(),
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(inventory(&exact), inventory(&baseline));
        assert_eq!(exact.threshold_diagnostics, baseline.threshold_diagnostics);
        assert_eq!(exact.failure_reason, baseline.failure_reason);
        assert_eq!(
            serde_yaml_ng::to_value(&exact.uncertain_products).unwrap(),
            serde_yaml_ng::to_value(&baseline.uncertain_products).unwrap()
        );

        let multiplier = u64::from(u32::MAX) + 1;
        let mut builder = MemoryCountBuilder::new(size, 1 << 20, 16).unwrap();
        for (&bits, &count) in counts.iter() {
            builder
                .add(
                    CanonicalKey::new(bits, size).unwrap(),
                    NonZeroU64::new(u64::from(count) * multiplier).unwrap(),
                )
                .unwrap();
        }
        let store = builder.finalize([0; 32]);
        let large = run(&FilteredKmerCounts::from_store(&store, 1));
        assert_eq!(large.records.len(), baseline.records.len());
        for (record, expected) in large.records.iter().zip(&baseline.records) {
            assert_eq!(record.seq(), expected.seq());
            assert!(
                record
                    .desc()
                    .unwrap()
                    .contains(&format!("kmer_count_min={}", 10 * multiplier))
            );
        }
        assert!(large.threshold_diagnostics[0].threshold > u64::from(u32::MAX));
        assert_eq!(compute_coverage_thresholds(u64::MAX, 3).last(), Some(&3));
        assert_eq!(graph::compute_mean(&[u64::MAX, u64::MAX]), u64::MAX as f64);
    }

    #[test]
    fn retained_sequence_results_release_single_target_admission_before_publication() {
        use crate::kmer::budget::{MEMORY_CLASSES, ResourceLimits};
        use std::sync::mpsc;
        use std::time::Duration;
        let (_, _, _, counts, params) = build_test_case();
        // Let the existing 18S/mismatch fixture finish graph construction;
        // this test isolates result ownership, not a graph-size ceiling.
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::TargetGraphs as usize] = 8 << 20;
        classes[MemoryClass::RetainedResults as usize] = 16 << 20;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 16 << 20,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let pool = Arc::new(
            memory::admission::GraphBufferPool::new(Arc::clone(&budget), 8 << 20).unwrap(),
        );
        let worker_pool = Arc::clone(&pool);
        let (sender, receiver) = mpsc::channel();
        let worker = std::thread::spawn(move || {
            for _ in 0..3 {
                let lease = worker_pool.acquire().unwrap();
                let result = do_pcr_with_graph_budget(
                    &counts.filtered_view(1),
                    "retained",
                    &params,
                    false,
                    "",
                    None,
                    100_000,
                    withheld_diagnostic_limits_for_gene(0, 1),
                    Some(lease.budget()),
                    Some(worker_pool.result_budget()),
                )
                .unwrap();
                drop(lease);
                sender.send(result).unwrap();
            }
        });
        let mut outcomes = Vec::new();
        for _ in 0..3 {
            outcomes.push(receiver.recv_timeout(Duration::from_secs(10)).unwrap());
        }
        worker.join().unwrap();
        assert!(outcomes.iter().all(|outcome| !outcome.records.is_empty()));
        assert_eq!(pool.summary().unwrap().peak_active_targets, 1);
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        assert!(budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize] > 0);
        for pair in outcomes.windows(2) {
            assert_eq!(pair[0].records, pair[1].records);
            assert_eq!(pair[0].sequence_annotations, pair[1].sequence_annotations);
        }
        drop(outcomes);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn retained_sequence_refusal_is_typed_and_drops_graph_and_result_owners() {
        use crate::kmer::budget::{
            MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits,
        };
        let (_, _, _, counts, params) = build_test_case();
        // Let the existing 18S/mismatch fixture finish graph construction;
        // this test isolates result ownership, not a graph-size ceiling.
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::TargetGraphs as usize] = 8 << 20;
        classes[MemoryClass::RetainedResults as usize] = 1;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 16 << 20,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let pool = memory::admission::GraphBufferPool::new(Arc::clone(&budget), 8 << 20).unwrap();
        let lease = pool.acquire().unwrap();
        let error = do_pcr_with_graph_budget(
            &counts.filtered_view(1),
            "retained",
            &params,
            false,
            "",
            None,
            100_000,
            withheld_diagnostic_limits_for_gene(0, 1),
            Some(lease.budget()),
            Some(pool.result_budget()),
        )
        .err()
        .unwrap();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::StageMemory(MemoryClass::RetainedResults)
        );
        assert_eq!(
            budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize],
            0
        );
        drop(lease);
        assert_eq!(pool.summary().unwrap().peak_active_targets, 1);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn owned_threshold_failure_summary_preserves_order_and_refuses_without_leaks() {
        use crate::kmer::budget::{MEMORY_CLASSES, ResourceLimits};
        let diagnostics = [
            PcrThresholdDiagnostic {
                threshold: 8,
                node_budget_reached: true,
                dfs_limit_reached: true,
                path_quota_reached: true,
                withheld_candidate_paths: 1,
                ..PcrThresholdDiagnostic::default()
            },
            PcrThresholdDiagnostic {
                threshold: 4,
                dfs_limit_reached: true,
                ..PcrThresholdDiagnostic::default()
            },
            PcrThresholdDiagnostic {
                threshold: 2,
                ..PcrThresholdDiagnostic::default()
            },
        ];
        for bytes in [1, 4096] {
            let mut classes = [0; MEMORY_CLASSES];
            classes[MemoryClass::RetainedResults as usize] = bytes;
            let budget = ResourceBudget::new(ResourceLimits {
                memory_bytes: bytes,
                process_overhead_bytes: 0,
                class_bytes: classes,
                temporary_disk_bytes: 0,
                temporary_open_files: 0,
            })
            .unwrap();
            match earlier_threshold_limit_summary_budgeted(&diagnostics, Some(&budget)) {
                Ok(summary) => {
                    assert_eq!(
                        summary.unwrap().as_str(),
                        "threshold 8: node budget reached, DFS state limit reached, path quota reached, candidate-local repeat uncertainty; threshold 4: DFS state limit reached"
                    );
                    assert_eq!(bytes, 4096);
                }
                Err(error) => {
                    assert_eq!(bytes, 1);
                    assert!(
                        error
                            .downcast_ref::<crate::kmer::budget::ResourceExhausted>()
                            .is_some()
                    );
                }
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        }
    }

    #[test]
    fn count_store_failures_propagate_from_primer_scan_and_graph_lookup() {
        use crate::kmer::identity::{CanonicalKey, KmerSize};
        use crate::kmer::store::{CountMetadata, CountStore, MemoryCountBuilder, MemoryCountStore};
        use std::num::NonZeroU64;

        struct FailingStore {
            store: MemoryCountStore,
            fail_visit: bool,
        }
        impl CountStore for FailingStore {
            fn metadata(&self) -> &CountMetadata {
                self.store.metadata()
            }
            fn lookup(&self, _: CanonicalKey) -> Result<Option<NonZeroU64>> {
                anyhow::bail!("injected exact lookup failure")
            }
            fn visit(
                &self,
                visitor: &mut dyn FnMut(CanonicalKey, NonZeroU64) -> Result<()>,
            ) -> Result<()> {
                if self.fail_visit {
                    anyhow::bail!("injected exact scan failure");
                }
                self.store.visit(visitor)
            }
        }
        let (_, k, _, counts, params) = build_test_case();
        for fail_visit in [true, false] {
            let size = KmerSize::new(k).unwrap();
            let mut builder = MemoryCountBuilder::new(size, 1 << 20, 16).unwrap();
            for (&bits, &count) in counts.iter() {
                builder
                    .add(
                        CanonicalKey::new(bits, size).unwrap(),
                        NonZeroU64::new(u64::from(count)).unwrap(),
                    )
                    .unwrap();
            }
            let store = FailingStore {
                store: builder.finalize([0; 32]),
                fail_visit,
            };
            let result = do_pcr(
                &FilteredKmerCounts::from_store(&store, 1),
                "failure",
                &params,
                false,
                "",
                None,
                100_000,
                withheld_diagnostic_limits_for_gene(0, 1),
            );
            let error = match result {
                Ok(_) => panic!("store failure became a PCR outcome"),
                Err(error) => error,
            };
            assert!(error.to_string().contains(if fail_visit {
                "scan failure"
            } else {
                "lookup failure"
            }));
        }
    }

    #[test]
    fn test_primer_preprocessing_steps() {
        let (_, _, _, kmer_counts, params) = build_test_case();
        let kmer_counts = kmer_counts.filtered_view(1);

        let reverse_variants =
            preprocess_primer(&params, PrimerDirection::Reverse, &kmer_counts.get_k()).unwrap();

        // There should be 991 variants of the reverse primer when r=2
        assert_eq!(reverse_variants.len(), 991);

        println!("Reverse variants: {:?}", reverse_variants);

        // Check for inclusion of original trimmed sequence TGCAGGTTCACCTAC
        let member = "TGCAGGTTCACCTAC".to_string();
        assert!(reverse_variants.contains(&member));

        // Check for inclusion of off by one variant GGCAGGTTCACCTAC
        let member = "GGCAGGTTCACCTAC".to_string();
        assert!(reverse_variants.contains(&member));

        let mut reverse_primer_kmers =
            get_kmers_from_primers(&reverse_variants, &kmer_counts, &params.min_count).unwrap();

        // Check for kmer
        assert_eq!(reverse_primer_kmers.len(), 1);

        reverse_primer_kmers =
            filter_primer_kmers(reverse_primer_kmers, DEFAULT_MAX_NUM_PRIMER_KMERS).unwrap();

        // Check for kmer after filtering
        assert_eq!(reverse_primer_kmers.len(), 1);
    }

    #[test]
    fn test_extension_steps() {
        let (_read_string, _k, _replicates, kmer_counts, _params) = build_test_case();
        let kmer_counts = kmer_counts.filtered_view(1);

        let seq = "TGATCCTGCCAGTATCATATG".to_string();
        let kmer: u64 = seq_to_kmer(&seq).unwrap();
        assert!(kmer_counts.get_canonical(&kmer).unwrap().is_some());
    }

    #[test]
    fn test_integration() {
        let (read_string, k, replicates, kmer_counts, params) = build_test_case();
        let min_count = 5;

        // Check the number of kmers (canonical only, no reverse complements stored)
        assert_eq!(kmer_counts.len(), read_string.len() - k + 1);

        // Check the total count of kmers
        assert_eq!(
            kmer_counts.get_n_kmers(),
            ((read_string.len() - k + 1) * replicates) as u64
        );

        let filtered = kmer_counts.filtered_view(1);

        let (forward_primer_kmers, reverse_primer_kmers) =
            primers::get_primer_kmers(&params, &filtered, None).unwrap();

        assert_eq!(forward_primer_kmers.len(), 1);
        assert_eq!(reverse_primer_kmers.len(), 1);

        let (seed_graph, node_lookup) = graph::create_seed_graph(
            &forward_primer_kmers,
            &reverse_primer_kmers,
            &filtered,
            None,
        )
        .unwrap();

        // Check the number of nodes in the seed graph
        // With unified primer handling, both primers use prefix (kmer >> 2).
        // The forward and reverse primer kmers produce distinct sub_kmers,
        // so we still expect 2 nodes.
        assert_eq!(seed_graph.node_count(), 2);
        assert_eq!(get_start_nodes(&seed_graph).len(), 1);
        assert_eq!(get_end_nodes(&seed_graph).len(), 1);

        let (
            mut graph_result,
            _node_lookup_final,
            _found,
            _extension_repeat_markers,
            _node_budget_reached,
        ) = graph::extend_graph(
            seed_graph,
            node_lookup,
            &filtered,
            &min_count,
            &params,
            0,
            graph::DEFAULT_MAX_NUM_NODES,
        )
        .unwrap();

        // Print the number of nodes and edges in the graph
        println!("There are {} nodes in the graph", graph_result.node_count());
        println!("There are {} edges in the graph", graph_result.edge_count());

        // With reverse extension, forward extension from the start node
        // and reverse extension from the end node should converge, producing
        // complete paths.
        let all_paths = get_assembly_paths(&graph_result, &filtered, &params, None);
        assert!(
            !all_paths.is_empty(),
            "Expected paths after reverse extension"
        );

        graph_result
            .prune_tips(&filtered.get_k(), DEFAULT_TIP_COVERAGE_FRACTION)
            .unwrap();
        graph_result.prune_reachability().unwrap();

        let all_paths = get_assembly_paths(&graph_result, &filtered, &params, None);
        assert!(!all_paths.is_empty(), "Expected paths after pruning");
    }
}

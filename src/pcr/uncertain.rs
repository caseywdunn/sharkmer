//! Separate, bounded output of unsupported graph-path hypotheses.
//!
//! This consumes the existing diagnostic collector after ordinary admission and
//! deduplication. It never adds a path, extends a search, or selects a threshold.

use bio::io::fasta;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

use super::{
    PcrThresholdDiagnostic, WithheldDiagnosticLimits, WithheldDiagnosticUsage,
    WithheldMarkerCoveredRun, WithheldMarkerOccurrence, WithheldPathDiagnostics,
    WithheldPathRecord,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum MarkerDetail {
    Complete,
    SummaryOnly,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProductSource {
    pub(crate) threshold: u32,
    pub(crate) marker_detail: MarkerDetail,
    pub(crate) node_visit_limit: usize,
    pub(crate) max_observed_node_visits: usize,
    pub(crate) omitted_self_loop_marker_occurrences: usize,
    pub(crate) cyclic_scc_marker_occurrences: usize,
    pub(crate) retained_collision_marker_occurrences: usize,
    pub(crate) marker_occurrences: Vec<WithheldMarkerOccurrence>,
    pub(crate) marker_covered_runs: Vec<WithheldMarkerCoveredRun>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProduct {
    pub(crate) id: String,
    pub(crate) description: String,
    pub(crate) oriented_sequence: String,
    pub(crate) sequence_sha256: String,
    /// Length of the emitted hypothesis, not an established biological length.
    pub(crate) sequence_length: usize,
    pub(crate) sources: Vec<UncertainProductSource>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainRetentionSummary {
    pub(crate) gene_usage_before_threshold: WithheldDiagnosticUsage,
    pub(crate) gene_usage_after_threshold: WithheldDiagnosticUsage,
    pub(crate) observed_withheld_paths: usize,
    pub(crate) retained_paths: usize,
    pub(crate) observed_not_retained_paths: usize,
    pub(crate) retained_sequence_bases: usize,
    pub(crate) retained_marker_occurrences: usize,
    pub(crate) retained_summary_only_paths: usize,
    pub(crate) omitted_detailed_marker_occurrences: usize,
    pub(crate) marker_detail_truncated: bool,
    pub(crate) dropped_oversize_sequence: usize,
    pub(crate) dropped_oversize_marker_set: usize,
    pub(crate) dropped_path_cap: usize,
    pub(crate) dropped_sequence_base_cap: usize,
    pub(crate) dropped_marker_cap: usize,
    pub(crate) node_visit_limit: usize,
    pub(crate) node_visit_skips: usize,
    pub(crate) retention_truncated: bool,
}

impl From<&WithheldPathDiagnostics> for UncertainRetentionSummary {
    fn from(value: &WithheldPathDiagnostics) -> Self {
        Self {
            gene_usage_before_threshold: value.gene_usage_before_threshold,
            gene_usage_after_threshold: value.gene_usage_after_threshold,
            observed_withheld_paths: value.observed_withheld_paths,
            retained_paths: value.retained_paths,
            observed_not_retained_paths: value.observed_not_retained_paths,
            retained_sequence_bases: value.retained_sequence_bases,
            retained_marker_occurrences: value.retained_marker_occurrences,
            retained_summary_only_paths: value.retained_summary_only_paths,
            omitted_detailed_marker_occurrences: value.omitted_detailed_marker_occurrences,
            marker_detail_truncated: value.marker_detail_truncated,
            dropped_oversize_sequence: value.dropped_oversize_sequence,
            dropped_oversize_marker_set: value.dropped_oversize_marker_set,
            dropped_path_cap: value.dropped_path_cap,
            dropped_sequence_base_cap: value.dropped_sequence_base_cap,
            dropped_marker_cap: value.dropped_marker_cap,
            node_visit_limit: value.node_visit_limit,
            node_visit_skips: value.node_visit_skips,
            retention_truncated: value.retention_truncated,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainThresholdSummary {
    /// The standard search counters are preserved, without diagnostic sequences.
    #[serde(flatten)]
    pub(crate) search: PcrThresholdDiagnostic,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) retention: Option<UncertainRetentionSummary>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProducts {
    pub(crate) schema_version: u32,
    pub(crate) purpose: String,
    pub(crate) coordinate_system: String,
    pub(crate) marker_span_scope: String,
    pub(crate) read_support_evaluation: String,
    pub(crate) repeat_length_evaluation: String,
    pub(crate) deduplication_policy: String,
    pub(crate) standard_equivalence_policy: String,
    pub(crate) limits: WithheldDiagnosticLimits,
    pub(crate) path_search_evaluated: bool,
    pub(crate) search_limit_scope: String,
    pub(crate) search_limits_reached: bool,
    pub(crate) retention_truncated: bool,
    pub(crate) marker_detail_truncated: bool,
    pub(crate) thresholds: Vec<UncertainThresholdSummary>,
    pub(crate) retained_path_records: usize,
    pub(crate) excluded_standard_equivalent_path_records: usize,
    pub(crate) exact_duplicate_path_records: usize,
    /// Counts oriented sequence hypotheses, not loci, haplotypes, or molecules.
    pub(crate) unique_products: usize,
    pub(crate) products: Vec<UncertainProduct>,
}

impl UncertainProducts {
    pub(super) fn empty(limits: WithheldDiagnosticLimits) -> Self {
        Self {
            schema_version: 2,
            purpose: "uncertain_repeat_structure_graph_path_hypotheses".to_string(),
            coordinate_system: "zero_based_half_open_in_emitted_oriented_sequence".to_string(),
            marker_span_scope: "marker_covered_positions_only_for_complete_detail_sources"
                .to_string(),
            read_support_evaluation: "not_evaluated".to_string(),
            repeat_length_evaluation: "unresolved".to_string(),
            deduplication_policy:
                "exact_oriented_sequence_only_reverse_complements_retained_as_separate_hypotheses"
                    .to_string(),
            standard_equivalence_policy:
                "exact_sequence_or_reverse_complement_of_final_standard_output".to_string(),
            limits,
            path_search_evaluated: false,
            search_limit_scope: "all_attempted_thresholds_for_gene".to_string(),
            search_limits_reached: false,
            retention_truncated: false,
            marker_detail_truncated: false,
            thresholds: Vec::new(),
            retained_path_records: 0,
            excluded_standard_equivalent_path_records: 0,
            exact_duplicate_path_records: 0,
            unique_products: 0,
            products: Vec::new(),
        }
    }
}

fn canonical_sequence(sequence: &[u8]) -> Vec<u8> {
    let reverse = bio::alphabets::dna::revcomp(sequence);
    if sequence <= reverse.as_slice() {
        sequence.to_vec()
    } else {
        reverse
    }
}

fn source(threshold: u32, path: WithheldPathRecord) -> UncertainProductSource {
    let marker_detail = if path.marker_occurrences.is_empty() {
        MarkerDetail::SummaryOnly
    } else {
        MarkerDetail::Complete
    };
    UncertainProductSource {
        threshold,
        marker_detail,
        node_visit_limit: path.node_visit_limit,
        max_observed_node_visits: path.max_observed_node_visits,
        omitted_self_loop_marker_occurrences: path.omitted_self_loop_marker_occurrences,
        cyclic_scc_marker_occurrences: path.cyclic_scc_marker_occurrences,
        retained_collision_marker_occurrences: path.retained_collision_marker_occurrences,
        marker_occurrences: path.marker_occurrences,
        marker_covered_runs: path.marker_covered_runs,
    }
}

/// Finalize only from actual in-range paths retained by the unchanged collector.
pub(super) fn collect(
    sample_name: &str,
    gene_name: &str,
    keep_hidden_diagnostics: bool,
    limits: WithheldDiagnosticLimits,
    standard_records: &[fasta::Record],
    threshold_diagnostics: &mut [PcrThresholdDiagnostic],
) -> UncertainProducts {
    let mut output = UncertainProducts::empty(limits);
    let standard_sequences: BTreeSet<Vec<u8>> = standard_records
        .iter()
        .map(|record| canonical_sequence(record.seq()))
        .collect();
    let mut products = BTreeMap::<String, UncertainProduct>::new();

    for diagnostic in threshold_diagnostics {
        let payload = diagnostic.withheld_path_diagnostics.take();
        let search = diagnostic.clone();
        if keep_hidden_diagnostics {
            diagnostic.withheld_path_diagnostics = payload.clone();
        }
        output.path_search_evaluated |= search.path_search_evaluated;
        output.search_limits_reached |=
            search.node_budget_reached || search.dfs_limit_reached || search.path_quota_reached;
        let retention = payload.as_ref().map(UncertainRetentionSummary::from);
        if let Some(ref retention) = retention {
            output.search_limits_reached |= retention.node_visit_skips > 0;
            output.retention_truncated |= retention.retention_truncated;
            output.marker_detail_truncated |= retention.marker_detail_truncated;
        }
        output
            .thresholds
            .push(UncertainThresholdSummary { search, retention });

        if let Some(payload) = payload {
            for path in payload.paths {
                output.retained_path_records += 1;
                if standard_sequences
                    .contains(&canonical_sequence(path.oriented_sequence.as_bytes()))
                {
                    output.excluded_standard_equivalent_path_records += 1;
                    continue;
                }
                if let Some(product) = products.get_mut(&path.oriented_sequence) {
                    output.exact_duplicate_path_records += 1;
                    product.sources.push(source(diagnostic.threshold, path));
                } else {
                    let sequence = path.oriented_sequence.clone();
                    let product = UncertainProduct {
                        id: String::new(),
                        description: String::new(),
                        oriented_sequence: sequence.clone(),
                        sequence_sha256: path.sequence_sha256.clone(),
                        sequence_length: path.sequence_length,
                        sources: vec![source(diagnostic.threshold, path)],
                    };
                    products.insert(sequence, product);
                }
            }
        }
    }

    for (index, mut product) in products.into_values().enumerate() {
        // Standard IDs always end in an integer. A terminal status suffix
        // avoids collisions even when a different gene ends in `.uncertain`.
        product.id = format!("{}_{}_{}.uncertain", sample_name, gene_name, index);
        let thresholds = product
            .sources
            .iter()
            .map(|source| source.threshold)
            .collect::<BTreeSet<_>>();
        let thresholds = thresholds
            .into_iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let mut causes = BTreeSet::new();
        let mut marker_detail = MarkerDetail::Complete;
        for source in &product.sources {
            if source.marker_detail == MarkerDetail::SummaryOnly {
                marker_detail = MarkerDetail::SummaryOnly;
            }
            if source.omitted_self_loop_marker_occurrences > 0 {
                causes.insert("pre_prune_omitted_self_loop_node");
            }
            if source.cyclic_scc_marker_occurrences > 0 {
                causes.insert("pre_prune_cyclic_scc_node");
            }
            if source.retained_collision_marker_occurrences > 0 {
                causes.insert("retained_collision_edge");
            }
        }
        product.description = format!(
            "assembly_status=uncertain uncertainty=repeat_structure repeat_length=unresolved read_support=not_evaluated marker_detail={} hypothesis_length={} sequence_sha256={} causes={} thresholds={} source_path_records={} enumeration=bounded limit_scope=gene search_limits_reached={} retention_truncated={}",
            if marker_detail == MarkerDetail::Complete {
                "complete"
            } else {
                "summary_only"
            },
            product.sequence_length,
            product.sequence_sha256,
            causes.into_iter().collect::<Vec<_>>().join(","),
            thresholds,
            product.sources.len(),
            output.search_limits_reached,
            output.retention_truncated,
        );
        output.products.push(product);
    }
    output.unique_products = output.products.len();
    debug_assert_eq!(
        output.retained_path_records,
        output.excluded_standard_equivalent_path_records
            + output.exact_duplicate_path_records
            + output.unique_products
    );
    debug_assert_eq!(
        output
            .products
            .iter()
            .map(|product| product.sources.len())
            .sum::<usize>(),
        output.retained_path_records - output.excluded_standard_equivalent_path_records
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcr::{WithheldMarkerCause, WithheldMarkerIdentityKind};
    use sha2::{Digest, Sha256};

    fn path(sequence: &str) -> WithheldPathRecord {
        let identity = sequence[..3].to_string();
        WithheldPathRecord {
            oriented_sequence: sequence.to_string(),
            sequence_sha256: format!("{:x}", Sha256::digest(sequence.as_bytes())),
            sequence_length: sequence.len(),
            node_visit_limit: 2,
            max_observed_node_visits: 1,
            omitted_self_loop_marker_occurrences: 0,
            cyclic_scc_marker_occurrences: 1,
            retained_collision_marker_occurrences: 0,
            marker_occurrences: vec![WithheldMarkerOccurrence {
                cause: WithheldMarkerCause::PrePruneCyclicSccNode,
                identity_kind: WithheldMarkerIdentityKind::OrientedNodeSequence,
                identity_sha256: format!("{:x}", Sha256::digest(identity.as_bytes())),
                identity_sequence: identity,
                start: 0,
                end: 3,
            }],
            marker_covered_runs: vec![WithheldMarkerCoveredRun { start: 0, end: 3 }],
        }
    }

    fn threshold(threshold: u32, paths: Vec<WithheldPathRecord>) -> PcrThresholdDiagnostic {
        PcrThresholdDiagnostic {
            threshold,
            connectivity_found: true,
            path_search_evaluated: true,
            completed_candidate_paths: paths.len(),
            withheld_candidate_paths: paths.len(),
            withheld_path_diagnostics: Some(WithheldPathDiagnostics {
                observed_withheld_paths: paths.len(),
                retained_paths: paths.len(),
                paths,
                ..WithheldPathDiagnostics::default()
            }),
            ..PcrThresholdDiagnostic::default()
        }
    }

    fn collect_test(
        diagnostics: &mut [PcrThresholdDiagnostic],
        standard: &[fasta::Record],
        keep: bool,
    ) -> UncertainProducts {
        collect(
            "sample",
            "gene",
            keep,
            WithheldDiagnosticLimits::default(),
            standard,
            diagnostics,
        )
    }

    #[test]
    fn exact_duplicates_merge_threshold_provenance_without_collapsing_repeat_lengths() {
        let mut diagnostics = vec![
            threshold(8, vec![path("ACGTA"), path("ACGTTA")]),
            threshold(4, vec![path("ACGTA")]),
        ];
        let output = collect_test(&mut diagnostics, &[], false);
        assert_eq!(output.unique_products, 2);
        assert_eq!(output.retained_path_records, 3);
        assert_eq!(output.exact_duplicate_path_records, 1);
        assert_eq!(output.products[0].sources.len(), 2);
        assert!(output.products[0].description.contains("thresholds=4,8"));
        assert!(
            diagnostics
                .iter()
                .all(|diagnostic| diagnostic.withheld_path_diagnostics.is_none())
        );
    }

    #[test]
    fn only_exact_or_reverse_complement_standard_sequences_are_excluded() {
        let mut diagnostics = vec![threshold(
            4,
            vec![path("ACGTA"), path("TACGT"), path("ACGTTA")],
        )];
        let standard = vec![fasta::Record::with_attrs("accepted", None, b"ACGTA")];
        let output = collect_test(&mut diagnostics, &standard, false);
        assert_eq!(output.excluded_standard_equivalent_path_records, 2);
        assert_eq!(output.unique_products, 1);
        assert_eq!(output.products[0].oriented_sequence, "ACGTTA");
        assert_eq!(
            output.retained_path_records,
            output.excluded_standard_equivalent_path_records
                + output.exact_duplicate_path_records
                + output.unique_products
        );
    }

    #[test]
    fn reverse_complement_uncertain_records_keep_their_own_orientation_and_coordinates() {
        let mut diagnostics = vec![threshold(4, vec![path("ACGTA"), path("TACGT")])];
        let output = collect_test(&mut diagnostics, &[], false);
        assert_eq!(output.unique_products, 2);
        for product in output.products {
            for marker in &product.sources[0].marker_occurrences {
                assert_eq!(
                    &product.oriented_sequence[marker.start..marker.end],
                    marker.identity_sequence
                );
            }
            assert!(product.description.contains("assembly_status=uncertain"));
            assert!(product.description.contains("repeat_length=unresolved"));
            assert!(product.description.contains("read_support=not_evaluated"));
        }
    }

    #[test]
    fn diagnostic_flag_keeps_original_payload_but_does_not_multiply_collection() {
        let original = threshold(4, vec![path("ACGTA")]);
        let mut kept = vec![original.clone()];
        let mut moved = vec![original.clone()];
        let with_diagnostics = collect_test(&mut kept, &[], true);
        let without_diagnostics = collect_test(&mut moved, &[], false);
        assert_eq!(with_diagnostics, without_diagnostics);
        assert_eq!(kept[0], original);
        assert!(moved[0].withheld_path_diagnostics.is_none());
    }

    #[test]
    fn caps_and_search_limits_remain_visible_when_no_records_were_retained() {
        let mut diagnostic = threshold(4, vec![]);
        diagnostic.node_budget_reached = true;
        diagnostic.dfs_limit_reached = true;
        diagnostic.path_quota_reached = true;
        let retention = diagnostic.withheld_path_diagnostics.as_mut().unwrap();
        retention.observed_withheld_paths = 1;
        retention.observed_not_retained_paths = 1;
        retention.dropped_path_cap = 1;
        retention.retention_truncated = true;
        retention.node_visit_skips = 2;
        let output = collect_test(&mut [diagnostic], &[], false);
        assert_eq!(output.unique_products, 0);
        assert!(output.search_limits_reached);
        assert!(output.retention_truncated);
        assert_eq!(
            output.thresholds[0]
                .retention
                .as_ref()
                .unwrap()
                .dropped_path_cap,
            1
        );
    }

    #[test]
    fn summary_provenance_preserves_sequence_and_labels_header() {
        let mut summary_path = path("ACGTA");
        summary_path.marker_occurrences.clear();
        summary_path.marker_covered_runs.clear();
        let mut diagnostic = threshold(2, vec![summary_path]);
        let retention = diagnostic.withheld_path_diagnostics.as_mut().unwrap();
        retention.retained_summary_only_paths = 1;
        retention.omitted_detailed_marker_occurrences = 1;
        retention.marker_detail_truncated = true;
        let output = collect_test(&mut [diagnostic], &[], false);
        assert_eq!(output.schema_version, 2);
        assert_eq!(output.unique_products, 1);
        assert!(output.marker_detail_truncated);
        assert_eq!(output.products[0].oriented_sequence, "ACGTA");
        assert_eq!(
            output.products[0].sources[0].marker_detail,
            MarkerDetail::SummaryOnly
        );
        assert!(
            output.products[0]
                .description
                .contains("marker_detail=summary_only")
        );
        assert!(
            output.products[0]
                .description
                .contains("assembly_status=uncertain")
        );
        assert!(!output.retention_truncated);
    }

    #[test]
    fn empty_collection_records_explicit_not_evaluated_status_and_limits() {
        let limits = crate::pcr::withheld_diagnostic_limits_for_gene(0, 1);
        let output = UncertainProducts::empty(limits);
        assert_eq!(output.limits, limits);
        assert!(!output.path_search_evaluated);
        assert_eq!(output.unique_products, 0);
        assert!(output.thresholds.is_empty());
    }
}

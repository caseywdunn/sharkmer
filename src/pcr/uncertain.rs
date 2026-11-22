//! Separate, bounded output of unsupported graph-path hypotheses.
//!
//! This consumes the existing diagnostic collector after ordinary admission and
//! deduplication. It never adds a path, extends a search, or selects a threshold.

use crate::kmer::budget::collections::{BudgetMap, BudgetVec, TryClone};
use crate::kmer::budget::text::BudgetString;
use crate::kmer::budget::{MemoryClass, ResourceBudget};
use anyhow::Result;
#[cfg(test)]
use bio::io::fasta;
use serde::Serialize;
use std::sync::Arc;

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

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProductSource {
    pub(crate) threshold: u64,
    pub(crate) marker_detail: MarkerDetail,
    pub(crate) node_visit_limit: usize,
    pub(crate) max_observed_node_visits: usize,
    pub(crate) omitted_self_loop_marker_occurrences: usize,
    pub(crate) cyclic_scc_marker_occurrences: usize,
    pub(crate) retained_collision_marker_occurrences: usize,
    pub(crate) marker_occurrences: BudgetVec<WithheldMarkerOccurrence>,
    pub(crate) marker_covered_runs: BudgetVec<WithheldMarkerCoveredRun>,
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProduct {
    pub(crate) id: BudgetString,
    pub(crate) description: BudgetString,
    pub(crate) assembly_status_causes: BudgetVec<&'static str>,
    pub(crate) oriented_sequence: BudgetString,
    pub(crate) sequence_sha256: BudgetString,
    /// Length of the emitted hypothesis, not an established biological length.
    pub(crate) sequence_length: usize,
    pub(crate) sources: BudgetVec<UncertainProductSource>,
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

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainThresholdSummary {
    /// The standard search counters are preserved, without diagnostic sequences.
    #[serde(flatten)]
    pub(crate) search: PcrThresholdDiagnostic,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) retention: Option<UncertainRetentionSummary>,
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, PartialEq, Eq, Serialize)]
pub(crate) struct UncertainProducts {
    pub(crate) schema_version: u32,
    pub(crate) purpose: &'static str,
    pub(crate) coordinate_system: &'static str,
    pub(crate) marker_span_scope: &'static str,
    pub(crate) read_support_evaluation: &'static str,
    pub(crate) repeat_length_evaluation: &'static str,
    pub(crate) deduplication_policy: &'static str,
    pub(crate) standard_equivalence_policy: &'static str,
    pub(crate) limits: WithheldDiagnosticLimits,
    pub(crate) path_search_evaluated: bool,
    pub(crate) search_limit_scope: &'static str,
    pub(crate) search_limits_reached: bool,
    pub(crate) retention_truncated: bool,
    pub(crate) marker_detail_truncated: bool,
    pub(crate) thresholds: BudgetVec<UncertainThresholdSummary>,
    pub(crate) retained_path_records: usize,
    pub(crate) excluded_standard_equivalent_path_records: usize,
    pub(crate) exact_duplicate_path_records: usize,
    /// Counts oriented sequence hypotheses, not loci, haplotypes, or molecules.
    pub(crate) unique_products: usize,
    pub(crate) products: BudgetVec<UncertainProduct>,
}

impl TryClone for MarkerDetail {
    fn try_clone(&self) -> Result<Self> {
        Ok(*self)
    }
}
impl TryClone for UncertainRetentionSummary {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
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
        })
    }
}
impl TryClone for UncertainThresholdSummary {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            search: TryClone::try_clone(&self.search)?,
            retention: TryClone::try_clone(&self.retention)?,
        })
    }
}
impl TryClone for UncertainProductSource {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            threshold: TryClone::try_clone(&self.threshold)?,
            marker_detail: TryClone::try_clone(&self.marker_detail)?,
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
impl TryClone for UncertainProduct {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            id: TryClone::try_clone(&self.id)?,
            description: TryClone::try_clone(&self.description)?,
            assembly_status_causes: TryClone::try_clone(&self.assembly_status_causes)?,
            oriented_sequence: TryClone::try_clone(&self.oriented_sequence)?,
            sequence_sha256: TryClone::try_clone(&self.sequence_sha256)?,
            sequence_length: TryClone::try_clone(&self.sequence_length)?,
            sources: TryClone::try_clone(&self.sources)?,
        })
    }
}
impl TryClone for UncertainProducts {
    fn try_clone(&self) -> Result<Self> {
        Ok(Self {
            schema_version: TryClone::try_clone(&self.schema_version)?,
            purpose: TryClone::try_clone(&self.purpose)?,
            coordinate_system: TryClone::try_clone(&self.coordinate_system)?,
            marker_span_scope: TryClone::try_clone(&self.marker_span_scope)?,
            read_support_evaluation: TryClone::try_clone(&self.read_support_evaluation)?,
            repeat_length_evaluation: TryClone::try_clone(&self.repeat_length_evaluation)?,
            deduplication_policy: TryClone::try_clone(&self.deduplication_policy)?,
            standard_equivalence_policy: TryClone::try_clone(&self.standard_equivalence_policy)?,
            limits: TryClone::try_clone(&self.limits)?,
            path_search_evaluated: TryClone::try_clone(&self.path_search_evaluated)?,
            search_limit_scope: TryClone::try_clone(&self.search_limit_scope)?,
            search_limits_reached: TryClone::try_clone(&self.search_limits_reached)?,
            retention_truncated: TryClone::try_clone(&self.retention_truncated)?,
            marker_detail_truncated: TryClone::try_clone(&self.marker_detail_truncated)?,
            thresholds: TryClone::try_clone(&self.thresholds)?,
            retained_path_records: TryClone::try_clone(&self.retained_path_records)?,
            excluded_standard_equivalent_path_records: TryClone::try_clone(
                &self.excluded_standard_equivalent_path_records,
            )?,
            exact_duplicate_path_records: TryClone::try_clone(&self.exact_duplicate_path_records)?,
            unique_products: TryClone::try_clone(&self.unique_products)?,
            products: TryClone::try_clone(&self.products)?,
        })
    }
}
impl UncertainProducts {
    pub(super) fn empty(
        limits: WithheldDiagnosticLimits,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Self {
        Self {
            schema_version: 4,
            purpose: "uncertain_repeat_structure_graph_path_hypotheses",
            coordinate_system: "zero_based_half_open_in_emitted_oriented_sequence",
            marker_span_scope: "marker_covered_positions_only_for_complete_detail_sources",
            read_support_evaluation: "not_evaluated",
            repeat_length_evaluation: "unresolved",
            deduplication_policy: "exact_oriented_sequence_only_reverse_complements_retained_as_separate_hypotheses",
            standard_equivalence_policy: "exact_sequence_or_reverse_complement_of_final_standard_output",
            limits,
            path_search_evaluated: false,
            search_limit_scope: "all_attempted_thresholds_for_gene",
            search_limits_reached: false,
            retention_truncated: false,
            marker_detail_truncated: false,
            thresholds: BudgetVec::in_class(budget, MemoryClass::RetainedResults),
            retained_path_records: 0,
            excluded_standard_equivalent_path_records: 0,
            exact_duplicate_path_records: 0,
            unique_products: 0,
            products: BudgetVec::in_class(budget, MemoryClass::RetainedResults),
        }
    }
}

fn canonical_sequence(
    sequence: &[u8],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<BudgetVec<u8>> {
    let reverse = || {
        sequence
            .iter()
            .rev()
            .map(|base| bio::alphabets::dna::complement(*base))
    };
    if reverse().cmp(sequence.iter().copied()).is_lt() {
        BudgetVec::collect_in_class(reverse(), budget, MemoryClass::RetainedResults)
    } else {
        BudgetVec::collect_in_class(
            sequence.iter().copied(),
            budget,
            MemoryClass::RetainedResults,
        )
    }
}

fn source(
    threshold: u64,
    path: WithheldPathRecord,
) -> (BudgetString, BudgetString, UncertainProductSource) {
    let marker_detail = if path.marker_occurrences.is_empty() {
        MarkerDetail::SummaryOnly
    } else {
        MarkerDetail::Complete
    };
    (
        path.oriented_sequence,
        path.sequence_sha256,
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
        },
    )
}

/// Finalize only paths admitted by the unchanged bounded diagnostic collector.
pub(super) fn collect(
    sample: &str,
    gene: &str,
    keep_hidden: bool,
    limits: WithheldDiagnosticLimits,
    standard_records: &[super::sequence::SequenceRecord],
    diagnostics: &mut [PcrThresholdDiagnostic],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<UncertainProducts> {
    let mut output = UncertainProducts::empty(limits, budget);
    let mut standard = BudgetMap::in_class(budget, MemoryClass::RetainedResults);
    for record in standard_records {
        standard.insert(canonical_sequence(record.seq(), budget)?, ())?;
    }
    let mut products: BudgetMap<BudgetString, UncertainProduct> =
        BudgetMap::in_class(budget, MemoryClass::RetainedResults);
    for diagnostic in diagnostics {
        let payload = diagnostic.withheld_path_diagnostics.take();
        let search = diagnostic.try_clone()?;
        if keep_hidden {
            diagnostic.withheld_path_diagnostics = TryClone::try_clone(&payload)?;
        }
        output.path_search_evaluated |= search.path_search_evaluated;
        output.search_limits_reached |=
            search.node_budget_reached || search.dfs_limit_reached || search.path_quota_reached;
        let retention = payload.as_ref().map(UncertainRetentionSummary::from);
        if let Some(retention) = &retention {
            output.search_limits_reached |= retention.node_visit_skips > 0;
            output.retention_truncated |= retention.retention_truncated;
            output.marker_detail_truncated |= retention.marker_detail_truncated;
        }
        output
            .thresholds
            .push(UncertainThresholdSummary { search, retention })?;
        if let Some(payload) = payload {
            for path in payload.paths {
                output.retained_path_records += 1;
                let canonical = canonical_sequence(path.oriented_sequence.as_bytes(), budget)?;
                if standard.get(canonical.as_ref()).is_some() {
                    output.excluded_standard_equivalent_path_records += 1;
                    continue;
                }
                if let Some(product) = products.get_mut(path.oriented_sequence.as_str()) {
                    output.exact_duplicate_path_records += 1;
                    let (_, _, source) = source(diagnostic.threshold, path);
                    product.sources.push(source)?;
                } else {
                    let key = path.oriented_sequence.try_clone()?;
                    let length = path.sequence_length;
                    let (sequence, sha, source) = source(diagnostic.threshold, path);
                    let product = UncertainProduct {
                        id: BudgetString::new(budget, MemoryClass::RetainedResults),
                        description: BudgetString::new(budget, MemoryClass::RetainedResults),
                        assembly_status_causes: BudgetVec::in_class(
                            budget,
                            MemoryClass::RetainedResults,
                        ),
                        oriented_sequence: sequence,
                        sequence_sha256: sha,
                        sequence_length: length,
                        sources: BudgetVec::collect_in_class(
                            [source],
                            budget,
                            MemoryClass::RetainedResults,
                        )?,
                    };
                    products.insert(key, product)?;
                }
            }
        }
    }
    // Full oriented sequence keys retain the old BTreeMap's lexical numbering;
    // the hash table is only a lookup structure, never sequence identity proof.
    let mut ordered = BudgetVec::collect_in_class(products, budget, MemoryClass::RetainedResults)?;
    ordered.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    for (index, (_, mut product)) in ordered.into_iter().enumerate() {
        product.id = BudgetString::format(
            format_args!("{sample}_{gene}_{index}.uncertain"),
            budget,
            MemoryClass::RetainedResults,
        )?;
        for (cause, present) in [
            (
                "high_coverage_repeat_suspected",
                product
                    .sources
                    .iter()
                    .any(|s| s.retained_collision_marker_occurrences > 0),
            ),
            (
                "homopolymer_length_unresolved",
                product
                    .sources
                    .iter()
                    .any(|s| s.omitted_self_loop_marker_occurrences > 0),
            ),
            (
                "repeat_structure_unresolved",
                product
                    .sources
                    .iter()
                    .any(|s| s.cyclic_scc_marker_occurrences > 0),
            ),
        ] {
            if present {
                product.assembly_status_causes.push(cause)?;
            }
        }
        let mut causes = BudgetString::new(budget, MemoryClass::RetainedResults);
        for (index, cause) in product.assembly_status_causes.iter().enumerate() {
            if index > 0 {
                causes.push_str(",")?;
            }
            causes.push_str(cause)?;
        }
        product.description = BudgetString::format(
            format_args!(
                "sample={} target={} sequence_index={} sequence_length={} assembly_status=warning assembly_status_causes={}",
                sample, gene, index, product.sequence_length, causes
            ),
            budget,
            MemoryClass::RetainedResults,
        )?;
        output.products.push(product)?;
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
            .map(|p| p.sources.len())
            .sum::<usize>(),
        output.retained_path_records - output.excluded_standard_equivalent_path_records
    );
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pcr::{WithheldMarkerCause, WithheldMarkerIdentityKind};
    use sha2::{Digest, Sha256};

    fn path(sequence: &str) -> WithheldPathRecord {
        let identity = sequence[..3].to_string();
        WithheldPathRecord {
            oriented_sequence: (sequence.to_string()).into(),
            sequence_sha256: (format!("{:x}", Sha256::digest(sequence.as_bytes()))).into(),
            sequence_length: sequence.len(),
            node_visit_limit: 2,
            max_observed_node_visits: 1,
            omitted_self_loop_marker_occurrences: 0,
            cyclic_scc_marker_occurrences: 1,
            retained_collision_marker_occurrences: 0,
            marker_occurrences: (vec![WithheldMarkerOccurrence {
                cause: WithheldMarkerCause::PrePruneCyclicSccNode,
                identity_kind: WithheldMarkerIdentityKind::OrientedNodeSequence,
                identity_sha256: (format!("{:x}", Sha256::digest(identity.as_bytes()))).into(),
                identity_sequence: (identity).into(),
                start: 0,
                end: 3,
            }])
            .into(),
            marker_covered_runs: (vec![WithheldMarkerCoveredRun { start: 0, end: 3 }]).into(),
        }
    }

    fn threshold(threshold: u64, paths: Vec<WithheldPathRecord>) -> PcrThresholdDiagnostic {
        PcrThresholdDiagnostic {
            threshold,
            connectivity_found: true,
            path_search_evaluated: true,
            completed_candidate_paths: paths.len(),
            withheld_candidate_paths: paths.len(),
            withheld_path_diagnostics: Some(WithheldPathDiagnostics {
                observed_withheld_paths: paths.len(),
                retained_paths: paths.len(),
                paths: paths.into(),
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
            &standard.iter().cloned().map(Into::into).collect::<Vec<_>>(),
            diagnostics,
            None,
        )
        .unwrap()
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
        assert_eq!(
            output.products[0]
                .sources
                .iter()
                .map(|source| source.threshold)
                .collect::<Vec<_>>(),
            vec![8, 4]
        );
        assert!(!output.products[0].description.contains("thresholds="));
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
            assert!(product.description.contains("assembly_status=warning"));
            assert!(
                product
                    .description
                    .contains("assembly_status_causes=repeat_structure_unresolved")
            );
            assert!(!product.description.contains("read_support="));
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
    fn summary_provenance_preserves_sequence_without_exposing_internal_detail_in_header() {
        let mut summary_path = path("ACGTA");
        summary_path.marker_occurrences.clear();
        summary_path.marker_covered_runs.clear();
        let mut diagnostic = threshold(2, vec![summary_path]);
        let retention = diagnostic.withheld_path_diagnostics.as_mut().unwrap();
        retention.retained_summary_only_paths = 1;
        retention.omitted_detailed_marker_occurrences = 1;
        retention.marker_detail_truncated = true;
        let output = collect_test(&mut [diagnostic], &[], false);
        assert_eq!(output.schema_version, 4);
        assert_eq!(output.unique_products, 1);
        assert!(output.marker_detail_truncated);
        assert_eq!(output.products[0].oriented_sequence, "ACGTA");
        assert_eq!(
            output.products[0].sources[0].marker_detail,
            MarkerDetail::SummaryOnly
        );
        assert!(!output.products[0].description.contains("marker_detail="));
        assert!(
            output.products[0]
                .description
                .contains("assembly_status=warning")
        );
        assert!(!output.retention_truncated);
    }

    #[test]
    fn warning_causes_merge_across_sources_with_stable_readable_names() {
        let mut first = path("ACGTA");
        first.omitted_self_loop_marker_occurrences = 1;
        first.retained_collision_marker_occurrences = 1;
        first.marker_occurrences.clear();
        first.marker_covered_runs.clear();
        let mut diagnostics = vec![threshold(8, vec![first]), threshold(4, vec![path("ACGTA")])];
        let output = collect_test(&mut diagnostics, &[], false);
        assert_eq!(output.products.len(), 1);
        let product = &output.products[0];
        assert_eq!(product.id, "sample_gene_0.uncertain");
        assert_eq!(product.oriented_sequence, "ACGTA");
        assert_eq!(product.sources.len(), 2);
        assert_eq!(product.sources[0].omitted_self_loop_marker_occurrences, 1);
        assert_eq!(product.sources[0].retained_collision_marker_occurrences, 1);
        assert_eq!(
            product.sources[1].marker_occurrences,
            path("ACGTA").marker_occurrences
        );
        assert_eq!(
            product.assembly_status_causes,
            vec![
                "high_coverage_repeat_suspected",
                "homopolymer_length_unresolved",
                "repeat_structure_unresolved",
            ]
        );
        assert_eq!(
            product.description,
            "sample=sample target=gene sequence_index=0 sequence_length=5 assembly_status=warning assembly_status_causes=high_coverage_repeat_suspected,homopolymer_length_unresolved,repeat_structure_unresolved"
        );
    }

    #[test]
    fn empty_collection_records_explicit_not_evaluated_status_and_limits() {
        let limits = crate::pcr::withheld_diagnostic_limits_for_gene(0, 1);
        let output = UncertainProducts::empty(limits, None);
        assert_eq!(output.limits, limits);
        assert!(!output.path_search_evaluated);
        assert_eq!(output.unique_products, 0);
        assert!(output.thresholds.is_empty());
    }
    fn warning_budget(bytes: u64) -> Arc<ResourceBudget> {
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

    fn owned_path(sequence: &str, budget: &Arc<ResourceBudget>) -> Result<WithheldPathRecord> {
        let mut template = path(sequence);
        let text =
            |value: &str| BudgetString::from_str(value, Some(budget), MemoryClass::RetainedResults);
        template.oriented_sequence = text(sequence)?;
        template.sequence_sha256 = BudgetString::format(
            format_args!("{:x}", Sha256::digest(sequence.as_bytes())),
            Some(budget),
            MemoryClass::RetainedResults,
        )?;
        let mut markers = BudgetVec::in_class(Some(budget), MemoryClass::RetainedResults);
        for marker in template.marker_occurrences {
            markers.push(WithheldMarkerOccurrence {
                identity_sequence: text(&marker.identity_sequence)?,
                identity_sha256: text(&marker.identity_sha256)?,
                ..marker
            })?;
        }
        template.marker_occurrences = markers;
        template.marker_covered_runs = BudgetVec::collect_in_class(
            template.marker_covered_runs,
            Some(budget),
            MemoryClass::RetainedResults,
        )?;
        Ok(template)
    }

    fn owned_thresholds(budget: &Arc<ResourceBudget>) -> Result<BudgetVec<PcrThresholdDiagnostic>> {
        let mut diagnostics = BudgetVec::in_class(Some(budget), MemoryClass::RetainedResults);
        for (threshold, sequences) in [(10, vec!["AACG", "ACCG"]), (5, vec!["AACG"])] {
            let mut diagnostic = super::tests::threshold(threshold, Vec::new());
            let payload = diagnostic.withheld_path_diagnostics.as_mut().unwrap();
            payload.paths = BudgetVec::in_class(Some(budget), MemoryClass::RetainedResults);
            for sequence in sequences {
                payload.paths.push(owned_path(sequence, budget)?)?;
            }
            payload.retained_paths = payload.paths.len();
            payload.observed_withheld_paths = payload.paths.len();
            diagnostic.completed_candidate_paths = payload.paths.len();
            diagnostic.withheld_candidate_paths = payload.paths.len();
            diagnostics.push(diagnostic)?;
        }
        Ok(diagnostics)
    }

    #[test]
    fn warning_collection_refusal_releases_copies_keys_sources_and_hidden_diagnostics() {
        for hidden in [false, true] {
            let (mut refused, mut succeeded) = (0, 0);
            for bytes in (1..=160).map(|n| n * 128).chain([1 << 20]) {
                let budget = warning_budget(bytes);
                let attempt = (|| -> Result<_> {
                    let mut diagnostics = owned_thresholds(&budget)?;
                    let output = collect(
                        "sample",
                        "gene",
                        hidden,
                        WithheldDiagnosticLimits::default(),
                        &[],
                        &mut diagnostics,
                        Some(&budget),
                    )?;
                    Ok((diagnostics, output))
                })();
                match attempt {
                    Ok((diagnostics, output)) => {
                        assert_eq!(output.retained_path_records, 3);
                        assert_eq!(output.exact_duplicate_path_records, 1);
                        assert_eq!(output.unique_products, 2);
                        assert_eq!(output.products[0].oriented_sequence, "AACG");
                        assert_eq!(output.products[0].sources.len(), 2);
                        assert_eq!(output.products[1].oriented_sequence, "ACCG");
                        assert_eq!(output.products[1].sources.len(), 1);
                        assert_eq!(
                            diagnostics
                                .iter()
                                .all(|d| d.withheld_path_diagnostics.is_some()),
                            hidden
                        );
                        assert_eq!(
                            budget.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
                            0
                        );
                        drop(diagnostics);
                        assert!(budget.snapshot().reserved_memory_bytes > 0);
                        drop(output);
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
                assert_eq!(
                    budget.snapshot().reserved_memory_bytes,
                    0,
                    "hidden {hidden}, budget {bytes}"
                );
            }
            assert!(refused > 10 && succeeded > 0);
        }
    }

    #[test]
    fn canonical_keys_match_byte_oracle_and_warning_hashes_never_replace_sequence_identity() {
        let alphabet = b"ACGTNRYacgt";
        let budget = warning_budget(1 << 20);
        for length in 0..33 {
            for offset in 0..alphabet.len() {
                let sequence: Vec<_> = (0..length)
                    .map(|position| alphabet[(position * 7 + offset) % alphabet.len()])
                    .collect();
                let reverse = bio::alphabets::dna::revcomp(&sequence);
                let expected = sequence.as_slice().min(reverse.as_slice());
                let canonical = canonical_sequence(&sequence, Some(&budget)).unwrap();
                assert_eq!(canonical.as_ref(), expected);
                let state = ahash::RandomState::new();
                assert_eq!(state.hash_one(&canonical), state.hash_one(expected));
                drop(canonical);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
        let mut diagnostics = owned_thresholds(&budget).unwrap();
        for diagnostic in &mut *diagnostics {
            for record in &mut *diagnostic.withheld_path_diagnostics.as_mut().unwrap().paths {
                record.sequence_sha256 = BudgetString::from_str(
                    "same deliberately forged digest",
                    Some(&budget),
                    MemoryClass::RetainedResults,
                )
                .unwrap();
            }
        }
        let output = collect(
            "sample",
            "gene",
            false,
            WithheldDiagnosticLimits::default(),
            &[],
            &mut diagnostics,
            Some(&budget),
        )
        .unwrap();
        assert_eq!(output.unique_products, 2);
        assert_eq!(output.exact_duplicate_path_records, 1);
        drop(diagnostics);
        drop(output);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}

"""Synthetic controls for the independent recount; no biological inputs."""
import copy
import unittest

from independent_outcome_review import historical_matches, retention, CAPS


class RecountControls(unittest.TestCase):
    def entry(self):
        threshold = dict(threshold=18, withheld_candidate_paths=2, path_search_evaluated=True,
                         node_budget_reached=False, dfs_limit_reached=False, path_quota_reached=False)
        usage = dict(paths=2, sequence_bases=10, marker_occurrences=2)
        kept = dict(gene_usage_before_threshold=dict.fromkeys(usage, 0), gene_usage_after_threshold=usage,
                    retained_paths=2, retained_sequence_bases=10, retained_marker_occurrences=2,
                    observed_withheld_paths=2, observed_not_retained_paths=0,
                    dropped_oversize_sequence=0, dropped_oversize_marker_set=0,
                    dropped_path_cap=0, dropped_sequence_base_cap=0, dropped_marker_cap=0,
                    retention_truncated=False, node_visit_limit=2, node_visit_skips=1)
        source = dict(threshold=18, node_visit_limit=2, marker_occurrences=[{}])
        return dict(gene_name="synthetic", n_products=0, threshold_diagnostics=[threshold],
                    uncertain_products=dict(limits=dict(**copy.deepcopy(CAPS),
                    allocated_run_share_for_gene={k: n // 21 + (0 < n % 21) for k, n in CAPS["run"].items()},
                    unused_share_is_not_redistributed=True),
                    thresholds=[dict(**threshold, retention=kept)], retained_path_records=2,
                    path_search_evaluated=True, retention_truncated=False, search_limits_reached=True,
                    excluded_standard_equivalent_path_records=0, exact_duplicate_path_records=1,
                    unique_products=1, products=[dict(sequence_length=5, sources=[source, copy.deepcopy(source)])]))

    def test_exact_reverse_and_near_are_distinct(self):
        rows = [dict(header="exact", sequence="AAAC"), dict(header="reverse", sequence="GTTT"),
                dict(header="near", sequence="AAAT")]
        self.assertEqual(historical_matches("AAAC", rows), dict(exact_oriented=["exact"], reverse_complement_only=["reverse"]))

    def test_palindrome_is_not_double_counted(self):
        self.assertEqual(historical_matches("ACGT", [dict(header="palindrome", sequence="ACGT")]),
                         dict(exact_oriented=["palindrome"], reverse_complement_only=[]))

    def test_retention_positive(self):
        result = retention(self.entry(), 0)
        self.assertEqual((result["retained_path_records"], result["exact_duplicate_path_records"]), (2, 1))

    def test_usage_discontinuity_and_missing_sources_rejected(self):
        for mutation in ["usage", "sources", "drops", "limit_flag"]:
            with self.subTest(mutation=mutation):
                entry = self.entry()
                value = entry["uncertain_products"]
                if mutation == "usage":
                    value["thresholds"][0]["retention"]["gene_usage_after_threshold"]["paths"] = 3
                elif mutation == "sources":
                    value["products"][0]["sources"].pop()
                elif mutation == "drops":
                    value["thresholds"][0]["retention"]["dropped_path_cap"] = 1
                else:
                    value["search_limits_reached"] = False
                with self.assertRaises(ValueError):
                    retention(entry, 0)


if __name__ == "__main__":
    unittest.main()

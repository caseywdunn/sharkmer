"""Independent reviewer controls; synthetic data and authenticated old artifacts only.

These supplement, rather than recount, the implementation author's controls.
They never invoke Sharkmer or read external biological FASTQs.
"""
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "inputs"))
import run_release_comparison as comparison
import extract_historical_identities as extractor


def record(sequence, label="standard", gene="p_g", scope="high_copy_candidate"):
    return comparison.sequence_record({"header": "id assembly_status=" + label,
                                       "sequence": sequence}, gene, label, scope)


def jobs():
    result = []
    for sample in range(10):
        for pair in (1, 2, 3):
            for role in ("baseline", "candidate"):
                rows = [record("AACG"), record("AAGC", gene="p_h", scope="deferred")]
                if role == "candidate":
                    rows = [record("CGTT", "uncertain"), record("AAGC", gene="p_h", scope="deferred"),
                            record("AGCC", "uncertain"), record("AGCC", "uncertain")]
                result.append({"status": "pass", "invocation": {"input": str(sample), "version": role, "pair_index": pair},
                               "outputs": {"records": rows, "counts": {"reads": 7}, "pcr_results": [{"threshold": 2}]},
                               "execution": {"metrics": {"wall_seconds": [1, 5, 9][pair - 1] if role == "baseline" else [2, 4, 6][pair - 1],
                                                           "peak_rss_bytes": 1024 if role == "baseline" else 2048}}})
    return result


class IndependentComparisonControls(unittest.TestCase):
    def test_released_oracle_requires_literal_order_class_multiplicity_and_gene(self):
        expected = {"accession": "synthetic", "panel": "p", "roles": {
            role: {"by_gene": {"g": [{"sequence": "AACG"}, {"sequence": "AAGC"}], "empty": []}}
            for role in ("baseline", "candidate")}}
        original = [record("AACG"), record("AAGC")]
        self.assertTrue(comparison.historical_identity_check(original, expected, "baseline")["ordered_standard_sequences_identical"])
        cases = [original[::-1], original[:1], original + original[:1],
                 [record("CGTT"), original[1]], [record("AACG", "uncertain"), original[1]],
                 [record("AACG", gene="p_other"), original[1]]]
        for changed in cases:
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                comparison.historical_identity_check(changed, expected, "baseline")

    def test_candidate_historical_standard_changes_are_reported_not_rejected(self):
        expected = {"accession": "synthetic", "panel": "p", "roles": {
            "candidate": {"by_gene": {"g": [{"sequence": "AACG"}]}}}}
        changed = [record("AACG", "uncertain"), record("AGCC")]
        result = comparison.historical_identity_check(changed, expected, "candidate")
        self.assertFalse(result["ordered_standard_sequences_identical"])
        self.assertEqual(result["historical_role"], "september_11_pre_inclusive_candidate")
        self.assertEqual(result["missing_literal_identities"], [{"gene": "p_g", "sha256": changed[0]["sha256"]}])
        self.assertEqual(result["gained_literal_identities"], [{"gene": "p_g", "sha256": changed[1]["sha256"]}])

    def test_archive_oracle_is_reproducible_from_authenticated_members(self):
        with tempfile.TemporaryDirectory() as temporary:
            regenerated = extractor.extract(HERE / "inputs/inventory.json", Path(temporary) / "oracle.json")
        recorded = json.loads((HERE / "inputs/historical_identities.json").read_text())
        self.assertEqual(regenerated, recorded)
        self.assertEqual(recorded["totals"]["baseline"], {"all_products": 86, "high_copy_products": 71, "deferred_products": 15})
        self.assertEqual(len(recorded["invocation_evidence"]), 60)

    def test_class_precedence_retains_competing_orientation_evidence(self):
        old = [record("AACG"), record("ACGT")]
        new = [record("AACG", "uncertain"), record("CGTT"), record("ACGT", "uncertain")]
        rows = {row["sha256"]: row for row in comparison.dispositions(old, new)}
        first = rows[old[0]["sha256"]]
        self.assertEqual(first["disposition"], "standard_rc")
        self.assertTrue(first["evidence"]["uncertain_exact"])
        self.assertTrue(first["evidence"]["standard_rc"])
        self.assertEqual(rows[old[1]["sha256"]]["disposition"], "uncertain_exact")
        self.assertEqual(comparison.dispositions([record("AACG")], [record("AACG", gene="other")])[0]["disposition"], "unobserved")

    def test_summary_counts_records_identities_scope_and_paired_changes_separately(self):
        summary = comparison.summarize(jobs())
        self.assertEqual(summary["totals"]["candidate"], {"standard": 10, "uncertain": 30})
        cell = summary["cells"][0]
        self.assertEqual(cell["roles"]["candidate"]["unique_counts"], {"standard": 1, "uncertain": 2})
        self.assertEqual(cell["roles"]["candidate"]["scope_unique_counts"]["deferred"], {"standard": 1, "uncertain": 0})
        self.assertEqual(cell["roles"]["candidate"]["standard_targets"], ["p_h"])
        self.assertEqual(cell["roles"]["candidate"]["available_targets"], ["p_g", "p_h"])
        self.assertEqual(len(cell["literal_gains"]), 2)
        self.assertEqual(sum(row["reverse_complement_in_baseline"] for row in cell["literal_gains"]), 1)
        self.assertEqual([pair["wall_seconds"]["absolute_delta"] for pair in cell["paired_resources"]], [1, -1, -3])
        self.assertEqual(cell["resource_changes"]["wall_seconds"]["absolute_delta"], -1)
        self.assertFalse(cell["resource_changes"]["wall_seconds"]["increase"])
        self.assertTrue(cell["resource_changes"]["peak_rss_bytes"]["priority_review"])

    def test_repeat_state_and_paired_input_counts_cannot_be_silently_selected(self):
        for field in ("counts", "pcr_results", "records"):
            changed = jobs()
            changed[-1]["outputs"][field] = {} if field == "counts" else []
            with self.subTest(field=field), self.assertRaises(ValueError):
                comparison.summarize(changed)
        changed = jobs()
        for job in changed:
            if job["invocation"]["version"] == "candidate":
                job["outputs"]["counts"]["reads"] = 8
        with self.assertRaises(ValueError):
            comparison.summarize(changed)

    def test_failed_or_repeated_pair_cannot_satisfy_sixty_run_total(self):
        changed = jobs()
        changed[-1]["status"] = "failed"
        with self.assertRaises(ValueError):
            comparison.summarize(changed)
        changed = jobs()
        changed[-1]["invocation"]["pair_index"] = 1
        with self.assertRaises(ValueError):
            comparison.summarize(changed)


if __name__ == "__main__":
    unittest.main()

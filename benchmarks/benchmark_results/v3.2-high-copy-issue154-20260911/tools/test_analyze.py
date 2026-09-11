import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SPEC = importlib.util.spec_from_file_location("analyze", Path(__file__).with_name("analyze.py"))
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def product(gene, sequence, status):
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    return {
        "sha256": digest,
        "length": len(sequence),
        "header": f"sample_{gene}_0 gene=panel_{gene} length={len(sequence)} kmer_count_median=6 kmer_count_min=2",
        "reference_match": {"status": status, "reference_id": "reference"},
    }


def result(invocation, genes):
    return {
        "signature": {"invocation": {**invocation, "sample_prefix": "panel_input_1M"}},
        "status": "complete",
        "timing_status": "complete",
        "classification_status": "complete",
        "execution": {
            "returncode": 0,
            "timed_out": False,
            "orphaned_process_group_cleaned": False,
            "measurement_complete": True,
            "measurement_matches_exit": True,
            "input_changed_during_invocation": False,
            "binary_changed_during_invocation": False,
            "panel_changed_during_invocation": False,
            "wall_time_s": 2.0 if invocation["version"] == "candidate" else 1.0,
            "gnu_time": {"exit_status": 0, "peak_rss_bytes": 200 if invocation["version"] == "candidate" else 100},
        },
        "metrics": {field: 10 for field in ANALYZER.COUNT_FIELDS},
        "genes": [{"gene": gene, "n_products": len(products), "products": products} for gene, products in genes.items()],
    }


def fixture(root, omit_result=False):
    protocol = {
        "target_metadata": {"panel": [
            {"gene": "high_lost", "scope": "high_copy_candidate"},
            {"gene": "deferred_lost", "scope": "deferred_or_unclassified"},
            {"gene": "high_gained", "scope": "high_copy_candidate"},
        ]},
        "priority": "high-copy recovery",
        "baseline_role": "baseline",
        "review_criteria": ["descriptive"],
    }
    schedule = []
    for pair_index in range(3):
        for version in ("baseline", "candidate"):
            invocation = {"invocation_id": f"{version}-{pair_index}", "cell": "cell", "panel": "panel", "input": "input", "depth": "1M", "pair_index": pair_index, "version": version}
            schedule.append(invocation)
            if omit_result and invocation["invocation_id"] == "candidate-2":
                continue
            genes = {
                "high_lost": [product("high_lost", "AAAA", "confirmed")],
                "deferred_lost": [product("deferred_lost", "CCCC", "unclassified")],
            } if version == "baseline" else {"high_gained": [product("high_gained", "GGGG", "candidate")]}
            write_json(root / "results" / f"{invocation['invocation_id']}.json", result(invocation, genes))
    receipt_path = root / "receipts" / "protocol.json"
    write_json(receipt_path, protocol)
    write_json(root / "provenance.json", {"protocol": protocol, "frozen_receipts": {"protocol": {"path": str(receipt_path), "sha256": ANALYZER.sha256_file(receipt_path)}}, "schedule": schedule})


class AnalyzeTest(unittest.TestCase):
    def test_finished_issue154_header_uses_panel_gene_prefix(self):
        product_from_finished_result = {
            "sha256": "4421359eb30276c48af559081691c91a74d132155b4192c72cd44eac885f4af4",
            "length": 352,
            "header": "insecta_SRR31887760_1000000_insecta_CO1_1_0 sample=insecta_SRR31887760_1000000 gene=insecta_CO1_1 product=0 length=352 kmer_count_mean=36.97 kmer_count_median=38 kmer_count_min=22 kmer_count_max=47 score=38.00",
            "reference_match": {"status": "insufficient_alignment", "expected_gene": "CO1_1"},
        }
        evidence = ANALYZER.product_evidence("insecta", "CO1_1", product_from_finished_result)
        self.assertEqual(evidence["kmer_count_min"], 22)
        unique, duplicates = ANALYZER.membership([evidence, evidence])
        self.assertEqual(len(unique), 1)
        self.assertEqual(len(duplicates), 1)
        product_from_finished_result["reference_match"] = {"status": "failed_run"}
        with self.assertRaisesRegex(ValueError, "classification"):
            ANALYZER.product_evidence("insecta", "CO1_1", product_from_finished_result)

    def test_complete_three_pairs_separate_scopes(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture(root)
            analysis = ANALYZER.analyze(root)
            cell = analysis["cells"][0]
            self.assertEqual(cell["repetitions"], 3)
            self.assertTrue(all(pair["counts_identical"] for pair in cell["pairs"]))
            self.assertEqual(len(cell["losses_by_scope"]["high_copy_candidate"]), 3)
            self.assertEqual(len(cell["losses_by_scope"]["deferred_or_unclassified"]), 3)
            self.assertEqual(len(cell["all_unfiltered_losses"]), 6)
            self.assertEqual(cell["pairs"][0]["gained_products"][0]["length"], 4)
            self.assertEqual(cell["pairs"][0]["gained_products"][0]["kmer_count_median"], "6")
            self.assertTrue(cell["same_version_stability"]["baseline"]["sequence_membership_stable"])

    def test_incomplete_frozen_schedule_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fixture(root, omit_result=True)
            with self.assertRaisesRegex(ValueError, "exactly match"):
                ANALYZER.analyze(root)


if __name__ == "__main__":
    unittest.main()

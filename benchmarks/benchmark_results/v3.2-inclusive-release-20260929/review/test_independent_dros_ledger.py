"""Independent ledger audit boundary controls; no assembler or biological scan.

Uses the author's explicit synthetic directory fixture and mocked publication
parser only to exercise the independent audit's source/standard linkage.
"""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import yaml

ROOT = Path(__file__).resolve().parents[4]
SCRATCH = ROOT / "tmp/inclusive-release-20260929/dros12s"
sys.path.insert(0, str(SCRATCH))
import audit_ledger
spec = importlib.util.spec_from_file_location("author_ledger_fixture", SCRATCH / "test_audit_ledger.py")
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


class IndependentLedgerControls(unittest.TestCase):
    def test_retained_remainder_requires_actual_standard_or_reverse_complement(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, paths, package = fixtures.LedgerControls().full_fixture(Path(temporary))
            for directory in args[:2]:
                path = directory / f"{args[2]}.stats.yaml"
                stats = yaml.safe_load(path.read_text())
                gene = stats["pcr_results"][0]
                removed = gene["uncertain_products"]["products"].pop(0)
                self.assertEqual(removed["oriented_sequence"], "C" * 495)
                gene["uncertain_products"]["excluded_standard_equivalent_path_records"] = 1
                path.write_text(yaml.safe_dump(stats))
            package.runner.parse_fasta_products = lambda *unused: [{"gene": "insecta_12S", "products": [{"sequence": "G" * 495}]}]
            with patch.dict("sys.modules", {"sharkmer_validate": package}):
                self.assertEqual(audit_ledger.audit(*args)["status"], "pass")
                package.runner.parse_fasta_products = lambda *unused: [{"gene": "insecta_12S", "products": [{"sequence": "A" * 495}]}]
                with self.assertRaisesRegex(ValueError, "not standard-equivalent"):
                    audit_ledger.audit(*args)

    def test_unchanged_ledger_does_not_mask_changed_normal_threshold_assessment(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, paths, package = fixtures.LedgerControls().full_fixture(Path(temporary))
            path = args[1] / f"{args[2]}.stats.yaml"
            stats = yaml.safe_load(path.read_text())
            stats["pcr_results"][0]["uncertain_products"]["thresholds"][0]["extra_search_state"] = "changed"
            path.write_text(yaml.safe_dump(stats))
            with patch.dict("sys.modules", {"sharkmer_validate": package}):
                with self.assertRaisesRegex(ValueError, "Normal result parity mismatch: pcr_results"):
                    audit_ledger.audit(*args)


if __name__ == "__main__":
    unittest.main()

"""Synthetic receipt-merge controls; only tiny temporary records are read."""

import copy
import json
from pathlib import Path
import tempfile
import unittest

import merge_availability as M
import restore_inputs as R
from test_restore_inputs import expected, record


class MergeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.samples = []
        for index in range(10):
            accession = f"SRR{index + 1}"
            path = self.base / f"{accession}.fastq"
            raw = record(accession)
            path.write_bytes(raw)
            self.samples.append({"accession": accession, "panel": "test", "expected_prefix": expected(raw),
                                 "current_file": {"path": str(path)}, "restoration_sources": []})
        self.inventory = self.base / "inventory.json"
        R.write_json(self.inventory, {"samples": self.samples, "panels": []})
        self.attempts = []
        for version in (1, 2):
            directory = self.base / f"v{version}"
            directory.mkdir()
            helper = M.HERE / ("restore_inputs.py" if version == 1 else "restore_inputs_v2.py")
            frozen = {"inventory": R.fingerprint(self.inventory), "helper": R.fingerprint(helper)}
            rows = []
            for index, sample in enumerate(self.samples):
                if version == 2 and index < 8:
                    continue
                if version == 1 and index >= 8:
                    receipt = {"schema_version": 1, "status": "failed", "accession": sample["accession"],
                               "panel_id": "test", "expected_prefix": sample["expected_prefix"],
                               "sources_before": frozen, "available_path": None, "error": "synthetic failure"}
                    R.write_json(directory / f"{sample['accession']}.restore.json", receipt)
                else:
                    receipt = R.restore_one(sample, directory, frozen)
                row = {key: receipt[key] for key in
                       ("accession", "panel_id", "status", "available_path", "expected_prefix")}
                row["receipt"] = R.fingerprint(directory / f"{sample['accession']}.restore.json")
                rows.append(row)
            path = directory / "availability.json"
            R.write_json(path, {"schema_version": 1, "status": "failed" if version == 1 else "pass",
                               "sources_before": frozen, "sources_after": frozen, "samples": rows})
            self.attempts.append(path)

    def merge(self):
        return M.merge(self.inventory, self.attempts, self.base / "merged.json")

    def test_eight_prior_successes_plus_two_retries(self):
        result = self.merge()
        self.assertEqual(result["status"], "pass")
        self.assertEqual(len(result["samples"]), 10)
        self.assertEqual([a["status"] for a in result["source_attempts"]], ["failed", "pass"])
        self.assertEqual(result["sources_before"], result["sources_after"])

    def test_missing_prefix_rejected(self):
        self.attempts = self.attempts[:1]
        with self.assertRaisesRegex(ValueError, "Missing verified"):
            self.merge()

    def test_duplicate_verified_prefix_rejected(self):
        self.attempts.append(self.attempts[1])
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            self.merge()

    def test_changed_receipt_rejected(self):
        receipt = self.attempts[0].parent / "SRR1.restore.json"
        receipt.write_text(receipt.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "Changed source/receipt"):
            self.merge()

    def test_changed_fastq_identity_rejected(self):
        path = Path(self.samples[0]["current_file"]["path"])
        path.write_bytes(path.read_bytes().replace(b"ACGT", b"TGCA"))
        with self.assertRaisesRegex(ValueError, "FASTQ file changed"):
            self.merge()

    def test_no_clobber(self):
        self.merge()
        before = (self.base / "merged.json").read_bytes()
        with self.assertRaisesRegex(ValueError, "exists"):
            self.merge()
        self.assertEqual((self.base / "merged.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()

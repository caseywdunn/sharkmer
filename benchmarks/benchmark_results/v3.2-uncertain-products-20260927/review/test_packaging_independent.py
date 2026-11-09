"""Independent packager assertions using the author's synthetic setup only.

These tests create tiny temporary archives; they never package real results.
The shared fixture is not independent biological or compiled-check evidence.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[4]
SETUP = ROOT / "tmp/uncertain-products-20260927/test_package_evidence.py"
SPEC = importlib.util.spec_from_file_location("packaging_author_fixture", SETUP)
AUTHOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUTHOR)
PACKAGE = AUTHOR.package


class IndependentPackagingChecks(unittest.TestCase):
    def fixture(self, root):
        return AUTHOR.PackagingControls().synthetic_final_evidence(root)

    def test_each_selected_evidence_category_requires_review_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            args, constants = self.fixture(Path(directory))
            review = json.loads(args.final_review.read_text())
            categories = [
                "selected/checks.json", "selected/baseline/build.json",
                "selected/candidate/snapshot.json", "comparison/freeze.json",
                "comparison/runs/assessment.json", "comparison/runs/job3.receipt.json",
                "comparison/runs/job3/output.yaml", "docs/result.md", "bundle/README.md",
            ]
            with patch.multiple(PACKAGE, **constants):
                for key in categories:
                    with self.subTest(key=key):
                        changed = json.loads(json.dumps(review))
                        del changed["reviewed_artifacts"][key]
                        args.final_review.write_bytes(PACKAGE.encode(changed))
                        with self.assertRaisesRegex(ValueError, "omits selected artifacts"):
                            PACKAGE.prepare(args)

    def test_changed_receipted_output_rejected_before_archiving(self):
        with tempfile.TemporaryDirectory() as directory:
            args, constants = self.fixture(Path(directory))
            (args.comparison_dir / "runs/job4/output.yaml").write_bytes(b"changed\n")
            with patch.multiple(PACKAGE, **constants):
                with self.assertRaisesRegex(ValueError, "output changed"):
                    PACKAGE.prepare(args)

    def test_duplicate_jobs_rejected_even_with_refreshed_review_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            args, constants = self.fixture(Path(directory))
            execution = args.comparison_dir / "runs/execution.json"
            value = json.loads(execution.read_text())
            value["jobs"][1]["id"] = value["jobs"][0]["id"]
            execution.write_bytes(PACKAGE.encode(value))
            review = json.loads(args.final_review.read_text())
            review["reviewed_artifacts"]["comparison/runs/execution.json"] = PACKAGE.fingerprint(execution)
            args.final_review.write_bytes(PACKAGE.encode(review))
            with patch.multiple(PACKAGE, **constants):
                with self.assertRaisesRegex(ValueError, "nine unique job receipts"):
                    PACKAGE.prepare(args)

    def test_compiled_candidate_copy_must_match_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            args, constants = self.fixture(Path(directory))
            (args.candidate_build.parent / "source/src/main.rs").write_bytes(b"changed source\n")
            with patch.multiple(PACKAGE, **constants):
                with self.assertRaisesRegex(ValueError, "Candidate source differs"):
                    PACKAGE.prepare(args)

    def test_two_archives_from_same_plan_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.txt"
            source.write_bytes(b"independent synthetic evidence\n")
            plan = {"members": {"evidence/source.txt": PACKAGE.fingerprint(source)}}
            first = PACKAGE.create(plan, root / "archive1")
            second = PACKAGE.create(plan, root / "archive2")
            self.assertEqual(first["archive"]["sha256"], second["archive"]["sha256"])
            self.assertEqual(first["verified_member_count"], 3)
            self.assertEqual(hashlib.sha256((root / "archive1/evidence.tar.gz").read_bytes()).hexdigest(), first["archive"]["sha256"])

    def test_symlinked_output_ancestor_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            real = root / "real"
            real.mkdir()
            (root / "alias").symlink_to(real, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "Symlinked output path"):
                PACKAGE.validate_output({"members": {}}, root / "alias/new_archive")


if __name__ == "__main__":
    unittest.main()

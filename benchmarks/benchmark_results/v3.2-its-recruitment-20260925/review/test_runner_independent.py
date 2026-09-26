#!/usr/bin/env python3
"""Independent runner checks; every executed input and command is synthetic.

Writes a compact, source-bound review receipt beside this script. No biological
data, archive, source_paths(), or production freeze operation is accessed.
"""

import contextlib
import datetime
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
RUNNER = HERE.parent / "run_experiment.py"
SPEC = importlib.util.spec_from_file_location("reviewed_experiment_runner", RUNNER)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def identity(path):
    raw = Path(path).read_bytes()
    return {"size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


class RunnerChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="its-runner-review-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.original_root = runner.ROOT
        runner.ROOT = self.root
        self.addCleanup(setattr, runner, "ROOT", self.original_root)
        self.fastq = self.root / "synthetic.fastq"
        self.fastq.write_text("@synthetic-only\nACGTACGT\n+\nIIIIIIII\n")
        self.source = self.root / "synthetic_source.txt"
        self.source.write_text("synthetic source content\n")
        self.output = self.root / "run"
        self.output.mkdir()
        self.freeze = self.root / "freeze.json"
        self.plan = {
            "schema_version": 1,
            "python": sys.version,
            "root": str(self.root),
            "run_dir": str(self.output),
            "fastq": str(self.fastq),
            "input_identity": identity(self.fastq),
            "sources": {self.source.name: identity(self.source)},
            "timeout_seconds_per_command": 3,
            "address_space_bytes": 256 * 1024 * 1024,
            "commands": [{"id": "synthetic", "argv": [sys.executable, "-c", "print('synthetic success')"]}],
        }

    def command(self, code):
        self.plan["commands"][0]["argv"] = [sys.executable, "-c", code]

    def execute(self, before=None):
        self.freeze.write_text(json.dumps(self.plan))
        frozen_identity = identity(self.freeze)
        if before:
            before()
        with contextlib.redirect_stdout(io.StringIO()):
            result = runner.execute(self.freeze)
        started = json.loads((self.output / "started.json").read_text())
        execution = json.loads((self.output / "execution.json").read_text())
        self.assertEqual(started["freeze"], frozen_identity)
        self.assertEqual(execution["freeze"], frozen_identity)
        self.assertEqual(execution["commands_planned"], [item["id"] for item in self.plan["commands"]])
        self.assertEqual(len(execution["commands"]), 1)
        receipt = execution["commands"][0]
        self.assertIn("started_utc", receipt)
        self.assertIn("finished_utc", receipt)
        self.assertGreaterEqual(receipt["elapsed_seconds"], 0)
        self.assertEqual(execution["status"], receipt["status"])
        if "receipt_publication_error" not in receipt:
            standalone = json.loads((self.output / "synthetic.receipt.json").read_text())
            self.assertEqual(standalone, receipt)
        for name, expected in receipt["outputs"].items():
            self.assertEqual(identity(self.output / name), expected)
        return result, receipt

    def test_success_and_durable_receipts(self):
        result, receipt = self.execute()
        self.assertEqual(result, 0)
        self.assertEqual(receipt["status"], "complete")
        self.assertEqual(receipt["exit_code"], 0)
        self.assertTrue(receipt["launched"])
        self.assertTrue(receipt["source_input_bookends_match"])
        self.assertEqual(receipt["input_before"], identity(self.fastq))
        self.assertEqual(receipt["input_after"], identity(self.fastq))
        self.assertEqual((self.output / "synthetic.stdout.txt").read_text(), "synthetic success\n")
        self.assertIn("synthetic.time.txt", receipt["outputs"])

    def test_nonzero_exit(self):
        self.command("import sys; print('intentional failure'); sys.exit(7)")
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertEqual(receipt["exit_code"], 7)
        self.assertEqual(receipt["status"], "failed")
        self.assertTrue(receipt["source_input_bookends_match"])

    def test_short_timeout(self):
        self.plan["timeout_seconds_per_command"] = 0.1
        self.command("import time; time.sleep(30)")
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertEqual(receipt["exit_code"], 124)
        self.assertTrue(receipt["timeout"])
        self.assertLess(receipt["elapsed_seconds"], 5)
        self.assertTrue(receipt["source_input_bookends_match"])

    def test_missing_executable(self):
        self.plan["commands"][0]["argv"] = [str(self.root / "definitely-not-an-executable")]
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        # GNU time launches successfully and reports its child's exec failure.
        self.assertEqual(receipt["exit_code"], 127)
        self.assertEqual(receipt["status"], "failed")
        self.assertTrue(receipt["launched"])

    def test_input_mutation_after_command(self):
        self.command(f"from pathlib import Path; Path({str(self.fastq)!r}).write_text('changed synthetic input')")
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertEqual(receipt["exit_code"], 0)
        self.assertFalse(receipt["source_input_bookends_match"])
        self.assertNotEqual(receipt["input_before"], receipt["input_after"])
        self.assertIn("Input changed", receipt["verification_error"])

    def test_source_mutation_after_command(self):
        self.command(f"from pathlib import Path; Path({str(self.source)!r}).write_text('changed synthetic source')")
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertEqual(receipt["exit_code"], 0)
        self.assertFalse(receipt["source_input_bookends_match"])
        self.assertIn("Source changed", receipt["verification_error"])

    def test_preverification_mismatch_stops_remaining_commands(self):
        self.plan["commands"].append({"id": "must_not_run", "argv": [sys.executable, "-c", "print('unexpected')"]})
        result, receipt = self.execute(before=lambda: self.fastq.write_text("changed before launch"))
        self.assertEqual(result, 1)
        self.assertFalse(receipt["launched"])
        self.assertIsNone(receipt["exit_code"])
        self.assertIn("Input changed", receipt["error"])
        self.assertFalse(receipt["source_input_bookends_match"])
        self.assertFalse(list(self.output.glob("must_not_run.*")))

    def check_no_clobber(self, suffix):
        existing = self.output / ("synthetic." + suffix)
        existing.write_text("original content must survive\n")
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertFalse(receipt["launched"])
        self.assertIn("Refusing to overwrite", receipt["error"])
        self.assertEqual(existing.read_text(), "original content must survive\n")
        self.assertFalse((self.output / "synthetic.stdout.txt").exists())
        return receipt

    def test_no_clobber_timing(self):
        self.check_no_clobber("time.txt")

    def test_no_clobber_output(self):
        self.check_no_clobber("json")

    def test_no_clobber_receipt_inline_failure_is_durable(self):
        receipt = self.check_no_clobber("receipt.json")
        self.assertIn("receipt_publication_error", receipt)

    def test_resource_limit_inherited(self):
        self.command("import json, resource; print(json.dumps(resource.getrlimit(resource.RLIMIT_AS)))")
        result, receipt = self.execute()
        self.assertEqual(result, 0)
        self.assertEqual(receipt["status"], "complete")
        observed = json.loads((self.output / "synthetic.stdout.txt").read_text())
        self.assertEqual(observed, [self.plan["address_space_bytes"]] * 2)

    def test_launch_setup_error_is_recorded(self):
        self.plan["address_space_bytes"] = "invalid synthetic limit"
        result, receipt = self.execute()
        self.assertEqual(result, 1)
        self.assertFalse(receipt["launched"])
        self.assertEqual(receipt["error_type"], "SubprocessError")
        self.assertIsNone(receipt["exit_code"])

    def test_repeat_execution_preserves_receipts(self):
        result, _ = self.execute()
        self.assertEqual(result, 0)
        before = {path.name: identity(path) for path in self.output.iterdir()}
        with self.assertRaises(FileExistsError):
            runner.execute(self.freeze)
        after = {path.name: identity(path) for path in self.output.iterdir()}
        self.assertEqual(before, after)


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RunnerChecks)
    names = [test.id() for test in suite]
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    payload = {
        "reviewed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "status": "pass" if result.wasSuccessful() else "fail",
        "scope": "Temporary synthetic files only; no biological data scanned or opened.",
        "runner": identity(RUNNER),
        "test_source": identity(Path(__file__)),
        "python": sys.version,
        "tests_run": result.testsRun,
        "test_ids": names,
        "failures": [{"test": test.id(), "traceback": detail} for test, detail in result.failures],
        "errors": [{"test": test.id(), "traceback": detail} for test, detail in result.errors],
        "limitations": "Tests verify published receipts and hashes, not storage recovery after power loss or hostile concurrent filesystem writes.",
    }
    (HERE / "runner_synthetic_results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    sys.exit(0 if result.wasSuccessful() else 1)

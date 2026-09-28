#!/usr/bin/env python3
"""Independent orchestration-gate controls; synthetic receipts/processes only."""

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import run_compiled_checks as checks
import run_allocated_comparison as allocated
from build_snapshot import fingerprint


class CheckReceiptGate(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)
        source = self.output / "source.rs"
        source.write_text("synthetic source identity")
        self.source_rows = {str(source): fingerprint(source)}
        self.runtime = {"interpreter": {"path": str(Path(sys.executable).resolve())}, "python": "synthetic-control"}
        tool_rows = {}
        for name in ("cargo", "rustc", "rustfmt", "cargo-fmt", "cargo-clippy", "clippy-driver"):
            path = self.output / name
            path.write_text(f"synthetic {name} identity")
            tool_rows[name] = fingerprint(path)
        steps = []
        for name, command in checks.command_plan(tool_rows["cargo"]["path"], self.runtime["interpreter"]["path"]):
            log = self.output / f"{name}.log"
            log.write_text("Synthetic receipt fixture: no compiled checks were executed.\n")
            steps.append({"name": name, "command": command, "cwd": str(checks.ROOT),
                          "status": "pass", "returncode": 0, "log": fingerprint(log)})
        self.document = {"status": "pass", "steps": steps, "sources_before": self.source_rows,
                         "sources_after": self.source_rows, "tools_before": tool_rows, "tools_after": tool_rows,
                         "cargo_configs_before": {}, "cargo_configs_after": {},
                         "runtime_before": self.runtime, "runtime_after": self.runtime,
                         "comparison_environment": checks.comparison_environment(os.environ)}

    def verify(self, document):
        path = self.output / "checks.json"
        path.write_text(json.dumps(document))
        with patch.object(checks, "sources", return_value=self.source_rows), \
             patch.object(checks, "runtime", return_value=self.runtime):
            return allocated.verify_checks(path)

    def test_valid_synthetic_gate_control(self):
        self.assertEqual(self.verify(self.document)["status"], "pass")

    def test_every_required_command_cannot_be_replaced_by_success(self):
        for index in range(7):
            with self.subTest(step=index):
                changed = copy.deepcopy(self.document)
                changed["steps"][index]["command"] = ["/usr/bin/true"]
                with self.assertRaises(ValueError):
                    self.verify(changed)

    def test_release_feature_and_all_targets_flags_are_required(self):
        for index, flag in ((1, "--all-targets"), (2, "--release"), (3, "--no-default-features"),
                            (4, "fxhashmap"), (6, "--verify-existing")):
            with self.subTest(index=index, flag=flag):
                changed = copy.deepcopy(self.document)
                changed["steps"][index]["command"].remove(flag)
                with self.assertRaises(ValueError):
                    self.verify(changed)

    def test_required_step_cwd_and_order_are_bound(self):
        for mutation in ("cwd", "order", "status", "returncode"):
            with self.subTest(mutation=mutation):
                changed = copy.deepcopy(self.document)
                if mutation == "order":
                    changed["steps"].reverse()
                else:
                    changed["steps"][0][mutation] = {"cwd": "/other", "status": "failed", "returncode": 1}[mutation]
                with self.assertRaises(ValueError):
                    self.verify(changed)

    def test_runtime_and_shared_compilation_environment_must_match(self):
        for mutation in ("runtime", "environment", "configuration"):
            with self.subTest(mutation=mutation):
                changed = copy.deepcopy(self.document)
                if mutation == "runtime":
                    changed["runtime_before"]["python"] = "different"
                elif mutation == "configuration":
                    changed["cargo_configs_after"] = {"unexpected": {}}
                else:
                    changed["comparison_environment"]["CARGO_BUILD_TARGET"] = "different-target"
                with self.assertRaises(ValueError):
                    self.verify(changed)

    def test_log_and_compiler_bytes_are_rechecked(self):
        for kind in ("log", "compiler"):
            with self.subTest(kind=kind):
                row = self.document["steps"][0]["log"] if kind == "log" else self.document["tools_before"]["rustc"]
                path = Path(row["path"])
                before = path.read_bytes()
                path.write_bytes(before + b"changed")
                with self.assertRaises(ValueError):
                    self.verify(self.document)
                path.write_bytes(before)


class FailureReceiptContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)

    def test_timeout_kills_descendant_and_keeps_receipt(self):
        sentinel = self.output / "descendant_survived"
        ready = self.output / "descendant_ready"
        descendant = "import pathlib,sys,time;pathlib.Path(sys.argv[1]).write_text('ready');time.sleep(0.5);pathlib.Path(sys.argv[2]).write_text('survived')"
        parent = "import subprocess,sys,time;subprocess.Popen([sys.executable,'-c',sys.argv[1],sys.argv[2],sys.argv[3]]);time.sleep(10)"
        with self.assertRaises(RuntimeError):
            allocated.call([sys.executable, "-c", parent, descendant, str(ready), str(sentinel)], self.output, "timeout", 0.2)
        receipt = json.loads((self.output / "timeout.orchestration.json").read_text())
        self.assertTrue(ready.is_file(), "Synthetic descendant did not start")
        self.assertTrue(receipt["process_group_killed"])
        self.assertEqual(receipt["status"], "failed")
        time.sleep(0.6)
        self.assertFalse(sentinel.exists())

    def test_spawn_failure_is_receipted_and_cannot_be_overwritten(self):
        with self.assertRaises(RuntimeError):
            allocated.call([str(self.output / "missing")], self.output, "spawn", 1)
        receipt = self.output / "spawn.orchestration.json"
        before = receipt.read_bytes()
        self.assertEqual(json.loads(before)["status"], "failed")
        with self.assertRaises(FileExistsError):
            allocated.call([sys.executable, "-c", "print('unexpected')"], self.output, "spawn", 1)
        self.assertEqual(receipt.read_bytes(), before)


if __name__ == "__main__":
    unittest.main(verbosity=2)

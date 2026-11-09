#!/usr/bin/env python3
"""Small process/gate controls for job scaffolding, without builds or real reads."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_compiled_checks as checks
import run_allocated_comparison as allocated


class Controls(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)

    def test_successful_step_is_receipted(self):
        result = checks.run_step('ok', [sys.executable, '-c', 'print("synthetic")'], self.output, dict(os.environ), self.output, timeout=3)
        self.assertEqual(result['status'], 'pass')
        self.assertEqual(result['returncode'], 0)
        self.assertEqual(result, json.loads((self.output / 'ok.json').read_text()))

    def test_failure_is_receipted_and_not_overwritten(self):
        command = [sys.executable, '-c', 'raise SystemExit(7)']
        result = checks.run_step('bad', command, self.output, dict(os.environ), self.output, timeout=3)
        self.assertEqual(result['status'], 'failed')
        self.assertEqual(result['returncode'], 7)
        before = (self.output / 'bad.json').read_bytes()
        with self.assertRaises(FileExistsError): checks.run_step('bad', command, self.output, dict(os.environ), self.output, timeout=3)
        self.assertEqual(before, (self.output / 'bad.json').read_bytes())

    def test_timeout_is_receipted(self):
        command = [sys.executable, '-c', 'import time; time.sleep(10)']
        result = checks.run_step('timeout', command, self.output, dict(os.environ), self.output, timeout=0.05)
        self.assertEqual(result['status'], 'failed')
        self.assertIn('TimeoutError', result['error'])
        self.assertLess(result['elapsed_seconds'], 2)

    def test_log_cap_fails_closed(self):
        command = [sys.executable, '-c', 'print("x" * 1000)']
        with patch.object(checks, 'MAX_LOG_BYTES', 64):
            result = checks.run_step('oversize', command, self.output, dict(os.environ), self.output, timeout=3)
        self.assertEqual(result['status'], 'failed')
        self.assertIn('termination threshold', result['error'])

    def test_outer_timeout_and_spawn_failure_have_receipts(self):
        for name, command in [('outer-timeout', [sys.executable, '-c', 'import time; time.sleep(10)']), ('missing-command', [str(self.output / 'nonexistent')])]:
            with self.subTest(name=name), self.assertRaises(RuntimeError):
                allocated.call(command, self.output, name, 0.05)
            receipt = json.loads((self.output / f'{name}.orchestration.json').read_text())
            self.assertEqual(receipt['status'], 'failed')
            self.assertEqual(receipt['command'], command)
            self.assertEqual(receipt['process_group_killed'], name == 'outer-timeout')

    def test_incomplete_checks_gate_rejected(self):
        path = self.output / 'checks.json'
        path.write_text(json.dumps({'status': 'pass', 'steps': []}))
        with self.assertRaises(ValueError): allocated.verify_checks(path)

    def test_check_inventory_prunes_generated_panel_outputs(self):
        root = self.output / 'repository'
        root.mkdir()
        meaningful = ['Cargo.toml', 'Cargo.lock', 'src/main.rs', 'panels/insecta.yaml', 'panels/reference_sources.json.gz', 'panels/examples/example.yaml', 'tests/fixtures/input.fastq', 'scripts/helper.py']
        for name in meaningful + ['panels/validation_runs/old.fastq', 'panels/validation_reports/old.html', 'panels/validation_results/old.json', 'tests/__pycache__/test.pyc']:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
        with patch.object(checks, 'ROOT', root): inventory = checks.sources()
        observed = {str(Path(name).relative_to(root)) for name in inventory if Path(name).is_relative_to(root)}
        self.assertEqual(observed, set(meaningful))

    def test_orchestration_refuses_missing_review_without_launch(self):
        output = self.output / 'attempt'
        argv = ['run_allocated_comparison.py', '--checks', str(self.output / 'checks'), '--baseline-build', str(self.output / 'baseline'), '--candidate-build', str(self.output / 'candidate'), '--review', str(self.output / 'review'), '--attempt-dir', str(output)]
        with patch.object(sys, 'argv', argv), patch.object(allocated, 'call') as launch:
            self.assertEqual(allocated.main(), 1)
            launch.assert_not_called()
        receipt = json.loads((output / 'orchestration.json').read_text())
        self.assertFalse(receipt['comparison_started'])
        self.assertEqual(receipt['status'], 'failed')


if __name__ == '__main__':
    unittest.main(verbosity=2)

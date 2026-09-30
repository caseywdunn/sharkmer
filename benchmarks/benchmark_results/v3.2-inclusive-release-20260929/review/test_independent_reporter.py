"""Independent reporter artifact-binding controls using author synthetic records.

No real invocation, external read, or measurement is performed. The deliberately
shared tiny payloads test fingerprint wiring, not biological capture validity.
"""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(RESULTS))
import report_comparison as reporter
spec = importlib.util.spec_from_file_location("author_reporting_fixtures", RESULTS / "test_report_comparison.py")
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


def prepared(base):
    execution, summary, oracle = fixtures.fixture()
    payload, log = base / "tiny.fasta", base / "tiny.log"
    payload.write_text(">synthetic\nACGA\n")
    log.write_text("synthetic log\n")
    receipts = []
    for job in execution["jobs"]:
        job["outputs"]["files"] = {payload.name: reporter.fingerprint(payload)}
        job["execution"]["files"] = {log.name: reporter.fingerprint(log)}
        path = base / (job["invocation"]["invocation_id"] + ".json")
        path.write_text(json.dumps(job))
        job["receipt"] = reporter.fingerprint(path)
        receipts.append(path)
    paths = [base / name for name in ("execution.json", "summary.json", "oracle.json")]
    paths[1].write_text(json.dumps(summary))
    paths[2].write_text(json.dumps(oracle))
    execution["summary"] = reporter.fingerprint(paths[1])
    execution["bindings"] = {str(paths[2]): reporter.fingerprint(paths[2])}
    paths[0].write_text(json.dumps(execution))
    return paths, payload, log, receipts


class IndependentReporterControls(unittest.TestCase):
    def test_source_bound_outputs_and_timing_logs_reject_postcapture_mutation(self):
        for changed in (None, "payload", "log"):
            with self.subTest(changed=changed), tempfile.TemporaryDirectory() as temporary:
                base = Path(temporary)
                paths, payload, log, _ = prepared(base)
                if changed:
                    (payload if changed == "payload" else log).write_text("mutated\n")
                result = reporter.generate(*paths, base / "report")
                self.assertEqual(result["status"], "pass" if changed is None else "failed")
                if changed:
                    self.assertIn("output/log changed", result["error"])

    def test_execution_cannot_replace_individual_receipt_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            paths, _, _, receipts = prepared(base)
            execution = json.loads(paths[0].read_text())
            execution["jobs"][0]["outputs"]["records"][0]["header"] = "replaced"
            paths[0].write_text(json.dumps(execution))
            result = reporter.generate(*paths, base / "report")
            self.assertEqual(result["status"], "failed")
            self.assertIn("receipt differs from execution", result["error"])

    def test_supplied_oracle_and_summary_must_match_execution_fingerprints(self):
        for index, expected in ((1, "supplied summary"), (2, "supplied historical oracle")):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as temporary:
                base = Path(temporary)
                paths, _, _, _ = prepared(base)
                paths[index].write_text(paths[index].read_text() + "\n")
                result = reporter.generate(*paths, base / "report")
                self.assertEqual(result["status"], "failed")
                self.assertIn(expected, result["error"])


if __name__ == "__main__":
    unittest.main()

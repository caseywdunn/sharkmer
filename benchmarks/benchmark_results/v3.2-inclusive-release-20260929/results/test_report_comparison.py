"""Author controls for posthoc reporting, using synthetic captured outcomes."""

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import report_comparison as R


def record(sequence="ACGA", gene="test_g", label="standard"):
    return {"gene": gene, "assembly_status": label, "scope": "high_copy_candidate", "sequence": sequence,
            "length": len(sequence), "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "reverse_complement_sha256": hashlib.sha256(R.reverse(sequence).encode()).hexdigest(), "header": "synthetic"}


def fixture():
    jobs, cells, expectations = [], [], []
    for index in range(10):
        accession = f"SRR{index + 1}"
        baseline = record()
        label = "uncertain" if index < 7 else "standard"
        candidate = record(label=label)
        disposition = label + "_exact"
        evidence = {name: [baseline["sha256"]] if name == disposition else []
                    for name in ("standard_exact", "standard_rc", "uncertain_exact", "uncertain_rc")}
        cell = {"accession": accession, "roles": {}, "literal_gains": [], "paired_resources": [],
                "released_identity_dispositions": [{"gene": "test_g", "sha256": baseline["sha256"],
                    "length": 4, "scope": "high_copy_candidate", "disposition": disposition, "evidence": evidence}],
                "resource_changes": {}}
        for role, row in (("baseline", baseline), ("candidate", candidate)):
            metrics = {"wall_seconds": [10, 11, 12] if role == "baseline" else [11, 12, 13],
                       "peak_rss_bytes": [1000, 1100, 1200] if role == "baseline" else [990, 1000, 1010]}
            count = {name: int(row["assembly_status"] == name) for name in R.LABELS}
            cell["roles"][role] = {"record_counts": count, "unique_counts": count,
                                  "metrics": {name: {"values": values, "median": values[1], "min": values[0], "max": values[2]}
                                              for name, values in metrics.items()}}
            for pair in (1, 2, 3):
                jobs.append({"status": "pass", "invocation": {"input": accession, "panel": "test", "version": role,
                            "pair_index": pair, "invocation_id": f"{accession}-{role}-{pair}"},
                            "outputs": {"records": [copy.deepcopy(row)], "counts": {"n_reads_read": 1}, "pcr_results": [], "files": {}},
                            "execution": {"status": "pass", "metrics": {name: values[pair - 1] for name, values in metrics.items()}, "files": {}}})
        for metric, priority in (("wall_seconds", 5), ("peak_rss_bytes", 2)):
            before, after = [cell["roles"][role]["metrics"][metric]["median"] for role in ("baseline", "candidate")]
            percent = 100 * (after / before - 1)
            cell["resource_changes"][metric] = {"absolute_delta": after - before, "percent": percent,
                                               "increase": after > before, "priority_review": percent > priority}
        for pair in (1, 2, 3):
            entry = {"pair": pair}
            for metric in ("wall_seconds", "peak_rss_bytes"):
                before, after = [cell["roles"][role]["metrics"][metric]["values"][pair - 1] for role in ("baseline", "candidate")]
                entry[metric] = {"baseline": before, "candidate": after, "absolute_delta": after - before, "ratio": after / before}
            cell["paired_resources"].append(entry)
        old_record = {"sequence": baseline["sequence"], "sequence_sha256": baseline["sha256"], "scope": baseline["scope"]}
        expectations.append({"accession": accession, "panel": "test", "roles": {
            "baseline": {"by_gene": {"g": [old_record]}},
            "candidate": {"by_gene": {"g": [] if index < 7 else [old_record]}}}})
        cells.append(cell)
    return ({"status": "pass", "jobs": jobs},
            {"status": "pass", "cells": cells, "totals": {"baseline": {"standard": 10, "uncertain": 0},
             "candidate": {"standard": 3, "uncertain": 7}}, "median_wall_sum_seconds": {"baseline": 110, "candidate": 120}},
            {"status": "pass", "sample_expectations": expectations})


class ReporterTests(unittest.TestCase):
    def test_duplicate_and_orientation_count_meanings(self):
        rows = [record(), record(), record("TCGT", label="uncertain"), record(gene="another")]
        value = R.counts(rows)
        self.assertEqual(value["total_records"], 4)
        self.assertEqual(value["total_unique_oriented"], 3)
        self.assertEqual(value["total_unique_strand_invariant"], 2)
        self.assertEqual(value["standard_unique_oriented"], 2)

    def test_disposition_priority_and_gains(self):
        baseline = [record(gene=f"g{index}") for index in range(5)]
        candidate = [record(gene="g0"), record("TCGT", "g1"),
                     record(gene="g2", label="uncertain"), record("TCGT", "g3", "uncertain")]
        dispositions, gains = R.identities(baseline, candidate)
        self.assertEqual([row["disposition"] for row in dispositions],
                         ["standard_exact", "standard_rc", "uncertain_exact", "uncertain_rc", "unobserved"])
        self.assertEqual(len(gains), 2)
        self.assertTrue(all(row["reverse_complement_in_baseline"] for row in gains))

    def test_same_literal_in_both_classes_keeps_standard_priority(self):
        dispositions, _ = R.identities([record()], [record(label="uncertain"), record()])
        self.assertEqual(dispositions[0]["disposition"], "standard_exact")
        self.assertTrue(dispositions[0]["evidence"]["uncertain_exact"])

    def test_resources_keep_paired_ratios_separate_from_medians(self):
        report = R.resources([10, 20, 30], [11, 18, 36], 5)
        self.assertEqual(report["paired_ratios"], [1.1, 0.9, 1.2])
        self.assertEqual(report["median_ratio"], 0.9)
        self.assertFalse(report["increase"])
        self.assertTrue(R.resources([10, 10, 10], [11, 11, 11], 5)["priority_review"])

    def test_zero_baseline_resources_preserve_undefined_relative_change(self):
        for after in ([0, 0, 0], [1, 2, 3]):
            report = R.resources([0, 0, 0], after, 5)
            self.assertIsNone(report["median_ratio"])
            self.assertIsNone(report["median_percent_change"])
            self.assertEqual(report["paired_ratios"], [None, None, None])
            self.assertFalse(report["priority_review"])
            self.assertEqual(report["increase"], bool(after[1]))

    def test_changed_historical_release_or_wall_sum_rejected(self):
        execution, summary, oracle = fixture()
        summary["median_wall_sum_seconds"]["baseline"] += 1
        with self.assertRaisesRegex(ValueError, "sum of sample median"):
            R.derive(execution, summary, oracle)
        execution, summary, oracle = fixture()
        oracle["sample_expectations"][0]["roles"]["baseline"]["by_gene"]["g"][0]["sequence"] = "TTTT"
        with self.assertRaisesRegex(ValueError, "ordered identities"):
            R.derive(execution, summary, oracle)

    def test_complete_synthetic_matrix_and_historical_subset(self):
        tables = R.derive(*fixture())
        self.assertEqual(len(tables["sample_counts"]), 60)
        self.assertEqual(len(tables["resources"]), 20)
        self.assertEqual(len(tables["paired_resources"]), 60)
        self.assertEqual(len(tables["historical_seven"]), 7)
        self.assertEqual({row["disposition"] for row in tables["historical_seven"]}, {"uncertain_exact"})

    def test_incomplete_or_duplicate_invocations_rejected(self):
        for mutation in (lambda jobs: jobs.pop(), lambda jobs: jobs.__setitem__(1, copy.deepcopy(jobs[0]))):
            execution, summary, oracle = fixture()
            mutation(execution["jobs"])
            with self.assertRaises(ValueError):
                R.derive(execution, summary, oracle)

    def test_changed_capture_or_summary_rejected(self):
        for mutation in (
            lambda execution, summary: execution["jobs"][0]["outputs"]["records"][0].update(sha256="0" * 64),
            lambda execution, summary: summary["cells"][0]["roles"]["candidate"]["record_counts"].update(uncertain=5),
            lambda execution, summary: summary["cells"][0]["released_identity_dispositions"][0].update(disposition="unobserved"),
            lambda execution, summary: summary["cells"][0]["resource_changes"]["wall_seconds"].update(priority_review=False),
        ):
            execution, summary, oracle = fixture()
            mutation(execution, summary)
            with self.assertRaises(ValueError):
                R.derive(execution, summary, oracle)

    def test_artifact_generation_and_no_clobber(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            execution, summary, oracle = fixture()
            for job in execution["jobs"]:
                path = base / (job["invocation"]["invocation_id"] + ".json")
                path.write_text(json.dumps(job))
                job["receipt"] = R.fingerprint(path)
            oracle_path, summary_path, execution_path = [base / name for name in ("oracle.json", "summary.json", "execution.json")]
            oracle_path.write_text(json.dumps(oracle))
            summary_path.write_text(json.dumps(summary))
            execution["summary"] = R.fingerprint(summary_path)
            execution["bindings"] = {str(oracle_path): R.fingerprint(oracle_path)}
            execution_path.write_text(json.dumps(execution))
            result = R.generate(execution_path, summary_path, oracle_path, base / "reports")
            self.assertEqual(result["status"], "pass", result)
            self.assertEqual(len(result["files"]), 6)
            self.assertEqual(result["sources_before"], result["sources_after"])
            with self.assertRaises(FileExistsError):
                R.generate(execution_path, summary_path, oracle_path, base / "reports")


if __name__ == "__main__":
    unittest.main()

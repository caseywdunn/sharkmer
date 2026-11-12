"""Protect the frozen engineering comparison's label and performance semantics."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "count_comparison", ROOT / "benchmarks/compare_count_backends.py"
)
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


class CountComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.protocol = json.loads(
            (
                ROOT
                / "benchmarks/benchmark_results/v4-performance-baseline-20261109/protocol.json"
            ).read_text()
        )

    def runs(self, ratio: float = 0.94, rss_ratio: float = 1.0) -> list[dict]:
        rows = []
        for sample in self.protocol["samples"]:
            for _ in range(3):
                for backend in ["baseline", "candidate"]:
                    rows.append(
                        {
                            "accession": sample["accession"],
                            "kind": "timed",
                            "backend": backend,
                            "metrics": {
                                "wall_s": 100
                                * (ratio if backend == "candidate" else 1),
                                "rss_kib": (1 << 20)
                                * (rss_ratio if backend == "candidate" else 1),
                            },
                            "allocator_peak_bytes": 1 << 30,
                            "stages": {"count_finalize_s": 0.1},
                        }
                    )
        return rows

    def test_historical_uncertainty_labels_normalize_without_erasing_warning_causes(
        self,
    ) -> None:
        sample = next(
            row for row in self.protocol["samples"] if row["accession"] == "SRR31887760"
        )
        expected = comparison.expected_records(sample)
        self.assertTrue(
            any(
                row["assembly_status"] == "uncertain"
                for row in sample["expected_records"]
            )
        )
        self.assertTrue(any(row["assembly_status"] == "warning" for row in expected))
        self.assertFalse(any(row["assembly_status"] == "uncertain" for row in expected))
        for actual, old in zip(
            expected,
            sorted(sample["expected_records"], key=lambda row: row["id"]),
            strict=True,
        ):
            self.assertEqual(
                {
                    key: value
                    for key, value in actual.items()
                    if key != "assembly_status"
                },
                {key: value for key, value in old.items() if key != "assembly_status"},
            )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "record.fasta").write_text(
                ">id sample=s gene=g length=3 assembly_status=warning assembly_status_causes=repeat_count_unresolved\nAAA\n"
            )
            warning = comparison.records(path)[0]
            self.assertEqual(warning["assembly_status"], "warning")
            self.assertEqual(
                warning["assembly_status_causes"], ["repeat_count_unresolved"]
            )

    def test_no_speedup_fails_aggregate_gate_even_when_individual_limits_pass(
        self,
    ) -> None:
        summary = comparison.summarize(self.runs(ratio=1.0), self.protocol)
        self.assertTrue(
            all(row["wall_passed"] and row["rss_passed"] for row in summary["samples"])
        )
        self.assertFalse(summary["performance_gates_passed"])
        self.assertTrue(
            comparison.summarize(self.runs(), self.protocol)["performance_gates_passed"]
        )

    def test_single_sample_wall_or_rss_regression_fails_despite_aggregate_improvement(
        self,
    ) -> None:
        for metric, multiplier in [("wall_s", 1.1), ("rss_kib", 1.03)]:
            rows = self.runs(ratio=0.8)
            for row in rows:
                if (
                    row["backend"] == "candidate"
                    and row["accession"] == self.protocol["samples"][0]["accession"]
                ):
                    row["metrics"][metric] = (
                        100 if metric == "wall_s" else 1 << 20
                    ) * multiplier
            self.assertFalse(
                comparison.summarize(rows, self.protocol)["performance_gates_passed"]
            )

    def test_warmups_are_excluded_from_timed_medians(self) -> None:
        rows = self.runs()
        expected = comparison.summarize(rows, self.protocol)
        rows += [
            {**row, "kind": "warmup", "metrics": {**row["metrics"], "wall_s": 10000}}
            for row in rows
        ]
        self.assertEqual(comparison.summarize(rows, self.protocol), expected)


if __name__ == "__main__":
    unittest.main()

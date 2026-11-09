#!/usr/bin/env python3
"""Independent, offline harness checks; no real read scan or assembler run.

These target capture invariants separately from the production metadata
validator. Its validation calls are mocked only in the capture unit controls;
the independent parser suite exercises real metadata and manifest validation.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import yaml

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))
import run_comparison as harness
import build_snapshot as builder
from sharkmer_validate import runner


class CaptureContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)
        self.accession = "SRR27962769"
        self.prefix = f"uncertain_{self.accession}"
        self.basename = f"{self.prefix}_insecta_test.fasta"
        (self.output / self.basename).write_text(">standard product=0\nACGT\n")
        (self.output / f"{self.prefix}.stats.yaml").write_text("synthetic: true\n")
        self.sample = {"path": "/synthetic/reads.fastq", "expected_prefix": {"records": 1, "bases": 4}}
        self.job = {"accession": self.accession, "role": "baseline_off",
                    "output": str(self.output), "argv": ["synthetic-assembler"]}
        self.stats = {
            "sample": self.prefix,
            "input_source": {"kind": "local_files", "inputs": [self.sample["path"]],
                             "paired": False, "max_reads": 1000000},
            "chunks": 0, "n_reads_read": 1, "n_bases_read": 4,
            "stage_timings": {"read_threading_input_s": None},
            "pcr_results": [{"gene_name": "insecta_test", "status": "success", "n_products": 1,
                             "product_lengths": [4], "output_file": self.basename,
                             "threshold_diagnostics": []}],
        }

    def capture(self):
        with patch.object(runner, "_parse_stats_yaml", return_value=self.stats), \
             patch.object(runner, "load_panel_yaml", return_value={"name": "insecta"}), \
             patch.object(runner, "panel_gene_names", return_value=["test"]), \
             patch.object(runner, "_validate_stats_manifest"), \
             patch.object(runner, "_validate_output_manifest"):
            return harness.capture(self.job, self.sample)

    def test_valid_capture_control(self):
        self.assertEqual(self.capture()["standard"][0]["sequence"], "ACGT")

    def test_per_gene_declared_count_and_lengths_are_not_replaced_by_total_comparison(self):
        for field, value in (("n_products", 2), ("product_lengths", [5]), ("product_lengths", [4, 4])):
            with self.subTest(field=field, value=value):
                original = copy.deepcopy(self.stats)
                self.stats["pcr_results"][0][field] = value
                with self.assertRaises(ValueError):
                    self.capture()
                self.stats = original

    def test_failed_gene_cannot_publish_standard_file(self):
        self.stats["pcr_results"][0].update(status="fail", n_products=0, product_lengths=[])
        with self.assertRaises(ValueError):
            self.capture()

    def test_extra_record_does_not_match_declared_count(self):
        with (self.output / self.basename).open("a") as output:
            output.write(">second product=1\nACGA\n")
        with self.assertRaises(ValueError):
            self.capture()

    def test_threading_observation_must_be_disabled(self):
        self.stats["stage_timings"]["read_threading_input_s"] = 0.0
        with self.assertRaises(ValueError):
            self.capture()

    def test_hidden_diagnostic_payload_is_forbidden_in_every_mode(self):
        self.stats["pcr_results"][0]["threshold_diagnostics"] = [{"withheld_path_diagnostics": {}}]
        with self.assertRaises(ValueError):
            self.capture()


class ParityContract(unittest.TestCase):
    def fixture(self):
        rows = []
        for accession, roles in harness.ORDER:
            for role in roles:
                gene = {"gene_name": "insecta_test", "status": "success", "n_products": 2,
                        "product_lengths": [4, 4], "output_file": "standard.fasta",
                        "threshold_diagnostics": [{"threshold": 2, "failure_reason": None}]}
                if role == "changed_on":
                    gene["uncertain_products"] = {"unique_products": 0, "products": []}
                rows.append({"accession": accession, "role": role, "outputs": {
                    "stats": {**{key: 1 for key in harness.COUNTS}, "pcr_results": [gene]},
                    "standard": [{"gene": "insecta_test", "header": "one", "sequence": "ACGT"},
                                 {"gene": "insecta_test", "header": "two", "sequence": "ACGA"}],
                    "uncertain": [],
                }})
        return rows

    def test_record_order_and_headers_are_part_of_standard_parity(self):
        for mutation in ("order", "header", "diagnostic", "lengths"):
            with self.subTest(mutation=mutation):
                rows = self.fixture()
                changed = rows[1]["outputs"]
                if mutation == "order":
                    changed["standard"].reverse()
                elif mutation == "header":
                    changed["standard"][0]["header"] += " modified=true"
                elif mutation == "lengths":
                    changed["stats"]["pcr_results"][0]["product_lengths"] = [3, 5]
                else:
                    changed["stats"]["pcr_results"][0]["threshold_diagnostics"][0]["failure_reason"] = "changed"
                self.assertEqual(harness.analyze(rows, [])["status"], "fail")

    def test_new_gene_assessment_alone_does_not_change_standard_parity(self):
        self.assertEqual(harness.analyze(self.fixture(), [])["status"], "pass")

    def test_optional_kmer_count_presence_and_value_are_part_of_parity(self):
        for key in ("n_multi_kmers", "n_singleton_kmers"):
            with self.subTest(key=key):
                rows = self.fixture()
                rows[1]["outputs"]["stats"][key] = 1
                self.assertEqual(harness.analyze(rows, [])["status"], "fail")
                for row in rows:
                    row["outputs"]["stats"][key] = 1
                self.assertEqual(harness.analyze(rows, [])["status"], "pass")
                rows[1]["outputs"]["stats"][key] = 2
                self.assertEqual(harness.analyze(rows, [])["status"], "fail")


class ProcessAndBuildContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)

    def test_hashing_detects_a_write_during_read(self):
        path = self.output / "source"
        path.write_bytes(b"before")
        original = hashlib.sha256

        class MutatingHash:
            def __init__(self):
                self.inner = original()
                self.changed = False

            def update(self, block):
                self.inner.update(block)
                if not self.changed:
                    self.changed = True
                    path.write_bytes(b"after changed")

            def hexdigest(self):
                return self.inner.hexdigest()

        with patch.object(builder.hashlib, "sha256", MutatingHash):
            with self.assertRaises(ValueError):
                builder.fingerprint(path)

    def test_source_symlink_is_not_a_frozen_regular_file(self):
        path = self.output / "source"
        path.write_text("original")
        link = self.output / "source_link"
        link.symlink_to(path)
        with self.assertRaises(ValueError):
            builder.fingerprint(link)

    def test_missing_compiler_probe_has_immutable_failure_receipt(self):
        source = self.output / "source"
        source.mkdir()
        snapshot = {"source": str(source), "files": {}, "exporter": builder.fingerprint(Path(builder.__file__))}
        (self.output / "snapshot.json").write_text(json.dumps(snapshot))
        with self.assertRaises(Exception):
            builder.build(self.output, self.output / "missing_cargo", self.output / "missing_rustc")
        receipt = self.output / "build.json"
        before = receipt.read_bytes()
        value = json.loads(before)
        self.assertEqual(value["status"], "failed")
        self.assertIn("error", value)
        self.assertEqual(value["source_before"], value["source_after"])
        with self.assertRaises(ValueError):
            builder.build(self.output, self.output / "missing_cargo", self.output / "missing_rustc")
        self.assertEqual(receipt.read_bytes(), before)

    def launch(self, wrong=None):
        cpus = sorted(os.sched_getaffinity(0))[:2]
        if not cpus:
            self.skipTest("No CPU available for synthetic enforcement control")
        cap = 512 * 1024 ** 2
        receipt = self.output / "enforcement.json"
        child = self.output / "child.json"
        code = "import json,os,resource,sys;open(sys.argv[1],'w').write(json.dumps({'pid':os.getpid(),'cpus':sorted(os.sched_getaffinity(0)),'as':list(resource.getrlimit(resource.RLIMIT_AS))}))"
        command = [sys.executable, "-c", code, str(child)]
        # Exercise actual process enforcement within this test's allocation;
        # the comparison-plan gate separately requires exactly two CPUs.
        expected_cpus = cpus if wrong != "cpus" else [cpus[0] + 1]
        argv = ["/usr/bin/prlimit", f"--as={cap}", "--", "/usr/bin/taskset", "-c", ",".join(map(str, cpus)),
                sys.executable, str(HERE / "launch_bounded.py"), "--receipt", str(receipt),
                "--expected-cpus", ",".join(map(str, expected_cpus)), "--expected-as", str(cap + int(wrong == "as")),
                "--", *command]
        result = subprocess.run(argv, capture_output=True, text=True, timeout=10)
        return result, receipt, child, command

    def test_observed_constraints_survive_exec_and_pid_is_preserved(self):
        result, receipt, child, command = self.launch()
        self.assertEqual(result.returncode, 0, result.stderr)
        observed, final = json.loads(receipt.read_text()), json.loads(child.read_text())
        self.assertEqual(observed["status"], "pass")
        self.assertEqual(observed["command"], command)
        self.assertEqual(observed["pid"], final["pid"])
        self.assertEqual(observed["affinity"], final["cpus"])
        self.assertEqual(observed["rlimit_as_soft_hard"], final["as"])

    def test_wrong_observed_constraints_prevent_exec(self):
        for wrong in ("cpus", "as"):
            with self.subTest(wrong=wrong):
                result, receipt, child, _ = self.launch(wrong)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(json.loads(receipt.read_text())["status"], "fail")
                self.assertFalse(child.exists())
                receipt.unlink()


class FrozenGateContract(unittest.TestCase):
    """Exercise freeze-shape/review coverage independently of building Rust."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)
        self.required = self.output / "required.py"
        self.extra = self.output / "extra_control.py"
        self.required.write_text("required source")
        self.extra.write_text("extra reviewed control")
        self.tools = {}
        for name in ("prlimit", "time", "taskset"):
            path = self.output / name
            path.write_text(name)
            self.tools[name] = path
        self.review = self.output / "review.json"
        reviewed = {"required": harness.fingerprint(self.required), "extra": harness.fingerprint(self.extra)}
        self.review.write_text(json.dumps({"status": "pass", "reviewed_sources": reviewed}))
        builds, receipts = {}, {}
        for role in ("baseline", "candidate"):
            builds[role] = {"binary": {"path": str(self.tools["time"])},
                            **{key: None for key in ("cargo", "rustc", "cargo_version", "rustc_version",
                                                   "compiler_sysroot", "comparison_environment", "cargo_configs")}}
            receipt = self.output / f"{role}.build.json"
            receipt.write_text(json.dumps(builds[role]))
            receipts[role] = harness.fingerprint(receipt)
        paths = [self.required, self.extra, self.review, *self.tools.values(), Path(sys.executable).resolve(),
                 *(Path(value["path"]) for value in receipts.values())]
        samples = {accession: {"path": str(self.output / accession)} for accession, _ in harness.ORDER}
        self.plan = {
            "schema_version": 1, "status": "frozen_before_execution", "samples": samples,
            "builds": builds, "build_receipts": receipts, "cpus": [7, 8],
            "runtime": {"python": sys.version, "interpreter": harness.fingerprint(Path(sys.executable).resolve()),
                        "yaml_version": yaml.__version__},
            "execution_environment": {"LC_ALL": "C"}, "review": harness.fingerprint(self.review),
            "reviewed_sources": reviewed, "sources": {str(path): harness.fingerprint(path) for path in paths},
            "limits": dict(harness.LIMITS), "output": str(self.output / "run"),
        }
        with patch.object(harness, "TOOLS", self.tools):
            self.plan["jobs"] = harness.commands(samples, builds, Path(self.plan["output"]), self.plan["cpus"])

    def check(self, plan):
        def build_info(path, role):
            return self.plan["builds"][role], {"files": {}}
        with patch.object(harness, "TOOLS", self.tools), \
             patch.object(harness, "build_info", side_effect=build_info), \
             patch.object(harness, "verify_sample_plan"), \
             patch.object(harness, "required_sources", return_value={self.required}), \
             patch.object(harness.os, "sched_getaffinity", return_value={7, 8}):
            harness.check_frozen(plan)

    def test_valid_review_coverage_control(self):
        self.check(self.plan)

    def test_source_coverage_cannot_omit_extra_reviewed_controls(self):
        for path in (self.required, self.extra):
            with self.subTest(path=path):
                plan = copy.deepcopy(self.plan)
                del plan["sources"][str(path)]
                with self.assertRaises(ValueError):
                    self.check(plan)

    def test_prepared_prefix_metadata_and_absolute_path_are_required(self):
        inputs = self.output / "inputs"
        inputs.mkdir()
        original = {"accession": "synthetic", "path": "historical_path", "expected_prefix": {
            "records": 1, "bases": 4, "size_bytes": 15, "sha256": "0" * 64}, "roles": ["calibration"]}
        (inputs / "read_inputs.json").write_text(json.dumps({"samples": [original]}))
        sample = {**original, "path": str(self.output / "reads"), "observed_prefix": original["expected_prefix"]}
        with patch.object(harness, "HERE", self.output):
            harness.verify_sample_plan({"synthetic": sample})
            for field, value in (("path", "relative_reads"), ("observed_prefix", {}),
                                 ("expected_prefix", {}), ("roles", ["validation"])):
                with self.subTest(field=field):
                    changed = copy.deepcopy(sample)
                    changed[field] = value
                    with self.assertRaises(ValueError):
                        harness.verify_sample_plan({"synthetic": changed})

    def test_changed_extra_reviewed_control_invalidates_plan(self):
        self.extra.write_text("changed control")
        with self.assertRaises(ValueError):
            self.check(self.plan)

    def test_invalid_frozen_shape_and_cpu_values_are_rejected(self):
        mutations = (("schema_version", True), ("schema_version", 2), ("status", "complete"),
                     ("sources", {}), ("reviewed_sources", {}), ("cpus", [7]),
                     ("cpus", [7, 7]), ("cpus", [True, 8]), ("cpus", [7, 9]),
                     ("samples", {}), ("build_receipts", {}), ("execution_environment", {}))
        for field, value in mutations:
            with self.subTest(field=field, value=value):
                plan = copy.deepcopy(self.plan)
                plan[field] = value
                with self.assertRaises(ValueError):
                    self.check(plan)

    def test_job_order_and_flags_cannot_change_after_freeze(self):
        for mutation in ("order", "flag", "time", "output"):
            with self.subTest(mutation=mutation):
                plan = copy.deepcopy(self.plan)
                if mutation == "order":
                    plan["jobs"].reverse()
                elif mutation == "flag":
                    plan["jobs"][0]["argv"].insert(-1, "--emit-uncertain")
                elif mutation == "time":
                    plan["limits"]["timeout_seconds"] += 1
                else:
                    plan["jobs"][0]["output"] += ".changed"
                with self.assertRaises(ValueError):
                    self.check(plan)


if __name__ == "__main__":
    unittest.main(verbosity=2)

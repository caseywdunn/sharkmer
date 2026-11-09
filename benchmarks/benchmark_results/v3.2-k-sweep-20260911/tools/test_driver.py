import copy
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parent
DRIVER_PATH = ROOT / "driver.py"
HIGH_COPY_PROTOCOL = Path("/tmp/sharkmer-high-copy-20260911/final-protocol.json")


def load_driver():
    specification = importlib.util.spec_from_file_location("k_sweep_driver", DRIVER_PATH)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def protocol_fixture(driver):
    source = json.loads(HIGH_COPY_PROTOCOL.read_text())
    samples = [sample for sample in source["samples"] if sample["panel"] == "insecta"]
    return {
        "schema_version": 1,
        "purpose": "Test preregistered k sweep",
        "allowed_k": [19, 23, 27, 31],
        "settings": {key: value for key, value in source["settings"].items() if key != "k"},
        "ordering": copy.deepcopy(driver.EXPECTED_ORDERING),
        "baseline_role": source["baseline_role"],
        "expected_commits": source["expected_commits"],
        "samples": samples,
        "target_metadata": {"insecta": source["target_metadata"]["insecta"]},
        "fixed_reference": {
            "execution": "/tmp/sharkmer-high-copy-20260911/final-execution",
            "analysis": "/tmp/sharkmer-high-copy-20260911/final-analysis",
            "protocol_sha256": driver.EXPECTED_DEPENDENCIES["high_copy_protocol"][1],
            "diagnostic_audit_sha256": "fcfee505da9fe0d7e2ef464546c530db0890adf6dbe8ca7c5bb4d0f5f9653194",
            "baseline_k": 19,
        },
        "limitations": ["Calibration only."],
    }


class FakeHelper:
    def __init__(self):
        self.validate_shared_stats = self.original_shared
        self.validate_legacy_outputs = lambda *arguments: ([], {})
        self.validate_current_outputs = self.original_current

    def original_shared(self, stats, invocation, panel, input_record, command, expected_version):
        if stats["kmer_length"] != 19:
            raise ValueError("original helper expected k19")
        return ["validated"]

    def original_current(self, output_dir, invocation, panel, stats, pcr_results, command, runner):
        if stats["kmer_length"] != 19:
            raise ValueError("original current validator expected k19")
        return ["genes"], {"complete": True}


class FakeRunner:
    def __init__(self):
        self.observed_k = None

    def _validate_stats_manifest(self, stats, sample, selected_k, genes, mode, expected_command):
        self.observed_k = selected_k
        if stats["kmer_length"] != selected_k:
            raise ValueError("manifest k mismatch")


class DriverTests(unittest.TestCase):
    def setUp(self):
        self.driver = load_driver()

    def test_protocol_allows_only_preregistered_k_and_exact_invariants(self):
        protocol = protocol_fixture(self.driver)
        self.assertEqual(self.driver.validate_protocol(protocol), [19, 23, 27, 31])
        invalid = copy.deepcopy(protocol)
        invalid["allowed_k"] = [19, 25]
        with self.assertRaisesRegex(ValueError, "allowed_k"):
            self.driver.validate_protocol(invalid)
        invalid = copy.deepcopy(protocol)
        invalid["settings"]["threads"] = 3
        with self.assertRaisesRegex(ValueError, "Non-k settings"):
            self.driver.validate_protocol(invalid)

    def test_registered_protocol_and_dependency_receipt_names_are_valid(self):
        registered = self.driver.load_json(ROOT / "discovery-protocol.json")
        self.assertEqual(self.driver.validate_protocol(registered), [19, 23, 27, 31])
        destinations = [
            f"dependency-{label}{source.suffix or '.receipt'}"
            for label, (source, unused_sha256) in self.driver.EXPECTED_DEPENDENCIES.items()
        ]
        self.assertEqual(len(destinations), len(set(destinations)))
        self.assertNotIn("driver.py", destinations)
        self.assertNotIn("protocol.json", destinations)
        helper = self.driver.load_module(
            "receipt_test_helper", self.driver.EXPECTED_DEPENDENCIES["helper"][0]
        )
        with TemporaryDirectory() as temporary:
            receipt_root = Path(temporary)
            protocol_path = receipt_root / "requested-protocol.json"
            builds_path = receipt_root / "requested-builds.json"
            protocol_path.write_text(json.dumps(registered))
            builds_path.write_text("{}\n")
            receipts = self.driver.freeze_execution_receipts(
                {"helper": helper, "protocol": registered},
                protocol_path,
                builds_path,
                receipt_root / "execution",
            )
            receipt_paths = [receipt["path"] for receipt in receipts.values()]
            self.assertEqual(len(receipt_paths), len(set(receipt_paths)))
            self.assertTrue(all(Path(path).is_file() for path in receipt_paths))

    def test_protocol_refuses_changed_sample_or_metadata(self):
        protocol = protocol_fixture(self.driver)
        protocol["samples"][0]["taxon"] = "changed"
        with self.assertRaisesRegex(ValueError, "Sample differs"):
            self.driver.validate_protocol(protocol)
        protocol = protocol_fixture(self.driver)
        protocol["target_metadata"]["insecta"][0]["scope"] = "changed"
        with self.assertRaisesRegex(ValueError, "Target metadata"):
            self.driver.validate_protocol(protocol)

    def test_master_schedule_rotation_and_version_order_are_exact(self):
        protocol = protocol_fixture(self.driver)
        schedule = self.driver.schedule_for_protocol(protocol)
        self.assertEqual(len(schedule), 72)
        blocks = [schedule[offset : offset + 2] for offset in range(0, len(schedule), 2)]
        first_k_by_sample = [blocks[index * 4][0]["k"] for index in range(3)]
        self.assertEqual(first_k_by_sample, [19, 23, 27])
        pair_two_offset = 3 * 4
        self.assertEqual(blocks[pair_two_offset][0]["k"], 23)
        for block in blocks:
            self.assertEqual(block[0]["k"], block[1]["k"])
            self.assertEqual({entry["version"] for entry in block}, {"baseline", "candidate"})
            invocation = block[0]
            allowed_index = protocol["allowed_k"].index(invocation["k"])
            baseline_first = (
                invocation["pair_index"] + invocation["cell_index"] + allowed_index
            ) % 2 == 1
            self.assertEqual(block[0]["version"] == "baseline", baseline_first)

    def test_fixed_k_schedule_never_mixes_k(self):
        protocol = protocol_fixture(self.driver)
        schedule = self.driver.fixed_k_schedule(
            self.driver.schedule_for_protocol(protocol), 27
        )
        self.assertEqual(len(schedule), 18)
        self.assertEqual({invocation["k"] for invocation in schedule}, {27})
        self.assertTrue(all("_k27_" in invocation["invocation_id"] for invocation in schedule))

    def test_dynamic_adapter_requires_actual_stats_k(self):
        helper = FakeHelper()
        runner = FakeRunner()
        legacy = object()
        self.driver.install_dynamic_adapters(helper, legacy, 27)
        stats = {
            "kmer_length": 27,
            "input_source": {
                "kind": "local_files",
                "inputs": ["/reads.fastq"],
                "paired": False,
                "max_reads": 100,
            },
        }
        invocation = {
            "version": "candidate",
            "depth": 100,
            "sample_prefix": "insecta_SRR1_100",
        }
        result = helper.validate_shared_stats(
            stats,
            invocation,
            {},
            {"path": "/reads.fastq"},
            ["sharkmer", "-k", "27"],
            "3.2.0-dev",
        )
        self.assertEqual(result, ["validated"])
        helper.validate_current_outputs(
            Path("/tmp"), invocation, {"stats_genes": set()}, stats, [], [], runner
        )
        self.assertEqual(runner.observed_k, 27)
        invalid = dict(stats)
        invalid["kmer_length"] = 19
        with self.assertRaisesRegex(ValueError, "selected sweep k"):
            helper.validate_shared_stats(
                invalid,
                invocation,
                {},
                {"path": "/reads.fastq"},
                [],
                "3.2.0-dev",
            )

    def test_result_k_checks_command_signature_and_raw_stats(self):
        invocation = {"invocation_id": "cell_k31_pair1_candidate", "k": 31}
        with TemporaryDirectory() as temporary:
            stats_path = Path(temporary) / "sample.stats.yaml"
            stats_path.write_text("kmer_length: 31\n")
            result = {
                "signature": {"settings": {"k": 31}, "invocation": invocation},
                "binary_command": ["sharkmer", "-k", "31"],
                "timing_status": "complete",
                "stats_path": str(stats_path),
            }
            self.driver.ensure_result_k(result, invocation, 31)
            result["binary_command"][-1] = "19"
            with self.assertRaisesRegex(ValueError, "command has wrong k"):
                self.driver.ensure_result_k(result, invocation, 31)
            result["binary_command"][-1] = "31"
            stats_path.write_text("kmer_length: 27\n")
            with self.assertRaisesRegex(ValueError, "stats has wrong k"):
                self.driver.ensure_result_k(result, invocation, 31)

    def test_classification_tool_receipt_binds_path_hash_and_version(self):
        with TemporaryDirectory() as temporary:
            tool = Path(temporary) / "blastn"
            tool.write_text("#!/bin/sh\necho blastn-test 1.0\n")
            tool.chmod(0o755)
            original_path = self.driver.os.environ.get("PATH", "")
            self.driver.os.environ["PATH"] = f"{temporary}:{original_path}"
            try:
                receipt = self.driver.tool_receipt("blastn")
            finally:
                self.driver.os.environ["PATH"] = original_path
            self.assertEqual(receipt["path"], str(tool.resolve()))
            self.assertEqual(receipt["sha256"], self.driver.sha256_file(tool))
            self.assertIn("blastn-test 1.0", receipt["version_stdout"])

    def test_measurement_index_receipt_detects_mutation(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text("{}\n")
            receipt = {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": self.driver.sha256_file(path),
            }
            self.assertTrue(self.driver.verify_receipt(receipt))
            path.write_text('{"changed": true}\n')
            self.assertFalse(self.driver.verify_receipt(receipt))

    def test_master_guard_requires_every_preregistered_k_before_classification(self):
        protocol = protocol_fixture(self.driver)
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt_paths = {}
            for label in ("driver", "protocol", "builds"):
                path = root / f"{label}.json"
                path.write_bytes(DRIVER_PATH.read_bytes() if label == "driver" else b"{}\n")
                receipt_paths[label] = {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": self.driver.sha256_file(path),
                }
            children = {}
            for selected_k in protocol["allowed_k"]:
                path = root / f"k{selected_k}" / "measurement-index.json"
                path.parent.mkdir()
                path.write_text("{}\n")
                children[f"k{selected_k}"] = {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": self.driver.sha256_file(path),
                }
            master_path = root / "master-schedule.json"
            master_path.write_text(json.dumps({
                "status": "complete",
                "schedule": self.driver.schedule_for_protocol(protocol),
                "receipts": receipt_paths,
                "child_measurement_indexes": children,
            }))
            self.driver.validate_master_measurement(protocol, master_path)
            del children["k31"]
            master_path.write_text(json.dumps({
                "status": "complete",
                "schedule": self.driver.schedule_for_protocol(protocol),
                "receipts": receipt_paths,
                "child_measurement_indexes": children,
            }))
            with self.assertRaisesRegex(ValueError, "child measurement set"):
                self.driver.validate_master_measurement(protocol, master_path)


if __name__ == "__main__":
    unittest.main()

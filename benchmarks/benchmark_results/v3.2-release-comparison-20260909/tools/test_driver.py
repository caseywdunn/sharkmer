import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path

import driver


def fastq_record(name, sequence):
    quality = "I" * len(sequence)
    return f"@{name}\n{sequence}\n+\n{quality}\n".encode()


class DriverTests(unittest.TestCase):
    def test_scan_fastq_records_prefix_bytes_bases_and_sha(self):
        content = fastq_record("one", "ACGT") + fastq_record("two", "AAAAAA") + fastq_record("three", "CC")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_bytes(content)
            observed = driver.scan_fastq(path, {1: {}, 2: {}, 5: {}})
        first = fastq_record("one", "ACGT")
        first_two = first + fastq_record("two", "AAAAAA")
        self.assertEqual(observed["total_records"], 3)
        self.assertEqual(observed["total_bases"], 12)
        self.assertEqual(observed["prefixes"][1], {"records": 1, "bases": 4, "size_bytes": len(first), "sha256": hashlib.sha256(first).hexdigest()})
        self.assertEqual(observed["prefixes"][2], {"records": 2, "bases": 10, "size_bytes": len(first_two), "sha256": hashlib.sha256(first_two).hexdigest()})
        self.assertEqual(observed["prefixes"][5], {"records": 3, "bases": 12, "size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})

    def test_scan_fastq_rejects_gzip_and_malformed_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            gzip_path = Path(temporary) / "reads.fastq"
            gzip_path.write_bytes(b"\x1f\x8bnot-really-gzip")
            with self.assertRaisesRegex(ValueError, "uncompressed"):
                driver.scan_fastq(gzip_path, {1: {}})
            malformed_path = Path(temporary) / "bad.fastq"
            malformed_path.write_bytes(b"@one\nACGT\n+\nIII\n")
            with self.assertRaisesRegex(ValueError, "length mismatch"):
                driver.scan_fastq(malformed_path, {1: {}})

    def test_schedule_has_three_primary_pairs_and_alternates_order(self):
        protocol = {
            "schedule": {"primary_depth": 1_000_000, "primary_pairs": 3, "deeper_pairs": 1},
            "order_seed": "release-comparison",
            "samples": [
                {"panel": "cnidaria", "input": "SRR1", "taxon": "Taxon", "depths": [1_000_000, 2_000_000]}
            ],
        }
        schedule = driver.build_schedule(protocol)
        primary = [entry for entry in schedule if entry["depth"] == 1_000_000]
        deeper = [entry for entry in schedule if entry["depth"] == 2_000_000]
        self.assertEqual(len(primary), 6)
        self.assertEqual(len(deeper), 2)
        orders = []
        for pair_index in (1, 2, 3):
            orders.append([entry["version"] for entry in primary if entry["pair_index"] == pair_index])
        self.assertEqual(orders[0], orders[2])
        self.assertEqual(orders[1], list(reversed(orders[0])))

    def test_schedule_places_all_primary_cells_before_deeper_cells(self):
        protocol = {
            "schedule": {"primary_depth": 1_000_000, "primary_pairs": 3, "deeper_pairs": 1},
            "order_seed": "release-comparison",
            "samples": [
                {"panel": "cnidaria", "input": "SRR1", "taxon": "One", "depths": [1_000_000, 2_000_000, 4_000_000]},
                {"panel": "insecta", "input": "SRR2", "taxon": "Two", "depths": [1_000_000, 2_000_000, 4_000_000]},
            ],
        }
        schedule = driver.build_schedule(protocol)
        cell_depths = []
        for entry in schedule:
            if not cell_depths or cell_depths[-1] != (entry["cell_index"], entry["depth"]):
                cell_depths.append((entry["cell_index"], entry["depth"]))
        self.assertEqual([depth for _, depth in cell_depths], [1_000_000, 1_000_000, 2_000_000, 2_000_000, 4_000_000, 4_000_000])

    def test_freeze_receipt_refuses_changed_resume_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "source.json"
            frozen = directory / "receipts" / "source.json"
            source.write_text('{"value": 1}\n')
            first = driver.freeze_receipt(source, frozen)
            second = driver.freeze_receipt(source, frozen)
            self.assertEqual(first, second)
            source.write_text('{"value": 2}\n')
            with self.assertRaisesRegex(ValueError, "Frozen receipt differs"):
                driver.freeze_receipt(source, frozen)

    def test_parse_fasta_reads_every_product_and_current_header(self):
        content = ">sample product=0 kmer_count_median=7\nACGT\n>sample product=1 kmer_count_median=11\nTTAA\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "products.fasta"
            path.write_text(content)
            products = driver.parse_fasta_file(path)
        self.assertEqual([product["product_index"] for product in products], [0, 1])
        self.assertEqual([product["kmer_count_median"] for product in products], [7, 11])
        self.assertEqual([product["sha256"] for product in products], [hashlib.sha256(b"ACGT").hexdigest(), hashlib.sha256(b"TTAA").hexdigest()])

    def test_parse_fasta_rejects_obsolete_or_duplicate_headers(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "products.fasta"
            path.write_text(">sample median 7\nACGT\n")
            with self.assertRaisesRegex(ValueError, "product index"):
                driver.parse_fasta_file(path)
            path.write_text(">sample product=0 kmer_count_median=7\nACGT\n>sample product=0 kmer_count_median=8\nTTAA\n")
            with self.assertRaisesRegex(ValueError, "duplicate product"):
                driver.parse_fasta_file(path)

    def test_legacy_adapter_derives_outputs_and_rejects_stale_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            invocation = {"sample_prefix": "panel_SRR1_1000000"}
            panel = {"output_prefix": "panel"}
            stats = {"sharkmer_version": "3.1.0"}
            pcr_results = [
                {"gene_name": "panel_gene", "status": "success", "n_products": 1, "product_lengths": [4], "output_file": None},
                {"gene_name": "panel_absent", "status": "fail", "n_products": 0, "failure_reason": "no path"},
            ]
            stats_path = output_dir / "panel_SRR1_1000000.stats.yaml"
            stats_path.write_text(json.dumps(stats))
            fasta_path = output_dir / "panel_SRR1_1000000_panel_gene.fasta"
            fasta_path.write_text(">sample product=0 kmer_count_median=7\nACGT\n")
            genes, evidence = driver.validate_legacy_outputs(output_dir, invocation, panel, stats, pcr_results)
            self.assertTrue(genes[0]["recovered"])
            self.assertFalse(evidence["manifest_available"])
            (output_dir / "stale.fasta").write_text(">stale\nACGT\n")
            with self.assertRaisesRegex(ValueError, "file set mismatch"):
                driver.validate_legacy_outputs(output_dir, invocation, panel, stats, pcr_results)

    def test_legacy_adapter_rejects_current_transaction_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            stats = {"sharkmer_version": "3.1.0", "run_id": "unexpected"}
            invocation = {"sample_prefix": "panel_SRR1_1000000"}
            with self.assertRaisesRegex(ValueError, "current-only"):
                driver.validate_legacy_outputs(output_dir, invocation, {"output_prefix": "panel"}, stats, [])

    def test_shared_stats_accepts_attested_source_exhaustion(self):
        invocation = {"sample_prefix": "panel_SRR1_2000000", "depth": 2_000_000}
        panel = {"stats_genes": {"panel_gene"}}
        input_record = {"prefixes": {2_000_000: {"records": 17, "bases": 68}}}
        command = ["/binary", "--max-reads", "2000000", "reads.fastq"]
        stats = {
            "sample": invocation["sample_prefix"],
            "sharkmer_version": "3.1.0",
            "kmer_length": 19,
            "chunks": 0,
            "command": " ".join(command),
            "n_reads_read": 17,
            "n_bases_read": 68,
            "n_subreads_ingested": 17,
            "n_bases_ingested": 68,
            "n_kmers": 34,
            "peak_memory_bytes": 100,
            "pcr_results": [
                {"gene_name": "panel_gene", "status": "success", "n_products": 1, "product_lengths": [40]}
            ],
        }
        results = driver.validate_shared_stats(stats, invocation, panel, input_record, command, "3.1.0")
        self.assertEqual(len(results), 1)
        stats["n_reads_read"] = 2_000_000
        with self.assertRaisesRegex(ValueError, "record count"):
            driver.validate_shared_stats(stats, invocation, panel, input_record, command, "3.1.0")

    def test_completed_resume_rehashes_logs_stats_and_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_root = Path(temporary)
            invocation_id = "000_panel_SRR1_1000000_pair1_baseline"
            attempt_id = "attempt"
            attempt_dir = output_root / "attempts" / invocation_id / attempt_id
            output_dir = attempt_dir / "output"
            output_dir.mkdir(parents=True)
            stdout_path = attempt_dir / "stdout.log"
            stderr_path = attempt_dir / "stderr.log"
            time_path = attempt_dir / "gnu-time.txt"
            stats_path = output_dir / "sample.stats.yaml"
            stdout_path.write_text("stdout\n")
            stderr_path.write_text("stderr\n")
            time_path.write_text("time\n")
            stats_path.write_text("stats\n")
            signature = {"frozen": True}
            result = {
                "status": "complete",
                "timing_status": "complete",
                "signature": signature,
                "attempt_id": attempt_id,
                "attempt_dir": str(attempt_dir),
                "logs": {
                    "stdout": driver.file_receipt(stdout_path),
                    "stderr": driver.file_receipt(stderr_path),
                    "gnu_time": driver.file_receipt(time_path),
                },
                "stats_path": str(stats_path),
                "stats_sha256": driver.sha256_file(stats_path),
                "raw_output_files": driver.raw_output_files(output_dir),
                "unvalidated_output_inventory": driver.inventory_untrusted_directory(output_dir),
            }
            self.assertEqual(
                driver.validate_preserved_result(result, signature, output_root, invocation_id), result
            )
            stdout_path.write_text("changed\n")
            with self.assertRaisesRegex(ValueError, "logs changed"):
                driver.validate_preserved_result(result, signature, output_root, invocation_id)

    def test_failed_resume_remains_primary_outcome(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_root = Path(temporary)
            invocation_id = "000_panel_SRR1_1000000_pair1_baseline"
            attempt_id = "attempt"
            attempt_dir = output_root / "attempts" / invocation_id / attempt_id
            output_dir = attempt_dir / "output"
            output_dir.mkdir(parents=True)
            stdout_path = attempt_dir / "stdout.log"
            stderr_path = attempt_dir / "stderr.log"
            stdout_path.write_text("stdout\n")
            stderr_path.write_text("failure\n")
            signature = {"frozen": True}
            result = {
                "status": "failed",
                "timing_status": "failed",
                "failure": "nonzero_exit",
                "signature": signature,
                "attempt_id": attempt_id,
                "attempt_dir": str(attempt_dir),
                "logs": {
                    "stdout": driver.file_receipt(stdout_path),
                    "stderr": driver.file_receipt(stderr_path),
                    "gnu_time": None,
                },
                "unvalidated_output_inventory": [],
            }
            self.assertEqual(
                driver.validate_preserved_result(result, signature, output_root, invocation_id), result
            )

    def test_classification_is_a_separate_post_timing_transition(self):
        class FakeBlast:
            @staticmethod
            def extract_references(panel_data):
                return []

            @staticmethod
            def blast_all_products(run_results, db_path, taxon, skip_blast, reference_genes):
                self.assertIsNone(db_path)
                self.assertEqual(taxon, "Taxon")
                self.assertFalse(skip_blast)
                self.assertEqual(reference_genes, set())
                run_results[0]["genes"][0]["products"][0]["reference_match"] = {"status": "no_reference"}

        with tempfile.TemporaryDirectory() as temporary:
            output_root = Path(temporary)
            invocation = {"invocation_id": "one", "taxon": "Taxon"}
            result_path = output_root / "results" / "one.json"
            timed_result = {
                "status": "timed_complete",
                "timing_status": "complete",
                "genes": [
                    {
                        "gene": "gene",
                        "products": [{"product_index": 0, "sequence": "ACGT", "sha256": hashlib.sha256(b"ACGT").hexdigest()}],
                    }
                ],
            }
            panel = {"data": {}}
            finalized = driver.finalize_classification(
                invocation, timed_result, panel, output_root, FakeBlast(), None
            )
            self.assertEqual(finalized["status"], "complete")
            self.assertEqual(finalized["classification_status"], "complete")
            self.assertNotIn("sequence", finalized["genes"][0]["products"][0])
            self.assertEqual(driver.load_json(result_path)["status"], "complete")

    def test_run_timed_preserves_nonzero_logs_and_timeout(self):
        allowed_cpu = min(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 0
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stdout_path = directory / "stdout.log"
            stderr_path = directory / "stderr.log"
            time_path = directory / "time.txt"
            failed = driver.run_timed(
                ["/bin/sh", "-c", "echo preserved-out; echo preserved-err >&2; exit 7"],
                stdout_path,
                stderr_path,
                time_path,
                5,
                1024 * 1024 * 1024,
                [allowed_cpu],
                {**os.environ, "LC_ALL": "C"},
            )
            self.assertEqual(failed["returncode"], 7)
            self.assertIn("preserved-out", stdout_path.read_text())
            self.assertIn("preserved-err", stderr_path.read_text())
            timed_out = driver.run_timed(
                ["/bin/sh", "-c", "sleep 5 & wait"],
                directory / "timeout-stdout.log",
                directory / "timeout-stderr.log",
                directory / "timeout-time.txt",
                0.05,
                1024 * 1024 * 1024,
                [allowed_cpu],
                {**os.environ, "LC_ALL": "C"},
            )
            self.assertTrue(timed_out["timed_out"])
            self.assertIsNotNone(timed_out["returncode"])

    def test_summary_keeps_null_legacy_metrics_and_compares_pairs(self):
        schedule = [
            {"invocation_id": "baseline", "cell": "panel/SRR1/1000000", "pair_index": 1, "version": "baseline"},
            {"invocation_id": "candidate", "cell": "panel/SRR1/1000000", "pair_index": 1, "version": "candidate"},
        ]
        metrics = {
            "n_reads_read": 10,
            "n_bases_read": 40,
            "n_subreads_ingested": 10,
            "n_bases_ingested": 40,
            "n_kmers": 22,
            "end_to_end_input_mbp_s": 0.00004,
            "end_to_end_ingested_kmers_s": 22,
        }
        product = {"product_index": 0, "length": 4, "sha256": hashlib.sha256(b"ACGT").hexdigest(), "reference_match": {"status": "confirmed_product"}}
        results = {
            "baseline": {"status": "complete", "metrics": {**metrics, "stage_timings": None}, "genes": [{"gene": "gene", "products": [product]}], "execution": {"wall_time_s": 2.0, "gnu_time": {"peak_rss_bytes": 100}}},
            "candidate": {"status": "complete", "metrics": {**metrics, "stage_timings": {"pipeline_total_s": 1.0}}, "genes": [{"gene": "gene", "products": [product]}], "execution": {"wall_time_s": 1.0, "gnu_time": {"peak_rss_bytes": 90}}},
        }
        summary = driver.summarize_results(schedule, results)
        comparison = summary["pair_comparisons"][0]
        self.assertTrue(comparison["counts_identical"])
        self.assertTrue(comparison["products_identical"])
        self.assertTrue(comparison["classifications_identical"])
        self.assertEqual(comparison["candidate_speedup_ratio"], 2.0)
        self.assertEqual(summary["failures"], [])

    def test_normalized_metrics_derives_end_to_end_rates(self):
        stats = {
            "n_reads_read": 10,
            "n_bases_read": 4_000_000,
            "n_subreads_ingested": 10,
            "n_bases_ingested": 4_000_000,
            "n_kmers": 6_000_000,
            "peak_memory_bytes": 100,
        }
        metrics = driver.normalized_metrics(stats, 2.0)
        self.assertEqual(metrics["end_to_end_input_mbp_s"], 2.0)
        self.assertEqual(metrics["end_to_end_ingested_kmers_s"], 3_000_000)

    def test_parse_gnu_time_uses_elapsed_and_rss(self):
        content = "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50\n\tMaximum resident set size (kbytes): 2048\n\tExit status: 0\n"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "time.txt"
            path.write_text(content)
            parsed = driver.parse_gnu_time(path)
        self.assertEqual(parsed["wall_time_s"], 62.5)
        self.assertEqual(parsed["peak_rss_bytes"], 2 * 1024 * 1024)
        self.assertEqual(parsed["exit_status"], 0)


if __name__ == "__main__":
    unittest.main()

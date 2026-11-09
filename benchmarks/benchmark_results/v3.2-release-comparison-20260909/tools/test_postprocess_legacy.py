import copy
import json
import tempfile
import unittest
from pathlib import Path

import driver
import postprocess_legacy
import yaml


def header(sample, gene, product, sequence, median, score=None):
    score = median if score is None else score
    return (
        f">{sample}_{gene}_{product} sample={sample} gene={gene} product={product} "
        f"length={len(sequence)} kmer_count_mean=67.25 kmer_count_median={median} "
        f"kmer_count_min=10 kmer_count_max=100 score={score}\n{sequence}\n"
    )


class LegacyPostprocessTests(unittest.TestCase):
    def test_parser_accepts_every_integer_and_half_median_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "products.fasta"
            path.write_text(
                header("sample", "panel_gene", 0, "ACGT", "67.5")
                + header("sample", "panel_gene", 1, "TTAA", "68")
            )
            products = postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")
        self.assertEqual([product["product_index"] for product in products], [0, 1])
        self.assertEqual([product["kmer_count_median"] for product in products], [67.5, 68])
        self.assertEqual([product["kmer_count_median_raw"] for product in products], ["67.5", "68"])

    def test_parser_preserves_independent_composite_score(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "products.fasta"
            content = header("sample", "panel_gene", 0, "ACGT", "6", score="0.11").replace(
                "kmer_count_min=10", "kmer_count_min=2"
            )
            path.write_text(content)
            product = postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")[0]
        self.assertEqual(product["kmer_count_median"], 6)
        self.assertEqual(product["score"], 0.11)
        self.assertEqual(product["score_raw"], "0.11")

    def test_parser_rejects_values_outside_legacy_median_domain(self):
        for invalid in ("67.25", "67.0", "-1", "01", "nan"):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "products.fasta"
                path.write_text(header("sample", "panel_gene", 0, "ACGT", invalid))
                with self.assertRaisesRegex(ValueError, "integer or .5"):
                    postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")

    def test_parser_rejects_header_identity_length_and_index_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "products.fasta"
            path.write_text(header("wrong", "panel_gene", 0, "ACGT", "67.5"))
            with self.assertRaisesRegex(ValueError, "sample or gene"):
                postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")
            path.write_text(header("sample", "panel_gene", 1, "ACGT", "67.5"))
            with self.assertRaisesRegex(ValueError, "contiguous"):
                postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")
            content = header("sample", "panel_gene", 0, "ACGT", "67.5").replace("length=4", "length=5")
            path.write_text(content)
            with self.assertRaisesRegex(ValueError, "header length"):
                postprocess_legacy.parse_legacy_fasta(path, "sample", "panel_gene")

    def test_only_exact_baseline_adapter_failure_is_correctable(self):
        result = {
            "status": "failed",
            "timing_status": "failed",
            "failure": "invalid_output",
            "error": postprocess_legacy.EXPECTED_LEGACY_ERROR + "/output/file.fasta",
            "signature": {"invocation": {"version": "baseline"}},
            "execution": {
                "returncode": 0,
                "timed_out": False,
                "orphaned_process_group_cleaned": False,
                "measurement_complete": True,
                "measurement_matches_exit": True,
                "input_changed_during_invocation": False,
                "binary_changed_during_invocation": False,
                "panel_changed_during_invocation": False,
            },
        }
        self.assertTrue(postprocess_legacy.should_correct_legacy_result(result))
        candidate = copy.deepcopy(result)
        candidate["signature"]["invocation"]["version"] = "candidate"
        self.assertFalse(postprocess_legacy.should_correct_legacy_result(candidate))
        other_failure = copy.deepcopy(result)
        other_failure["error"] = "ValueError: unrelated invalid output"
        self.assertFalse(postprocess_legacy.should_correct_legacy_result(other_failure))
        nonzero = copy.deepcopy(result)
        nonzero["execution"]["returncode"] = 1
        self.assertFalse(postprocess_legacy.should_correct_legacy_result(nonzero))

    def test_normalization_rederives_stats_files_lengths_and_lineage(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            attempt_dir = directory / "attempt"
            output_dir = attempt_dir / "output"
            output_dir.mkdir(parents=True)
            sample = "panel_SRR1_100"
            gene = "panel_gene"
            command = ["/binary", "--max-reads", "100", "reads.fastq"]
            stats = {
                "sample": sample,
                "sharkmer_version": "3.1.0",
                "kmer_length": 19,
                "chunks": 0,
                "command": " ".join(command),
                "n_reads_read": 100,
                "n_bases_read": 400,
                "n_subreads_ingested": 100,
                "n_bases_ingested": 400,
                "n_kmers": 200,
                "peak_memory_bytes": 1024,
                "pcr_results": [
                    {
                        "gene_name": gene,
                        "status": "success",
                        "n_products": 1,
                        "product_lengths": [4],
                        "output_file": None,
                    }
                ],
            }
            stats_path = output_dir / f"{sample}.stats.yaml"
            stats_path.write_text(yaml.safe_dump(stats, sort_keys=False))
            fasta_path = output_dir / f"{sample}_{gene}.fasta"
            fasta_path.write_text(header(sample, gene, 0, "ACGT", "67.5"))
            result_path = directory / "original.json"
            original = {
                "status": "failed",
                "timing_status": "failed",
                "failure": "invalid_output",
                "error": postprocess_legacy.EXPECTED_LEGACY_ERROR + str(fasta_path),
                "attempt_dir": str(attempt_dir),
                "binary_command": command,
                "execution": {"wall_time_s": 2.0},
                "signature": {"binary_sha256": "a" * 64},
            }
            result_path.write_text(json.dumps(original))
            invocation = {"sample_prefix": sample, "depth": 100}
            panel = {"output_prefix": "panel", "stats_genes": {gene}}
            input_record = {"prefixes": {100: {"records": 100, "bases": 400}}}
            build = {"binary_sha256": "a" * 64}
            normalized = postprocess_legacy.normalize_legacy_adapter_failure(
                result_path,
                original,
                invocation,
                panel,
                input_record,
                build,
                driver,
            )
        self.assertEqual(normalized["status"], "normalized_legacy_pending_classification")
        self.assertEqual(normalized["timing_status"], "complete")
        self.assertEqual(normalized["genes"][0]["products"][0]["kmer_count_median"], 67.5)
        self.assertEqual(normalized["genes"][0]["products"][0]["kmer_count_median_raw"], "67.5")
        self.assertEqual(normalized["lineage"]["original_failure"], "invalid_output")
        self.assertEqual(normalized["metrics"]["end_to_end_ingested_kmers_s"], 100)

    def test_normalization_rejects_unexpected_extra_legacy_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output_dir = directory / "attempt" / "output"
            output_dir.mkdir(parents=True)
            sample = "panel_SRR1_100"
            gene = "panel_gene"
            command = ["/binary"]
            stats = {
                "sample": sample,
                "sharkmer_version": "3.1.0",
                "kmer_length": 19,
                "chunks": 0,
                "command": "/binary",
                "n_reads_read": 100,
                "n_bases_read": 400,
                "n_subreads_ingested": 100,
                "n_bases_ingested": 400,
                "n_kmers": 200,
                "peak_memory_bytes": 1024,
                "pcr_results": [
                    {"gene_name": gene, "status": "success", "n_products": 1, "product_lengths": [4]}
                ],
            }
            (output_dir / f"{sample}.stats.yaml").write_text(yaml.safe_dump(stats))
            (output_dir / f"{sample}_{gene}.fasta").write_text(header(sample, gene, 0, "ACGT", "67.5"))
            (output_dir / "stale.fasta").write_text(header(sample, gene, 0, "ACGT", "67.5"))
            result = {
                "attempt_dir": str(directory / "attempt"),
                "binary_command": command,
                "execution": {"wall_time_s": 2.0},
                "signature": {"binary_sha256": "a" * 64},
            }
            with self.assertRaisesRegex(ValueError, "file set mismatch"):
                postprocess_legacy.validate_correctable_legacy_result(
                    result,
                    {"sample_prefix": sample, "depth": 100},
                    {"output_prefix": "panel", "stats_genes": {gene}},
                    {"prefixes": {100: {"records": 100, "bases": 400}}},
                    {"binary_sha256": "a" * 64},
                    driver,
                )

    def test_completed_legacy_result_is_rederived_from_raw_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output_dir = directory / "attempt" / "output"
            output_dir.mkdir(parents=True)
            sample = "panel_SRR1_100"
            gene = "panel_gene"
            command = ["/binary"]
            stats = {
                "sample": sample,
                "sharkmer_version": "3.1.0",
                "kmer_length": 19,
                "chunks": 0,
                "command": "/binary",
                "n_reads_read": 100,
                "n_bases_read": 400,
                "n_subreads_ingested": 100,
                "n_bases_ingested": 400,
                "n_kmers": 200,
                "peak_memory_bytes": 1024,
                "pcr_results": [
                    {"gene_name": gene, "status": "success", "n_products": 1, "product_lengths": [4]}
                ],
            }
            stats_path = output_dir / f"{sample}.stats.yaml"
            stats_path.write_text(yaml.safe_dump(stats))
            (output_dir / f"{sample}_{gene}.fasta").write_text(header(sample, gene, 0, "ACGT", "68"))
            invocation = {
                "invocation_id": "one",
                "version": "baseline",
                "sample_prefix": sample,
                "depth": 100,
            }
            panel = {"output_prefix": "panel", "stats_genes": {gene}}
            input_record = {"prefixes": {100: {"records": 100, "bases": 400}}}
            parsed_products = driver.parse_fasta_file(output_dir / f"{sample}_{gene}.fasta")
            parsed_products[0].pop("sequence")
            parsed_products[0]["reference_match"] = {"status": "no_reference"}
            result = {
                "status": "complete",
                "timing_status": "complete",
                "attempt_dir": str(directory / "attempt"),
                "binary_command": command,
                "execution": {"wall_time_s": 2.0},
                "completion": {
                    "completion_evidence": "exit_zero_validated_stats_fasta_exact_fresh_directory",
                    "manifest_available": False,
                },
                "stats_path": str(stats_path),
                "stats_sha256": driver.sha256_file(stats_path),
                "metrics": driver.normalized_metrics(stats, 2.0),
                "genes": [
                    {"gene": "gene", "recovered": True, "n_products": 1, "products": parsed_products}
                ],
            }
            postprocess_legacy.revalidate_completed_result(
                result,
                invocation,
                panel,
                input_record,
                {"binary_version": "sharkmer 3.1.0 (source)"},
                driver,
                None,
            )
            result["genes"][0]["products"][0]["length"] = 5
            with self.assertRaisesRegex(ValueError, "genes differ"):
                postprocess_legacy.revalidate_completed_result(
                    result,
                    invocation,
                    panel,
                    input_record,
                    {"binary_version": "sharkmer 3.1.0 (source)"},
                    driver,
                    None,
                )


if __name__ == "__main__":
    unittest.main()

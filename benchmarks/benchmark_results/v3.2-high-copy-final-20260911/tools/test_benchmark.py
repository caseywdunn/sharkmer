import importlib.util
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parent
FROZEN_ROOT = Path("/tmp/sharkmer-release-comparison")
BENCHMARK_PATH = ROOT / "benchmark.py"
HELPER_PATH = FROZEN_ROOT / "execution" / "receipts" / "driver.py"
LEGACY_PATH = FROZEN_ROOT / "postprocess_legacy.py"


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def legacy_header(sample, gene, median="67.5"):
    return (
        f">{sample}_{gene}_0 sample={sample} gene={gene} product=0 length=4 "
        f"kmer_count_mean=67.5 kmer_count_median={median} "
        "kmer_count_min=10 kmer_count_max=100 score=0.11\nACGT\n"
    )


def invocation(version):
    return {
        "version": version,
        "sample_prefix": "panel_SRR1_100",
        "depth": 100,
    }


def panel():
    return {"output_prefix": "panel", "stats_genes": {"panel_gene"}}


def pcr_results():
    return [
        {
            "gene_name": "panel_gene",
            "status": "success",
            "n_products": 1,
            "product_lengths": [4],
        }
    ]


def input_record():
    return {
        "path": "/frozen/reads.fastq",
        "prefixes": {100: {"records": 100, "bases": 400}},
    }


def shared_stats(version):
    command = ["/frozen/sharkmer", "--max-reads", "100", "/frozen/reads.fastq"]
    return {
        "sample": "panel_SRR1_100",
        "sharkmer_version": version,
        "kmer_length": 19,
        "chunks": 0,
        "command": " ".join(command),
        "n_reads_read": 100,
        "n_bases_read": 400,
        "n_subreads_ingested": 100,
        "n_bases_ingested": 400,
        "n_kmers": 200,
        "peak_memory_bytes": 1024,
        "pcr_results": pcr_results(),
    }, command


class BenchmarkAdapterTests(unittest.TestCase):
    def configured_helper(self):
        benchmark = load_module("high_copy_benchmark", BENCHMARK_PATH)
        helper = load_module("high_copy_frozen_helper", HELPER_PATH)
        legacy = load_module("high_copy_legacy_parser", LEGACY_PATH)
        benchmark.install_adapters(helper, legacy)
        return helper

    def write_legacy_output(self, directory, header):
        sample = invocation("baseline")["sample_prefix"]
        (directory / f"{sample}.stats.yaml").write_text("sample: panel_SRR1_100\n")
        fasta_path = directory / f"{sample}_panel_gene.fasta"
        fasta_path.write_text(header)
        return fasta_path

    def test_legacy_half_median_uses_corrected_parser_and_original_file_checks(self):
        helper = self.configured_helper()
        with TemporaryDirectory() as temporary:
            output_directory = Path(temporary)
            self.write_legacy_output(output_directory, legacy_header("panel_SRR1_100", "panel_gene"))
            genes, completion = helper.validate_legacy_outputs(
                output_directory,
                invocation("baseline"),
                panel(),
                {"sharkmer_version": "3.1.0"},
                pcr_results(),
            )
        self.assertEqual(genes[0]["products"][0]["kmer_count_median"], 67.5)
        self.assertEqual(genes[0]["products"][0]["kmer_count_median_raw"], "67.5")
        self.assertEqual(completion["manifest_available"], False)

    def test_legacy_adapter_refuses_wrong_header_identity_and_extra_file(self):
        helper = self.configured_helper()
        with TemporaryDirectory() as temporary:
            output_directory = Path(temporary)
            self.write_legacy_output(output_directory, legacy_header("wrong", "panel_gene"))
            with self.assertRaisesRegex(ValueError, "sample or gene"):
                helper.validate_legacy_outputs(
                    output_directory,
                    invocation("baseline"),
                    panel(),
                    {"sharkmer_version": "3.1.0"},
                    pcr_results(),
                )
            self.write_legacy_output(output_directory, legacy_header("panel_SRR1_100", "wrong_gene"))
            with self.assertRaisesRegex(ValueError, "sample or gene"):
                helper.validate_legacy_outputs(
                    output_directory,
                    invocation("baseline"),
                    panel(),
                    {"sharkmer_version": "3.1.0"},
                    pcr_results(),
                )
            self.write_legacy_output(output_directory, legacy_header("panel_SRR1_100", "panel_gene"))
            (output_directory / "foreign.fasta").write_text(legacy_header("panel_SRR1_100", "panel_gene"))
            with self.assertRaisesRegex(ValueError, "file set mismatch"):
                helper.validate_legacy_outputs(
                    output_directory,
                    invocation("baseline"),
                    panel(),
                    {"sharkmer_version": "3.1.0"},
                    pcr_results(),
                )

    def test_legacy_adapter_restores_current_parser_after_error(self):
        helper = self.configured_helper()
        original_parser = helper.parse_fasta_file
        with TemporaryDirectory() as temporary:
            output_directory = Path(temporary)
            fasta_path = self.write_legacy_output(output_directory, legacy_header("wrong", "panel_gene"))
            with self.assertRaisesRegex(ValueError, "sample or gene"):
                helper.validate_legacy_outputs(
                    output_directory,
                    invocation("baseline"),
                    panel(),
                    {"sharkmer_version": "3.1.0"},
                    pcr_results(),
                )
            self.assertIs(helper.parse_fasta_file, original_parser)
            with self.assertRaisesRegex(ValueError, "current kmer median"):
                helper.parse_fasta_file(fasta_path)

    def test_current_input_source_is_exact_and_baseline_is_exempt(self):
        helper = self.configured_helper()
        baseline_stats, command = shared_stats("3.1.0")
        baseline_stats["input_source"] = {"kind": "remote", "inputs": [], "paired": True, "max_reads": 1}
        helper.validate_shared_stats(
            baseline_stats,
            invocation("baseline"),
            panel(),
            input_record(),
            command,
            "3.1.0",
        )

        candidate_stats, command = shared_stats("3.2.0-dev")
        candidate_stats["input_source"] = {
            "kind": "local_files",
            "inputs": ["/frozen/reads.fastq"],
            "paired": False,
            "max_reads": 100,
        }
        helper.validate_shared_stats(
            candidate_stats,
            invocation("candidate"),
            panel(),
            input_record(),
            command,
            "3.2.0-dev",
        )

        for field, invalid_value in (
            ("kind", "remote"),
            ("inputs", ["/other/reads.fastq"]),
            ("paired", True),
            ("max_reads", 99),
        ):
            with self.subTest(field=field):
                invalid_stats, command = shared_stats("3.2.0-dev")
                invalid_stats["input_source"] = {
                    "kind": "local_files",
                    "inputs": ["/frozen/reads.fastq"],
                    "paired": False,
                    "max_reads": 100,
                }
                invalid_stats["input_source"][field] = invalid_value
                with self.assertRaisesRegex(ValueError, "input source differs"):
                    helper.validate_shared_stats(
                        invalid_stats,
                        invocation("candidate"),
                        panel(),
                        input_record(),
                        command,
                        "3.2.0-dev",
                    )


if __name__ == "__main__":
    unittest.main()

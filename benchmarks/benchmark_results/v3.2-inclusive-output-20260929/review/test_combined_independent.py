"""Independent consumer controls for inclusive delivery; synthetic data only.

Fixture construction is borrowed explicitly from the implementation author's
tests/test_uncertain_output.py. The scenarios and assertions below were written
independently; borrowed fixtures are not independent biological evidence.
"""
import copy
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location("inclusive_author_fixtures", ROOT / "tests/test_uncertain_output.py")
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)
runner = fixtures.runner
uncertain = fixtures.uncertain


def two_targets(directory, second_standard=True):
    first = fixtures.combined_transaction(directory, gene="A")
    second = fixtures.combined_transaction(directory, gene="A.uncertain", standard=second_standard)
    first["pcr_results"].extend(second["pcr_results"])
    fixtures.publish(directory, first, [row["output_file"] for row in first["pcr_results"]])
    return first


class IndependentCombinedControls(unittest.TestCase):
    def test_suffix_named_target_has_disjoint_record_ids_and_standard_selection(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = two_targets(directory)
            runner._validate_stats_manifest(stats, "sample", 19, {"A", "A.uncertain"}, "end-to-end")
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual([(row["gene"], row["sequences"]) for row in parsed],
                             [("A", ["ACGT"]), ("A.uncertain", ["ACGT"])])
            ids = [line.split()[0][1:] for row in stats["pcr_results"]
                   for line in (directory / row["output_file"]).read_text().splitlines()
                   if line.startswith(">")]
            self.assertEqual(len(ids), 4)
            self.assertEqual(len(set(ids)), 4)
            self.assertIn("sample_A.uncertain_0", ids)
            self.assertIn("sample_A_0.uncertain", ids)

    def test_requested_standard_file_does_not_bypass_other_uncertain_only_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = two_targets(directory, second_standard=False)
            path = directory / "sample_A.uncertain.fasta"
            path.write_text(path.read_text().replace("assembly_status=uncertain", "assembly_status=standard", 1))
            fixtures.publish(directory, stats, [row["output_file"] for row in stats["pcr_results"]])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory, ["sample_A.fasta"])

    def test_class_partition_cannot_hide_duplicate_id(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = fixtures.combined_transaction(directory)
            path = directory / "sample_gene.fasta"
            path.write_text(path.read_text().replace(">sample_gene_0 assembly_status=standard",
                                                    ">sample_gene_0.uncertain assembly_status=standard", 1))
            fixtures.publish(directory, stats, [path.name])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_outer_schema2_and_inner_assessment_schema1_are_distinct(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = fixtures.combined_transaction(directory)
            self.assertEqual(len(runner.parse_fasta_products("sample", directory)), 1)
            stats["pcr_results"][0]["uncertain_products"]["schema_version"] = 2
            fixtures.publish(directory, stats, ["sample_gene.fasta"])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_relabeling_stats_as_legacy_cannot_admit_manifest2(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = fixtures.combined_transaction(directory, disabled=True)
            for key in ("product_output_schema_version", "product_output_layout"):
                del stats[key]
            for key in ("n_standard_products", "n_uncertain_products"):
                del stats["pcr_results"][0][key]
            manifest = fixtures.publish(directory, stats, ["sample_gene.fasta"])
            manifest["schema_version"] = 2
            (directory / "sample.manifest.yaml").write_text(fixtures.yaml.safe_dump(manifest))
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_unowned_legacy_sidecar_is_ignored_but_cannot_join_current_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = fixtures.combined_transaction(directory, disabled=True)
            sidecar = directory / "sample_gene.uncertain.fasta"
            sidecar.write_text(">foreign assembly_status=uncertain\nAAAA\n")
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual(parsed[0]["sequences"], ["ACGT"])
            fixtures.publish(directory, stats, ["sample_gene.fasta", sidecar.name])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_quoted_duplicate_and_empty_status_are_not_classes(self):
        for header in ("id assembly_status=standard 'assembly_status=uncertain'",
                       "id assembly_status=standard assembly_status=",
                       "id assembly_status=''", "id assembly_status=STANDARD"):
            with self.subTest(header=header), self.assertRaises(ValueError):
                uncertain.assembly_status(header)


if __name__ == "__main__":
    unittest.main()

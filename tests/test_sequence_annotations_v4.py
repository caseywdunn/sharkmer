#!/usr/bin/env python3
"""Schema-4 field migration, historical reading and per-record evidence joins."""

import copy
import tempfile
import unittest
from pathlib import Path

from test_sequence_annotations import warning_transaction
from test_uncertain_output import (
    collect_amplicons_from_runs,
    combined_transaction,
    publish,
    runner,
    transaction,
)


def annotation_transaction(
    directory: Path,
    standard: bool = True,
    empty: bool = False,
    disabled: bool = False,
    prefix: str = "sample",
    gene: str = "gene",
) -> dict:
    """Emit schema 4 while keeping historical fixture factories unchanged."""
    stats = warning_transaction(directory, standard, empty, disabled, prefix, gene)
    stats["product_output_schema_version"] = 4
    entry = stats["pcr_results"][0]
    rows, contents = [], ""
    if standard:
        rows.append(
            {
                "sequence_id": f"{prefix}_{gene}_0",
                "sequence_index": 0,
                "sequence_length": 4,
                "assembly_rank_score": 7.0,
            }
        )
        contents = (
            f">{prefix}_{gene}_0 sample={prefix} target={gene} sequence_index=0 "
            "sequence_length=4 kmer_count_mean=7.00 kmer_count_median=7 "
            "kmer_count_min=7 kmer_count_max=7\nACGT\n"
        )
    value = entry.get("uncertain_products")
    if value is not None:
        value["schema_version"] = 4
        for index, product in enumerate(value["products"]):
            product["description"] = (
                f"sample={prefix} target={gene} sequence_index={index} "
                f"sequence_length={product['sequence_length']} assembly_status=warning "
                "assembly_status_causes=repeat_structure_unresolved"
            )
            rows.append(
                {
                    "sequence_id": product["id"],
                    "sequence_index": index,
                    "sequence_length": product["sequence_length"],
                }
            )
            contents += f">{product['id']} {product['description']}\n{product['oriented_sequence']}\n"
    if rows:
        entry["sequence_annotations"] = rows
    filename = entry.get("output_file")
    if filename:
        (directory / filename).write_text(contents)
    publish(directory, stats, [filename] if filename else [])
    return stats


class SequenceAnnotationMigrationTests(unittest.TestCase):
    def test_all_historical_formats_and_schema4_preserve_ordinary_recovery(self):
        for factory in (
            transaction,
            combined_transaction,
            warning_transaction,
            annotation_transaction,
        ):
            for ordinary in (False, True):
                with (
                    self.subTest(factory=factory.__name__, ordinary=ordinary),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    directory = Path(temporary)
                    stats = factory(directory, standard=ordinary)
                    runner._validate_stats_manifest(
                        stats, "sample", 19, {"gene"}, "end-to-end"
                    )
                    parsed = runner.parse_fasta_products("sample", directory)
                    self.assertEqual(
                        [seq for row in parsed for seq in row["sequences"]],
                        ["ACGT"] if ordinary else [],
                    )

    def test_schema4_mixed_warning_only_empty_and_opt_out_have_exact_joins(self):
        for ordinary in (False, True):
            for empty, disabled in ((False, False), (True, False), (True, True)):
                with (
                    self.subTest(ordinary=ordinary, empty=empty, disabled=disabled),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    directory = Path(temporary)
                    stats = annotation_transaction(directory, ordinary, empty, disabled)
                    parsed = runner.parse_fasta_products("sample", directory)
                    entry = stats["pcr_results"][0]
                    rows = entry.get("sequence_annotations", [])
                    self.assertEqual(len(rows), entry["n_products"])
                    self.assertEqual(
                        sum("assembly_rank_score" in row for row in rows), int(ordinary)
                    )
                    if ordinary:
                        record = parsed[0]["products"][0]
                        self.assertIn("sequence_index=0", record["header"])
                        self.assertNotIn("score=", record["header"])
                        self.assertEqual(
                            rows[0]["sequence_id"], record["header"].split()[0]
                        )
                        self.assertEqual(rows[0]["assembly_rank_score"], 7.0)
                    self.assertEqual(
                        sum(row["n_products"] for row in parsed), int(ordinary)
                    )

    def test_bad_evidence_joins_and_invented_warning_scores_fail_with_valid_receipts(
        self,
    ):
        cases = [
            (0, "sequence_id", "other"),
            (0, "sequence_index", 1),
            (0, "sequence_index", True),
            (0, "sequence_length", 5),
            (0, "assembly_rank_score", "7"),
            (0, "assembly_rank_score", True),
            (0, "assembly_rank_score", float("nan")),
            (0, "assembly_rank_score", float("inf")),
            (0, "assembly_rank_score", -1),
            (1, "assembly_rank_score", 0),
        ]
        for index, field, value in cases:
            with (
                self.subTest(index=index, field=field, value=value),
                tempfile.TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                stats = annotation_transaction(directory)
                stats["pcr_results"][0]["sequence_annotations"][index][field] = value
                publish(directory, stats, ["sample_gene.fasta"])
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)
        for mutation in ("absent", "missing_score", "duplicate", "reordered"):
            with (
                self.subTest(mutation=mutation),
                tempfile.TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                stats = annotation_transaction(directory)
                entry = stats["pcr_results"][0]
                if mutation == "absent":
                    del entry["sequence_annotations"]
                elif mutation == "missing_score":
                    del entry["sequence_annotations"][0]["assembly_rank_score"]
                elif mutation == "duplicate":
                    entry["sequence_annotations"][1] = copy.deepcopy(
                        entry["sequence_annotations"][0]
                    )
                else:
                    entry["sequence_annotations"].reverse()
                publish(directory, stats, ["sample_gene.fasta"])
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_ordinary_header_tampering_and_old_new_alias_mixtures_fail(self):
        for old, new in [
            ("sample=sample", "sample=other"),
            ("target=gene", "target=other"),
            ("sequence_index=0", "sequence_index=1"),
            ("sequence_length=4", "sequence_length=5"),
            ("sequence_index=0", "sequence_index=0 product=0"),
            ("target=gene", "target=gene gene=gene"),
            ("sequence_length=4", "sequence_length=4 length=4"),
            ("kmer_count_mean=7.00", "kmer_count_mean=7.00 score=7.00"),
            ("target=gene", "target=gene read_support=evaluated"),
        ]:
            with self.subTest(new=new), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = annotation_transaction(directory)
                path = directory / "sample_gene.fasta"
                path.write_text(path.read_text().replace(old, new, 1))
                publish(directory, stats, [path.name])
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_standalone_reader_uses_new_index_and_preserves_identifier_with_equals(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "sample_gene.fasta"
            path.write_text(
                ">sample=x sample=sample target=gene sequence_index=7 sequence_length=4\nACGT\n"
                ">warning target=gene sequence_index=0 sequence_length=5 assembly_status=warning "
                "assembly_status_causes=repeat_structure_unresolved\nACGTT\n"
            )
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual(parsed[0]["sequences"], ["ACGT"])
            self.assertEqual(parsed[0]["products"][0]["product_index"], 7)

    def test_bootstrap_never_promotes_schema4_warning_hypotheses(self):
        panel = {
            "name": "panel",
            "primers": [{"gene": "gene", "forward": "ACGT", "reverse": "TGCA"}],
            "validation": {
                "samples": [{"accession": "A", "taxon": "Taxon", "max_reads": [1000]}]
            },
        }
        for ordinary in (False, True):
            with (
                self.subTest(ordinary=ordinary),
                tempfile.TemporaryDirectory() as temporary,
            ):
                directory = Path(temporary)
                annotation_transaction(
                    directory, standard=ordinary, prefix="panel_A_1k", gene="panel_gene"
                )
                rows = collect_amplicons_from_runs(panel, directory, "panel")
                self.assertEqual(
                    [row["sequence"] for row in rows], ["ACGT"] if ordinary else []
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)

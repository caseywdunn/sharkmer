#!/usr/bin/env python3
"""Schema-3 warning annotations and legacy compatibility regression controls."""

import copy
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import yaml

from test_uncertain_output import (
    collect_amplicons_from_runs,
    combined_transaction,
    publish,
    runner,
    transaction,
    uncertain,
)


def warning_transaction(
    directory,
    standard=True,
    empty=False,
    disabled=False,
    prefix="sample",
    gene="gene",
    marker_detail="complete",
):
    """Produce the new wire format without altering the legacy fixture helpers."""
    stats = combined_transaction(
        directory, standard, empty, disabled, prefix=prefix, gene=gene
    )
    stats["product_output_schema_version"] = 3
    stats["product_output_layout"] = "combined_warnings"
    entry = stats["pcr_results"][0]
    value = entry.get("uncertain_products")
    contents = (
        f">{prefix}_{gene}_0 product=0 length=4 kmer_count_median=7 score=7.00\nACGT\n"
        if standard
        else ""
    )
    if value is not None:
        value["schema_version"] = 3
        value["marker_span_scope"] = "marker_covered_positions_only_for_complete_detail_sources"
        value["marker_detail_truncated"] = marker_detail == "summary_only" and not empty
        for product in value["products"]:
            product["assembly_status_causes"] = ["repeat_structure_unresolved"]
            product["description"] = (
                f"sample={prefix} gene={gene} product=0 length={product['sequence_length']} "
                "assembly_status=warning assembly_status_causes=repeat_structure_unresolved"
            )
            for source in product["sources"]:
                source["marker_detail"] = marker_detail
                if marker_detail == "summary_only":
                    source["marker_occurrences"] = []
                    source["marker_covered_runs"] = []
            contents += (
                f">{product['id']} {product['description']}\n{product['oriented_sequence']}\n"
            )
    filename = entry.get("output_file")
    if filename:
        (directory / filename).write_text(contents)
    publish(directory, stats, [filename] if filename else [])
    return stats


class SequenceAnnotationTests(unittest.TestCase):
    def test_legacy_and_new_formats_preserve_standard_recovery(self):
        for factory in (transaction, combined_transaction, warning_transaction):
            for standard in (False, True):
                with self.subTest(factory=factory.__name__, standard=standard):
                    with tempfile.TemporaryDirectory() as temporary:
                        directory = Path(temporary)
                        stats = factory(directory, standard=standard)
                        runner._validate_stats_manifest(
                            stats, "sample", 19, {"gene"}, "end-to-end"
                        )
                        parsed = runner.parse_fasta_products("sample", directory)
                        self.assertEqual(
                            [sequence for row in parsed for sequence in row["sequences"]],
                            ["ACGT"] if standard else [],
                        )

    def test_new_mixed_warning_only_empty_and_opt_out_outputs(self):
        for standard in (False, True):
            for empty, disabled in ((False, False), (True, False), (True, True)):
                for detail in ("complete", "summary_only"):
                    with self.subTest(
                        standard=standard, empty=empty, disabled=disabled, detail=detail
                    ), tempfile.TemporaryDirectory() as temporary:
                        directory = Path(temporary)
                        stats = warning_transaction(
                            directory, standard, empty, disabled, marker_detail=detail
                        )
                        runner._validate_stats_manifest(
                            stats, "sample", 19, {"gene"}, "end-to-end"
                        )
                        parsed = runner.parse_fasta_products("sample", directory)
                        self.assertEqual(sum(row["n_products"] for row in parsed), int(standard))
                        if standard:
                            record = parsed[0]["products"][0]
                            self.assertEqual(record["product_index"], 0)
                            self.assertEqual(record["kmer_count_median"], 7)
                            self.assertIn("score=7.00", record["header"])
                            self.assertNotIn("assembly_status", record["header"])
                        if not empty:
                            header = stats["pcr_results"][0]["uncertain_products"]["products"][0][
                                "description"
                            ]
                            self.assertNotIn("read_support=", header)
                            self.assertNotIn("retention_truncated=", header)

    def test_missing_or_invalid_status_does_not_fall_back_between_schemas(self):
        for schema in (2, 3):
            invalid = (
                "assembly_status=other",
                "assembly_status=standard assembly_status=warning",
                "assembly_status=standard assembly_status=standard",
                "assembly_status_causes=repeat_structure_unresolved",
                "assembly_status",
            )
            invalid += ("", "assembly_status=warning") if schema == 2 else (
                "assembly_status=standard", "assembly_status=uncertain"
            )
            for annotation in invalid:
                with self.subTest(schema=schema, annotation=annotation):
                    with self.assertRaises(ValueError):
                        uncertain.assembly_status(f"record {annotation}", schema)
        self.assertEqual(uncertain.assembly_status("record length=4", 3), "standard")

    def test_warning_header_tampering_rejected_with_recomputed_receipts(self):
        mutations = (
            ("assembly_status=warning", ""),
            ("assembly_status=warning", "assembly_status=standard"),
            ("assembly_status=warning", "assembly_status=uncertain"),
            ("assembly_status=warning", "assembly_status=unsupported"),
            ("assembly_status=warning", "assembly_status=warning assembly_status=warning"),
            ("assembly_status_causes=repeat_structure_unresolved", ""),
            ("assembly_status_causes=repeat_structure_unresolved", "assembly_status_causes="),
            ("assembly_status_causes=repeat_structure_unresolved", "assembly_status_causes=unknown"),
            (
                "assembly_status_causes=repeat_structure_unresolved",
                "assembly_status_causes=repeat_structure_unresolved,repeat_structure_unresolved",
            ),
            ("length=5", "length=6"),
            ("sample=sample", "sample=other"),
            ("gene=gene", "gene=other"),
            ("product=0", "product=1"),
            ("assembly_status=warning", "assembly_status=warning read_support=evaluated"),
            ("assembly_status=warning", "assembly_status=warning uncertainty=repeat_structure"),
        )
        for old, new in mutations:
            with self.subTest(old=old, new=new), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = warning_transaction(directory)
                product = stats["pcr_results"][0]["uncertain_products"]["products"][0]
                original = product["description"]
                product["description"] = original.replace(old, new)
                path = directory / "sample_gene.fasta"
                path.write_text(path.read_text().replace(original, product["description"]))
                publish(directory, stats, [path.name])
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_assessment_and_manifest_versions_must_match_output_contract(self):
        for mutation in ("assessment", "manifest", "layout", "schema", "downgrade", "missing"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = warning_transaction(directory)
                if mutation == "assessment":
                    stats["pcr_results"][0]["uncertain_products"]["schema_version"] = 2
                elif mutation == "layout":
                    stats["product_output_layout"] = "combined_classified"
                elif mutation == "schema":
                    stats["product_output_schema_version"] = 4
                elif mutation == "downgrade":
                    stats["product_output_schema_version"] = 2
                    stats["product_output_layout"] = "combined_classified"
                elif mutation == "missing":
                    del stats["product_output_schema_version"]
                manifest = publish(directory, stats, ["sample_gene.fasta"])
                if mutation == "manifest":
                    manifest["schema_version"] = 2
                    (directory / "sample.manifest.yaml").write_text(yaml.safe_dump(manifest))
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_read_support_cannot_be_claimed_on_unassessed_ordinary_record(self):
        with self.assertRaises(ValueError):
            uncertain.assembly_status("record length=4 read_support=evaluated", 3)

    def test_semantic_causes_are_bound_to_internal_source_evidence(self):
        groups = [(cause,) for cause in uncertain.CAUSE_COUNTS]
        groups.append(tuple(uncertain.CAUSE_COUNTS))
        for causes in groups:
            with self.subTest(causes=causes), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = warning_transaction(directory, marker_detail="summary_only")
                entry = stats["pcr_results"][0]
                product = entry["uncertain_products"]["products"][0]
                source = product["sources"][0]
                for cause, field in uncertain.CAUSE_COUNTS.items():
                    source[field] = int(cause in causes)
                expected = sorted(uncertain.WARNING_CAUSES[cause] for cause in causes)
                product["assembly_status_causes"] = expected
                product["description"] = (
                    "sample=sample gene=gene product=0 length=5 assembly_status=warning "
                    f"assembly_status_causes={','.join(expected)}"
                )
                self.assertEqual(
                    uncertain.validate_combined_metadata(entry, "sample", 3),
                    "sample_gene.fasta",
                )
                for wrong in ([], ["unknown"], expected + expected, list(reversed(expected))):
                    if wrong == expected:
                        continue
                    bad = copy.deepcopy(entry)
                    bad["uncertain_products"]["products"][0]["assembly_status_causes"] = wrong
                    with self.assertRaises(ValueError):
                        uncertain.validate_combined_metadata(bad, "sample", 3)

    def test_standalone_renamed_warning_records_are_never_promoted(self):
        ordinary = ">ordinary product=0\nACGT\n"
        for status in ('assembly_status=warning', 'assembly_status="warning"'):
            warning = (
                f">renamed {status} assembly_status_causes=repeat_structure_unresolved\nACGTT\n"
            )
            for content in (warning, ordinary + warning, warning + ordinary):
                with self.subTest(content=content), tempfile.TemporaryDirectory() as temporary:
                    directory = Path(temporary)
                    (directory / "sample_gene.fasta").write_text(content)
                    parsed = runner.parse_fasta_products("sample", directory)
                    observed = [sequence for row in parsed for sequence in row["sequences"]]
                    self.assertEqual(observed, ["ACGT"] if ordinary in content else [])
                    with self.assertRaises(ValueError):
                        runner.parse_fasta_products("sample", directory, ["sample_gene.fasta"])

    def test_standalone_unknown_duplicate_and_orphaned_metadata_rejected(self):
        cases = (
            "assembly_status=unknown",
            "assembly_status=warning",
            "assembly_status=warning assembly_status_causes=unknown",
            "assembly_status=warning assembly_status=standard",
            "assembly_status_causes=repeat_structure_unresolved",
            "assembly_status=standard causes=pre_prune_cyclic_scc_node",
            "assembly_status=warning assembly_status_causes=repeat_structure_unresolved "
            "assembly_status_causes=repeat_structure_unresolved",
        )
        for annotation in cases:
            with self.subTest(annotation=annotation), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                (directory / "sample_gene.fasta").write_text(f">candidate {annotation}\nACGTT\n")
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_legacy_identifier_with_equals_is_not_an_annotation(self):
        ordinary = ">sample=x sample=s gene=g product=0 length=4\nACGT\n"
        warning = (
            ">sample=y sample=s gene=g product=1 length=5 assembly_status=warning "
            "assembly_status_causes=repeat_structure_unresolved\nACGTT\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "sample_gene.fasta"
            path.write_text(ordinary + warning)
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual(parsed[0]["sequences"], ["ACGT"])
            self.assertTrue(parsed[0]["products"][0]["header"].startswith("sample=x "))
            path.write_text(ordinary.replace("sample=s", "sample=s sample=t"))
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_warning_annotation_in_identifier_position_is_rejected(self):
        identifiers = (
            "assembly_status=warning", 'assembly_status="warning"',
            "'assembly_status=warning'", "assembly_status=uncertain",
            '"assembly_status"=uncertain', "'assembly_status'=uncertain",
            '"assembly_status"="warning"',
            "assembly_status_causes=repeat_structure_unresolved",
            '"assembly_status_causes"=repeat_structure_unresolved',
            "causes=pre_prune_cyclic_scc_node", "repeat_length=unresolved",
            "uncertainty=repeat_structure",
            "'assembly_status=warning", 'assembly_status="warning',
        )
        for identifier in identifiers:
            with self.subTest(identifier=identifier), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                (directory / "sample_gene.fasta").write_text(f">{identifier} length=5\nACGTT\n")
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_bootstrap_excludes_warning_records(self):
        panel = {
            "name": "panel",
            "primers": [{"gene": "gene", "forward": "ACGT", "reverse": "TGCA"}],
            "validation": {
                "samples": [{"accession": "A", "taxon": "Taxon", "max_reads": [1000]}]
            },
        }
        for standard in (False, True):
            with self.subTest(standard=standard), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                warning_transaction(
                    directory, standard=standard, prefix="panel_A_1k", gene="panel_gene"
                )
                rows = collect_amplicons_from_runs(panel, directory, "panel")
                self.assertEqual([row["sequence"] for row in rows], ["ACGT"] if standard else [])

    def test_current_run_recovery_and_bootstrap_use_only_unflagged_partition(self):
        for standard in (False, True):
            with self.subTest(standard=standard), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                panel = {
                    "name": "panel",
                    "primers": [{"gene": "gene", "forward": "ACGT", "reverse": "TGCA"}],
                    "validation": {
                        "samples": [{"accession": "A", "taxon": "Taxon", "max_reads": [1000]}]
                    },
                }
                panel_path = directory / "panel.yaml"
                panel_path.write_text(yaml.safe_dump(panel))
                fastq = directory / "reads.fastq"
                fastq.write_text("@read\nACGT\n+\nIIII\n")

                def simulated_process(command, resource_path):
                    output_dir = Path(command[command.index("-o") + 1])
                    prefix = command[command.index("-s") + 1]
                    stats = warning_transaction(
                        output_dir, standard=standard, prefix=prefix, gene="panel_gene"
                    )
                    stats.update(
                        command=" ".join(command),
                        n_reads_read=1,
                        n_bases_read=4,
                        input_source={
                            "kind": "local_files", "inputs": [str(fastq)],
                            "paired": False, "max_reads": 1000,
                        },
                    )
                    publish(output_dir, stats, [f"{prefix}_panel_gene.fasta"])
                    return subprocess.CompletedProcess(command, 0, "", ""), None, "synthetic"

                with mock.patch.object(runner, "_run_with_rss", side_effect=simulated_process):
                    run = runner.run_sharkmer(
                        panel_path, "panel", "A", 1000, directory / "runs",
                        executable=directory / "unused", input_path=fastq,
                    )
                self.assertTrue(run["success"], run.get("failure"))
                self.assertEqual(run["genes"][0]["recovered"], standard)
                self.assertEqual(run["genes"][0]["n_products"], int(standard))
                rows = collect_amplicons_from_runs(
                    panel, directory / "runs", "panel", {("A", 1000): run}
                )
                self.assertEqual([row["sequence"] for row in rows], ["ACGT"] if standard else [])


if __name__ == "__main__":
    unittest.main(verbosity=2)

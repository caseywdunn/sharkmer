import sys
import unittest
import copy
import json
import random
import shutil
import tempfile
from dataclasses import asdict
from pathlib import Path
from unittest import mock

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from sharkmer_validate import blast_references, report, reference_provenance
from sharkmer_validate.reference_targets import (
    logical_gene_name,
    reference_target_names,
    target_logical_genes,
)


def reviewed_its_panel(name="insecta"):
    if name == "c_elegans":
        first_forward = "TACACACCGCCCGTCGCTATCC"
        second_forward = "GTAGGTGAACCTGCAGCTGGATCA"
        reverse = "ACTCGCCGTTACTAAGG"
    else:
        first_forward = "TACACACCGCCCGTCGCTACTA"
        second_forward = "GTAGGTGAACCTGCAGAAGGATCA"
        reverse = "ACTCGCCGTTACTRRGG"
    return {
        "name": name,
        "primers": [
            {
                "gene": "ITS",
                "index": 1,
                "forward_seq": first_forward,
                "reverse_seq": reverse,
            },
            {
                "gene": "ITS",
                "index": 2,
                "forward_seq": second_forward,
                "reverse_seq": reverse,
            },
        ],
    }


def fixture_xml(targets):
    hits = []
    for target in targets:
        hits.append(
            f"<Hit><Hit_def>{target}|Test_taxon|TEST.1</Hit_def><Hit_len>100</Hit_len>"
            "<Hit_hsps><Hsp><Hsp_bit-score>200</Hsp_bit-score><Hsp_identity>100</Hsp_identity>"
            "<Hsp_align-len>100</Hsp_align-len><Hsp_gaps>0</Hsp_gaps>"
            "<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>100</Hsp_query-to>"
            "<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>100</Hsp_hit-to>"
            "</Hsp></Hit_hsps></Hit>"
        )
    return (
        "<BlastOutput><BlastOutput_iterations><Iteration><Iteration_query-len>100</Iteration_query-len>"
        "<Iteration_hits>" + "".join(hits)
        + "</Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"
    )


class ReferenceTargetTests(unittest.TestCase):
    def classify(self, expected, targets):
        return blast_references._parse_blast_xml(
            fixture_xml(targets), expected, "Test taxon",
            fixture_provenance="explicit_test_fixture",
        )

    def test_reviewed_same_gene_indices_share_support_not_primer_region_truth(self):
        for gene in ("16S", "18S", "28S", "CO1", "CO2"):
            with self.subTest(gene=gene):
                match = self.classify(f"{gene}_2", [f"{gene}_1", f"{gene}_2"])
                self.assertEqual(match.status, "gene_supported_expected_taxon")
                self.assertEqual(match.expected_gene, f"{gene}_2")
                self.assertEqual(match.matched_gene, f"{gene}_1")
                self.assertEqual(match.expected_logical_gene, gene.casefold())
                self.assertEqual(match.matched_logical_gene, gene.casefold())
                self.assertEqual(match.primer_region_support, "not_established")

    def test_distinct_genes_and_regions_are_not_merged(self):
        for expected, other in (
            ("18S_1", "28S_1"),
            ("ITS_1", "ITS_2"),
            ("mt-1404-3947", "mt-3734-6739"),
            ("atpB-rbcL", "atpB"),
            ("custom_1", "custom_2"),
            ("CO1_1", "CO1_12"),
        ):
            with self.subTest(expected=expected, other=other):
                self.assertEqual(self.classify(expected, [expected, other]).status, "ambiguous_gene")
                self.assertEqual(self.classify(expected, [other]).status, "wrong_gene")

    def test_gene_equivalence_is_case_insensitive_and_explicit(self):
        self.assertEqual(logical_gene_name("cO1_1"), "co1")
        self.assertEqual(logical_gene_name("ITS_1"), "its_1")
        self.assertEqual(logical_gene_name("18S-V9"), "18s-v9")
        self.assertEqual(logical_gene_name("psbA-trnH"), "psba-trnh")

    def test_panel_mapping_preserves_region_identity(self):
        mapping = target_logical_genes({
            "name": "external_panel",
            "primers": [
                {"gene": "mt", "region": "1404-3947"},
                {"gene": "mt", "region": "3734-6739"},
                {"gene": "atpB", "region": "rbcL"},
            ],
        })
        self.assertEqual(mapping["mt-1404-3947"], "mt-1404-3947")
        self.assertEqual(mapping["mt-3734-6739"], "mt-3734-6739")
        self.assertEqual(mapping["atpB-rbcL"], "atpb-rbcl")

    def test_report_scores_gene_support_from_another_reviewed_primer_index(self):
        match = asdict(self.classify("18S_2", ["18S_1"]))
        self.assertEqual(
            report._score_gene(True, "18S_2", "Test taxon", match, {"18S_1": {"Test taxon"}}),
            "+**",
        )

    def test_product_evaluation_uses_logical_gene_availability(self):
        runs = [{"success": True, "genes": [{"gene": "18S_2", "products": [{"sequence": "ACGT"}]}]}]
        match = self.classify("18S_2", ["18S_1"])
        manifest = {
            "reference_audit": {
                "target_logical_genes": {"18S_2": "18s"},
                "verified": [{"logical_gene": "18s"}],
            }
        }
        with mock.patch.object(
            blast_references, "_load_database_manifest", return_value=(manifest, {})
        ), mock.patch.object(
            blast_references, "blast_against_references", return_value=match
        ) as command:
            blast_references.blast_all_products(
                runs,
                Path("db"),
                "Test taxon",
                reference_genes={"18S_1"},
                target_mapping={"18S_2": "18s"},
            )
        command.assert_called_once()
        self.assertEqual(runs[0]["genes"][0]["products"][0]["reference_match"]["target_support"], "supported")

    def test_reviewed_panel_its_targets_share_only_logical_gene_support(self):
        for panel_name in ("insecta", "cnidaria", "c_elegans"):
            with self.subTest(panel=panel_name):
                panel = reviewed_its_panel(panel_name)
                mapping = target_logical_genes(panel)
                self.assertEqual(mapping["ITS_1"], "its_rdna_cluster")
                self.assertEqual(mapping["ITS_2"], "its_rdna_cluster")
                self.assertEqual(
                    reference_target_names(panel, {"ITS_1"}),
                    {"ITS_1", "ITS_2"},
                )
                metadata = {
                    "ITS_1|Test_taxon|TEST.1": {
                        "gene": "ITS_1",
                        "logical_gene": "its_rdna_cluster",
                        "taxon": "Test taxon",
                        "accession": "TEST.1",
                        "length": 100,
                        "provenance": {"status": "explicit_test_fixture"},
                    }
                }
                match = blast_references._parse_blast_xml(
                    fixture_xml(["ITS_1"]),
                    "ITS_2",
                    "Test taxon",
                    reference_metadata=metadata,
                    fixture_provenance="explicit_test_fixture",
                    expected_logical_gene=mapping["ITS_2"],
                )
                self.assertEqual(match.status, "gene_supported_expected_taxon")
                self.assertEqual(match.expected_gene, "ITS_2")
                self.assertEqual(match.matched_gene, "ITS_1")
                self.assertEqual(match.primer_region_support, "not_established")

    def test_registered_its_signatures_match_repository_panels(self):
        panels_directory = Path(__file__).resolve().parents[1] / "panels"
        for panel_name in ("insecta", "cnidaria", "c_elegans"):
            with self.subTest(panel=panel_name):
                panel = yaml.safe_load(
                    (panels_directory / f"{panel_name}.yaml").read_text()
                )
                mapping = target_logical_genes(panel)
                self.assertEqual(mapping["ITS_1"], "its_rdna_cluster")
                self.assertEqual(mapping["ITS_2"], "its_rdna_cluster")

    def test_its_grouping_requires_reviewed_panel_and_exact_primer_context(self):
        panel = reviewed_its_panel()
        mutated = copy.deepcopy(panel)
        mutated["primers"][1]["forward_seq"] += "A"
        external = copy.deepcopy(panel)
        external["name"] = "external_panel"
        for candidate in (mutated, external):
            with self.subTest(panel=candidate["name"]):
                mapping = target_logical_genes(candidate)
                self.assertEqual(mapping["ITS_1"], "its_1")
                self.assertEqual(mapping["ITS_2"], "its_2")
                self.assertNotEqual(mapping["ITS_1"], mapping["ITS_2"])

    def test_context_receipt_changes_after_primer_mutation(self):
        panel = reviewed_its_panel()
        mutated = copy.deepcopy(panel)
        mutated["primers"][1]["reverse_seq"] += "A"
        empty_audit = {
            "verified": [],
            "excluded": [],
            "verified_count": 0,
            "excluded_count": 0,
            "total_count": 0,
            "catalog": {"status": "verified_local_snapshot"},
            "audit_sha256": "a" * 64,
        }
        with mock.patch.object(
            reference_provenance, "audit_references", return_value=empty_audit
        ):
            reviewed = blast_references.reference_checksums(panel)
            changed = blast_references.reference_checksums(mutated)
        self.assertNotEqual(
            reviewed["target_logical_genes_sha256"],
            changed["target_logical_genes_sha256"],
        )
        self.assertEqual(
            reviewed["target_logical_genes"]["ITS_2"], "its_rdna_cluster"
        )
        self.assertEqual(changed["target_logical_genes"]["ITS_2"], "its_2")

    def test_stale_database_mapping_is_not_used_after_primer_mutation(self):
        runs = [{
            "success": True,
            "genes": [{
                "gene": "ITS_2",
                "products": [{"product_index": 0, "sequence": "A" * 100}],
            }],
        }]
        manifest = {
            "reference_audit": {
                "target_logical_genes": {
                    "ITS_1": "its_rdna_cluster",
                    "ITS_2": "its_rdna_cluster",
                },
                "verified": [{"logical_gene": "its_rdna_cluster"}],
            }
        }
        with mock.patch.object(
            blast_references, "_load_database_manifest", return_value=(manifest, {})
        ), mock.patch.object(blast_references, "blast_against_references") as command:
            blast_references.blast_all_products(
                runs,
                Path("db"),
                "Test taxon",
                target_mapping={"ITS_1": "its_1", "ITS_2": "its_2"},
            )
        match = runs[0]["genes"][0]["products"][0]["reference_match"]
        self.assertEqual(match["status"], "failed_run")
        self.assertIn("differs from current panel context", match["error"])
        command.assert_not_called()

    def test_metadata_audit_rejects_logical_target_inconsistent_with_context(self):
        metadata = {
            "reference_000000": {
                "gene": "ITS_1",
                "logical_gene": "its_1",
                "taxon": "Test taxon",
                "accession": "TEST.1",
                "length": 100,
                "sha256": "a" * 64,
                "provenance": {},
            }
        }
        manifest = {
            "reference_audit": {
                "target_logical_genes": {"ITS_1": "its_rdna_cluster"},
                "verified": [{
                    "gene": "ITS_1",
                    "logical_gene": "its_1",
                    "taxon": "Test taxon",
                    "accession": "TEST.1",
                    "length": 100,
                    "sha256": "a" * 64,
                    "provenance": None,
                    "source_taxid": None,
                    "gene_assignment_basis": None,
                    "contains_ambiguity": None,
                }],
            }
        }
        with self.assertRaisesRegex(ValueError, "differs from panel context"):
            blast_references._validate_metadata_audit(metadata, manifest)

    def test_report_uses_bound_logical_target_without_claiming_primer_support(self):
        match = asdict(
            blast_references.RefBlastResult(
                "gene_supported_expected_taxon",
                "ITS_2",
                "Test taxon",
                matched_gene="ITS_1",
                matched_taxon="Test taxon",
                target_support="supported",
                taxon_support="supported",
                expected_logical_gene="its_rdna_cluster",
                matched_logical_gene="its_rdna_cluster",
            )
        )
        self.assertEqual(
            report._score_gene(
                True,
                "ITS_2",
                "Test taxon",
                match,
                {"ITS_1": {"Test taxon"}, "ITS_2": {"Test taxon"}},
            ),
            "+**",
        )
        result = {
            "provenance": {
                "references": {
                    "verified_count": 1,
                    "excluded_count": 0,
                    "catalog": {"status": "verified_local_snapshot"},
                    "target_logical_genes": {
                        "ITS_1": "its_rdna_cluster",
                        "ITS_2": "its_rdna_cluster",
                    },
                }
            }
        }
        rendered = "\n".join(report._reference_provenance_summary(result))
        self.assertIn("Reviewed logical target groups", rendered)
        self.assertIn("Primer-region support**: not established", rendered)
        details = report._reference_details(
            {
                "samples": [{
                    "accession": "TEST",
                    "taxon": "Test taxon",
                    "depths": [{
                        "success": True,
                        "max_reads": 100,
                        "genes": [{
                            "gene": "ITS_2",
                            "logical_gene": "its_rdna_cluster",
                            "recovered": True,
                            "products": [{
                                "product_index": 0,
                                "reference_match": match,
                            }],
                        }],
                    }],
                }]
            },
            ["ITS_2"],
        )
        rendered_details = "\n".join(details)
        self.assertIn("Matched reference target", rendered_details)
        self.assertIn("| ITS_2 | its_rdna_cluster | 0 |", rendered_details)
        self.assertIn("| ITS_1 | its_rdna_cluster | Test taxon |", rendered_details)
        self.assertIn("| not_established | not_established | not_evaluated |", rendered_details)

    @unittest.skipUnless(shutil.which("blastn") and shutil.which("makeblastdb"), "BLAST+ not installed")
    def test_real_database_build_and_query_binds_logical_gene_metadata(self):
        generator = random.Random(1703)
        sequence = "".join(generator.choice("ACGT") for _ in range(240))
        digest = reference_provenance.sequence_digest(sequence)
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            catalog_path = temporary / "catalog.json"
            record = {
                "accession_version": "SYNTHETIC.1", "sequence": sequence,
                "sequence_sha256": digest, "organism": "Test taxon", "taxid": 1,
                "topology": "linear", "url": "https://example.invalid/synthetic-test",
                "retrieved_at": "2026-09-13T00:00:00Z",
            }
            catalog_path.write_text(json.dumps({"schema_version": 1, "records": {"SYNTHETIC.1": record}}))
            provenance = {
                "schema_version": 1, "kind": "public_record_region", "accession_version": "SYNTHETIC.1",
                "source_sequence_sha256": digest, "source_length": len(sequence), "start": 0,
                "end": len(sequence), "strand": "+", "wraps_origin": False, "sequence_sha256": digest,
            }
            reference = {"taxon": "Test taxon", "accession": "SYNTHETIC.1", "sequence": sequence, "provenance": provenance}
            panel = {"references": [{"gene": "18S_1", "sequences": [reference]}]}
            database = blast_references.build_reference_db(panel, temporary, catalog_path)
            self.assertIsNotNone(database)
            result = blast_references.blast_against_references(sequence, database, "18S_2", "Test taxon")
            self.assertEqual(result.status, "gene_supported_expected_taxon", result.error)
            self.assertEqual(result.matched_logical_gene, "18s")


if __name__ == "__main__":
    unittest.main()

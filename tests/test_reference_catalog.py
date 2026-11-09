import copy
import gzip
import json
import sys
import tempfile
import unittest
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from sharkmer_validate import reference_provenance as provenance
import bootstrap_references


class ReferenceCatalogTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name) / "catalog.json.gz"
        self.source = "ACGTRYSWKMBDHVN"
        self.record = {
            "accession_version": "TEST_123.1",
            "sequence": self.source,
            "sequence_sha256": provenance.sequence_digest(self.source),
            "organism": "Test taxon",
            "taxid": 1,
            "topology": "linear",
            "url": "https://www.ncbi.nlm.nih.gov/nuccore/TEST_123.1",
            "retrieved_at": "2026-09-13T00:00:00Z",
        }
        self.catalog = {"schema_version": 1, "records": {"TEST_123.1": self.record}}
        self.reference = {
            "accession": "TEST_123.1",
            "taxon": "Test taxon",
            "sequence": self.source[2:8],
            "provenance": {
                "schema_version": 1,
                "kind": "public_record_region",
                "accession_version": "TEST_123.1",
                "source_sequence_sha256": self.record["sequence_sha256"],
                "source_length": len(self.source),
                "start": 2,
                "end": 8,
                "strand": "+",
                "wraps_origin": False,
                "sequence_sha256": provenance.sequence_digest(self.source[2:8]),
            },
        }

    def write_catalog(self, content=None):
        serialized = json.dumps(self.catalog).encode() if content is None else content
        self.path.write_bytes(gzip.compress(serialized, mtime=0))

    def audit(self, reference=None):
        panel = {"references": [{"gene": "target", "sequences": [reference or self.reference]}]}
        return provenance.audit_references(panel, self.path)

    def test_exact_public_region_is_verified_without_imputing_iupac(self):
        self.write_catalog()
        result = self.audit()
        self.assertEqual(result["verified_count"], 1)
        self.assertEqual(result["excluded_count"], 0)
        self.assertEqual(result["verified"][0]["sequence"], self.source[2:8])
        self.assertTrue(result["verified"][0]["contains_ambiguity"])
        self.assertEqual(result["catalog"]["status"], "verified_local_snapshot")

    def test_reverse_complement_and_coordinates_are_exact(self):
        self.reference["provenance"]["strand"] = "-"
        sequence = provenance.reverse_complement(self.source[2:8])
        self.reference["sequence"] = sequence
        self.reference["provenance"]["sequence_sha256"] = provenance.sequence_digest(sequence)
        self.write_catalog()
        self.assertEqual(self.audit()["verified_count"], 1)
        self.reference["provenance"]["start"] = 3
        self.assertEqual(self.audit()["excluded_count"], 1)

    def test_circular_origin_requires_topology_and_at_most_one_traversal(self):
        self.record["topology"] = "circular"
        region = self.source[10:] + self.source[:3]
        self.reference["sequence"] = region
        self.reference["provenance"].update(
            start=10, end=3, wraps_origin=True,
            sequence_sha256=provenance.sequence_digest(region),
        )
        self.write_catalog()
        self.assertEqual(self.audit()["verified_count"], 1)
        self.record["topology"] = "linear"
        self.write_catalog()
        self.assertEqual(self.audit()["excluded_count"], 1)
        self.record["topology"] = "circular"
        self.reference["provenance"].update(start=1, end=3)
        self.write_catalog()
        self.assertEqual(self.audit()["excluded_count"], 1)

    def test_provenance_failures_are_excluded_not_treated_as_alleles(self):
        self.write_catalog()
        mutations = [
            ("accession_version", "TEST_123"),
            ("accession_version", "TEST_123.2"),
            ("kind", "sharkmer_product"),
            ("schema_version", True),
            ("source_sequence_sha256", "0" * 64),
            ("sequence_sha256", "0" * 64),
            ("source_length", len(self.source) + 1),
            ("start", -1),
            ("start", True),
            ("start", len(self.source)),
            ("end", len(self.source) + 1),
            ("strand", "unknown"),
            ("wraps_origin", 0),
        ]
        for field, value in mutations:
            with self.subTest(field=field, value=value):
                reference = copy.deepcopy(self.reference)
                reference["provenance"][field] = value
                result = self.audit(reference)
                self.assertEqual(result["verified_count"], 0)
                self.assertEqual(result["excluded_count"], 1)
                self.assertTrue(result["excluded"][0]["reason"])

    def test_unverified_source_and_taxon_mismatch_are_explicit(self):
        self.write_catalog()
        reference = copy.deepcopy(self.reference)
        reference.pop("provenance")
        self.assertIn("missing public-record", self.audit(reference)["excluded"][0]["reason"])
        reference = copy.deepcopy(self.reference)
        reference["taxon"] = "Different sample taxon"
        self.assertIn("source organism", self.audit(reference)["excluded"][0]["reason"])

    def test_missing_catalog_never_falls_back_to_embedded_sequence(self):
        result = self.audit()
        self.assertEqual(result["verified_count"], 0)
        self.assertEqual(result["excluded_count"], 1)
        self.assertEqual(result["catalog"]["status"], "unavailable")

    def test_corrupt_source_checksum_excludes_catalog(self):
        self.record["sequence_sha256"] = "0" * 64
        self.write_catalog()
        self.assertEqual(self.audit()["verified_count"], 0)

    def test_duplicate_catalog_keys_are_rejected(self):
        self.write_catalog(b'{"schema_version":1,"records":{},"records":{}}')
        self.assertIn("duplicate catalog key", self.audit()["catalog"]["error"])

    def test_corrupt_gzip_is_reported_as_unavailable(self):
        self.path.write_bytes(bytes.fromhex("1f8b0800000000000203ff"))
        result = self.audit()
        self.assertEqual(result["catalog"]["status"], "unavailable")
        self.assertEqual(result["excluded_count"], 1)

    def test_malformed_topology_is_reported_as_unavailable(self):
        self.record["topology"] = []
        self.write_catalog()
        self.assertIn("source topology", self.audit()["catalog"]["error"])

    def test_known_annotation_conflict_cannot_be_promoted_by_exact_source_match(self):
        self.record["annotation_conflicts"] = {"tArGet": "public annotation conflicts with target label"}
        self.write_catalog()
        self.assertEqual(provenance.verify_reference(self.reference, self.catalog)["sequence"], self.reference["sequence"])
        for gene_name in ("target", "target_1", "target-V9_2", "TARGET", "TaRgeT_1"):
            panel = {"references": [{"gene": gene_name, "sequences": [self.reference]}]}
            result = provenance.audit_references(panel, self.path)
            self.assertEqual(result["verified_count"], 0)
            self.assertIn("source region verified but gene annotation excluded", result["excluded"][0]["reason"])

    def test_malformed_annotation_conflicts_are_reported(self):
        self.record["annotation_conflicts"] = []
        self.write_catalog()
        self.assertIn("annotation conflicts", self.audit()["catalog"]["error"])

    def test_iupac_compatible_substitution_is_not_literal_provenance(self):
        self.reference["sequence"] = self.reference["sequence"].replace("R", "A")
        self.reference["provenance"]["sequence_sha256"] = provenance.sequence_digest(self.reference["sequence"])
        self.write_catalog()
        self.assertEqual(self.audit()["verified_count"], 0)

    def test_source_catalog_symlink_is_not_trusted(self):
        target = Path(self.temporary.name) / "real.json"
        target.write_text(json.dumps(self.catalog))
        self.path.symlink_to(target)
        self.assertEqual(self.audit()["verified_count"], 0)

    def test_oversize_compressed_catalog_is_bounded(self):
        from unittest import mock

        self.write_catalog(b" " * 3000)
        with mock.patch.object(provenance, "MAX_CATALOG_BYTES", 2000):
            self.assertIn("size limit", self.audit()["catalog"]["error"])

    def test_empty_reference_collection_still_records_catalog_receipt(self):
        self.write_catalog()
        result = provenance.audit_references({}, self.path)
        self.assertEqual(result["total_count"], 0)
        self.assertIsNotNone(result["catalog"]["sha256"])

    def test_bootstrap_candidates_are_not_published_as_references(self):
        import yaml

        output = yaml.safe_load(bootstrap_references.format_references_yaml({"target": [self.reference]}))
        self.assertNotIn("references", output)
        self.assertEqual(output["purpose"], "unverified_reference_candidates")
        panel = Path(self.temporary.name) / "panel.yaml"
        panel.write_text("unchanged\n")
        with self.assertRaisesRegex(ValueError, "publication is disabled"):
            bootstrap_references.write_references_to_panel(panel, {})
        self.assertEqual(panel.read_text(), "unchanged\n")

    def test_bootstrap_write_fails_before_input_or_network_work(self):
        result = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts/bootstrap_references.py"),
             str(Path(self.temporary.name) / "missing.yaml"), "--write"],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--write is disabled", result.stderr)


if __name__ == "__main__":
    unittest.main()

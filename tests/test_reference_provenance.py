import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.sharkmer_validate.reference_provenance import extract_region


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "audit_panel_references.py"
SCRIPT_SPEC = importlib.util.spec_from_file_location("audit_panel_references", SCRIPT_PATH)
audit_panel_references = importlib.util.module_from_spec(SCRIPT_SPEC)
SCRIPT_SPEC.loader.exec_module(audit_panel_references)


def write_panel(path, sequence, accession="ABC123"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "references:\n"
        "  - gene: target\n"
        "    sequences:\n"
        f"      - taxon: Test organism\n        accession: {accession}\n        sequence: {sequence}\n"
    )


def write_archive(path, sequence, accession_version="ABC123.1", topology="linear"):
    path.mkdir(parents=True, exist_ok=True)
    metadata = {
        "result": {
            "uids": ["1"],
            "1": {
                "accessionversion": accession_version,
                "slen": len(sequence),
                "organism": "Test organism",
                "taxid": 1,
                "topology": topology,
                "title": "Test source",
            },
        },
    }
    (path / "metadata.json").write_text(json.dumps(metadata))
    (path / "records.fasta").write_text(f">{accession_version} test\n{sequence}\n")


def write_verified_panel(path, provenance):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "references:\n"
        "  - gene: target\n"
        "    sequences:\n"
        "      - taxon: Test organism\n"
        "        accession: ABC123.1\n"
        "        sequence: CGTR\n"
        "        provenance:\n"
        "          schema_version: 1\n"
        "          kind: public_record_region\n"
        "          accession_version: ABC123.1\n"
        f"          source_sequence_sha256: {provenance['source_sequence_sha256']}\n"
        f"          source_length: {provenance['source_length']}\n"
        f"          start: {provenance['start']}\n"
        f"          end: {provenance['end']}\n"
        f"          strand: '{provenance['strand']}'\n"
        "          wraps_origin: false\n"
        f"          sequence_sha256: {provenance['sequence_sha256']}\n"
    )


class ReferenceProvenanceTests(unittest.TestCase):
    def build_catalog(self, sequence, source_sequence, topology="linear"):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panels_root = root / "panels"
            write_panel(panels_root / "test.yaml", sequence)
            write_panel(panels_root / "examples/reference.yaml", sequence)
            archive = root / "archive"
            write_archive(archive, source_sequence, topology=topology)
            return audit_panel_references.build_catalog(
                panels_root,
                archive,
                1024,
                4096,
                "https://example.test/efetch",
                "2026-09-13T00:00:00+00:00",
            )

    def test_exact_subsequence_records_one_based_coordinates(self):
        catalog = self.build_catalog("CGTR", "AACGTRTT")
        entry = catalog["entries"][0]
        self.assertEqual(entry["disposition"], "exact_subsequence")
        self.assertEqual(entry["provenance"]["strand"], "+")
        self.assertEqual(entry["provenance"]["start"], 2)
        self.assertEqual(entry["provenance"]["end"], 6)
        self.assertEqual(entry["provenance"]["sequence_sha256"], audit_panel_references.sha256_text("CGTR"))

    def test_iupac_is_literal_not_wildcard(self):
        catalog = self.build_catalog("CGTR", "AACGTATT")
        self.assertEqual(catalog["entries"][0]["disposition"], "no_exact_literal_match")

    def test_nonunique_coordinates_preserve_alternatives(self):
        catalog = self.build_catalog("ACGA", "AACGATTTACGA")
        entry = catalog["entries"][0]
        self.assertEqual(entry["disposition"], "exact_subsequence")
        self.assertTrue(entry["nonunique_placement"])
        self.assertEqual(len(entry["alternative_coordinates"]), 2)
        self.assertEqual(entry["provenance"]["start"], 1)

    def test_circular_wrap_is_explicit(self):
        catalog = self.build_catalog("CA", "AAGC", topology="circular")
        entry = catalog["entries"][0]
        self.assertEqual(entry["disposition"], "circular_wrap_subsequence")
        self.assertEqual(entry["provenance"]["start"], 3)
        self.assertEqual(entry["provenance"]["end"], 1)
        self.assertTrue(entry["provenance"]["wraps_origin"])

    def test_shared_extractor_rebuilds_recorded_region(self):
        catalog = self.build_catalog("CGTR", "AACGTRTT")
        entry = catalog["entries"][0]
        record = catalog["records"][entry["provenance"]["accession_version"]]
        self.assertEqual(extract_region(record, entry["provenance"]), "CGTR")

    def test_verify_existing_accepts_complete_catalog(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panels_root = root / "panels"
            write_panel(panels_root / "test.yaml", "CGTR")
            write_panel(panels_root / "examples/reference.yaml", "CGTR")
            archive = root / "archive"
            write_archive(archive, "AACGTRTT")
            catalog = audit_panel_references.build_catalog(
                panels_root,
                archive,
                1024,
                4096,
                "https://example.test/efetch",
                "2026-09-13T00:00:00+00:00",
            )
            catalog_path = root / "catalog.json"
            catalog_path.write_text(json.dumps(catalog))
            provenance = catalog["entries"][0]["provenance"]
            write_verified_panel(panels_root / "test.yaml", provenance)
            write_verified_panel(panels_root / "examples/reference.yaml", provenance)
            verification = audit_panel_references.verify_existing(panels_root, catalog_path)
        self.assertEqual(verification["failure_count"], 0)
        self.assertEqual(verification["verified_references"], 2)

    def test_verify_existing_cli_accepts_complete_catalog(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panels_root = root / "panels"
            write_panel(panels_root / "test.yaml", "CGTR")
            write_panel(panels_root / "examples/reference.yaml", "CGTR")
            archive = root / "archive"
            write_archive(archive, "AACGTRTT")
            catalog = audit_panel_references.build_catalog(
                panels_root,
                archive,
                1024,
                4096,
                "https://example.test/efetch",
                "2026-09-13T00:00:00+00:00",
            )
            provenance = catalog["entries"][0]["provenance"]
            write_verified_panel(panels_root / "test.yaml", provenance)
            write_verified_panel(panels_root / "examples/reference.yaml", provenance)
            catalog_path = root / "catalog.json"
            catalog_path.write_text(json.dumps(catalog))
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--verify-existing",
                    "--panels-root",
                    str(panels_root),
                    "--reference-catalog",
                    str(catalog_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["failure_count"], 0)

    def test_oversize_sources_are_not_read(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            panels_root = root / "panels"
            write_panel(panels_root / "test.yaml", "ACGT")
            write_panel(panels_root / "examples/reference.yaml", "ACGT")
            archive = root / "archive"
            write_archive(archive, "AACGTT")
            catalog = audit_panel_references.build_catalog(
                panels_root,
                archive,
                4,
                4096,
                "https://example.test/efetch",
                "2026-09-13T00:00:00+00:00",
            )
        self.assertEqual(catalog["entries"][0]["disposition"], "skipped_oversize_source")


if __name__ == "__main__":
    unittest.main()

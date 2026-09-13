import sys
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from sharkmer_validate import reference_provenance


def panel_references(panel_name, gene):
    panel = yaml.safe_load((ROOT / "panels" / f"{panel_name}.yaml").read_text())
    return [reference for group in panel.get("references", []) if group["gene"] == gene for reference in group["sequences"]]


class PublicReplacementTests(unittest.TestCase):
    def test_all_active_references_reconstruct_from_pinned_sources(self):
        paths = sorted((ROOT / "panels").glob("*.yaml")) + [ROOT / "panels/examples/reference.yaml"]
        for path in paths:
            with self.subTest(panel=path.name):
                audit = reference_provenance.audit_references(yaml.safe_load(path.read_text()))
                self.assertEqual(audit["excluded"], [])
                self.assertEqual(audit["verified_count"], audit["total_count"])

    def test_public_16s_replacements_do_not_reintroduce_muts(self):
        references = panel_references("cnidaria", "16S")
        self.assertIn("LC467070.1", {reference["accession"] for reference in references})
        self.assertNotIn("XM_047005986.1", {reference["accession"] for reference in references})
        catalog, _receipt = reference_provenance.load_catalog()
        self.assertIn("16s", {gene.casefold() for gene in catalog["records"]["XM_047005986.1"]["annotation_conflicts"]})

    def test_rhopilema_cox1_uses_reverse_origin_spanning_gene(self):
        reference = next(reference for reference in panel_references("cnidaria", "CO1") if reference["accession"] == "NC_035741.1")
        provenance = reference["provenance"]
        self.assertEqual((provenance["start"], provenance["end"], provenance["strand"], provenance["wraps_origin"]), (15807, 1542, "-", True))
        self.assertEqual(len(reference["sequence"]), 1590)

    def test_genomic_h3_replacement_is_shared_not_fabricated_alleles(self):
        references = [reference for reference in panel_references("c_elegans", "H3") if reference["accession"] == "Z98866.1"]
        self.assertEqual(len(references), 1)
        self.assertEqual(len(references[0]["sequence"]), 433)
        self.assertEqual(references[0]["provenance"]["source_length"], 125590)
        self.assertNotEqual(references[0]["accession"], "NC_003281.10")

    def test_full_public_genes_replace_short_benchmark_comparators(self):
        gryllus = next(reference for reference in panel_references("insecta", "CO1_1") if reference["accession"] == "PP230540.1")
        heliconius = next(reference for reference in panel_references("insecta", "ND1") if reference["accession"] == "NC_024741.1")
        self.assertEqual(len(gryllus["sequence"]), 1531)
        self.assertEqual(len(heliconius["sequence"]), 939)

    def test_unresolved_gene_conflicts_remain_unfilled(self):
        self.assertEqual(panel_references("angiospermae", "trnV-atpE"), [])
        for gene in ("ITS_1", "ITS_2"):
            self.assertNotIn("Agalma elegans", {reference["taxon"] for reference in panel_references("cnidaria", gene)})


if __name__ == "__main__":
    unittest.main()

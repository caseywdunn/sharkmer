import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SPEC = importlib.util.spec_from_file_location("analyze_sweep", Path(__file__).with_name("analyze_sweep.py"))
ANALYZER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYZER)


LOSS_GENES = ("12S", "16S_2", "CO1_1", "ND1", "ITS_2", "12S", "12S")


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def product(gene, value):
    return {
        "sha256": digest(value),
        "length": 100,
        "reference_match": {"status": "confirmed_product", "on_target": True, "matched_gene": gene},
    }


def result(invocation, product_values, kmer_count):
    grouped = {}
    for gene, value in product_values:
        grouped.setdefault(gene, []).append(product(gene, value))
    return {
        "signature": {"invocation": invocation},
        "status": "complete",
        "classification_status": "complete",
        "execution": {"returncode": 0, "timed_out": False, "wall_time_s": 3.0, "gnu_time": {"peak_rss_bytes": 1000}},
        "metrics": {
            "n_reads_read": 4,
            "n_bases_read": 400,
            "n_subreads_ingested": 4,
            "n_bases_ingested": 400,
            "n_kmers": kmer_count,
        },
        "threshold_diagnostics": [{"threshold": 2, "outcome": "valid"}],
        "genes": [{"gene": gene, "products": values} for gene, values in grouped.items()],
    }


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def protocol(targets):
    return {
        "baseline_role": "baseline",
        "expected_commits": {"baseline": "released", "candidate": "candidate"},
        "samples": [{"panel": "insecta", "input": "sample", "depth": 100, "pairs": 3}],
        "target_metadata": {"insecta": [{"gene": gene, "scope": "high_copy_candidate"} for gene in targets]},
    }


def write_execution(root, frozen_protocol, kmer_length, records):
    schedule = []
    for pair_index, role, values in records:
        invocation = {
            "invocation_id": f"k{kmer_length}_{pair_index}_{role}",
            "cell": "insecta/sample/100",
            "panel": "insecta",
            "pair_index": pair_index,
            "version": role,
            "kmer_length": kmer_length,
        }
        schedule.append(invocation)
        output = result(invocation, values, 1000 - kmer_length)
        stats = {
            "kmer_length": kmer_length,
            "pcr_results": [
                {
                    "gene_name": f"insecta_{gene_result['gene']}",
                    **({"threshold_diagnostics": [{"threshold": 2}]} if role == "candidate" else {}),
                }
                for gene_result in output["genes"]
            ],
        }
        stats_path = root / "attempts" / invocation["invocation_id"] / "output" / f"{invocation['invocation_id']}.stats.yaml"
        write_json(stats_path, stats)
        stats_digest = hashlib.sha256(stats_path.read_bytes()).hexdigest()
        output["stats_path"] = str(stats_path)
        output["stats_sha256"] = stats_digest
        output["raw_output_files"] = [{"path": stats_path.name, "sha256": stats_digest}]
        output["signature"] = {
            "invocation": invocation,
            "source_commit": frozen_protocol["expected_commits"][role],
            "binary_sha256": f"binary-{role}",
            "input_sha256": "input",
            "input_subset": {"sha256": "prefix"},
            "panel_sha256": "panel",
            "settings": {"k": kmer_length},
        }
        write_json(root / "results" / f"{invocation['invocation_id']}.json", output)
    write_json(root / "provenance.json", {"protocol": frozen_protocol, "schedule": schedule})
    write_json(root / "receipts" / "protocol.json", frozen_protocol)


class SweepTest(unittest.TestCase):
    def test_changed_same_gene_is_visible_but_not_exact_retained(self):
        anchor, changed = product("ITS_2", "anchor"), product("ITS_2", "changed")
        anchor["gene"] = "ITS_2"
        changed["gene"] = "ITS_2"
        comparison = ANALYZER.comparison([anchor], [changed])
        self.assertEqual(comparison["exact_retained"], [])
        self.assertEqual(len(comparison["lost"]), 1)
        self.assertEqual(len(comparison["changed_endpoint_or_sequence_unresolved"]), 1)

    def test_nested_sweep_requires_exact_every_replicate_restoration(self):
        loss_values = [(gene, f"loss-{index}") for index, gene in enumerate(LOSS_GENES)]
        targets = {gene for gene, _ in loss_values} | {"retained"}
        reference_protocol = protocol(targets)
        its_digest = digest("loss-4")
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            fixed_root = temporary / "fixed"
            write_execution(
                fixed_root,
                reference_protocol,
                19,
                [
                    record
                    for pair_index in (1, 2, 3)
                    for record in (
                        (pair_index, "baseline", loss_values + [("retained", "stable")]),
                        (pair_index, "candidate", [("retained", "stable")]),
                    )
                ],
            )
            diagnostic_path = fixed_root / "analysis" / "diagnostic-audit.json"
            write_json(diagnostic_path, {"status": "frozen"})
            discovery_protocol = {
                **protocol(targets),
                "allowed_k": [19, 23],
                "fixed_reference": {
                    "execution": str(fixed_root),
                    "analysis": str(fixed_root / "analysis"),
                    "protocol_sha256": hashlib.sha256((fixed_root / "receipts" / "protocol.json").read_bytes()).hexdigest(),
                    "diagnostic_audit_sha256": hashlib.sha256(diagnostic_path.read_bytes()).hexdigest(),
                    "baseline_k": 19,
                    "expected_lost_high_copy_products": 7,
                    "expected_candidate_high_copy_products_all_ten_samples": 1,
                    "exact_its2_sha256": its_digest,
                },
            }
            root = temporary / "discovery"
            for kmer_length in (19, 23):
                records = []
                for pair_index in (1, 2, 3):
                    candidate = [("retained", "stable")]
                    if kmer_length == 23:
                        candidate.append(("ITS_2", "loss-4"))
                    records.extend([(pair_index, "baseline", [("retained", "stable")]), (pair_index, "candidate", candidate)])
                write_execution(root / f"k{kmer_length}", discovery_protocol, kmer_length, records)
            analysis = ANALYZER.analyze(root)
        option = next(option for option in analysis["confirmation_options_ranked"] if option["kmer_length"] == 23)
        self.assertTrue(option["eligible"])
        self.assertTrue(option["exact_its2_ak281180_restored_every_replicate"])
        self.assertEqual(option["exact_restored_known_losses_every_replicate"], [("insecta/sample/100", "ITS_2", its_digest)])
        candidate_counts = analysis["cells"][1]["pairs"][0]["candidate"]["counts"]
        baseline_counts = analysis["cells"][0]["pairs"][0]["baseline"]["counts"]
        self.assertNotEqual(candidate_counts["n_kmers"], baseline_counts["n_kmers"])

    def test_selection_uses_intersection_not_union_across_replicates(self):
        restored = ("ITS_2", digest("its"))
        cells = [{
            "cell": "insecta/sample/100",
            "repeat_sequence_and_classification_stable": {"baseline": True, "candidate": True},
            "pairs": [
                {
                    "same_k_count_parity": True,
                    "candidate_exact_restored_known_losses": [restored],
                    "candidate_vs_k19_candidate": {"lost": []},
                    "candidate_vs_k19_released": {"exact_retained": []},
                },
                {
                    "same_k_count_parity": True,
                    "candidate_exact_restored_known_losses": [],
                    "candidate_vs_k19_candidate": {"lost": []},
                    "candidate_vs_k19_released": {"exact_retained": []},
                },
            ],
        }]
        option = ANALYZER.option_summary(23, cells, restored, True, True)
        self.assertEqual(option["exact_restored_known_losses_every_replicate"], [])
        self.assertFalse(option["eligible"])

    def test_three_sample_selection_allows_gryllus_only_rescue_but_blocks_instability_or_confirmed_loss(self):
        restored = ("ITS_2", digest("gryl-its"))

        def cell(name, restorations, confirmed_loss=False):
            pairs = []
            for recovery in restorations:
                lost = [
                    {
                        "gene": "CO1_1",
                        "sha256": digest("confirmed-anchor"),
                        "reference_match": {"status": "confirmed_product"},
                    }
                ] if confirmed_loss else []
                pairs.append({
                    "same_k_count_parity": True,
                    "candidate_exact_restored_known_losses": recovery,
                    "candidate_vs_k19_candidate": {"lost": lost},
                    "candidate_vs_k19_released": {"exact_retained": []},
                })
            return {
                "cell": name,
                "repeat_sequence_and_classification_stable": {"baseline": True, "candidate": True},
                "pairs": pairs,
            }

        stable = [
            cell("insecta/drosophila/100", [[], [], []]),
            cell("insecta/heliconius/100", [[], [], []]),
            cell("insecta/gryllus/100", [[restored], [restored], [restored]]),
        ]
        option = ANALYZER.option_summary(23, stable, restored, True, True)
        self.assertTrue(option["eligible"])
        self.assertEqual(option["exact_restored_known_losses_every_replicate"], [("insecta/gryllus/100", *restored)])
        self.assertFalse(ANALYZER.option_summary(23, stable, restored, True, False)["eligible"])
        unstable = stable[:-1] + [cell("insecta/gryllus/100", [[restored], [restored], []])]
        self.assertFalse(ANALYZER.option_summary(23, unstable, restored, True, True)["eligible"])
        blocked = stable[:-1] + [cell("insecta/gryllus/100", [[restored], [restored], [restored]], confirmed_loss=True)]
        self.assertFalse(ANALYZER.option_summary(23, blocked, restored, True, True)["eligible"])


if __name__ == "__main__":
    unittest.main()

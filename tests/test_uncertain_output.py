#!/usr/bin/env python3
"""CI copy of independently developed parser/publication controls; synthetic data only.

The later copy in tests/test_uncertain_output.py is a CI regression copy of
these controls, not additional independent evidence.
"""

import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from sharkmer_validate import runner, uncertain


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def assessment(sequence="ACGTT", empty=False):
    cap = {"paths": 32, "sequence_bases": 524288, "marker_occurrences": 4096}
    description = (
        "assembly_status=uncertain uncertainty=repeat_structure "
        "repeat_length=unresolved read_support=not_evaluated "
        f"hypothesis_length={len(sequence)} sequence_sha256={digest(sequence)} "
        "causes=pre_prune_cyclic_scc_node thresholds=2 source_path_records=1 "
        "enumeration=bounded limit_scope=gene search_limits_reached=false retention_truncated=false"
    )
    product = {
        "id": "sample_gene.uncertain_0", "description": description,
        "oriented_sequence": sequence, "sequence_sha256": digest(sequence),
        "sequence_length": len(sequence), "sources": [{
            "threshold": 2, "node_visit_limit": 2, "max_observed_node_visits": 1,
            "omitted_self_loop_marker_occurrences": 0,
            "cyclic_scc_marker_occurrences": 1,
            "retained_collision_marker_occurrences": 0,
            "marker_occurrences": [{
                "cause": "pre_prune_cyclic_scc_node",
                "identity_kind": "oriented_node_sequence",
                "identity_sequence": sequence[:3], "identity_sha256": digest(sequence[:3]),
                "start": 0, "end": 3,
            }], "marker_covered_runs": [{"start": 0, "end": 3}],
        }],
    }
    value = {
        "schema_version": 1, "n_products": int(not empty),
        "product_lengths": [] if empty else [len(sequence)],
        "purpose": "uncertain_repeat_structure_graph_path_hypotheses",
        "coordinate_system": "zero_based_half_open_in_emitted_oriented_sequence",
        "marker_span_scope": "marker_covered_positions_not_complete_ambiguity_or_bridge_ready_intervals",
        "read_support_evaluation": "not_evaluated", "repeat_length_evaluation": "unresolved",
        "deduplication_policy": "exact_oriented_sequence_only_reverse_complements_retained_as_separate_hypotheses",
        "standard_equivalence_policy": "exact_sequence_or_reverse_complement_of_final_standard_output",
        "limits": {"threshold": cap, "gene": cap, "run": cap,
                   "allocated_run_share_for_gene": cap, "unused_share_is_not_redistributed": True},
        "path_search_evaluated": not empty, "search_limits_reached": False,
        "search_limit_scope": "all_attempted_thresholds_for_gene",
        "retention_truncated": False, "thresholds": [],
        "retained_path_records": int(not empty),
        "excluded_standard_equivalent_path_records": 0,
        "exact_duplicate_path_records": 0, "unique_products": int(not empty),
        "products": [] if empty else [product],
    }
    if not empty:
        value["output_file"] = "sample_gene.uncertain.fasta"
    return value


def transaction(directory, standard=True, empty=False):
    uncertain = assessment(empty=empty)
    entry = {"gene_name": "gene", "status": "success" if standard else "fail",
             "n_products": int(standard), "uncertain_products": uncertain}
    files = {}
    if standard:
        entry.update(output_file="sample_gene.fasta", product_lengths=[4])
        files[entry["output_file"]] = ">standard product=0\nACGT\n"
    else:
        entry["failure_reason"] = "repeat structure unresolved"
    if not empty:
        item = uncertain["products"][0]
        files[uncertain["output_file"]] = f'>{item["id"]} {item["description"]}\n{item["oriented_sequence"]}\n'
    stats = {"sharkmer_version": "3.2.0-dev", "sample": "sample", "kmer_length": 19,
             "command": "synthetic fixture", "run_id": "independent-control",
             "run_status": "complete", "output_manifest": "sample.manifest.yaml",
             "pcr_results": [entry]}
    for name, text in files.items():
        (directory / name).write_text(text)
    publish(directory, stats, list(files))
    return stats


def publish(directory, stats, files):
    (directory / "sample.stats.yaml").write_text(yaml.safe_dump(stats))
    manifest = {"schema_version": 1, "producer": "sharkmer", "sample": "sample",
                "run_id": stats["run_id"], "status": "complete", "files": [
                    {"path": name, "sha256": runner._sha256_file(directory / name)}
                    for name in ["sample.stats.yaml", *files]
                ]}
    (directory / "sample.manifest.yaml").write_text(yaml.safe_dump(manifest))
    return manifest


class UncertainIsolationTests(unittest.TestCase):
    def test_manifest_standard_and_uncertain_are_separate(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory)
            runner._validate_stats_manifest(stats, "sample", 19, {"gene"}, "end-to-end")
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual([item["gene"] for item in parsed], ["gene"])
            self.assertEqual(parsed[0]["sequences"], ["ACGT"])
            self.assertEqual(parsed[0]["n_products"], 1)

    def test_uncertain_only_is_zero_standard_products(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory, standard=False)
            runner._validate_stats_manifest(stats, "sample", 19, {"gene"}, "end-to-end")
            self.assertEqual(runner.parse_fasta_products("sample", directory), [])

    def test_enabled_empty_has_no_uncertain_fasta_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            transaction(directory, standard=False, empty=True)
            self.assertEqual(runner.parse_fasta_products("sample", directory), [])

    def test_uncertain_output_cannot_be_explicitly_requested_as_standard(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            transaction(directory)
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory, ["sample_gene.uncertain.fasta"])

    def test_changed_uncertain_sequence_invalidates_entire_transaction(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            transaction(directory, standard=False)
            with (directory / "sample_gene.uncertain.fasta").open("a") as output:
                output.write("A\n")
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_missing_uncertain_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory)
            publish(directory, stats, ["sample_gene.fasta"])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_extra_uncertain_receipt_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory, empty=True)
            (directory / "sample_gene.uncertain.fasta").write_text(">unknown assembly_status=uncertain\nACGTT\n")
            publish(directory, stats, ["sample_gene.fasta", "sample_gene.uncertain.fasta"])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_wrong_gene_uncertain_path_rejected_even_with_valid_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory)
            source = directory / "sample_gene.uncertain.fasta"
            source.rename(directory / "sample_other.uncertain.fasta")
            stats["pcr_results"][0]["uncertain_products"]["output_file"] = "sample_other.uncertain.fasta"
            publish(directory, stats, ["sample_gene.fasta", "sample_other.uncertain.fasta"])
            with self.assertRaises(ValueError):
                runner.parse_fasta_products("sample", directory)

    def test_renamed_or_mixed_uncertain_headers_are_never_promoted(self):
        standard = ">standard product=0\nACGT\n"
        uncertain = ">candidate assembly_status=uncertain uncertainty=repeat_structure\nACGTT\n"
        for contents in (uncertain, uncertain + standard, standard + uncertain):
            with self.subTest(contents=contents), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                (directory / "sample_gene.fasta").write_text(contents)
                try:
                    parsed = runner.parse_fasta_products("sample", directory)
                except ValueError:
                    continue  # Rejecting mixed standard/uncertain output is fail-closed.
                sequences = [sequence for gene in parsed for sequence in gene["sequences"]]
                self.assertNotIn("ACGTT", sequences)

    def test_manifest_does_not_override_uncertain_header_guard(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stats = transaction(directory, empty=True)
            (directory / "sample_gene.fasta").write_text(">mislabelled product=0 assembly_status=uncertain\nACGT\n")
            publish(directory, stats, ["sample_gene.fasta"])
            try:
                parsed = runner.parse_fasta_products("sample", directory)
            except ValueError:
                return
            self.assertFalse(any(gene["n_products"] for gene in parsed))

    def test_legacy_standard_gene_with_uncertain_suffix_remains_standard(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "sample_gene.uncertain.fasta").write_text(">standard product=0\nACGT\n")
            parsed = runner.parse_fasta_products("sample", directory)
            self.assertEqual([gene["gene"] for gene in parsed], ["gene.uncertain"])
            self.assertEqual(parsed[0]["sequences"], ["ACGT"])

    def test_quoted_and_duplicate_uncertainty_labels_cannot_be_promoted(self):
        labels = (
            'assembly_status="uncertain"',
            "'assembly_status=uncertain'",
            'assembly_status=uncertain assembly_status=standard',
            'assembly_status="uncertain" assembly_status=standard',
        )
        for label in labels:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                (directory / "sample_gene.fasta").write_text(f">candidate {label}\nACGTT\n")
                try:
                    parsed = runner.parse_fasta_products("sample", directory)
                except ValueError:
                    continue
                self.assertFalse(any(gene["n_products"] for gene in parsed))

    def test_valid_hash_receipt_does_not_replace_fasta_metadata_binding(self):
        for mutation in ("sequence", "description", "extra_record", "missing_record", "id"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = transaction(directory)
                path = directory / "sample_gene.uncertain.fasta"
                original = path.read_text()
                if mutation == "sequence":
                    modified = original.replace("\nACGTT\n", "\nACGTA\n")
                elif mutation == "description":
                    modified = original.replace("read_support=not_evaluated", "read_support=confirmed")
                elif mutation == "extra_record":
                    modified = original + original
                elif mutation == "missing_record":
                    modified = ""
                else:
                    modified = original.replace("sample_gene.uncertain_0", "other_id")
                path.write_text(modified)
                publish(directory, stats, ["sample_gene.fasta", path.name])
                with self.assertRaises(ValueError):
                    runner.parse_fasta_products("sample", directory)

    def test_invalid_scalar_metadata_is_rejected(self):
        mutations = (
            ("schema_version", True), ("schema_version", 2),
            ("n_products", True), ("n_products", -1), ("unique_products", 2),
            ("retained_path_records", 2), ("exact_duplicate_path_records", 1),
            ("excluded_standard_equivalent_path_records", True),
            ("product_lengths", [True]), ("product_lengths", [6]),
            ("read_support_evaluation", "confirmed"),
            ("repeat_length_evaluation", "resolved"),
            ("search_limits_reached", "false"),
            ("coordinate_system", "one_based"), ("products", None),
        )
        for key, value in mutations:
            with self.subTest(key=key, value=value):
                entry = {"gene_name": "gene", "uncertain_products": assessment()}
                entry["uncertain_products"][key] = value
                with self.assertRaises(ValueError):
                    uncertain.validate_metadata(entry, "sample")

    def test_invalid_source_local_metadata_is_rejected(self):
        mutations = (
            ("product", "sequence_sha256", "0" * 64),
            ("product", "sequence_length", True),
            ("product", "oriented_sequence", "ACGTN"),
            ("product", "sources", []),
            ("source", "threshold", True), ("source", "threshold", 0),
            ("source", "cyclic_scc_marker_occurrences", 0),
            ("source", "marker_occurrences", []),
            ("source", "node_visit_limit", 0),
            ("source", "node_visit_limit", True),
            ("source", "node_visit_limit", "2"),
            ("source", "max_observed_node_visits", 0),
            ("source", "max_observed_node_visits", 3),
            ("source", "max_observed_node_visits", True),
            ("source", "marker_covered_runs", None),
            ("source", "marker_covered_runs", []),
            ("source", "marker_covered_runs", [{"start": False, "end": 3}]),
            ("source", "marker_covered_runs", [{"start": 0, "end": 2}]),
            ("source", "marker_covered_runs", [{"start": 0, "end": 4}]),
            ("source", "marker_covered_runs", [{"start": 0, "end": 3}, {"start": 4, "end": 5}]),
            ("marker", "cause", "verified_repeat_copy_number"),
            ("marker", "identity_kind", "oriented_edge_sequence"),
            ("marker", "identity_kind", None),
            ("marker", "start", True), ("marker", "start", 3),
            ("marker", "end", 6), ("marker", "identity_sequence", "AAA"),
            ("marker", "identity_sha256", "0" * 64),
        )
        for level, key, value in mutations:
            with self.subTest(level=level, key=key, value=value):
                entry = {"gene_name": "gene", "uncertain_products": assessment()}
                product = entry["uncertain_products"]["products"][0]
                target = product
                if level in ("source", "marker"):
                    target = product["sources"][0]
                if level == "marker":
                    target = target["marker_occurrences"][0]
                target[key] = value
                with self.assertRaises(ValueError):
                    uncertain.validate_metadata(entry, "sample")

    def test_marker_union_covers_overlap_adjacency_disjoint_and_duplicates(self):
        cases = (
            ([(0, 2), (2, 4), (4, 5)], [(0, 5)]),
            ([(1, 4), (0, 2), (2, 3)], [(0, 4)]),
            ([(3, 5), (0, 1)], [(0, 1), (3, 5)]),
            ([(0, 3), (0, 3)], [(0, 3)]),
        )
        for spans, runs in cases:
            with self.subTest(spans=spans):
                entry = {"gene_name": "gene", "uncertain_products": assessment()}
                product = entry["uncertain_products"]["products"][0]
                source = product["sources"][0]
                sequence = product["oriented_sequence"]
                source["marker_occurrences"] = [{
                    "cause": "pre_prune_cyclic_scc_node", "identity_kind": "oriented_node_sequence",
                    "start": start, "end": end, "identity_sequence": sequence[start:end],
                    "identity_sha256": digest(sequence[start:end]),
                } for start, end in spans]
                source["cyclic_scc_marker_occurrences"] = len(spans)
                source["marker_covered_runs"] = [{"start": start, "end": end} for start, end in runs]
                self.assertEqual(uncertain.validate_metadata(entry, "sample"), "sample_gene.uncertain.fasta")

    def test_each_cause_has_its_expected_identity_kind(self):
        cases = (
            ("pre_prune_omitted_self_loop_node", "omitted_self_loop_marker_occurrences", "oriented_node_sequence"),
            ("pre_prune_cyclic_scc_node", "cyclic_scc_marker_occurrences", "oriented_node_sequence"),
            ("retained_collision_edge", "retained_collision_marker_occurrences", "oriented_edge_sequence"),
        )
        for cause, count_field, identity_kind in cases:
            with self.subTest(cause=cause):
                entry = {"gene_name": "gene", "uncertain_products": assessment()}
                product = entry["uncertain_products"]["products"][0]
                source = product["sources"][0]
                source["cyclic_scc_marker_occurrences"] = 0
                source[count_field] = 1
                source["marker_occurrences"][0].update(cause=cause, identity_kind=identity_kind)
                product["description"] = product["description"].replace("causes=pre_prune_cyclic_scc_node", f"causes={cause}")
                self.assertEqual(uncertain.validate_metadata(entry, "sample"), "sample_gene.uncertain.fasta")
                source["marker_occurrences"][0]["identity_kind"] = (
                    "oriented_node_sequence" if identity_kind == "oriented_edge_sequence" else "oriented_edge_sequence"
                )
                with self.assertRaises(ValueError):
                    uncertain.validate_metadata(entry, "sample")

    def test_header_cannot_overclaim_or_erase_uncertainty(self):
        mutations = (
            ("assembly_status=uncertain", "assembly_status=standard"),
            ("read_support=not_evaluated", "read_support=confirmed"),
            ("repeat_length=unresolved", "repeat_length=resolved"),
            ("thresholds=2", "thresholds=4"),
            ("source_path_records=1", "source_path_records=3"),
            ("causes=pre_prune_cyclic_scc_node", "causes=retained_collision_edge"),
            ("retention_truncated=false", "retention_truncated=false retention_truncated=true"),
        )
        for old, replacement in mutations:
            with self.subTest(replacement=replacement):
                entry = {"gene_name": "gene", "uncertain_products": assessment()}
                product = entry["uncertain_products"]["products"][0]
                product["description"] = product["description"].replace(old, replacement)
                with self.assertRaises(ValueError):
                    uncertain.validate_metadata(entry, "sample")


if __name__ == "__main__":
    unittest.main(verbosity=2)

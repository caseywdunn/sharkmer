"""Synthetic controls for experimental exact graph histories, never real assays."""

import copy
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch


SOURCE = Path(__file__).resolve().parents[1] / "scripts" / "audit_read_histories.py"
SPEC = importlib.util.spec_from_file_location("audit_read_histories", SOURCE)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def dna(length, seed):
    generator = random.Random(seed)
    return bytes(generator.choice(b"ACGT") for _ in range(length))


LEFT, CYCLE, RIGHT = dna(21, 17), dna(21, 23), dna(21, 41)
SHORT = LEFT + CYCLE * 2 + RIGHT
LONG = LEFT + CYCLE * 3 + RIGHT


def graph_document(sequences, k=19, omitted=()):
    edges = {sequence[start:start + k] for sequence in sequences for start in range(len(sequence) - k + 1)} - set(omitted)
    nodes = {node for edge in edges for node in (edge[:-1], edge[1:])}
    outgoing = {node: [] for node in nodes}
    incoming = {node: [] for node in nodes}
    for edge in edges:
        outgoing[edge[:-1]].append(edge[1:])
        incoming[edge[1:]].append(edge[:-1])
    labels = {node: None for node in nodes}
    for ordinal, component in enumerate(audit.cyclic_components(nodes, outgoing, incoming)):
        labels.update({node: ordinal for node in component})
    return {"schema_version": 1, "target": {"kmer_length": k}, "graph": {
        "phase": "synthetic", "node_count": len(nodes), "edge_count": len(edges),
        "nodes": [{"node_id": node.decode(), "oriented_sequence": node.decode(),
                   "incoming_degree": len(incoming[node]), "outgoing_degree": len(outgoing[node]),
                   "cyclic_scc_id": labels[node]} for node in sorted(nodes)],
        "edges": [{"edge_id": edge.decode(), "oriented_sequence": edge.decode(),
                   "source_node_id": edge[:-1].decode(), "target_node_id": edge[1:].decode(),
                   "source_cyclic_scc_id": labels[edge[:-1]], "target_cyclic_scc_id": labels[edge[1:]],
                   "count": 10} for edge in sorted(edges)]}}


def manifest(sequences, regions=None):
    regions = ["region"] * len(sequences) if regions is None else regions
    by_region = {}
    for ordinal, (sequence, region) in enumerate(zip(sequences, regions)):
        by_region.setdefault(region, []).append({"id": "class_" + str(ordinal), "sequence": sequence.decode(),
                                               "sequence_sha256": audit.digest(sequence), "kind": "synthetic"})
    return {"schema_version": 1, "regions": [{"id": region, "classes": values} for region, values in by_region.items()]}


def record(sequence, graph, classes, ordinal=1, read_id=None, quality=None):
    quality = b"I" * len(sequence) if quality is None else quality
    histories = audit.read_histories(sequence, quality, graph)
    return {"ordinal": ordinal, "read_id": str(ordinal) if read_id is None else read_id,
            "canonical_sequence_sha256": audit.digest(min(sequence, audit.reverse_complement(sequence))),
            "histories": histories, "observations": audit.observations_for_read(sequence, quality, classes, histories)}


def marginal_counts(runs):
    """Current threading counts an edge/local pair once per read, not per visit."""
    edges, links = Counter(), Counter()
    for run in runs:
        edges.update(set(run["edges"]))
        links.update({(event["incoming_edge"], event["outgoing_edge"]) for event in run["branch_events"]})
    return edges, links


class Histories(unittest.TestCase):
    def setUp(self):
        self.graph = audit.prepare_graph(graph_document([SHORT, LONG]))
        self.classes = audit.prepare_manifest(manifest([SHORT, LONG]), self.graph)

    def records(self, sequence, reverse=False):
        template = b"AC" + sequence + b"GT"
        reads = [template[offset:offset + len(sequence) + 2] for offset in range(3)]
        return [record(audit.reverse_complement(read) if reverse else read, self.graph, self.classes, index + 1) for index, read in enumerate(reads)]

    def test_three_positioned_ids_and_sequence_groups_pass(self):
        summary = audit.summarize_classes(self.classes, self.records(LONG))[1]
        self.assertTrue(summary["experimental_local_eligibility"])
        self.assertEqual(len(summary["witness"]), 3)
        self.assertEqual(summary["distinct_footprints"], 3)
        self.assertEqual(len(summary["canonical_sequence_groups"]), 3)

    def test_reverse_complement_normalizes_footprints(self):
        forward = audit.summarize_classes(self.classes, self.records(LONG))[1]
        reverse = audit.summarize_classes(self.classes, self.records(LONG, True))[1]
        self.assertEqual(forward["canonical_sequence_groups"], reverse["canonical_sequence_groups"])
        self.assertEqual([x["footprint"] for x in forward["witness"]], [x["footprint"] for x in reverse["witness"]])
        self.assertTrue(reverse["experimental_local_eligibility"])

    def test_more_than_two_cycle_visits_preserve_occurrences(self):
        longer = LEFT + CYCLE * 4 + RIGHT
        pairs = [audit.read_histories(sequence, b"I" * len(sequence), self.graph)["+"]["q20"]["runs"][0] for sequence in (LONG, longer)]
        self.assertEqual(set(pairs[0]["edges"]), set(pairs[1]["edges"]))
        local = lambda run: {(event["incoming_edge"], event["outgoing_edge"]) for event in run["branch_events"]}
        self.assertEqual(local(pairs[0]), local(pairs[1]))
        self.assertEqual(marginal_counts([pairs[0]]), marginal_counts([pairs[1]]))
        self.assertNotEqual(pairs[0]["edges"], pairs[1]["edges"])
        self.assertEqual(pairs[1]["sequence"], longer.decode())
        self.assertGreater(pairs[1]["repeated_edge_occurrences"], pairs[0]["repeated_edge_occurrences"])
        self.assertGreater(len(pairs[1]["cyclic_node_visits"]), len(pairs[0]["cyclic_node_visits"]))
        for visit in pairs[1]["cyclic_node_visits"]:
            start = visit["read_node_start"]
            self.assertEqual(visit["node_id"], longer[start:start + 18].decode())

    def test_unequal_repeat_mixture_preserves_both_classes(self):
        reads = self.records(SHORT) + [record(SHORT, self.graph, self.classes, 4)]
        reads += [dict(item, ordinal=item["ordinal"] + 10, read_id="long" + item["read_id"]) for item in self.records(LONG)]
        summaries = audit.summarize_classes(self.classes, reads)
        self.assertTrue(all(item["experimental_local_eligibility"] for item in summaries))
        self.assertEqual([item["qualifying_records"] for item in summaries], [4, 3])

    def test_invalid_base_and_low_quality_split_runs(self):
        for ambiguity in (b"N", b"R", b"X"):
            sequence = SHORT[:40] + ambiguity + SHORT[41:]
            history = audit.read_histories(sequence, b"I" * len(sequence), self.graph)["+"]["q20"]
            self.assertGreaterEqual(len(history["runs"]), 2)
            self.assertTrue(any("invalid_base" in value["reasons"] for value in history["breaks"]))
            self.assertFalse(any(run["start"] <= 40 < run["end"] for run in history["runs"]))
        quality = b"I" * 40 + b"!" + b"I" * (len(SHORT) - 41)
        history = record(SHORT, self.graph, self.classes, quality=quality)
        self.assertFalse(history["observations"][0]["qualifying_local_witness"])
        self.assertEqual(len(history["histories"]["+"]["unfiltered"]["runs"]), 1)
        self.assertGreaterEqual(len(history["histories"]["+"]["q20"]["runs"]), 2)

    def test_missing_graph_kmer_splits_without_stitching(self):
        edge = SHORT[10:29]
        graph = audit.prepare_graph(graph_document([SHORT, LONG], omitted=[edge]))
        classes = audit.prepare_manifest(manifest([SHORT]), graph)
        history = record(SHORT, graph, classes)
        self.assertEqual(len(history["histories"]["+"]["q20"]["runs"]), 2)
        self.assertFalse(history["observations"][0]["qualifying_local_witness"])
        self.assertFalse(classes[0]["representations"][0]["represented"])

    def test_substitution_and_indel_do_not_repair_path(self):
        substitution = b"A" if SHORT[40:41] != b"A" else b"C"
        for changed in (SHORT[:40] + substitution + SHORT[41:], SHORT[:40] + b"A" + SHORT[40:], SHORT[:40] + SHORT[41:]):
            observed = record(changed, self.graph, self.classes)
            self.assertFalse(observed["observations"])
            self.assertFalse(any(run["sequence"] == SHORT.decode() for run in observed["histories"]["+"]["q20"]["runs"]))

    def test_low_quality_endpoint_never_clipped(self):
        observed = record(SHORT, self.graph, self.classes, quality=b"I" * (len(SHORT) - 1) + b",")
        self.assertEqual(observed["observations"][0]["physical_interval"], [0, len(SHORT)])
        self.assertFalse(observed["observations"][0]["qualifying_local_witness"])
        self.assertEqual(observed["observations"][0]["footprint"], [0, len(SHORT)])

    def test_duplicate_rc_and_quality_variants_do_not_multiply(self):
        reads = [record(LONG, self.graph, self.classes, 1),
                 record(LONG, self.graph, self.classes, 2, quality=b"H" * len(LONG)),
                 record(audit.reverse_complement(LONG), self.graph, self.classes, 3)]
        summary = audit.summarize_classes(self.classes, reads)[1]
        self.assertEqual(summary["qualifying_records"], 3)
        self.assertEqual(summary["distinct_footprints"], 1)
        self.assertEqual(len(summary["canonical_sequence_groups"]), 1)
        self.assertFalse(summary["experimental_local_eligibility"])

    def test_iupac_outside_class_preserves_rc_sequence_group(self):
        sequence = b"RKMBD" + LONG + b"HVSYN"
        reverse = audit.reverse_complement(sequence)
        self.assertEqual(reverse[:5], b"NRSBD")
        forward_record = record(sequence, self.graph, self.classes)
        reverse_record = record(reverse, self.graph, self.classes, 2)
        self.assertEqual(forward_record["canonical_sequence_sha256"], reverse_record["canonical_sequence_sha256"])
        self.assertEqual(forward_record["observations"][0]["footprint"], reverse_record["observations"][0]["footprint"])
        self.assertTrue(forward_record["observations"][0]["qualifying_local_witness"])
        self.assertTrue(reverse_record["observations"][0]["qualifying_local_witness"])

    def test_reused_ids_and_uncoupled_sets_fail(self):
        reads = self.records(LONG)
        for item in reads:
            item["read_id"] = "same_id"
        self.assertFalse(audit.summarize_classes(self.classes, reads)[1]["experimental_local_eligibility"])
        links = {"a": {(0, 150)}, "b": {(0, 150)}, "c": {(1, 151), (2, 152)}}
        self.assertEqual(len(audit.corroboration_witness(links)), 2)

    def test_low_quality_competing_placement_remains_ambiguous(self):
        sequence = SHORT + b"N" + SHORT
        quality = b"I" * len(SHORT) + b"I" + b"!" * len(SHORT)
        observed = record(sequence, self.graph, self.classes, quality=quality)
        self.assertEqual(len(observed["observations"]), 2)
        self.assertFalse(any(value["unique_unfiltered_placement"] for value in observed["observations"]))
        self.assertFalse(any(value["qualifying_local_witness"] for value in observed["observations"]))

    def test_class_tie_and_palindrome_remain_ambiguous(self):
        classes = audit.prepare_manifest(manifest([SHORT, SHORT]), self.graph)
        self.assertFalse(any(item["qualifying_local_witness"] for item in record(SHORT, self.graph, classes)["observations"]))
        half = dna(25, 77)
        palindrome = half + audit.reverse_complement(half)
        graph = audit.prepare_graph(graph_document([palindrome]))
        classes = audit.prepare_manifest(manifest([palindrome]), graph)
        observed = record(palindrome, graph, classes)
        self.assertEqual(len(observed["observations"]), 2)
        self.assertFalse(any(item["qualifying_local_witness"] for item in observed["observations"]))

    def test_two_region_phase_histories_separate_equal_marginals(self):
        left, middle, right = dna(25, 82), dna(60, 90), dna(25, 99)
        templates = [left + a + middle + b + right for a, b in ((b"A", b"A"), (b"C", b"C"), (b"A", b"C"), (b"C", b"A"))]
        graph = audit.prepare_graph(graph_document(templates))
        runs = [audit.read_histories(sequence, b"I" * len(sequence), graph)["+"]["q20"]["runs"][0] for sequence in templates]
        edges = lambda pool: set().union(*(set(run["edges"]) for run in pool))
        links = lambda pool: {(event["incoming_edge"], event["outgoing_edge"]) for run in pool for event in run["branch_events"]}
        self.assertEqual(edges(runs[:2]), edges(runs[2:]))
        self.assertEqual(links(runs[:2]), links(runs[2:]))
        self.assertEqual(marginal_counts(runs[:2]), marginal_counts(runs[2:]))
        self.assertNotEqual({run["sequence"] for run in runs[:2]}, {run["sequence"] for run in runs[2:]})

    def test_omitted_self_loop_and_partial_repeat_have_no_witness(self):
        sequence = LEFT + b"A" * 45 + RIGHT
        graph = audit.prepare_graph(graph_document([sequence], omitted=[b"A" * 19]))
        classes = audit.prepare_manifest(manifest([sequence]), graph)
        self.assertFalse(record(sequence, graph, classes)["observations"][0]["qualifying_local_witness"])
        self.assertFalse(record(LONG[10:80], self.graph, self.classes)["observations"])

    def test_shuffled_records_preserve_class_summary(self):
        reads = self.records(LONG)
        self.assertEqual(audit.summarize_classes(self.classes, reads), audit.summarize_classes(self.classes, list(reversed(reads))))

    def test_graph_labels_degrees_and_endpoints_fail_closed(self):
        for mutation in ("scc", "degree", "endpoint", "duplicate"):
            graph = graph_document([SHORT])
            if mutation == "scc":
                graph["graph"]["nodes"][0]["cyclic_scc_id"] = 999
            elif mutation == "degree":
                graph["graph"]["nodes"][0]["outgoing_degree"] += 1
            elif mutation == "endpoint":
                graph["graph"]["edges"][0]["target_node_id"] = "A" * 18
            else:
                graph["graph"]["edges"].append(copy.deepcopy(graph["graph"]["edges"][0]))
            with self.assertRaises(ValueError):
                audit.prepare_graph(graph)

    def test_manifest_checksum_and_bounds_fail_closed(self):
        document = manifest([SHORT])
        document["regions"][0]["classes"][0]["sequence_sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            audit.prepare_manifest(document, self.graph)
        with patch.object(audit, "MAX_CLASSES", 1):
            with self.assertRaises(ValueError):
                audit.prepare_manifest(manifest([SHORT, LONG]), self.graph)

    def write_fastq(self, directory, sequences):
        raw = b"".join(b"@read" + str(index).encode() + b"\n" + sequence + b"\n+\n" + b"I" * len(sequence) + b"\n" for index, sequence in enumerate(sequences))
        path = Path(directory) / "reads.fastq"
        path.write_bytes(raw)
        expected = {"records": len(sequences), "bases": sum(map(len, sequences)), "size_bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}
        return path, expected

    def test_streaming_prefix_and_ledger_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.write_fastq(directory, [SHORT, LONG])
            report = audit.scan_fastq(path, expected, self.graph, self.classes)
            self.assertEqual(report["consumed_prefix"], expected)
            self.assertEqual(report["counts"]["retained_records"], 2)
            self.assertEqual(report["staged_fastq_before"], report["staged_fastq_after"])
            self.assertEqual("".join(record["raw_fastq"] for record in report["records"]), path.read_text())

    def test_completely_unrepresented_class_is_still_reported(self):
        sequence = dna(70, 123)
        classes = audit.prepare_manifest(manifest([sequence]), self.graph)
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.write_fastq(directory, [sequence])
            report = audit.scan_fastq(path, expected, self.graph, classes)
        self.assertEqual(report["classes"][0]["exact_records"], 1)
        self.assertFalse(report["classes"][0]["experimental_local_eligibility"])
        self.assertFalse(report["records"][0]["histories"]["+"]["q20"]["runs"])

    def test_caps_and_identity_fail_without_partial_success(self):
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.write_fastq(directory, [SHORT, LONG])
            for name, value in (("MAX_HISTORY_RECORDS", 1), ("MAX_TOTAL_EDGE_VISITS", 1), ("MAX_READ_LENGTH", 30), ("MAX_LINE_BYTES", 20)):
                with patch.object(audit, name, value), self.assertRaises(ValueError):
                    audit.scan_fastq(path, expected, self.graph, self.classes)
            with self.assertRaises(ValueError):
                audit.scan_fastq(path, {**expected, "sha256": "0" * 64}, self.graph, self.classes)
            expected["records"] += 1
            with self.assertRaises(ValueError):
                audit.scan_fastq(path, expected, self.graph, self.classes)

    def test_input_change_during_scan_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.write_fastq(directory, [SHORT])
            before = audit.file_identity(path)
            after = {**before, "ctime_ns": before["ctime_ns"] + 1}
            with patch.object(audit, "file_identity", side_effect=[before, after]), self.assertRaises(ValueError):
                audit.scan_fastq(path, expected, self.graph, self.classes)

    def test_graph_change_during_scan_fails_before_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            reads, expected = self.write_fastq(directory, [SHORT])
            graph_path, classes_path = directory / "graph.json", directory / "classes.json"
            expected_path, output = directory / "prefix.json", directory / "output.json"
            graph_path.write_text(json.dumps(graph_document([SHORT, LONG])))
            classes_path.write_text(json.dumps(manifest([SHORT, LONG])))
            expected_path.write_text(json.dumps(expected))

            def modified_scan(*args):
                graph_path.write_text(graph_path.read_text() + "\n")
                return {}

            arguments = [str(SOURCE), "--graph", str(graph_path), "--manifest", str(classes_path),
                         "--expected-prefix", str(expected_path), "--fastq", str(reads), "--output", str(output)]
            with patch("sys.argv", arguments), patch.object(audit, "scan_fastq", side_effect=modified_scan), self.assertRaises(ValueError):
                audit.main()
            self.assertFalse(output.exists())

    def test_malformed_quality_and_sequence_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, expected = self.write_fastq(directory, [SHORT])
            raw = path.read_bytes()
            for changed in (raw.replace(b"\n+\n", b"\n-\n"), raw[:-2] + b" \n", raw.replace(SHORT, b"!" + SHORT[1:])):
                path.write_bytes(changed)
                with self.assertRaises(ValueError):
                    audit.scan_fastq(path, expected, self.graph, self.classes)

    def test_publish_is_no_clobber_and_byte_bounded(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.json"
            audit.publish_json(output, {"status": "complete"})
            with self.assertRaises(FileExistsError):
                audit.publish_json(output, {"status": "changed"})
            self.assertEqual(json.loads(output.read_text()), {"status": "complete"})
            second = Path(directory) / "second.json"
            with patch.object(audit, "MAX_OUTPUT_BYTES", 1), self.assertRaises(ValueError):
                audit.publish_json(second, {"status": "complete"})
            self.assertFalse(second.exists())
            exact_bytes = len(json.dumps({}, sort_keys=True, separators=(",", ":")).encode())
            with patch.object(audit, "MAX_OUTPUT_BYTES", exact_bytes), self.assertRaises(ValueError):
                audit.publish_json(second, {})
            self.assertFalse(second.exists())


if __name__ == "__main__":
    unittest.main()

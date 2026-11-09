import copy
import hashlib
import importlib.util
import random
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("audit_read_bridges", ROOT / "scripts/audit_read_bridges.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def make_candidate(identifier, sequence):
    return {"id": identifier, "sequence": sequence, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()}


class ReadBridgeTests(unittest.TestCase):
    def setUp(self):
        generator = random.Random(813)
        self.left = "".join(generator.choices("ACGT", k=60))
        self.right = "".join(generator.choices("ACGT", k=60))
        self.true = self.left + "A" * 40 + self.right
        self.short = self.left + "A" * 18 + self.right
        self.manifest = {
            "schema_version": 1,
            "candidates": [make_candidate("true", self.true), make_candidate("short", self.short)],
            "events": [
                {"id": "true_repeat", "candidate_id": "true", "start": 60, "end": 100, "comparison_group": "repeat"},
                {"id": "short_repeat", "candidate_id": "short", "start": 60, "end": 78, "comparison_group": "repeat"},
            ],
        }

    def scan(self, reads, manifest=None, max_substitutions=0, modify_expected=None):
        contents = b"".join(f"@{identifier}\n{sequence}\n+\n{quality}\n".encode() for identifier, sequence, quality in reads)
        expected = {"records": len(reads), "bases": sum(len(sequence) for _identifier, sequence, _quality in reads), "size_bytes": len(contents), "sha256": hashlib.sha256(contents).hexdigest()}
        if modify_expected:
            expected.update(modify_expected)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reads.fastq"
            path.write_bytes(contents)
            return audit.scan_fastq(path, expected, manifest or self.manifest, max_substitutions)

    def reads_from(self, template):
        return [(f"read{start}", template[start:-5], "I" * (len(template) - start - 5)) for start in (5, 10, 15)]

    def test_long_repeat_bridge_does_not_support_collapsed_template(self):
        result = self.scan(self.reads_from(self.true))
        self.assertEqual(result["events"][0]["local_bridge_status"], "corroborated_local_bridge")
        self.assertEqual(result["events"][1]["records"], 0)
        self.assertEqual(result["events"][1]["local_bridge_status"], "unresolved")

    def test_clean_short_repeat_and_unequal_mixture_remain_visible(self):
        reads = self.reads_from(self.short) + self.reads_from(self.true)[:1]
        result = self.scan(reads)
        self.assertEqual(result["events"][1]["local_bridge_status"], "corroborated_local_bridge")
        self.assertEqual(result["events"][0]["records"], 1)
        self.assertEqual(result["events"][0]["local_bridge_status"], "unresolved")

    def test_identical_copies_and_reverse_complements_do_not_multiply_footprints(self):
        sequence = self.true[5:-5]
        reverse = audit.reverse_complement(sequence.encode()).decode()
        reads = [(f"copy{index}", sequence if index % 2 else reverse, "I" * len(sequence)) for index in range(12)]
        event = self.scan(reads)["events"][0]
        self.assertEqual(event["records"], 12)
        self.assertEqual(event["distinct_footprints"], 1)
        self.assertEqual(event["distinct_sequences"], 1)
        self.assertEqual(event["local_bridge_status"], "unresolved")

    def test_reused_identifier_cannot_supply_corroboration(self):
        reads = [("same_id", sequence, quality) for _identifier, sequence, quality in self.reads_from(self.true)]
        event = self.scan(reads)["events"][0]
        self.assertEqual(event["distinct_footprints"], 3)
        self.assertEqual(event["distinct_ids"], 1)
        self.assertEqual(event["local_bridge_status"], "unresolved")

    def test_alternative_or_repeated_template_placements_stay_ambiguous(self):
        manifest = copy.deepcopy(self.manifest)
        manifest["candidates"].append(make_candidate("copy", self.true))
        event = self.scan(self.reads_from(self.true), manifest)["events"][0]
        self.assertEqual(event["records"], 3)
        self.assertEqual(event["unique_records"], 0)
        self.assertEqual(event["local_bridge_status"], "unresolved")
        manifest["candidates"] = [make_candidate("true", self.true + self.true), make_candidate("short", self.short)]
        event = self.scan(self.reads_from(self.true), manifest)["events"][0]
        self.assertFalse(event["unique_template_placement"])

    def test_distinct_ids_and_footprints_require_a_distinct_matching_witness(self):
        original = self.reads_from(self.true)
        reads = [("same_id", sequence, quality) for _identifier, sequence, quality in original]
        reads.extend([(identifier, original[0][1], original[0][2]) for identifier in ("extra1", "extra2")])
        event = self.scan(reads)["events"][0]
        self.assertEqual(event["distinct_footprints"], 3)
        self.assertEqual(event["distinct_ids"], 3)
        self.assertEqual(len(event["corroboration_witness"]), 2)
        self.assertEqual(event["local_bridge_status"], "unresolved")

    def test_complete_flanks_are_required_individually(self):
        for start, end in ((20, 100), (60, len(self.true) - 20)):
            manifest = copy.deepcopy(self.manifest)
            manifest["events"][0].update(start=start, end=end)
            event = self.scan(self.reads_from(self.true), manifest)["events"][0]
            self.assertEqual(event["callability"], "missing_complete_flanks")
            self.assertEqual(event["records"], 0)

    def test_one_substitution_is_separate_and_never_accepts_n_or_an_indel(self):
        sequence = self.true[5:-5]
        mutation = "C" if sequence[75] != "C" else "G"
        changed = sequence[:75] + mutation + sequence[76:]
        reads = [("changed", changed, "I" * len(changed))]
        self.assertEqual(self.scan(reads)["events"][0]["records"], 0)
        event = self.scan(reads, max_substitutions=1)["events"][0]
        self.assertEqual(event["one_substitution_records"], 1)
        self.assertEqual(event["exact_records"], 0)
        for changed in (sequence[:75] + "N" + sequence[76:], sequence[:75] + "A" + sequence[75:]):
            self.assertEqual(self.scan([("altered", changed, "I" * len(changed))], max_substitutions=1)["events"][0]["records"], 0)

    def test_low_quality_inside_bridge_is_not_joined_around(self):
        sequence = self.true[5:-5]
        quality = "I" * 75 + "!" + "I" * (len(sequence) - 76)
        self.assertEqual(self.scan([("poor", sequence, quality)])["events"][0]["records"], 0)

    def test_one_substitution_observed_window_competes_with_unassayed_alternative(self):
        sequence = self.left + self.right
        alternate = list(sequence)
        for offset in (54, 65):
            alternate[offset] = "C" if alternate[offset] != "C" else "G"
        observed = sequence[:54] + alternate[54] + sequence[55:]
        manifest = {"schema_version": 1, "candidates": [make_candidate("focal", sequence), make_candidate("alternate", "".join(alternate))], "events": [{"id": "focal", "candidate_id": "focal", "start": 60, "end": 61}]}
        event = self.scan([("tie", observed, "I" * len(observed))], manifest, 1)["events"][0]
        self.assertEqual(event["records"], 1)
        self.assertEqual(event["unique_records"], 0)
        self.assertEqual(event["ambiguous_records"], 1)

    def test_palindromic_window_and_repeated_read_occurrences_are_not_unique(self):
        half = self.left[:21]
        palindrome = half + audit.reverse_complement(half.encode()).decode()
        manifest = {"schema_version": 1, "candidates": [make_candidate("palindrome", palindrome)], "events": [{"id": "palindrome", "candidate_id": "palindrome", "start": 21, "end": 21}]}
        event = self.scan([("palindrome", palindrome, "I" * len(palindrome))], manifest)["events"][0]
        self.assertEqual(event["records"], 1)
        self.assertEqual(event["unique_records"], 0)
        repeated = self.true + self.true
        event = self.scan([("repeated", repeated, "I" * len(repeated))])["events"][0]
        self.assertEqual(event["records"], 1)
        self.assertEqual(event["unique_records"], 0)

    def test_disconnected_reads_do_not_bridge_a_repeat(self):
        sequence = self.left + "A" * 180 + self.right
        manifest = {"schema_version": 1, "candidates": [make_candidate("long", sequence)], "events": [{"id": "long", "candidate_id": "long", "start": 60, "end": 240}]}
        reads = [("fragment/1", sequence[:100], "I" * 100), ("fragment/2", sequence[-100:], "I" * 100)]
        event = self.scan(reads, manifest)["events"][0]
        self.assertEqual(event["records"], 0)
        self.assertTrue(event["longer_than_every_read"])
        self.assertEqual(event["local_bridge_status"], "unresolved")

    def test_recombinant_with_only_marginal_support_is_unresolved(self):
        sequence = self.left + "ACGT" * 20 + self.right
        first = self.left + "ACGT" * 20 + "C" * 60
        second = "T" * 60 + "ACGT" * 20 + self.right
        manifest = {"schema_version": 1, "candidates": [make_candidate("recombinant", sequence), make_candidate("first", first), make_candidate("second", second)], "events": [{"id": "complete", "candidate_id": "recombinant", "start": 60, "end": 140}]}
        reads = [("left", first[:140], "I" * 140), ("right", second[60:], "I" * 140)]
        event = self.scan(reads, manifest)["events"][0]
        self.assertEqual(event["records"], 0)
        self.assertEqual(event["local_bridge_status"], "unresolved")

    def test_tied_patterns_on_one_read_never_depend_on_event_order(self):
        sequence = self.true + self.short
        reads = [("conflicting", sequence, "I" * len(sequence))]
        forward = self.scan(reads)
        manifest = copy.deepcopy(self.manifest)
        manifest["events"].reverse()
        reverse = self.scan(reads, manifest)
        self.assertTrue(all(event["unique_records"] == 0 for event in forward["events"] + reverse["events"]))

    def test_changed_candidate_or_prefix_and_malformed_fastq_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "Consumed-prefix"):
            self.scan(self.reads_from(self.true), modify_expected={"sha256": "0" * 64})
        manifest = copy.deepcopy(self.manifest)
        manifest["candidates"][0]["sequence_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "Candidate sequence checksum"):
            self.scan(self.reads_from(self.true), manifest)
        with self.assertRaisesRegex(ValueError, "Malformed FASTQ"):
            self.scan([("truncated", self.true, "I")])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "truncated.fastq"
            path.write_bytes(b"@header\nACGT\n+\n")
            with self.assertRaisesRegex(ValueError, "Truncated FASTQ"):
                audit.scan_fastq(path, {"records": 1}, self.manifest)


if __name__ == "__main__":
    unittest.main()

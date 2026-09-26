import importlib.util
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("audit_its_recruitment", ROOT / "scripts/audit_its_recruitment.py")
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def entry(identifier, sequence, **extra):
    return {"id": identifier, "sequence": sequence.decode(), "sequence_sha256": audit.digest(sequence), **extra}


def manifest(sequence, others=(), decoys=()):
    return {"schema_version": 1, "regions": [{"id": "region", "classes": [entry("historical", sequence)] + list(others)}], "decoys": list(decoys)}


def oracle_distance(first, second):
    """Independent unbanded Wagner-Fischer matrix, without production helpers."""
    matrix = [list(range(len(second) + 1))]
    for row, base in enumerate(first, 1):
        matrix.append([row] + [0] * len(second))
        for column, other in enumerate(second, 1):
            matrix[row][column] = min(matrix[row - 1][column] + 1,
                                      matrix[row][column - 1] + 1,
                                      matrix[row - 1][column - 1] + (base != other))
    return matrix[-1][-1]


def fastq_bytes(reads):
    return b"".join(b"@" + identifier + b"\n" + sequence + b"\n+\n" + quality + b"\n"
                    for identifier, sequence, quality in reads)


def expectation(raw, reads):
    return {"records": len(reads), "bases": sum(len(read[1]) for read in reads),
            "size_bytes": len(raw), "sha256": audit.digest(raw)}


class RecruitmentTests(unittest.TestCase):
    def setUp(self):
        randomizer = random.Random(417)
        self.sequence = bytes(randomizer.choice(b"ACGT") for _ in range(61))

    def observe(self, read, quality=None, specification=None):
        regions, decoys = audit.prepare_manifest(specification or manifest(self.sequence))
        return audit.observations_for_read(read, quality or b"I" * len(read), regions, decoys)

    def scan(self, reads, specification=None):
        raw = fastq_bytes(reads)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            path.write_bytes(raw)
            return audit.scan_fastq(path, expectation(raw, reads), specification or manifest(self.sequence))

    def test_distance_and_traceback(self):
        self.assertEqual(audit.edit_alignment(b"ACGT", b"ACGT", 4)["distance"], 0)
        for observed, counts in [(b"ACTT", (1, 0, 0)), (b"ACGGT", (0, 1, 0)), (b"ACT", (0, 0, 1))]:
            alignment = audit.edit_alignment(b"ACGT", observed, 4)
            self.assertEqual(tuple(alignment[key] for key in ("substitutions", "insertions", "deletions")), counts)
        self.assertTrue(audit.edit_alignment(b"AAA", b"AA", 2)["tied_optimal_alignments"])
        self.assertIsNone(audit.edit_alignment(b"AAAA", b"CCCC", 2))

    def test_independent_distance_oracle(self):
        rng = random.Random(109)
        for _ in range(250):
            first = bytes(rng.choice(b"ACGT") for _ in range(rng.randrange(12)))
            second = bytes(rng.choice(b"ACGT") for _ in range(rng.randrange(12)))
            distance = oracle_distance(first, second)
            for maximum in (0, 2, 4):
                expected = distance if distance <= maximum else None
                self.assertEqual(audit.edit_distance(first, second, maximum), expected)
                observed = audit.edit_alignment(first, second, maximum)
                self.assertEqual(None if observed is None else observed["distance"], expected)
                if observed:
                    self.assertEqual(sum(observed[key] for key in ("substitutions", "insertions", "deletions")), distance)

    def test_seed_hits_match_exhaustive_oracle(self):
        rng = random.Random(819)
        pattern = self.sequence[:21]
        cases = [pattern, pattern[1:], b"A" + pattern, pattern[:7] + b"T" + pattern[7:],
                 pattern[:7] + pattern[8:14] + b"G" + pattern[14:]]
        for _ in range(20):
            altered = bytearray(pattern)
            for _ in range(2):
                operation = rng.randrange(3)
                position = rng.randrange(len(altered))
                if operation == 0:
                    altered[position] = rng.choice(b"ACGT")
                elif operation == 1:
                    altered.insert(position, rng.choice(b"ACGT"))
                else:
                    del altered[position]
            cases.append(bytes(altered))
        for case in cases:
            sequence = b"GCT" + case + b"ATC"
            expected = []
            for start in range(len(sequence)):
                for end in range(start + 19, min(len(sequence), start + 23) + 1):
                    distance = oracle_distance(pattern, sequence[start:end])
                    if distance <= 2:
                        expected.append({"start": start, "end": end, "distance": distance})
            self.assertEqual(audit.all_anchor_hits(sequence, pattern), expected)

    def test_four_edit_context_seeds_against_exhaustive_oracle(self):
        pattern = self.sequence[:25]
        cases = [pattern, pattern[:3] + b"AA" + pattern[3:-2], pattern[1:7] + b"G" + pattern[7:] + b"AA"]
        for case in cases:
            sequence = b"TCG" + case + b"GCA"
            expected = []
            for start in range(len(sequence)):
                for end in range(start + 21, min(len(sequence), start + 29) + 1):
                    distance = oracle_distance(pattern, sequence[start:end])
                    if distance <= 4:
                        expected.append({"start": start, "end": end, "distance": distance})
            self.assertEqual(audit.all_anchor_hits(sequence, pattern, 4), expected)

    def test_exact_wins_over_approximate_flank_shifts(self):
        result = self.observe(b"GCTAG" + self.sequence + b"ATCGT")[0]
        self.assertGreater(len(result["pairs"]), 1)
        self.assertEqual(result["status"], "qualified_unique_local_alignment")
        self.assertEqual(result["best_q20"][0]["distance"], 0)
        self.assertEqual(result["best_q20"][0]["read_interval"], [5, 66])

    def test_reverse_orientation(self):
        result = self.observe(audit.reverse_complement(b"GCTAG" + self.sequence + b"ATCGT"))[0]
        self.assertEqual(result["best_q20"][0]["strand"], "-")
        self.assertEqual(result["best_q20"][0]["oriented_interval"], [5, 66])

    def test_multiple_substitutions_and_indels(self):
        altered = bytearray(self.sequence)
        for position in (25, 29, 35):
            altered[position] = ord("A") if altered[position] != ord("A") else ord("C")
        winner = self.observe(bytes(altered))[0]["best_q20"][0]
        self.assertEqual(winner["distance"], 3)
        inserted = self.sequence[:29] + b"AA" + self.sequence[29:]
        self.assertEqual(self.observe(inserted)[0]["best_q20"][0]["insertions"], 2)
        deleted = self.sequence[:29] + self.sequence[31:]
        self.assertEqual(self.observe(deleted)[0]["best_q20"][0]["deletions"], 2)

    def test_low_quality_retained_but_not_qualified(self):
        quality = bytearray(b"I" * len(self.sequence))
        quality[30] = ord("!")
        result = self.observe(self.sequence, bytes(quality))[0]
        self.assertTrue(result["best_unfiltered"])
        self.assertFalse(result["best_q20"])

    def test_quality_boundary_reverse_and_qualified_competitor(self):
        quality = bytearray(b"5" * len(self.sequence))  # Phred33+20.
        quality[30] = ord("4")
        self.assertTrue(self.observe(self.sequence, b"5" * len(self.sequence))[0]["best_q20"])
        self.assertFalse(self.observe(audit.reverse_complement(self.sequence), bytes(quality[::-1]))[0]["best_q20"])
        altered = bytearray(self.sequence)
        altered[30] = ord("A") if altered[30] != ord("A") else ord("C")
        read = self.sequence + b"NNNNN" + bytes(altered)
        quality = b"!" * len(self.sequence) + b"IIIII" + b"I" * len(altered)
        result = self.observe(read, quality)[0]
        self.assertEqual(result["best_unfiltered"][0]["distance"], 0)
        self.assertEqual(result["best_q20"][0]["distance"], 1)
        self.assertEqual(result["best_q20"][0]["read_interval"], [66, 127])

    def test_non_acgt_inner_does_not_discard_flanks(self):
        read = self.sequence[:30] + b"N" + self.sequence[31:]
        result = self.observe(read)[0]
        self.assertTrue(result["anchors"])
        self.assertTrue(result["pairs"])
        self.assertFalse(result["best_q20"])

    def test_invalid_low_quality_and_outside_distance_summaries(self):
        invalid = self.sequence[:30] + b"N" + self.sequence[31:]
        distant = self.sequence[:21] + b"A" * 70 + self.sequence[-21:]
        reads = [(b"invalid", invalid, b"I" * len(invalid)),
                 (b"low", self.sequence, b"!" * len(self.sequence)),
                 (b"distant", distant, b"I" * len(distant))]
        summary = self.scan(reads)["regions"][0]
        self.assertEqual(summary["non_acgt_only_pair_records"], 1)
        self.assertEqual(summary["q20_pair_records"], 1)
        self.assertEqual(summary["low_quality_only_pair_records"], 1)
        self.assertEqual(summary["low_quality_alignment_records"], 1)
        self.assertEqual(summary["low_quality_only_alignment_records"], 1)
        self.assertEqual(summary["unaligned_pair_records"], 1)

    def test_flank_and_span_edit_limits_and_flank_indels(self):
        for observed in (self.sequence[:7] + b"AA" + self.sequence[7:],
                         self.sequence[:7] + self.sequence[9:],
                         self.sequence[:-9] + b"AA" + self.sequence[-9:],
                         self.sequence[:-9] + self.sequence[-7:]):
            result = self.observe(observed)[0]
            self.assertTrue(result["best_q20"])
            self.assertEqual(result["best_q20"][0]["distance"], 2)
        four = bytearray(self.sequence)
        for position in (22, 26, 30, 34):
            four[position] = ord("A") if four[position] != ord("A") else ord("C")
        self.assertEqual(oracle_distance(self.sequence, four), 4)
        self.assertEqual(self.observe(bytes(four))[0]["best_q20"][0]["distance"], 4)
        five = self.sequence[:30] + b"AAAAA" + self.sequence[30:]
        self.assertEqual(oracle_distance(self.sequence, five), 5)
        self.assertFalse(self.observe(five)[0]["best_q20"])
        pattern = self.sequence[:21]
        for count in (2, 3):
            altered = bytearray(pattern)
            for position in (4, 10, 16)[:count]:
                altered[position] = ord("A") if altered[position] != ord("A") else ord("C")
            self.assertEqual(oracle_distance(pattern, altered), count)
            exact_interval = [hit for hit in audit.all_anchor_hits(bytes(altered), pattern)
                              if hit["start"] == 0 and hit["end"] == len(altered)]
            self.assertEqual(bool(exact_interval), count <= 2)

    def test_single_flank_and_nonspanning_distance_retained(self):
        result = self.observe(self.sequence[:21] + b"N" * 40)[0]
        self.assertTrue(any(anchor["side"] == "left" for anchor in result["anchors"]))
        self.assertFalse(result["pairs"])
        read = self.sequence[:21] + b"A" * 70 + self.sequence[-21:]
        result = self.observe(read)[0]
        self.assertTrue(result["pairs"])
        self.assertFalse(result["best_q20"])

    def test_disconnected_support_and_wrong_flank_order(self):
        left, right = self.sequence[:21], self.sequence[-21:]
        reads = [(b"left", left, b"I" * len(left)), (b"right", right, b"I" * len(right)),
                 (b"wrong", right + b"N" * 30 + left, b"I" * 72)]
        result = self.scan(reads)
        summary = result["regions"][0]
        self.assertEqual(summary["left_only_records"], 1)
        self.assertEqual(summary["right_only_records"], 1)
        self.assertEqual(summary["both_flanks_no_ordered_pair_records"], 1)
        self.assertEqual(summary["ordered_pair_records"], 0)
        self.assertEqual(summary["no_ordered_pair_records"], 3)

    def test_class_tie_is_ambiguous(self):
        first, second, read = map(bytearray, (self.sequence, self.sequence, self.sequence))
        first[30], second[30], read[30] = ord("A"), ord("C"), ord("G")
        result = self.observe(bytes(read), specification=manifest(bytes(first), [entry("alternative", bytes(second))]))[0]
        self.assertEqual(result["status"], "qualified_ambiguous_local_alignment")
        self.assertEqual(len(result["best_q20"]), 2)

    def test_repeat_placement_tie_and_pair_cap(self):
        result = self.observe(self.sequence + b"GGGGG" + self.sequence)[0]
        self.assertEqual(result["status"], "qualified_ambiguous_local_alignment")
        self.assertEqual(len(result["best_q20"]), 2)
        self.assertGreater(len(result["pairs"]), 256)
        with patch.object(audit, "MAX_PAIRS_PER_READ", 1):
            with self.assertRaisesRegex(ValueError, "pair cap"):
                self.observe(self.sequence)

    def test_repeat_copy_classes_do_not_collapse(self):
        alternative = self.sequence[:30] + self.sequence[21:42] * 3 + self.sequence[30:]
        specification = manifest(self.sequence, [entry("extra_cycle", alternative)])
        result = self.observe(alternative, specification=specification)[0]
        self.assertEqual(result["best_q20"][0]["class_id"], "extra_cycle")
        self.assertEqual(result["best_q20"][0]["distance"], 0)

    def test_unequal_mixture_and_duplicate_identifiers(self):
        extra = self.sequence[:30] + self.sequence[21:42] + self.sequence[30:]
        reads = [(b"duplicate", self.sequence, b"I" * len(self.sequence))] * 3
        reads += [(b"extra", extra, b"I" * len(extra))]
        result = self.scan(reads, manifest(self.sequence, [entry("extra_cycle", extra)]))
        self.assertEqual([record["ordinal"] for record in result["read_evidence"]], [1, 2, 3, 4])
        summary = result["regions"][0]
        self.assertEqual(summary["distinct_recruited_read_ids"], 2)
        self.assertEqual(summary["classes"]["historical"]["unique_records"], 3)
        self.assertEqual(summary["classes"]["historical"]["distinct_read_ids"], 1)
        self.assertEqual(summary["classes"]["extra_cycle"]["unique_records"], 1)

    def test_palindrome_keeps_both_orientations(self):
        half = self.sequence[:30]
        palindrome = half + audit.reverse_complement(half)
        result = self.observe(palindrome, specification=manifest(palindrome))[0]
        self.assertEqual(result["status"], "qualified_ambiguous_local_alignment")
        self.assertEqual({item["strand"] for item in result["best_q20"]}, {"+", "-"})
        self.assertTrue(all(item["read_interval"] == [0, 60] for item in result["best_q20"]))

    def test_decoy_both_strands_self_and_other_contexts(self):
        decoy = self.sequence + b"GGG" + audit.reverse_complement(self.sequence)
        specification = manifest(self.sequence, decoys=[entry("decoy", decoy)])
        specification["regions"][0]["classes"][0]["self_context"] = [{"decoy_id": "decoy", "strand": "+", "start": 0, "end": 61}]
        result = self.observe(self.sequence, specification=specification)[0]
        hits = result["best_q20"][0]["decoy_placements"]
        self.assertTrue(any(hit["exact_self_interval"] for hit in hits))
        self.assertTrue(any(hit["strand"] == "-" and hit["distance"] == 0 and not hit["overlaps_self_context"] for hit in hits))
        self.assertTrue(any(hit["overlaps_self_context"] and not hit["exact_self_interval"] for hit in hits))

    def test_shared_cross_locus_string_remains_descriptive(self):
        specification = manifest(self.sequence, decoys=[entry("other_locus", b"GGA" + self.sequence + b"AAA")])
        specification["regions"].append({"id": "second_locus", "classes": [entry("shared", self.sequence)]})
        results = self.observe(self.sequence, specification=specification)
        self.assertEqual(len(results), 2)
        for result in results:
            self.assertEqual(result["status"], "qualified_unique_local_alignment")
            exact = [hit for hit in result["best_q20"][0]["decoy_placements"] if hit["distance"] == 0]
            self.assertEqual(len(exact), 1)
            self.assertFalse(exact[0]["overlaps_self_context"])

    def test_manifest_validation(self):
        specification = manifest(self.sequence)
        specification["regions"][0]["classes"].append(entry("same", self.sequence))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            audit.prepare_manifest(specification)
        specification = manifest(self.sequence)
        specification["regions"][0]["classes"][0]["sequence_sha256"] = "bad"
        with self.assertRaisesRegex(ValueError, "hash"):
            audit.prepare_manifest(specification)

    def test_all_declared_evidence_bounds_fail_closed(self):
        specification = manifest(self.sequence, decoys=[entry("decoy", self.sequence)])
        specification["regions"][0]["classes"][0]["self_context"] = [
            {"decoy_id": "decoy", "strand": "+", "start": 0, "end": 61}]
        manifest_limits = [("MAX_MANIFEST_BYTES", 10), ("MAX_REGIONS", 0),
                           ("MAX_CLASSES_PER_REGION", 0), ("MAX_TOTAL_CLASSES", 0),
                           ("MAX_DECOYS", 0), ("MAX_DECOY_LENGTH", 1),
                           ("MAX_TOTAL_DECOY_BASES", 1), ("MAX_SELF_CONTEXTS_PER_CLASS", 0),
                           ("MAX_IDENTIFIER_BYTES", 1)]
        for constant, maximum in manifest_limits:
            with self.subTest(constant=constant), patch.object(audit, constant, maximum):
                with self.assertRaisesRegex(ValueError, "cap"):
                    audit.prepare_manifest(specification)
        observation_limits = [("MAX_READ_LENGTH", 1), ("MAX_FLANK_PLACEMENTS_PER_READ", 0),
                              ("MAX_PAIRS_PER_READ", 0), ("MAX_DECOY_PLACEMENTS", 0),
                              ("MAX_DECOY_PLACEMENTS_PER_READ", 0), ("MAX_RECORD_EVIDENCE_BYTES", 1)]
        for constant, maximum in observation_limits:
            with self.subTest(constant=constant), patch.object(audit, constant, maximum):
                with self.assertRaisesRegex(ValueError, "cap"):
                    self.observe(self.sequence, specification=specification)
        reads = [(b"read", self.sequence, b"I" * len(self.sequence))]
        for constant in ("MAX_OUTPUT_BYTES", "MAX_FASTQ_LINE_BYTES"):
            with self.subTest(constant=constant), patch.object(audit, constant, 1):
                with self.assertRaisesRegex(ValueError, "cap"):
                    self.scan(reads)

    def test_malformed_fastq_and_input_mutation(self):
        reads = [(b"one", self.sequence, b"I" * len(self.sequence))]
        raw = fastq_bytes(reads)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            for malformed in (raw[:-5], raw.replace(b"@one", b"bad!"), raw.replace(b"\n+\n", b"\n+other\n"),
                              raw.replace(b"I", b" "), raw[:-1], raw.replace(b"@one", b"@\xff")):
                path.write_bytes(malformed)
                with self.subTest(raw=malformed[:10]), self.assertRaises(ValueError):
                    audit.scan_fastq(path, expectation(raw, reads), manifest(self.sequence))
            path.write_bytes(raw)
            real_observations = audit.observations_for_read
            def mutate(*arguments):
                observed = real_observations(*arguments)
                path.write_bytes(raw + b"trailing-mutation")
                return observed
            with patch.object(audit, "observations_for_read", side_effect=mutate):
                with self.assertRaisesRegex(ValueError, "changed"):
                    audit.scan_fastq(path, expectation(raw, reads), manifest(self.sequence))

    def test_staged_identity_includes_unconsumed_suffix(self):
        reads = [(b"one", self.sequence, b"I" * len(self.sequence))]
        raw = fastq_bytes(reads)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            full = raw + b"ignored suffix outside the declared prefix\n"
            path.write_bytes(full)
            result = audit.scan_fastq(path, expectation(raw, reads), manifest(self.sequence))
            self.assertEqual(result["input_subset"]["sha256"], audit.digest(raw))
            self.assertEqual(result["staged_file_before"]["sha256"], audit.digest(full))
            self.assertEqual(result["staged_file_before"], result["staged_file_after"])

    def test_cli_provenance_no_clobber_and_failure_receipt(self):
        reads = [(b"one", self.sequence, b"I" * len(self.sequence))]
        raw = fastq_bytes(reads)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq, hypotheses, prefix, output = [root / name for name in ("reads.fastq", "manifest.json", "prefix.json", "result.json")]
            fastq.write_bytes(raw)
            hypotheses.write_bytes(audit.encoded_json(manifest(self.sequence)))
            prefix.write_bytes(audit.encoded_json(expectation(raw, reads)))
            arguments = [str(ROOT / "scripts/audit_its_recruitment.py"), "--fastq", str(fastq), "--manifest", str(hypotheses),
                         "--expected-prefix", str(prefix), "--output", str(output)]
            with patch.object(sys, "argv", arguments):
                audit.main()
                original = output.read_bytes()
                result = json.loads(original)
                self.assertEqual(result["provenance"]["argv"], arguments)
                self.assertEqual(result["provenance"]["scanner_sha256"], audit.digest(Path(audit.__file__).read_bytes()))
                self.assertEqual(result["provenance"]["manifest_sha256"], audit.digest(hypotheses.read_bytes()))
                with self.assertRaisesRegex(ValueError, "overwrite"):
                    audit.main()
                self.assertEqual(output.read_bytes(), original)
            failed = root / "failed.json"
            prefix.write_bytes(audit.encoded_json(dict(expectation(raw, reads), sha256="0" * 64)))
            with patch.object(sys, "argv", arguments[:-1] + [str(failed)]):
                with self.assertRaisesRegex(ValueError, "identity"):
                    audit.main()
                self.assertFalse(failed.exists())
                receipt = failed.with_name(failed.name + ".failure.json")
                self.assertEqual(json.loads(receipt.read_bytes())["status"], "failed")
                original_failure = receipt.read_bytes()
                with self.assertRaisesRegex(ValueError, "overwrite"):
                    audit.main()
                self.assertEqual(receipt.read_bytes(), original_failure)

    def test_source_change_blocks_publication(self):
        reads = [(b"one", self.sequence, b"I" * len(self.sequence))]
        raw = fastq_bytes(reads)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq, hypotheses, prefix, output = [root / name for name in ("reads.fastq", "manifest.json", "prefix.json", "result.json")]
            fastq.write_bytes(raw)
            hypotheses.write_bytes(audit.encoded_json(manifest(self.sequence)))
            prefix.write_bytes(audit.encoded_json(expectation(raw, reads)))
            arguments = [audit.__file__, "--fastq", str(fastq), "--manifest", str(hypotheses),
                         "--expected-prefix", str(prefix), "--output", str(output)]
            real_read = Path.read_bytes
            source_calls = []
            def altered_source(path):
                result = real_read(path)
                if path == Path(audit.__file__):
                    source_calls.append(path)
                    if len(source_calls) > 1:
                        result += b"\n# simulated mutation\n"
                return result
            with patch.object(sys, "argv", arguments), patch.object(Path, "read_bytes", altered_source):
                with self.assertRaisesRegex(ValueError, "source changed"):
                    audit.main()
            self.assertFalse(output.exists())
            self.assertTrue(output.with_name(output.name + ".failure.json").exists())

    def test_atomic_publication_and_output_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            path.write_bytes(b"original")
            with self.assertRaises(FileExistsError):
                audit.publish_json(path, {"replacement": True})
            self.assertEqual(path.read_bytes(), b"original")
            absent = Path(directory) / "oversize.json"
            with self.assertRaisesRegex(ValueError, "cap"):
                audit.publish_json(absent, {"payload": "too large"}, maximum=1)
            self.assertFalse(absent.exists())
            self.assertEqual(sorted(item.name for item in Path(directory).iterdir()), ["result.json"])

    def test_prefix_identity_counts_and_recruited_cap(self):
        raw = b"@one\n" + self.sequence + b"\n+\n" + b"I" * len(self.sequence) + b"\n"
        expected = {"records": 1, "bases": len(self.sequence), "size_bytes": len(raw), "sha256": audit.digest(raw)}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq"
            path.write_bytes(raw)
            result = audit.scan_fastq(path, expected, manifest(self.sequence))
            self.assertEqual(result["regions"][0]["unique_records"], 1)
            self.assertEqual(result["regions"][0]["classes"]["historical"]["distinct_read_ids"], 1)
            with self.assertRaisesRegex(ValueError, "identity"):
                audit.scan_fastq(path, dict(expected, sha256="wrong"), manifest(self.sequence))
            with patch.object(audit, "MAX_RECRUITED_RECORDS", 0):
                with self.assertRaisesRegex(ValueError, "Recruited-record cap"):
                    audit.scan_fastq(path, expected, manifest(self.sequence))


if __name__ == "__main__":
    unittest.main()

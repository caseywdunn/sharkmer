import unittest
import tempfile
from pathlib import Path

from triage_saved_gains import (
    analyze, bind_inputs, canonical_kmers, differing_window, parse_fasta,
    reverse_complement, sha256_bytes, write_results,
)


class SavedGainControls(unittest.TestCase):
    def test_reverse_complement_sets_and_multiplicity(self):
        sequence = "AACCGTACGTTAGCCATGCA"
        self.assertEqual(canonical_kmers(sequence, 5),
                         canonical_kmers(reverse_complement(sequence), 5))
        self.assertEqual(canonical_kmers("AAAAA", 3), canonical_kmers("AAAAAA", 3))

    def test_windows_reconstruct_and_do_not_overlap(self):
        for first, second in [("AACGT", "AATGT"), ("AAACT", "AACT"),
                              ("AAAA", "AAAAA"), ("ACGT", "ACGT"),
                              ("ACGT", "AGCT"), ("A", "T")]:
            with self.subTest(first=first, second=second):
                row = differing_window(first, second)
                start, end = row["query_start"], row["query_end"]
                other_start, other_end = row["neighbor_start"], row["neighbor_end"]
                self.assertLessEqual(start, end)
                self.assertLessEqual(other_start, other_end)
                self.assertEqual(first[:start] + row["neighbor_allele"] + first[end:], second)
                self.assertEqual(second[:other_start] + row["query_allele"] + second[other_end:], first)
        insertion = differing_window("AAACT", "AACT")
        self.assertEqual((insertion["query_allele"], insertion["neighbor_allele"]), ("A", ""))

    def test_long_window_is_not_edit_distance(self):
        row = differing_window("AT" + "C" * 100 + "GA", "AG" + "C" * 100 + "TA")
        self.assertEqual(row["max_allele_length"], 102)

    def test_binding_and_hash_mismatch(self):
        row = {"accession": "sample", "gene": "target", "length": 4,
               "sha256": sha256_bytes(b"ACGT")}
        mapping = row | {"query_id": "q000"}
        self.assertEqual(bind_inputs([row], [mapping], {"q000": "ACGT"}, 1)[0]["sequence"], "ACGT")
        with self.assertRaisesRegex(ValueError, "length/hash"):
            bind_inputs([row], [mapping], {"q000": "ACGA"}, 1)
        with self.assertRaisesRegex(ValueError, "one-to-one"):
            bind_inputs([row], [mapping], {"q001": "ACGT"}, 1)

    def test_duplicate_fasta_id_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            parse_fasta(">q000\nACGT\n>q000\nACTG\n")

    def test_existing_results_are_not_overwritten_or_partially_created(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "existing.tsv").write_bytes(b"saved evidence")
            with self.assertRaises(FileExistsError):
                write_results(output, {"new.tsv": b"new", "existing.tsv": b"replacement"})
            self.assertFalse((output / "new.tsv").exists())
            self.assertEqual((output / "existing.tsv").read_bytes(), b"saved evidence")

    def test_nearest_stays_within_target_and_preserves_ties(self):
        records = []
        for index, (gene, sequence) in enumerate([
                ("a", "ACGT"), ("a", "ACGA"), ("a", "ACGC"), ("b", "ACGT")]):
            records.append({"query_id": str(index), "accession": "s", "gene": gene,
                            "sha256": sha256_bytes(sequence.encode()),
                            "length": len(sequence), "sequence": sequence})
        _, windows = analyze(records)
        self.assertEqual(len(windows), 3)
        self.assertEqual(windows[0]["equal_score_neighbors_before_hash_tiebreak"], 2)
        self.assertNotEqual(windows[0]["neighbor_query_id"], "3")


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression controls for inference mistakes, using no biological truth labels."""
import unittest
import tempfile
from pathlib import Path

import analyze_linkage as audit


def read(sequence, ordinal=1, quality=None):
    return {"sequence":sequence,"quality":quality or "I"*len(sequence),
            "ordinal":ordinal,"header":f"@synthetic.{ordinal}"}


class LinkageControls(unittest.TestCase):
    def setUp(self):
        self.span = audit.make_span("TTCGCC", "TCTGT", "T")

    def test_separate_local_reads_do_not_link(self):
        reads = [read(self.span[:60]),read(self.span[-60:],2)]
        self.assertTrue(audit.exact_matches(reads,self.span[:45]))
        self.assertTrue(audit.exact_matches(reads,self.span[-43:]))
        self.assertEqual(audit.exact_matches(reads,self.span),[])
        self.assertEqual(audit.anchored_spans(reads),[])

    def test_opposite_strand_links_with_correct_original_coordinates(self):
        synthetic = read("GG"+audit.rc(self.span)+"AAA")
        matches = audit.exact_matches([synthetic],self.span)
        self.assertEqual(len(matches),1)
        self.assertEqual(matches[0]["strand"],"-")
        self.assertEqual(matches[0]["original_read_span"],[2,2+len(self.span)])

    def test_low_quality_between_high_quality_windows_prevents_link(self):
        quality = list("I"*len(self.span));quality[65]="!"
        synthetic=read(self.span,quality="".join(quality))
        self.assertTrue(audit.exact_matches([synthetic],self.span[:45]))
        self.assertTrue(audit.exact_matches([synthetic],self.span[-43:]))
        self.assertEqual(audit.exact_matches([synthetic],self.span),[])
        self.assertFalse(audit.anchored_spans([synthetic])[0]["passes_q20_entire_span"])

    def test_intervening_mismatch_is_retained_and_not_linked(self):
        alternative=audit.make_span("TTCGCC","ACGGG","T")
        observed=read(alternative)
        self.assertEqual(audit.exact_matches([observed],self.span),[])
        self.assertEqual(audit.anchored_spans([observed])[0]["sequence"],alternative)
        windows=[{"sequence":self.span[:45],"event_id":"a"}]
        placed=audit.candidate_placements(observed,{"sequence":self.span},windows)
        self.assertTrue(placed)
        self.assertFalse(placed[0]["exact_overlap"])
        self.assertTrue(placed[0]["mismatches"])
        self.assertLess(placed[0]["contiguous_q20_exact_intervals"][0]["candidate_span"][1],len(self.span))

    def test_competing_offsets_are_preserved(self):
        synthetic=read("GACCTAGT")
        candidate={"sequence":"AAAGACCTAGTCCCCGACCTAGTTTT"}
        placed=audit.candidate_placements(synthetic,candidate,[{"sequence":"GACCTAGT","event_id":"x"}])
        self.assertEqual([p["candidate_start_offset"] for p in placed if p["strand"]=="+"],[3,15])

    def test_repeated_anchor_combinations_are_preserved(self):
        synthetic=read("ACGTTTACGTTGCA")
        spans=audit.anchored_spans([synthetic],"ACGT","TGCA")
        self.assertEqual(len([s for s in spans if s["strand"]=="+"]),2)

    def test_novel_intervening_block_is_not_forced_into_candidate(self):
        novel=audit.make_span("CTCGCCC","TCGGG","T")
        self.assertEqual(audit.parse_blocks(novel)["C"],"TCGGG")
        self.assertEqual(audit.anchored_spans([read(novel)])[0]["sequence"],novel)

    def test_indel_literal_strings_do_not_get_forced_equal(self):
        longer=audit.make_span("TTCGCCC","TCTGT","T")
        self.assertEqual(audit.exact_matches([read(longer)],self.span),[])
        self.assertEqual(audit.parse_blocks(audit.anchored_spans([read(longer)])[0]["sequence"])["A"],"TTCGCCC")

    def test_quality_on_reverse_strand_is_reversed_together_with_bases(self):
        sequence="AAA"+audit.rc(self.span)+"GG"
        quality=list("I"*len(sequence));quality[3+65]="!"
        self.assertEqual(audit.exact_matches([read(sequence,quality="".join(quality))],self.span),[])

    def test_variable_core_with_partial_flanks_is_separate_from_full_span(self):
        synthetic=read(self.span[:-18])
        matches=audit.core_matches([synthetic],"TTCGCC","TCTGT","T")
        self.assertEqual(len(matches),1)
        self.assertEqual(matches[0]["consecutive_exact_q20_left_flank_bases"],21)
        self.assertEqual(matches[0]["consecutive_exact_q20_right_flank_bases"],3)
        self.assertEqual(audit.exact_matches([synthetic],self.span),[])

    def test_saved_outputs_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            directory=Path(temp)
            audit.check_output_paths(directory)
            (directory/"report.json").write_text("saved witness")
            with self.assertRaises(FileExistsError):
                audit.check_output_paths(directory)
            self.assertEqual((directory/"report.json").read_text(),"saved witness")


if __name__=="__main__":
    unittest.main()

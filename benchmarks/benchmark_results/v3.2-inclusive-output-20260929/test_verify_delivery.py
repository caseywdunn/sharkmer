"""Mutation controls for the output-regression comparison's allowed changes."""
import copy
import subprocess
import unittest

import verify_delivery as delivery


class ComparisonControls(unittest.TestCase):
    def test_documented_awk_filters_preserve_wrapped_sequences(self):
        text = '>s assembly_status=standard\nACG\nT\n>u assembly_status=uncertain\nTT\nAAA\n'
        for label, expected in [('standard', [{'header': 's assembly_status=standard', 'sequence': 'ACGT'}]),
                                ('uncertain', [{'header': 'u assembly_status=uncertain', 'sequence': 'TTAAA'}])]:
            output = subprocess.check_output(['awk', f'/^>/ {{keep = / assembly_status={label}( |$)/}} keep'],
                                             input=text, text=True)
            self.assertEqual(delivery.fasta_text(output), expected)

    def test_class_annotations_cannot_hide_missing_duplicate_or_unknown_status(self):
        for header in (
            'seq product=0',
            'seq assembly_status=unknown',
            'seq assembly_status=standard assembly_status=uncertain',
            'seq assembly_status=standard assembly_status=standard',
        ):
            with self.subTest(header=header), self.assertRaises(ValueError):
                delivery.class_of(header)
        original = 'sample_gene_0 product=0 length=12 score=5.00'
        self.assertEqual(delivery.old_standard_header(original + ' assembly_status=standard'), original)
        changed = original.replace('score=5.00', 'score=4.00') + ' assembly_status=standard'
        self.assertNotEqual(delivery.old_standard_header(changed), original)

    def test_only_uncertainty_id_and_delivery_path_changes_are_ignored(self):
        original = {'output_file': 'sample_gene.uncertain.fasta', 'retained_path_records': 1,
                    'products': [{'id': 'sample_gene.uncertain_0', 'sequence_sha256': 'abc',
                                  'oriented_sequence': 'ACGT', 'sources': [{'threshold': 2,
                                      'marker_occurrences': [{'start': 0, 'end': 2}]}]}]}
        delivery_change = copy.deepcopy(original)
        delivery_change['output_file'] = 'sample_gene.fasta'
        delivery_change['products'][0]['id'] = 'sample_gene_0.uncertain'
        expected = delivery.assessment_identity(original)
        self.assertEqual(delivery.assessment_identity(delivery_change), expected)
        for mutation in ('sequence', 'threshold', 'marker', 'retention'):
            changed = copy.deepcopy(delivery_change)
            if mutation == 'sequence':
                changed['products'][0]['oriented_sequence'] = 'ACGA'
            elif mutation == 'threshold':
                changed['products'][0]['sources'][0]['threshold'] = 3
            elif mutation == 'marker':
                changed['products'][0]['sources'][0]['marker_occurrences'][0]['start'] = 1
            else:
                changed['retained_path_records'] = 2
            with self.subTest(mutation=mutation):
                self.assertNotEqual(delivery.assessment_identity(changed), expected)


if __name__ == '__main__':
    unittest.main()

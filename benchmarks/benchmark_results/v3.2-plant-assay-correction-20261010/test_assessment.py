"""Controls for silent status, multiplicity and synthetic-truth overclaims."""
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('plant_comparison', HERE / 'run_comparison.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def job(variant, records, scope='real', status='pass'):
    entries = []
    for gene in ('control_gene', runner.CHANGED_GENE):
        products = [record for record in records if record['gene'] == gene]
        entries.append({'gene_name': gene, 'status': 'success' if products else 'fail',
                        'n_products': len(products), 'n_standard_products': sum(r['assembly_status'] == 'standard' for r in products),
                        'n_uncertain_products': sum(r['assembly_status'] == 'uncertain' for r in products),
                        'failure_reason': None if products else 'no connectivity',
                        'threshold_diagnostics': [{'threshold': 2, 'dfs_limit_reached': False}],
                        'output_file': f'{variant}_{gene}.fasta' if products else None,
                        'uncertain_products': {'search_limits_reached': False, 'retention_truncated': False,
                                               'products': [{'marker_detail_truncated': False}]}})
    return {'variant': variant, 'sample_id': 'control', 'scope': scope, 'status': status,
            'outputs': {'records': records, 'pcr_results': entries,
                        'counts': {key: 100 for key in runner.capture.COUNTS}}, 'expected_corrected_sha256': 'expected'}


def record(gene='control_gene', sha='shared', status='standard'):
    return {'gene': gene, 'sha256': sha, 'assembly_status': status, 'length': 909}


class AssessmentControls(unittest.TestCase):
    genes = {'control_gene': 'control', runner.CHANGED_GENE: 'changed'}

    def test_same_sequence_status_switch_is_control_change(self):
        rows = runner.assess([job('before', [record()]), job('after', [record(status='uncertain')])], self.genes, [])
        self.assertEqual(rows[0]['status'], 'control_target_changes_require_review')

    def test_duplicate_count_change_is_not_hidden_by_set_comparison(self):
        rows = runner.assess([job('before', [record()]), job('after', [record(), record()])], self.genes, [])
        self.assertFalse(rows[0]['other_nine_targets_unchanged'])

    def test_reference_hash_at_wrong_target_is_not_synthetic_success(self):
        rows = runner.assess([job('before', [], 'synthetic'), job('after', [record(sha='expected')], 'synthetic')], self.genes, [])
        self.assertEqual(rows[0]['status'], 'source_reconstruction_not_observed')

    def test_exact_uncertain_source_product_is_reconstruction_not_sample_phase(self):
        rows = runner.assess([job('before', [], 'synthetic'),
                              job('after', [record(runner.CHANGED_GENE, 'expected', 'uncertain')], 'synthetic')],
                             self.genes, [{'accession': 'source.1', 'provenance': {'sequence_sha256': 'expected'}}])
        self.assertTrue(rows[0]['exact_corrected_source_emitted_after'])
        self.assertEqual(rows[0]['target_reference_matches'][0]['sample_full_phase'], 'not_established')

    def test_failed_execution_is_not_zero_products(self):
        rows = runner.assess([job('before', []), job('after', [], status='failed')], self.genes, [])
        self.assertEqual(rows[0]['status'], 'incomplete_execution')

    def test_synthetic_control_change_requires_review_despite_source_match(self):
        rows = runner.assess([job('before', [record()], 'synthetic'),
                              job('after', [record(runner.CHANGED_GENE, 'expected')], 'synthetic')], self.genes, [])
        self.assertTrue(rows[0]['exact_corrected_source_emitted_after'])
        self.assertEqual(rows[0]['status'], 'control_target_changes_require_review')

    def test_each_aggregate_count_change_requires_review(self):
        for key in runner.capture.COUNTS:
            with self.subTest(key=key):
                before, after = job('before', []), job('after', [])
                after['outputs']['counts'][key] += 1
                rows = runner.assess([before, after], self.genes, [])
                self.assertEqual(rows[0]['status'], 'aggregate_counts_require_review')
                self.assertTrue(rows[0]['other_nine_targets_unchanged'])

    def test_zero_product_failure_or_threshold_change_requires_review(self):
        for field, value in [('failure_reason', 'DFS limit'),
                             ('threshold_diagnostics', [{'threshold': 2, 'dfs_limit_reached': True}])]:
            with self.subTest(field=field):
                before, after = job('before', []), job('after', [])
                after['outputs']['pcr_results'][0][field] = value
                rows = runner.assess([before, after], self.genes, [])
                self.assertTrue(rows[0]['other_nine_literal_outputs_unchanged'])
                self.assertFalse(rows[0]['other_nine_diagnostics_unchanged'])
                self.assertEqual(rows[0]['status'], 'control_target_changes_require_review')

    def test_uncertainty_search_retention_and_nested_marker_flags_are_compared(self):
        for flag in ('search_limits_reached', 'retention_truncated', 'marker_detail_truncated'):
            with self.subTest(flag=flag):
                before, after = job('before', []), job('after', [])
                uncertainty = after['outputs']['pcr_results'][0]['uncertain_products']
                if flag == 'marker_detail_truncated':
                    uncertainty['products'][0][flag] = True
                else:
                    uncertainty[flag] = True
                rows = runner.assess([before, after], self.genes, [])
                self.assertEqual(rows[0]['status'], 'control_target_changes_require_review')

    def test_run_specific_output_filenames_do_not_make_false_control_change(self):
        rows = runner.assess([job('before', [record()]), job('after', [record()])], self.genes, [])
        self.assertEqual(rows[0]['status'], 'pass')
        self.assertEqual(rows[0]['partial_reference_assessment'], 'pending_post_run')


if __name__ == '__main__':
    unittest.main()

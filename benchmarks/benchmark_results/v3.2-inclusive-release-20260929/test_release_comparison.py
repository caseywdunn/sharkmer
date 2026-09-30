import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import importlib.util

import run_release_comparison as comparison


def record(sequence, label='standard', gene='p_g', scope='high_copy_candidate'):
    return comparison.sequence_record({'header': 'id assembly_status=' + label, 'sequence': sequence}, gene, label, scope)


class ReleaseComparisonControls(unittest.TestCase):
    def test_released_identity_gate_is_mandatory_and_candidate_changes_are_reported(self):
        expectation = {'panel': 'p', 'accession': 'SRR1', 'roles': {role: {'by_gene': {'g': [{'sequence': 'AACG'}], 'empty': []}} for role in ('baseline', 'candidate')}}
        rows = [record('AACG')]
        self.assertTrue(comparison.historical_identity_check(rows, expectation, 'baseline')['ordered_standard_sequences_identical'])
        for invalid in ([record('AACG', 'uncertain')], [record('CGTT')], [record('AACG'), record('AACG')], [record('AACG', gene='p_wrong')]):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                comparison.historical_identity_check(invalid, expectation, 'baseline')
        result = comparison.historical_identity_check([record('AGGG')], expectation, 'candidate')
        self.assertFalse(result['ordered_standard_sequences_identical'])
        self.assertEqual(len(result['missing_literal_identities']), 1)
        self.assertEqual(len(result['gained_literal_identities']), 1)
        self.assertTrue(comparison.historical_identity_check(rows + [record('AGGG', 'uncertain')], expectation, 'candidate')['ordered_standard_sequences_identical'])

    def test_dispositions_preserve_classes_orientation_and_gene(self):
        old = [record(sequence) for sequence in ('AACG', 'AAGC', 'ACCC', 'AGGG', 'AACA')]
        new = [record('AACG'), record('GCTT'), record('ACCC', 'uncertain'), record('CCCT', 'uncertain'), record('AACA', gene='p_other')]
        actual = {row['sha256']: row['disposition'] for row in comparison.dispositions(old, new)}
        self.assertEqual([actual[row['sha256']] for row in old], ['standard_exact', 'standard_rc', 'uncertain_exact', 'uncertain_rc', 'unobserved'])
        self.assertEqual(len(comparison.dispositions(old + [old[0]], new)), 5)

    def test_classification_rejects_ambiguous_headers_and_invalid_sequence(self):
        for header in ('id', 'id assembly_status=unknown', 'id assembly_status=standard assembly_status=uncertain'):
            with self.assertRaises(ValueError):
                comparison.status_of(header)
        self.assertEqual(comparison.status_of('id assembly_status=uncertain'), 'uncertain')
        with self.assertRaises(ValueError):
            record('ACGN')

    def test_legacy_half_median_and_file_stats_consistency(self):
        import yaml
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            command = ['sharkmer', 'reads.fastq']
            stats = {'sample': 's', 'kmer_length': 19, 'chunks': 0, 'command': ' '.join(command), 'sharkmer_version': '3.1.0',
                     'n_reads_read': 1, 'n_bases_read': 4, 'n_subreads_ingested': 1, 'n_bases_ingested': 4, 'n_kmers': 0,
                     'peak_memory_bytes': 1000, 'pcr_results': [{'gene_name': 'p_g', 'status': 'success', 'n_products': 1, 'product_lengths': [4]}]}
            (directory / 's.stats.yaml').write_text(yaml.safe_dump(stats))
            path = directory / 's_p_g.fasta'
            path.write_text('>s_p_g_0 sample=s gene=p_g product=0 length=4 kmer_count_mean=2.5 kmer_count_median=2.5 kmer_count_min=2 kmer_count_max=3 score=0.1\nAACG\n')
            arguments = (directory, {'sample_prefix': 's', 'version': 'baseline'}, {'expected_prefix': {'records': 1, 'bases': 4}}, command, {'p_g': 'high_copy_candidate'})
            actual = comparison.capture(*arguments)
            self.assertEqual(actual['records'][0]['sequence'], 'AACG')
            path.write_text(path.read_text().replace('gene=p_g', 'gene=wrong'))
            with self.assertRaises(ValueError):
                comparison.capture(*arguments)
            path.write_text(path.read_text().replace('gene=wrong', 'gene=p_g'))
            (directory / 'foreign.fasta').write_text('>foreign\nACGT\n')
            with self.assertRaises(ValueError):
                comparison.capture(*arguments)

    def test_summary_requires_all_repeats_and_does_not_multiply_identities(self):
        jobs = []
        for sample in range(10):
            for pair in (1, 2, 3):
                for role in ('baseline', 'candidate'):
                    records = [record('AACG')]
                    if role == 'candidate':
                        records += [record('AGGG', 'uncertain'), record('AGGG', 'uncertain')]
                    jobs.append({'status': 'pass', 'invocation': {'input': str(sample), 'version': role, 'pair_index': pair},
                                 'outputs': {'records': records, 'counts': {'n_reads_read': 1}, 'pcr_results': []},
                                 'execution': {'metrics': {'wall_seconds': 2 if role == 'baseline' else 3, 'peak_rss_bytes': 1000}}})
        result = comparison.summarize(jobs)
        self.assertEqual(result['totals']['candidate'], {'standard': 10, 'uncertain': 20})
        self.assertEqual(result['cells'][0]['roles']['candidate']['unique_counts']['uncertain'], 1)
        self.assertEqual(len(result['cells'][0]['literal_gains']), 1)
        self.assertEqual(result['cells'][0]['resource_changes']['wall_seconds']['percent'], 50)
        self.assertEqual(result['cells'][0]['paired_resources'][0]['wall_seconds']['ratio'], 1.5)
        with self.assertRaises(ValueError):
            comparison.summarize(jobs[:-1])
        changed = copy.deepcopy(jobs)
        changed[0]['outputs']['records'][0]['sequence'] = 'AAAA'
        with self.assertRaises(ValueError):
            comparison.summarize(changed)
        changed = copy.deepcopy(jobs)
        for job in changed:
            if job['invocation']['version'] == 'candidate':
                job['outputs']['counts']['n_reads_read'] = 2
        with self.assertRaises(ValueError):
            comparison.summarize(changed)

    def test_measured_child_and_constraints_are_bound(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'job'
            invocation = {'directory': str(directory), 'outputs': str(directory / 'outputs'),
                          'command': [str(Path(sys.executable).resolve()), '-c', 'print("controlled child")']}
            cpus = sorted(os.sched_getaffinity(0))
            result = comparison.run_timed(invocation, cpus, dict(os.environ, LC_ALL='C'))
            self.assertEqual(result['status'], 'pass')
            observed = json.loads((directory / 'process_constraints.json').read_text())
            self.assertEqual(observed['command'][0], '/usr/bin/time')
            self.assertEqual(observed['command'][-3:], invocation['command'])
            self.assertEqual((directory / 'stdout.log').read_text().strip(), 'controlled child')

    def test_failed_and_timed_out_children_are_not_measurements(self):
        for code, timeout in [('raise SystemExit(3)', 10), ('import time; time.sleep(30)', 0.1)]:
            with self.subTest(code=code), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary) / 'job'
                invocation = {'directory': str(directory), 'outputs': str(directory / 'outputs'),
                              'command': [str(Path(sys.executable).resolve()), '-c', code]}
                with mock.patch.object(comparison, 'TIMEOUT', timeout):
                    result = comparison.run_timed(invocation, sorted(os.sched_getaffinity(0)), dict(os.environ, LC_ALL='C'))
                self.assertEqual(result['status'], 'failed')
                self.assertNotIn('metrics', result)
                if timeout < 1:
                    self.assertTrue(result['timed_out'])
                else:
                    self.assertEqual(result['returncode'], 3)

    def test_candidate_capture_checks_mixed_uncertain_only_empty_and_opt_out(self):
        # Borrow publication fixtures explicitly from the output implementation's tests.
        spec = importlib.util.spec_from_file_location('inclusive_release_fixtures', comparison.ROOT / 'tests/test_uncertain_output.py')
        fixtures = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fixtures)
        for standard, empty, disabled in [(True, False, False), (False, False, False), (False, True, False), (True, False, True)]:
            with self.subTest(standard=standard, empty=empty, disabled=disabled), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                stats = fixtures.combined_transaction(directory, standard=standard, empty=empty, disabled=disabled)
                command = ['sharkmer', 'reads.fastq']
                stats.update(command=' '.join(command), chunks=0, n_reads_read=1, n_bases_read=4,
                             n_subreads_ingested=1, n_bases_ingested=4, n_kmers=0, peak_memory_bytes=1000,
                             input_source={'kind': 'local_files', 'inputs': ['reads.fastq'], 'paired': False, 'max_reads': 1000000},
                             stage_timings={'read_threading_input_s': None})
                files = [row['output_file'] for row in stats['pcr_results'] if row.get('output_file')]
                fixtures.publish(directory, stats, files)
                arguments = (directory, {'sample_prefix': 'sample', 'version': 'candidate'},
                             {'expected_prefix': {'records': 1, 'bases': 4}, 'available_path': 'reads.fastq'}, command, {'gene': 'high_copy_candidate'})
                captured = comparison.capture(*arguments)
                self.assertEqual(sum(row['assembly_status'] == 'standard' for row in captured['records']), int(standard))
                self.assertEqual(sum(row['assembly_status'] == 'uncertain' for row in captured['records']), int(not empty and not disabled))
                if not empty and not disabled:
                    path = directory / 'sample_gene.fasta'
                    path.write_text(path.read_text().replace('assembly_status=uncertain', 'assembly_status=standard'))
                    fixtures.publish(directory, stats, files)
                    with self.assertRaises(ValueError):
                        comparison.capture(*arguments)


if __name__ == '__main__':
    unittest.main()

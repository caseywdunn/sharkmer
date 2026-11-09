#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path('/tmp/sharkmer-withheld-diagnostics-20260913/run_three_way.py')
SPECIFICATION = importlib.util.spec_from_file_location('run_three_way', MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPECIFICATION)
assert SPECIFICATION.loader is not None
SPECIFICATION.loader.exec_module(MODULE)


class HarnessTests(unittest.TestCase):
    def test_protocol_has_twelve_distinct_jobs(self) -> None:
        protocol = json.loads(MODULE.PROTOCOL_PATH.read_text())
        jobs = MODULE.build_jobs(protocol)
        self.assertEqual(len(jobs), 12)
        self.assertEqual(len({(job['accession'], job['role']) for job in jobs}), 12)
        self.assertEqual([job['role'] for job in jobs[:3]], ['baseline_off', 'changed_off', 'changed_on'])
        self.assertEqual([job['role'] for job in jobs[9:]], ['changed_on_focused'] * 3)

    def test_stable_receipt_ignores_capture_time_only(self) -> None:
        before = {'hashes': {'scanner': 'one'}, 'diagnostic_args': ['--diagnostic'], 'inputs': {'input': 'one'}, 'reviewed_files': {'source': 'one'}, 'captured_at': 'first'}
        after = {**before, 'captured_at': 'second'}
        self.assertEqual(MODULE.stable_receipt(before), MODULE.stable_receipt(after))

    def test_prefix_detects_content_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq = Path(temporary) / 'reads.fastq'
            fastq.write_text('@read\nACGT\n+\n!!!!\n')
            observed = MODULE.inspect_prefix(fastq, 1)
            self.assertEqual(observed['sha256'], MODULE.inspect_prefix(fastq, 1)['sha256'])
            fastq.write_text('@read\nACGA\n+\n!!!!\n')
            self.assertNotEqual(observed['sha256'], MODULE.inspect_prefix(fastq, 1)['sha256'])

    def test_optional_payload_only_is_stripped(self) -> None:
        baseline = {'threshold': 3, 'nested': {'withheld_path_diagnostics': {'nodes': 1}}}
        changed = {'threshold': 3, 'nested': {'withheld_path_diagnostics': {'nodes': 2}}}
        different = {'threshold': 4, 'nested': {'withheld_path_diagnostics': {'nodes': 2}}}
        self.assertEqual(MODULE.strip_optional_payload(baseline, 'withheld_path_diagnostics'), MODULE.strip_optional_payload(changed, 'withheld_path_diagnostics'))
        self.assertNotEqual(MODULE.strip_optional_payload(baseline, 'withheld_path_diagnostics'), MODULE.strip_optional_payload(different, 'withheld_path_diagnostics'))

    def test_full_comparator_rejects_nonpayload_change(self) -> None:
        outputs = {'stats': {'n_reads_read': 1, 'n_bases_read': 2, 'n_subreads_ingested': 1, 'n_bases_ingested': 2, 'n_kmers': 3, 'pcr_results': []}, 'fasta_records_ordered': []}
        changed = {'stats': {**outputs['stats'], 'n_kmers': 4}, 'fasta_records_ordered': []}
        records = [
            {'accession': 'sample', 'role': 'baseline_off', 'focused_diagnostic_only': False, 'outputs': outputs},
            {'accession': 'sample', 'role': 'changed_off', 'focused_diagnostic_only': False, 'outputs': changed},
            {'accession': 'sample', 'role': 'changed_on', 'focused_diagnostic_only': False, 'outputs': outputs},
        ]
        result = MODULE.compare_full_roles(records, 'withheld_path_diagnostics')
        self.assertFalse(result[0]['passed'])
        self.assertIn('baseline_vs_changed_off:n_kmers', result[0]['failures'])

    def test_full_comparator_rejects_gene_status_and_failure_changes(self) -> None:
        stats = {'n_reads_read': 1, 'n_bases_read': 2, 'n_subreads_ingested': 1, 'n_bases_ingested': 2, 'n_kmers': 3, 'pcr_results': [{'gene_name': 'gene', 'status': 'fail', 'n_products': 0, 'failure_reason': 'first'}]}
        changed = {**stats, 'pcr_results': [{'gene_name': 'gene', 'status': 'success', 'n_products': 1, 'output_file': 'sample_gene.fasta'}]}
        records = [
            {'accession': 'sample', 'role': 'baseline_off', 'focused_diagnostic_only': False, 'outputs': {'stats': stats, 'fasta_records_ordered': []}},
            {'accession': 'sample', 'role': 'changed_off', 'focused_diagnostic_only': False, 'outputs': {'stats': changed, 'fasta_records_ordered': []}},
            {'accession': 'sample', 'role': 'changed_on', 'focused_diagnostic_only': False, 'outputs': {'stats': stats, 'fasta_records_ordered': []}},
        ]
        result = MODULE.compare_full_roles(records, 'withheld_path_diagnostics')
        self.assertIn('baseline_vs_changed_off:pcr_results', result[0]['failures'])

    def test_diagnostic_payload_rejects_wrong_placement_and_off_mode(self) -> None:
        wrong = {'withheld_path_diagnostics': [], 'pcr_results': []}
        with self.assertRaises(ValueError):
            MODULE.validate_diagnostic_payload(wrong, 'withheld_path_diagnostics', True)
        off = {'pcr_results': [{'threshold_diagnostics': [{'withheld_path_diagnostics': []}]}]}
        with self.assertRaises(ValueError):
            MODULE.validate_diagnostic_payload(off, 'withheld_path_diagnostics', False)

    def test_diagnostic_payload_rejects_nested_threshold_field(self) -> None:
        nested = {'pcr_results': [{'threshold_diagnostics': [{'nested': {'withheld_path_diagnostics': []}}]}]}
        with self.assertRaises(ValueError):
            MODULE.validate_diagnostic_payload(nested, 'withheld_path_diagnostics', True)

    def test_missing_current_manifest_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                MODULE.capture_outputs(
                    Path(temporary), 'sample', Path('/home/claude/repos/sharkmer/panels/insecta.yaml'), 19,
                    ['sharkmer'], Path('/tmp/input.fastq'), {'records': 1, 'bases': 1, 'size_bytes': 1, 'sha256': 'one'}, False,
                    'withheld_path_diagnostics',
                )


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python3
"""Synthetic controls only: no biological inputs or assembler invocation."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import subprocess
import unittest
from unittest.mock import patch

import run_comparison as run
import build_snapshot as builder


def identity(data, records=1, bases=4):
    return {'records': records, 'bases': bases, 'size_bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def fixtures():
    records = []
    for accession, roles in run.ORDER:
        for role in roles:
            gene = {'gene_name': 'insecta_test', 'status': 'success', 'n_products': 1, 'output_file': 'same.fasta', 'threshold_diagnostics': [{'threshold': 2, 'accepted_paths': 1}]}
            if role == 'changed_on':
                gene['uncertain_products'] = {'limits': {'paths_per_gene': 12}, 'retained_path_records': 2, 'unique_products': 1, 'search_limits_reached': True, 'retention_truncated': True}
            records.append({'accession': accession, 'role': role, 'outputs': {'stats': {**{key: 1 for key in run.COUNTS}, 'pcr_results': [gene]}, 'standard': [{'gene': 'insecta_test', 'header': 'standard', 'sequence': 'ACGT'}], 'uncertain': [{'gene': 'insecta_test', 'header': 'uncertain', 'sequence': 'AATC'}] if role == 'changed_on' else []}})
    losses = [{'accession': run.ORDER[0][0], 'gene': 'test', 'length': 4, 'sequence': 'AATC', 'sequence_sha256': hashlib.sha256(b'AATC').hexdigest()}]
    return records, losses


class Controls(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_prefix_exact(self):
        data = b'@a\nACGT\n+\nIIII\n'
        path = self.root / 'reads'
        path.write_bytes(data)
        self.assertEqual(run.inspect_prefix(path, identity(data)), identity(data))

    def test_prefix_rejects_truncated_malformed_extra_and_identity(self):
        good = b'@a\nACGT\n+\nIIII\n'
        for data in (good[:-1], good.replace(b'+', b'-'), good + good, good.replace(b'ACGT', b'ACZT'), good.replace(b'IIII', b'II I'), good.replace(b'ACGT', b'AGGT')):
            with self.subTest(data=data):
                path = self.root / 'reads'
                path.write_bytes(data)
                with self.assertRaises(ValueError):
                    run.inspect_prefix(path, identity(good))

    def test_no_clobber_json(self):
        path = self.root / 'receipt'
        run.write_new(path, {'one': 1})
        before = path.read_bytes()
        with self.assertRaises(FileExistsError):
            run.write_new(path, {'two': 2})
        self.assertEqual(path.read_bytes(), before)

    def test_fingerprint_mutation(self):
        path = self.root / 'source'
        path.write_text('before')
        row = run.fingerprint(path)
        path.write_text('after')
        with self.assertRaises(ValueError):
            run.verify(row)

    def test_review_requires_pass_nonempty_and_coverage(self):
        path, review = self.root / 'source', self.root / 'review'
        path.write_text('code')
        for document in ({'status': 'pass', 'reviewed_sources': {}}, {'status': 'fail', 'reviewed_sources': {'s': run.fingerprint(path)}}, {'status': 'pass', 'reviewed_sources': {'s': run.fingerprint(path)}}):
            review.write_text(json.dumps(document))
            with self.assertRaises(ValueError):
                run.review_sources(review, {self.root / 'missing'})
        self.assertTrue(run.review_sources(review, {path}))

    def test_exact_nine_commands_flags_order_and_limits(self):
        samples = {accession: {'path': str(self.root / accession)} for accession, roles in run.ORDER}
        builds = {role: {'binary': {'path': str(self.root / role)}} for role in ('baseline', 'candidate')}
        jobs = run.commands(samples, builds, self.root, [4, 7])
        self.assertEqual([(job['accession'], job['role']) for job in jobs], [(acc, role) for acc, roles in run.ORDER for role in roles])
        for job in jobs:
            self.assertEqual('--emit-uncertain' in job['argv'], job['role'] == 'changed_on')
            self.assertNotIn('--read-threading', job['argv'])
            self.assertNotIn('--diagnose-withheld-paths', job['argv'])
            self.assertIn('--as=42949672960', job['command'])
            self.assertIn('4,7', job['command'])

    def test_parity_and_exact_separate(self):
        records, losses = fixtures()
        output = run.analyze(records, losses)
        self.assertEqual(output['status'], 'pass')
        match = output['historical_losses'][0]['matches']['changed_on']
        self.assertEqual(match['standard']['exact_oriented'], [])
        self.assertEqual(match['uncertain']['exact_oriented'], ['uncertain'])
        self.assertTrue(output['historical_losses'][0]['assessment']['retention_truncated'])

    def test_reverse_complement_separate(self):
        records, losses = fixtures()
        losses[0].update(sequence='GATT', sequence_sha256=hashlib.sha256(b'GATT').hexdigest())
        match = run.analyze(records, losses)['historical_losses'][0]['matches']['changed_on']['uncertain']
        self.assertEqual(match, {'exact_oriented': [], 'reverse_complement_only': ['uncertain']})

    def test_approximate_not_exact(self):
        records, losses = fixtures()
        losses[0].update(sequence='AATT', sequence_sha256=hashlib.sha256(b'AATT').hexdigest())
        match = run.analyze(records, losses)['historical_losses'][0]['matches']['changed_on']['uncertain']
        self.assertEqual(match, {'exact_oriented': [], 'reverse_complement_only': []})

    def test_each_preexisting_parity_dimension_matters(self):
        for field in (*run.COUNTS, 'header', 'sequence', 'diagnostic', 'status'):
            with self.subTest(field=field):
                records, losses = fixtures()
                output = records[1]['outputs']
                if field in run.COUNTS:
                    output['stats'][field] += 1
                elif field in ('header', 'sequence'):
                    output['standard'][0][field] += 'A'
                elif field == 'diagnostic':
                    output['stats']['pcr_results'][0]['threshold_diagnostics'][0]['accepted_paths'] += 1
                else:
                    output['stats']['pcr_results'][0]['status'] = 'fail'
                self.assertEqual(run.analyze(records, losses)['status'], 'fail')

    def test_only_direct_new_payload_removed(self):
        value = {'pcr_results': [{'uncertain_products': {}, 'threshold_diagnostics': [{'uncertain_products': 'unexpected'}]}]}
        self.assertEqual(run.standard_stats(value), [{'threshold_diagnostics': [{'uncertain_products': 'unexpected'}]}])

    def test_missing_and_duplicate_jobs_rejected(self):
        records, losses = fixtures()
        for altered in (records[:-1], records + records[:1]):
            with self.assertRaises(ValueError):
                run.analyze(altered, losses)

    def test_loss_identity_corruption_rejected(self):
        records, losses = fixtures()
        losses[0]['sequence_sha256'] = '0' * 64
        with self.assertRaises(ValueError):
            run.analyze(records, losses)

    def test_frozen_commands_tamper_rejected(self):
        import yaml
        source = self.root / 'source'
        source.write_text('reviewed source')
        review = self.root / 'review'
        reviewed = {'source': run.fingerprint(source)}
        review.write_text(json.dumps({'status': 'pass', 'reviewed_sources': reviewed}))
        samples = {accession: {'path': str(self.root / accession)} for accession, roles in run.ORDER}
        builds = {role: {'binary': {'path': str(self.root / role)}} for role in ('baseline', 'candidate')}
        receipts = {}
        for role in builds:
            path = self.root / f'{role}.build.json'
            path.write_text('{}')
            receipts[role] = run.fingerprint(path)
        cpus = [0, 1]
        paths = {source, review, Path(sys.executable).resolve(), *run.TOOLS.values(), *(Path(row['path']) for row in receipts.values())}
        frozen = {'schema_version': 1, 'status': 'frozen_before_execution', 'review': run.fingerprint(review), 'reviewed_sources': reviewed, 'build_receipts': receipts, 'sources': {str(path): run.fingerprint(path) for path in paths}, 'limits': run.LIMITS, 'samples': samples, 'builds': builds, 'output': str(self.root), 'cpus': cpus, 'runtime': {'python': sys.version, 'interpreter': run.fingerprint(Path(sys.executable).resolve()), 'yaml_version': yaml.__version__}, 'execution_environment': {'LC_ALL': 'C'}}
        frozen['jobs'] = run.commands(samples, builds, self.root, cpus)
        with patch.object(run.os, 'sched_getaffinity', return_value={0, 1}), patch.object(run, 'build_info', side_effect=lambda path, role: (builds[role], {})), patch.object(run, 'required_sources', return_value={source}), patch.object(run, 'same_build_configuration'), patch.object(run, 'verify_sample_plan'):
            run.check_frozen(frozen)
            for alteration in ('command', 'empty_sources', 'missing_source', 'empty_review', 'schema', 'duplicate_cpu', 'boolean_cpu', 'extra_sample'):
                changed = copy.deepcopy(frozen)
                if alteration == 'command': changed['jobs'][0]['argv'].insert(-1, '--emit-uncertain')
                elif alteration == 'empty_sources': changed['sources'] = {}
                elif alteration == 'missing_source': del changed['sources'][str(source)]
                elif alteration == 'empty_review': changed['reviewed_sources'] = {}
                elif alteration == 'schema': changed['schema_version'] = True
                elif alteration == 'duplicate_cpu': changed['cpus'] = [0, 0]
                elif alteration == 'boolean_cpu': changed['cpus'] = [False, 1]
                elif alteration == 'extra_sample': changed['samples']['extra'] = {}
                with self.subTest(alteration=alteration), self.assertRaises(ValueError):
                    run.check_frozen(changed)

    def test_preparation_receipts_reject_derived_input_mutation(self):
        inputs = self.root / 'inputs'
        inputs.mkdir()
        preparer = self.root / 'prepare_inputs.py'
        preparer.write_text('original preparer')
        outputs = {}
        for index in range(12):
            path = inputs / str(index)
            path.write_text(str(index))
            outputs[str(index)] = {key: value for key, value in run.fingerprint(path).items() if key != 'path'}
        receipt = {'status': 'verified', 'outputs': outputs, 'sources': {}, 'preparation_source': run.fingerprint(preparer)}
        (inputs / 'provenance.json').write_text(json.dumps(receipt))
        with patch.object(run, 'HERE', self.root):
            self.assertEqual(len(run.prepared_sources()), 12)
            (inputs / '3').write_text('changed')
            with self.assertRaises(ValueError):
                run.prepared_sources()

    def test_optional_kmer_count_presence_and_values_are_compared(self):
        for key in ('n_multi_kmers', 'n_singleton_kmers'):
            records, losses = fixtures()
            records[1]['outputs']['stats'][key] = 0
            self.assertEqual(run.analyze(records, losses)['status'], 'fail')
            for record in records: record['outputs']['stats'][key] = 0
            self.assertEqual(run.analyze(records, losses)['status'], 'pass')
            records[1]['outputs']['stats'][key] = 1
            self.assertEqual(run.analyze(records, losses)['status'], 'fail')

    def test_symlink_source_is_rejected(self):
        source, link = self.root / 'source', self.root / 'link'
        source.write_text('source')
        link.symlink_to(source)
        with self.assertRaises(ValueError): run.fingerprint(link)

    def test_missing_compiler_creates_failed_build_receipt(self):
        output = self.root / 'build'
        source = output / 'source'
        source.mkdir(parents=True)
        (source / 'Cargo.toml').write_text('synthetic')
        snapshot = {'source': str(source), 'files': builder.tree(source), 'exporter': run.fingerprint(Path(builder.__file__))}
        (output / 'snapshot.json').write_text(json.dumps(snapshot))
        with self.assertRaises(RuntimeError): builder.build(output, self.root / 'missing-cargo', self.root / 'missing-rustc')
        receipt = json.loads((output / 'build.json').read_text())
        self.assertEqual(receipt['status'], 'failed')
        self.assertIn('FileNotFoundError', receipt['error'])
        self.assertEqual(receipt['source_before'], receipt['source_after'])

    def test_candidate_export_prunes_only_generated_panel_outputs_and_caches(self):
        root = self.root / 'repository'
        root.mkdir()
        meaningful = ['Cargo.toml', 'Cargo.lock', 'src/main.rs', 'panels/insecta.yaml', 'panels/reference_sources.json.gz', 'panels/examples/example.yaml', 'tests/uncertain_products.rs', 'tests/fixtures/input.fastq', 'tests/validation_results/meaningful_fixture.txt']
        excluded = ['panels/validation_runs/old.fastq', 'panels/validation_results/old.json', 'panels/validation_reports/old.html', 'tests/__pycache__/test.pyc']
        for name in meaningful + excluded:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
        with patch.object(builder, 'ROOT', root), patch.object(builder.subprocess, 'check_output', return_value='synthetic-head\n'):
            result = builder.export('candidate', self.root / 'export')
        self.assertEqual(set(result['files']), set(meaningful))

    def test_launcher_observes_enforced_limits_and_rejects_wrong_expectation(self):
        cpu = min(os.sched_getaffinity(0))
        cap = 256 * 1024**2
        for wrong in (False, True):
            receipt = self.root / f'launch-{wrong}.json'
            command = [str(run.TOOLS['prlimit']), f'--as={cap}', '--', str(run.TOOLS['taskset']), '-c', str(cpu), str(Path(sys.executable).resolve()), str(run.HERE / 'launch_bounded.py'), '--receipt', str(receipt), '--expected-cpus', str(cpu), '--expected-as', str(cap + int(wrong)), '--', str(Path(sys.executable).resolve()), '-c', 'print("synthetic")']
            result = subprocess.run(command, capture_output=True, text=True, timeout=5)
            observed = json.loads(receipt.read_text())
            self.assertEqual(observed['affinity'], [cpu])
            self.assertEqual(observed['rlimit_as_soft_hard'], [cap, cap])
            self.assertEqual(observed['status'], 'fail' if wrong else 'pass')
            self.assertEqual(result.returncode == 0, not wrong)

    def test_external_freeze_pin_rejected_before_output(self):
        path = self.root / 'freeze'
        path.write_text('{}')
        with self.assertRaises(ValueError):
            run.execute(argparse.Namespace(freeze=path, freeze_sha256='0' * 64))

    def test_execute_failed_gate_has_receipt_and_no_clobber(self):
        path, output = self.root / 'freeze', self.root / 'output'
        path.write_text(json.dumps({'output': str(output)}))
        args = argparse.Namespace(freeze=path, freeze_sha256=run.fingerprint(path)['sha256'])
        with patch.object(run, 'check_frozen', side_effect=ValueError('gate failed')):
            self.assertEqual(run.execute(args), 1)
        self.assertEqual(json.loads((output / 'execution.json').read_text())['status'], 'failed')
        with self.assertRaises(FileExistsError):
            run.execute(args)

    def test_failed_process_receipt_and_input_bookends(self):
        path = self.root / 'reads'
        data = b'@a\nACGT\n+\nIIII\n'
        path.write_bytes(data)
        sample = {'path': str(path), 'expected_prefix': identity(data)}
        job = {'output': str(self.root / 'out'), 'command': [sys.executable, '-c', 'raise SystemExit(9)']}
        with patch.object(run, 'check_frozen'):
            receipt = run.run_job(job, sample, {})
        self.assertEqual(receipt['returncode'], 9)
        self.assertEqual(receipt['status'], 'failed')
        self.assertEqual(receipt['prefix_before'], receipt['prefix_after'])

    def test_process_timeout_kills_group_and_preserves_receipt(self):
        path = self.root / 'reads'
        data = b'@a\nACGT\n+\nIIII\n'
        path.write_bytes(data)
        sample = {'path': str(path), 'expected_prefix': identity(data)}
        job = {'output': str(self.root / 'out'), 'command': [sys.executable, '-c', 'import time; time.sleep(20)']}
        with patch.object(run, 'check_frozen'), patch.dict(run.LIMITS, {'timeout_seconds': 0.05}):
            receipt = run.run_job(job, sample, {})
        self.assertTrue(receipt['timed_out'])
        self.assertEqual(receipt['status'], 'failed')
        self.assertLess(receipt['wall_seconds'], 5)


if __name__ == '__main__':
    unittest.main(verbosity=2)

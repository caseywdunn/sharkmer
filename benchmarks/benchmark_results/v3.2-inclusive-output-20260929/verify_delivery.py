#!/usr/bin/env python3
"""Check/build the inclusive interface, then compare six runs with saved outcomes."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
PREVIOUS = ROOT / 'benchmarks/benchmark_results/v3.2-uncertain-products-20260927'
PRIOR_RUNS = ROOT / 'tmp/uncertain-products-20260927/comparison-v1/runs'
sys.path.insert(0, str(PREVIOUS))
import run_compiled_checks as checks
from build_snapshot import fingerprint, write_new
from run_comparison import COUNTS, fasta, inspect_prefix, verify
sys.path.insert(0, str(ROOT / 'scripts'))
from sharkmer_validate import runner


def run(command, cwd, log, environment, timeout):
    result = {'command': command, 'cwd': str(cwd), 'status': 'failed'}
    started = time.monotonic()
    with log.open('x') as stream:
        process = subprocess.Popen(command, cwd=cwd, env=environment,
                                   stdout=stream, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            result['returncode'] = process.wait(timeout=timeout)
            result['status'] = 'pass' if result['returncode'] == 0 else 'failed'
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            result.update(timed_out=True, returncode=process.returncode)
    result.update(wall_seconds=time.monotonic() - started, log=fingerprint(log))
    return result


def require_pass(result):
    if result['status'] != 'pass':
        raise ValueError(f"Command failed; preserve its log: {result['log']['path']}")


def class_of(header):
    values = [token.split('=', 1)[1] for token in shlex.split(header)
              if token.startswith('assembly_status=')]
    if len(values) != 1 or values[0] not in ('standard', 'uncertain'):
        raise ValueError('Each new record must have exactly one valid assembly status')
    return values[0]


def old_standard_header(header):
    if class_of(header) != 'standard':
        raise ValueError('Expected a standard record')
    tokens = header.split(' ')
    tokens.remove('assembly_status=standard')
    return ' '.join(tokens)


def assessment_identity(assessment):
    value = copy.deepcopy(assessment)
    value.pop('output_file', None)
    for product in value['products']:
        product['id'] = product['sequence_sha256']
    return value


def fasta_text(text):
    records = []
    for line in text.splitlines():
        if line.startswith('>'):
            records.append({'header': line[1:], 'sequence': ''})
        elif line:
            if not records:
                raise ValueError('Sequence before header in filtered FASTA')
            records[-1]['sequence'] += line.strip()
    if any(not record['sequence'] for record in records):
        raise ValueError('Empty filtered FASTA record')
    return records


def capture(directory, sample, command, prior, include):
    name = f"uncertain_{sample['accession']}"
    stats = runner._parse_stats_yaml(directory / f'{name}.stats.yaml')
    genes = {entry['gene_name'] for entry in prior['stats']['pcr_results']}
    runner._validate_stats_manifest(stats, name, 19, genes, 'end-to-end',
                                    expected_command=command)
    manifest = runner._validate_output_manifest(stats, name, directory)
    assert manifest['schema_version'] == 2
    assert stats['product_output_schema_version'] == 2
    assert stats['product_output_layout'] == 'combined_classified'
    assert len(stats['pcr_results']) == len(genes) == 21
    assert stats['chunks'] == 0
    assert stats['input_source'] == {'kind': 'local_files', 'inputs': [sample['file']['path']],
                                     'paired': False, 'max_reads': 1000000}
    assert stats['stage_timings']['read_threading_input_s'] is None
    for key in (*COUNTS, 'n_multi_kmers', 'n_singleton_kmers'):
        assert (key in stats) == (key in prior['stats'])
        assert stats.get(key) == prior['stats'].get(key), key
    old_genes = {entry['gene_name']: entry for entry in prior['stats']['pcr_results']}
    standard, uncertain, gene_rows = [], [], []
    for entry in stats['pcr_results']:
        gene = entry['gene_name']
        old = old_genes[gene]
        records = fasta(directory / entry['output_file']) if entry.get('output_file') else []
        standard_records = [row for row in records if class_of(row['header']) == 'standard']
        uncertain_records = [row for row in records if class_of(row['header']) == 'uncertain']
        assert records == standard_records + uncertain_records
        assert entry['n_products'] == len(records)
        assert entry.get('product_lengths', []) == [len(row['sequence']) for row in records]
        assert entry['n_standard_products'] == len(standard_records) == old['n_products']
        assert entry['n_uncertain_products'] == len(uncertain_records)
        assert bool(entry.get('output_file')) == bool(records)
        assert entry['status'] == ('success' if standard_records else 'uncertain' if records else 'fail')
        assert entry.get('failure_reason') == old.get('failure_reason')
        assert entry.get('threshold_diagnostics', []) == old.get('threshold_diagnostics', [])
        old_standard = [row for row in prior['standard'] if row['gene'] == gene]
        new_standard = [{'gene': gene, 'header': old_standard_header(row['header']),
                         'sequence': row['sequence']} for row in standard_records]
        assert new_standard == old_standard, f'{gene}: original standard records changed'
        old_uncertain = [row for row in prior['uncertain'] if row['gene'] == gene] if include else []
        assert [row['sequence'] for row in uncertain_records] == [row['sequence'] for row in old_uncertain]
        assert [row['header'].split(' ', 1)[1] for row in uncertain_records] == [
            row['header'].split(' ', 1)[1] for row in old_uncertain]
        if include:
            assert assessment_identity(entry['uncertain_products']) == assessment_identity(old['uncertain_products'])
        else:
            assert 'uncertain_products' not in entry
        if records:
            path = directory / entry['output_file']
            assert path.name == f'{name}_{gene}.fasta'
            for label, expected in [('standard', standard_records), ('uncertain', uncertain_records)]:
                filtered = subprocess.check_output(['awk', f'/^>/ {{keep = / assembly_status={label}( |$)/}} keep', str(path)], text=True)
                assert fasta_text(filtered) == expected, f'{gene}: documented shell filter differs'
        standard.extend({'gene': gene, **row} for row in standard_records)
        uncertain.extend({'gene': gene, **row} for row in uncertain_records)
        gene_rows.append({'gene': gene, 'status': entry['status'], 'standard': len(standard_records),
                          'uncertain': len(uncertain_records), 'total': len(records)})
    ids = [row['header'].split()[0] for row in standard + uncertain]
    assert len(ids) == len(set(ids)), 'Cross-target output IDs collided'
    consumed = runner.parse_fasta_products(name, directory)
    consumed_sequences = [(row['gene'], record['sequence']) for row in consumed for record in row['products']]
    assert consumed_sequences == [(row['gene'], row['sequence']) for row in standard]
    assert runner._validate_output_manifest(stats, name, directory) == manifest
    return {'standard': standard, 'uncertain': uncertain, 'genes': gene_rows,
            'stats': fingerprint(directory / f'{name}.stats.yaml'),
            'manifest': fingerprint(directory / f'{name}.manifest.yaml'),
            'files': {path.name: fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()}}


def main():
    if not __debug__:
        raise SystemExit('Verification must run without Python assertion optimization')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'schema_version': 1, 'status': 'failed', 'jobs': [], 'claims_performance_equivalence': False}
    try:
        cpus = sorted(os.sched_getaffinity(0))[:2]
        assert len(cpus) == 2
        os.sched_setaffinity(0, cpus)
        environment = dict(os.environ)
        rust_bin = Path(environment['RUSTUP_HOME']) / 'toolchains/stable-x86_64-unknown-linux-gnu/bin'
        cargo, rustc = str(rust_bin / 'cargo'), str(rust_bin / 'rustc')
        receipt['source_before'] = checks.sources()
        extra = [Path(__file__).resolve(), HERE / 'verify_delivery.sbatch', HERE / 'test_verify_delivery.py', ROOT / 'AGENTS.md',
                 ROOT / 'dev_docs/INCLUSIVE_OUTPUT.md', PREVIOUS / 'inputs/current_prefix_availability.json',
                 PREVIOUS / 'inputs/insecta.yaml', PREVIOUS / 'inputs/historical_losses.json',
                 PRIOR_RUNS / 'assessment.json', PRIOR_RUNS / 'execution.json']
        receipt['extra_sources'] = {str(path): fingerprint(path) for path in extra}
        receipt['runtime'] = checks.runtime()
        receipt['affinity'] = cpus
        receipt['slurm_job_id'] = os.environ.get('SLURM_JOB_ID')
        step = run([str(Path(sys.executable).resolve()), '-m', 'unittest', 'discover',
                    '-s', str(HERE), '-p', 'test_verify_delivery.py', '-v'], ROOT,
                   output / 'comparison_controls.log', environment, 60)
        receipt['comparison_controls'] = step
        require_pass(step)
        checked = output / 'checks'
        step = run([str(Path(sys.executable).resolve()), str(PREVIOUS / 'run_compiled_checks.py'),
                    '--cargo', cargo, '--rustc', rustc, '--output', str(checked)], ROOT,
                   output / 'checks.log', environment, 14400)
        receipt['checks_step'] = step
        require_pass(step)
        check_receipt = json.loads((checked / 'checks.json').read_text())
        assert check_receipt['status'] == 'pass'
        assert check_receipt['sources_before'] == receipt['source_before'] == checks.sources()
        receipt['checks'] = fingerprint(checked / 'checks.json')
        environment.update(PATH=str(rust_bin) + os.pathsep + environment.get('PATH', ''),
                           RUSTC=rustc, CARGO_TARGET_DIR=str(checked / 'target'), CARGO_BUILD_JOBS='2')
        for key in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER'):
            environment.pop(key, None)
        step = run([cargo, 'build', '--locked', '--offline', '--release', '-j', '2'], ROOT,
                   output / 'build.log', environment, 3600)
        receipt['build_step'] = step
        require_pass(step)
        assert checks.sources() == receipt['source_before']
        binary = output / 'sharkmer'
        shutil.copy2(checked / 'target/release/sharkmer', binary)
        binary.chmod(0o555)
        receipt['binary'] = fingerprint(binary)
        receipt['compiler'] = check_receipt['tools_before']
        for tool in receipt['compiler'].values():
            verify(tool)
        samples = json.loads((PREVIOUS / 'inputs/current_prefix_availability.json').read_text())['samples']
        prior_execution = json.loads((PRIOR_RUNS / 'execution.json').read_text())
        assert prior_execution['status'] == 'complete'
        prior_receipts = {row['id']: row['receipt'] for row in prior_execution['jobs']}
        for sample in samples:
            accession = sample['accession']
            old_path = PRIOR_RUNS / f'{accession}.changed_on.receipt.json'
            verify(prior_receipts[f'{accession}.changed_on'])
            old = json.loads(old_path.read_text())
            assert old['status'] == 'complete'
            for file in old['outputs']['files'].values():
                verify(file)
            for mode in ('default', 'standard_only'):
                directory = output / 'runs' / accession / mode
                directory.mkdir(parents=True)
                job = {'accession': accession, 'mode': mode, 'status': 'failed',
                       'prior_receipt': fingerprint(old_path)}
                receipt['jobs'].append(job)
                try:
                    job['prefix_before'] = inspect_prefix(Path(sample['file']['path']), sample['expected_prefix'])
                    command = [str(binary), '-k', '19', '-t', '2', '--chunks', '0', '--max-reads', '1000000',
                               '-o', str(directory), '-s', f'uncertain_{accession}',
                               '--pcr-panel-file', str(PREVIOUS / 'inputs/insecta.yaml')]
                    if mode == 'standard_only':
                        command.append('--no-emit-uncertain')
                    command.append(sample['file']['path'])
                    job['command'] = command
                    constraints = directory / 'process_constraints.json'
                    bounded = ['/usr/bin/prlimit', '--as=42949672960', '--',
                               str(Path(sys.executable).resolve()), str(PREVIOUS / 'launch_bounded.py'),
                               '--receipt', str(constraints), '--expected-cpus', ','.join(map(str, cpus)),
                               '--expected-as', '42949672960', '--', *command]
                    job['execution'] = run(bounded, ROOT, directory / 'run.log', environment, 1800)
                    require_pass(job['execution'])
                    observed_constraints = json.loads(constraints.read_text())
                    assert observed_constraints['status'] == 'pass'
                    assert observed_constraints['command'] == command
                    assert observed_constraints['affinity'] == cpus
                    assert observed_constraints['rlimit_as_soft_hard'] == [42949672960, 42949672960]
                    job['outputs'] = capture(directory, sample, command, old['outputs'], mode == 'default')
                    job['prefix_after'] = inspect_prefix(Path(sample['file']['path']), sample['expected_prefix'])
                    assert job['prefix_before'] == job['prefix_after']
                    job['status'] = 'pass'
                finally:
                    job_path = output / f'{accession}.{mode}.json'
                    write_new(job_path, job)
                    job['receipt'] = fingerprint(job_path)
                for file in old['outputs']['files'].values():
                    verify(file)
                verify(job['prior_receipt'])
        losses = json.loads((PREVIOUS / 'inputs/historical_losses.json').read_text())
        if isinstance(losses, dict):
            losses = losses['rows']
        receipt['historical_losses'] = []
        for loss in losses:
            rows = next(row['outputs']['uncertain'] for row in receipt['jobs']
                        if row['accession'] == loss['accession'] and row['mode'] == 'default')
            matches = [row['header'] for row in rows if row['gene'] == 'insecta_' + loss['gene']
                       and hashlib.sha256(row['sequence'].encode()).hexdigest() == loss['sequence_sha256']]
            receipt['historical_losses'].append({key: loss[key] for key in ('accession', 'gene', 'length', 'sequence_sha256')} | {'exact_matches': matches})
        assert sum(bool(row['exact_matches']) for row in receipt['historical_losses']) == 5
        receipt['status'] = 'pass'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        try:
            receipt['source_after'] = checks.sources()
            assert receipt.get('source_before') == receipt['source_after']
            for row in receipt.get('extra_sources', {}).values():
                verify(row)
            if 'binary' in receipt:
                verify(receipt['binary'])
            for tool in receipt.get('compiler', {}).values():
                verify(tool)
            if 'checks' in receipt:
                verify(receipt['checks'])
                checked_receipt = json.loads(Path(receipt['checks']['path']).read_text())
                assert checks.cargo_configs(ROOT, os.environ) == checked_receipt['cargo_configs_after']
            for job in receipt['jobs']:
                if 'receipt' in job:
                    verify(job['receipt'])
                for row in job.get('outputs', {}).get('files', {}).values():
                    verify(row)
                if 'prior_receipt' in job:
                    verify(job['prior_receipt'])
                    previous = json.loads(Path(job['prior_receipt']['path']).read_text())
                    for row in previous['outputs']['files'].values():
                        verify(row)
            assert checks.runtime() == receipt['runtime']
        except Exception as error:
            receipt.update(status='failed', bookend_error=repr(error))
        write_new(output / 'verification.json', receipt)
    print(json.dumps({'status': receipt['status'], 'output': str(output), 'error': receipt.get('error')}))
    return 0 if receipt['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

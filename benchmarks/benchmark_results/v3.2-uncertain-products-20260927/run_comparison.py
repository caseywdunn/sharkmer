#!/usr/bin/env python3
"""Freeze, execute and assess the preregistered nine-run uncertain-product comparison."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

from build_snapshot import BASELINE, fingerprint, write_new, tree

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
PROTOCOL = ROOT / 'dev_docs/UNCERTAIN_PRODUCTS.md'
COUNTS = ('n_reads_read', 'n_bases_read', 'n_subreads_ingested', 'n_bases_ingested', 'n_kmers')
ORDER = [('SRR27962769', ['baseline_off', 'changed_off', 'changed_on']), ('SRR31887760', ['changed_on', 'baseline_off', 'changed_off']), ('SRR1057608', ['changed_off', 'changed_on', 'baseline_off'])]
LIMITS = {'address_space_bytes': 40 * 1024**3, 'timeout_seconds': 1800, 'threads': 2, 'k': 19, 'chunks': 0, 'max_reads': 1000000}
TOOLS = {name: Path('/usr/bin', name).resolve() for name in ('time', 'taskset', 'prlimit')}


def verify(row):
    if fingerprint(row['path']) != row:
        raise ValueError(f"Fingerprint changed: {row['path']}")


def inspect_prefix(path, expected):
    h = hashlib.sha256()
    bases = size = 0
    with Path(path).open('rb') as handle:
        for ordinal in range(expected['records']):
            lines = [handle.readline(4097) for unused in range(4)]
            if any(not line.endswith(b'\n') or len(line) > 4096 for line in lines):
                raise ValueError(f'Invalid/truncated FASTQ at record {ordinal + 1}')
            seq, qual = lines[1].rstrip(b'\r\n'), lines[3].rstrip(b'\r\n')
            if not lines[0].startswith(b'@') or not lines[2].startswith(b'+') or not seq or len(seq) != len(qual):
                raise ValueError(f'Malformed FASTQ at record {ordinal + 1}')
            if any(base not in b'ACGTNacgtn' for base in seq) or any(q < 33 or q > 126 for q in qual):
                raise ValueError(f'Invalid sequence/quality at record {ordinal + 1}')
            for line in lines:
                h.update(line)
                size += len(line)
            bases += len(seq)
        if handle.read(1):
            raise ValueError('Staged input must contain only the exact consumed prefix')
    observed = {'records': expected['records'], 'bases': bases, 'size_bytes': size, 'sha256': h.hexdigest()}
    if observed != expected:
        raise ValueError(f'Consumed prefix differs: {path}')
    return observed


def review_sources(review_path, required):
    review = json.loads(review_path.read_text())
    if review.get('status') != 'pass' or not review.get('reviewed_sources'):
        raise ValueError('Independent passing review with nonempty reviewed_sources is required')
    reviewed = review['reviewed_sources']
    paths = set()
    for row in reviewed.values():
        verify(row)
        paths.add(str(Path(row['path']).resolve()))
    missing = {str(path.resolve()) for path in required} - paths
    if missing:
        raise ValueError(f'Review does not bind required sources: {sorted(missing)}')
    return reviewed


def prepared_sources():
    """Recheck derivation receipts instead of trusting modified derived inputs."""
    provenance = json.loads((HERE / 'inputs/provenance.json').read_text())
    if provenance.get('status') != 'verified' or len(provenance.get('outputs', {})) != 12:
        raise ValueError('Input preparation receipt is incomplete')
    sources = set()
    for base, rows in ((HERE / 'inputs', provenance['outputs']), (ROOT, provenance['sources'])):
        for name, expected in rows.items():
            path = base / name
            actual = fingerprint(path)
            if any(actual[key] != expected[key] for key in ('sha256', 'size_bytes')):
                raise ValueError(f'Prepared input/source changed: {path}')
            sources.add(path)
    actual = fingerprint(HERE / 'prepare_inputs.py')
    if any(actual[key] != provenance['preparation_source'][key] for key in ('sha256', 'size_bytes')):
        raise ValueError('Input preparer differs from preparation receipt')
    return sources


def build_info(path, role):
    document = json.loads(path.read_text())
    if document.get('status') != 'complete':
        raise ValueError(f'Build incomplete: {path}')
    verify(document['binary'])
    verify(document['snapshot'])
    for key in ('builder', 'builder_after', 'cargo', 'rustc'):
        verify(document[key])
    if document['builder'] != document['builder_after'] or document.get('returncode') != 0:
        raise ValueError('Builder changed or build did not succeed')
    command = [document['cargo']['path'], 'build', '--locked', '--offline', '--release', '--message-format=json-render-diagnostics', '-j', '2']
    if document['command'] != command or document['artifact']['profile'].get('opt_level') != '3' or document['artifact']['profile'].get('test') is not False:
        raise ValueError('Build command/profile differs from preregistered release build')
    for row in document['cargo_configs'].values():
        verify(row)
    snapshot = json.loads(Path(document['snapshot']['path']).read_text())
    verify(snapshot['exporter'])
    if snapshot['role'] != role or document['source_before'] != document['source_after'] or tree(Path(snapshot['source'])) != snapshot['files'] or snapshot['files'] != document['source_before']:
        raise ValueError('Build source binding differs')
    if role == 'baseline' and snapshot['provenance'].get('commit') != BASELINE:
        raise ValueError('Wrong baseline commit')
    if role == 'baseline':
        verify(snapshot['provenance']['archive'])
    if role == 'candidate':
        for name, row in snapshot['files'].items():
            current = fingerprint(ROOT / name)
            if any(current[key] != row[key] for key in ('sha256', 'size_bytes')):
                raise ValueError(f'Candidate compiled inputs differ from reviewed tree: {name}')
    return document, snapshot


def required_sources(snapshot):
    required = {PROTOCOL, Path(__file__), HERE / 'build_snapshot.py', HERE / 'test_run_comparison.py', HERE / 'prepare_inputs.py', HERE / 'restore_prefix.py', HERE / 'launch_bounded.py'}
    required.update(ROOT / name for name in snapshot['files'])
    required.update((ROOT / 'scripts/sharkmer_validate').rglob('*.py'))
    required.update(path for path in (HERE / 'inputs').rglob('*') if path.is_file())
    return required


def same_build_configuration(builds):
    for key in ('cargo', 'rustc', 'cargo_version', 'rustc_version', 'compiler_sysroot', 'comparison_environment', 'cargo_configs'):
        if builds['baseline'][key] != builds['candidate'][key]:
            raise ValueError(f'Build toolchain/configuration differs: {key}')


def verify_sample_plan(samples):
    original = json.loads((HERE / 'inputs/read_inputs.json').read_text())['samples']
    if set(samples) != {row['accession'] for row in original}:
        raise ValueError('Frozen sample set differs from prepared inputs')
    for row in original:
        sample = samples[row['accession']]
        if any(sample.get(key) != value for key, value in row.items() if key != 'path') or sample.get('observed_prefix') != row['expected_prefix']:
            raise ValueError('Frozen consumed identity differs from prepared inputs')
        if not Path(sample['path']).is_absolute():
            raise ValueError('Frozen read path must be absolute')


def commands(samples, builds, output, cpus):
    jobs = []
    for accession, roles in ORDER:
        sample = samples[accession]
        for role in roles:
            directory = output / accession / role
            binary = builds['baseline' if role == 'baseline_off' else 'candidate']['binary']['path']
            argv = [binary, '-k', '19', '-t', '2', '--chunks', '0', '--max-reads', '1000000', '-o', str(directory), '-s', f'uncertain_{accession}', '--pcr-panel-file', str(HERE / 'inputs/insecta.yaml')]
            if role == 'changed_on':
                argv.append('--emit-uncertain')
            argv.append(sample['path'])
            cpu_text = ','.join(map(str, cpus))
            command = [str(TOOLS['prlimit']), f"--as={LIMITS['address_space_bytes']}", '--', str(TOOLS['time']), '-v', '-o', str(directory / 'time.txt'), str(TOOLS['taskset']), '-c', cpu_text, str(Path(sys.executable).resolve()), str(HERE / 'launch_bounded.py'), '--receipt', str(directory / 'process_constraints.json'), '--expected-cpus', cpu_text, '--expected-as', str(LIMITS['address_space_bytes']), '--', *argv]
            jobs.append({'id': f'{accession}.{role}', 'accession': accession, 'role': role, 'output': str(directory), 'argv': argv, 'command': command})
    return jobs


def freeze(args):
    if args.output.exists() or args.run_dir.exists():
        raise ValueError('Freeze/output directory already exists')
    prepared = prepared_sources()
    manifest_path = HERE / 'inputs/read_inputs.json'
    manifest = json.loads(manifest_path.read_text())
    overrides = {}
    for item in args.read:
        accession, path = item.split('=', 1)
        if accession in overrides:
            raise ValueError('Duplicate read override')
        overrides[accession] = path
    samples = {row['accession']: dict(row) for row in manifest['samples']}
    if set(samples) != {row[0] for row in ORDER} or set(overrides) - set(samples):
        raise ValueError('Unexpected sample set')
    for accession, sample in samples.items():
        sample['path'] = str((ROOT / overrides.get(accession, sample['path'])).resolve())
        sample['observed_prefix'] = inspect_prefix(sample['path'], sample['expected_prefix'])
    builds, snapshots = {}, {}
    for role, path in [('baseline', args.baseline_build), ('candidate', args.candidate_build)]:
        builds[role], snapshots[role] = build_info(path.resolve(), role)
    same_build_configuration(builds)
    required = required_sources(snapshots['candidate'])
    reviewed = review_sources(args.review.resolve(), required)
    cpus = sorted(os.sched_getaffinity(0))[:2]
    if len(cpus) != 2:
        raise ValueError('Two allocated CPUs required')
    paths = required | prepared | {args.review.resolve(), args.baseline_build.resolve(), args.candidate_build.resolve(), Path(sys.executable).resolve(), *TOOLS.values()}
    paths.update(Path(row['path']) for row in reviewed.values())
    paths.update(path for path in (HERE / 'inputs').rglob('*') if path.is_file())
    for build, snapshot in zip(builds.values(), snapshots.values()):
        paths.update(Path(build[key]['path']) for key in ('binary', 'snapshot', 'cargo', 'rustc'))
        paths.update(Path(row['path']) for row in snapshot['files'].values())
        paths.update(Path(row['path']) for row in build['cargo_configs'].values())
        if 'archive' in snapshot['provenance']:
            paths.add(Path(snapshot['provenance']['archive']['path']))
    # Preserve every restoration attempt, including failed network attempts.
    paths.update((ROOT / 'tmp/uncertain-products-20260927/reads').glob('*.restore.json'))
    import yaml
    paths.update(path for path in Path(yaml.__file__).resolve().parent.rglob('*') if path.is_file() and path.suffix in ('.py', '.so'))
    frozen = {'schema_version': 1, 'status': 'frozen_before_execution', 'limits': LIMITS, 'cpus': cpus, 'samples': samples, 'builds': builds, 'review': fingerprint(args.review.resolve()), 'build_receipts': {'baseline': fingerprint(args.baseline_build.resolve()), 'candidate': fingerprint(args.candidate_build.resolve())}, 'reviewed_sources': reviewed, 'sources': {str(path.absolute()): fingerprint(path) for path in sorted(paths)}, 'output': str(args.run_dir.resolve()), 'created_unix': time.time(), 'freeze_argv': sys.argv, 'runtime': {'python': sys.version, 'interpreter': fingerprint(Path(sys.executable).resolve()), 'yaml_version': yaml.__version__}, 'execution_environment': {'LC_ALL': 'C'}}
    frozen['jobs'] = commands(samples, builds, args.run_dir.resolve(), cpus)
    write_new(args.output, frozen)
    print(json.dumps(fingerprint(args.output)))


def fasta(path):
    rows = []
    header = None
    sequence = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if header is not None:
                rows.append({'header': header, 'sequence': ''.join(sequence)})
            header, sequence = line[1:], []
        elif header is not None:
            sequence.append(line)
        elif line:
            raise ValueError('Sequence before FASTA header')
    if header is not None:
        rows.append({'header': header, 'sequence': ''.join(sequence)})
    if any(not row['sequence'] for row in rows):
        raise ValueError('Empty FASTA record')
    return rows


def capture(job, sample):
    sys.path.insert(0, str(ROOT / 'scripts'))
    from sharkmer_validate import runner
    directory = Path(job['output'])
    name = f"uncertain_{job['accession']}"
    stats = runner._parse_stats_yaml(directory / f'{name}.stats.yaml')
    panel = runner.load_panel_yaml(HERE / 'inputs/insecta.yaml')
    genes = {f"{panel['name']}_{gene}" for gene in runner.panel_gene_names(panel)}
    runner._validate_stats_manifest(stats, name, 19, genes, 'end-to-end', expected_command=job['argv'])
    runner._validate_output_manifest(stats, name, directory)
    if len(stats['pcr_results']) != len(genes):
        raise ValueError('Duplicate/surplus gene rows')
    expected = sample['expected_prefix']
    source = stats.get('input_source', {})
    if source != {'kind': 'local_files', 'inputs': [sample['path']], 'paired': False, 'max_reads': 1000000} or stats.get('chunks') != 0:
        raise ValueError('Unexpected input or chunks settings')
    if stats.get('n_reads_read') != expected['records'] or stats.get('n_bases_read') != expected['bases']:
        raise ValueError('Consumed input counts differ')
    if stats.get('stage_timings', {}).get('read_threading_input_s') is not None:
        raise ValueError('Read threading was unexpectedly evaluated')
    standard, uncertain = [], []
    for gene in stats['pcr_results']:
        assessment = gene.get('uncertain_products')
        if (assessment is not None) != (job['role'] == 'changed_on'):
            raise ValueError('Uncertain assessment presence differs from flag')
        if any('withheld_path_diagnostics' in item for item in gene.get('threshold_diagnostics', [])):
            raise ValueError('Hidden diagnostics unexpectedly enabled')
        ordinary = fasta(directory / gene['output_file']) if gene.get('output_file') else []
        if len(ordinary) != gene['n_products'] or [len(row['sequence']) for row in ordinary] != gene.get('product_lengths', []):
            raise ValueError('Ordinary FASTA count/length differs from gene stats')
        if gene['status'] == 'fail' and (ordinary or gene.get('output_file')):
            raise ValueError('Failed gene published ordinary output')
        if gene['n_products'] == 0 and gene.get('output_file'):
            raise ValueError('Zero-product gene declared ordinary output')
        standard.extend({'gene': gene['gene_name'], **row} for row in ordinary)
        if assessment and assessment.get('output_file'):
            uncertain.extend({'gene': gene['gene_name'], **row} for row in fasta(directory / assessment['output_file']))
    # Exercise the normal consumer's exclusion of uncertain records as well.
    products = runner.parse_fasta_products(name, directory)
    if sum(len(row['products']) for row in products) != len(standard):
        raise ValueError('Standard consumer count differs from explicit ordinary files')
    return {'stats': stats, 'standard': standard, 'uncertain': uncertain, 'files': {path.name: fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()}}


def standard_stats(stats):
    # Remove only the new optional gene assessment; retain every pre-existing field.
    return [{key: value for key, value in gene.items() if key != 'uncertain_products'} for gene in stats['pcr_results']]


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def analyze(records, losses):
    by_sample = {}
    for record in records:
        roles = by_sample.setdefault(record['accession'], {})
        if record['role'] in roles:
            raise ValueError('Duplicate job record')
        roles[record['role']] = record['outputs']
    if set(by_sample) != {row[0] for row in ORDER} or any(set(roles) != {'baseline_off', 'changed_off', 'changed_on'} for roles in by_sample.values()):
        raise ValueError('Expected all nine distinct runs')
    comparisons = []
    for accession, roles in by_sample.items():
        failures = []
        for left, right in [('baseline_off', 'changed_off'), ('changed_off', 'changed_on')]:
            a, b = roles[left], roles[right]
            for key in COUNTS:
                if key not in a['stats'] or key not in b['stats'] or a['stats'][key] != b['stats'][key]:
                    failures.append(f'{left}/{right}:{key}')
            for key in ('n_multi_kmers', 'n_singleton_kmers'):
                if (key in a['stats']) != (key in b['stats']) or a['stats'].get(key) != b['stats'].get(key):
                    failures.append(f'{left}/{right}:{key}')
            if a['standard'] != b['standard']:
                failures.append(f'{left}/{right}:ordered_standard_fasta')
            if standard_stats(a['stats']) != standard_stats(b['stats']):
                failures.append(f'{left}/{right}:gene_outcomes_and_threshold_diagnostics')
        comparisons.append({'accession': accession, 'passed': not failures, 'failures': failures, 'standard_counts': {role: len(value['standard']) for role, value in roles.items()}, 'uncertain_count': len(roles['changed_on']['uncertain'])})
    loss_rows = []
    for loss in losses:
        seq, gene = loss['sequence'], 'insecta_' + loss['gene']
        if hashlib.sha256(seq.encode()).hexdigest() != loss['sequence_sha256'] or len(seq) != loss['length']:
            raise ValueError('Historical sequence identity differs')
        roles = by_sample[loss['accession']]
        matches = {}
        for role, output in roles.items():
            matches[role] = {}
            for kind in ('standard', 'uncertain'):
                rows = [row for row in output[kind] if row['gene'] == gene]
                matches[role][kind] = {'exact_oriented': [row['header'] for row in rows if row['sequence'] == seq], 'reverse_complement_only': [row['header'] for row in rows if row['sequence'] == rc(seq) and row['sequence'] != seq]}
        entries = [entry for entry in roles['changed_on']['stats']['pcr_results'] if entry['gene_name'] == gene]
        if len(entries) != 1:
            raise ValueError('Historical target missing or duplicated')
        assessment = entries[0]['uncertain_products']
        loss_rows.append({key: value for key, value in loss.items() if key not in ('sequence', 'historical_header')} | {'matches': matches, 'assessment': assessment, 'non_observation_scope': 'Unchanged full-panel search and static 21-target run-budget shares (12–13 retained path records per target before other caps); bounded non-observation does not establish path or biological absence. Prior focused diagnostics used different per-target shares.'})
    return {'schema_version': 1, 'status': 'pass' if all(row['passed'] for row in comparisons) else 'fail', 'comparisons': comparisons, 'historical_losses': loss_rows, 'interpretation': 'Uncertain sequence availability is separate from standard recovery and does not validate repeat structure. Runtime/RSS are descriptive single-trio observations.'}


def check_frozen(frozen):
    import yaml
    if type(frozen.get('schema_version')) is not int or frozen['schema_version'] != 1 or frozen.get('status') != 'frozen_before_execution' or not frozen.get('sources') or not frozen.get('reviewed_sources'):
        raise ValueError('Frozen schema/review/source coverage is incomplete')
    if set(frozen['samples']) != {row[0] for row in ORDER} or set(frozen['builds']) != {'baseline', 'candidate'} or set(frozen['build_receipts']) != {'baseline', 'candidate'}:
        raise ValueError('Frozen sample/build roles differ')
    verify_sample_plan(frozen['samples'])
    if len(frozen['cpus']) != 2 or any(type(cpu) is not int for cpu in frozen['cpus']) or len(set(frozen['cpus'])) != 2 or not set(frozen['cpus']).issubset(os.sched_getaffinity(0)):
        raise ValueError('Frozen two-CPU affinity unavailable or invalid')
    if frozen['runtime'] != {'python': sys.version, 'interpreter': fingerprint(Path(sys.executable).resolve()), 'yaml_version': yaml.__version__} or frozen.get('execution_environment') != {'LC_ALL': 'C'}:
        raise ValueError('Execution interpreter differs from frozen runtime')
    for name, row in frozen['sources'].items():
        if name != row['path']:
            raise ValueError('Frozen source key differs from its path')
        verify(row)
    for role, receipt in frozen['build_receipts'].items():
        verify(receipt)
        build, snapshot = build_info(Path(receipt['path']), role)
        if build != frozen['builds'][role]:
            raise ValueError('Build receipt differs from frozen build')
        if role == 'candidate':
            required = required_sources(snapshot)
    same_build_configuration(frozen['builds'])
    verify(frozen['review'])
    reviewed = review_sources(Path(frozen['review']['path']), required)
    if reviewed != frozen['reviewed_sources']:
        raise ValueError('Independent review differs from frozen review')
    mandatory = {str(path.absolute()) for path in required} | {row['path'] for row in reviewed.values()} | {row['path'] for row in frozen['build_receipts'].values()} | {frozen['review']['path'], str(Path(sys.executable).resolve()), *map(str, TOOLS.values())}
    if mandatory - set(frozen['sources']):
        raise ValueError('Mandatory frozen source identities missing')
    if frozen['limits'] != LIMITS or frozen['jobs'] != commands(frozen['samples'], frozen['builds'], Path(frozen['output']), frozen['cpus']):
        raise ValueError('Frozen command plan differs from protocol')


def run_job(job, sample, frozen):
    directory = Path(job['output'])
    directory.mkdir(parents=True, exist_ok=False)
    receipt = dict(job, status='failed', affinity_before=sorted(os.sched_getaffinity(0)))
    try:
        check_frozen(frozen)
        receipt['prefix_before'] = inspect_prefix(sample['path'], sample['expected_prefix'])
        # Prefix verification immediately above also prewarms every consumed byte.
        started = time.monotonic()
        receipt['started_unix'] = time.time()
        with (directory / 'stdout.txt').open('x') as stdout, (directory / 'stderr.txt').open('x') as stderr:
            process = subprocess.Popen(job['command'], stdout=stdout, stderr=stderr, start_new_session=True, env=dict(os.environ, LC_ALL='C'))
            try:
                receipt['returncode'] = process.wait(timeout=LIMITS['timeout_seconds'])
                receipt['timed_out'] = False
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                receipt.update(returncode=process.returncode, timed_out=True)
        receipt['wall_seconds'] = time.monotonic() - started
        receipt['prefix_after'] = inspect_prefix(sample['path'], sample['expected_prefix'])
        check_frozen(frozen)
        if receipt['returncode'] != 0 or receipt['timed_out']:
            raise ValueError('Command failed or timed out')
        enforced = json.loads((directory / 'process_constraints.json').read_text())
        if enforced.get('status') != 'pass' or enforced.get('command') != job['argv'] or enforced.get('affinity') != frozen['cpus'] or enforced.get('rlimit_as_soft_hard') != [LIMITS['address_space_bytes']] * 2:
            raise ValueError('Actual process constraints differ from frozen command')
        receipt['outputs'] = capture(job, sample)
        match = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', (directory / 'time.txt').read_text())
        if not match:
            raise ValueError('RSS measurement missing')
        receipt['maximum_rss_bytes'] = int(match[1]) * 1024
        receipt['enforcement'] = {'observed_before_exec': enforced, 'timeout_seconds': LIMITS['timeout_seconds'], 'sequential': True}
        receipt['status'] = 'complete'
    except Exception as error:
        receipt['error'] = repr(error)
    return receipt


def execute(args):
    before = fingerprint(args.freeze)
    if before['sha256'] != args.freeze_sha256:
        raise ValueError('Externally pinned freeze hash differs')
    frozen = json.loads(args.freeze.read_text())
    output = Path(frozen['output'])
    output.mkdir(parents=True, exist_ok=False)
    execution = {'schema_version': 1, 'status': 'failed', 'freeze_before': before, 'argv': sys.argv, 'jobs': []}
    try:
        check_frozen(frozen)
        for job in frozen['jobs']:
            receipt = run_job(job, frozen['samples'][job['accession']], frozen)
            receipt_path = output / (job['id'] + '.receipt.json')
            write_new(receipt_path, receipt)
            execution['jobs'].append({'id': job['id'], 'receipt': fingerprint(receipt_path), 'status': receipt['status']})
            if receipt['status'] != 'complete':
                raise ValueError(f"Job failed: {job['id']}")
        records = [json.loads(Path(item['receipt']['path']).read_text()) for item in execution['jobs']]
        result = analyze(records, json.loads((HERE / 'inputs/historical_losses.json').read_text())['rows'])
        for record in records:
            for row in record['outputs']['files'].values():
                verify(row)
        for item in execution['jobs']:
            verify(item['receipt'])
        for sample in frozen['samples'].values():
            inspect_prefix(sample['path'], sample['expected_prefix'])
        check_frozen(frozen)
        execution['freeze_after'] = fingerprint(args.freeze)
        if execution['freeze_after'] != before:
            raise ValueError('Freeze changed during execution')
        write_new(output / 'assessment.json', result)
        execution['assessment'] = fingerprint(output / 'assessment.json')
        if result['status'] != 'pass':
            raise ValueError('Standard parity failed')
        execution['status'] = 'complete'
    except Exception as error:
        execution['error'] = repr(error)
    finally:
        write_new(output / 'execution.json', execution)
    return 0 if execution['status'] == 'complete' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    plan = sub.add_parser('freeze')
    plan.add_argument('--baseline-build', type=Path, required=True)
    plan.add_argument('--candidate-build', type=Path, required=True)
    plan.add_argument('--review', type=Path, required=True)
    plan.add_argument('--read', action='append', default=[], help='ACCESSION=path override; frozen identities remain unchanged')
    plan.add_argument('--run-dir', type=Path, required=True)
    plan.add_argument('--output', type=Path, required=True)
    run = sub.add_parser('execute')
    run.add_argument('--freeze', type=Path, required=True)
    run.add_argument('--freeze-sha256', required=True)
    args = parser.parse_args()
    if args.action == 'freeze':
        freeze(args)
        return 0
    return execute(args)


if __name__ == '__main__':
    raise SystemExit(main())

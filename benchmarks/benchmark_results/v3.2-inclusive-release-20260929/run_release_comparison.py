#!/usr/bin/env python3
"""Run the frozen ten-sample release comparison with class-aware outputs."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OLD = ROOT / 'benchmarks/benchmark_results/v3.2-uncertain-products-20260927'
LEGACY = ROOT / 'benchmarks/benchmark_results/v3.2-release-comparison-20260909/tools'
sys.path.insert(0, str(OLD))
from build_snapshot import fingerprint, write_new, tree, cargo_configs
from run_comparison import inspect_prefix, fasta, verify, COUNTS
from run_compiled_checks import runtime
sys.path.insert(0, str(LEGACY))
from postprocess_legacy import parse_legacy_fasta
sys.path.insert(0, str(ROOT / 'scripts'))
from sharkmer_validate import runner
from build_pair import COMMITS

ADDRESS_SPACE = 40 * 1024**3
TIMEOUT = 1800


def require(value, message):
    if not value:
        raise ValueError(message)


def status_of(header):
    statuses = [token.split('=', 1)[1] for token in shlex.split(header) if token.startswith('assembly_status=')]
    require(len(statuses) == 1 and statuses[0] in ('standard', 'uncertain'), 'Missing, duplicate or unknown assembly status')
    return statuses[0]


def sequence_record(row, gene, status, scope):
    sequence = row['sequence']
    require(bool(sequence) and set(sequence) <= set('ACGT'), 'Invalid output sequence')
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return {'gene': gene, 'assembly_status': status, 'scope': scope, 'header': row['header'],
            'sequence': sequence, 'length': len(sequence),
            'sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'reverse_complement_sha256': hashlib.sha256(reverse.encode()).hexdigest()}


def capture(directory, invocation, sample, command, scopes):
    name = invocation['sample_prefix']
    stats_path = directory / f'{name}.stats.yaml'
    stats = runner._parse_stats_yaml(stats_path)
    expected = sample['expected_prefix']
    genes = set(scopes)
    require(stats.get('sample') == name and stats.get('kmer_length') == 19 and stats.get('chunks') == 0, 'Stats identity/settings differ')
    require(stats.get('command') == ' '.join(command), 'Stats command differs')
    require(stats.get('sharkmer_version') == ('3.1.0' if invocation['version'] == 'baseline' else '3.2.0-dev'), 'Wrong binary version')
    require(stats.get('n_reads_read') == expected['records'] and stats.get('n_bases_read') == expected['bases'], 'Consumed records/bases differ')
    for key in (*COUNTS, 'peak_memory_bytes'):
        require(type(stats.get(key)) is int and stats[key] >= 0, f'Invalid {key}')
    entries = stats.get('pcr_results', [])
    require(len(entries) == len(genes) and {entry.get('gene_name') for entry in entries} == genes, 'Wrong/duplicate target rows')
    manifest = None
    if invocation['version'] == 'candidate':
        runner._validate_stats_manifest(stats, name, 19, genes, 'end-to-end', expected_command=command)
        manifest = runner._validate_output_manifest(stats, name, directory)
        require(manifest['schema_version'] == 2 and stats.get('product_output_schema_version') == 2
                and stats.get('product_output_layout') == 'combined_classified', 'Wrong candidate schema')
        require(stats['input_source'] == {'kind': 'local_files', 'inputs': [sample['available_path']],
                                         'paired': False, 'max_reads': 1000000}, 'Wrong input provenance')
        require(stats['stage_timings']['read_threading_input_s'] is None, 'Threading unexpectedly enabled')
    else:
        require(not any(key in stats for key in ('run_id', 'run_status', 'output_manifest', 'input_source', 'stage_timings')), 'Legacy schema unexpectedly changed')
    records, expected_files = [], {stats_path.name}
    for entry in entries:
        gene = entry['gene_name']
        count = entry['n_products']
        require(type(count) is int and count >= 0, 'Invalid product count')
        path = directory / f'{name}_{gene}.fasta'
        if invocation['version'] == 'baseline':
            require(entry['status'] == ('success' if count else 'fail'), 'Legacy count/status mismatch')
            rows = parse_legacy_fasta(path, name, gene) if count else []
            classified = [('standard', row) for row in rows]
        else:
            require(bool(entry.get('output_file')) == bool(count), 'Candidate output locator/count differs')
            if count:
                require(entry['output_file'] == path.name, 'Wrong target filename')
            rows = fasta(path) if count else []
            classified = [(status_of(row['header']), row) for row in rows]
            labels = [label for label, row in classified]
            require(labels == sorted(labels, key=lambda label: label != 'standard'), 'Mixed class ordering differs')
            require(entry['n_standard_products'] == labels.count('standard') and entry['n_uncertain_products'] == labels.count('uncertain'), 'Class subtotals differ')
        require(count == len(rows) and entry.get('product_lengths', []) == [len(row['sequence']) for row in rows], 'FASTA counts/lengths differ')
        if rows:
            expected_files.add(path.name)
        records.extend(sequence_record(row, gene, label, scopes[gene]) for label, row in classified)
    require({path.name for path in directory.glob('*.fasta')} == {name for name in expected_files if name.endswith('.fasta')}, 'Unexpected/missing FASTA')
    if invocation['version'] == 'baseline':
        require({path.name for path in directory.iterdir()} == expected_files, 'Unexpected legacy file')
    else:
        consumed = runner.parse_fasta_products(name, directory)
        selected = {row['gene']: [product['sequence'] for product in row['products']] for row in consumed}
        expected_selection = {gene: [row['sequence'] for row in records if row['gene'] == gene and row['assembly_status'] == 'standard'] for gene in genes}
        require({gene: seqs for gene, seqs in selected.items() if seqs} == {gene: seqs for gene, seqs in expected_selection.items() if seqs}, 'Standard-only consumer differs')
        require(runner._validate_output_manifest(stats, name, directory) == manifest, 'Manifest changed during parsing')
    return {'records': records, 'counts': {key: stats[key] for key in COUNTS},
            'allocator_peak_bytes': stats['peak_memory_bytes'], 'pcr_results': entries,
            'files': {path.name: fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()}}


def physical_cpus():
    selected, seen, topology = [], set(), []
    for cpu in sorted(os.sched_getaffinity(0)):
        directory = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        pair = tuple(int((directory / field).read_text()) for field in ('physical_package_id', 'core_id'))
        topology.append({'cpu': cpu, 'socket': pair[0], 'core': pair[1]})
        if pair not in seen and len(selected) < 2:
            seen.add(pair)
            selected.append(cpu)
    require(len(selected) == 2, 'Need two distinct physical cores')
    os.sched_setaffinity(0, selected)
    return selected, topology


def verify_builds(path, bindings):
    document = json.loads(path.read_text())
    require(document.get('status') == 'pass' and set(document['versions']) == set(COMMITS), 'Incomplete paired builds')
    bindings[str(path)] = fingerprint(path)
    builds = {}
    for key in ('wrapper', 'builder'):
        bindings[document[key]['path']] = document[key]
    for role, version in document['versions'].items():
        require(version['commit'] == COMMITS[role], 'Wrong source commit')
        for key in ('snapshot', 'build', 'binary'):
            verify(version[key])
            bindings[version[key]['path']] = version[key]
        snapshot = json.loads(Path(version['snapshot']['path']).read_text())
        build = json.loads(Path(version['build']['path']).read_text())
        require(build['status'] == 'complete' and snapshot['provenance']['commit'] == version['commit'], 'Build source identity differs')
        require(snapshot['files'] == build['source_before'] == build['source_after'] == tree(Path(snapshot['source'])), 'Build source changed')
        require(build['binary'] == version['binary'] and build['returncode'] == 0, 'Executable differs from completed build')
        require(set(build['artifact']['features']) == {'default', 'ahashmap'}, 'Wrong hash backend')
        require(build['command'] == [build['cargo']['path'], 'build', '--locked', '--offline', '--release', '--message-format=json-render-diagnostics', '-j', '2'], 'Wrong build flags')
        require(build['cargo_configs'] == cargo_configs(Path(snapshot['source']), os.environ), 'Cargo configuration changed')
        for row in [snapshot['provenance']['archive'], snapshot['exporter'], *snapshot['files'].values(), build['cargo'], build['rustc'], *build['cargo_configs'].values()]:
            bindings[row['path']] = row
        builds[role] = build
    for key in ('cargo', 'rustc', 'cargo_version', 'rustc_version', 'comparison_environment'):
        require(builds['baseline'][key] == builds['candidate'][key], f'Build configurations differ: {key}')
    return document['versions']


def run_timed(invocation, cpus, environment):
    directory = Path(invocation['directory'])
    directory.mkdir(parents=True, exist_ok=False)
    Path(invocation['outputs']).mkdir()
    timing = directory / 'time.json'
    constraints = directory / 'process_constraints.json'
    timed = ['/usr/bin/time', '-o', str(timing), '-f', '{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}', '--', *invocation['command']]
    command = ['/usr/bin/prlimit', f'--as={ADDRESS_SPACE}', '--', str(Path(sys.executable).resolve()),
               str(OLD / 'launch_bounded.py'), '--receipt', str(constraints), '--expected-cpus', ','.join(map(str, cpus)),
               '--expected-as', str(ADDRESS_SPACE), '--', *timed]
    result = {'status': 'failed', 'command': command, 'timing_boundary': 'GNU time directly measures Sharkmer; constraint observer executes before timing'}
    start = time.monotonic()
    with (directory / 'stdout.log').open('x') as stdout, (directory / 'stderr.log').open('x') as stderr:
        process = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            result['returncode'] = process.wait(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            result.update(timed_out=True, returncode=process.returncode)
    result['observer_wall_seconds'] = time.monotonic() - start
    result['files'] = {path.name: fingerprint(path) for path in sorted(directory.iterdir()) if path.is_file()}
    if result['returncode'] != 0:
        return result
    observed = json.loads(constraints.read_text())
    require(observed['status'] == 'pass' and observed['command'] == timed and observed['affinity'] == cpus
            and observed['rlimit_as_soft_hard'] == [ADDRESS_SPACE] * 2, 'Actual process constraints differ')
    measured = json.loads(timing.read_text())
    require(measured['exit_code'] == 0 and all(math.isfinite(measured[key]) and measured[key] >= 0 for key in ('wall_seconds', 'user_seconds', 'system_seconds', 'max_rss_kib')), 'Invalid GNU time result')
    measured['peak_rss_bytes'] = measured['max_rss_kib'] * 1024
    result.update(status='pass', metrics=measured)
    return result


def unique(records):
    return {(row['gene'], row['sha256']): row for row in records}


def historical_identity_check(records, expectation, role):
    """Require release reproduction; expose historical candidate differences."""
    expected = {expectation['panel'] + '_' + gene: [row['sequence'] for row in rows]
                for gene, rows in expectation['roles'][role]['by_gene'].items()}
    observed = {gene: [] for gene in expected}
    for row in records:
        if row['assembly_status'] == 'standard':
            observed.setdefault(row['gene'], []).append(row['sequence'])
    same = expected == observed
    if role == 'baseline':
        require(same, f"Released sequence reproduction differs for {expectation['accession']}")
    old = {(gene, hashlib.sha256(sequence.encode()).hexdigest()) for gene, sequences in expected.items() for sequence in sequences}
    new = {(gene, hashlib.sha256(sequence.encode()).hexdigest()) for gene, sequences in observed.items() for sequence in sequences}
    return {'ordered_standard_sequences_identical': same,
            'historical_role': 'released_v3.1.0' if role == 'baseline' else 'september_11_pre_inclusive_candidate',
            'missing_literal_identities': [{'gene': gene, 'sha256': digest} for gene, digest in sorted(old - new)],
            'gained_literal_identities': [{'gene': gene, 'sha256': digest} for gene, digest in sorted(new - old)]}


def dispositions(baseline, candidate):
    rows = []
    for key, original in sorted(unique(baseline).items()):
        evidence = {}
        for label in ('standard', 'uncertain'):
            same_gene = [row for row in candidate if row['gene'] == original['gene'] and row['assembly_status'] == label]
            for orientation, digest in [('exact', original['sha256']), ('rc', original['reverse_complement_sha256'])]:
                evidence[f'{label}_{orientation}'] = [row['sha256'] for row in same_gene if row['sha256'] == digest]
        disposition = next((label for label in ('standard_exact', 'standard_rc', 'uncertain_exact', 'uncertain_rc') if evidence[label]), 'unobserved')
        rows.append({'gene': key[0], 'sha256': key[1], 'length': original['length'], 'scope': original['scope'],
                     'disposition': disposition, 'evidence': evidence})
    return rows


def summarize(jobs):
    require(len(jobs) == 60 and all(job['status'] == 'pass' for job in jobs), 'Need all sixty completed invocations')
    result = {'status': 'pass', 'cells': [], 'biological_support': 'unassessed in this identity/resource comparison'}
    for accession in dict.fromkeys(job['invocation']['input'] for job in jobs):
        cell = {'accession': accession, 'roles': {}, 'paired_resources': []}
        role_jobs = {}
        for role in ('baseline', 'candidate'):
            selected = [job for job in jobs if job['invocation']['input'] == accession and job['invocation']['version'] == role]
            require(len(selected) == 3 and {job['invocation']['pair_index'] for job in selected} == {1, 2, 3}, 'Wrong repetitions')
            for key in ('records', 'counts', 'pcr_results'):
                require(all(job['outputs'][key] == selected[0]['outputs'][key] for job in selected), f'{accession}/{role}: unstable {key}')
            records = selected[0]['outputs']['records']
            counts = {label: len([row for row in records if row['assembly_status'] == label]) for label in ('standard', 'uncertain')}
            scopes = {}
            for scope in sorted({row['scope'] for row in records}):
                scopes[scope] = {label: len(unique([row for row in records if row['scope'] == scope and row['assembly_status'] == label])) for label in ('standard', 'uncertain')}
            cell['roles'][role] = {'record_counts': counts, 'unique_counts': {label: len(unique([row for row in records if row['assembly_status'] == label])) for label in counts},
                                   'scope_unique_counts': scopes, 'standard_targets': sorted({row['gene'] for row in records if row['assembly_status'] == 'standard'}),
                                   'available_targets': sorted({row['gene'] for row in records}), 'metrics': {}}
            for metric in ('wall_seconds', 'peak_rss_bytes'):
                values = [job['execution']['metrics'][metric] for job in sorted(selected, key=lambda job: job['invocation']['pair_index'])]
                cell['roles'][role]['metrics'][metric] = {'values': values, 'median': statistics.median(values), 'min': min(values), 'max': max(values)}
            role_jobs[role] = {job['invocation']['pair_index']: job for job in selected}
        require(role_jobs['baseline'][1]['outputs']['counts'] == role_jobs['candidate'][1]['outputs']['counts'], f'{accession}: aggregate count mismatch')
        old = role_jobs['baseline'][1]['outputs']['records']
        new = role_jobs['candidate'][1]['outputs']['records']
        cell['released_identity_dispositions'] = dispositions(old, new)
        old_keys = unique(old)
        cell['literal_gains'] = [{'gene': row['gene'], 'sha256': row['sha256'], 'length': row['length'], 'assembly_status': row['assembly_status'], 'scope': row['scope'],
                                 'reverse_complement_in_baseline': (row['gene'], row['reverse_complement_sha256']) in old_keys}
                                for row in {(row['gene'], row['sha256'], row['assembly_status']): row for row in new}.values()
                                if (row['gene'], row['sha256']) not in old_keys]
        cell['resource_changes'] = {}
        for metric, threshold in [('wall_seconds', 5), ('peak_rss_bytes', 2)]:
            before, after = [cell['roles'][role]['metrics'][metric]['median'] for role in ('baseline', 'candidate')]
            percent = 100 * (after / before - 1) if before else None
            cell['resource_changes'][metric] = {'absolute_delta': after - before, 'percent': percent,
                                              'increase': after > before, 'priority_review': percent is not None and percent > threshold}
        for pair in (1, 2, 3):
            paired = {'pair': pair}
            for metric in ('wall_seconds', 'peak_rss_bytes'):
                before, after = [role_jobs[role][pair]['execution']['metrics'][metric] for role in ('baseline', 'candidate')]
                paired[metric] = {'baseline': before, 'candidate': after, 'absolute_delta': after - before,
                                  'ratio': after / before if before else None}
            cell['paired_resources'].append(paired)
        result['cells'].append(cell)
    require(len(result['cells']) == 10, 'Wrong sample count')
    result['totals'] = {role: {label: sum(cell['roles'][role]['record_counts'][label] for cell in result['cells']) for label in ('standard', 'uncertain')} for role in ('baseline', 'candidate')}
    result['median_wall_sum_seconds'] = {role: sum(cell['roles'][role]['metrics']['wall_seconds']['median'] for cell in result['cells']) for role in ('baseline', 'candidate')}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--restored', type=Path, required=True)
    parser.add_argument('--builds', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    result = {'status': 'failed', 'jobs': [], 'bindings': {}, 'slurm_job_id': os.environ.get('SLURM_JOB_ID')}
    try:
        cpus, topology = physical_cpus()
        result.update(cpus=cpus, allowed_topology=topology, machine=platform.uname()._asdict(), runtime=runtime())
        bindings = result['bindings']
        paths = [Path(__file__), HERE / 'build_pair.py', HERE / 'build_pair.sbatch', HERE / 'run_release_comparison.sbatch', HERE / 'test_release_comparison.py',
                 ROOT / 'dev_docs/INCLUSIVE_RELEASE_PROTOCOL.md', OLD / 'build_snapshot.py', OLD / 'run_comparison.py', OLD / 'run_compiled_checks.py',
                 OLD / 'launch_bounded.py', LEGACY / 'postprocess_legacy.py', args.inventory.resolve(), args.restored.resolve(),
                 Path('/usr/bin/time'), Path('/usr/bin/prlimit'), ROOT / 'benchmarks/benchmark_results/v3.2-high-copy-final-20260911/provenance.json',
                 ROOT / 'tests/test_uncertain_output.py', ROOT / 'scripts/bootstrap_from_runs.py', HERE / 'inputs/historical_identities.json']
        paths.extend((HERE / 'inputs').glob('*.py'))
        paths.extend((HERE / 'inputs').glob('*.sbatch'))
        paths.extend((ROOT / 'scripts/sharkmer_validate').glob('*.py'))
        for path in paths:
            bindings[str(path.resolve())] = fingerprint(path.resolve())
        builds = verify_builds(args.builds.resolve(), bindings)
        inventory = json.loads(args.inventory.read_text())
        restored = json.loads(args.restored.read_text())
        require(restored['status'] == 'pass' and restored['complete_inventory_selected'], 'Inputs are not restored/verified')
        samples = {row['accession']: row for row in restored['samples']}
        require(set(samples) == {row['accession'] for row in inventory['samples']} and len(samples) == 10, 'Wrong restored inputs')
        for row in inventory['samples']:
            sample = samples[row['accession']]
            require(sample['status'] == 'verified' and sample['expected_prefix'] == row['expected_prefix'] and sample['panel_id'] == row['panel'], 'Restored input differs')
            verify(sample['receipt'])
            bindings[sample['receipt']['path']] = sample['receipt']
        historical = json.loads((ROOT / 'benchmarks/benchmark_results/v3.2-high-copy-final-20260911/provenance.json').read_text())
        require(inventory['historical_schedule'] == historical['schedule'] and inventory['historical_protocol'] == historical['protocol'], 'Historical scope changed')
        oracle = json.loads((HERE / 'inputs/historical_identities.json').read_text())
        require(oracle['status'] == 'pass' and oracle['baseline_commit'] == COMMITS['baseline'], 'Historical identity oracle is incomplete')
        expectations = {row['accession']: row for row in oracle['sample_expectations']}
        require(set(expectations) == set(samples) and len(oracle['sample_expectations']) == 10, 'Historical oracle sample set differs')
        for row in oracle['sources_before'].values():
            verify(row)
            bindings[row['path']] = row
        require(oracle['sources_before'] == oracle['sources_after'], 'Historical identity derivation source changed')
        for expectation in expectations.values():
            for role in ('baseline', 'candidate'):
                for products in expectation['roles'][role]['by_gene'].values():
                    for product in products:
                        require(product['sequence_sha256'] == hashlib.sha256(product['sequence'].encode()).hexdigest()
                                and product['length'] == len(product['sequence']), 'Historical oracle sequence identity differs')
        panels = {row['name']: row['assembly_panel'] for row in inventory['panels']}
        for row in panels.values():
            verify(row)
            bindings[row['path']] = row
        scopes = {panel: {panel + '_' + row['gene']: row['scope'] for row in rows} for panel, rows in historical['protocol']['target_metadata'].items()}
        plan = []
        for original in inventory['historical_schedule']:
            invocation = dict(original)
            invocation['sample_prefix'] = f"{original['panel']}_{original['input']}_{original['depth']}"
            directory = output / 'runs' / original['invocation_id']
            invocation.update(directory=str(directory), outputs=str(directory / 'outputs'))
            invocation['command'] = [builds[original['version']]['binary']['path'], '-k', '19', '-t', '2', '--chunks', '0', '--max-reads', '1000000',
                                     '-o', invocation['outputs'] + '/', '-s', invocation['sample_prefix'], '--pcr-panel-file', panels[original['panel']]['path'], samples[original['input']]['available_path']]
            plan.append(invocation)
        require(len(plan) == 60 and len({row['invocation_id'] for row in plan}) == 60, 'Invalid frozen schedule')
        for row in bindings.values():
            verify(row)
        write_new(output / 'plan.json', {'schedule': plan, 'bindings': bindings, 'cpus': cpus, 'topology': topology, 'runtime': result['runtime']})
        result['plan'] = fingerprint(output / 'plan.json')
        environment = dict(os.environ, LC_ALL='C')
        for invocation in plan:
            sample = samples[invocation['input']]
            job = {'status': 'failed', 'invocation': invocation}
            result['jobs'].append(job)
            try:
                prefix = Path(sample['available_path'])
                job['prefix_before'] = inspect_prefix(prefix, sample['expected_prefix'])
                job['input_file'] = fingerprint(prefix)
                verify(builds[invocation['version']]['binary'])
                verify(panels[invocation['panel']])
                job['execution'] = run_timed(invocation, cpus, environment)
                require(job['execution']['status'] == 'pass', f"Invocation failed: {invocation['invocation_id']}")
                job['outputs'] = capture(Path(invocation['outputs']), invocation, sample, invocation['command'], scopes[invocation['panel']])
                job['historical_identity_check'] = historical_identity_check(job['outputs']['records'], expectations[invocation['input']], invocation['version'])
                verify(job['input_file'])
                verify(builds[invocation['version']]['binary'])
                verify(panels[invocation['panel']])
                job['status'] = 'pass'
            finally:
                path = output / f"{invocation['invocation_id']}.json"
                write_new(path, job)
                job['receipt'] = fingerprint(path)
            print(json.dumps({'completed': len(result['jobs']), 'invocation': invocation['invocation_id'], 'status': job['status']}), flush=True)
        summary = summarize(result['jobs'])
        write_new(output / 'summary.json', summary)
        result['summary'] = fingerprint(output / 'summary.json')
        result['status'] = 'pass'
    except Exception as error:
        result['error'] = repr(error)
    finally:
        try:
            for row in result['bindings'].values():
                verify(row)
            if 'plan' in result:
                verify(result['plan'])
            if 'summary' in result:
                verify(result['summary'])
            for job in result['jobs']:
                for key in ('input_file', 'receipt'):
                    if key in job:
                        verify(job[key])
                for row in job.get('outputs', {}).get('files', {}).values():
                    verify(row)
                for row in job.get('execution', {}).get('files', {}).values():
                    verify(row)
            require(runtime() == result['runtime'], 'Runtime changed during comparison')
        except Exception as error:
            result.update(status='failed', bookend_error=repr(error))
        write_new(output / 'execution.json', result)
    print(json.dumps({'status': result['status'], 'error': result.get('error'), 'output': str(output)}))
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

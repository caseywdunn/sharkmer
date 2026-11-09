#!/usr/bin/env python3
"""Run reviewed plant panel snapshots using one frozen binary and bounded jobs."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SHARED = ROOT / 'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/run_release_comparison.py'
sys.path.insert(0, str(SHARED.parent))
spec = importlib.util.spec_from_file_location('strict_capture', SHARED)
capture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(capture)
CHANGED_GENE = 'angiospermae_trnV-atpE'


def dump_new(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True) + '\n')


def require(value, message):
    if not value:
        raise ValueError(message)


def identities(records, gene):
    return Counter((r['sha256'], r['assembly_status']) for r in records if r['gene'] == gene)


def identity_rows(counter):
    return [{'sha256': key[0], 'assembly_status': key[1], 'count': count} for key, count in sorted(counter.items())]


def target_diagnostics(entry):
    """Compare outcomes without run-specific output filenames or descriptions."""
    outcome_keys = ('status', 'n_products', 'n_standard_products', 'n_uncertain_products',
                    'failure_reason', 'threshold_diagnostics')
    projection = {key: entry.get(key) for key in outcome_keys}
    flags = []
    flag_keys = {'search_limits_reached', 'retention_truncated', 'marker_detail_truncated', 'marker_detail',
                 'path_search_evaluated', 'search_limit_scope'}

    def walk(value, path):
        if isinstance(value, dict):
            for key in sorted(value):
                if key in flag_keys:
                    flags.append({'path': path + [key], 'value': value[key]})
                walk(value[key], path + [key])
        elif isinstance(value, list):
            for index, item in enumerate(value):
                walk(item, path + [index])

    walk(entry.get('uncertain_products'), ['uncertain_products'])
    projection['uncertainty_limit_flags'] = flags
    return projection


def assess(jobs, genes, reference_rows):
    reference_by_sha = {}
    for row in reference_rows:
        reference_by_sha.setdefault(row['provenance']['sequence_sha256'], []).append(row['accession'])
    comparisons = []
    pairs = {(j['scope'], j['sample_id']) for j in jobs}
    for scope, sample_id in sorted(pairs):
        pair = {j['variant']: j for j in jobs if (j['scope'], j['sample_id']) == (scope, sample_id)}
        row = {'scope': scope, 'sample_id': sample_id, 'targets': []}
        if set(pair) != {'before', 'after'} or any(j['status'] != 'pass' for j in pair.values()):
            row['status'] = 'incomplete_execution'
            comparisons.append(row)
            continue
        counts = {variant: {key: pair[variant]['outputs']['counts'][key] for key in capture.COUNTS}
                  for variant in ('before', 'after')}
        row['aggregate_counts'] = counts
        row['aggregate_counts_unchanged'] = counts['before'] == counts['after']
        diagnostics = {variant: {entry['gene_name']: target_diagnostics(entry)
                                  for entry in pair[variant]['outputs']['pcr_results']}
                       for variant in ('before', 'after')}
        require(all(set(diagnostics[variant]) == set(genes) for variant in diagnostics),
                'Captured target diagnostics missing or unexpected')
        for gene in sorted(genes):
            old, new = (identities(pair[variant]['outputs']['records'], gene) for variant in ('before', 'after'))
            row['targets'].append({'gene': gene, 'assay_changed': gene == CHANGED_GENE,
                                   'before': identity_rows(old), 'after': identity_rows(new),
                                   'lost_identities': identity_rows(old - new), 'gained_identities': identity_rows(new - old),
                                   'literal_identity_status_multiplicity_unchanged': old == new,
                                   'diagnostics_before': diagnostics['before'][gene],
                                   'diagnostics_after': diagnostics['after'][gene],
                                   'outcome_and_search_diagnostics_unchanged': diagnostics['before'][gene] == diagnostics['after'][gene]})
        controls = [r for r in row['targets'] if not r['assay_changed']]
        row['other_nine_literal_outputs_unchanged'] = all(r['literal_identity_status_multiplicity_unchanged'] for r in controls)
        row['other_nine_diagnostics_unchanged'] = all(r['outcome_and_search_diagnostics_unchanged'] for r in controls)
        row['other_nine_targets_unchanged'] = row['other_nine_literal_outputs_unchanged'] and row['other_nine_diagnostics_unchanged']
        row['partial_reference_assessment'] = 'pending_post_run'
        row['target_reference_matches'] = []
        for variant, job in pair.items():
            for record in job['outputs']['records']:
                if record['gene'] == CHANGED_GENE:
                    row['target_reference_matches'].append({'variant': variant, 'sha256': record['sha256'],
                                                           'length': record['length'], 'assembly_status': record['assembly_status'],
                                                           'exact_full_reference_accessions': reference_by_sha.get(record['sha256'], []),
                                                           'sample_full_phase': 'not_established'})
        if scope == 'synthetic':
            expected_sha = pair['after']['expected_corrected_sha256']
            row['expected_corrected_sha256'] = expected_sha
            row['exact_corrected_source_emitted_before'] = any(r['gene'] == CHANGED_GENE and r['sha256'] == expected_sha for r in pair['before']['outputs']['records'])
            row['exact_corrected_source_emitted_after'] = any(r['gene'] == CHANGED_GENE and r['sha256'] == expected_sha for r in pair['after']['outputs']['records'])
            if not row['exact_corrected_source_emitted_after']:
                row['status'] = 'source_reconstruction_not_observed'
            else:
                row['status'] = 'pass' if row['other_nine_targets_unchanged'] else 'control_target_changes_require_review'
        else:
            row['status'] = 'pass' if row['other_nine_targets_unchanged'] else 'control_target_changes_require_review'
        if not row['aggregate_counts_unchanged']:
            row['status'] = 'aggregate_counts_require_review'
        comparisons.append(row)
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--scope', choices=('synthetic', 'real', 'all'), default='all')
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    require(not args.output.exists(), 'Output must be a fresh directory')
    args.output = args.output.resolve()
    args.binary = args.binary.resolve()
    proposal = HERE / 'proposal-v2'
    protocol_path = HERE / 'protocol-v3.json'
    protocol = json.loads(protocol_path.read_text())
    mapping = json.loads((proposal / 'mapping.json').read_text())
    for name, expected in json.loads((proposal / 'checksums.json').read_text()).items():
        require(hashlib.sha256((proposal / name).read_bytes()).hexdigest() == expected, f'Proposal changed: {name}')
    panels = {variant: capture.fingerprint(proposal / f'angiospermae-{variant}.yaml') for variant in ('before', 'after')}
    panel = yaml.safe_load((proposal / 'angiospermae-after.yaml').read_text())
    genes = {f'angiospermae_{p["gene"]}': 'assay_correction_evaluation' for p in panel['primers']}
    require(len(genes) == 10, 'Expected complete ten-target panel')
    binary = capture.fingerprint(args.binary)
    available = json.loads((ROOT / 'tmp/inclusive-release-20260929/availability.json').read_text())
    selected = [r for r in available['samples'] if r['accession'] in protocol['real']['sample_accessions']]
    require(len(selected) == 2 and all(r['status'] == 'verified' for r in selected), 'Real input receipts missing')
    samples = []
    if args.scope in ('synthetic', 'all'):
        for row in mapping['synthetic_inputs']:
            path = (proposal / row['fastq']).resolve()
            sample = {'scope': 'synthetic', 'sample_id': row['source_accession'], 'available_path': str(path),
                      'expected_prefix': {'records': row['read_count'], 'bases': row['read_count'] * row['read_length'],
                                          'size_bytes': path.stat().st_size, 'sha256': row['fastq_sha256']},
                      'expected_corrected_sha256': row['expected_corrected_sha256']}
            samples.append(sample)
    if args.scope in ('real', 'all'):
        for row in selected:
            require(row['expected_prefix'] == protocol['real']['prefixes'][row['accession']], 'Protocol/input prefix mismatch')
            samples.append({'scope': 'real', 'sample_id': row['accession'], 'available_path': row['available_path'],
                            'expected_prefix': row['expected_prefix']})
    for sample in samples:
        # These staged files consist of exactly the previously validated prefix.
        # Matching all bytes binds the prior FASTQ syntax/count/base validation.
        sample['input_file'] = capture.fingerprint(Path(sample['available_path']))
        require(all(sample['input_file'][k] == sample['expected_prefix'][k] for k in ('sha256', 'size_bytes')), 'Input bytes differ')
    cpus = [min(os.sched_getaffinity(0))]
    if not args.plan_only:
        require(bool(os.environ.get('SLURM_JOB_ID')), 'Execution requires a Slurm allocation')
        os.sched_setaffinity(0, cpus)
    bindings = {str(Path(__file__).resolve()): capture.fingerprint(Path(__file__).resolve()),
                str(SHARED): capture.fingerprint(SHARED),
                str(protocol_path): capture.fingerprint(protocol_path),
                str(proposal / 'mapping.json'): capture.fingerprint(proposal / 'mapping.json'),
                str(proposal / 'checksums.json'): capture.fingerprint(proposal / 'checksums.json')}
    for module in list(sys.modules.values()):
        source = getattr(module, '__file__', None)
        if source:
            path = Path(source).resolve()
            if path.is_file() and ROOT in path.parents and path.suffix == '.py':
                bindings[str(path)] = capture.fingerprint(path)
    launcher = ROOT / 'benchmarks/benchmark_results/v3.2-uncertain-products-20260927/launch_bounded.py'
    bindings[str(launcher)] = capture.fingerprint(launcher)
    jobs = []
    for sample in samples:
        for variant in ('before', 'after'):
            sample_prefix = f'{sample["scope"]}_{sample["sample_id"].replace(".", "_")}_{variant}'
            directory = args.output / sample_prefix
            outputs = directory / 'outputs'
            command = [str(args.binary), '-k', '19', '-t', '1', '--chunks', '0', '--max-reads', '1000000',
                       '-o', str(outputs) + '/', '-s', sample_prefix, '--pcr-panel-file', panels[variant]['path'], sample['available_path']]
            jobs.append({**sample, 'variant': variant, 'sample_prefix': sample_prefix, 'directory': str(directory),
                         'outputs_directory': str(outputs), 'command': command,
                         'address_space_bytes': protocol[sample['scope']]['address_space_limit_bytes'], 'status': 'planned'})
    args.output.mkdir(parents=True)
    plan = {'scope': args.scope, 'protocol': protocol_path.name, 'binary': binary, 'panels': panels,
            'bindings': bindings, 'jobs': jobs, 'cpus': cpus, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'timing_scope': 'Bounded functionality checks; not a paired performance benchmark.'}
    dump_new(args.output / 'plan.json', plan)
    if args.plan_only:
        print(json.dumps({'status': 'plan_only', 'jobs': len(jobs), 'output': str(args.output)}))
        return
    environment = dict(os.environ, LC_ALL='C')
    for job in jobs:
        try:
            for fingerprint in [binary, panels[job['variant']], job['input_file'], *bindings.values()]:
                capture.verify(fingerprint)
            capture.ADDRESS_SPACE = job['address_space_bytes']
            capture.TIMEOUT = protocol['execution']['timeout_seconds_per_invocation']
            invocation = {'directory': job['directory'], 'outputs': job['outputs_directory'],
                          'command': job['command'], 'sample_prefix': job['sample_prefix'], 'version': 'candidate'}
            job['execution'] = capture.run_timed(invocation, cpus, environment)
            require(job['execution']['status'] == 'pass', 'Sharkmer execution failed or exceeded bound')
            job['outputs'] = capture.capture(Path(job['outputs_directory']), invocation, job, job['command'], genes)
            for fingerprint in [binary, panels[job['variant']], job['input_file'], *bindings.values()]:
                capture.verify(fingerprint)
            job['status'] = 'pass'
        except Exception as error:
            job['status'] = 'failed'
            job['error'] = f'{type(error).__name__}: {error}'
        dump_new(Path(job['directory']) / 'receipt.json' if Path(job['directory']).exists() else args.output / f'{job["sample_prefix"]}.failed.json', job)
        print(json.dumps({'sample': job['sample_id'], 'variant': job['variant'], 'status': job['status'],
                          'wall_seconds': job.get('execution', {}).get('metrics', {}).get('wall_seconds'), 'error': job.get('error')}), flush=True)
    comparisons = assess(jobs, genes, mapping['references'])
    summary = {'status': 'pass' if all(j['status'] == 'pass' for j in jobs) and all(r['status'] == 'pass' for r in comparisons) else 'review_required',
               'jobs': jobs, 'comparisons': comparisons,
               'limits': ['Exact reference support is separate from full sample phasing.', 'Synthetic source reconstruction is not real-sample validation.',
                          'Only the trnV-atpE assay changes; report any changes to the other nine targets.', 'No release or performance acceptance is implied.']}
    dump_new(args.output / 'summary.json', summary)
    print(json.dumps({'status': summary['status'], 'summary': str(args.output / 'summary.json')}))
    if summary['status'] != 'pass':
        raise SystemExit(1)


if __name__ == '__main__':
    main()

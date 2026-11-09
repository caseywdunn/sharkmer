#!/usr/bin/env python3
"""Receipt-bound three-way parity harness. Refuses execution without --execute."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.insert(0, '/home/claude/repos/sharkmer/scripts')
from sharkmer_validate import runner as validation_runner


ROOT = Path('/tmp/sharkmer-withheld-diagnostics-20260913')
PROTOCOL_PATH = ROOT / 'protocol.json'
TIME_BINARY = Path('/usr/bin/time')
TASKSET_BINARY = Path('/usr/bin/taskset')
PRLIMIT_BINARY = Path('/usr/bin/prlimit')


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open('rb') as handle:
        while block := handle.read(1024 * 1024):
            hasher.update(block)
    return hasher.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')


def write_new_json(path: Path, value: object) -> None:
    serialized = (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o664)
    with os.fdopen(descriptor, 'wb') as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())


def inspect_prefix(path: Path, expected_records: int) -> dict[str, int | str]:
    hasher = hashlib.sha256()
    bases = 0
    with path.open('rb') as handle:
        for ordinal in range(1, expected_records + 1):
            header = handle.readline()
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not all((header, sequence, separator, quality)):
                raise ValueError(f'{path}: EOF before record {ordinal}')
            if not header.startswith(b'@') or not separator.startswith(b'+'):
                raise ValueError(f'{path}: malformed four-line FASTQ at record {ordinal}')
            sequence_bases = sequence.rstrip(b'\r\n')
            quality_bases = quality.rstrip(b'\r\n')
            if len(sequence_bases) != len(quality_bases):
                raise ValueError(f'{path}: sequence/quality length mismatch at record {ordinal}')
            for line in (header, sequence, separator, quality):
                hasher.update(line)
            bases += len(sequence_bases)
        return {'records': expected_records, 'bases': bases, 'size_bytes': handle.tell(), 'sha256': hasher.hexdigest()}


def prewarm(path: Path, records: int) -> None:
    with path.open('rb') as handle:
        for unused in range(records):
            for unused in range(4):
                if not handle.readline():
                    raise ValueError(f'{path}: EOF during prefix prewarm')


def load_reviewed_snapshot(path: Path) -> dict[str, str]:
    snapshot = json.loads(path.read_text())
    files = snapshot.get('files')
    if not isinstance(files, dict) or not files:
        raise ValueError('Reviewed source snapshot requires a nonempty files mapping')
    for name, expected in files.items():
        if not isinstance(name, str) or not isinstance(expected, str) or digest(Path(name)) != expected:
            raise ValueError(f'Reviewed source snapshot differs: {name}')
    return files


def source_receipt(protocol: dict, binaries: dict[str, Path], diagnostic_args: list[str], inputs: dict, changed_snapshot: Path, changed_build_receipt: Path) -> dict:
    panel = Path(protocol['panel']['path'])
    source_paths = {
        'protocol': PROTOCOL_PATH,
        'runner': Path(__file__),
        'panel': panel,
        'baseline_binary': binaries['baseline_off'],
        'changed_binary': binaries['changed_off'],
        'interpreter': Path(sys.executable).resolve(),
        'time_binary': TIME_BINARY,
        'taskset_binary': TASKSET_BINARY,
        'prlimit_binary': PRLIMIT_BINARY,
        'validator': Path('/home/claude/repos/sharkmer/scripts/sharkmer_validate/runner.py'),
        'baseline_build_receipt': Path(protocol['baseline_binary']['build_receipt']),
        'changed_source_snapshot': changed_snapshot,
        'changed_build_receipt': changed_build_receipt,
        'analyzer': ROOT / 'analyze_withheld_results.py',
        'analyzer_tests': ROOT / 'test_analyze_withheld_results.py',
        'historical_localization': ROOT / 'historical_localization.json',
        'historical_map': ROOT / 'historical_lost_retained_map.json',
        'protocol_snapshot': ROOT / 'WITHHELD_PATH_DIAGNOSTICS.snapshot.md',
    }
    focused = protocol['focused_diagnostic_only']
    source_paths['focused_manifest'] = Path(focused['manifest'])
    source_paths['focused_identity_receipt'] = Path(focused['identity_receipt'])
    for row in focused['runs_after_full_panel']:
        source_paths[f"focused_panel_{row['accession']}"] = Path(row['panel'])
    hashes = {name: digest(path) for name, path in source_paths.items()}
    if hashes['panel'] != protocol['panel']['sha256_before_run']:
        raise ValueError('Pinned panel hash differs from protocol')
    if hashes['baseline_binary'] != protocol['baseline_binary']['sha256']:
        raise ValueError('Baseline binary hash differs from protocol')
    if digest(binaries['changed_on']) != hashes['changed_binary']:
        raise ValueError('Changed off/on paths must have identical binary bytes')
    if hashes['focused_manifest'] != focused['manifest_sha256']:
        raise ValueError('Focused-panel manifest hash differs from protocol')
    if hashes['focused_identity_receipt'] != focused['identity_receipt_sha256']:
        raise ValueError('Focused-panel identity receipt hash differs from protocol')
    for row in focused['runs_after_full_panel']:
        if hashes[f"focused_panel_{row['accession']}"] != row['panel_sha256']:
            raise ValueError(f"Focused panel hash differs for {row['accession']}")
    baseline_receipt = json.loads(Path(protocol['baseline_binary']['build_receipt']).read_text())
    changed_receipt = json.loads(changed_build_receipt.read_text())
    if baseline_receipt.get('frozen_binary_sha256') != hashes['baseline_binary']:
        raise ValueError('Baseline build receipt does not bind frozen binary')
    if changed_receipt.get('binary_sha256') != hashes['changed_binary']:
        raise ValueError('Changed build receipt does not bind changed binary')
    reviewed_files = load_reviewed_snapshot(changed_snapshot)
    return {'hashes': hashes, 'diagnostic_args': diagnostic_args, 'inputs': inputs, 'reviewed_files': reviewed_files, 'captured_at': utc_now()}


def stable_receipt(receipt: dict) -> dict:
    return {key: receipt[key] for key in ('hashes', 'diagnostic_args', 'inputs', 'reviewed_files')}


def validate_inputs(protocol: dict) -> dict:
    verified = {}
    for accession, sample in protocol['samples'].items():
        path = Path(sample['input'])
        full_hash = digest(path)
        if full_hash != sample['full_sha256']:
            raise ValueError(f'{accession}: full input hash changed')
        prefix = inspect_prefix(path, sample['prefix']['records'])
        if prefix != sample['prefix']:
            raise ValueError(f'{accession}: consumed input prefix changed')
        verified[accession] = {'path': str(path), 'full_sha256': full_hash, 'prefix': prefix}
    return verified


def parse_rss(time_log: Path) -> int | None:
    for line in time_log.read_text(errors='replace').splitlines():
        match = re.fullmatch(r'\s*Maximum resident set size \(kbytes\):\s*(\d+)\s*', line)
        if match:
            return int(match.group(1)) * 1024
    return None


def parse_fasta(path: Path) -> list[dict[str, str]]:
    records = []
    header = None
    sequence_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if header is not None:
                sequence = ''.join(sequence_parts)
                records.append({'header': header, 'sequence': sequence, 'sha256': hashlib.sha256(sequence.encode()).hexdigest()})
            header = line[1:]
            sequence_parts = []
        elif header is not None:
            sequence_parts.append(line)
        elif line:
            raise ValueError(f'{path}: content before FASTA header')
    if header is not None:
        sequence = ''.join(sequence_parts)
        records.append({'header': header, 'sequence': sequence, 'sha256': hashlib.sha256(sequence.encode()).hexdigest()})
    return records


def expected_gene_names(panel_path: Path) -> set[str]:
    panel = validation_runner.load_panel_yaml(panel_path)
    return {f"{panel['name']}_{gene}" for gene in validation_runner.panel_gene_names(panel)}


def strip_optional_payload(value: object, key: str) -> object:
    if isinstance(value, dict):
        return {name: strip_optional_payload(item, key) for name, item in value.items() if name != key}
    if isinstance(value, list):
        return [strip_optional_payload(item, key) for item in value]
    return value


def diagnostic_payload_paths(value: object, key: str, path: tuple[str, ...] = ()) -> list[tuple[str, ...]]:
    paths = []
    if isinstance(value, dict):
        for name, item in value.items():
            item_path = (*path, str(name))
            if name == key:
                paths.append(item_path)
            paths.extend(diagnostic_payload_paths(item, key, item_path))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            paths.extend(diagnostic_payload_paths(item, key, (*path, str(index))))
    return paths


def validate_diagnostic_payload(stats: dict, key: str, enabled: bool) -> None:
    paths = diagnostic_payload_paths(stats, key)
    for path in paths:
        if len(path) != 5 or path[0] != 'pcr_results' or not path[1].isdigit() or path[2] != 'threshold_diagnostics' or not path[3].isdigit() or path[4] != key:
            raise ValueError(f'Diagnostic payload appears outside threshold diagnostics: {path}')
    if not enabled and paths:
        raise ValueError('Diagnostic payload present while opt-in flag is absent')
    if enabled and not paths:
        raise ValueError('Diagnostic payload absent while opt-in flag is present')


def build_jobs(protocol: dict) -> list[dict]:
    jobs = [
        {'accession': accession, 'role': role, 'binary_role': role, 'panel': protocol['panel']['path'], 'focused_diagnostic_only': False}
        for accession, sample in protocol['samples'].items() for role in sample['roles']
    ]
    jobs.extend({
        'accession': row['accession'], 'role': 'changed_on_focused', 'binary_role': 'changed_on', 'panel': row['panel'],
        'focused_diagnostic_only': True, 'selected_targets': row['targets'],
    } for row in protocol['focused_diagnostic_only']['runs_after_full_panel'])
    if len(jobs) != 12 or len({(job['accession'], job['role']) for job in jobs}) != 12:
        raise ValueError('Protocol must define twelve distinct jobs')
    return jobs


def capture_outputs(output_directory: Path, sample_name: str, panel_path: Path, expected_k: int, expected_command: list[str], input_path: Path, prefix: dict, diagnostic_enabled: bool, diagnostic_key: str) -> dict:
    stats_path = output_directory / f'{sample_name}.stats.yaml'
    stats = validation_runner._parse_stats_yaml(stats_path)
    validation_runner._validate_stats_manifest(stats, sample_name, expected_k, expected_gene_names(panel_path), 'end-to-end', expected_command=expected_command)
    validation_runner._validate_output_manifest(stats, sample_name, output_directory)
    source_plan = stats.get('input_source')
    if not isinstance(source_plan, dict):
        raise ValueError('Stats manifest lacks input_source')
    observed_inputs = [str(Path(item).resolve()) for item in source_plan.get('inputs', [])]
    if source_plan.get('kind') != 'local_files' or observed_inputs != [str(input_path.resolve())]:
        raise ValueError('Stats input_source does not bind the local staged input')
    if source_plan.get('paired') is not False or source_plan.get('max_reads') != prefix['records'] or stats.get('chunks') != 0:
        raise ValueError('Stats input_source or chunks differs from protocol')
    if stats.get('n_reads_read') != prefix['records'] or stats.get('n_bases_read') != prefix['bases']:
        raise ValueError('Stats consumed count differs from pinned prefix')
    validate_diagnostic_payload(stats, diagnostic_key, diagnostic_enabled)
    fasta_records = []
    for product in validation_runner.parse_fasta_products(sample_name, output_directory):
        for record in product['products']:
            fasta_records.append({'gene': product['gene'], 'path': product['output_file'], **record, 'sha256': hashlib.sha256(record['sequence'].encode()).hexdigest()})
    return {
        'stats_path': str(stats_path),
        'stats_sha256': digest(stats_path),
        'stats': stats,
        'fasta_records_ordered': fasta_records,
        'fasta_files': {path.name: digest(path) for path in sorted(output_directory.glob('*.fasta'))},
    }


def threshold_diagnostics(stats: dict, targets: set[str] | None = None) -> dict:
    diagnostics = {}
    for entry in stats['pcr_results']:
        gene = entry['gene_name']
        if targets is None or gene in targets:
            diagnostics[gene] = entry.get('threshold_diagnostics')
    return diagnostics


def pcr_results(stats: dict, targets: set[str] | None = None) -> list[dict]:
    entries = stats['pcr_results']
    return entries if targets is None else [entry for entry in entries if entry['gene_name'] in targets]


def compare_values(left: object, right: object, label: str) -> list[str]:
    return [] if left == right else [label]


def compare_full_roles(records: list[dict], payload_key: str) -> list[dict]:
    count_fields = ('n_reads_read', 'n_bases_read', 'n_subreads_ingested', 'n_bases_ingested', 'n_kmers')
    by_accession: dict[str, dict[str, dict]] = {}
    for record in records:
        if not record['focused_diagnostic_only']:
            by_accession.setdefault(record['accession'], {})[record['role']] = record
    comparisons = []
    for accession, roles in sorted(by_accession.items()):
        baseline = roles['baseline_off']['outputs']
        changed_off = roles['changed_off']['outputs']
        changed_on = roles['changed_on']['outputs']
        failures = []
        for name, left, right in (('baseline_vs_changed_off', baseline, changed_off), ('changed_off_vs_changed_on', changed_off, changed_on)):
            for field in count_fields:
                failures.extend(f'{name}:{field}' for ignored in compare_values(left['stats'].get(field), right['stats'].get(field), field))
            failures.extend(f'{name}:ordered_fasta_records' for ignored in compare_values(left['fasta_records_ordered'], right['fasta_records_ordered'], 'ordered_fasta_records'))
            failures.extend(f'{name}:pcr_results' for ignored in compare_values(
                strip_optional_payload(pcr_results(left['stats']), payload_key),
                strip_optional_payload(pcr_results(right['stats']), payload_key),
                'pcr_results',
            ))
            left_diagnostics = strip_optional_payload(threshold_diagnostics(left['stats']), payload_key)
            right_diagnostics = strip_optional_payload(threshold_diagnostics(right['stats']), payload_key)
            failures.extend(f'{name}:threshold_diagnostics' for ignored in compare_values(left_diagnostics, right_diagnostics, 'threshold_diagnostics'))
        comparisons.append({'kind': 'full_panel', 'accession': accession, 'passed': not failures, 'failures': failures})
    return comparisons


def compare_focused_runs(records: list[dict], payload_key: str) -> list[dict]:
    full = {record['accession']: record for record in records if record['role'] == 'changed_on' and not record['focused_diagnostic_only']}
    comparisons = []
    for record in records:
        if not record['focused_diagnostic_only']:
            continue
        targets = {f"insecta_{target}" for target in record['selected_targets']}
        focused = record['outputs']
        full_record = full[record['accession']]['outputs']
        failures = []
        for field in ('n_reads_read', 'n_bases_read', 'n_subreads_ingested', 'n_bases_ingested', 'n_kmers'):
            failures.extend(f'count:{field}' for ignored in compare_values(focused['stats'].get(field), full_record['stats'].get(field), field))
        focused_products = [item for item in focused['fasta_records_ordered'] if item['gene'] in targets]
        full_products = [item for item in full_record['fasta_records_ordered'] if item['gene'] in targets]
        failures.extend('selected_ordered_fasta_records' for ignored in compare_values(focused_products, full_products, 'selected_ordered_fasta_records'))
        failures.extend('selected_pcr_results' for ignored in compare_values(
            strip_optional_payload(pcr_results(focused['stats'], targets), payload_key),
            strip_optional_payload(pcr_results(full_record['stats'], targets), payload_key),
            'selected_pcr_results',
        ))
        focused_diagnostics = strip_optional_payload(threshold_diagnostics(focused['stats'], targets), payload_key)
        full_diagnostics = strip_optional_payload(threshold_diagnostics(full_record['stats'], targets), payload_key)
        failures.extend('selected_threshold_diagnostics' for ignored in compare_values(focused_diagnostics, full_diagnostics, 'selected_threshold_diagnostics'))
        comparisons.append({'kind': 'focused_vs_full', 'accession': record['accession'], 'targets': sorted(targets), 'passed': not failures, 'failures': failures})
    return comparisons


def make_command(binary: Path, protocol: dict, sample: dict, panel_path: str, output_directory: Path, sample_name: str, diagnostic_args: list[str]) -> list[str]:
    common = protocol['common_invocation']
    return [
        str(PRLIMIT_BINARY), f"--as={common['address_space_limit_bytes']}", '--', str(TIME_BINARY), '-v', '-o', str(output_directory / 'time.txt'),
        str(TASKSET_BINARY), '-c', common['cpu_list'], str(binary), '-k', str(common['k']), '-t', str(common['threads']),
        '--chunks', str(common['chunks']), '--max-reads', str(common['max_reads']), '-o', str(output_directory), '-s', sample_name,
        '--pcr-panel-file', panel_path, *diagnostic_args, sample['input'],
    ]


def run_job(job: dict, protocol: dict, binaries: dict[str, Path], diagnostic_args: list[str], output_root: Path, source_before: dict) -> dict:
    output_directory = output_root / job['accession'] / job['role']
    if output_directory.exists() or output_directory.is_symlink():
        raise ValueError(f'Refusing existing output directory: {output_directory}')
    output_directory.mkdir(parents=True)
    sample = protocol['samples'][job['accession']]
    input_path = Path(sample['input'])
    prefix_before = inspect_prefix(input_path, sample['prefix']['records'])
    if prefix_before != sample['prefix']:
        raise ValueError(f"{job['accession']}: input prefix changed before command")
    prewarm(input_path, sample['prefix']['records'])
    sample_name = f"withheld_{job['accession']}"
    command = make_command(
        binaries[job['binary_role']], protocol, sample, job['panel'], output_directory,
        sample_name, diagnostic_args if job['binary_role'] == 'changed_on' else [],
    )
    sharkmer_command = command[command.index(str(binaries[job['binary_role']])):]
    started_at = utc_now()
    started = time.monotonic()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=protocol['common_invocation']['timeout_seconds'])
        returncode = process.returncode
        timed_out = False
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        returncode = None
        timed_out = True
    (output_directory / 'stdout.txt').write_text(stdout)
    (output_directory / 'stderr.txt').write_text(stderr)
    record = {
        'accession': job['accession'], 'role': job['role'], 'binary_role': job['binary_role'], 'command': command, 'started_at': started_at, 'finished_at': utc_now(),
        'wall_time_seconds': time.monotonic() - started, 'returncode': returncode, 'timed_out': timed_out,
        'source_before': source_before, 'input_prefix_before': prefix_before, 'focused_diagnostic_only': job['focused_diagnostic_only'],
        'selected_targets': job.get('selected_targets'), 'stdout_sha256': hashlib.sha256(stdout.encode()).hexdigest(),
        'stderr_sha256': hashlib.sha256(stderr.encode()).hexdigest(), 'output_directory': str(output_directory),
    }
    try:
        record['input_prefix_after'] = inspect_prefix(input_path, sample['prefix']['records'])
        if record['input_prefix_after'] != prefix_before:
            raise ValueError('Consumed input prefix changed during command')
        if returncode == 0 and not timed_out:
            record['maximum_rss_bytes'] = parse_rss(output_directory / 'time.txt')
            record['outputs'] = capture_outputs(
                output_directory, sample_name, Path(job['panel']), protocol['common_invocation']['k'], sharkmer_command, input_path, sample['prefix'],
                job['binary_role'] == 'changed_on', protocol['comparison_contract']['preregistered_payload_key'],
            )
            record['validated_output'] = True
        else:
            record['validated_output'] = False
    except (OSError, ValueError, yaml.YAMLError) as error:
        record['validated_output'] = False
        record['validation_error'] = str(error)
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--changed', type=Path, required=True)
    parser.add_argument('--changed-source-snapshot', type=Path, required=True)
    parser.add_argument('--changed-build-receipt', type=Path, required=True)
    parser.add_argument('--diagnostic-arg', action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    arguments = parser.parse_args()
    if not arguments.execute:
        print('Refusing to execute without --execute; protocol and harness are pre-run only.')
        return 2
    if arguments.diagnostic_arg != ['--diagnose-withheld-paths']:
        raise ValueError('Changed-on runs require exactly one reviewed argument: --diagnose-withheld-paths')
    if arguments.output.exists() or arguments.output.is_symlink():
        raise ValueError(f'Refusing existing output root: {arguments.output}')
    protocol = json.loads(PROTOCOL_PATH.read_text())
    binaries = {'baseline_off': arguments.baseline, 'changed_off': arguments.changed, 'changed_on': arguments.changed}
    input_before = validate_inputs(protocol)
    source_before = source_receipt(protocol, binaries, arguments.diagnostic_arg, input_before, arguments.changed_source_snapshot, arguments.changed_build_receipt)
    arguments.output.mkdir(parents=True)
    write_new_json(arguments.output / 'source_before.json', source_before)
    jobs = build_jobs(protocol)
    records = []
    for job in jobs:
        receipt_path = arguments.output / f"{job['accession']}.{job['role']}.receipt.json"
        try:
            record = run_job(job, protocol, binaries, arguments.diagnostic_arg, arguments.output, source_before)
        except (OSError, ValueError, yaml.YAMLError) as error:
            record = {
                'accession': job['accession'], 'role': job['role'], 'binary_role': job['binary_role'],
                'returncode': None, 'timed_out': False, 'validated_output': False, 'execution_error': str(error),
            }
        write_new_json(receipt_path, record)
        records.append(record)
        if record['returncode'] != 0 or record['timed_out'] or not record.get('validated_output'):
            break
    after_validation_error = None
    try:
        input_after = validate_inputs(protocol)
        source_after = source_receipt(protocol, binaries, arguments.diagnostic_arg, input_after, arguments.changed_source_snapshot, arguments.changed_build_receipt)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        input_after = None
        source_after = None
        after_validation_error = str(error)
    comparisons = []
    if after_validation_error is None and len(records) == 12 and all(record.get('validated_output') for record in records):
        payload_key = protocol['comparison_contract']['preregistered_payload_key']
        comparisons = compare_full_roles(records, payload_key) + compare_focused_runs(records, payload_key)
        write_new_json(arguments.output / 'comparison.json', {'schema_version': 1, 'comparisons': comparisons})
    write_new_json(arguments.output / 'execution.json', {
        'schema_version': 1, 'jobs': records, 'source_before': source_before, 'source_after': source_after,
        'after_validation_error': after_validation_error,
        'source_changed_during_execution': after_validation_error is not None or stable_receipt(source_before) != stable_receipt(source_after),
        'comparisons': comparisons,
    })
    return 0 if after_validation_error is None and len(records) == 12 and stable_receipt(source_before) == stable_receipt(source_after) and all(record['returncode'] == 0 and record.get('validated_output') for record in records) and all(comparison['passed'] for comparison in comparisons) else 1


if __name__ == '__main__':
    raise SystemExit(main())

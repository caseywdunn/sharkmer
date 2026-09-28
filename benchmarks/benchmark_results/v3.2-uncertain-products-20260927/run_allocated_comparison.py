#!/usr/bin/env python3
"""After compiled/final-review clearance, freeze and execute in one allocation."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from build_snapshot import fingerprint, write_new
import run_comparison as comparison
import run_compiled_checks as checks

HERE = Path(__file__).resolve().parent
EXPECTED_CHECKS = ['format', 'clippy-default', 'release-tests-default', 'clippy-fxhashmap', 'release-tests-fxhashmap', 'python-tests', 'reference-audit']


def verify_checks(path):
    before = fingerprint(path)
    document = json.loads(path.read_text())
    if document.get('status') != 'pass' or [row.get('name') for row in document.get('steps', [])] != EXPECTED_CHECKS:
        raise ValueError('All seven compiled/Python/reference checks must pass')
    if document['sources_before'] != document['sources_after'] or document['sources_before'] != checks.sources():
        raise ValueError('Check source differs from current source')
    if document['tools_before'] != document['tools_after'] or document['cargo_configs_before'] != document['cargo_configs_after']:
        raise ValueError('Tools/configuration changed during checks')
    if document['runtime_before'] != document['runtime_after'] or document['runtime_before'] != checks.runtime():
        raise ValueError('Check Python runtime differs from current runtime')
    if document['comparison_environment'] != checks.comparison_environment(os.environ):
        raise ValueError('Shared compilation environment differs from checked environment')
    plan = checks.command_plan(document['tools_before']['cargo']['path'], document['runtime_before']['interpreter']['path'])
    for row, (name, command) in zip(document['steps'], plan):
        if row['name'] != name or row['command'] != command or row['cwd'] != str(checks.ROOT):
            raise ValueError('Check command/cwd differs from required plan')
        if row.get('status') != 'pass' or row.get('returncode') != 0:
            raise ValueError('A required check did not succeed')
        comparison.verify(row['log'])
    for mapping in ('tools_before', 'cargo_configs_before', 'sources_before'):
        for row in document[mapping].values():
            comparison.verify(row)
    comparison.verify(before)
    return document


def call(command, output, name, timeout):
    log = output / f'{name}.log'
    started = time.monotonic()
    row = {'command': command, 'status': 'failed', 'timeout_seconds': timeout, 'process_group_killed': False}
    try:
        with log.open('x') as stream:
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            row['pid'] = process.pid
            try:
                row['returncode'] = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                row.update(process_group_killed=True, returncode=process.returncode)
                raise TimeoutError(f'{name} exceeded its outer orchestration deadline')
        if row['returncode']:
            raise ValueError(f'{name} exited with {row["returncode"]}; preserve this attempt')
        row['status'] = 'complete'
    except Exception as error:
        row['error'] = repr(error)
    finally:
        row['elapsed_seconds'] = time.monotonic() - started
        if log.is_file():
            row['log'] = fingerprint(log)
        write_new(output / f'{name}.orchestration.json', row)
    if row['status'] != 'complete':
        raise RuntimeError(row['error'])
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', type=Path, required=True)
    parser.add_argument('--baseline-build', type=Path, required=True)
    parser.add_argument('--candidate-build', type=Path, required=True)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--attempt-dir', type=Path, required=True)
    args = parser.parse_args()
    output = args.attempt_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'schema_version': 1, 'status': 'failed', 'argv': sys.argv, 'comparison_started': False}
    try:
        required = {Path(__file__).resolve(), HERE / 'run_compiled_checks.py', HERE / 'test_orchestration.py', HERE / 'checks.sbatch', HERE / 'build.sbatch', HERE / 'comparison.sbatch', args.checks.resolve(), args.baseline_build.resolve(), args.candidate_build.resolve(), HERE / 'inputs/current_prefix_availability.json'}
        comparison.review_sources(args.review.resolve(), required)
        checked = verify_checks(args.checks.resolve())
        for path, role in ((args.baseline_build, 'baseline'), (args.candidate_build, 'candidate')):
            build, unused_snapshot = comparison.build_info(path.resolve(), role)
            if build['cargo'] != checked['tools_before']['cargo'] or build['rustc'] != checked['tools_before']['rustc']:
                raise ValueError('Build and check compiler identities differ')
            if build['comparison_environment'] != checked['comparison_environment'] or build['cargo_configs'] != checked['cargo_configs_before']:
                raise ValueError('Build and check configuration/environment differ')
        receipt['prerequisites_before'] = {str(path): fingerprint(path) for path in sorted(required | {args.review.resolve()})}
        availability_path = HERE / 'inputs/current_prefix_availability.json'
        available = json.loads(availability_path.read_text())
        if available.get('status') != 'verified_bytes' or {row['accession'] for row in available['samples']} != {row[0] for row in comparison.ORDER} or len(available['samples']) != 3:
            raise ValueError('Verified three-prefix availability receipt is required')
        command = [str(Path(sys.executable).resolve()), str(HERE / 'run_comparison.py'), 'freeze', '--baseline-build', str(args.baseline_build.resolve()), '--candidate-build', str(args.candidate_build.resolve()), '--review', str(args.review.resolve()), '--run-dir', str(output / 'runs'), '--output', str(output / 'freeze.json')]
        for row in available['samples']:
            command.extend(['--read', f"{row['accession']}={row['file']['path']}"])
        receipt['freeze_call'] = call(command, output, 'freeze', 600)
        frozen = fingerprint(output / 'freeze.json')
        receipt['freeze'] = frozen
        receipt['comparison_started'] = True
        command = [str(Path(sys.executable).resolve()), str(HERE / 'run_comparison.py'), 'execute', '--freeze', frozen['path'], '--freeze-sha256', frozen['sha256']]
        receipt['execution_call'] = call(command, output, 'execute', 18800)
        for row in receipt['prerequisites_before'].values():
            comparison.verify(row)
        comparison.verify(frozen)
        execution = json.loads((output / 'runs/execution.json').read_text())
        if execution.get('status') != 'complete':
            raise ValueError('Execution lacks a complete success receipt')
        receipt['execution'] = fingerprint(output / 'runs/execution.json')
        receipt['status'] = 'complete'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        write_new(output / 'orchestration.json', receipt)
    return 0 if receipt['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())

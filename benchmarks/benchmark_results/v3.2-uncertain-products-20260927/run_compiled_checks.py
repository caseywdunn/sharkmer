#!/usr/bin/env python3
"""Run source-bound compilation/validation checks; never run the comparison inputs."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from build_snapshot import GENERATED_PANEL_DIRS, cargo_configs, fingerprint, write_new

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
STEP_TIMEOUT_SECONDS = 3600
MAX_LOG_BYTES = 64 * 1024**2


def sources():
    paths = {ROOT / 'Cargo.toml', ROOT / 'Cargo.lock', Path(__file__).resolve(), Path(__file__).resolve().with_name('build_snapshot.py')}
    for name in ('src', 'panels', 'tests', 'scripts'):
        for directory, child_dirs, files in os.walk(ROOT / name):
            child_dirs[:] = [child for child in child_dirs if child != '__pycache__' and not (name == 'panels' and child in GENERATED_PANEL_DIRS)]
            paths.update(Path(directory) / filename for filename in files if not filename.endswith('.pyc'))
    paths.update(HERE.glob('*.py'))
    paths.update(HERE.glob('*.sbatch'))
    return {str(path): fingerprint(path) for path in sorted(paths)}


def runtime():
    import yaml
    paths = [path for path in Path(yaml.__file__).resolve().parent.rglob('*') if path.is_file() and path.suffix in ('.py', '.so')]
    return {'python': sys.version, 'interpreter': fingerprint(Path(sys.executable).resolve()), 'yaml_version': yaml.__version__, 'yaml_sources': {str(path): fingerprint(path) for path in sorted(paths)}}


def comparison_environment(environment):
    # This is the shared subset required by the independently reviewed build receipt.
    return {key: environment.get(key) for key in ('CARGO_HOME', 'RUSTUP_HOME', 'RUSTUP_TOOLCHAIN', 'CARGO_BUILD_TARGET', 'CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'LDFLAGS')}


def command_plan(cargo, python):
    commands = [('format', [cargo, 'fmt', '--all', '--check'])]
    for backend, feature_args in [('default', []), ('fxhashmap', ['--no-default-features', '--features', 'fxhashmap'])]:
        commands.append((f'clippy-{backend}', [cargo, 'clippy', '--locked', '--offline', '--all-targets', *feature_args, '--', '-D', 'warnings']))
        commands.append((f'release-tests-{backend}', [cargo, 'test', '--locked', '--offline', '--release', *feature_args]))
    commands.extend([('python-tests', [python, '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_*.py']), ('reference-audit', [python, 'scripts/audit_panel_references.py', '--verify-existing'])])
    return commands


def run_step(name, command, cwd, environment, output, timeout=STEP_TIMEOUT_SECONDS):
    log = output / f'{name}.log'
    row = {'name': name, 'command': command, 'cwd': str(cwd), 'status': 'failed', 'timeout_seconds': timeout, 'polled_log_termination_threshold_bytes': MAX_LOG_BYTES, 'log_poll_interval_seconds': 1, 'started_unix': time.time()}
    started = time.monotonic()
    try:
        with log.open('x') as stream:
            process = subprocess.Popen(command, cwd=cwd, env=environment, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                while True:
                    try:
                        process.wait(timeout=min(1, timeout))
                        break
                    except subprocess.TimeoutExpired:
                        if time.monotonic() - started > timeout:
                            raise TimeoutError('Compiled-check step exceeded timeout')
                        if log.stat().st_size > MAX_LOG_BYTES:
                            raise ValueError('Compiled-check log exceeds polled termination threshold')
                if log.stat().st_size > MAX_LOG_BYTES:
                    raise ValueError('Compiled-check log exceeds polled termination threshold')
            except BaseException:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                raise
        row['returncode'] = process.returncode
        if process.returncode != 0:
            raise ValueError(f'Check exited with status {process.returncode}')
        row['status'] = 'pass'
    except Exception as error:
        row['error'] = repr(error)
    finally:
        row['elapsed_seconds'] = time.monotonic() - started
        if log.is_file():
            row['log'] = fingerprint(log)
        write_new(output / f'{name}.json', row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cargo', type=Path, required=True)
    parser.add_argument('--rustc', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'schema_version': 1, 'status': 'failed', 'argv': sys.argv, 'steps': [], 'comparison_input_runs_performed': False}
    try:
        receipt['sources_before'] = sources()
        receipt['runtime_before'] = runtime()
        environment = dict(os.environ)
        rejected = [key for key, value in environment.items() if value and (key == 'RUSTC_BOOTSTRAP' or key.startswith('CARGO_PROFILE_') or key == 'CARGO_BUILD_RUSTFLAGS' or (key.startswith('CARGO_TARGET_') and key.endswith('_RUSTFLAGS')))]
        if rejected:
            raise ValueError(f'Uncontrolled compiler/profile environment: {sorted(rejected)}')
        for key in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER'):
            environment.pop(key, None)
        sysroot = Path(subprocess.check_output([str(args.rustc.absolute()), '--print', 'sysroot'], text=True, timeout=30).strip()).resolve()
        tools = {name: (sysroot / 'bin' / name).resolve() for name in ('cargo', 'rustc', 'rustfmt', 'cargo-fmt', 'cargo-clippy', 'clippy-driver')}
        receipt['tools_before'] = {name: fingerprint(path) for name, path in tools.items()}
        requested_version = subprocess.check_output([str(args.cargo.absolute()), '--version'], text=True, timeout=30).strip()
        actual_version = subprocess.check_output([str(tools['cargo']), '--version'], text=True, timeout=30).strip()
        if requested_version != actual_version:
            raise ValueError('Requested cargo and compiler sysroot select different toolchains')
        receipt['cargo_version'] = actual_version
        receipt['rustc_version'] = subprocess.check_output([str(tools['rustc']), '-vV'], text=True, timeout=30).strip()
        cpus = sorted(os.sched_getaffinity(0))[:2]
        if len(cpus) != 2:
            raise ValueError('Compiled checks require a two-CPU allocation')
        os.sched_setaffinity(0, cpus)
        receipt['affinity'] = sorted(os.sched_getaffinity(0))
        receipt['slurm'] = {key: os.environ.get(key) for key in ('SLURM_JOB_ID', 'SLURM_JOB_CPUS_PER_NODE', 'SLURM_CPUS_PER_TASK', 'SLURM_JOB_NODELIST')}
        environment.update(PATH=str(sysroot / 'bin') + os.pathsep + environment.get('PATH', ''), RUSTC=str(tools['rustc']), RUSTFMT=str(tools['rustfmt']), CARGO_TARGET_DIR=str(output / 'target'), CARGO_BUILD_JOBS='2', LC_ALL='C')
        receipt['cargo_configs_before'] = cargo_configs(ROOT, environment)
        receipt['comparison_environment'] = comparison_environment(environment)
        receipt['environment'] = {key: environment.get(key) for key in ('CARGO_HOME', 'RUSTUP_HOME', 'RUSTUP_TOOLCHAIN', 'RUSTC', 'RUSTFMT', 'CARGO_TARGET_DIR', 'CARGO_BUILD_JOBS', 'RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER', 'RUSTC_BOOTSTRAP', 'LC_ALL')}
        commands = command_plan(str(tools['cargo']), str(Path(sys.executable).resolve()))
        for name, command in commands:
            row = run_step(name, command, ROOT, environment, output)
            receipt['steps'].append(row)
            if row['status'] != 'pass':
                raise ValueError(f'Check failed: {name}')
            if sources() != receipt['sources_before']:
                raise ValueError('Reviewed source changed during compiled checks')
        receipt['status'] = 'pass'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        try:
            receipt['sources_after'] = sources()
            if receipt.get('sources_before') != receipt['sources_after']:
                raise ValueError('Source changed during checks')
            receipt['runtime_after'] = runtime()
            if receipt.get('runtime_before') != receipt['runtime_after']:
                raise ValueError('Python runtime changed during checks')
            if 'tools_before' in receipt:
                receipt['tools_after'] = {name: fingerprint(Path(row['path'])) for name, row in receipt['tools_before'].items()}
                if receipt['tools_after'] != receipt['tools_before']:
                    raise ValueError('Toolchain changed during checks')
            if 'cargo_configs_before' in receipt:
                receipt['cargo_configs_after'] = cargo_configs(ROOT, os.environ)
                if receipt['cargo_configs_before'] != receipt['cargo_configs_after']:
                    raise ValueError('Cargo configuration changed during checks')
        except Exception as error:
            receipt.update(status='failed', bookend_error=repr(error))
        write_new(output / 'checks.json', receipt)
    return 0 if receipt['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

#!/usr/bin/env python3
"""Export immutable compiled inputs; optionally build with a caller-selected toolchain."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[3]
BASELINE = '956e771437a2410a7bcbb1504b033c61dd1e11b6'
GENERATED_PANEL_DIRS = {'validation_runs', 'validation_results', 'validation_reports'}


def fingerprint(path):
    path = Path(path).absolute()
    initial = path.lstat()
    if not stat.S_ISREG(initial.st_mode):
        raise ValueError(f'Expected a regular file, not a symlink: {path}')
    identity = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
    h = hashlib.sha256()
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW), 'rb') as stream:
        if identity(os.fstat(stream.fileno())) != identity(initial):
            raise ValueError(f'File changed while opening: {path}')
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
        if identity(os.fstat(stream.fileno())) != identity(initial):
            raise ValueError(f'File changed while hashing: {path}')
    if identity(path.lstat()) != identity(initial):
        raise ValueError(f'File changed after hashing: {path}')
    return {'path': str(path), 'sha256': h.hexdigest(), 'size_bytes': initial.st_size}


def write_new(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def tree(path):
    return {str(item.relative_to(path)): fingerprint(item) for item in sorted(path.rglob('*')) if item.is_file()}


def export(role, output):
    output.mkdir(parents=True, exist_ok=False)
    source = output / 'source'
    source.mkdir()
    if role == 'baseline':
        archive = output / 'source.tar'
        subprocess.run(['git', 'archive', '--format=tar', f'--output={archive}', BASELINE], cwd=ROOT, check=True)
        with tarfile.open(archive) as bundle:
            bundle.extractall(source, filter='data')
        provenance = {'commit': BASELINE, 'archive': fingerprint(archive)}
    else:
        # Includes every current compile-time panel and test module. Do not copy target/.git.
        for name in ('Cargo.toml', 'Cargo.lock', 'src', 'panels', 'tests'):
            original = ROOT / name
            if original.is_dir():
                ignored = ('__pycache__', '*.pyc', *(sorted(GENERATED_PANEL_DIRS) if name == 'panels' else ()))
                shutil.copytree(original, source / name, ignore=shutil.ignore_patterns(*ignored))
            else:
                shutil.copy2(original, source / name)
        if (ROOT / 'build.rs').exists():
            shutil.copy2(ROOT / 'build.rs', source / 'build.rs')
        provenance = {'working_tree_compiled_inputs': True, 'repository_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()}
    snapshot = {'schema_version': 1, 'role': role, 'source': str(source), 'files': tree(source), 'provenance': provenance, 'exporter': fingerprint(__file__)}
    write_new(output / 'snapshot.json', snapshot)
    return snapshot


def cargo_configs(source, environment):
    cargo_home = Path(environment.get('CARGO_HOME', str(Path.home() / '.cargo')))
    locations = [cargo_home, *[parent / '.cargo' for parent in (source, *source.parents)]]
    return {str(path.resolve()): fingerprint(path) for location in locations for name in ('config', 'config.toml') if (path := location / name).is_file()}


def build(output, cargo, rustc):
    if any((output / name).exists() for name in ('build.json', 'cargo.jsonl', 'cargo.stderr', 'target', 'sharkmer')):
        raise ValueError('Build attempt already exists; preserve it and export a fresh snapshot')
    source = None
    receipt = {'schema_version': 1, 'status': 'failed', 'builder': fingerprint(__file__), 'requested_cargo': str(cargo), 'requested_rustc': str(rustc)}
    try:
        snapshot = json.loads((output / 'snapshot.json').read_text())
        source = Path(snapshot['source'])
        receipt.update(snapshot=fingerprint(output / 'snapshot.json'), source_before=tree(source))
        if receipt['source_before'] != snapshot['files'] or fingerprint(snapshot['exporter']['path']) != snapshot['exporter']:
            raise ValueError('Snapshot or exporter changed before build; export a fresh attempt')
        environment = dict(os.environ)
        rejected = [key for key, value in environment.items() if value and (key == 'RUSTC_BOOTSTRAP' or key.startswith('CARGO_PROFILE_') or key == 'CARGO_BUILD_RUSTFLAGS' or (key.startswith('CARGO_TARGET_') and key.endswith('_RUSTFLAGS')))]
        if rejected:
            raise ValueError(f'Uncontrolled compiler/profile environment: {sorted(rejected)}')
        for key in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER'):
            environment.pop(key, None)
        # Resolve rustup proxies before building, and invoke concrete tools thereafter.
        sysroot = Path(subprocess.check_output([str(rustc), '--print', 'sysroot'], text=True, timeout=30).strip()).resolve()
        effective_rustc, effective_cargo = (sysroot / 'bin/rustc').resolve(), (sysroot / 'bin/cargo').resolve()
        if not effective_rustc.is_file() or not effective_cargo.is_file():
            raise ValueError('Compiler sysroot does not contain concrete cargo/rustc')
        requested_cargo_version = subprocess.check_output([str(cargo), '--version'], text=True, timeout=30).strip()
        actual_cargo_version = subprocess.check_output([str(effective_cargo), '--version'], text=True, timeout=30).strip()
        if requested_cargo_version != actual_cargo_version:
            raise ValueError('Requested cargo and compiler sysroot select different toolchains')
        environment['CARGO_TARGET_DIR'] = str(output / 'target')
        environment['RUSTC'] = str(effective_rustc)
        command = [str(effective_cargo), 'build', '--locked', '--offline', '--release', '--message-format=json-render-diagnostics', '-j', '2']
        receipt.update(command=command, cwd=str(source), cargo=fingerprint(effective_cargo), rustc=fingerprint(effective_rustc), cargo_version=actual_cargo_version, rustc_version=subprocess.check_output([str(effective_rustc), '-vV'], text=True, timeout=30).strip(), compiler_sysroot=str(sysroot), cargo_configs=cargo_configs(source, environment), comparison_environment={key: environment.get(key) for key in ('CARGO_HOME', 'RUSTUP_HOME', 'RUSTUP_TOOLCHAIN', 'CARGO_BUILD_TARGET', 'CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'LDFLAGS')}, environment={key: environment.get(key) for key in ('CARGO_TARGET_DIR', 'RUSTC', 'RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER')})
        with (output / 'cargo.jsonl').open('x') as stdout, (output / 'cargo.stderr').open('x') as stderr:
            process = subprocess.Popen(command, cwd=source, env=environment, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                process.wait(timeout=3600)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                receipt['timed_out'] = True
                raise
        receipt['returncode'] = process.returncode
        if process.returncode:
            raise ValueError(f'cargo failed: {process.returncode}')
        artifacts = [json.loads(line) for line in (output / 'cargo.jsonl').read_text().splitlines() if line.startswith('{')]
        artifacts = [item for item in artifacts if item.get('reason') == 'compiler-artifact' and item.get('executable') and item.get('target', {}).get('name') == 'sharkmer']
        if len(artifacts) != 1:
            raise ValueError('Expected exactly one sharkmer artifact')
        shutil.copy2(artifacts[0]['executable'], output / 'sharkmer')
        (output / 'sharkmer').chmod(0o555)
        receipt.update(binary=fingerprint(output / 'sharkmer'), artifact=artifacts[0], binary_version=subprocess.check_output([str(output / 'sharkmer'), '--version'], text=True).strip())
        receipt['status'] = 'complete'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        try:
            receipt['builder_after'] = fingerprint(__file__)
            if source is not None:
                receipt['source_after'] = tree(source)
            if receipt['builder'] != receipt['builder_after'] or receipt.get('source_before') != receipt.get('source_after'):
                raise ValueError('Builder/source changed during build')
            for key in ('cargo', 'rustc'):
                if key in receipt and fingerprint(receipt[key]['path']) != receipt[key]:
                    raise ValueError(f'{key} changed during build')
            if source is not None and 'cargo_configs' in receipt and cargo_configs(source, os.environ) != receipt['cargo_configs']:
                raise ValueError('Cargo configuration changed during build')
        except Exception as error:
            receipt.update(status='failed', bookend_error=repr(error))
        write_new(output / 'build.json', receipt)
    if receipt['status'] != 'complete':
        raise RuntimeError(receipt.get('error', 'Build failed'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['export', 'build'])
    parser.add_argument('--role', choices=['baseline', 'candidate'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cargo', type=Path)
    parser.add_argument('--rustc', type=Path)
    args = parser.parse_args()
    if args.action == 'export':
        if not args.role:
            parser.error('export requires --role')
        export(args.role, args.output.resolve())
    else:
        if not args.cargo or not args.rustc:
            parser.error('build requires explicit --cargo and --rustc paths')
        # Preserve cargo/rustc proxy basenames when the supplied tools are symlinks.
        build(args.output.resolve(), args.cargo.absolute(), args.rustc.absolute())


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Build pristine released and inclusive-development sources with one toolchain."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[3]
HELPERS = ROOT / 'benchmarks/benchmark_results/v3.2-uncertain-products-20260927'
sys.path.insert(0, str(HELPERS))
import build_snapshot as builder

COMMITS = {'baseline': '5a664680c91ad59f32b8c2a847b8fef37f34a0ae',
           'candidate': 'ed04bdec5178d055b47cc8085135ee5ae2600998'}


def main():
    if not __debug__:
        raise SystemExit('Build verification must not disable assertions')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'failed', 'versions': {}, 'wrapper': builder.fingerprint(__file__),
               'builder': builder.fingerprint(HELPERS / 'build_snapshot.py'),
               'slurm_job_id': os.environ.get('SLURM_JOB_ID')}
    try:
        tool_bin = Path(os.environ['RUSTUP_HOME']) / 'toolchains/stable-x86_64-unknown-linux-gnu/bin'
        for role, revision in COMMITS.items():
            commit = subprocess.check_output(['git', 'rev-parse', revision + '^{commit}'], cwd=ROOT, text=True).strip()
            if role == 'baseline':
                assert subprocess.check_output(['git', 'rev-parse', 'v3.1.0^{commit}'], cwd=ROOT, text=True).strip() == commit
            directory = output / role
            directory.mkdir()
            source = directory / 'source'
            source.mkdir()
            archive = directory / 'source.tar'
            tracked = set(subprocess.check_output(['git', 'ls-tree', '--name-only', commit], cwd=ROOT, text=True).splitlines())
            paths = [name for name in ('Cargo.toml', 'Cargo.lock', 'src', 'panels', 'tests', 'build.rs', '.cargo') if name in tracked]
            subprocess.run(['git', 'archive', '--format=tar', f'--output={archive}', commit, *paths], cwd=ROOT, check=True)
            with tarfile.open(archive) as stream:
                stream.extractall(source, filter='data')
            snapshot = {'schema_version': 1, 'role': role, 'source': str(source),
                        'files': builder.tree(source), 'exporter': builder.fingerprint(__file__),
                        'provenance': {'commit': commit, 'archive': builder.fingerprint(archive),
                                       'git_tree': subprocess.check_output(['git', 'rev-parse', commit + '^{tree}'], cwd=ROOT, text=True).strip(),
                                       'exported_paths': paths, 'pristine_git_export': True}}
            builder.write_new(directory / 'snapshot.json', snapshot)
            builder.build(directory, tool_bin / 'cargo', tool_bin / 'rustc')
            built = json.loads((directory / 'build.json').read_text())
            assert built['status'] == 'complete'
            assert set(built['artifact']['features']) == {'default', 'ahashmap'}
            assert built['artifact']['profile']['opt_level'] == '3'
            assert not built['artifact']['profile']['debug_assertions']
            receipt['versions'][role] = {'commit': commit, 'snapshot': builder.fingerprint(directory / 'snapshot.json'),
                                         'build': builder.fingerprint(directory / 'build.json'), 'binary': built['binary']}
        baseline, candidate = [json.loads(Path(receipt['versions'][role]['build']['path']).read_text()) for role in ('baseline', 'candidate')]
        for key in ('cargo', 'rustc', 'cargo_version', 'rustc_version', 'comparison_environment'):
            assert baseline[key] == candidate[key], key
        receipt['status'] = 'pass'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        try:
            assert receipt['wrapper'] == builder.fingerprint(__file__)
            assert receipt['builder'] == builder.fingerprint(HELPERS / 'build_snapshot.py')
            for version in receipt['versions'].values():
                for key in ('snapshot', 'build', 'binary'):
                    assert version[key] == builder.fingerprint(version[key]['path'])
        except Exception as error:
            receipt.update(status='failed', bookend_error=repr(error))
        builder.write_new(output / 'builds.json', receipt)
    print(json.dumps({'status': receipt['status'], 'error': receipt.get('error'), 'output': str(output)}))
    return 0 if receipt['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

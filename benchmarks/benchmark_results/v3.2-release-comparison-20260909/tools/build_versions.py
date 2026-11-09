import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone


REPOSITORY = Path('/home/claude/repos/sharkmer')
WORKSPACE = Path('/tmp/sharkmer-release-comparison')
VERSIONS = {
    'baseline': '5a664680c91ad59f32b8c2a847b8fef37f34a0ae',
    'candidate': 'ba64f573048b6c19a528278028810af5bd475b81',
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def source_tree_hash(directory):
    digest = hashlib.sha256()
    for path in sorted(directory.rglob('*')):
        if path.is_symlink():
            value = 'symlink:' + os.readlink(path)
        elif path.is_file():
            value = sha256(path)
        else:
            continue
        digest.update(str(path.relative_to(directory)).encode() + b'\0' + value.encode() + b'\n')
    return digest.hexdigest()


def capture(command, cwd=REPOSITORY):
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def main():
    WORKSPACE.mkdir(parents=True, exist_ok=True)
    if (WORKSPACE / 'builds.json').exists():
        raise RuntimeError('Build manifest already exists; do not overwrite immutable build provenance')
    rustc = capture(['rustc', '-vV'])
    cargo = capture(['cargo', '--version'])
    builds = {'schema_version': 1, 'created_at': datetime.now(timezone.utc).isoformat(), 'versions': {}}
    for version, commit in VERSIONS.items():
        source = WORKSPACE / 'sources' / version
        source.mkdir(parents=True, exist_ok=False)
        archive = WORKSPACE / f'{version}-source.tar'
        subprocess.run(['git', 'archive', '--format=tar', f'--output={archive}', commit], cwd=REPOSITORY, check=True)
        with tarfile.open(archive) as source_archive:
            source_archive.extractall(source, filter='data')
        initial_tree = source_tree_hash(source)
        target = WORKSPACE / 'targets' / version
        command = ['cargo', 'build', '--locked', '--offline', '--release', '--message-format=json-render-diagnostics']
        environment = dict(os.environ)
        for name in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER'):
            environment.pop(name, None)
        environment['CARGO_TARGET_DIR'] = str(target)
        output_path = WORKSPACE / f'{version}-cargo.jsonl'
        error_path = WORKSPACE / f'{version}-cargo.stderr'
        print(f'Building {version} from pristine archive {commit}', flush=True)
        with output_path.open('w') as output, error_path.open('w') as error:
            subprocess.run(command, cwd=source, env=environment, stdout=output, stderr=error, check=True)
        artifacts = []
        for line in output_path.read_text().splitlines():
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get('reason') == 'compiler-artifact' and event.get('executable') and event.get('target', {}).get('name') == 'sharkmer':
                artifacts.append(event)
        if len(artifacts) != 1:
            raise RuntimeError(f'Expected one sharkmer artifact, found {len(artifacts)}')
        artifact = artifacts[0]
        final_tree = source_tree_hash(source)
        if initial_tree != final_tree:
            raise RuntimeError(f'{version} export changed during build')
        binary = WORKSPACE / 'binaries' / f'sharkmer-{version}'
        binary.parent.mkdir(exist_ok=True)
        shutil.copy2(artifact['executable'], binary)
        binary.chmod(0o555)
        builds['versions'][version] = {
            'commit': commit,
            'git_tree_oid': capture(['git', 'rev-parse', f'{commit}^{{tree}}']),
            'release_tag': 'v3.1.0' if version == 'baseline' else None,
            'release_tag_commit': capture(['git', 'rev-parse', 'v3.1.0^{commit}']) if version == 'baseline' else None,
            'source_export': str(source),
            'source_archive': str(archive),
            'source_archive_sha256': sha256(archive),
            'source_tree_sha256_before_build': initial_tree,
            'source_tree_sha256_after_build': final_tree,
            'source_export_pristine_after_build': True,
            'binary_path': str(binary),
            'binary_sha256': sha256(binary),
            'binary_version': capture([str(binary), '--version'], source),
            'cargo_lock_sha256': sha256(source / 'Cargo.lock'),
            'build_command': command,
            'build_cwd': str(source),
            'build_environment': {'CARGO_TARGET_DIR': str(target), 'RUSTFLAGS': None, 'CARGO_ENCODED_RUSTFLAGS': None, 'RUSTC_WRAPPER': None, 'RUSTC_WORKSPACE_WRAPPER': None},
            'cargo_artifact_features': artifact['features'],
            'cargo_artifact_profile': artifact['profile'],
            'rustc': rustc,
            'cargo': cargo,
            'source_binding': 'built_from_verified_clean_git_archive_with_locked_dependencies',
        }
        (WORKSPACE / 'build-progress.json').write_text(json.dumps(builds, indent=2) + '\n')
        print(f"Built {version}: {builds['versions'][version]['binary_sha256']}", flush=True)
    (WORKSPACE / 'builds.json').write_text(json.dumps(builds, indent=2) + '\n')


if __name__ == '__main__':
    main()

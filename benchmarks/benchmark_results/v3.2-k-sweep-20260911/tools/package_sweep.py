import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def digest(path):
    checksum = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--directory', action='append', required=True)
    arguments = parser.parse_args()
    workspace = arguments.workspace.resolve()
    destination = arguments.destination
    if destination.exists():
        raise ValueError('Preserve earlier packages rather than overwrite')
    sources = [workspace / name for name in arguments.directory]
    if len({source.name for source in sources}) != len(sources):
        raise ValueError('Archive directory names collide')
    if any(not source.is_dir() or source.is_symlink() for source in sources):
        raise ValueError('Evidence directory missing or symlinked')
    destination.mkdir(parents=True)
    inventory = []
    with (destination / 'evidence.tar.gz').open('xb') as output:
        with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w|') as archive:
                for source in sources:
                    for path in sorted(source.rglob('*')):
                        if path.is_symlink():
                            raise ValueError(f'Unexpected evidence symlink: {path}')
                        if not path.is_file() or '__pycache__' in path.parts:
                            continue
                        relative = f'{source.name}/{path.relative_to(source).as_posix()}'
                        inventory.append({'path': relative, 'bytes': path.stat().st_size, 'sha256': digest(path)})
                        information = archive.gettarinfo(str(path), arcname=relative)
                        information.uid = information.gid = 0
                        information.uname = information.gname = ''
                        information.mtime = 0
                        information.mode = 0o644
                        with path.open('rb') as stream:
                            archive.addfile(information, stream)
    with (destination / 'ARCHIVE_CONTENTS.json').open('x') as output:
        json.dump(inventory, output, indent=2)
        output.write('\n')
    tools = destination / 'tools'
    tools.mkdir()
    for source in sorted(workspace.iterdir()):
        if not source.is_file() or source.suffix not in {'.py', '.json', '.md', '.log'}:
            continue
        target = tools if source.suffix == '.py' else destination
        shutil.copyfile(source, target / source.name)
    for name in ['report.md']:
        source = workspace / 'discovery-analysis' / name
        if source.is_file():
            shutil.copyfile(source, destination / name)
    files = sorted(path for path in destination.rglob('*') if path.is_file())
    with (destination / 'SHA256SUMS').open('x') as output:
        output.write(''.join(f'{digest(path)}  {path.relative_to(destination).as_posix()}\n' for path in files))
    print(f'Packaged {len(inventory)} original evidence files in {destination}')


if __name__ == '__main__':
    main()

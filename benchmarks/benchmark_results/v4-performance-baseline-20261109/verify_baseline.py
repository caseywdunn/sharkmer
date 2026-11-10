#!/usr/bin/env python3
"""Reconcile raw baseline outputs, bound source, and every frozen case."""
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tarfile

import yaml

BUNDLE = Path(__file__).resolve().parent
ROOT = BUNDLE.parents[2]
RUN = ROOT / 'tmp/v4-performance-baseline-20261109/run-v1'
SOURCE = ROOT / 'tmp/v4-performance-baseline-20261109/source'


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def inventory(directory):
    result = []
    for path in sorted(directory.glob('*.fasta')):
        lines = path.read_text().splitlines()
        current = None
        sequence = []
        def finish():
            if current is None:
                return
            fields = dict(token.split('=', 1) for token in current.split()[1:])
            bases = ''.join(sequence)
            status = fields.get('assembly_status', 'standard')
            assert status in ('standard', 'warning')
            assert len(bases) == int(fields['length'])
            result.append({'id': current.split()[0], 'gene': fields['gene'],
                           'length': len(bases), 'sha256': hashlib.sha256(bases.encode()).hexdigest(),
                           'assembly_status': status})
        for line in lines:
            if line.startswith('>'):
                finish()
                current, sequence = line[1:], []
            else:
                assert current is not None
                sequence.append(line)
        finish()
    return sorted(result, key=lambda row: row['id'])


def main():
    protocol = json.loads((BUNDLE / 'protocol.json').read_text())
    result = json.loads((RUN / 'result.json').read_text())
    assert result['status'] == 'passed'
    assert result['protocol_sha256'] == digest(BUNDLE / 'protocol.json')
    assert result['binary_sha256'] == digest(result['binary']) == '0d57521d5c6060df3e954d366b72fdccb0b3e40f9532162b400d62b14a59109c'
    source_hashes = {}
    for path in sorted(SOURCE.rglob('*')):
        if path.is_file():
            relative = path.relative_to(SOURCE).as_posix()
            expected = subprocess.check_output(['git', 'rev-parse', f'v3.2.1:{relative}'], cwd=ROOT, text=True).strip()
            observed = subprocess.check_output(['git', 'hash-object', str(path)], cwd=ROOT, text=True).strip()
            assert expected == observed, relative
            source_hashes[relative] = digest(path)
    (BUNDLE / 'source-bindings.json').write_text(json.dumps(source_hashes, indent=2) + '\n')
    summary = []
    seen = set()
    for sample in protocol['samples']:
        assert digest(sample['input']['path']) == sample['input']['sha256']
        assert digest(sample['panel_path']) == sample['panel_sha256']
        rows = [row for row in result['runs'] if row['accession'] == sample['accession']]
        assert sorted(row['repeat'] for row in rows) == [0, 1, 2, 3]
        expected_records = sorted([{**row, 'assembly_status': 'warning' if row['assembly_status'] == 'uncertain' else row['assembly_status']} for row in sample['expected_records']], key=lambda row: row['id'])
        timed = []
        for row in rows:
            identity = (row['accession'], row['repeat'])
            assert identity not in seen
            seen.add(identity)
            work = RUN / row['accession'] / str(row['repeat'])
            outputs = work / 'outputs'
            stats_path = outputs / f'{sample["sample"]}.stats.yaml'
            stats = yaml.safe_load(stats_path.read_text())
            manifest = yaml.safe_load((outputs / f'{sample["sample"]}.manifest.yaml').read_text())
            assert stats['run_status'] == manifest['status'] == 'complete'
            assert stats['run_id'] == manifest['run_id']
            assert {key: stats[key] for key in sample['expected_counts']} == sample['expected_counts']
            assert digest(stats_path) == row['stats_sha256']
            assert row['stage_timings'] == stats['stage_timings']
            assert row['metrics'] == json.loads((work / 'time.json').read_text())
            assert row['timed'] == (row['repeat'] > 0)
            assert inventory(outputs) == row['records'] == expected_records
            expected_command = [result['binary'], '-k', '19', '-t', '2', '--chunks', '0',
                                '--max-reads', str(sample['max_reads']), '-s', sample['sample'],
                                '-o', str(outputs), '--pcr-panel-file', sample['panel_path'], sample['input']['path']]
            assert row['command'] == expected_command
            for file in manifest['files']:
                assert Path(file['path']).name == file['path']
                assert digest(outputs / file['path']) == file['sha256']
            if row['timed']:
                timed.append(row['metrics'])
        summary.append({'accession': sample['accession'], 'panel': sample['panel'],
                        'records': len(expected_records),
                        'wall_median_s': statistics.median(row['wall_s'] for row in timed),
                        'wall_range_s': [min(row['wall_s'] for row in timed), max(row['wall_s'] for row in timed)],
                        'rss_median_kib': statistics.median(row['rss_kib'] for row in timed)})
    assert len(seen) == len(result['runs']) == 40
    (BUNDLE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    members = {path.relative_to(RUN).as_posix(): {'sha256': digest(path), 'bytes': path.stat().st_size}
               for path in sorted(RUN.rglob('*')) if path.is_file()}
    (BUNDLE / 'artifact-manifest.json').write_text(json.dumps(members, indent=2) + '\n')
    archive = BUNDLE / 'raw-results-v1.tar.gz'
    with tarfile.open(archive, 'w:gz', compresslevel=6) as output:
        output.add(RUN, arcname='run-v1')
    receipt = {'status': 'passed', 'invocations': 40, 'samples': 10,
               'distinct_baseline_records': sum(row['records'] for row in summary),
               'protocol_sha256': result['protocol_sha256'], 'binary_sha256': result['binary_sha256'],
               'source_file_count': len(source_hashes), 'raw_member_count': len(members),
               'raw_bytes': sum(member['bytes'] for member in members.values()),
               'archive_sha256': digest(archive), 'archive_bytes': archive.stat().st_size,
               'scope': 'engineering calibration baseline; no candidate performance/biological improvement claim'}
    (BUNDLE / 'verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()

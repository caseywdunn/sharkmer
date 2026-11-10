#!/usr/bin/env python3
"""Run frozen v3.2.1 calibration; later candidates require paired comparisons."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import time

import yaml


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def records(directory):
    rows = []
    for fasta in sorted(directory.glob('*.fasta')):
        for record in fasta.read_text().split('>')[1:]:
            header, *sequence = record.splitlines()
            fields = dict(part.split('=', 1) for part in header.split()[1:])
            sequence = ''.join(sequence)
            status = fields.get('assembly_status', 'standard')
            rows.append({'id': header.split()[0], 'gene': fields.get('target', fields.get('gene')),
                         'length': len(sequence), 'sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                         'assembly_status': 'warning' if status in ('warning', 'uncertain') else status})
            assert int(fields.get('sequence_length', fields.get('length'))) == len(sequence)
    return sorted(rows, key=lambda row: row['id'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    bundle = Path(__file__).resolve().parent
    protocol_file = bundle / 'protocol.json'
    protocol_sha = digest(protocol_file)
    assert protocol_sha == (bundle / 'protocol.sha256').read_text().split()[0]
    protocol = json.loads(protocol_file.read_text())
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    binary_sha = digest(binary)
    started = time.monotonic()
    receipt = {'status': 'running', 'protocol_sha256': protocol_sha,
               'binary_sha256': binary_sha, 'binary': str(binary),
               'baseline_commit': protocol['baseline_commit'],
               'started_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'runs': []}

    def save():
        path = output / 'result.json.tmp'
        path.write_text(json.dumps(receipt, indent=2) + '\n')
        path.replace(output / 'result.json')

    try:
        version = subprocess.check_output([str(binary), '--version'], text=True).strip()
        assert version == 'sharkmer 3.2.1 (https://github.com/caseywdunn/sharkmer)'
        receipt['version'] = version
        # This verifies the complete retained prefixes before timing, without
        # changing system caches. All runs use the same warm-cache protocol.
        for sample in protocol['samples']:
            source = sample['input']
            assert Path(source['path']).stat().st_size == source['size_bytes']
            assert digest(source['path']) == source['sha256'], source['path']
            assert digest(sample['panel_path']) == sample['panel_sha256']
        save()
        for sample in protocol['samples']:
            expected = sorted([{**row, 'assembly_status': 'warning' if row['assembly_status'] in ('warning', 'uncertain') else row['assembly_status']} for row in sample['expected_records']], key=lambda row: row['id'])
            for repeat in range(protocol['timed_replicates'] + protocol['warmup_replicates']):
                work = output / sample['accession'] / str(repeat)
                work.mkdir(parents=True)
                destination = work / 'outputs'
                metrics_file = work / 'time.json'
                command = [str(binary), '-k', str(protocol['k']), '-t', str(protocol['threads']),
                           '--chunks', str(protocol['chunks']), '--max-reads', str(sample['max_reads']),
                           '-s', sample['sample'], '-o', str(destination),
                           '--pcr-panel-file', sample['panel_path'], sample['input']['path']]
                timed = ['/usr/bin/time', '-f', '{"wall_s":%e,"user_s":%U,"system_s":%S,"rss_kib":%M}', '-o', str(metrics_file), *command]
                with (work / 'stdout.log').open('wb') as stdout, (work / 'stderr.log').open('wb') as stderr:
                    result = subprocess.run(timed, stdout=stdout, stderr=stderr,
                                            timeout=protocol['timeout_seconds'])
                assert result.returncode == 0, f'{sample["accession"]}: exit {result.returncode}'
                stats = yaml.safe_load((destination / f'{sample["sample"]}.stats.yaml').read_text())
                manifest = yaml.safe_load((destination / f'{sample["sample"]}.manifest.yaml').read_text())
                assert stats['run_status'] == manifest['status'] == 'complete'
                assert stats['run_id'] == manifest['run_id']
                for file in manifest['files']:
                    assert digest(destination / file['path']) == file['sha256']
                assert {key: stats[key] for key in sample['expected_counts']} == sample['expected_counts']
                observed = records(destination)
                assert observed == expected, f'{sample["accession"]}: sequence inventory differs'
                row = {'accession': sample['accession'], 'repeat': repeat,
                       'timed': repeat >= protocol['warmup_replicates'], 'command': command,
                       'metrics': json.loads(metrics_file.read_text()), 'stage_timings': stats['stage_timings'],
                       'records': observed, 'stats_sha256': digest(destination / f'{sample["sample"]}.stats.yaml')}
                receipt['runs'].append(row)
                save()
                print(sample['accession'], repeat, row['metrics'], flush=True)
        assert digest(binary) == binary_sha
        assert digest(protocol_file) == protocol_sha
        receipt['status'] = 'passed'
    except Exception as error:
        receipt['status'] = 'failed'
        receipt['error'] = str(error)
        raise
    finally:
        receipt['elapsed_s'] = time.monotonic() - started
        receipt['finished_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        save()


if __name__ == '__main__':
    main()

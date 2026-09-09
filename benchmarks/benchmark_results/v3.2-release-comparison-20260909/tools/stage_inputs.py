import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import time
import urllib.parse
import urllib.request

import yaml


REPOSITORY = Path('/home/claude/repos/sharkmer')
WORKSPACE = Path('/tmp/sharkmer-release-comparison')
PRIOR_RESULTS = Path('/tmp/sharkmer-v3.2-review/calibration-final/results')
LIMITS = (1000000, 2000000, 4000000, 8000000)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sequence_length(sequence):
    return len(sequence.rstrip(b'\r\n'))


def record(stream):
    header = stream.readline()
    if not header:
        return None
    sequence = stream.readline()
    separator = stream.readline()
    quality = stream.readline()
    if not header.startswith(b'@') or not separator.startswith(b'+') or not quality or sequence_length(sequence) != sequence_length(quality):
        raise ValueError('Malformed or truncated four-line FASTQ record during staging')
    return header + sequence + separator + quality, sequence_length(sequence)


def prior_selections():
    selections = {}
    for path in sorted(PRIOR_RESULTS.glob('*.yaml')):
        result = yaml.safe_load(path.read_text())
        for sample in result['samples']:
            depth = sample['depths'][0]
            selections[sample['accession']] = {
                'taxon': sample['taxon'],
                'selected_cache_entries': depth['input']['selected_cache_entries'],
                'n_reads_read': depth['run_stats']['n_reads_read'],
                'n_bases_read': depth['run_stats']['n_bases_read'],
            }
    return selections


def cached_prefix(selection):
    digest = hashlib.sha256()
    records = 0
    bases = 0
    size = 0
    for entry in selection['selected_cache_entries']:
        path = Path(entry['path'])
        if sha256(path) != entry['sha256']:
            raise ValueError(f'Historical input checksum mismatch: {path}')
        with gzip.open(path, 'rb') as source:
            while records < 1000000:
                value = record(source)
                if value is None:
                    break
                payload, length = value
                digest.update(payload)
                records += 1
                bases += length
                size += len(payload)
        if records >= 1000000:
            break
    if records != selection['n_reads_read'] or bases != selection['n_bases_read']:
        raise ValueError('Cached historical prefix disagrees with observed record/base counts')
    return {'sha256': digest.hexdigest(), 'records': records, 'bases': bases, 'size_bytes': size}


def metadata(accession):
    directory = WORKSPACE / 'ena'
    directory.mkdir(exist_ok=True)
    path = directory / f'{accession}.tsv'
    query = urllib.parse.urlencode({
        'accession': accession, 'result': 'read_run',
        'fields': 'run_accession,fastq_ftp,scientific_name,read_count,fastq_md5,fastq_bytes',
        'format': 'tsv',
    })
    url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + query
    if not path.exists():
        request = urllib.request.Request(url, headers={'User-Agent': 'Sharkmer-release-comparison/1'})
        with urllib.request.urlopen(request, timeout=120) as response:
            body = response.read()
        path.write_bytes(body)
    rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'))
    if len(rows) != 1 or rows[0].get('run_accession') != accession or not rows[0].get('fastq_ftp'):
        raise ValueError(f'Unexpected ENA metadata for {accession}')
    urls = []
    for location in rows[0]['fastq_ftp'].split(';'):
        parsed = urllib.parse.urlsplit(location if '://' in location else 'https://' + location)
        urls.append(urllib.parse.urlunsplit(('https', parsed.netloc, parsed.path, parsed.query, '')))
    return urls, {'query_url': url, 'path': str(path), 'sha256': sha256(path), 'fields': rows[0]}


def url_identity(url):
    parsed = urllib.parse.urlsplit(url)
    return parsed.netloc, parsed.path


def stage_sample(entry, selection):
    accession = entry['accession']
    directory = WORKSPACE / 'inputs'
    directory.mkdir(exist_ok=True)
    final_path = directory / f'{accession}.fastq'
    manifest_path = directory / f'{accession}.json'
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text())
        if final_path.exists() and sha256(final_path) == existing['sha256']:
            print(f'Reusing staged {accession}', flush=True)
            return existing
        raise ValueError(f'Staged input changed: {accession}')
    if final_path.exists():
        raise ValueError(f'Unreceipted staged output exists: {final_path}')
    requested = max(entry['max_reads'])
    urls, ena = metadata(accession)
    historical = cached_prefix(selection)
    cached = {url_identity(item['url']): item for item in selection['selected_cache_entries']}
    observed_order = [url_identity(item['url']) for item in selection['selected_cache_entries']]
    if observed_order != [url_identity(url) for url in urls[:len(observed_order)]]:
        raise ValueError(f'ENA URL order changed since historical comparison: {accession}')
    temporary = directory / f'{accession}.fastq.part'
    for attempt in range(1, 4):
        digest = hashlib.sha256()
        prefixes = {}
        records = 0
        bases = 0
        size = 0
        consumed = []
        started = time.monotonic()
        try:
            with temporary.open('wb', buffering=1024 * 1024) as output:
                for url in urls:
                    if records >= requested:
                        break
                    available = cached.get(url_identity(url))
                    use_cache = available and (available['complete'] or available['n_reads'] >= requested - records)
                    if use_cache:
                        transport = Path(available['path']).open('rb')
                        transport_provenance = {'kind': 'verified_historical_cache', 'path': available['path'], 'compressed_sha256': available['sha256']}
                    else:
                        request = urllib.request.Request(url, headers={'User-Agent': 'Sharkmer-release-comparison/1'})
                        transport = urllib.request.urlopen(request, timeout=120)
                        transport_provenance = {'kind': 'https_gzip_prefix', 'status': transport.status, 'headers': dict(transport.headers), 'full_upstream_checksum_verified': False}
                    before = records
                    exhausted = False
                    with transport, gzip.GzipFile(fileobj=transport, mode='rb') as decompressed:
                        while records < requested:
                            value = record(decompressed)
                            if value is None:
                                exhausted = True
                                break
                            payload, length = value
                            output.write(payload)
                            digest.update(payload)
                            records += 1
                            bases += length
                            size += len(payload)
                            if records in LIMITS:
                                prefixes[str(records)] = {'records': records, 'bases': bases, 'size_bytes': size, 'sha256': digest.hexdigest()}
                                print(f'{accession}: {records:,} records staged ({time.monotonic() - started:.1f}s)', flush=True)
                    consumed.append({'url': url, 'records': records - before, 'eof_verified': exhausted, 'transport': transport_provenance})
                output.flush()
                os.fsync(output.fileno())
            for depth in LIMITS:
                if depth <= requested and str(depth) not in prefixes:
                    prefixes[str(depth)] = {'records': records, 'bases': bases, 'size_bytes': size, 'sha256': digest.hexdigest()}
            if prefixes['1000000'] != historical:
                raise ValueError(f'{accession} new staging differs from historical 1M prefix')
            if sha256(temporary) != digest.hexdigest() or temporary.stat().st_size != size:
                raise ValueError(f'Staged disk integrity failure for {accession}')
            temporary.rename(final_path)
            final_path.chmod(0o444)
            result = {
                'id': accession, 'accession': accession, 'path': str(final_path),
                'sha256': digest.hexdigest(), 'size_bytes': size,
                'total_records': records, 'total_bases': bases,
                'requested_max_records': requested,
                'source_urls_ordered': urls,
                'source_consumption': consumed,
                'ena_metadata': ena,
                'prefixes': prefixes,
                'historical_1m_prefix': historical,
                'historical_1m_prefix_identical': True,
                'prepared_at': datetime.now(timezone.utc).isoformat(),
            }
            manifest_path.write_text(json.dumps(result, indent=2) + '\n')
            print(f'Completed {accession}: {records:,} records, {size / 1024**3:.2f} GiB', flush=True)
            return result
        except Exception as error:
            print(f'Staging {accession} attempt {attempt} failed: {error}', flush=True)
            if isinstance(error, ValueError) or attempt == 3:
                raise
            time.sleep(5 * attempt)
    raise RuntimeError('Staging retry loop exhausted')


def main():
    config = yaml.safe_load((REPOSITORY / 'benchmarks/benchmark.yaml').read_text())
    selections = prior_selections()
    inputs = []
    for entry in config['samples']:
        inputs.append(stage_sample(entry, selections[entry['accession']]))
    result = {
        'schema_version': 1,
        'preparation': {
            'description': 'Immutable record-aligned uncompressed FASTQ; each 1M prefix verified against earlier checksum-verified calibration cache.',
            'command': ['python', str(Path(__file__).resolve())],
            'tool_versions': {'python': platform.python_version()},
            'record_order': 'ENA URL order, sequential, unpaired',
            'limits': 'first8M for six configured sweep samples; first1M or EOF for remaining samples',
            'remote_integrity_note': 'Staged full/prefix SHA-256 verifies comparison identity. Partial remote downloads do not verify full upstream compressed MD5.',
        },
        'inputs': inputs,
    }
    (WORKSPACE / 'inputs.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()

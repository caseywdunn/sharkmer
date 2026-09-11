import collections
import hashlib
import json
from pathlib import Path

import yaml


panel_path = Path('/tmp/sharkmer-release-comparison/sources/candidate/panels/insecta.yaml')
panel_sha256 = hashlib.sha256(panel_path.read_bytes()).hexdigest()
assert panel_sha256 == '043a4d224f8186d58459287c54bec7935d5f4921a5087c14a2136cfac42df874'
panel = yaml.safe_load(panel_path.read_text())
references = [
    reference
    for group in panel['references']
    if group['gene'] == 'ITS_2'
    for reference in group['sequences']
    if reference['accession'] == 'AK281180'
]
assert len(references) == 1
sequence = references[0]['sequence']
assert hashlib.sha256(sequence.encode()).hexdigest() == '23163ca76f97f5d1685f0fcf4bda6ffc798bdf58c63953416833e78a480af68c'
probes = []
for kmer_length in [19, 23, 27, 31]:
    node_length = kmer_length - 1
    positions = collections.defaultdict(list)
    for start in range(len(sequence) - node_length + 1):
        positions[sequence[start:start + node_length]].append(start + 1)
    repeated = [
        {'sequence': node, 'positions_1_based': starts}
        for node, starts in sorted(positions.items())
        if len(starts) > 1
    ]
    probes.append({
        'k': kmer_length,
        'node_length': node_length,
        'distinct_nodes': len(positions),
        'repeated_distinct_nodes': len(repeated),
        'repeated_nodes': repeated,
    })
result = {
    'purpose': 'Supplementary reference-only substring diagnostic, not Sharkmer graph execution or read evidence',
    'reference': 'AK281180',
    'length': len(sequence),
    'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
    'panel': {'path': str(panel_path), 'sha256': panel_sha256},
    'orientation': 'Sense-strand substrings only; no canonicalization, other templates, sequencing errors, coverage filtering, or primer discovery',
    'interpretation': 'Repeated (k-1)-mers can identify intrinsic merges in a single-reference graph; their absence does not rule out cycles in the actual read-derived graph',
    'probes': probes,
}
with Path(__file__).with_name('reference-repeat-probe.json').open('x') as output:
    json.dump(result, output, indent=2)
    output.write('\n')
print(json.dumps([{key: value for key, value in probe.items() if key != 'repeated_nodes'} for probe in probes]))

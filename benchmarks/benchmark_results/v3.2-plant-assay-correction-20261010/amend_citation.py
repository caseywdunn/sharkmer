#!/usr/bin/env python3
"""Preserve proposal-v1 and make a citation-only, hash-bound proposal-v2."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    source = HERE / 'proposal-v1'
    output = HERE / 'proposal-v2'
    if output.exists():
        raise SystemExit('proposal-v2 must not already exist')
    checksums = json.loads((source / 'checksums.json').read_text())
    for name, expected in checksums.items():
        assert sha((source / name).read_bytes()) == expected, name
    output.mkdir()
    for name in checksums:
        if name != 'mapping.json':
            (output / name).write_bytes((source / name).read_bytes())
    panel_path = output / 'angiospermae-after.yaml'
    panel = panel_path.read_text()
    assert panel.count('Kress et al 2005 (forward/locus)') == 1
    panel = panel.replace('Kress et al 2005 (forward/locus)', 'Kress et al 2005 (legacy locus citation)')
    panel_path.write_text(panel)
    mapping = json.loads((source / 'mapping.json').read_text())
    mapping['after_panel_sha256'] = sha(panel.encode())
    mapping['citation_amendment'] = {
        'producer_sha256': sha(Path(__file__).read_bytes()),
        'source_mapping_sha256': sha((source / 'mapping.json').read_bytes()),
        'source_checksums_sha256': sha((source / 'checksums.json').read_bytes()),
        'reason': 'Do not imply that the original Kress supplemental forward oligo table was independently verified. Mapping, assay, reference sequences and synthetic inputs are unchanged.'}
    (output / 'mapping.json').write_text(json.dumps(mapping, indent=2, sort_keys=True) + '\n')
    (output / 'checksums.json').write_text(json.dumps({p.name: sha(p.read_bytes()) for p in sorted(output.iterdir())}, indent=2, sort_keys=True) + '\n')
    print(output)


if __name__ == '__main__':
    main()

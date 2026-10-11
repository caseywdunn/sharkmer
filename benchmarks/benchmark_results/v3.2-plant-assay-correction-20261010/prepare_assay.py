#!/usr/bin/env python3
"""Prepare reviewable plant-assay mapping, panel snapshots and synthetic inputs.

Reads only the pre-existing pinned catalog/annotated records; never edits them or
the production panel. Outputs require a fresh directory.
"""
import argparse
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PREVIOUS = ROOT / 'benchmarks/benchmark_results/v3.2-target-coverage-20261009/reference'
sys.path.insert(0, str(ROOT / 'scripts'))
from sharkmer_validate import reference_provenance as provenance

FORWARD = 'GTGTAAACGAGTTGCTCTACCA'
OLD_REVERSE = 'CGACATTTGCACATTTAGATGCTAC'
NEW_REVERSE = 'GGAGAGCAATTTGAAGAAATG'
SOURCE_IDS = ('NC_008326.1', 'NC_056221.1', 'MN990625.1')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--before-panel', type=Path, default=ROOT / 'panels/angiospermae.yaml')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must not already exist')
    before = args.before_panel.read_bytes()
    panel = yaml.safe_load(before)
    assay = next(p for p in panel['primers'] if p['gene'] == 'trnV-atpE')
    assert panel['panel_version'] == '1.0.2'
    assert assay['forward_seq'] == FORWARD and assay['reverse_seq'] == OLD_REVERSE
    assert [assay[k] for k in ('min_length', 'max_length', 'trim', 'mismatches')] == [400, 1000, 15, 2]
    catalog, catalog_receipt = provenance.load_catalog()
    spec = importlib.util.spec_from_file_location('frozen_mapping', PREVIOUS / 'map_assay.py')
    mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mapper)
    rows, references, synthetic = [], [], []
    args.output.mkdir(parents=True)
    (args.output / 'angiospermae-before.yaml').write_bytes(before)
    for accession in SOURCE_IDS:
        record = catalog['records'][accession]
        raw = (PREVIOUS / 'assay-v1' / f'{accession}.gb').read_bytes()
        text = raw.decode()
        sequence = ''.join(re.findall('[acgt]+', text.split('ORIGIN', 1)[1].split('//', 1)[0])).upper()
        assert sequence == record['sequence']
        genes = mapper.genes_from_gb(text)
        row = {'accession': accession, 'organism': record['organism'],
               'source_length': len(sequence), 'source_sequence_sha256': record['sequence_sha256'],
               'raw_genbank_sha256': digest(raw), 'raw_genbank_path': str((PREVIOUS / 'assay-v1' / f'{accession}.gb').relative_to(ROOT)),
               'boundary_annotations': [g for g in genes if g['gene'] in ('trnV-UAC', 'trnM-CAU', 'atpE', 'atpB')],
               'assays': {}}
        for label, reverse in [('legacy_S1022', OLD_REVERSE), ('corrected_6R', NEW_REVERSE),
                               ('rejected_internal_4R', 'GAACCAATGACTCCCGCCGTATG')]:
            row['assays'][label] = {
                'forward_full': FORWARD, 'reverse_full': reverse,
                'full_pairs': mapper.map_record(sequence, genes, FORWARD, reverse),
                'effective_pairs': mapper.map_record(sequence, genes, FORWARD[-15:], reverse[-15:])}
        corrected = [p for p in row['assays']['corrected_6R']['effective_pairs'] if p['within_400_1000_bound']]
        legacy = [p for p in row['assays']['legacy_S1022']['effective_pairs'] if p['within_400_1000_bound']]
        assert len(corrected) == 1 and not legacy
        pair = corrected[0]
        assert pair['strand'] == '+' and not pair['wraps_origin']
        start, end = pair['forward']['start'], pair['reverse_binding']['end']
        region = sequence[start:end]
        assert digest(region.encode()) == pair['amplicon_sha256']
        row['corrected_reference_interval'] = {'start': start, 'end': end, 'strand': '+', 'length': len(region), 'sha256': digest(region.encode())}
        if accession != 'MN990625.1':
            reference = {'accession': accession, 'taxon': record['organism'], 'sequence': region,
                         'provenance': {'schema_version': 1, 'kind': 'public_record_region',
                                        'accession_version': accession, 'start': start, 'end': end,
                                        'strand': '+', 'wraps_origin': False,
                                        'source_length': len(sequence), 'source_sequence_sha256': record['sequence_sha256'],
                                        'sequence_sha256': digest(region.encode())}}
            provenance.verify_reference(reference, catalog)
            references.append(reference)
            old_pair = next(p for p in row['assays']['legacy_S1022']['full_pairs'] if p['strand'] == '+' and not p['wraps_origin'])
            synth_start = old_pair['forward']['start'] - 100
            synth_end = old_pair['reverse_binding']['end'] + 100
            source_region = sequence[synth_start:synth_end]
            fastq = args.output / f'{accession}.synthetic.fastq'
            count = 0
            with fastq.open('x') as stream:
                # All 151-base windows, both orientations. Every source 19-mer
                # occurs at least twice, including the two outermost windows.
                for offset in range(len(source_region) - 150):
                    read = source_region[offset:offset + 151]
                    for strand, oriented in [('+', read), ('-', mapper.rc(read))]:
                        stream.write(f'@{accession}:{synth_start + offset}:{strand}\n{oriented}\n+\n' + 'I' * len(oriented) + '\n')
                        count += 1
            expected = args.output / f'{accession}.expected-corrected.fasta'
            expected.write_text(f'>{accession}|trnV-atpE|source={start}:{end}:+\n{region}\n')
            synthetic.append({'source_accession': accession, 'source_start': synth_start, 'source_end': synth_end,
                              'source_interval_sha256': digest(source_region.encode()), 'fastq': fastq.name,
                              'fastq_sha256': digest(fastq.read_bytes()), 'read_count': count, 'read_length': 151,
                              'quality': 40, 'expected_corrected_fasta': expected.name,
                              'expected_corrected_sha256': digest(region.encode()),
                              'scope': 'Error-free public-source synthetic reads; not real-sample validation.'})
        rows.append(row)
    assert rows[0]['corrected_reference_interval']['sha256'] == rows[2]['corrected_reference_interval']['sha256']
    after = before.decode().replace('panel_version: 1.0.2', 'panel_version: 2.0.0', 1)
    after = after.replace('changelog:\n', 'changelog:\n- panel_version: 2.0.0\n  date: 2026-10-10\n  sharkmer_version: "3.2.0-dev"\n  changes: "Replace trnV-atpE reverse primer S1022 with published atpE/trnV_6R and add exact public references. Target name is retained, but its sequence boundaries change incompatibly with panel 1.x. The extended marker includes the complete atpE gene and ends at the overlapping atpE/atpB boundary. Other primers and assay bounds are unchanged."\n', 1)
    after = after.replace(f'reverse_seq: "{OLD_REVERSE}"', f'reverse_seq: "{NEW_REVERSE}"', 1)
    old_metadata = '  citation: "Kress et al 2005. https://doi.org/10.1073/pnas.0503123102"\n  notes: "trnV5f; S1022"'
    new_metadata = '  citation: "Kress et al 2005 (forward/locus). https://doi.org/10.1073/pnas.0503123102; Tank and Donoghue 2010 Appendix 2 (atpE/trnV_6R, attributed to Bremer et al 2002). https://donoghuelab.yale.edu/sites/default/files/190_tank_systbot10.pdf"\n  notes: "Legacy trnV5f forward; atpE/trnV_6R reverse. Extended trnV-to-atpE marker includes trnM and complete atpE; reverse 3-prime end overlaps the atpE start and its 5-prime portion lies in overlapping atpB. Panel 2.0.0 changes this assay boundary; 400-1000 bp bounds, trim 15 and two mismatches per primer remain unchanged."'
    assert old_metadata in after
    after = after.replace(old_metadata, new_metadata, 1)
    references_yaml = yaml.safe_dump([{'gene': 'trnV-atpE', 'sequences': references}], sort_keys=False, width=100000)
    after = after.replace('references:\n', 'references:\n' + references_yaml, 1)
    (args.output / 'angiospermae-after.yaml').write_text(after)
    after_panel = yaml.safe_load(after)
    for old, new in zip(panel['primers'], after_panel['primers']):
        if old['gene'] != 'trnV-atpE':
            assert old == new
    assert after_panel['validation'] == panel['validation']
    assert not provenance.audit_references(after_panel)['excluded']
    receipt = {'schema_version': 1, 'coordinate_system': 'zero_based_half_open_forward_source',
               'producer_sha256': digest(Path(__file__).read_bytes()),
               'mapping_helper_sha256': digest((PREVIOUS / 'map_assay.py').read_bytes()),
               'catalog': catalog_receipt, 'before_panel_sha256': digest(before),
               'after_panel_sha256': digest(after.encode()), 'mismatches_each_primer': 2,
               'effective_primer_policy': '3-prime 15 bases; k19, min(trim,k-1)', 'length_bounds': [400, 1000],
               'records': rows, 'references': references, 'synthetic_inputs': synthetic,
               'duplicate_reference_omitted': 'MN990625.1 corrected interval is identical to NC_008326.1',
               'interpretation': 'Published 6R replaces the internal atpB S1022 endpoint. The extended marker contains complete atpE and trnM; 6R overlaps the atpE start by 3 bp and atpB by all 21 bp. This is not a pure intergenic-spacer assay.',
               'literature': [{'url': 'https://donoghuelab.yale.edu/sites/default/files/190_tank_systbot10.pdf',
                               'location': 'Appendix 2, printed page 441; atpE/trnV IGS and atpB gene entries',
                               'facts': ['atpE/trnV_6R is GGAGAGCAATTTGAAGAAATG, attributed to Bremer et al. 2002.',
                                         'S1022 is CGACATTTGCACATTTAGATGCTAC, listed among atpB primers.']},
                              {'url': 'https://pubmed.ncbi.nlm.nih.gov/12144762/',
                               'doi': '10.1016/S1055-7903(02)00240-3',
                               'facts': ['Bremer et al. describe the trnV, trnM and intervening region toward atpE.']}],
               'limits': ['Hamming source-coordinate mapping is not wet-lab validation.',
                          'Related-species Acer opalus is a reference proxy for Acer monspessulanum.',
                          'Public-source truth does not establish a real sample haplotype.',
                          'The original Kress supplemental oligo table has not been independently verified.']}
    write_json(args.output / 'mapping.json', receipt)
    write_json(args.output / 'checksums.json', {p.name: digest(p.read_bytes()) for p in sorted(args.output.iterdir()) if p.is_file()})
    print(json.dumps({'output': str(args.output), 'records': len(rows), 'references': len(references), 'synthetic_reads': [s['read_count'] for s in synthetic]}))


if __name__ == '__main__':
    main()

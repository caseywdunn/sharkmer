#!/usr/bin/env python3
"""Independent proposal verification; no biological input scan."""
import hashlib
import importlib.util
import json
from pathlib import Path
import yaml

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
HELPER=ROOT/'benchmarks/benchmark_results/v3.2-target-coverage-20261009/review/check_references.py'
spec=importlib.util.spec_from_file_location('independent_mapping_review',HELPER)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def load(p):return json.loads(p.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()

def main():
    directory=HERE/'proposal-v1';report=load(directory/'mapping.json')
    for name,digest in load(directory/'checksums.json').items():assert sha((directory/name).read_bytes())==digest
    before=yaml.safe_load((directory/'angiospermae-before.yaml').read_text())
    after=yaml.safe_load((directory/'angiospermae-after.yaml').read_text())
    assert before['panel_version']=='1.0.2' and after['panel_version']=='2.0.0'
    assert len(before['primers'])==len(after['primers'])==10
    for old,new in zip(before['primers'],after['primers']):
        if old['gene']=='trnV-atpE':
            assert new['reverse_seq']=='GGAGAGCAATTTGAAGAAATG'
            assert {k:v for k,v in old.items() if k not in ('reverse_seq','citation','notes')}=={k:v for k,v in new.items() if k not in ('reverse_seq','citation','notes')}
        else:assert old==new
    assert after['validation']==before['validation']
    assert after['references'][1:]==before['references']
    assert after['references'][0]=={'gene':'trnV-atpE','sequences':report['references']}
    checks=[];sources={}
    for record in report['records']:
        path=ROOT/record['raw_genbank_path'];raw=path.read_bytes()
        assert sha(raw)==record['raw_genbank_sha256']
        text=raw.decode();sequence=''.join(audit.re.findall('[acgt]+',text.split('ORIGIN',1)[1].split('//',1)[0])).upper()
        assert sha(sequence.encode())==record['source_sequence_sha256']
        sources[record['accession']]=sequence
        features=audit.features(text)
        assert record['boundary_annotations']==[g for g in features if g['gene'] in ('trnV-UAC','trnM-CAU','atpE','atpB')]
        for label,assay in record['assays'].items():
            for mode in ('full','effective'):
                fwd=assay['forward_full'];rev=assay['reverse_full']
                if mode=='effective':fwd,rev=fwd[-15:],rev[-15:]
                expected=[];n=len(sequence)
                for strand,s in [('+',sequence),('-',audit.rc(sequence))]:
                    left=audit.brute_hamming(s+s[:len(fwd)-1],fwd);right=audit.brute_hamming(s+s,audit.rc(rev))
                    for a,fd in left:
                        for z,rd in right:
                            length=z+len(rev)-a
                            if z>=a+len(fwd) and length<=n:expected.append((strand,a,z,length,tuple(fd),tuple(rd)))
                pairs=assay[mode+'_pairs']
                actual=[(p['strand'],p['forward']['start'],p['reverse_binding']['start'],p['amplicon_length'],tuple(p['forward']['mismatch_positions']),tuple(p['reverse_binding']['mismatch_positions'])) for p in pairs]
                assert sorted(actual)==sorted(expected)
                eligible=[p for p in pairs if 400<=p['amplicon_length']<=1000]
                assert all(p['within_400_1000_bound']==(400<=p['amplicon_length']<=1000) for p in pairs)
                assert len(eligible)==(1 if label=='corrected_6R' else 0)
                checks.append({'accession':record['accession'],'assay':label,'primer_mode':mode,'all_pairs':len(pairs),'eligible_lengths':[p['amplicon_length'] for p in eligible]})
        interval=record['corrected_reference_interval'];fragment=sequence[interval['start']:interval['end']]
        assert len(fragment)==interval['length'] and sha(fragment.encode())==interval['sha256']
    for reference in report['references']:
        p=reference['provenance'];s=sources[reference['accession']]
        assert reference['sequence']==s[p['start']:p['end']] and sha(reference['sequence'].encode())==p['sequence_sha256']
        assert sha(s.encode())==p['source_sequence_sha256'] and len(s)==p['source_length']
    for synthetic in report['synthetic_inputs']:
        s=sources[synthetic['source_accession']][synthetic['source_start']:synthetic['source_end']]
        assert sha(s.encode())==synthetic['source_interval_sha256']
        data=(directory/synthetic['fastq']).read_bytes();assert sha(data)==synthetic['fastq_sha256']
        lines=data.decode().splitlines();assert len(lines)//4==synthetic['read_count']
        expected=[word for offset in range(len(s)-150) for word in (s[offset:offset+151],audit.rc(s[offset:offset+151]))]
        assert lines[1::4]==expected and all(q=='I'*151 for q in lines[3::4]) and all(p=='+' for p in lines[2::4])
        expected_fasta=(directory/synthetic['expected_corrected_fasta']).read_text().splitlines()[1]
        assert sha(expected_fasta.encode())==synthetic['expected_corrected_sha256']
    result={'status':'approved_for_bounded_validation','mapping_checks':checks,'synthetic_inputs_verified':2,'unique_references_verified':2,
            'interpretation':'The corrected published6R primer has one eligible public-source pair per plastome; full physical endpoints include complete trnM and atpE. Its3-prime terminal3bp overlap atpE; all21bp overlap atpB, including18bp exclusively atpB. LegacyS1022 and rejected4R have no eligible400–1000bp pair.',
            'conditions':['Compare emitted corrected synthetic sequences literally to frozen909/929bp expected references; standard vs uncertain outcomes must be explicit.',
                          'Legacy/no-template negatives do not authorize suppressing any bounded uncertainty; classify unexpected standard products rather than treating process success as target validation.',
                          'Compare both real plant samples before/after with full panels; preserve every other target identity/status/count and distinguish intentional assay-boundary changes from retention regression.',
                          'Use legacy locus citation wording for Kress2005 unless original supplemental primer table is independently verified.'],
            'artifacts':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in [directory/'mapping.json',directory/'checksums.json',directory/'angiospermae-after.yaml',Path(__file__),HELPER]}}
    with (HERE/'independent_proposal_review.json').open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'mapping_checks':len(checks),'synthetic_inputs_verified':2}))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Independent literal/brute checks of the bounded annotation assessments."""
import gzip
import hashlib
import json
from pathlib import Path
import random
import re
import map_assay

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def brute(sequence,primer):
    found=[]
    for start in range(len(sequence)-len(primer)+1):
        mismatches=0
        for a,b in zip(sequence[start:start+len(primer)],primer):
            mismatches+=a!=b
            if mismatches>2:break
        if mismatches<=2 and set(sequence[start:start+len(primer)])<=set('ACGT'):found.append((start,mismatches))
    return found


def main():
    rng=random.Random(23079);controls=0
    for n in [15,22,25]:
        for _ in range(20):
            sequence=''.join(rng.choices('ACGTN',k=150));primer=''.join(rng.choices('ACGT',k=n))
            selftest=sequence+primer+map_assay.rc(primer)
            assert [(r['start'],r['mismatches']) for r in map_assay.hits(selftest,primer)]==brute(selftest,primer)
            controls+=1
    with gzip.open(ROOT/'panels/reference_sources.json.gz','rt') as f:catalog=json.load(f)
    mapping=json.load(open(HERE/'assay-v1/mapping.json'));annotation=[]
    for record in mapping['records']:
        source=catalog['records'][record['accession']]
        assert source['sequence_sha256']==record['source_sha256']
        raw=HERE/'assay-v1'/(record['accession']+'.gb');assert digest(raw)==record['raw_genbank_sha256']
        # Independent full-source unseeded Hamming scan for all orientations.
        for kind,forward,reverse in [('full_primer_pairs',mapping['forward_full'],mapping['reverse_full']),
                                     ('effective_primer_pairs',mapping['effective_forward'],mapping['effective_reverse'])]:
            pairs=[]
            for strand,sequence in [('+',source['sequence']),('-',map_assay.rc(source['sequence']))]:
                n=len(sequence);doubled=sequence+sequence
                starts=[x for x in brute(sequence+sequence[:len(forward)-1],forward) if x[0]<n]
                ends=[x for x in brute(sequence+sequence[:len(reverse)-1],map_assay.rc(reverse)) if x[0]<n]
                for start,fm in starts:
                    for end,rm in ends:
                        for actual_end in [end,end+n]:
                            span=actual_end+len(reverse)-start
                            if actual_end>=start+len(forward) and span<=n:pairs.append((strand,start,actual_end,span,fm,rm))
            produced=[(p['strand'],p['forward']['start'],p['reverse_binding']['start'],p['amplicon_length'],p['forward']['mismatches'],p['reverse_binding']['mismatches']) for p in record[kind]]
            assert sorted(pairs)==sorted(produced)
            assert len(pairs)==1 and pairs[0][3]>1000
            annotation.append({'accession':record['accession'],'assay':kind,'independent_unseeded_pairs':pairs})
    report=json.load(open(HERE/'assessment-v2/report.json'))
    for path,expected in report['source_receipts'].items():assert digest(Path(path))==expected
    exact_rows=[]
    for candidate in report['candidates']:
        if candidate['exact_full_sequence_placements']:
            seq=next(c['sequence'] for c in json.load(open(HERE.parent/'matrix/novel_target_inventory.json')) if c['sha256']==candidate['sha256'])
            for placement in candidate['exact_full_sequence_placements']:
                public=catalog['records'][placement['accession']]['sequence'][placement['source_start']:placement['source_end']]
                assert (public if placement['strand']=='+' else map_assay.rc(public))==seq
                assert placement['annotated_target_boundaries_agree']
                assert {g['gene'] for g in placement['left_endpoint_genes']}=={'atpB'}
                assert {g['gene'] for g in placement['right_endpoint_genes']}=={'rbcL'}
            assert candidate['annotated_exact_target_support'] and candidate['full_sequence_reference_support']=='exact_full_length_same_taxon_public_sequence'
            exact_rows.append(candidate['sha256'])
    assert len(report['candidates'])==154 and len(exact_rows)==1
    result={'status':'pass','randomized_optimized_vs_brute_mapping_controls':controls,'public_source_assay_pair_checks':annotation,
            'candidate_rows_checked':154,'exact_full_annotated_target_candidates':exact_rows,
            'sources_verified':report['source_receipts'],
            'output_hashes':{str(p):digest(p) for p in [HERE/'assessment-v2/report.json',HERE/'assessment-v2/reference_overlay.json',HERE/'assessment-v2/per_candidate.tsv',HERE/'assay-v1/mapping.json']},
            'verifier_sha256':digest(Path(__file__))}
    with (HERE/'assessment-v2/verification.json').open('x') as f:f.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'pass','randomized_controls':controls,'source_assay_cases':len(annotation),'candidate_rows_checked':154,'exact_full_target_candidates':len(exact_rows)}))


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Frozen public plastome feature/primer mapping; no assay parameter changes."""
import hashlib
import json
import re
import tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
ARCHIVE=ROOT/'benchmarks/benchmark_results/v3.2-public-replacements-20260913'
FWD='GTGTAAACGAGTTGCTCTACCA'
REV='CGACATTTGCACATTTAGATGCTAC'


def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def hits(sequence,primer,mismatches=2):
    # Any <=2-mismatch alignment has at least one of three disjoint exact seeds.
    candidates=set();n=len(primer);segments=[(n*i//3,n*(i+1)//3) for i in range(3)]
    for a,b in segments:
        offset=0
        while True:
            index=sequence.find(primer[a:b],offset)
            if index<0:break
            if 0<=index-a<=len(sequence)-n:candidates.add(index-a)
            offset=index+1
    result=[]
    for start in sorted(candidates):
        observed=sequence[start:start+n]
        differences=[i for i,(a,b) in enumerate(zip(observed,primer)) if a!=b]
        if len(differences)<=mismatches and set(observed)<=set('ACGT'):
            result.append({'start':start,'end':start+n,'observed':observed,'mismatches':len(differences),'mismatch_positions':differences})
    return result


def genes_from_gb(text):
    features=[]
    for block in re.split(r'\n(?=     \S)',text.split('FEATURES',1)[1].split('ORIGIN',1)[0]):
        lines=block.splitlines()
        if not lines or not lines[0].startswith('     gene '):continue
        gene=re.search(r'/gene="([^"]+)"',block)
        if not gene:continue
        location=lines[0][21:].strip()
        coords=[int(v) for v in re.findall(r'\d+',location)]
        if len(coords)==2:
            features.append({'gene':gene.group(1),'start':min(coords)-1,'end':max(coords),'strand':'-' if 'complement' in location else '+','genbank_location':location})
    return features


def map_record(sequence,features,forward,reverse):
    n=len(sequence);rows=[]
    for strand,oriented in [('+',sequence),('-',rc(sequence))]:
        # Explicit circular wrap support, all starts exactly once per orientation.
        doubled=oriented+oriented
        starts=[h for h in hits(doubled,forward) if h['start']<n]
        ends=hits(doubled,rc(reverse))
        for f in starts:
            for r in ends:
                length=r['end']-f['start']
                if r['start']<f['end'] or length>n:continue
                def original(a,b):
                    return [a%n,b%n] if strand=='+' else [(n-b)%n,(n-a)%n]
                def overlaps(h):
                    a,b=original(h['start'],h['end'])
                    return [g for g in features if a<b and a<g['end'] and b>g['start']]
                rows.append({'strand':strand,'forward':{**f,'source_span':original(f['start'],f['end']),'genes':overlaps(f)},
                             'reverse_binding':{**r,'source_span':original(r['start'],r['end']),'genes':overlaps(r)},
                             'amplicon_length':length,'within_400_1000_bound':400<=length<=1000,
                             'wraps_origin':r['end']>n,'amplicon_sha256':hashlib.sha256(doubled[f['start']:r['end']].encode()).hexdigest()})
    return rows


def main():
    out=HERE/'assay-v1';out.mkdir()
    members={m['path']:m for m in json.load(open(ARCHIVE/'ARCHIVE_CONTENTS.json'))['members']}
    rows=[]
    with tarfile.open(ARCHIVE/'evidence.tar.gz') as archive:
        for short in ['NC_008326','NC_056221','MN990625']:
            name=f'plants_nematode/raw/genbank/{short}.gb';data=archive.extractfile(name).read();expected=members[name]
            assert hashlib.sha256(data).hexdigest()==expected['sha256'] and len(data)==expected['size_bytes']
            text=data.decode();sequence=''.join(re.findall('[acgt]+',text.split('ORIGIN',1)[1].split('//',1)[0])).upper()
            accession=re.search(r'^VERSION\s+(\S+)',text,re.M).group(1)
            organism=re.search(r'^  ORGANISM\s+(.+)',text,re.M).group(1)
            (out/f'{accession}.gb').write_bytes(data)
            genes=genes_from_gb(text)
            rows.append({'accession':accession,'organism':organism,'source_length':len(sequence),'source_sha256':hashlib.sha256(sequence.encode()).hexdigest(),
                         'raw_genbank_sha256':expected['sha256'],'archive_member':name,
                         'boundary_gene_annotations':[g for g in genes if g['gene'] in ['trnV-UAC','atpE','atpB']],
                         'full_primer_pairs':map_record(sequence,genes,FWD,REV),
                         'effective_primer_pairs':map_record(sequence,genes,FWD[-15:],REV[-15:])})
    report={'schema_version':1,'status':'completed','coordinate_system':'zero_based_half_open_forward_source_unless_strand_specified',
            'forward_full':FWD,'reverse_full':REV,'effective_forward':FWD[-15:],'effective_reverse':REV[-15:],
            'mismatches_each_primer':2,'trim_bases':15,'kmer_length':19,'configured_length_bounds':[400,1000],
            'records':rows,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact source-coordinate assay mapping under Hamming mismatches only; this is not a Sharkmer graph run or sample callability test.'}
    with (out/'mapping.json').open('x') as stream:stream.write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps([{ 'accession':r['accession'],'full_pairs':len(r['full_primer_pairs']),
                       'effective_pairs':len(r['effective_primer_pairs']),
                       'eligible_effective_pairs':sum(x['within_400_1000_bound'] for x in r['effective_primer_pairs'])} for r in rows]))


if __name__=='__main__':main()

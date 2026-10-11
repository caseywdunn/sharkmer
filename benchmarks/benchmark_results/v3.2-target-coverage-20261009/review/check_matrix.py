#!/usr/bin/env python3
"""Independent panel, receipt, archived FASTA and identity-set reconciliation."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import tarfile

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parents[2]

def digest(data):
    return hashlib.sha256(data).hexdigest()

def load(path):
    return json.loads(path.read_text())

def fasta(data):
    output=[]
    header=None
    sequence=[]
    for line in data.decode().splitlines():
        if line.startswith('>'):
            if header is not None:
                output.append((header,''.join(sequence)))
            header=line[1:]
            sequence=[]
        elif line.strip():
            assert header is not None
            sequence.append(line.strip())
    if header is not None:
        output.append((header,''.join(sequence)))
    return output

def key(row):
    return tuple(row[k] for k in ('accession','gene','sha256'))

def panel_genes(path):
    body=path.read_text().split('\nprimers:\n',1)[1]
    body=re.split(r'(?m)^[A-Za-z_][^:\n]*:',body,maxsplit=1)[0]
    rows=[]
    for block in re.split(r'(?m)^- gene: ',body)[1:]:
        gene=block.splitlines()[0].strip().strip('"\'')
        index=re.search(r'(?m)^  index: (\S+)',block)
        region=re.search(r'(?m)^  region: (\S+)',block)
        if region:
            gene+='-'+region.group(1).strip('"\'')
        elif index:
            gene+='_'+index.group(1).strip('"\'')
        rows.append(path.stem+'_'+gene)
    assert len(rows)==len(set(rows))
    return set(rows)

def main():
    matrix=BASE/'matrix'
    validation=load(matrix/'validation.json')
    for source,expected in validation['source_files'].items():
        data=Path(source).read_bytes()
        assert digest(data)==expected['sha256'] and len(data)==expected['size_bytes'],source
    for name,expected in validation['outputs'].items():
        data=(matrix/name).read_bytes()
        assert digest(data)==expected['sha256'] and len(data)==expected['size_bytes'],name
    current=load(matrix/'sequence_inventory.json')
    historical=load(matrix/'historical_inventory.json')
    novel=load(matrix/'novel_target_inventory.json')
    rows=load(matrix/'target_matrix.json')
    inventories=[current,historical,novel]
    for inv in inventories:
        assert len(inv)==len({key(r) for r in inv})
        for r in inv:
            assert digest(r['sequence'].encode())==r['sha256']
            assert len(r['sequence'])==r['length']
    current_map={key(r):r for r in current}
    old_map={key(r):r for r in historical}
    row_map={(r['accession'],r['gene']):r for r in rows}
    assert len(row_map)==len(rows)
    observed_current={}
    observed_rows=set()
    status_counts=Counter()
    zero_reasons=Counter()
    pcr_outcomes=[]
    for run in validation['current_runs']:
        accession,panel=run['accession'],run['panel']
        receipt=load(Path(run['receipt']))
        assert receipt['status']==run['status']=='pass'
        assert receipt['prefix_before']==run['prefix']
        expected_genes=panel_genes(ROOT/'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/panels'/(panel+'.yaml'))
        pcr=receipt['outputs']['pcr_results']
        assert {p['gene_name'] for p in pcr}==expected_genes
        assert len(pcr)==len(expected_genes)
        fasta_records=[]
        for name,info in receipt['outputs']['files'].items():
            if not name.endswith('.fasta'):
                continue
            data=Path(info['path']).read_bytes()
            assert digest(data)==info['sha256'] and len(data)==info['size_bytes']
            fasta_records.extend(fasta(data))
        actual_records=receipt['outputs']['records']
        assert Counter(fasta_records)==Counter((r['header'],r['sequence']) for r in actual_records)
        for r in actual_records:
            k=(accession,r['gene'],r['sha256'])
            assert k not in observed_current
            observed_current[k]=r
            inv=current_map[k]
            assert all(inv[f]==r[f] for f in ('header','sequence','length','assembly_status','scope'))
            assert inv['source']==run['receipt']
        for p in pcr:
            k=(accession,p['gene_name'])
            observed_rows.add(k)
            row=row_map[k]
            selected=[r for r in actual_records if r['gene']==p['gene_name']]
            counts=Counter(r['assembly_status'] for r in selected)
            assert len(selected)==p['n_products']==row['current_total']
            assert counts['standard']==p['n_standard_products']==row['current_standard']
            assert counts['uncertain']==p['n_uncertain_products']==row['current_uncertain']
            assert row['attempted'] and row['run_status']=='pass'
            assert row['zero_output_current']==(not selected)
            status_counts[p['status']]+=1
            if not selected:
                zero_reasons[p.get('failure_reason')]+=1
            pcr_outcomes.append({'accession':accession,'gene':p['gene_name'],'target_status':p['status'],
                                 'failure_reason':p.get('failure_reason'),'current_total':len(selected)})
    assert observed_rows==set(row_map)
    assert set(observed_current)==set(current_map)
    authenticated=load(ROOT/'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/historical_identities.json')
    old_expectations={}
    for sample in authenticated['sample_expectations']:
        for gene,records in sample['roles']['baseline']['by_gene'].items():
            for r in records:
                old_expectations[(sample['accession'],sample['panel']+'_'+gene,r['sequence_sha256'])]=r
    assert set(old_expectations)==set(old_map)
    archives={}
    try:
        for source,expected in validation['archive_members'].items():
            archive,member=source.split('::',1)
            if archive not in archives:
                archives[archive]=tarfile.open(archive,'r:gz')
            data=archives[archive].extractfile(member).read()
            assert digest(data)==expected['sha256'] and len(data)==expected['size_bytes']
            emitted=fasta(data)
            inventory_records=[r for r in historical if r['source']==source]
            assert Counter(emitted)==Counter((r['header'],r['sequence']) for r in inventory_records)
        for k,r in old_map.items():
            e=old_expectations[k]
            assert r['source'] in validation['archive_members']
            assert r['sequence']==e['sequence'] and r['header']==e['description']
            assert r['source'].split('::',1)[1]==e['archived_member']
    finally:
        for archive in archives.values():
            archive.close()
    for k,r in current_map.items():
        retained=k in old_map
        assert r['historical_retained']==retained and r['new_vs_release']==(not retained)
        assert r['newly_represented_target']==row_map[k[:2]]['newly_represented_target']
    for k,r in old_map.items():
        assert r['retained_in_current']==(k in current_map)
    for k,r in row_map.items():
        old={identity[2] for identity in old_map if identity[:2]==k}
        new={identity[2] for identity in current_map if identity[:2]==k}
        for name,values in (('v3_1',old),('current',new),('retained',old&new),('new',new-old),('lost',old-new)):
            assert set(r[name+'_identities'])==values
        assert r['v3_1_products']==len(old)
        for name,values in (('retained',old&new),('new',new-old),('lost',old-new)):
            assert r[name+'_count']==len(values)
        assert r['newly_represented_target']==(not old and bool(new))
        assert r['whole_product_validation']=='not_established'
    assert {key(r) for r in novel}=={key(r) for r in current if r['newly_represented_target']}
    for name,inv in zip(('sequence_inventory','historical_inventory','novel_target_inventory'),inventories):
        emitted=fasta((matrix/(name+'.fasta')).read_bytes())
        assert Counter((h.split()[0],s) for h,s in emitted)==Counter(('|'.join(key(r)),r['sequence']) for r in inv)
    counts={'target_rows':len(rows),'historical_products':len(historical),'current_products':len(current),
            'retained_exact_products':len(set(old_map)&set(current_map)),'lost_exact_products':len(set(old_map)-set(current_map)),
            'new_exact_products':len(set(current_map)-set(old_map)),
            'standard_products':sum(r['assembly_status']=='standard' for r in current),
            'uncertain_products':sum(r['assembly_status']=='uncertain' for r in current),
            'historical_represented_targets':sum(bool(r['v3_1_products']) for r in rows),
            'current_represented_targets':sum(bool(r['current_total']) for r in rows),
            'newly_represented_targets':sum(r['newly_represented_target'] for r in rows),
            'novel_target_products':len(novel),'zero_output_current_targets':sum(r['zero_output_current'] for r in rows)}
    assert counts==validation['checks']
    result={'status':'pass','counts':counts,'source_files_verified':len(validation['source_files']),
            'authenticated_archived_fasta_members_verified':len(validation['archive_members']),
            'actual_per_target_status_counts':dict(status_counts),'zero_output_failure_reasons':dict(zero_reasons),
            'target_outcomes':pcr_outcomes,'validation_sha256':digest((matrix/'validation.json').read_bytes()),
            'reviewer_source_sha256':digest(Path(__file__).read_bytes()),
            'scope':'Exactly ten matched samples and their complete five frozen panels; 118 attempted sample-target rows. Literal oriented full-sequence identities, with reverse complements retained separately. Not an overlap claim for the three absent bacterial/metagenomic samples.',
            'limits':['Run status pass denotes command completion. All 44 zero-output targets have actual target status fail and bounded-search/connectivity reasons; zeros do not establish biological absence.',
                      'Sequence retention and emission counts do not establish reference identity, local read support or whole-product validity.']}
    with (HERE/'matrix_review.json').open('x') as out:
        out.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','counts','actual_per_target_status_counts','zero_output_failure_reasons')}))

if __name__=='__main__':
    main()

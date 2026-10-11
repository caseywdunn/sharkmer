#!/usr/bin/env python3
"""Independent frozen-query, HSP, feature and effective-primer audit."""
from collections import Counter,defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import yaml

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parents[2]
OLD=ROOT/'benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/blast-gains-v3'
NEW=BASE/'reference/blast-missing-v1'
COMP=str.maketrans('ACGT','TGCA')

def load(path): return json.loads(path.read_text())
def sha(data): return hashlib.sha256(data).hexdigest()
def rc(sequence): return sequence.translate(COMP)[::-1]
def key(row): return tuple(row[k] for k in ('accession','gene','sha256'))

def fasta(path):
    result={}
    for record in path.read_text().split('>')[1:]:
        lines=record.splitlines();identifier=lines[0].split()[0]
        assert identifier not in result
        result[identifier]=''.join(lines[1:])
    return result

def raw_hits(path):
    groups=defaultdict(list)
    for row in csv.reader(path.open(),delimiter='\t'):
        if row:
            groups[row[0]].append({'query_id':row[0],'subject_id':row[1],'identity_pct':float(row[2]),
                                  'length':int(row[3]),'query_length':int(row[4]),'qstart':int(row[5]),'qend':int(row[6]),
                                  'sstart':int(row[8]),'send':int(row[9])})
    return groups

def brute_hamming(sequence,primer):
    hits=[]
    for start in range(len(sequence)-len(primer)+1):
        differences=[]
        for i,base in enumerate(primer):
            if sequence[start+i]!=base:
                differences.append(i)
                if len(differences)>2:break
        else:
            hits.append((start,differences))
    return hits

def features(text):
    results=[]
    body=text.split('FEATURES',1)[1].split('ORIGIN',1)[0]
    for match in re.finditer(r'(?ms)^     gene\s+([^\n]+)\n(.*?)(?=^     \S|\Z)',body):
        name=re.search(r'/gene="([^"]+)"',match.group(2))
        span=re.fullmatch(r'(complement\()?([0-9]+)\.\.([0-9]+)\)?',match.group(1).strip())
        if name and span:
            results.append({'gene':name.group(1),'start':int(span.group(2))-1,'end':int(span.group(3)),
                            'strand':'-' if span.group(1) else '+','genbank_location':match.group(1).strip()})
    return results

def main():
    report_path=BASE/'reference/assessment-v2/report.json'
    report=load(report_path)
    inventory=load(BASE/'matrix/novel_target_inventory.json')
    inventories={key(r):r for r in inventory}
    rows={key(r):r for r in report['candidates']}
    assert len(rows)==len(inventories)==154 and set(rows)==set(inventories)
    for path,digest in report['source_receipts'].items():assert sha(Path(path).read_bytes())==digest,path
    execution=load(NEW/'execution.json')
    assert execution['status']=='pass' and execution['sources_unchanged']
    assert execution['source_hashes_before']==execution['source_hashes_after']
    for path,digest in execution['source_hashes_after'].items():assert sha(Path(path).read_bytes())==digest,path
    for name,receipt in execution['outputs'].items():assert sha((NEW/name).read_bytes())==receipt['sha256']
    old_summary=load(OLD/'summary.json')
    for receipt in [old_summary['query_fasta'],old_summary['source_fasta'],old_summary['source_hits']]+[r[k] for r in old_summary['panel'].values() for k in ('hits','references')]:
        assert sha((ROOT/receipt['path']).read_bytes())==receipt['sha256']
    catalog_path=ROOT/'panels/reference_sources.json.gz'
    assert sha(catalog_path.read_bytes())==old_summary['catalog']['sha256']
    catalog=json.loads(gzip.decompress(catalog_path.read_bytes()))['records']
    assert len(catalog)==121
    source_metadata=load(OLD/'public_sources_metadata.json')
    queries={'saved':fasta(ROOT/old_summary['query_fasta']['path']),'new':fasta(NEW/'public_queries.fasta')}
    new_map={key(r):r for r in load(NEW/'queries.json')}
    old_map={r['query_id']:r for r in load(OLD/'per_sequence.json')}
    panels={}
    panel_meta={}
    for panel in ('angiospermae','insecta'):
        current=yaml.safe_load((ROOT/'panels'/(panel+'.yaml')).read_text())
        frozen=yaml.safe_load((ROOT/'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/panels'/(panel+'.yaml')).read_text())
        assert current['primers']==frozen['primers']
        panels[panel]=frozen
        refs=[(block['gene'],r) for block in current['references'] for r in block['sequences']]
        text=''.join(f'>reference_{i:03d}\n{r["sequence"]}\n' for i,(_,r) in enumerate(refs))
        assert sha(text.encode())==old_summary['panel'][panel]['references']['sha256']
        panel_meta[panel]={f'reference_{i:03d}':{'gene':gene,**r} for i,(gene,r) in enumerate(refs)}
    hits={}
    for namespace,folder in [('saved',OLD),('new',NEW)]:
        for scope in ('public','angiospermae','insecta'):
            name=('public_sources_hits.tsv' if scope=='public' else scope+'_reference_hits.tsv') if namespace=='saved' else scope+'_hits.tsv'
            hits[namespace,scope]=raw_hits(folder/name)
    mapping_path=BASE/'reference/assay-v1/mapping.json'
    mapping=load(mapping_path)
    primer=next(p for p in panels['angiospermae']['primers'] if p['gene']=='trnV-atpE')
    assert [primer['min_length'],primer['max_length']]==mapping['configured_length_bounds']==[400,1000]
    assert primer['trim']==mapping['trim_bases']==15 and mapping['kmer_length']==19
    assert primer['mismatches']==mapping['mismatches_each_primer']==2
    assert primer['forward_seq']==mapping['forward_full'] and primer['reverse_seq']==mapping['reverse_full']
    assert primer['forward_seq'][-15:]==mapping['effective_forward'] and primer['reverse_seq'][-15:]==mapping['effective_reverse']
    source_features={}
    mapping_checks=[]
    archive_index=load(ROOT/'benchmarks/benchmark_results/v3.2-public-replacements-20260913/ARCHIVE_CONTENTS.json')
    archive_members={r['path']:r for r in archive_index['members']}
    for r in mapping['records']:
        path=BASE/'reference/assay-v1'/(r['accession']+'.gb')
        data=path.read_bytes();text=data.decode()
        assert sha(data)==r['raw_genbank_sha256']==archive_members[r['archive_member']]['sha256']
        sequence=''.join(re.findall('[acgt]+',text.split('ORIGIN',1)[1].split('//',1)[0])).upper()
        assert sha(sequence.encode())==r['source_sha256']==catalog[r['accession']]['sequence_sha256']
        assert sequence==catalog[r['accession']]['sequence']
        annotations=features(text);source_features[r['accession']]=annotations
        assert [g for g in annotations if g['gene'] in ('trnV-UAC','atpE','atpB')]==r['boundary_gene_annotations']
        for mode in ('full','effective'):
            fwd=mapping['forward_full'] if mode=='full' else mapping['effective_forward']
            rev=mapping['reverse_full'] if mode=='full' else mapping['effective_reverse']
            expected=[]
            n=len(sequence)
            for strand,oriented in [('+',sequence),('-',rc(sequence))]:
                # Exhaustive independent Hamming scan, with one circular traversal.
                left=brute_hamming(oriented+oriented[:len(fwd)-1],fwd)
                right=brute_hamming(oriented+oriented,rc(rev))
                for fs,fd in left:
                    for rs,rd in right:
                        length=rs+len(rev)-fs
                        if rs>=fs+len(fwd) and length<=n:
                            expected.append((strand,fs,rs,length,tuple(fd),tuple(rd)))
            observed=[(p['strand'],p['forward']['start'],p['reverse_binding']['start'],p['amplicon_length'],tuple(p['forward']['mismatch_positions']),tuple(p['reverse_binding']['mismatch_positions'])) for p in r[mode+'_primer_pairs']]
            assert sorted(expected)==sorted(observed)
            for p in r[mode+'_primer_pairs']:
                assert p['strand']=='+' and not p['wraps_origin']
                start,end=p['forward']['start'],p['reverse_binding']['end']
                assert sha(sequence[start:end].encode())==p['amplicon_sha256']
                assert p['within_400_1000_bound']==(400<=end-start<=1000)
                for side in ('forward','reverse_binding'):
                    h=p[side];a,z=h['source_span']
                    assert sequence[a:z]==h['observed']
                    assert h['genes']==[g for g in annotations if a<g['end'] and z>g['start']]
            mapping_checks.append({'accession':r['accession'],'primer_mode':mode,'pairs':len(expected),'lengths':[p['amplicon_length'] for p in r[mode+'_primer_pairs']]})
    used=Counter()
    positive=[]
    for k,row in rows.items():
        candidate=inventories[k];namespace=row['query_namespace'];qid=row['query_id'];used[namespace]+=1
        assert queries[namespace][qid]==candidate['sequence']
        assert sha(candidate['sequence'].encode())==candidate['sha256']
        assert key(old_map[qid])==k if namespace=='saved' else key(new_map[k])==k
        public=hits[namespace,'public'].get(qid,[]);panel=hits[namespace,candidate['panel']].get(qid,[])
        logical=re.sub(r'_\d+$','',candidate['target'])
        expected=[h for h in panel if re.sub(r'_\d+$','',panel_meta[candidate['panel']][h['subject_id']]['gene'])==logical]
        same_taxon=[h for h in public if source_metadata[h['subject_id']]['taxon']==candidate['taxon']]
        strong=lambda h:h['identity_pct']>=90 and (abs(h['qend']-h['qstart'])+1)*100>=90*h['query_length']
        assert row['public_hsp_count']==len(public) and row['panel_hsp_count']==len(panel)
        assert row['expected_target_hsp_count']==len(expected) and row['same_taxon_hsp_count']==len(same_taxon)
        assert row['strong_public_single_hsp_count']==sum(strong(h) for h in public)
        assert row['strong_expected_target_single_hsp_count']==sum(strong(h) for h in expected)
        exact=[accession for accession,r in catalog.items() if r['organism']==candidate['taxon'] and (candidate['sequence'] in r['sequence'] or rc(candidate['sequence']) in r['sequence'])]
        assert row['exact_full_sequence_same_taxon_sources']==exact
        for p in row['exact_full_sequence_placements']:
            seq=catalog[p['accession']]['sequence'][p['source_start']:p['source_end']]
            assert seq==(candidate['sequence'] if p['strand']=='+' else rc(candidate['sequence']))
            ann=source_features[p['accession']]
            left=[g for g in ann if g['start']<=p['source_start']<g['end']]
            right=[g for g in ann if g['start']<p['source_end']<=g['end']]
            assert p['left_endpoint_genes']==left and p['right_endpoint_genes']==right
            assert p['annotated_target_boundaries_agree']==(candidate['target']=='atpB-rbcL' and any(g['gene']=='atpB' for g in left) and any(g['gene']=='rbcL' for g in right))
        annotated=any(p['annotated_target_boundaries_agree'] for p in row['exact_full_sequence_placements'])
        assert row['annotated_exact_target_support']==annotated
        if annotated:
            assert row['target_identity']=='supported_by_exact_public_sequence_and_independent_feature_boundaries'
            assert row['full_sequence_reference_support']=='exact_full_length_same_taxon_public_sequence'
            positive.append({'accession':candidate['accession'],'gene':candidate['gene'],'sha256':candidate['sha256'],'placements':row['exact_full_sequence_placements']})
        else:
            assert row['target_identity']=='not_established' and not any(strong(h) for h in public+panel)
        if candidate['target']=='trnV-atpE':
            assert not same_taxon and not row['trnV_mapped_assay_overlap_hsps']
            spans={r['accession']:(r['effective_primer_pairs'][0]['forward']['start'],r['effective_primer_pairs'][0]['reverse_binding']['end']) for r in mapping['records']}
            for h in public:
                accession=source_metadata[h['subject_id']]['accession']
                if accession in spans:
                    a,z=spans[accession];hs,he=min(h['sstart'],h['send'])-1,max(h['sstart'],h['send'])
                    assert he<=a or hs>=z
        assert row['whole_product_biological_validation']=='not_established'
    assert used=={'saved':113,'new':41} and len(positive)==1
    overlay=load(BASE/'reference/assessment-v2/reference_overlay.json')
    assert {key(r) for r in overlay['sequence_evidence']}==set(rows)
    for r in overlay['sequence_evidence']:
        row=rows[key(r)]
        assert r['reference_evidence_status']==row['reference_evidence_status']
        assert all(row[k]==v for k,v in r['reference_evidence_details'].items())
    result={'status':'pass','candidate_scope':154,'query_reuse':dict(used),'mapping_checks':mapping_checks,
            'annotated_exact_target_support':positive,
            'source_review':'Production src/pcr/primers.rs:191 and :218 retain the last min(trim,k-1) bases; k19 and trim15 use the checked 15-base primers. Coordinate mapping uses substitutions only and does not re-run graph search.',
            'conclusions':['One 812 bp Liriodendron atpB-rbcL candidate matches two same-species public plastomes exactly with independently checked atpB/rbcL endpoint annotations.',
                           'All other 153 candidates lack a qualifying 90%-identity/90%-query-coverage public or expected-target panel HSP in the bounded frozen searches; absence of a hit is not a false-sequence classification.',
                           'The three mapped plastomes each have exactly one effective primer pair, at 1359/1379/1359 bp, outside configured 400–1000 bp. Reverse primer is inside annotated atpB; atpE lies within the interval.',
                           'Zero strong expected-target panel HSPs for the positive atpB-rbcL sequence coexist with exact full annotated genomic support; these are distinct evidence routes.'],
            'limitations':['Pinned public sources and panel regions are bounded. Reference matches establish sequence relationship, not sample whole-product linkage or biological truth.',
                          'Effective-primer coordinate mapping is a reference assay-configuration concern; it does not demonstrate a production assembly invariant failure or justify deleting any emitted hypothesis.'],
            'artifacts':{str(p.relative_to(BASE)):{'sha256':sha(p.read_bytes()),'bytes':p.stat().st_size} for p in [report_path,mapping_path,BASE/'reference/assessment-v2/reference_overlay.json',NEW/'execution.json']},
            'reviewer_source_sha256':sha(Path(__file__).read_bytes())}
    with (HERE/'reference_review.json').open('x') as stream:stream.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'pass','candidate_scope':154,'query_reuse':dict(used),'mapping_checks':mapping_checks,'annotated_exact_target_support_candidates':len(positive)}))

if __name__=='__main__':main()

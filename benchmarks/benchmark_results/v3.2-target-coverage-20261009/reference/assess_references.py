#!/usr/bin/env python3
"""Join frozen old/new reference searches for all six new targets."""
import csv
import hashlib
import importlib.util
import json
from collections import defaultdict
from pathlib import Path
import sys
import yaml
import map_assay

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OLD=ROOT/'benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/blast-gains-v3'
NEW=HERE/'blast-missing-v1'
sys.path.insert(0,str(ROOT))
from scripts.sharkmer_validate import reference_provenance
from scripts.sharkmer_validate.reference_targets import logical_gene_name
from scripts.sharkmer_validate.primer_analysis import IUPAC_SETS
from scripts.sharkmer_validate.runner import derive_gene_name


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(path,sha):
    if digest(path)!=sha:raise ValueError(f'Source hash mismatch {path}')


def hit_rows(path):
    with path.open() as f:
        for row in csv.reader(f,delimiter='\t'):
            if not row:continue
            q,s,pid,length,qlen,qs,qe,slen,ss,se,e,b,gaps,qcovs=row
            yield {'query_id':q,'subject_id':s,'identity_pct':float(pid),'alignment_length':int(length),
                   'query_length':int(qlen),'qstart':int(qs),'qend':int(qe),'subject_length':int(slen),
                   'sstart':int(ss),'send':int(se),'evalue':float(e),'bitscore':float(b),'gaps':int(gaps),
                   'hsp_query_coverage_pct':round(100*(abs(int(qe)-int(qs))+1)/int(qlen),3)}


def best(hits):return max(hits,key=lambda h:(h['hsp_query_coverage_pct'],h['identity_pct'],h['bitscore']),default=None)


def main():
    out=HERE/'assessment-v2';out.mkdir(exist_ok=True)
    if any(out.iterdir()):raise FileExistsError('Refusing to overwrite assessment outputs')
    receipt=json.load(open(NEW/'execution.json'));assert receipt['status']=='pass' and receipt['sources_unchanged']
    for path,expected in receipt['source_hashes_after'].items():verify(Path(path),expected)
    for path,expected in receipt['outputs'].items():verify(NEW/path,expected['sha256'])
    summary=json.load(open(OLD/'summary.json'));inventory_path=HERE.parent/'matrix/novel_target_inventory.json'
    inventory=json.load(open(inventory_path));source_meta=json.load(open(OLD/'public_sources_metadata.json'))
    catalog,catalog_receipt=reference_provenance.load_catalog()
    assert catalog_receipt['sha256']==summary['catalog']['sha256']
    for item in [summary['source_hits'],summary['query_fasta']]:verify(ROOT/item['path'],item['sha256'])
    oldrows={x['query_id']:x for x in json.load(open(OLD/'per_sequence.json'))}
    newrows={(x['accession'],x['gene'],x['sha256']):x for x in json.load(open(NEW/'queries.json'))}
    panel_meta={};panels={}
    for panel_name in ['angiospermae','insecta']:
        panel=yaml.safe_load((ROOT/'panels'/f'{panel_name}.yaml').read_text());panels[panel_name]=panel
        refs=reference_provenance.verified_references(panel)
        fasta=''.join(f'>reference_{i:03d}\n'+ref['sequence']+'\n' for i,ref in enumerate(refs))
        assert hashlib.sha256(fasta.encode()).hexdigest()==summary['panel'][panel_name]['references']['sha256']
        verify(ROOT/summary['panel'][panel_name]['hits']['path'],summary['panel'][panel_name]['hits']['sha256'])
        panel_meta[panel_name]={f'reference_{i:03d}':{'gene':r['gene_name'],'logical_gene':logical_gene_name(r['gene_name'],panel),
                                                   'taxon':r['taxon'],'accession':r['accession']} for i,r in enumerate(refs)}
    hits=defaultdict(list)
    for namespace,folder in [('saved',OLD),('new',NEW)]:
        for scope in ['public','angiospermae','insecta']:
            name=('public_sources_hits.tsv' if scope=='public' else scope+'_reference_hits.tsv') if namespace=='saved' else scope+'_hits.tsv'
            for hit in hit_rows(folder/name):
                meta=source_meta[hit['subject_id']] if scope=='public' else panel_meta[scope][hit['subject_id']]
                hits[(namespace,scope,hit['query_id'])].append({**hit,**meta})
    assay=json.load(open(HERE/'assay-v1/mapping.json'))
    feature_context={}
    for source in assay['records']:
        path=HERE/'assay-v1'/(source['accession']+'.gb')
        verify(path,source['raw_genbank_sha256'])
        assert catalog['records'][source['accession']]['sequence_sha256']==source['source_sha256']
        feature_context[source['accession']]=map_assay.genes_from_gb(path.read_text())
    mapped_spans={r['accession']:(r['effective_primer_pairs'][0]['forward']['source_span'][0],r['effective_primer_pairs'][0]['reverse_binding']['source_span'][1]) for r in assay['records']}
    rows=[]
    for c in inventory:
        key=(c['accession'],c['gene'],c['sha256']);qid=c['existing_blast_query_id'];namespace='saved' if qid else 'new'
        if qid:
            prior=oldrows[qid];assert (prior['accession'],prior['gene'],prior['sha256'])==key
        else:qid=newrows[key]['query_id']
        public=hits[(namespace,'public',qid)];panel=hits[(namespace,c['panel'],qid)]
        if namespace=='saved':assert len(public)==prior['public_hit_count'] and len(panel)==prior['panel_hit_count']
        assert all(h['query_length']==c['length'] for h in public+panel)
        logical=logical_gene_name(c['target'],panels[c['panel']]);expected=[h for h in panel if h['logical_gene']==logical]
        same_taxon=[h for h in public if h['taxon']==c['taxon']]
        strong=lambda h:h['identity_pct']>=90 and h['hsp_query_coverage_pct']>=90
        exact=[accession for accession,r in catalog['records'].items() if r['organism']==c['taxon'] and
               (c['sequence'] in r['sequence'] or reference_provenance.reverse_complement(c['sequence']) in r['sequence'])]
        exact_placements=[]
        for accession in exact:
            sequence=catalog['records'][accession]['sequence'];features=feature_context.get(accession,[])
            for strand,query in [('+',c['sequence']),('-',reference_provenance.reverse_complement(c['sequence']))]:
                offset=0
                while True:
                    start=sequence.find(query,offset)
                    if start<0:break
                    end=start+len(query);offset=start+1
                    left=[g for g in features if g['start']<=start<g['end']]
                    right=[g for g in features if g['start']<end<=g['end']]
                    supports=(c['target']=='atpB-rbcL' and any(g['gene']=='atpB' for g in left) and any(g['gene']=='rbcL' for g in right))
                    exact_placements.append({'accession':accession,'source_start':start,'source_end':end,'strand':strand,
                                             'left_endpoint_genes':left,'right_endpoint_genes':right,
                                             'annotated_target_boundaries_agree':supports})
        annotated_exact=any(p['annotated_target_boundaries_agree'] for p in exact_placements)
        primer=next(p for p in panels[c['panel']]['primers'] if derive_gene_name(p)==c['target'])
        trim=min(primer.get('trim',15),18);fwd=primer['forward_seq'][-trim:];rev=reference_provenance.reverse_complement(primer['reverse_seq'][-trim:])
        primer_mismatches={'forward':sum(b not in IUPAC_SETS[a] for a,b in zip(fwd,c['sequence'][:len(fwd)])),
                           'reverse_binding':sum(b not in IUPAC_SETS[a] for a,b in zip(rev,c['sequence'][-len(rev):])),
                           'trim':trim,'allowed_each':primer['mismatches']}
        locus=[]
        if c['target']=='trnV-atpE':
            for h in public:
                if h['accession'] in mapped_spans:
                    a,b=mapped_spans[h['accession']];hs,he=min(h['sstart'],h['send'])-1,max(h['sstart'],h['send'])
                    if hs<b and he>a:locus.append({**h,'overlap_with_mapped_assay_bases':min(b,he)-max(a,hs)})
        status=('exact_full_public_sequence_and_annotated_target_support' if annotated_exact else
                'expected_target_reference_support' if any(map(strong,expected)) else
                'partial_expected_target_hit_full_sequence_unresolved' if expected else
                'non_target_or_ambiguous_reference_hits' if public or panel else 'no_significant_hit_in_pinned_catalog')
        rows.append({'accession':c['accession'],'gene':c['gene'],'sha256':c['sha256'],'length':c['length'],
                     'taxon':c['taxon'],'target':c['target'],'query_namespace':namespace,'query_id':qid,
                     'public_hsp_count':len(public),'panel_hsp_count':len(panel),'expected_target_hsp_count':len(expected),
                     'same_taxon_hsp_count':len(same_taxon),'best_public_hsp':best(public),'best_expected_target_hsp':best(expected),
                     'best_same_taxon_hsp':best(same_taxon),'best_panel_hsp':best(panel),
                     'strong_public_single_hsp_count':sum(map(strong,public)),
                     'strong_expected_target_single_hsp_count':sum(map(strong,expected)),
                     'exact_full_sequence_same_taxon_sources':exact,
                     'exact_full_sequence_placements':exact_placements,'annotated_exact_target_support':annotated_exact,
                     'primer_endpoint_mismatches':primer_mismatches,
                     'trnV_mapped_assay_overlap_hsps':locus,
                     'reference_evidence_status':status,'target_identity': ('supported_by_exact_public_sequence_and_independent_feature_boundaries' if annotated_exact else
                               'supported_by_qualifying_reference_alignment' if any(map(strong,expected)) else 'not_established'),
                     'full_sequence_reference_support':('exact_full_length_same_taxon_public_sequence' if exact else
                                                        'qualifying_alignment' if any(map(strong,expected)) else 'not_established'),
                     'whole_product_biological_validation':'not_established'})
    groups={}
    for row in rows:
        group=groups.setdefault(row['accession']+':'+row['gene'],{'count':0,'public_hit_any':0,'expected_target_hit_any':0,'strong_expected_target':0,'strong_any_public':0,'exact_full_same_taxon':0,'annotated_exact_target_support':0,'same_taxon_hit_any':0,'mapped_trnV_assay_overlap_any':0})
        group['count']+=1
        for field,condition in [('public_hit_any',row['public_hsp_count']),('expected_target_hit_any',row['expected_target_hsp_count']),
                                ('strong_expected_target',row['strong_expected_target_single_hsp_count']),('strong_any_public',row['strong_public_single_hsp_count']),
                                ('exact_full_same_taxon',row['exact_full_sequence_same_taxon_sources']),('same_taxon_hit_any',row['same_taxon_hsp_count']),
                                ('annotated_exact_target_support',row['annotated_exact_target_support']),
                                ('mapped_trnV_assay_overlap_any',row['trnV_mapped_assay_overlap_hsps'])]:group[field]+=bool(condition)
    report={'status':'completed','candidates':rows,'by_target':groups,'total':len(rows),'reused_saved_queries':113,'new_queries':41,
            'single_hsp_reference_criterion':'querycoverage>=90pct AND identity>=90pct; this is reference similarity, not sample biological truth',
            'limitations':['Pinned121records and61panelregions are bounded and not an exhaustive public search.',
                           'No significant hit does not prove a wrong sequence; partial or wrong-target hits do not validate a named target.',
                           'Related-taxonomy and short public references do not establish full sample haplotypes or repeat phase.',
                           'trnV-atpE has no accepted panel reference; independent effective-assay coordinates are evaluated separately.'],
            'source_receipts':{str(p):digest(p) for p in [inventory_path,Path(__file__),HERE/'map_assay.py',OLD/'summary.json',OLD/'review.json',OLD/'per_sequence.json',NEW/'execution.json',HERE/'assay-v1/mapping.json']}}
    with (out/'report.json').open('x') as f:f.write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    overlay={'sequence_evidence':[{'accession':r['accession'],'gene':r['gene'],'sha256':r['sha256'],
                                 'reference_evidence_status':r['reference_evidence_status'],
                                 'reference_evidence_links':[str(out/'report.json')],
                                 'reference_evidence_details':{k:r[k] for k in ['query_namespace','query_id','target_identity','full_sequence_reference_support','expected_target_hsp_count','same_taxon_hsp_count','trnV_mapped_assay_overlap_hsps','exact_full_sequence_placements','annotated_exact_target_support','primer_endpoint_mismatches']}} for r in rows]}
    with (out/'reference_overlay.json').open('x') as f:f.write(json.dumps(overlay,indent=2,sort_keys=True)+'\n')
    fields=['accession','gene','sha256','length','query_namespace','query_id','public_hsp_count','panel_hsp_count','expected_target_hsp_count','same_taxon_hsp_count','strong_expected_target_single_hsp_count','reference_evidence_status','target_identity','full_sequence_reference_support']
    with (out/'per_candidate.tsv').open('x') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows({k:r[k] for k in fields} for r in rows)
    print(json.dumps(groups,indent=2))


if __name__=='__main__':main()

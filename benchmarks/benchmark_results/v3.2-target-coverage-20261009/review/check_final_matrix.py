#!/usr/bin/env python3
"""Independently audit scoped evidence joins without new biological scans."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tarfile

HERE=Path(__file__).resolve().parent
BASE=HERE.parent

def load(path):return json.loads(path.read_text())
def sha(data):return hashlib.sha256(data).hexdigest()
def key(row):return tuple(row[k] for k in ('accession','gene','sha256'))

def main():
    matrix=BASE/'matrix';final=matrix/'final-v1';validation=load(final/'validation.json')
    for source,expected in validation['source_files'].items():
        data=Path(source).read_bytes()
        assert sha(data)==expected['sha256'] and len(data)==expected['size_bytes']
    for name,expected in validation['outputs'].items():
        data=(final/name).read_bytes()
        assert sha(data)==expected['sha256'] and len(data)==expected['size_bytes']
    archives={}
    try:
        for source,expected in validation['archive_members'].items():
            archive,member=source.split('::',1)
            if archive not in archives:archives[archive]=tarfile.open(archive,'r:gz')
            data=archives[archive].extractfile(member).read()
            assert sha(data)==expected['sha256'] and len(data)==expected['size_bytes']
    finally:
        for archive in archives.values():archive.close()
    base_rows={(r['accession'],r['gene']):r for r in load(matrix/'target_matrix.json')}
    rows={(r['accession'],r['gene']):r for r in load(final/'target_matrix.json')}
    base_current={key(r):r for r in load(matrix/'sequence_inventory.json')}
    current={key(r):r for r in load(final/'sequence_inventory.json')}
    assert set(base_rows)==set(rows) and len(rows)==118
    assert set(base_current)==set(current) and len(current)==326
    old_join={key(r):r for r in load(matrix/'available-v1/sequence_inventory.json')}
    reviewed={}
    for name in ('matrix_review.json','read_assay_review.json','read_export_review.json','reference_review.json','annotation_review.json'):
        path=HERE/name
        assert load(path)['status']=='pass'
        reviewed[name]=sha(path.read_bytes())
    overlay_paths=[BASE/'reads/results-v1/evidence_overlay.json',matrix/'control_evidence_overlay.json',
                   BASE/'reference/assessment-v2/reference_overlay.json',BASE/'reference/annotation_overlay.json']
    overlays=[(str(p.resolve()),load(p)) for p in overlay_paths]
    controls={key(r):r for s in load(BASE/'reads/results-v1/summary.json')['samples'] for r in s['positive_controls']}
    for r in load(matrix/'control_evidence_overlay.json')['sequence_evidence']:
        assert r['read_evidence_details']==controls[key(r)]
    assert len(controls)==8
    assessed=Counter()
    for k,row in current.items():
        for field,value in base_current[k].items():
            if field not in ('read_evidence_status','read_evidence_links','reference_evidence_status','reference_evidence_links'):
                assert row[field]==value
        assert row['whole_product_validation']=='not_established'
        for kind in ('read','reference'):
            details=row[kind+'_evidence_details']
            expected=old_join[k][kind+'_evidence_details'].copy()
            for source,overlay in overlays:
                for item in overlay.get('sequence_evidence',[]):
                    if key(item)==k and kind+'_evidence_status' in item:
                        expected.append({'status':item[kind+'_evidence_status'],'links':item[kind+'_evidence_links'],
                                         'details':item[kind+'_evidence_details'],'evidence_source':source,'whole_product_truth':'not_established'})
            def normalized(entries):
                return [{**e,'evidence_source':str(Path(e['evidence_source']).resolve()) if Path(e['evidence_source']).exists() else e['evidence_source']} for e in entries]
            assert normalized(details)==normalized(expected)
            assert row[kind+'_evidence_status']==(';'.join(sorted({e['status'] for e in details})) or 'unassessed')
            assert row[kind+'_evidence_links']==sorted({p for e in details for p in e['links']})
            assessed[kind]+=bool(details)
    outcomes={}
    for run in load(matrix/'validation.json')['current_runs']:
        for r in load(Path(run['receipt']))['outputs']['pcr_results']:
            outcomes[run['accession'],r['gene_name']]=r
    for k,row in rows.items():
        for field,value in base_rows[k].items():
            if field not in ('read_evidence_status','read_evidence_links','reference_evidence_status','reference_evidence_links'):
                assert row[field]==value
        selected=[r for identity,r in current.items() if identity[:2]==k]
        result=outcomes[k]
        assert row['target_assembly_status']==result['status'] and row['target_failure_reason']==result.get('failure_reason')
        flags=sorted({f for step in result.get('threshold_diagnostics',[]) for f in ('dfs_limit_reached','node_budget_reached','path_quota_reached','max_length_reached') if step.get(f)})
        assert row['search_limit_flags']==flags
        uncertain=result.get('uncertain_products',{})
        assert row['search_limits_reached']==(bool(flags) or uncertain.get('search_limits_reached',False))
        assert row['retention_truncated']==uncertain.get('retention_truncated',False)
        assert row['marker_detail_truncated']==uncertain.get('marker_detail_truncated',False)
        for kind in ('read','reference'):
            counted=[r for r in selected if r[kind+'_evidence_details']]
            assert row[kind+'_assessed_candidates']==len(counted)
            assert row[kind+'_unassessed_candidates']==len(selected)-len(counted)
            assert row[kind+'_status_counts']==dict(Counter(r[kind+'_evidence_status'] for r in counted))
            expected_status='not_applicable_no_emitted_product' if not selected else 'unassessed' if not counted else 'all_emitted_candidates_have_scoped_assessments' if len(counted)==len(selected) else 'subset_of_emitted_candidates_has_scoped_assessments'
            assert row[kind+'_evidence_status']==expected_status
            expected_notes=[item for _,overlay in overlays for item in overlay.get('target_evidence',[]) if (item['accession'],item['gene'])==k and kind+'_evidence_status' in item]
            assert row[kind+'_target_notes']==expected_notes
        supported=sum(any(e['details'].get('target_support')=='supported' for e in r['reference_evidence_details']) for r in selected)
        assert row['historical_gene_supported_products']==supported
    assert assessed=={'read':228,'reference':313}
    assert sum(r['historical_gene_supported_products'] for r in rows.values())==51
    assert Counter(r['target_assembly_status'] for r in rows.values())=={'success':64,'uncertain':10,'fail':44}
    assert {key(r) for r in load(final/'novel_target_inventory.json')}=={k for k,r in current.items() if r['newly_represented_target']}
    compact=list(csv.DictReader((final/'target_matrix_compact.csv').open()))
    assert len(compact)==118
    compact_keys={(r['accession'],r['target']) for r in compact}
    assert compact_keys=={(r['accession'],r['target']) for r in rows.values()}
    for exported in compact:
        row=next(r for r in rows.values() if (r['accession'],r['target'])==(exported['accession'],exported['target']))
        for field,text in exported.items():
            value=row[field]
            assert (json.loads(text)==value if isinstance(value,dict) else text==('' if value is None else str(value))),field
    result={'status':'pass','target_rows':118,'current_identity_count':326,'read_assessed_candidates':228,'read_unassessed_candidates':98,
            'reference_assessed_candidates':313,'reference_unassessed_candidates':13,'historical_gene_supported_identities_preserved':51,
            'source_receipts_verified':len(validation['source_files']),'archive_member_receipts_verified':len(validation['archive_members']),
            'findings':'Frozen counts and all exact identities unchanged. Every available old evidence entry and every new full-identity overlay is preserved; each target count/status/limit flag and scoped assessment count agrees with source.',
            'caveats':['Assessed is not supported. Reference and local read evidence remain separate; target summaries do not validate untested products.',
                       'Only one of the154 new-target candidates has exact full annotated public target support plus local Q20 tiles; this does not establish a single full sample molecule.',
                       'All44 zero-output targets are actual assembly failures with explicit reasons, not biological negative detections.',
                       'All326 whole_product_validation fields remain not_established. trnV-atpE annotation/length concern is a target note with unchanged emitted identities.'],
            'reviewed_component_receipts':reviewed,'validation_sha256':sha((final/'validation.json').read_bytes()),
            'reviewer_source_sha256':sha(Path(__file__).read_bytes())}
    with (HERE/'final_matrix_review.json').open('x') as stream:stream.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','target_rows','current_identity_count','read_assessed_candidates','reference_assessed_candidates','source_receipts_verified')}))

if __name__=='__main__':main()

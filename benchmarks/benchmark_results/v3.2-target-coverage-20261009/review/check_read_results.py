#!/usr/bin/env python3
"""Independent saved-result and selected-ledger audit; no FASTQ rescan."""
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
COMP=str.maketrans('ACGT','TGCA')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def load(path):
    return json.loads(path.read_text())

def fingerprint(path):
    data=path.read_bytes()
    return {'sha256':sha(data),'bytes':len(data)}

def rc(sequence):
    return sequence.translate(COMP)[::-1]

def oriented(sequence):
    reverse=rc(sequence)
    return [('both',sequence)] if reverse==sequence else [('+',sequence),('-',reverse)]

def positions(sequence,pattern):
    start=sequence.find(pattern)
    while start!=-1:
        yield start
        start=sequence.find(pattern,start+1)

def main():
    out=BASE/'reads/results-v1'
    summary=load(out/'summary.json')
    candidate_rows=load(out/'per_candidate.json')
    candidate_lookup={(r['accession'],r['gene'],r['sha256']):r for r in candidate_rows}
    overlay=load(out/'evidence_overlay.json')
    assert len(candidate_rows)==len(candidate_lookup)==len(overlay['sequence_evidence'])==154
    for row in overlay['sequence_evidence']:
        key=tuple(row[k] for k in ('accession','gene','sha256'))
        assert row['read_evidence_details']==candidate_lookup[key]
        assert row['read_evidence_status']==candidate_lookup[key]['local_read_status']
    artifacts={}
    for name,expected in load(out/'checksums.json').items():
        assert fingerprint(out/name)==expected
        artifacts[str((out/name).relative_to(BASE))]=expected
    for name,expected in summary['source_files'].items():
        assert fingerprint(Path(name))==expected
    verified_rows=[]
    sample_results=[]
    for manifest_path in sorted((BASE/'reads/manifests-v2').glob('*.json')):
        manifest=load(manifest_path)
        accession=manifest['accession']
        run=BASE/'reads/run-v1'/accession
        receipt=load(run/'execution.json')
        raw=load(run/'raw.json')
        metrics=load(run/'time.json')
        assert receipt['status']=='pass' and raw['status']=='completed'
        assert receipt['returncode']==metrics['exit_code']==0
        assert receipt['source_files_before']==receipt['source_files_after'] and receipt['sources_unchanged']
        assert len(receipt['cpu_affinity'])==1
        assert metrics['wall_seconds']<600 and metrics['max_rss_kib']<2*1024*1024
        for name,expected in receipt['source_files_before'].items():
            assert fingerprint(Path(name))==expected
        for name,expected in receipt['outputs'].items():
            assert fingerprint(run/name)==expected
            artifacts[str((run/name).relative_to(BASE))]=expected
        artifacts[str((run/'execution.json').relative_to(BASE))]=fingerprint(run/'execution.json')
        assert raw['verified_prefix']==manifest['expected_prefix']
        histogram={int(k):v for k,v in raw['read_length_histogram'].items()}
        assert sum(histogram.values())==raw['verified_prefix']['records']==1000000
        assert sum(k*v for k,v in histogram.items())==raw['verified_prefix']['bases']
        events=manifest['events']
        counts=raw['events']
        assert [e['id'] for e in events]==[e['id'] for e in counts]
        candidates={c['id']:c for c in manifest['candidates']}
        index_by_id={e['id']:i for i,e in enumerate(events)}
        for event,count in zip(events,counts):
            pattern=event['sequence']
            expected=[]
            for candidate in manifest['candidates']:
                for strand,sequence in oriented(candidate['sequence']):
                    for start in positions(sequence,pattern):
                        expected.append({'candidate_id':candidate['id'],'strand':strand,'start':start,'end':start+len(pattern)})
            assert count['candidate_placements']==expected
            assert count['eligible_length_records']==sum(v for k,v in histogram.items() if k>=len(pattern))
            assert count['callability']==('assayed' if count['eligible_length_records'] else 'longer_than_every_read')
            assert 0<=count['q20_records']<=count['raw_records']<=count['eligible_length_records']
            assert count['raw_only_records']==count['raw_records']-count['q20_records']
            assert count['raw_occurrences']>=count['raw_records'] and count['q20_occurrences']>=count['q20_records']
            assert count['q20_occurrences']<=count['raw_occurrences']
        retained_q20=Counter()
        retained_raw_only=Counter()
        ledger=raw['read_evidence']
        assert len({r['ordinal'] for r in ledger})==len(ledger)
        assert [r['ordinal'] for r in ledger]==sorted(r['ordinal'] for r in ledger)
        for read in ledger:
            assert 1<=read['ordinal']<=1000000 and len(read['sequence'])==len(read['quality'])
            expected=[]
            for i,event in enumerate(events):
                matches=[]
                for strand,word in oriented(event['sequence']):
                    for start in positions(read['sequence'],word):
                        end=start+len(word)
                        qmin=min(ord(q)-33 for q in read['quality'][start:end])
                        matches.append({'strand':strand,'start':start,'end':end,'minimum_phred':qmin,'q20':qmin>=20})
                if matches:
                    expected.append({'event_index':i,'matches':matches})
                    if any(m['q20'] for m in matches):
                        retained_q20[i]+=1
                    else:
                        retained_raw_only[i]+=1
            assert read['observations']==expected
        ledger_summary=raw['ledger']
        assert ledger_summary['stored_reads']==len(ledger)<=ledger_summary['limits']['reads']
        assert ledger_summary['stored_bases']==sum(len(r['sequence']) for r in ledger)<=ledger_summary['limits']['bases']
        assert ledger_summary['stored_match_positions']==sum(len(o['matches']) for r in ledger for o in r['observations'])<=ledger_summary['limits']['matches']
        assert not ledger_summary['cap_omitted_requests_by_reason']
        for i,count in enumerate(counts):
            assert count['retained_q20_records']==retained_q20[i]<=count['q20_records']
            assert count['retained_raw_only_records']==retained_raw_only[i]<=count['raw_only_records']
            assert count['example_ledger_complete_for_event']==(retained_q20[i]+retained_raw_only[i]==count['raw_records'])
            assert retained_q20[i]>=min(count['q20_records'],8)
            assert retained_raw_only[i]>=min(count['raw_only_records'],2)
            assert count['ledger_cap_omitted_requests']==0
        sample=next(s for s in summary['samples'] if s['accession']==accession)
        assert sample['ledger']==ledger_summary and sample['verified_prefix']==raw['verified_prefix']
        assert sample['read_length_histogram']==raw['read_length_histogram'] and sample['metrics']==metrics
        control_ids=set(manifest['positive_control_ids'])
        for identifier in manifest['focused_candidate_ids']+manifest['positive_control_ids']:
            candidate=candidates[identifier]
            rows=[(i,m) for i,e in enumerate(events) for m in e['memberships'] if m['candidate_id']==identifier]
            tiles=[(i,m) for i,m in rows if m['role'].endswith('_tile')]
            tile_indices={i for i,_ in tiles}
            endpoints={i for i,m in rows if m['role']=='nearest_difference_endpoint'}
            observed={i for i in tile_indices|endpoints if counts[i]['q20_records']}
            qcovered={j for i,m in tiles if counts[i]['q20_records'] for j in range(m['start'],m['end'])}
            rcovered={j for i,m in tiles if counts[i]['raw_records'] for j in range(m['start'],m['end'])}
            key=tuple(candidate[k] for k in ('gene','sha256'))
            row=next(c for c in sample['positive_controls'] if (c['gene'],c['sha256'])==key) if identifier in control_ids else candidate_lookup[(accession,*key)]
            expected={'scheduled_tiles':len(tiles),'distinct_tile_strings':len(tile_indices),
                      'q20_observed_tile_strings':sum(bool(counts[i]['q20_records']) for i in tile_indices),
                      'raw_observed_tile_strings':sum(bool(counts[i]['raw_records']) for i in tile_indices),
                      'q20_observed_tile_span_union_bases':len(qcovered),'raw_observed_tile_span_union_bases':len(rcovered),
                      'distinct_nearest_endpoint_strings':len(endpoints),
                      'q20_observed_nearest_endpoint_strings':sum(bool(counts[i]['q20_records']) for i in endpoints),
                      'raw_observed_nearest_endpoint_strings':sum(bool(counts[i]['raw_records']) for i in endpoints),
                      'overlength_nearest_endpoint_strings':sum(counts[i]['eligible_length_records']==0 for i in endpoints),
                      'incomplete_flank_endpoint_memberships':sum(x['candidate_id']==identifier for x in manifest['uncallable_endpoints']),
                      'q20_observed_strings_unique_to_candidate_in_finite_context':sum({p['candidate_id'] for p in counts[i]['candidate_placements']}=={identifier} for i in observed),
                      'local_read_status':'local_exact_q20_occurrences_observed' if observed else 'raw_only_local_exact_occurrences' if any(counts[i]['raw_records'] for i in tile_indices|endpoints) else 'no_tested_local_exact_occurrence',
                      'whole_product_validation':'not_established'}
            assert all(row[k]==v for k,v in expected.items()),identifier
            if identifier in control_ids:
                assert row['positive_control_pass'] and expected['q20_observed_tile_strings']>0
            else:
                verified_rows.append(row)
        sample_results.append({'accession':accession,'events_verified':len(events),'all_candidate_placements_verified':sum(len(c['candidate_placements']) for c in counts),
                               'selected_ledger_reads_independently_checked':len(ledger),'saved_match_positions_independently_checked':ledger_summary['stored_match_positions'],
                               'positive_controls_pass':len(control_ids),'global_ledger_caps_reached':False,
                               'callable_events':sum(c['eligible_length_records']>0 for c in counts),'metrics':metrics})
    assert len(verified_rows)==154
    assert sum(r['local_read_status']=='local_exact_q20_occurrences_observed' for r in verified_rows)==1
    assert sum(r['local_read_status']=='no_tested_local_exact_occurrence' for r in verified_rows)==153
    result={'status':'pass','samples':sample_results,'focused_candidates':154,'positive_controls_pass':8,
            'findings':'One Liriodendron atpB-rbcL candidate has all scheduled 61 bp tiles observed at Q20; other 153 candidates have no scheduled tile or callable nearest-contrast exact string in raw or Q20 complete counts.',
            'complete_count_evidence':'Complete counts are accepted from hash-bound successful prefix scans and independently controlled scanner semantics. This review does not rescan FASTQ; selected ledgers cannot establish negative findings.',
            'limitations':['Saved examples are selected and capped; per-event completeness is verified separately. No global cap was reached.',
                          'All 833 emitted events have at least some eligible reads, but 22 missing-complete-flank endpoint memberships remain explicitly uncallable.',
                          'Observed tiles may be shared or overlapping; covered-base unions are marginal coverage and do not establish a single full-length molecule.',
                          'No exact local observation is not biological absence; reference relationship is assessed separately.'],
            'artifacts':artifacts,'reviewer_source_sha256':sha(Path(__file__).read_bytes()),
            'protocol_review_sha256':sha((HERE/'protocol_review_v2.json').read_bytes()),
            'independent_synthetic_controls_sha256':sha((HERE/'exact_independent_controls_v2.json').read_bytes())}
    with (HERE/'read_assay_review.json').open('x') as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','focused_candidates','positive_controls_pass','samples')}))

if __name__=='__main__':
    main()

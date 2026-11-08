import hashlib
import json
from pathlib import Path
import shlex
from collections import Counter

repo = Path.cwd()
out = repo / 'benchmarks/benchmark_results/v3.2-final-check-20261108/inputs/input_inventory.json'
out.parent.mkdir(parents=True, exist_ok=True)
source_bindings = {}
def fp(path):
    p = Path(path).resolve()
    if p.stat().st_size > 50_000_000:
        raise ValueError(f'Refusing large file hash in lightweight oracle preparation: {p}')
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'size_bytes': p.stat().st_size}
def read(path):
    p = Path(path).resolve()
    source_bindings[str(p.relative_to(repo))] = fp(p)
    return json.loads(p.read_text())
def identities(records):
    return Counter((r['gene'],r['length'],r['sha256'],r['assembly_status']) for r in records)
def outputs(d):
    result = {'counts': d['counts'], 'records': [], 'targets': []}
    for r in d['records']:
        assert len(r['sequence']) == r['length']
        assert hashlib.sha256(r['sequence'].encode()).hexdigest() == r['sha256']
        record = {k:r[k] for k in ('gene','length','sha256','assembly_status')}
        record['id'] = r['header'].split(maxsplit=1)[0]
        if r['assembly_status'] == 'standard':
            fields = dict(t.split('=',1) for t in shlex.split(r['header'])[1:] if '=' in t)
            fields.pop('assembly_status',None)
            record['ordinary_fields'] = fields
        result['records'].append(record)
    for t in d['pcr_results']:
        result['targets'].append({
            'gene_name': t['gene_name'], 'status': t['status'], 'n_products': t['n_products'],
            'n_standard_products': t.get('n_standard_products',t['n_products']),
            'n_uncertain_products': t.get('n_uncertain_products',0),
            'product_lengths': t.get('product_lengths',[])
        })
    return result

availability_path = repo/'tmp/inclusive-release-20260929/availability.json'
availability = read(availability_path)
historical_inputs = read(repo/'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/inventory.json')
historical_panels = {row['name']: row['assembly_panel'] for row in historical_inputs['panels']}
for row in historical_panels.values():
    assert fp(row['path']) == row
plant_checksums = read(repo/'benchmarks/benchmark_results/v3.2-plant-assay-correction-20261010/proposal-v2/checksums.json')
candidate_inventory = read(repo/'benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/final-v1/sequence_inventory.json')
release_inventory = read(repo/'benchmarks/benchmark_results/v3.2-target-coverage-20261009/matrix/historical_inventory.json')
current_inventory = read(repo/'benchmarks/benchmark_results/v3.2-plant-assay-correction-20261010/assessment-v1/derived_combined_inventory.json')
plant_summary = read(repo/'benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/results-v1/plant/summary.json')
assert len(candidate_inventory)==326 and len(release_inventory)==86 and len(current_inventory)==287
samples = {}
for i,a in enumerate(availability['samples']):
    accession = a['accession']
    sample_prefix = f"{a['panel_id']}_{accession}_1000000"
    prefix = f'{i:03}_{sample_prefix}_pair1_'
    baseline_path=repo/'tmp/inclusive-release-20260929/comparison-v1'/f'{prefix}baseline.json'
    candidate_path=repo/'tmp/sequence-retention-20261001/full-v3'/f'{prefix}candidate.json'
    baseline=read(baseline_path); candidate=read(candidate_path)
    assert baseline['status']=='pass' and candidate['status']=='pass'
    assert baseline['prefix_before']==candidate['prefix_before']==a['expected_prefix']
    expected_input={'path':a['available_path'], 'sha256':a['expected_prefix']['sha256'], 'size_bytes':a['expected_prefix']['size_bytes']}
    assert baseline['input_file']==candidate['input_file']==expected_input
    assert Path(expected_input['path']).is_file() and Path(expected_input['path']).stat().st_size==expected_input['size_bytes']
    assert baseline['invocation']['sample_prefix']==candidate['invocation']['sample_prefix']==sample_prefix
    command=candidate['invocation']['command']
    panel_path=command[command.index('--pcr-panel-file')+1]
    assert panel_path==baseline['invocation']['command'][baseline['invocation']['command'].index('--pcr-panel-file')+1]
    assert fp(panel_path) == historical_panels[a['panel_id']]
    expected={'baseline':outputs(baseline['outputs']),'candidate':outputs(candidate['outputs'])}
    assert expected['baseline']['counts']==expected['candidate']['counts']
    ci=[r for r in candidate_inventory if r['accession']==accession]
    ri=[dict(r,assembly_status='standard') for r in release_inventory if r['accession']==accession]
    assert identities(expected['candidate']['records'])==identities(ci)
    assert identities(expected['baseline']['records'])==identities(ri)
    samples[accession]={'accession':accession,'panel_id':a['panel_id'],'taxon':candidate['invocation']['taxon'],
        'assay_scope':'historical_unchanged_panel','sample_prefix':sample_prefix,'max_reads':1000000,
        'input':expected_input,'expected_prefix':a['expected_prefix'],'panel':fp(panel_path),
        'expected':expected,'source_receipts':{'baseline':fp(baseline_path),'candidate':fp(candidate_path)}}
for job in plant_summary['jobs']:
    if job['scope']!='real' or not job['sample_prefix'].endswith('_after'):
        continue
    accession=job['sample_id']; key=accession+'_corrected'
    command=job['command']; panel_path=command[command.index('--pcr-panel-file')+1]
    assert fp(panel_path)['sha256']==fp(repo/'panels/angiospermae.yaml')['sha256']
    assert fp(panel_path)['sha256']==plant_checksums['angiospermae-after.yaml']
    assert job['expected_prefix']==samples[accession]['expected_prefix']
    assert job['input_file']==samples[accession]['input']
    expected={'candidate':outputs(job['outputs'])}
    ci=[r for r in current_inventory if r['accession']==accession]
    assert identities(expected['candidate']['records'])==identities(ci)
    samples[key]={'accession':accession,'panel_id':'angiospermae','taxon':samples[accession]['taxon'],
        'assay_scope':'corrected_plant_panel_2.0.0','sample_prefix':job['sample_prefix'],'max_reads':1000000,
        'input':job['input_file'],'expected_prefix':job['expected_prefix'],'panel':fp(panel_path),
        'expected':expected,'source_receipts':{'candidate':fp(repo/'benchmarks/benchmark_results/v3.2-prerelease-checks-20261010/results-v1/plant/summary.json')},
        'source_job_sample_prefix':job['sample_prefix'],
        'current_production_panel':fp(repo/'panels/angiospermae.yaml')}
# The derived present-panel inventory must also agree for all eight unchanged samples.
for accession,s in samples.items():
    if accession.endswith('_corrected') or accession in ('SRR25378184','ERR14009273'):
        continue
    assert identities(s['expected']['candidate']['records'])==identities([r for r in current_inventory if r['accession']==accession])

builds=read(repo/'tmp/inclusive-release-20260929/build-v1/builds.json')
baseline_build=read(builds['versions']['baseline']['build']['path'])
baseline_snapshot=read(builds['versions']['baseline']['snapshot']['path'])
assert builds['versions']['baseline']['commit']=='5a664680c91ad59f32b8c2a847b8fef37f34a0ae'
assert fp(builds['versions']['baseline']['binary']['path'])==builds['versions']['baseline']['binary']
assert fp(builds['versions']['baseline']['build']['path'])==builds['versions']['baseline']['build']
assert fp(builds['versions']['baseline']['snapshot']['path'])==builds['versions']['baseline']['snapshot']
candidate_validation_path=repo/'benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/validation/result.json'
validation=read(candidate_validation_path)
assert validation['status']=='passed' and validation['source_before']==validation['source_after'] and validation['source_unchanged']
assert all(hashlib.sha256((repo/p).read_bytes()).hexdigest()==h for p,h in validation['source_after'].items())
candidate_binary=fp(repo/'tmp/annotation-conventions-20261108/run-v3/sharkmer')
assert candidate_binary['sha256']==validation['candidate_sha256']
result={
    'schema_version':1,
    'purpose':'Frozen input and independently cross-checked sequence/class oracles for final v3.2 regression check; not biological truth.',
    'preparation_scope':'Small receipts/panels/oracles and executables hashed; large FASTQ files only checked for existence and size. Runner must verify full input fingerprints and prefix counts before execution.',
    'record_order':'Preserved from saved receipt outputs, in per-target pcr_results order and emitted FASTA order; classes normalize schema-3 warning to legacy uncertain for comparisons.',
    'panel_provenance':{'historical_panels': historical_panels,
        'historical_inventory': fp(repo/'benchmarks/benchmark_results/v3.2-inclusive-release-20260929/inputs/inventory.json'),
        'corrected_panel_sha256': plant_checksums['angiospermae-after.yaml'],
        'corrected_proposal_checksums': fp(repo/'benchmarks/benchmark_results/v3.2-plant-assay-correction-20261010/proposal-v2/checksums.json')},
    'baseline': dict(builds['versions']['baseline'],cargo_version=baseline_build['cargo_version'],rustc_version=baseline_build['rustc_version'],build_command=baseline_build['command']),
    'candidate':{'binary':candidate_binary,'validation_receipt':fp(candidate_validation_path),'source_bindings':validation['source_after']},
    'samples':samples,
    'verification':{'historical_samples':10,'corrected_plant_samples':2,'historical_candidate_records':326,'historical_release_records':86,'current_combined_records':287,
        'cross_checks':['All historical receipt sequence/class identities equal tracked 326-record candidate and 86-record release inventories.','Two corrected plant receipt inventories plus eight unchanged sample inventories equal tracked 287-record combined inventory.','Every saved sequence literal hashes to its recorded SHA-256 and has expected length.','All five historical panels equal their frozen original input inventory fingerprints; corrected panel equals the reviewed proposal checksum and current production panel.','All ten FASTQ paths exist with recorded byte sizes; no large input hash or assembly run performed.','Saved candidate source hashes match current production/test sources and its binary hash matches passed frozen validation.','Baseline binary/build/snapshot fingerprints agree with preserved v3.1 build index.']},
    'source_bindings':source_bindings}
out.write_text(json.dumps(result,indent=2)+'\n')
print('Wrote',out.relative_to(repo),out.stat().st_size,'bytes')
for key,s in samples.items():
 print(key,s['sample_prefix'],{role:len(v['records']) for role,v in s['expected'].items()})
print('SHA256',hashlib.sha256(out.read_bytes()).hexdigest())

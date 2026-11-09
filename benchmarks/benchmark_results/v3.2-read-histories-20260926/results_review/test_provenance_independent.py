import importlib.util,copy,json,random,tempfile,sys,contextlib,io,itertools
from pathlib import Path
from unittest.mock import patch
root=Path(__file__).resolve().parents[4]
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
check=load('histories_post_auditor',root/'benchmarks/benchmark_results/v3.2-read-histories-20260926/results_review/provenance_ledger_review.py')
fixtures=load('history_fixture_helpers',root/'tests/test_read_histories.py');scanner=fixtures.audit
runner=load('history_runner',root/'benchmarks/benchmark_results/v3.2-read-histories-20260926/run_experiment.py')
rng=random.Random(620)
for case in range(200):
 links={str(i):{(j,j+10) for j in range(4) if rng.randrange(2)} for i in range(4)}
 pairs=[(i,j) for i,values in links.items() for j in values];expected=0
 for n in (1,2,3):
  if any(len({p[0] for p in combo})==n and len({p[1] for p in combo})==n for combo in itertools.combinations(pairs,n)):expected=n
 assert check.matching_size_three(links)==expected
with tempfile.TemporaryDirectory(prefix='history-post-auditor-synthetic-') as d:
 tmp=Path(d);evidence=tmp/'evidence';inputs=evidence/'inputs';inputs.mkdir(parents=True);scripts=tmp/'scripts';scripts.mkdir();source=scripts/'audit_read_histories.py';source.write_bytes((root/'scripts/audit_read_histories.py').read_bytes())
 short,long=fixtures.SHORT,fixtures.LONG
 seqs=[b'A'*i+short+b'C'*(3-i) for i in range(3)]+[b'G'*i+long+b'T'*(3-i) for i in range(3)]
 seqs += [scanner.reverse_complement(sequence) for sequence in seqs[:3]]
 seqs += [short[:40]+b'N'+short[41:],short,fixtures.dna(150,671)]
 qualities=[b'I'*len(s) for s in seqs];qualities[-2]=b'I'*40+b'!'+b'I'*(len(short)-41)
 raw=b''.join(b'@r'+str(i).encode()+b'\n'+s+b'\n+\n'+q+b'\n' for i,(s,q) in enumerate(zip(seqs,qualities),1));fastq=tmp/'synthetic.fastq';fastq.write_bytes(raw)
 expected={'records':len(seqs),'bases':sum(map(len,seqs)),'size_bytes':len(raw),'sha256':check.sha(raw)}
 manifest=fixtures.manifest([short,long]);(inputs/'its_manifest.json').write_text(json.dumps(manifest));(inputs/'expected_prefix.json').write_text(json.dumps(expected))
 documents=[fixtures.graph_document([short,long]),fixtures.graph_document([short,long],omitted=[short[10:29]])]
 for name,document in zip(('pre','post'),documents):(inputs/f'graph_{name}.json').write_text(json.dumps(document))
 graph_sha={name:check.identity(inputs/f'graph_{name}.json')['sha256'] for name in ('pre','post')}
 run=tmp/'run';run.mkdir();commands=[]
 for phase in ('pre','post'):
  commands.append({'id':'histories_'+phase,'argv':[sys.executable,str(source),'--fastq',str(fastq),'--graph',str(inputs/f'graph_{phase}.json'),'--manifest',str(inputs/'its_manifest.json'),'--expected-prefix',str(inputs/'expected_prefix.json'),'--output',str(run/f'histories_{phase}.json')]})
 sources={str(path.relative_to(tmp)):runner.fingerprint(path) for path in [source,*inputs.iterdir()]};snap=run/'source_snapshot'
 for relative in sources:
  path=snap/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes((tmp/relative).read_bytes())
 plan={'schema_version':1,'python':sys.version,'root':str(tmp),'run_dir':str(run),'fastq':str(fastq),'input_identity':runner.fingerprint(fastq),'sources':sources,'source_snapshot_root':str(snap),'timeout_seconds_per_command':10,'address_space_bytes':1024*1024*1024,'commands':commands}
 freeze=tmp/'freeze.json';freeze.write_text(json.dumps(plan))
 with patch.object(runner,'ROOT',tmp),contextlib.redirect_stdout(io.StringIO()):
  code=runner.execute(freeze)
 if code:
  print((run/'execution.json').read_text());print((run/'histories_pre.stderr.txt').read_text());raise AssertionError('synthetic scanner execution failed')
 reports=[json.loads((run/f'histories_{phase}.json').read_text()) for phase in ('pre','post')]
 archived={'q20_exact_records':{row['id']:row['q20_exact_records'] for row in reports[0]['classes']}}
 with patch.object(check,'ROOT',tmp),patch.object(check,'EVIDENCE',evidence),patch.object(check,'EXPECTED_PREFIX',expected),patch.object(check,'GRAPH_SHA',graph_sha),patch.object(check,'archived_exact_counts',return_value=archived):
  review=check.audit(freeze);assert review['status']=='pass'
  models=[check.graph_model(doc) for doc in documents];classes=[{**item,'region_id':region['id']} for region in manifest['regions'] for item in region['classes']]
  represented=[{**item,'representations':check.representations(item['sequence'],models[0])} for item in classes]
  rejected=0
  for mutation in ('run_start','edge_repeat','cycle_position','raw_binding','ledger_omission','eligibility','out_of_range_ordinal','duplicate_class_summary'):
   changed=copy.deepcopy(reports);record=changed[0]['records'][0]
   if mutation=='run_start':record['histories']['+']['q20']['runs'][0]['start']+=1
   elif mutation=='edge_repeat':record['histories']['+']['q20']['runs'][0]['repeated_edge_occurrences']+=1
   elif mutation=='cycle_position':record['histories']['+']['q20']['runs'][0]['cyclic_node_visits'][0]['read_node_start']+=1
   elif mutation=='raw_binding':record['raw_fastq']=record['raw_fastq'].replace('@r1','@changed')
   elif mutation=='ledger_omission':changed[0]['records'].pop(0)
   elif mutation=='eligibility':changed[0]['classes'][0]['experimental_local_eligibility']=not changed[0]['classes'][0]['experimental_local_eligibility']
   elif mutation=='out_of_range_ordinal':record['ordinal']=expected['records']+1
   elif mutation=='duplicate_class_summary':changed[0]['classes'].append(copy.deepcopy(changed[0]['classes'][0]))
   try:
    if mutation in ('raw_binding','ledger_omission'):check.prefix_binding(fastq,changed,models,classes)
    else:check.validate_report(changed[0],models[0],represented)
    raise AssertionError('mutation accepted '+mutation)
   except ValueError:rejected+=1
  assert rejected==8
  path=run/'histories_pre.json';path.write_bytes(path.read_bytes()+b' ')
  try:check.audit(freeze);raise AssertionError('outputhash mutation accepted')
  except ValueError as error:assert 'Output hash differs' in str(error)
print(json.dumps({'status':'pass','synthetic_end_to_end_graphs':2,'synthetic_prefix_records':len(seqs),'independent_matching_networks':200,'mutation_rejections':9,'checker':check.identity(Path(check.__file__))}))
m=json.loads((root/'benchmarks/benchmark_results/v3.2-read-histories-20260926/inputs/its_manifest.json').read_text());cs=[{**c,'region_id':r['id']} for r in m['regions'] for c in r['classes']];a=check.archived_exact_counts(cs);assert set(a['q20_exact_records'].values())=={0};print(json.dumps({'archived_exact_controls_verified':a}))

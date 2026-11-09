#!/usr/bin/env python3
"""Bounded BLAST of 41 previously unassessed new-target identities."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OLD=ROOT/'benchmarks/benchmark_results/v3.2-sequence-first-retention-20261001/blast-gains-v3'
BLAST=Path('/nfs/roberts/project/pi_cwd7/cwd7/ycrc_conda/envs/sharkmer-bench/bin/blastn')


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output=HERE/'blast-missing-v1';output.mkdir()
    inventory=HERE.parent/'matrix/novel_target_inventory.json'
    candidates=json.load(open(inventory));old=json.load(open(OLD/'per_sequence.json'))
    bykey={(x['accession'],x['gene'],x['sha256']):x for x in old}
    missing=[]
    for c in candidates:
        assert hashlib.sha256(c['sequence'].encode()).hexdigest()==c['sha256']
        oldrow=bykey.get((c['accession'],c['gene'],c['sha256']))
        if oldrow:assert oldrow['query_id']==c['existing_blast_query_id']
        else:missing.append({**c,'query_id':f'n{len(missing):03d}'})
    assert len(candidates)==154 and len(missing)==41
    with (output/'queries.json').open('x') as f:f.write(json.dumps(missing,indent=2,sort_keys=True)+'\n')
    query_paths={}
    for scope in ['public','angiospermae','insecta']:
        query_paths[scope]=output/(scope+'_queries.fasta')
        with query_paths[scope].open('x') as f:
            for c in missing:
                if scope=='public' or c['panel']==scope:f.write('>'+c['query_id']+'\n'+c['sequence']+'\n')
    summary=json.load(open(OLD/'summary.json'))
    for item in [summary['catalog'],summary['source_fasta'],summary['source_hits'],summary['query_fasta']]:
        assert digest(ROOT/item['path'])==item['sha256']
    sources=[inventory,OLD/'summary.json',OLD/'per_sequence.json',OLD/'public_sources_metadata.json',BLAST,Path(__file__)]
    sources+=list(OLD.glob('*_db.*'))+list(query_paths.values())+[output/'queries.json']
    before={str(p):digest(p) for p in sources}
    assert len(os.sched_getaffinity(0))==1
    version=subprocess.check_output([str(BLAST),'-version'],text=True)
    assert version.startswith('blastn: 2.17.0+')
    deadline=time.monotonic()+600;commands=[];results=[]
    fields='qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore gaps qcovs'
    for scope in query_paths:
        db=OLD/('public_sources_db' if scope=='public' else scope+'_references_db')
        remaining=max(1,int(deadline-time.monotonic()))
        command=['/usr/bin/timeout','--kill-after=5s',str(remaining)+'s','/usr/bin/time','-o',str(output/(scope+'_time.json')),
                 '-f','{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}',
                 '/usr/bin/prlimit','--as=2147483648','--',str(BLAST),'-query',str(query_paths[scope]),'-db',str(db),
                 '-task','blastn','-dust','no','-evalue','1e-10','-max_target_seqs','121','-max_hsps','20',
                 '-num_threads','1','-outfmt','6 '+fields,'-out',str(output/(scope+'_hits.tsv'))]
        commands.append(command)
        with (output/(scope+'_stdout.log')).open('x') as out,(output/(scope+'_stderr.log')).open('x') as err:
            process=subprocess.run(command,stdout=out,stderr=err)
        results.append(process.returncode)
        if process.returncode:break
    after={str(p):digest(p) for p in sources}
    receipt={'status':'pass' if len(results)==3 and not any(results) and before==after else 'fail',
             'commands':commands,'returncodes':results,'source_hashes_before':before,'source_hashes_after':after,
             'sources_unchanged':before==after,'cpu_affinity':sorted(os.sched_getaffinity(0)),
             'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'wall_budget_seconds':600,'address_space_limit_bytes':2147483648,
             'blastn_version':version,'inventory_candidates':154,'reused_queries':113,'new_queries':41,
             'outputs':{p.name:{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(output.iterdir()) if p.is_file()}}
    with (output/'execution.json').open('x') as f:f.write(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:receipt[key] for key in ['status','reused_queries','new_queries','returncodes']}))
    if receipt['status']!='pass':sys.exit(1)


if __name__=='__main__':main()

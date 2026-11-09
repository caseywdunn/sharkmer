#!/usr/bin/env python3
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if len(os.sched_getaffinity(0))!=1 or os.environ.get('SLURM_JOB_ID')!='28875555':
        raise ValueError('Expected existing job28875555 with one CPU')
    output=HERE/'run-v1';output.mkdir()
    sources=[HERE/name for name in ['count_k19.py','test_count_k19.py','run_check.py']]
    sources += [HERE.parent/name for name in ['manifest.json','expected_prefix.json','run-v1/execution.json']]
    before={str(path):digest(path) for path in sources}
    test_command=[sys.executable,'-m','unittest','discover','-s',str(HERE),'-p','test_count_k19.py','-v']
    tests=subprocess.run(test_command,text=True,capture_output=True)
    with (output/'tests.log').open('x') as stream:stream.write(tests.stdout+tests.stderr)
    if tests.returncode:raise ValueError('Synthetic k-mer controls failed')
    fastq=HERE.parents[3]/'tmp/uncertain-products-20260927/reads/SRR31887760.first1M_R1.retry1.fastq'
    command=['/usr/bin/timeout','--kill-after=5s','600s','/usr/bin/time','-o',str(output/'time.json'),
             '-f','{"wall_seconds":%e,"user_seconds":%U,"system_seconds":%S,"max_rss_kib":%M,"exit_code":%x}',
             '/usr/bin/prlimit','--as=2147483648','--',sys.executable,str(HERE/'count_k19.py'),
             '--fastq',str(fastq),'--output',str(output/'results')]
    start=time.monotonic();started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (output/'stdout.log').open('x') as out,(output/'stderr.log').open('x') as err:
        completed=subprocess.run(command,stdout=out,stderr=err)
    after={str(path):digest(path) for path in sources}
    receipt={'status':'pass' if completed.returncode==0 and before==after else 'fail',
             'returncode':completed.returncode,'command':command,'test_command':test_command,'test_exit_code':tests.returncode,
             'started_utc':started,'observer_elapsed_seconds':time.monotonic()-start,
             'source_hashes_before':before,'source_hashes_after':after,'sources_unchanged':before==after,
             'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'cpu_affinity':sorted(os.sched_getaffinity(0)),
             'address_space_limit_bytes':2147483648,'time_limit_seconds':600,
             'outputs':{str(path.relative_to(output)):{'sha256':digest(path),'bytes':path.stat().st_size}
                        for path in sorted(output.rglob('*')) if path.is_file()}}
    with (output/'execution.json').open('x') as stream:stream.write(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:receipt[key] for key in ['status','returncode','observer_elapsed_seconds','sources_unchanged']}))
    if receipt['status']!='pass':raise SystemExit(1)


if __name__=='__main__':main()

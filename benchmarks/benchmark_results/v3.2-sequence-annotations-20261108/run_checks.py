import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

repo = Path.cwd()
destination = Path(sys.argv[1]).resolve()
destination.mkdir(parents=True, exist_ok=False)
python = sys.executable
target = Path(os.environ['CARGO_TARGET_DIR'])
candidate = destination / 'sharkmer'
source_paths = sorted({p for folder in ('src', 'scripts', 'tests')
                       for p in (repo / folder).rglob('*')
                       if p.suffix in ('.rs', '.py')} |
                      {repo / 'Cargo.toml', repo / 'Cargo.lock'})

def bindings():
    return {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_paths}

result = {'started': time.time(), 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
          'source_before': bindings(), 'gates': [], 'status': 'running'}

def save():
    (destination / 'result.json').write_text(json.dumps(result, indent=2) + '\n')

gates = [
    ('fmt', ['cargo', 'fmt', '--all', '--check']),
    ('test_default', ['cargo', 'test', '--release', '--locked', '--offline']),
    ('clippy_default', ['cargo', 'clippy', '--locked', '--offline', '--', '-D', 'warnings']),
    ('test_fxhash', ['cargo', 'test', '--release', '--locked', '--offline', '--no-default-features', '--features', 'fxhashmap']),
    ('clippy_fxhash', ['cargo', 'clippy', '--locked', '--offline', '--no-default-features', '--features', 'fxhashmap', '--', '-D', 'warnings']),
    ('python', [python, '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_*.py']),
    ('references', [python, 'scripts/audit_panel_references.py', '--verify-existing']),
    ('parity', [python, 'benchmarks/benchmark_results/v3.2-sequence-annotations-20261108/parity.py', '--candidate', str(candidate), '--output', str(destination / 'parity')]),
]
save()
for name, command in gates:
    started = time.time()
    print(f'START {name}', flush=True)
    with (destination / f'{name}.log').open('w') as stream:
        try:
            outcome = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                     timeout=900, check=False)
            code = outcome.returncode
        except subprocess.TimeoutExpired:
            code = 124
    if name == 'test_default' and code == 0:
        shutil.copy2(target / 'release/sharkmer', candidate)
        result['candidate_sha256'] = hashlib.sha256(candidate.read_bytes()).hexdigest()
    result['gates'].append({'name': name, 'command': command,
                            'exit_code': code, 'elapsed_seconds': time.time() - started})
    save()
    print(f'END {name} exit={code}', flush=True)
result['source_after'] = bindings()
result['source_unchanged'] = result['source_before'] == result['source_after']
result['finished'] = time.time()
result['status'] = 'passed' if result['source_unchanged'] and all(
    gate['exit_code'] == 0 for gate in result['gates']) else 'failed'
save()
print(result['status'], flush=True)
raise SystemExit(0 if result['status'] == 'passed' else 1)

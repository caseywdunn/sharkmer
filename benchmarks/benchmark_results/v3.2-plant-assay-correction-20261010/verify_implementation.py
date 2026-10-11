#!/usr/bin/env python3
"""Freeze focused implementation checks in a fresh receipt directory."""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    output = HERE / 'implementation-v1'
    if output.exists():
        raise SystemExit('implementation-v1 must not already exist')
    output.mkdir()
    owned = [ROOT / 'panels/angiospermae.yaml', ROOT / 'tests/test_public_reference_replacements.py', ROOT / 'CHANGELOG.md', Path(__file__)]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in owned}
    assert (ROOT / 'panels/angiospermae.yaml').read_bytes() == (HERE / 'proposal-v2/angiospermae-after.yaml').read_bytes()
    command = [sys.executable, '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_public_reference_replacements.py', '-v']
    start = time.monotonic()
    result = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60)
    elapsed = time.monotonic() - start
    (output / 'tests.log').write_text(result.stdout)
    assert result.returncode == 0, result.stdout
    controls = []
    for script, arguments, directory in [
        ('prepare_assay.py', ['--output', str(HERE / 'proposal-v1')], HERE / 'proposal-v1'),
        ('amend_citation.py', [], HERE / 'proposal-v2')]:
        before = {p.name: digest(p) for p in directory.iterdir() if p.is_file()}
        refusal = subprocess.run([sys.executable, str(HERE / script), *arguments], cwd=ROOT, capture_output=True, text=True, timeout=60)
        after = {p.name: digest(p) for p in directory.iterdir() if p.is_file()}
        assert refusal.returncode != 0 and 'must not already exist' in refusal.stderr and before == after
        controls.append({'script': script, 'exit_code': refusal.returncode, 'existing_output_unchanged': True})
    assert hashes == {str(p.relative_to(ROOT)): digest(p) for p in owned}
    receipt = {'status': 'pass', 'command': command, 'exit_code': result.returncode, 'elapsed_seconds': elapsed,
               'source_hashes': hashes, 'tests_log_sha256': digest(output / 'tests.log'),
               'production_panel_equals_reviewed_proposal_v2': True, 'overwrite_refusal_controls': controls,
               'scope': 'Focused panel/reference regression checks; full CLI reconstruction and real plant evaluation have separate execution receipts.'}
    (output / 'verification.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()

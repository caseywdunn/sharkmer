#!/usr/bin/env python3
"""Diagnostic-only warm-cache probes after a failed test fixture; not final clearance."""
import argparse
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from build_snapshot import fingerprint, write_new
import run_compiled_checks as checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cargo', type=Path, required=True)
    parser.add_argument('--rustc', type=Path, required=True)
    parser.add_argument('--warm-target', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'failed', 'purpose': 'diagnostic_only_final_full_checks_still_required', 'warm_target': str(args.warm_target.resolve()), 'steps': [], 'probe_source': fingerprint(__file__)}
    try:
        before = checks.sources()
        receipt['sources_before'] = before
        receipt['cargo'] = fingerprint(args.cargo.resolve())
        receipt['rustc'] = fingerprint(args.rustc.resolve())
        environment = dict(os.environ, RUSTC=str(args.rustc.resolve()), CARGO_TARGET_DIR=str(args.warm_target.resolve()), CARGO_BUILD_JOBS='2', LC_ALL='C')
        for key in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'RUSTC_WRAPPER', 'RUSTC_WORKSPACE_WRAPPER'):
            environment.pop(key, None)
        cargo = str(args.cargo.resolve())
        commands = [('cli-unit', [cargo, 'test', '--locked', '--offline', '--release', '--bin', 'sharkmer', 'cli::tests::uncertain_output_name_collisions_fail_before_input_resolution_only_when_enabled', '--', '--exact']), ('uncertain-integration', [cargo, 'test', '--locked', '--offline', '--release', '--test', 'uncertain_products'])]
        for name, command in commands:
            row = checks.run_step(name, command, checks.ROOT, environment, output)
            receipt['steps'].append(row)
            if row['status'] != 'pass':
                raise ValueError(f'Probe failed: {name}')
        receipt['sources_after'] = checks.sources()
        if before != receipt['sources_after']:
            raise ValueError('Source changed during targeted probe')
        receipt['status'] = 'pass'
    except Exception as error:
        receipt['error'] = repr(error)
    finally:
        write_new(output / 'probe.json', receipt)
    return 0 if receipt['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

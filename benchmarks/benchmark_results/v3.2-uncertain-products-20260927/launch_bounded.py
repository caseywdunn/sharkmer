#!/usr/bin/env python3
"""Observe enforced process constraints immediately before replacing self with sharkmer."""
import argparse
import json
import os
from pathlib import Path
import resource
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--expected-cpus', required=True)
    parser.add_argument('--expected-as', type=int, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    expected_cpus = sorted(int(item) for item in args.expected_cpus.split(','))
    cpus = sorted(os.sched_getaffinity(0))
    limits = list(resource.getrlimit(resource.RLIMIT_AS))
    passed = bool(command) and cpus == expected_cpus and limits == [args.expected_as] * 2
    receipt = {'status': 'pass' if passed else 'fail', 'pid': os.getpid(), 'unix_time': time.time(), 'affinity': cpus, 'rlimit_as_soft_hard': limits, 'command': command}
    with args.receipt.open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    if not passed:
        raise ValueError('Observed process limits differ from the frozen plan')
    os.execv(command[0], command)


if __name__ == '__main__':
    main()

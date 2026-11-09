#!/usr/bin/env python3
"""Compare optimized exact evidence with direct scanner on adversarial reads."""

import hashlib
import json
import random
import time
from pathlib import Path

from scan_prefiltered import ExactPrefilter, HERE, load_scanner


def main():
    scanner = load_scanner()
    events = scanner.prepare_events(json.loads((HERE / "manifest.json").read_text()), 0)
    direct = scanner.observations_for_read
    adapter = ExactPrefilter(scanner)
    cases = [(b"N" * 101, b"I" * 101)]
    # Both orientations, overlapping/repeated seeds, whole-window Q20 boundary,
    # a failed seed retained as a false-positive filter hit, and outside-window
    # low qualities. Patterns >101 are tested with synthetic longer reads only.
    for event in events:
        for _, pattern in event["patterns"]:
            cases.extend([(b"NN" + pattern + b"NN", b"!" * 2 + b"5" * len(pattern) + b"!" * 2),
                          (pattern, b"4" + b"I" * (len(pattern) - 1)),
                          (pattern[:19] + b"N" + pattern[20:], b"I" * len(pattern))])
    randomizer = random.Random(20261009)
    cases.extend((bytes(randomizer.choice(b"ACGTN") for _ in range(101)), b"I" * 101) for _ in range(1000))
    start = time.perf_counter()
    for sequence, quality in cases:
        if adapter(sequence, quality, events, 0) != direct(sequence, quality, events, 0):
            raise ValueError("Optimized and direct observations differ")
    equality_seconds = time.perf_counter() - start
    try:
        adapter(b"A" * 101, b"I" * 101, events, 1)
    except ValueError:
        pass
    else:
        raise ValueError("Nonexact assay was accepted")
    # A synthetic random-read runtime sample, without reading experimental data.
    random_reads = [bytes(randomizer.choice(b"ACGT") for _ in range(101)) for _ in range(2000)]
    start = time.perf_counter()
    for sequence in random_reads:
        adapter(sequence, b"I" * 101, events, 0)
    sample_seconds = time.perf_counter() - start
    result = {"status": "pass", "exact_observation_equality_cases": len(cases),
              "equality_seconds": equality_seconds, "synthetic_random_reads": len(random_reads),
              "synthetic_random_seconds": sample_seconds,
              "random_only_1M_projection_seconds": sample_seconds * 500,
              "projection_limitation": "Excludes FASTQ IO, parsing, hashing and experimental matching density; execution remains bounded by 600 seconds",
              "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                  [HERE / "scan_prefiltered.py", Path(__file__), HERE / "manifest.json"]}}
    with (HERE / "prefilter_controls.json").open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()

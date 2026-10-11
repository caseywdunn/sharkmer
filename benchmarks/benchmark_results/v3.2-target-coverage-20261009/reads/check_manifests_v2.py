#!/usr/bin/env python3
"""Validate frozen real patterns on synthetic reads without sampling real inputs."""

import argparse
import json
from pathlib import Path
import random
import time

from scan_exact_v2 import Matcher, prepare_placements, rc, sha


HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifests", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = sorted(args.manifests.glob("*.json"))
    if len(paths) != 4:
        raise ValueError("Expected four frozen sample manifests")
    results = []
    for path in paths:
        manifest = json.loads(path.read_text())
        events = manifest["events"]
        matcher = Matcher(events)
        prepare_placements(manifest)
        candidates = {x["id"]: x for x in manifest["candidates"]}
        tile_spans = {x: [] for x in manifest["focused_candidate_ids"] + manifest["positive_control_ids"]}
        cases = 0
        for event in events:
            for membership in event["memberships"]:
                candidate = candidates[membership["candidate_id"]]
                piece = candidate["sequence"][membership["start"]:membership["end"]].encode()
                if min(piece, rc(piece)).decode() != event["sequence"]:
                    raise ValueError("Membership does not reproduce pattern")
                if membership["role"] in ("candidate_backbone_tile", "positive_control_tile"):
                    tile_spans[membership["candidate_id"]].append((membership["start"], membership["end"]))
            pattern = event["sequence"].encode()
            for oriented in (pattern, rc(pattern)):
                sequence = b"NN" + oriented + b"NN"
                for quality in (b"!!" + b"5" * len(oriented) + b"!!", b"!!4" + b"I" * (len(oriented) - 1) + b"!!"):
                    fast = matcher.observations(sequence, quality)
                    direct = matcher.observations(sequence, quality, prefilter=False)
                    if fast != direct:
                        raise ValueError("Seeded matching omits exact observations")
                    cases += 1
        for identifier, spans in tile_spans.items():
            covered = {i for start, end in spans for i in range(start, end)}
            if covered != set(range(len(candidates[identifier]["sequence"]))):
                raise ValueError("Tile schedule does not cover entire candidate")
        randomizer = random.Random(20261009)
        sequences = [bytes(randomizer.choice(b"ACGT") for _ in range(151)) for _ in range(2000)]
        quality = b"I" * 151
        started = time.monotonic()
        for sequence in sequences:
            matcher.observations(sequence, quality)
        elapsed = time.monotonic() - started
        results.append({"accession": manifest["accession"], "events": len(events),
                        "synthetic_pattern_orientation_quality_cases": cases,
                        "candidate_tile_schedules_verified": len(tile_spans),
                        "synthetic_151base_reads": len(sequences), "random_matcher_seconds": elapsed,
                        "random_only_million_read_projection_seconds": elapsed * 500,
                        "projection_limit": "Excludes real match density, FASTQ parsing, quality checks, hashing and output; hard execution limits remain necessary"})
    source_paths = paths + [Path(__file__), HERE / "scan_exact_v2.py", HERE / "prepare_assays_v2.py", HERE / "run_assay_v2.py"]
    result = {"status": "pass", "samples": results,
              "sources": {str(p.resolve()): sha(p.read_bytes()) for p in source_paths}}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result["samples"]))


if __name__ == "__main__":
    main()

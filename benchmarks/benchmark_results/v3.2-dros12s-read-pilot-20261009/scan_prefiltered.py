#!/usr/bin/env python3
"""Exact-only seed prefilter around the unchanged read bridge scanner.

Every exact window occurrence necessarily contains its first 19 bases in
the same orientation. The prefilter only omits events for which neither
orientation's prefix seed occurs. Matching/quality/placement semantics are
delegated to the archived scanner implementation, without approximation.
"""

import argparse
import hashlib
import importlib.util
import json
import re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCANNER = ROOT / "scripts/audit_read_bridges.py"


def digest(value):
    return hashlib.sha256(value).hexdigest()


def load_scanner():
    spec = importlib.util.spec_from_file_location("unchanged_scanner", SCANNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExactPrefilter:
    def __init__(self, scanner):
        self.scanner = scanner
        self.original = scanner.observations_for_read
        self.events = None
        self.raw_counts = defaultdict(int)
        self.prefiltered_records = 0
        self.records = 0

    def bind(self, events):
        self.events = events
        self.seed_events = defaultdict(set)
        for i, event in enumerate(events):
            for _, pattern in event["patterns"]:
                if len(pattern) < 19:
                    raise ValueError("Prefilter requires at least 19 bases")
                self.seed_events[pattern[:19]].add(i)
        # Every alternative is a literal equal-length seed. A lookahead retains
        # overlapping seed occurrences and returns the union of eligible events.
        self.regex = re.compile(b"(?=(" + b"|".join(sorted(self.seed_events)) + b"))")

    def __call__(self, sequence, quality, events, max_substitutions=0):
        if max_substitutions != 0:
            raise ValueError("Seed filtering is valid only for exact matching")
        if events is not self.events:
            self.bind(events)
        self.records += 1
        indices = set()
        for hit in self.regex.finditer(sequence):
            indices.update(self.seed_events[hit[1]])
        if not indices:
            return []
        self.prefiltered_records += 1
        selected = [events[i] for i in sorted(indices)]
        # Supplemental exact occurrence without a base-quality threshold.
        # This does not change, expand, or replace the primary Q20 assay.
        for event in selected:
            if any(pattern in sequence for _, pattern in event["patterns"]):
                self.raw_counts[event["id"]] += 1
        return self.original(sequence, quality, selected, max_substitutions=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing output overwrite")
    paths = [SCANNER, Path(__file__), HERE / "manifest.json", HERE / "expected_prefix.json", HERE / "derivation.json"]
    frozen = {str(p): digest(p.read_bytes()) for p in paths}
    scanner = load_scanner()
    adapter = ExactPrefilter(scanner)
    scanner.observations_for_read = adapter
    result = scanner.scan_fastq(args.fastq, json.loads((HERE / "expected_prefix.json").read_text()),
                                json.loads((HERE / "manifest.json").read_text()), 0)
    for event in result["events"]:
        event["unfiltered_exact_records"] = adapter.raw_counts[event["id"]]
        if event["unfiltered_exact_records"] < event["records"]:
            raise ValueError("Filtered count exceeds raw count")
        event["assay_callability"] = "longer_than_every_read" if event["longer_than_every_read"] else event["callability"]
    result["prefilter"] = {"seed_length": 19, "distinct_oriented_seeds": len(adapter.seed_events),
                           "records_examined": adapter.records, "records_with_any_seed": adapter.prefiltered_records,
                           "exact_only": True, "semantic_delegate": str(SCANNER)}
    if frozen != {str(p): digest(p.read_bytes()) for p in paths}:
        raise ValueError("Assay sources changed")
    result["provenance"] = frozen
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"events": len(result["events"]), "matching_q20_records": len(result["read_evidence"]), "prefilter": result["prefilter"]}))


if __name__ == "__main__":
    main()

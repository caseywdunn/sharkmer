#!/usr/bin/env python3
"""Bounded exact local-string counts with a separately capped example ledger."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re


COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
MIN_QUALITY = 20
SEED_LENGTH = 19
LEDGER_READS = 5000
LEDGER_BASES = 1000000
LEDGER_MATCHES = 100000
EXAMPLES_Q20 = 8
EXAMPLES_RAW_ONLY = 2
MAX_READ_LENGTH = 10000


def sha(data):
    return hashlib.sha256(data).hexdigest()


def rc(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def positions(sequence, pattern):
    start = sequence.find(pattern)
    while start >= 0:
        yield start
        start = sequence.find(pattern, start + 1)


def oriented_patterns(pattern):
    other = rc(pattern)
    return [("both", pattern)] if other == pattern else [("+", pattern), ("-", other)]


def seed_trie_regex(seeds):
    """Factor equal-length literal seeds without changing overlapping matches."""
    root = {}
    for seed in seeds:
        node = root
        for letter in seed:
            node = node.setdefault(letter, {})
    def emit(node):
        branches = [bytes([letter]) + emit(child) for letter, child in sorted(node.items())]
        if not branches:
            return b""
        if len(branches) == 1:
            return branches[0]
        if all(len(branch) == 1 for branch in branches):
            return b"[" + b"".join(branches) + b"]"
        return b"(?:" + b"|".join(branches) + b")"
    return re.compile(b"(?=(" + emit(root) + b"))")


class Matcher:
    def __init__(self, events):
        self.events = events
        self.patterns = []
        seeds = defaultdict(set)
        for index, event in enumerate(events):
            pattern = event["sequence"].encode("ascii")
            if len(pattern) < SEED_LENGTH or set(pattern) - set(b"ACGT"):
                raise ValueError("Only literal patterns of at least 19 bases are allowed")
            if sha(pattern) != event["sha256"]:
                raise ValueError("Pattern hash differs")
            oriented = oriented_patterns(pattern)
            self.patterns.append(oriented)
            for _, text in oriented:
                seeds[text[:SEED_LENGTH]].add(index)
        self.seeds = dict(seeds)
        self.regex = seed_trie_regex(seeds)

    def observations(self, sequence, quality, prefilter=True):
        if len(sequence) != len(quality):
            raise ValueError("Sequence/quality length mismatch")
        if prefilter:
            selected = set()
            for match in self.regex.finditer(sequence):
                selected.update(self.seeds[match[1]])
        else:
            selected = range(len(self.events))
        result = []
        for index in sorted(selected):
            matches = []
            for strand, pattern in self.patterns[index]:
                for start in positions(sequence, pattern):
                    end = start + len(pattern)
                    minimum = min(quality[start:end]) - 33
                    matches.append({"strand": strand, "start": start, "end": end,
                                    "minimum_phred": minimum, "q20": minimum >= MIN_QUALITY})
            if matches:
                result.append({"event_index": index, "matches": matches})
        return result


def prepare_placements(manifest):
    candidates = manifest["candidates"]
    if len({x["id"] for x in candidates}) != len(candidates):
        raise ValueError("Duplicate candidate IDs")
    for candidate in candidates:
        seq = candidate["sequence"].encode("ascii")
        if not seq or set(seq) - set(b"ACGT") or sha(seq) != candidate["sha256"]:
            raise ValueError("Candidate sequence identity invalid")
    placements = []
    for event in manifest["events"]:
        pattern = event["sequence"].encode("ascii")
        items = []
        for candidate in candidates:
            sequence = candidate["sequence"].encode("ascii")
            for strand, text in oriented_patterns(sequence):
                for start in positions(text, pattern):
                    items.append({"candidate_id": candidate["id"], "strand": strand,
                                  "start": start, "end": start + len(pattern)})
        if not items:
            raise ValueError("Event lacks an exact candidate placement")
        placements.append(items)
    return placements


def summarize_records(records, events, limits=None):
    """Shared streaming aggregation; complete counts are independent of ledger caps."""
    limits = limits or {"reads": LEDGER_READS, "bases": LEDGER_BASES, "matches": LEDGER_MATCHES}
    matcher = Matcher(events)
    counts = [{"raw_records": 0, "q20_records": 0, "raw_occurrences": 0,
               "q20_occurrences": 0, "raw_only_records": 0,
               "retained_q20_records": 0, "retained_raw_only_records": 0,
               "ledger_cap_omitted_requests": 0} for _ in events]
    ledger, stored_bases, stored_matches = [], 0, 0
    cap_reasons = Counter()
    read_lengths = Counter()
    matching_records = 0
    for ordinal, header, sequence, quality in records:
        read_lengths[len(sequence)] += 1
        observations = matcher.observations(sequence, quality)
        if not observations:
            continue
        matching_records += 1
        requested = []
        for observation in observations:
            i = observation["event_index"]
            matches = observation["matches"]
            q20_count = sum(x["q20"] for x in matches)
            count = counts[i]
            count["raw_records"] += 1
            count["q20_records"] += bool(q20_count)
            count["raw_only_records"] += not bool(q20_count)
            count["raw_occurrences"] += len(matches)
            count["q20_occurrences"] += q20_count
            if (q20_count and count["retained_q20_records"] < EXAMPLES_Q20) or (
                    not q20_count and count["retained_raw_only_records"] < EXAMPLES_RAW_ONLY):
                requested.append(i)
        if not requested:
            continue
        new_matches = sum(len(x["matches"]) for x in observations)
        reasons = []
        if len(ledger) >= limits["reads"]:
            reasons.append("reads")
        if stored_bases + len(sequence) > limits["bases"]:
            reasons.append("bases")
        if stored_matches + new_matches > limits["matches"]:
            reasons.append("matches")
        if reasons:
            cap_reasons.update(reasons)
            for i in requested:
                counts[i]["ledger_cap_omitted_requests"] += 1
            continue
        ledger.append({"ordinal": ordinal, "header": header.decode("ascii"),
                       "sequence": sequence.decode("ascii"), "quality": quality.decode("ascii"),
                       "observations": observations})
        stored_bases += len(sequence)
        stored_matches += new_matches
        for observation in observations:
            key = "retained_q20_records" if any(x["q20"] for x in observation["matches"]) else "retained_raw_only_records"
            counts[observation["event_index"]][key] += 1
    summaries = []
    for event, count in zip(events, counts):
        length = len(event["sequence"])
        eligible = sum(n for size, n in read_lengths.items() if size >= length)
        summaries.append({"id": event["id"], **count, "eligible_length_records": eligible,
                          "callability": "assayed" if eligible else "longer_than_every_read",
                          "example_ledger_complete_for_event":
                          count["raw_records"] == count["retained_q20_records"] + count["retained_raw_only_records"]})
    return {"events": summaries, "read_evidence": ledger,
            "records_with_any_exact_pattern": matching_records,
            "read_length_histogram": dict(sorted(read_lengths.items())),
            "ledger": {"selection": "first 8 Q20 and first 2 raw-only records per event, with shared-read bonus and global caps",
                       "limits": limits, "stored_reads": len(ledger), "stored_bases": stored_bases,
                       "stored_match_positions": stored_matches, "cap_omitted_requests_by_reason": dict(cap_reasons),
                       "complete_counts_despite_example_caps": True,
                       "no_molecule_or_footprint_uniqueness_inferred": True},
            "prefilter": {"seed_length": SEED_LENGTH, "oriented_seeds": len(matcher.seeds)}}


def scan(path, manifest):
    expected = manifest["expected_prefix"]
    path = Path(path)
    if path.is_symlink() or not path.is_file() or not 0 < expected["records"] <= 1000000:
        raise ValueError("Invalid input or prefix bound")
    before = path.stat()
    digest = hashlib.sha256()
    observed = {"records": 0, "bases": 0, "size_bytes": 0}

    def records():
        with path.open("rb") as stream:
            for ordinal in range(1, expected["records"] + 1):
                lines = [stream.readline() for _ in range(4)]
                if not all(lines):
                    raise ValueError("Truncated input")
                header, sequence, separator, quality = [x.rstrip(b"\r\n") for x in lines]
                if (not header.startswith(b"@") or not separator.startswith(b"+") or
                        not 0 < len(sequence) <= MAX_READ_LENGTH or len(sequence) != len(quality) or
                        min(quality) < 33 or max(quality) > 126):
                    raise ValueError("Malformed FASTQ or unsupported length/quality")
                content = b"".join(lines)
                digest.update(content)
                observed["records"] += 1
                observed["bases"] += len(sequence)
                observed["size_bytes"] += len(content)
                yield ordinal, header, sequence.upper(), quality

    placements = prepare_placements(manifest)
    result = summarize_records(records(), manifest["events"])
    observed["sha256"] = digest.hexdigest()
    if observed != expected:
        raise ValueError("Consumed input identity differs")
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns")
    if any(getattr(before, k) != getattr(after, k) for k in fields):
        raise ValueError("Input changed during scan")
    for summary, places in zip(result["events"], placements):
        summary["candidate_placements"] = places
    result.update({"status": "completed", "schema_version": 1, "accession": manifest["accession"],
                   "verified_prefix": observed, "minimum_phred": MIN_QUALITY,
                   "match_policy": "exact, both orientations, overlapping occurrences; one record counted once per event",
                   "palindromes": "one physical occurrence with strand=both",
                   "interpretation": "Local occurrences only; no full-product linkage, locus identity or biological validation is implied."})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    frozen = {str(p.resolve()): sha(p.read_bytes()) for p in (args.manifest, Path(__file__))}
    manifest = json.loads(args.manifest.read_text())
    result = scan(manifest["fastq_path"], manifest)
    if frozen != {p: sha(Path(p).read_bytes()) for p in frozen}:
        raise ValueError("Manifest or scanner changed")
    result["source_hashes"] = frozen
    with args.output.open("x") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps({"accession": result["accession"], "events": len(result["events"]),
                      "matching_records": result["records_with_any_exact_pattern"], "ledger": result["ledger"]}))


if __name__ == "__main__":
    main()

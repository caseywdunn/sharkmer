#!/usr/bin/env python3
"""Audit literal local read bridges without asserting full-amplicon truth."""

import argparse
import hashlib
import json
from pathlib import Path


COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
FLANK_BASES = 21
MIN_QUALITY = 20
MAX_MATCHING_RECORDS = 100000


def digest(contents):
    return hashlib.sha256(contents).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def positions(sequence, pattern):
    offset = sequence.find(pattern)
    while offset >= 0:
        yield offset
        offset = sequence.find(pattern, offset + 1)


def window_hits(sequence, pattern, max_substitutions):
    if len(sequence) < len(pattern):
        return []
    if max_substitutions == 0:
        return [(offset, 0) for offset in positions(sequence, pattern)]
    midpoint = len(pattern) // 2
    starts = set(positions(sequence, pattern[:midpoint]))
    starts.update(offset - midpoint for offset in positions(sequence, pattern[midpoint:]))
    matches = []
    for start in sorted(starts):
        if start < 0 or start + len(pattern) > len(sequence):
            continue
        observed = sequence[start:start + len(pattern)]
        if set(observed) - set(b"ACGT"):
            continue
        substitutions = sum(expected != actual for expected, actual in zip(pattern, observed))
        if substitutions <= max_substitutions:
            matches.append((start, substitutions))
    return matches


def prepare_events(manifest, max_substitutions=0):
    if max_substitutions not in (0, 1):
        raise ValueError("Only exact or separately reported one-substitution assays are supported")
    if manifest.get("schema_version") != 1 or not manifest.get("candidates") or not manifest.get("events"):
        raise ValueError("Expected version 1 candidate/event manifest")
    candidates = {}
    for candidate in manifest["candidates"]:
        identifier = candidate["id"]
        sequence = candidate["sequence"].encode("ascii")
        if identifier in candidates or not sequence or set(sequence) - set(b"ACGT"):
            raise ValueError("Duplicate candidate or nonliteral candidate sequence")
        if digest(sequence) != candidate["sequence_sha256"]:
            raise ValueError("Candidate sequence checksum differs")
        candidates[identifier] = sequence
    events = []
    identifiers = set()
    for event in manifest["events"]:
        identifier = event["id"]
        if identifier in identifiers:
            raise ValueError("Duplicate event identifier")
        identifiers.add(identifier)
        candidate = candidates[event["candidate_id"]]
        start, end = event["start"], event["end"]
        if type(start) is not int or type(end) is not int or not 0 <= start <= end <= len(candidate):
            raise ValueError("Event coordinates are invalid")
        window_start, window_end = start - FLANK_BASES, end + FLANK_BASES
        prepared = {**event, "window_start": window_start, "window_end": window_end}
        prepared["_candidate_sequences"] = candidates
        if window_start < 0 or window_end > len(candidate):
            prepared.update({"callability": "missing_complete_flanks", "patterns": []})
        else:
            pattern = candidate[window_start:window_end]
            template_placements = []
            for candidate_id, sequence in candidates.items():
                for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                    for offset, substitutions in window_hits(oriented, pattern, max_substitutions):
                        template_placements.append({"candidate_id": candidate_id, "strand": strand, "start": offset, "substitutions": substitutions})
            prepared.update({
                "callability": "assayed",
                "window_sequence": pattern.decode("ascii"),
                "template_placements": template_placements,
                "unique_template_placement": len(template_placements) == 1,
                "patterns": [("+", pattern), ("-", reverse_complement(pattern))],
            })
        events.append(prepared)
    return events


def observations_for_read(sequence, quality, events, max_substitutions=0):
    observations = []
    for event in events:
        matches = []
        for strand, pattern in event["patterns"]:
            for start, substitutions in window_hits(sequence, pattern, max_substitutions):
                end = start + len(pattern)
                if min(quality[start:end]) < MIN_QUALITY + 33:
                    continue
                oriented_start = start if strand == "+" else len(sequence) - end
                projected_start = event["window_start"] - oriented_start
                observed_window = sequence[start:end]
                if strand == "-":
                    observed_window = reverse_complement(observed_window)
                placements = []
                for candidate_id, candidate_sequence in event["_candidate_sequences"].items():
                    for template_strand, oriented in (("+", candidate_sequence), ("-", reverse_complement(candidate_sequence))):
                        for offset, differences in window_hits(oriented, observed_window, max_substitutions):
                            placements.append({"candidate_id": candidate_id, "strand": template_strand, "start": offset, "substitutions": differences})
                matches.append({
                    "strand": strand, "read_start": start, "read_end": end,
                    "substitutions": substitutions,
                    "projected_footprint": [projected_start, projected_start + len(sequence)],
                    "observed_template_placements": placements,
                })
        if matches:
            observations.append({
                "event_id": event["id"], "matches": matches,
                "candidate_id": event["candidate_id"],
                "comparison_group": event.get("comparison_group", event["candidate_id"]),
                "unique": len(matches) == 1 and len(matches[0]["observed_template_placements"]) == 1,
            })
    competing = {}
    for observation in observations:
        competing.setdefault(observation["comparison_group"], set()).add(observation["candidate_id"])
    for observation in observations:
        if len(competing[observation["comparison_group"]]) > 1:
            observation["unique"] = False
    return observations


def corroboration_witness(links):
    owners = {}

    def assign(read_id, visited):
        for footprint in sorted(links[read_id]):
            if footprint in visited:
                continue
            visited.add(footprint)
            if footprint not in owners or assign(owners[footprint], visited):
                owners[footprint] = read_id
                return True
        return False

    for read_id in sorted(links):
        assign(read_id, set())
        if len(owners) >= 3:
            break
    return [{"read_id": read_id, "footprint": list(footprint)} for footprint, read_id in sorted(owners.items())]


def scan_fastq(path, expected, manifest, max_substitutions=0):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Input must be a regular FASTQ file")
    records = expected["records"]
    if type(records) is not int or records <= 0:
        raise ValueError("An explicit positive consumed-prefix record count is required")
    events = prepare_events(manifest, max_substitutions)
    accumulators = {
        event["id"]: {"records": 0, "unique_records": 0, "ambiguous_records": 0, "exact_records": 0, "one_substitution_records": 0, "ids": set(), "sequences": set(), "footprints": set(), "strands": set(), "links": {}}
        for event in events
    }
    prefix_digest = hashlib.sha256()
    size_bytes = 0
    bases = 0
    maximum_read_length = 0
    ledger = []
    before = path.stat()
    with path.open("rb") as stream:
        for ordinal in range(1, records + 1):
            lines = [stream.readline() for _field in range(4)]
            if any(not line for line in lines):
                raise ValueError(f"Truncated FASTQ before consumed record {ordinal}")
            header, sequence, separator, quality = [line.rstrip(b"\r\n") for line in lines]
            if not header.startswith(b"@") or not separator.startswith(b"+") or not sequence or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record {ordinal}")
            if min(quality) < 33 or max(quality) > 126:
                raise ValueError(f"Unsupported quality encoding at record {ordinal}")
            contents = b"".join(lines)
            prefix_digest.update(contents)
            size_bytes += len(contents)
            bases += len(sequence)
            maximum_read_length = max(maximum_read_length, len(sequence))
            sequence = sequence.upper()
            observations = observations_for_read(sequence, quality, events, max_substitutions)
            if not observations:
                continue
            if len(ledger) >= MAX_MATCHING_RECORDS:
                raise ValueError("Matching-record evidence limit exceeded; no successful partial audit")
            read_id = header.split()[0].decode("ascii")
            sequence_digest = digest(min(sequence, reverse_complement(sequence)))
            ledger.append({"ordinal": ordinal, "header": header.decode("ascii"), "sequence": sequence.decode("ascii"), "quality": quality.decode("ascii"), "observations": observations})
            for observation in observations:
                accumulator = accumulators[observation["event_id"]]
                accumulator["records"] += 1
                accumulator["exact_records"] += any(match["substitutions"] == 0 for match in observation["matches"])
                accumulator["one_substitution_records"] += any(match["substitutions"] == 1 for match in observation["matches"])
                if observation["unique"]:
                    accumulator["unique_records"] += 1
                    accumulator["ids"].add(read_id)
                    accumulator["sequences"].add(sequence_digest)
                    accumulator["footprints"].add(tuple(observation["matches"][0]["projected_footprint"]))
                    accumulator["links"].setdefault(read_id, set()).add(tuple(observation["matches"][0]["projected_footprint"]))
                    accumulator["strands"].add(observation["matches"][0]["strand"])
                else:
                    accumulator["ambiguous_records"] += 1
    observed = {"records": records, "bases": bases, "size_bytes": size_bytes, "sha256": prefix_digest.hexdigest()}
    if any(observed[key] != expected[key] for key in observed):
        raise ValueError("Consumed-prefix identity differs; no read-support result published")
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("Input changed during audit")
    summaries = []
    for event in events:
        accumulator = accumulators[event["id"]]
        summary = {key: value for key, value in event.items() if key not in {"patterns", "_candidate_sequences"}}
        summary.update({key: value for key, value in accumulator.items() if not isinstance(value, set) and key != "links"})
        summary.update({"distinct_ids": len(accumulator["ids"]), "distinct_sequences": len(accumulator["sequences"]), "distinct_footprints": len(accumulator["footprints"]), "footprints": sorted(accumulator["footprints"]), "strands": sorted(accumulator["strands"]), "haplotype_truth": "not_established"})
        summary["corroboration_witness"] = corroboration_witness(accumulator["links"])
        corroborated = len(summary["corroboration_witness"]) >= 3
        summary["local_bridge_status"] = "corroborated_local_bridge" if corroborated else "unresolved"
        summary["longer_than_every_read"] = event["window_end"] - event["window_start"] > maximum_read_length
        summaries.append(summary)
    return {"schema_version": 1, "input_path": str(path.resolve()), "input_subset": observed, "maximum_read_length": maximum_read_length, "max_substitutions": max_substitutions, "minimum_quality": MIN_QUALITY, "flank_bases": FLANK_BASES, "paired": False, "independence": "Distinct observed R1 footprints, not UMI-certified molecules", "events": summaries, "read_evidence": ledger}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", type=Path, required=True)
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--expected-prefix", type=Path, required=True)
    parser.add_argument("--max-substitutions", type=int, choices=(0, 1), default=0)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists() or arguments.output.is_symlink():
        raise ValueError("Refusing to overwrite an audit output")
    event_bytes = arguments.events.read_bytes()
    prefix_bytes = arguments.expected_prefix.read_bytes()
    script_digest = digest(Path(__file__).read_bytes())
    result = scan_fastq(arguments.fastq, json.loads(prefix_bytes), json.loads(event_bytes), arguments.max_substitutions)
    if arguments.events.read_bytes() != event_bytes or arguments.expected_prefix.read_bytes() != prefix_bytes or digest(Path(__file__).read_bytes()) != script_digest:
        raise ValueError("Audit specification or scanner changed during evaluation")
    result["provenance"] = {"event_manifest_sha256": digest(event_bytes), "expected_prefix_manifest_sha256": digest(prefix_bytes), "scanner_sha256": script_digest}
    with arguments.output.open("x") as output:
        output.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(arguments.output), "events": len(result["events"]), "matching_records": len(result["read_evidence"])}))


if __name__ == "__main__":
    main()

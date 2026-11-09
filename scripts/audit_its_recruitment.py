#!/usr/bin/env python3
"""Bounded, standalone flank recruitment; alignments do not establish haplotypes."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time


FLANK_BASES = 21
FLANK_EDITS = 2
SPAN_EDITS = 4
MIN_QUALITY = 20
MAX_RECRUITED_RECORDS = 100000
MAX_PAIRS_PER_READ = 4096
MAX_DECOY_PLACEMENTS = 4096
MAX_DECOY_PLACEMENTS_PER_READ = 16384
MAX_MANIFEST_BYTES = 1024 * 1024
MAX_REGIONS = 16
MAX_CLASSES_PER_REGION = 16
MAX_TOTAL_CLASSES = 64
MAX_DECOYS = 32
MAX_DECOY_LENGTH = 10000
MAX_TOTAL_DECOY_BASES = 100000
MAX_SELF_CONTEXTS_PER_CLASS = 64
MAX_IDENTIFIER_BYTES = 256
MAX_READ_LENGTH = 1000
MAX_FASTQ_LINE_BYTES = 4096
MAX_FLANK_PLACEMENTS_PER_READ = 4096
MAX_RECORD_EVIDENCE_BYTES = 16 * 1024 * 1024
MAX_OUTPUT_BYTES = 1024 * 1024 * 1024
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
BASES = frozenset(b"ACGT")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded_json(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("ascii")


def bounded_json(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(MAX_MANIFEST_BYTES + 1)
    if len(raw) > MAX_MANIFEST_BYTES:
        raise ValueError("JSON manifest byte cap exceeded")
    return raw


def file_identity(path):
    """Hash the complete staged file, separately from its consumed prefix."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Input must be a regular FASTQ file")
    before = path.stat()
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, field) != getattr(after, field) for field in fields):
        raise ValueError("Input changed during staged-file hashing")
    return {"sha256": hasher.hexdigest(), "size_bytes": after.st_size,
            "device": after.st_dev, "inode": after.st_ino,
            "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns}


def publish_json(path, value, maximum=None):
    """Publish a complete fsynced file using an atomic, no-clobber hard link."""
    maximum = MAX_OUTPUT_BYTES if maximum is None else maximum
    raw = encoded_json(value)
    if len(raw) > maximum:
        raise ValueError("Serialized output byte cap exceeded; no successful partial result")
    descriptor, temporary = tempfile.mkstemp(prefix=".its-audit-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(raw)
            output.flush()
            os.fsync(output.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def valid_identifier(identifier):
    if not isinstance(identifier, str) or not identifier or len(identifier.encode("utf8")) > MAX_IDENTIFIER_BYTES:
        raise ValueError("Invalid identifier or identifier byte cap exceeded")
    return identifier


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def positions(sequence, pattern):
    start = sequence.find(pattern)
    while start >= 0:
        yield start
        start = sequence.find(pattern, start + 1)


def edit_alignment(reference, observed, maximum):
    """Unit-cost global alignment; preserve a tie flag and one deterministic path.

    I consumes an observed base, D a reference base; counts therefore describe
    the observed read relative to the hypothesis. A tie never selects a class.
    """
    if abs(len(reference) - len(observed)) > maximum:
        return None
    limit = maximum + 1
    scores = {(0, 0): 0}
    ways = {(0, 0): 1}
    previous = {}
    for i in range(len(reference) + 1):
        for j in range(max(0, i - maximum), min(len(observed), i + maximum) + 1):
            if i == j == 0:
                continue
            choices = []
            if i and j:
                operation = "=" if reference[i - 1] == observed[j - 1] else "X"
                choices.append((scores.get((i - 1, j - 1), limit) + (operation == "X"), (i - 1, j - 1), operation))
            if j:
                choices.append((scores.get((i, j - 1), limit) + 1, (i, j - 1), "I"))
            if i:
                choices.append((scores.get((i - 1, j), limit) + 1, (i - 1, j), "D"))
            best = min(choice[0] for choice in choices)
            if best > maximum:
                continue
            tied = [choice for choice in choices if choice[0] == best]
            scores[i, j] = best
            ways[i, j] = min(2, sum(ways.get(choice[1], 0) for choice in tied))
            previous[i, j] = tied[0][1:]
    endpoint = len(reference), len(observed)
    if endpoint not in scores:
        return None
    operations = []
    cursor = endpoint
    while cursor != (0, 0):
        cursor, operation = previous[cursor]
        operations.append(operation)
    operations.reverse()
    cigar = []
    for operation in operations:
        if cigar and cigar[-1][1] == operation:
            cigar[-1][0] += 1
        else:
            cigar.append([1, operation])
    return {"distance": scores[endpoint], "substitutions": operations.count("X"),
            "insertions": operations.count("I"), "deletions": operations.count("D"),
            "cigar": "".join(f"{count}{operation}" for count, operation in cigar),
            "tied_optimal_alignments": ways[endpoint] > 1}


def edit_distance(reference, observed, maximum):
    """Allocation-light banded distance used for seed verification."""
    if abs(len(reference) - len(observed)) > maximum:
        return None
    limit = maximum + 1
    previous = list(range(len(observed) + 1))
    for i, base in enumerate(reference, 1):
        current = [limit] * (len(observed) + 1)
        if i <= maximum:
            current[0] = i
        for j in range(max(1, i - maximum), min(len(observed), i + maximum) + 1):
            current[j] = min(previous[j - 1] + (base != observed[j - 1]),
                             previous[j] + 1, current[j - 1] + 1)
        if min(current) > maximum:
            return None
        previous = current
    distance = previous[-1]
    return distance if distance <= maximum else None


def all_anchor_hits(sequence, pattern, maximum=FLANK_EDITS, hit_limit=None):
    """All bounded edit-distance substring hits, with exhaustive seed checking.

    maximum+1 disjoint blocks guarantee an intact block under maximum edits.
    A seed position implies a start within +/-maximum; check every legal length.
    """
    if len(pattern) < maximum + 1:
        raise ValueError("Pattern too short for the frozen seed scheme")
    starts = set()
    for block in range(maximum + 1):
        offset = block * len(pattern) // (maximum + 1)
        end = (block + 1) * len(pattern) // (maximum + 1)
        for hit in positions(sequence, pattern[offset:end]):
            starts.update(range(max(0, hit - offset - maximum), min(len(sequence), hit - offset + maximum) + 1))
    hits = []
    for start in sorted(starts):
        for size in range(max(1, len(pattern) - maximum), len(pattern) + maximum + 1):
            end = start + size
            if end > len(sequence):
                continue
            observed = sequence[start:end]
            if set(observed) - BASES:
                continue
            distance = edit_distance(pattern, observed, maximum)
            if distance is not None:
                hits.append({"start": start, "end": end, "distance": distance})
                if hit_limit is not None and len(hits) > hit_limit:
                    raise ValueError("Alignment placement cap exceeded; no successful partial result")
    return hits


def validated_sequence(item):
    sequence = item["sequence"].encode("ascii")
    if not sequence or set(sequence) - BASES or digest(sequence) != item["sequence_sha256"]:
        raise ValueError("Nonliteral sequence or mismatched sequence hash")
    return sequence


def prepare_manifest(manifest):
    if len(encoded_json(manifest)) > MAX_MANIFEST_BYTES:
        raise ValueError("JSON manifest byte cap exceeded")
    if manifest.get("schema_version") != 1 or not manifest.get("regions"):
        raise ValueError("Expected a version-1 nonempty region manifest")
    if len(manifest["regions"]) > MAX_REGIONS:
        raise ValueError("Region cap exceeded")
    regions = []
    seen_regions = set()
    class_count = 0
    for region in manifest["regions"]:
        valid_identifier(region["id"])
        if region["id"] in seen_regions or not region.get("classes"):
            raise ValueError("Duplicate region ID or empty classes")
        class_count += len(region["classes"])
        if len(region["classes"]) > MAX_CLASSES_PER_REGION or class_count > MAX_TOTAL_CLASSES:
            raise ValueError("Manifest class cap exceeded")
        seen_regions.add(region["id"])
        classes = []
        seen_classes, seen_sequences = set(), set()
        for item in region["classes"]:
            valid_identifier(item["id"])
            if len(item.get("self_context", [])) > MAX_SELF_CONTEXTS_PER_CLASS:
                raise ValueError("Self-context interval cap exceeded")
            sequence = validated_sequence(item)
            if len(sequence) < 2 * FLANK_BASES or len(sequence) > 150:
                raise ValueError("Class requires two complete flanks and <=150 bases")
            if item["id"] in seen_classes or sequence in seen_sequences:
                raise ValueError("Duplicate class identity or spelling in region")
            seen_classes.add(item["id"])
            seen_sequences.add(sequence)
            classes.append({**item, "_sequence": sequence})
        regions.append({**region, "classes": classes,
                        "_left": sorted({item["_sequence"][:FLANK_BASES] for item in classes}),
                        "_right": sorted({item["_sequence"][-FLANK_BASES:] for item in classes})})
    decoys = []
    seen_decoys = set()
    total_decoy_bases = 0
    if len(manifest.get("decoys", [])) > MAX_DECOYS:
        raise ValueError("Decoy cap exceeded")
    for decoy in manifest.get("decoys", []):
        valid_identifier(decoy["id"])
        if decoy["id"] in seen_decoys:
            raise ValueError("Duplicate decoy ID")
        seen_decoys.add(decoy["id"])
        sequence = validated_sequence(decoy)
        total_decoy_bases += len(sequence)
        if len(sequence) > MAX_DECOY_LENGTH or total_decoy_bases > MAX_TOTAL_DECOY_BASES:
            raise ValueError("Decoy base cap exceeded")
        decoys.append({**decoy, "_sequence": sequence})
    for region in regions:
        for item in region["classes"]:
            for context in item.get("self_context", []):
                if context["decoy_id"] not in seen_decoys or context["strand"] not in ("+", "-"):
                    raise ValueError("Invalid self-context identity/orientation")
                decoy = next(value for value in decoys if value["id"] == context["decoy_id"])
                start, end = context["start"], context["end"]
                if type(start) is not int or type(end) is not int or not 0 <= start < end <= len(decoy["_sequence"]):
                    raise ValueError("Invalid self-context interval")
                expected = decoy["_sequence"][start:end]
                if context["strand"] == "-":
                    expected = reverse_complement(expected)
                if expected != item["_sequence"]:
                    raise ValueError("Self-context does not spell the complete class")
    return regions, decoys


def physical_interval(start, end, length, strand):
    return [start, end] if strand == "+" else [length - end, length - start]


def decoy_placements(observed, decoys, item):
    placements = []
    for decoy in decoys:
        sequence = decoy["_sequence"]
        for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
            for hit in all_anchor_hits(oriented, observed, SPAN_EDITS, MAX_DECOY_PLACEMENTS):
                interval = physical_interval(hit["start"], hit["end"], len(sequence), strand)
                contexts = [context for context in item.get("self_context", [])
                            if context["decoy_id"] == decoy["id"] and context["strand"] == strand]
                placements.append({"decoy_id": decoy["id"], "strand": strand,
                                   "interval": interval, "distance": hit["distance"],
                                   "overlaps_self_context": any(interval[0] < context["end"] and context["start"] < interval[1] for context in contexts),
                                   "exact_self_interval": any(interval == [context["start"], context["end"]] for context in contexts)})
                if len(placements) > MAX_DECOY_PLACEMENTS:
                    raise ValueError("Decoy placement cap exceeded; no successful partial result")
    return placements


def observations_for_read(sequence, quality, regions, decoys):
    if not sequence or len(sequence) != len(quality) or len(sequence) > MAX_READ_LENGTH:
        raise ValueError("Invalid read length or read length cap exceeded")
    if min(quality) < 33 or max(quality) > 126:
        raise ValueError("Unsupported quality encoding")
    observations = []
    total_pairs = 0
    total_anchors = 0
    total_decoy_placements = 0
    orientations = [("+", sequence, quality), ("-", reverse_complement(sequence), quality[::-1])]
    cache = {}
    for region in regions:
        anchors, pairs = [], []
        for strand, oriented, oriented_quality in orientations:
            sides = {}
            for side in ("left", "right"):
                hits = {}
                for pattern in region["_" + side]:
                    key = strand, pattern
                    if key not in cache:
                        cache[key] = all_anchor_hits(oriented, pattern, hit_limit=MAX_FLANK_PLACEMENTS_PER_READ)
                    for hit in cache[key]:
                        identity = hit["start"], hit["end"]
                        record = hits.setdefault(identity, {"strand": strand, "side": side,
                            "start": hit["start"], "end": hit["end"],
                            "read_interval": physical_interval(hit["start"], hit["end"], len(sequence), strand),
                            "q20": min(oriented_quality[hit["start"]:hit["end"]]) >= MIN_QUALITY + 33,
                            "flank_matches": []})
                        record["flank_matches"].append({"sequence": pattern.decode("ascii"), "distance": hit["distance"]})
                sides[side] = list(hits.values())
                total_anchors += len(sides[side])
                if total_anchors > MAX_FLANK_PLACEMENTS_PER_READ:
                    raise ValueError("Per-read flank placement cap exceeded; no successful partial result")
                anchors.extend(sides[side])
            for left in sides["left"]:
                for right in sides["right"]:
                    if left["end"] > right["start"]:
                        continue
                    total_pairs += 1
                    if total_pairs > MAX_PAIRS_PER_READ:
                        raise ValueError("Per-read ordered flank-pair cap exceeded; no successful partial result")
                    span = oriented[left["start"]:right["end"]]
                    literal = not (set(span) - BASES)
                    matches = []
                    if literal:
                        for item in region["classes"]:
                            alignment = edit_alignment(item["_sequence"], span, SPAN_EDITS)
                            if alignment is not None:
                                matches.append({"class_id": item["id"], **alignment})
                    pairs.append({"strand": strand, "start": left["start"], "end": right["end"],
                                  "read_interval": physical_interval(left["start"], right["end"], len(sequence), strand),
                                  "left_interval": [left["start"], left["end"]],
                                  "right_interval": [right["start"], right["end"]],
                                  "observed_span": span.decode("ascii"), "span_length": len(span),
                                  "inter_anchor_bases": right["start"] - left["end"],
                                  "q20": min(oriented_quality[left["start"]:right["end"]]) >= MIN_QUALITY + 33,
                                  "literal_acgt": literal, "class_matches": matches})
        if not anchors:
            continue
        # Deduplicate the same class/physical span reached through alternative
        # approximate flank endpoints, while preserving every original pair.
        alternatives = {}
        qualified_alternatives = {}
        for pair_index, pair in enumerate(pairs):
            for match in pair["class_matches"]:
                identity = match["class_id"], pair["strand"], tuple(pair["read_interval"])
                contender = {**match, "strand": pair["strand"], "read_interval": pair["read_interval"],
                             "oriented_interval": [pair["start"], pair["end"]], "pair_indices": []}
                for mapping in ([alternatives, qualified_alternatives] if pair["q20"] else [alternatives]):
                    stored = mapping.setdefault(identity, dict(contender, pair_indices=[]))
                    stored["pair_indices"].append(pair_index)
        def select(mapping):
            if not mapping:
                return []
            best = min(value["distance"] for value in mapping.values())
            return [value for value in mapping.values() if value["distance"] == best]
        winners = select(alternatives)
        qualified = select(qualified_alternatives)
        # Decoy checks are descriptive scope controls for every minimal contender,
        # including a qualified minimum if the unfiltered minimum is low quality.
        decoy_cache = {}
        for contender in winners + qualified:
            pair = pairs[contender["pair_indices"][0]]
            key = contender["class_id"], pair["observed_span"]
            if key not in decoy_cache:
                item = next(item for item in region["classes"] if item["id"] == contender["class_id"])
                decoy_cache[key] = decoy_placements(pair["observed_span"].encode("ascii"), decoys, item)
                total_decoy_placements += len(decoy_cache[key])
                if total_decoy_placements > MAX_DECOY_PLACEMENTS_PER_READ:
                    raise ValueError("Per-read decoy placement cap exceeded; no successful partial result")
            contender["decoy_placements"] = decoy_cache[key]
        status = "qualified_unique_local_alignment" if len(qualified) == 1 else "qualified_ambiguous_local_alignment" if qualified else "no_qualified_local_alignment"
        present = {anchor["side"] for anchor in anchors}
        recruitment = "both_flanks" if len(present) == 2 else next(iter(present)) + "_only"
        observations.append({"region_id": region["id"], "anchors": anchors, "pairs": pairs,
                             "recruitment": recruitment, "no_ordered_pair": not bool(pairs),
                             "best_unfiltered": winners, "best_q20": qualified, "status": status})
        if len(encoded_json(observations)) > MAX_RECORD_EVIDENCE_BYTES:
            raise ValueError("Per-record serialized evidence cap exceeded; no successful partial result")
    return observations


def scan_fastq(path, expected, manifest, staged_before=None):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Input must be a regular FASTQ file")
    records = expected.get("records")
    if type(records) is not int or records <= 0:
        raise ValueError("An explicit positive prefix record count is required")
    if any(type(expected.get(key)) is not int or expected[key] <= 0 for key in ("bases", "size_bytes")):
        raise ValueError("Positive expected prefix base and byte counts are required")
    if not isinstance(expected.get("sha256"), str) or len(expected["sha256"]) != 64 or any(base not in "0123456789abcdef" for base in expected["sha256"]):
        raise ValueError("Expected prefix identity requires a lowercase SHA-256")
    regions, decoys = prepare_manifest(manifest)
    accumulators = {region["id"]: {"recruited_records": 0, "left_records": 0, "right_records": 0,
        "left_q20_records": 0, "right_q20_records": 0, "both_flanks_records": 0,
        "left_only_records": 0, "right_only_records": 0, "no_ordered_pair_records": 0,
        "both_flanks_no_ordered_pair_records": 0, "ids": set(),
        "ordered_pair_records": 0, "q20_pair_records": 0, "non_acgt_pair_records": 0,
        "literal_pair_records": 0, "non_acgt_only_pair_records": 0,
        "low_quality_pair_records": 0, "low_quality_only_pair_records": 0,
        "low_quality_alignment_records": 0, "low_quality_only_alignment_records": 0,
        "outside_distance_pair_records": 0, "unaligned_pair_records": 0,
        "unique_records": 0, "ambiguous_records": 0,
        "classes": {item["id"]: {"unique_records": 0, "exact_records": 0,
            "substitution_only_records": 0, "indel_records": 0,
            "tied_traceback_records": 0, "ids": set(), "offsets": set(), "strands": set(), "spans": set()}
            for item in region["classes"]}} for region in regions}
    before = file_identity(path) if staged_before is None else staged_before
    prefix = hashlib.sha256()
    size_bytes = bases = maximum_read_length = 0
    ledger = []
    ledger_bytes = 0
    with path.open("rb") as stream:
        for ordinal in range(1, records + 1):
            lines = [stream.readline(MAX_FASTQ_LINE_BYTES + 1) for _ in range(4)]
            if any(len(line) > MAX_FASTQ_LINE_BYTES for line in lines):
                raise ValueError(f"FASTQ line byte cap exceeded at record {ordinal}")
            if any(not line for line in lines):
                raise ValueError(f"Truncated FASTQ at record {ordinal}")
            if any(not line.endswith(b"\n") for line in lines):
                raise ValueError(f"Unterminated FASTQ line at record {ordinal}")
            header, sequence, separator, quality = [line[:-2] if line.endswith(b"\r\n") else line[:-1] for line in lines]
            if not header.startswith(b"@") or not header[1:].strip() or not separator.startswith(b"+") or not sequence or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record {ordinal}")
            if len(separator) > 1 and separator[1:] != header[1:]:
                raise ValueError(f"FASTQ separator/header mismatch at record {ordinal}")
            if any(value < 33 or value > 126 for value in sequence) or any(value < 32 or value > 126 for value in header):
                raise ValueError(f"Nonprintable or non-ASCII FASTQ at record {ordinal}")
            if len(sequence) > MAX_READ_LENGTH:
                raise ValueError(f"Read length cap exceeded at record {ordinal}")
            if min(quality) < 33 or max(quality) > 126:
                raise ValueError(f"Unsupported quality encoding at record {ordinal}")
            prefix.update(b"".join(lines))
            size_bytes += sum(map(len, lines))
            bases += len(sequence)
            maximum_read_length = max(maximum_read_length, len(sequence))
            sequence = sequence.upper()
            observations = observations_for_read(sequence, quality, regions, decoys)
            if not observations:
                continue
            if len(ledger) >= MAX_RECRUITED_RECORDS:
                raise ValueError("Recruited-record cap exceeded; no successful partial result")
            read_id = header.split()[0].decode("ascii")
            record = {"ordinal": ordinal, "header": header.decode("ascii"), "read_id": read_id,
                           "sequence": sequence.decode("ascii"), "quality": quality.decode("ascii"),
                           "observations": observations}
            record_bytes = len(encoded_json(record))
            if record_bytes > MAX_RECORD_EVIDENCE_BYTES:
                raise ValueError("Per-record serialized evidence cap exceeded; no successful partial result")
            ledger_bytes += record_bytes
            if ledger_bytes > MAX_OUTPUT_BYTES:
                raise ValueError("Serialized output byte cap exceeded; no successful partial result")
            ledger.append(record)
            for observation in observations:
                accumulator = accumulators[observation["region_id"]]
                accumulator["recruited_records"] += 1
                accumulator["ids"].add(read_id)
                present = set()
                for side in ("left", "right"):
                    anchors = [anchor for anchor in observation["anchors"] if anchor["side"] == side]
                    accumulator[side + "_records"] += bool(anchors)
                    accumulator[side + "_q20_records"] += any(anchor["q20"] for anchor in anchors)
                    if anchors:
                        present.add(side)
                accumulator["both_flanks_records"] += len(present) == 2
                for side in ("left", "right"):
                    accumulator[side + "_only_records"] += present == {side}
                pairs = observation["pairs"]
                literal_pairs = [pair for pair in pairs if pair["literal_acgt"]]
                accumulator["no_ordered_pair_records"] += not pairs
                accumulator["both_flanks_no_ordered_pair_records"] += len(present) == 2 and not pairs
                accumulator["ordered_pair_records"] += bool(pairs)
                accumulator["literal_pair_records"] += bool(literal_pairs)
                accumulator["q20_pair_records"] += any(pair["q20"] for pair in literal_pairs)
                accumulator["non_acgt_pair_records"] += any(not pair["literal_acgt"] for pair in pairs)
                accumulator["non_acgt_only_pair_records"] += bool(pairs) and not literal_pairs
                accumulator["low_quality_pair_records"] += any(not pair["q20"] for pair in literal_pairs)
                accumulator["low_quality_only_pair_records"] += bool(literal_pairs) and not any(pair["q20"] for pair in literal_pairs)
                low_quality_alignment = any(pair["class_matches"] and not pair["q20"] for pair in literal_pairs)
                accumulator["low_quality_alignment_records"] += low_quality_alignment
                accumulator["low_quality_only_alignment_records"] += low_quality_alignment and not observation["best_q20"]
                accumulator["outside_distance_pair_records"] += any(not pair["class_matches"] for pair in literal_pairs)
                accumulator["unaligned_pair_records"] += bool(literal_pairs) and not observation["best_unfiltered"]
                winners = observation["best_q20"]
                if len(winners) > 1:
                    accumulator["ambiguous_records"] += 1
                elif winners:
                    accumulator["unique_records"] += 1
                    winner = winners[0]
                    result = accumulator["classes"][winner["class_id"]]
                    result["unique_records"] += 1
                    result["exact_records"] += winner["distance"] == 0
                    result["substitution_only_records"] += winner["distance"] > 0 and winner["insertions"] + winner["deletions"] == 0
                    result["indel_records"] += winner["insertions"] + winner["deletions"] > 0
                    result["tied_traceback_records"] += winner["tied_optimal_alignments"]
                    result["ids"].add(read_id)
                    result["offsets"].add(tuple(winner["oriented_interval"]))
                    result["strands"].add(winner["strand"])
                    result["spans"].add(pairs[winner["pair_indices"][0]]["observed_span"])
    observed = {"records": records, "bases": bases, "size_bytes": size_bytes, "sha256": prefix.hexdigest()}
    if any(observed[key] != expected[key] for key in observed):
        raise ValueError("Consumed-prefix identity differs; no read-support result published")
    after = file_identity(path)
    if before != after:
        raise ValueError("Input changed during audit")
    summaries = []
    for region in regions:
        summary = {"region_id": region["id"], **accumulators[region["id"]]}
        summary["distinct_recruited_read_ids"] = len(summary.pop("ids"))
        for result in summary["classes"].values():
            result["distinct_read_ids"] = len(result.pop("ids"))
            result["distinct_oriented_span_offsets"] = len(result["offsets"])
            result["oriented_span_offsets"] = sorted(result.pop("offsets"))
            result["strands"] = sorted(result["strands"])
            result["distinct_observed_spans"] = len(result.pop("spans"))
        summary["status"] = ("qualified_local_alignments" if summary["unique_records"] else
                             "ambiguous_local_alignments" if summary["ambiguous_records"] else
                             "recruited_without_qualified_alignment" if summary["recruited_records"] else
                             "no_recruited_flanks")
        summaries.append(summary)
    result = {"schema_version": 1, "input_subset": observed,
            "staged_file_before": before, "staged_file_after": after,
            "maximum_read_length": maximum_read_length,
            "parameters": {"flank_bases": FLANK_BASES, "flank_edits": FLANK_EDITS,
                           "span_edits": SPAN_EDITS, "minimum_quality": MIN_QUALITY,
                           "max_pairs_per_read": MAX_PAIRS_PER_READ, "max_recruited_records": MAX_RECRUITED_RECORDS,
                           **{name.lower(): value for name, value in globals().items() if name.startswith("MAX_")}},
            "limitations": ["Finite hypothesis/decoy universe, not genome-wide uniqueness",
                            "Gapped local alignment does not establish repeat count or full-product phase",
                            "Q20 applies to observed span bases; deleted hypothesis bases have no observed quality",
                            "Distinct R1 IDs and oriented span offsets are not UMI-certified molecules or exact gapped genomic footprints",
                            "Edit-type summaries describe the deterministic representative optimal traceback; tied_traceback flag preserves uncertainty",
                            "Per-region record categories overlap; q20_pair_records and low-quality pair categories require literal ACGT spans"],
            "regions": summaries, "read_evidence": ledger}
    if len(encoded_json(result)) > MAX_OUTPUT_BYTES:
        raise ValueError("Serialized output byte cap exceeded; no successful partial result")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    failure_path = arguments.output.with_name(arguments.output.name + ".failure.json")
    for path in (arguments.output, failure_path):
        if path.exists() or path.is_symlink():
            raise ValueError("Refusing to overwrite an audit output or failure receipt")
    started = time.monotonic()
    provenance = {"argv": sys.argv, "fastq": str(arguments.fastq.resolve()),
                  "manifest": str(arguments.manifest.resolve()),
                  "expected_prefix": str(arguments.expected_prefix.resolve())}
    try:
        manifest_bytes = bounded_json(arguments.manifest)
        prefix_bytes = bounded_json(arguments.expected_prefix)
        source_bytes = Path(__file__).read_bytes()
        provenance.update({"manifest_sha256": digest(manifest_bytes), "expected_prefix_manifest_sha256": digest(prefix_bytes),
                           "scanner_sha256": digest(source_bytes)})
        staged_before = file_identity(arguments.fastq)
        provenance["staged_file_before"] = staged_before
        result = scan_fastq(arguments.fastq, json.loads(prefix_bytes), json.loads(manifest_bytes), staged_before)
        if bounded_json(arguments.manifest) != manifest_bytes or bounded_json(arguments.expected_prefix) != prefix_bytes or Path(__file__).read_bytes() != source_bytes:
            raise ValueError("Manifest or source changed during evaluation")
        result["provenance"] = provenance
        result["status"] = "completed"
        result["elapsed_seconds"] = time.monotonic() - started
        publish_json(arguments.output, result)
    except Exception as error:
        publish_json(failure_path, {"status": "failed", "error_type": type(error).__name__, "error": str(error),
                                    "elapsed_seconds": time.monotonic() - started,
                                    "provenance": provenance}, maximum=MAX_MANIFEST_BYTES)
        raise
    print(json.dumps({"output": str(arguments.output), "recruited_records": len(result["read_evidence"])}))


if __name__ == "__main__":
    main()

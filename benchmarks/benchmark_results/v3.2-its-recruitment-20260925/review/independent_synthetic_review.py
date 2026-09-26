#!/usr/bin/env python3
"""Independent synthetic oracle; no biological FASTQ input is read.

The oracle uses an unbanded rectangular DP rather than production seed/band
logic. All reviewed real sequence strings come from the predeclared hypothesis
manifest, never from read outcomes. Run from the repository root.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import time

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "scripts/audit_its_recruitment.py"
MANIFEST = Path(__file__).resolve().parents[1] / "inputs/manifest.json"
spec = importlib.util.spec_from_file_location("reviewed_recruitment", SOURCE)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def oracle(reference, observed):
    score = [[0] * (len(observed) + 1) for _ in range(len(reference) + 1)]
    paths = [[1] * (len(observed) + 1) for _ in range(len(reference) + 1)]
    for i in range(len(reference) + 1):
        score[i][0] = i
    for j in range(len(observed) + 1):
        score[0][j] = j
    for i in range(1, len(reference) + 1):
        for j in range(1, len(observed) + 1):
            candidates = [(score[i-1][j-1] + (reference[i-1] != observed[j-1]), paths[i-1][j-1]),
                          (score[i-1][j] + 1, paths[i-1][j]),
                          (score[i][j-1] + 1, paths[i][j-1])]
            best = min(value[0] for value in candidates)
            score[i][j] = best
            paths[i][j] = min(2, sum(value[1] for value in candidates if value[0] == best))
    return score[-1][-1], paths[-1][-1] > 1


def exhaustive(sequence, pattern, maximum):
    result = []
    for start in range(len(sequence)):
        for size in range(max(1, len(pattern) - maximum), len(pattern) + maximum + 1):
            end = start + size
            observed = sequence[start:end]
            if end > len(sequence) or any(base not in b"ACGT" for base in observed):
                continue
            distance, _ = oracle(pattern, observed)
            if distance <= maximum:
                result.append({"start": start, "end": end, "distance": distance})
    return result


def rc(sequence):
    return bytes({65:84, 84:65, 67:71, 71:67}[base] for base in reversed(sequence))


def mutate(sequence, positions):
    result = bytearray(sequence)
    for index in positions:
        result[index] = b"ACGT"[(b"ACGT".index(result[index]) + 1) % 4]
    return bytes(result)


def main():
    started = time.monotonic()
    source_start = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    rng = random.Random(20260925)
    counts = {"distance_traceback_comparisons": 0, "exhaustive_seed_comparisons": 0}
    for number in range(450):
        n, m = rng.randrange(0, 25), rng.randrange(0, 25)
        reference = bytes(rng.choice(b"ACGT") for _ in range(n))
        if number % 2:
            observed = reference
            for _ in range(rng.randrange(6)):
                if not observed:
                    observed = b"A"
                where = rng.randrange(len(observed))
                operation = rng.randrange(3)
                observed = (mutate(observed, [where]) if operation == 0 else
                            observed[:where] + b"G" + observed[where:] if operation == 1 else
                            observed[:where] + observed[where+1:])
        else:
            observed = bytes(rng.choice(b"ACGT") for _ in range(m))
        expected_distance, expected_tie = oracle(reference, observed)
        for maximum in (0, 1, 2, 4):
            expected = expected_distance if expected_distance <= maximum else None
            assert audit.edit_distance(reference, observed, maximum) == expected
            alignment = audit.edit_alignment(reference, observed, maximum)
            assert (alignment["distance"] if alignment else None) == expected
            if alignment:
                assert alignment["tied_optimal_alignments"] == expected_tie
                assert sum(alignment[key] for key in ("substitutions", "insertions", "deletions")) == expected_distance
            counts["distance_traceback_comparisons"] += 1
    for maximum, length, repeats in ((2, 21, 70), (4, 25, 30)):
        for number in range(repeats):
            pattern = (b"A" * length if number == 0 else
                       (b"ACGT" * length)[:length] if number == 1 else
                       bytes(rng.choice(b"ACGT") for _ in range(length)))
            altered = pattern
            for _ in range(number % (maximum + 2)):
                where = rng.randrange(len(altered))
                op = rng.randrange(3)
                altered = (mutate(altered, [where]) if op == 0 else
                           altered[:where] + b"T" + altered[where:] if op == 1 else
                           altered[:where] + altered[where+1:])
            sequence = b"GCA" + altered + b"TGC"
            if number % 11 == 10:
                sequence = sequence[:7] + b"N" + sequence[8:]
            assert audit.all_anchor_hits(sequence, pattern, maximum) == exhaustive(sequence, pattern, maximum)
            counts["exhaustive_seed_comparisons"] += 1
    regions, decoys = audit.prepare_manifest(manifest)
    fixtures = []
    assert audit.MAX_PAIRS_PER_READ == 4096, "Review must retain preregistered synthetic-supported pair cap"
    for region in manifest["regions"]:
        for item in region["classes"]:
            sequence = item["sequence"].encode()
            variants = {"exact": sequence, "reverse_complement": rc(sequence),
                        "two_substitutions": mutate(sequence, [24, 28]),
                        "four_substitutions": mutate(sequence, [23, 26, 30, 34]),
                        "flank_substitutions_limit": mutate(sequence, [4, 13]),
                        "flank_insertion": sequence[:9] + b"A" + sequence[9:],
                        "flank_deletion": sequence[:9] + sequence[10:],
                        "interior_insertion": sequence[:29] + b"AC" + sequence[29:],
                        "interior_deletion": sequence[:29] + sequence[31:]}
            for label, read in variants.items():
                # Class variants exercise real frozen flank complexity. Only exact
                # and reverse-complement fixtures invoke the full decoy collection.
                selected_decoys = decoys if label in ("exact", "reverse_complement") else []
                observations = audit.observations_for_read(read, b"I" * len(read), regions, selected_decoys)
                own = next(observation for observation in observations if observation["region_id"] == region["id"])
                expected_sequence = rc(read) if label == "reverse_complement" else read
                expected_distance = oracle(sequence, expected_sequence)[0]
                best = own["best_q20"]
                assert best, (item["id"], label, "no qualified alignment")
                assert any(winner["class_id"] == item["id"] for winner in best), (item["id"], label, "class lost")
                assert min(winner["distance"] for winner in best) <= expected_distance
                if selected_decoys:
                    expected_exact_contexts = set()
                    for decoy in manifest["decoys"]:
                        raw = decoy["sequence"].encode()
                        for orientation, search_string in (("+", raw), ("-", rc(raw))):
                            for position in range(len(raw)-len(sequence)+1):
                                if search_string[position:position+len(sequence)] == sequence:
                                    interval = ((position, position+len(sequence)) if orientation == "+" else
                                                (len(raw)-position-len(sequence), len(raw)-position))
                                    expected_exact_contexts.add((decoy["id"], orientation, interval))
                    actual_exact_contexts = {(value["decoy_id"], value["strand"], tuple(value["interval"]))
                                             for value in best[0]["decoy_placements"] if value["distance"] == 0}
                    assert actual_exact_contexts == expected_exact_contexts
                    for context in item.get("self_context", []):
                        assert any(value["exact_self_interval"] and value["distance"] == 0 and
                                   value["decoy_id"] == context["decoy_id"] and value["strand"] == context["strand"]
                                   for value in best[0]["decoy_placements"])
                pair_count = sum(len(observation["pairs"]) for observation in observations)
                assert pair_count <= 256, "Individual frozen variants also fit the original proposed bound"
                fixtures.append({"class_id": item["id"], "variant": label, "read_length": len(read),
                                 "total_pairs": pair_count, "status": own["status"],
                                 "best_distance": best[0]["distance"], "oracle_whole_read_distance": expected_distance,
                                 "best_contenders": len(best),
                                 "decoy_placements": sum(len(winner.get("decoy_placements", [])) for winner in best)})
    repeated_fixtures = []
    for region in manifest["regions"]:
        item = region["classes"][0]
        sequence = item["sequence"].encode()
        read = sequence + b"ACGTA" + sequence
        assert len(read) <= 150
        observations = audit.observations_for_read(read, b"I" * len(read), regions, [])
        own = next(observation for observation in observations if observation["region_id"] == region["id"])
        best = own["best_q20"]
        assert own["status"] == "qualified_ambiguous_local_alignment"
        assert len(best) == 2 and all(winner["distance"] == 0 for winner in best)
        assert {tuple(winner["read_interval"]) for winner in best} == {(0, len(sequence)), (len(sequence)+5, len(read))}
        pairs = sum(len(observation["pairs"]) for observation in observations)
        assert pairs > 256 and pairs <= 4096
        repeated_fixtures.append({"class_id": item["id"], "read_length": len(read),
                                  "total_pairs": pairs, "best_exact_placements": len(best),
                                  "status": own["status"]})
    source_end = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    result = {"status": "passed", "biological_reads_examined": 0,
              "scanner_sha256_start": source_start, "scanner_sha256_end": source_end,
              "source_unchanged_during_review": source_start == source_end,
              "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
              "oracle_seed": 20260925, **counts, "frozen_class_fixtures": fixtures,
              "repeated_placement_fixtures": repeated_fixtures,
              "pair_cap": audit.MAX_PAIRS_PER_READ,
              "maximum_frozen_fixture_pairs": max(value["total_pairs"] for value in fixtures),
              "elapsed_seconds": time.monotonic() - started}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

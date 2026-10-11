#!/usr/bin/env python3
"""Bounded, exact linkage audit of saved positive reads; no FASTQ rescan.

Coordinates are zero-based half-open. A, C, B name literal blocks, not
uniquely located biological mutations. Every witness is one contiguous read.
"""
import argparse
import csv
import hashlib
import itertools
import json
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
PILOT = BASE.parent / "v3.2-its-gain-pilot-20261009"
LEFT = "GCTTGCTCCTATACGCTCGCG"
MID_A = "TCCGGGCGGGCCGGGTGGAAATGACCGTCGCCGGCGCGTTCT"
MID_B = "GCGCGTCGGCGTCGTA"
RIGHT = "TTTTCGGCTAACTCGTCTTTT"
MIN_Q = 20
RC = str.maketrans("ACGTN", "TGCAN")


def sha(value):
    return hashlib.sha256(value).hexdigest()


def receipt(path):
    return {"path": str(path), "sha256": sha(path.read_bytes()), "bytes": path.stat().st_size}


def rc(sequence):
    return sequence.translate(RC)[::-1]


def positions(sequence, word):
    start = 0
    while True:
        start = sequence.find(word, start)
        if start < 0:
            return
        yield start
        start += 1


def orientations(read):
    yield "+", read["sequence"], read["quality"]
    yield "-", rc(read["sequence"]), read["quality"][::-1]


def witness(read, strand, start, end, quality):
    n = len(read["sequence"])
    return {"ordinal": read["ordinal"], "read_id": read["header"].split()[0],
            "strand": strand, "oriented_read_span": [start, end],
            "original_read_span": [start, end] if strand == "+" else [n-end, n-start],
            "minimum_phred": min(ord(q)-33 for q in quality[start:end])}


def exact_matches(reads, sequence):
    result = []
    for read in reads:
        for strand, oriented, quality in orientations(read):
            for start in positions(oriented, sequence):
                item = witness(read, strand, start, start+len(sequence), quality)
                if item["minimum_phred"] >= MIN_Q:
                    result.append(item)
    return result


def anchored_spans(reads, left=LEFT, right=RIGHT):
    """Keep every ordered placement, both strands, including low-Q spans.

    Sequence between anchors is unconstrained: candidate-incompatible strings
    remain visible. The anchors themselves do not certify target specificity.
    """
    result = []
    for read in reads:
        for strand, sequence, quality in orientations(read):
            for start in positions(sequence, left):
                for right_start in positions(sequence, right):
                    if right_start >= start + len(left):
                        end = right_start + len(right)
                        item = witness(read, strand, start, end, quality)
                        item["sequence"] = sequence[start:end]
                        item["passes_q20_entire_span"] = item["minimum_phred"] >= MIN_Q
                        result.append(item)
    return result


def parse_blocks(sequence):
    expression = (re.escape(LEFT) + "([ACGT]+?)" + re.escape(MID_A)
                  + "([ACGT]+?)" + re.escape(MID_B) + "([ACGT]+?)" + re.escape(RIGHT))
    match = re.fullmatch(expression, sequence)
    return dict(zip(("A", "C", "B"), match.groups())) if match else None


def make_span(a, c, b):
    return LEFT + a + MID_A + c + MID_B + b + RIGHT


def segment_sequences(a, c, b):
    return {"A": LEFT+a+MID_A[:21],
            "C": MID_A[-21:]+c+MID_B,
            "B": MID_B+b+RIGHT,
            "A_C": LEFT+a+MID_A+c+MID_B,
            "C_B": MID_A[-21:]+c+MID_B+b+RIGHT,
            "A_C_B": make_span(a, c, b),
            "A_C_B_variable_core": a+MID_A+c+MID_B+b}


def core_matches(reads, a, c, b):
    """Partial-flank co-occurrence is distinct from the full flanked assay."""
    matches=exact_matches(reads,a+MID_A+c+MID_B+b)
    lookup={r["ordinal"]:r for r in reads}
    for match in matches:
        read=lookup[match["ordinal"]]
        sequence=read["sequence"] if match["strand"]=="+" else rc(read["sequence"])
        quality=read["quality"] if match["strand"]=="+" else read["quality"][::-1]
        start,end=match["oriented_read_span"]
        left=right=0
        while left<len(LEFT) and start-left>0 and sequence[start-left-1]==LEFT[-left-1] and ord(quality[start-left-1])-33>=MIN_Q:
            left+=1
        while right<len(RIGHT) and end+right<len(sequence) and sequence[end+right]==RIGHT[right] and ord(quality[end+right])-33>=MIN_Q:
            right+=1
        match["consecutive_exact_q20_left_flank_bases"]=left
        match["consecutive_exact_q20_right_flank_bases"]=right
        match["variable_core_length"]=end-start
    return matches


def ordinal_ids(matches):
    return sorted({m["ordinal"] for m in matches})


def candidate_placements(read, candidate, windows):
    """All offsets from exact Q20 event seeds, with no forced indel alignment.

    A seed is evidence of a local placement. A conflicting base outside that
    seed remains a conflict; a maximal compatible interval never jumps it.
    Competing offsets and strands are retained instead of picking a best one.
    """
    result = []
    template = candidate["sequence"]
    for strand, sequence, quality in orientations(read):
        offsets = {}
        for window in windows:
            for rpos in positions(sequence, window["sequence"]):
                rend = rpos+len(window["sequence"])
                if min(ord(q)-33 for q in quality[rpos:rend]) < MIN_Q:
                    continue
                for cpos in positions(template, window["sequence"]):
                    offsets.setdefault(cpos-rpos, []).append((rpos, rend, window["event_id"]))
        for offset, seeds in sorted(offsets.items()):
            rstart, rend = max(0, -offset), min(len(sequence), len(template)-offset)
            mismatches = [{"oriented_read_position": p, "candidate_position": p+offset,
                           "observed": sequence[p], "candidate": template[p+offset],
                           "phred": ord(quality[p])-33}
                          for p in range(rstart, rend) if sequence[p] != template[p+offset]]
            intervals = set()
            for left, right, _ in seeds:
                while left > rstart and sequence[left-1] == template[left-1+offset] and ord(quality[left-1])-33 >= MIN_Q:
                    left -= 1
                while right < rend and sequence[right] == template[right+offset] and ord(quality[right])-33 >= MIN_Q:
                    right += 1
                intervals.add((left, right))
            result.append({"strand": strand, "candidate_start_offset": offset,
                           "candidate_overlap": [rstart+offset, rend+offset],
                           "oriented_read_overlap": [rstart, rend],
                           "exact_overlap": not mismatches,
                           "full_read_within_candidate": rstart == 0 and rend == len(sequence),
                           "overlap_minimum_phred": min(ord(q)-33 for q in quality[rstart:rend]),
                           "mismatches": mismatches,
                           "seed_event_ids": sorted({seed[2] for seed in seeds}),
                           "contiguous_q20_exact_intervals": [
                               {"oriented_read_span": [a,b], "candidate_span": [a+offset,b+offset]}
                               for a,b in sorted(intervals)]})
    return result


def verify_sources():
    execution = json.loads((PILOT/"run-v1/execution.json").read_text())
    if execution["status"] != "pass" or not execution["sources_unchanged"]:
        raise ValueError("Unattested pilot execution")
    for name, expected in execution["outputs"].items():
        actual = receipt(PILOT/"run-v1"/name)
        if any(actual[k] != expected[k] for k in ("sha256", "bytes")):
            raise ValueError(f"Changed archived output: {name}")
    for name, expected in execution["source_hashes_after"].items():
        if sha(Path(name).read_bytes()) != expected:
            raise ValueError(f"Changed pilot input source: {name}")
    derivation = json.loads((PILOT/"derivation.json").read_text())
    for name, expected in derivation["result_files"].items():
        if sha((PILOT/name).read_bytes()) != expected["sha256"]:
            raise ValueError(f"Changed derived pilot input: {name}")
    return derivation


def analyze():
    derivation = verify_sources()
    raw = json.loads((PILOT/"run-v1/raw.json").read_text())
    manifest = json.loads((PILOT/"manifest.json").read_text())
    all_reads = raw["read_evidence"]
    reads = [r for r in all_reads if any(o["event_id"].startswith("its_local_") for o in r["observations"])]
    if len(reads) != 25 or raw["minimum_quality"] != MIN_Q or raw["max_substitutions"] != 0:
        raise ValueError("Pilot input ledger or assay changed")
    # Recompute all eight ITS and both control windows, preserving the 27/0 control.
    controls = []
    for event in raw["events"]:
        ids = ordinal_ids(exact_matches(all_reads, event["window_sequence"]))
        if len(ids) != event["records"]:
            raise ValueError(f"Archived event count failed: {event['id']}")
        if not event["id"].startswith("its_local_"):
            controls.append({"event_id": event["id"], "records": len(ids), "ordinals": ids})
    if [c["records"] for c in controls] != [27,0]:
        raise ValueError("Positive and historical longer-string controls changed")
    context = [c for c in manifest["candidates"] if c["id"].startswith("q") or
               any(s.get("gene") == "ITS_2" for s in c.get("source_candidates", []))]
    candidates = []
    for candidate in context:
        spans = [(start, end+len(RIGHT)) for start in positions(candidate["sequence"], LEFT)
                 for end in positions(candidate["sequence"], RIGHT) if end >= start+len(LEFT)]
        if len(spans) != 1:
            raise ValueError("Candidate block framing is ambiguous; explicit review required")
        start, end = spans[0]
        blocks = parse_blocks(candidate["sequence"][start:end])
        if blocks is None or sha(candidate["sequence"].encode()) != candidate["sequence_sha256"]:
            raise ValueError("Candidate frame or sequence digest invalid")
        candidates.append({**candidate, "blocks": blocks, "span": [start,end]})
    segment_cache = {}
    per_candidate = []
    for candidate in candidates:
        blocks = candidate["blocks"]
        segments = segment_sequences(blocks["A"], blocks["C"], blocks["B"])
        support = {}
        for name, sequence in segments.items():
            if sequence not in segment_cache:
                segment_cache[sequence] = exact_matches(reads, sequence)
            support[name] = ordinal_ids(segment_cache[sequence])
        per_candidate.append({"candidate_id": candidate["id"],
                              "role": "gained_uncertain" if candidate["id"].startswith("q") else "historical_ITS",
                              "sequence_sha256": candidate["sequence_sha256"],
                              "length": len(candidate["sequence"]), "blocks": blocks,
                              "A_C_B_candidate_span": candidate["span"],
                              "A_C_B_sequence": segments["A_C_B"],
                              "segment_read_ordinals": support,
                              "disposition": ("full_flanked_local_A_C_B_string_observed" if support["A_C_B"] else
                                              "variable_core_observed_with_partial_flanks" if support["A_C_B_variable_core"] else
                                              "variable_core_unobserved_in_saved_ledger"),
                              "full_candidate_phase": "unresolved"})
    span_observations = anchored_spans(reads)
    groups = []
    for sequence in sorted({m["sequence"] for m in span_observations}):
        matches = [m for m in span_observations if m["sequence"] == sequence]
        highq = [m for m in matches if m["passes_q20_entire_span"]]
        ids = [c["candidate_id"] for c in per_candidate if c["A_C_B_sequence"] == sequence]
        groups.append({"sequence": sequence, "length": len(sequence), "sequence_sha256": sha(sequence.encode()),
                       "blocks": parse_blocks(sequence), "candidate_ids": ids,
                       "q20_read_ordinals": ordinal_ids(highq), "all_placements": matches,
                       "in_saved_candidate_context": bool(ids)})
    # Enumerate finite combinations of candidate block strings explicitly as
    # analysis constructs. These are not new emitted products or inferred alleles.
    combinations = []
    for a,c,b in itertools.product(*(sorted({x["blocks"][k] for x in candidates}) for k in ("A","C","B"))):
        sequence = make_span(a,c,b)
        combinations.append({"blocks": {"A":a,"C":c,"B":b}, "sequence": sequence,
                             "candidate_ids": [x["candidate_id"] for x in per_candidate if x["A_C_B_sequence"] == sequence],
                             "q20_read_ordinals": ordinal_ids(exact_matches(reads, sequence)),
                             "variable_core_placements":core_matches(reads,a,c,b)})
    per_read = []
    for read in reads:
        placements = {c["id"]: candidate_placements(read,c,derivation["windows"]) for c in candidates}
        placements = {k:v for k,v in placements.items() if v}
        per_read.append({"ordinal": read["ordinal"], "header": read["header"],
                         "sequence": read["sequence"], "quality": read["quality"],
                         "sequence_sha256": sha(read["sequence"].encode()),
                         "original_event_ids": [o["event_id"] for o in read["observations"]],
                         "candidate_placements": placements,
                         "complete_A_C_B_placements": [m for m in span_observations if m["ordinal"] == read["ordinal"]]})
    gained = [c for c in per_candidate if c["role"] == "gained_uncertain"]
    return {"schema_version":1, "status":"completed_saved_read_linkage_audit", "minimum_quality_entire_span":MIN_Q,
            "coordinate_system":"zero_based_half_open", "retained_ITS_reads":len(reads),
            "candidate_context_count":len(candidates), "candidate_gains":len(gained),
            "anchors":{"left":LEFT,"between_A_C":MID_A,"between_C_B":MID_B,"right":RIGHT},
            "blocks_description": {"A":"overlapping original events 01–06, grouped into three literal strings",
                                   "C":"intervening variable string; required when evaluating A-to-B linkage",
                                   "B":"original events 07–08 omit the alternate C context; rechecked with both contexts"},
            "control_counts_reproduced":controls, "global_gain_shared_prefix":derivation["global_gain_shared_prefix"],
            "global_variable_extent_with_21bp_flanks_lengths":derivation["global_variable_extent_with_21bp_flanks_lengths"],
            "maximum_saved_read_length":max(len(r["sequence"]) for r in reads),
            "gains_with_complete_local_span_observed":sum(bool(c["segment_read_ordinals"]["A_C_B"]) for c in gained),
            "gains_with_variable_core_observed":sum(bool(c["segment_read_ordinals"]["A_C_B_variable_core"]) for c in gained),
            "gains_with_variable_core_unobserved":sum(not c["segment_read_ordinals"]["A_C_B_variable_core"] for c in gained),
            "saved_exact_q20_interval_coordinate_bounds":[
                min(i["candidate_span"][0] for r in per_read for ps in r["candidate_placements"].values() for p in ps for i in p["contiguous_q20_exact_intervals"]),
                max(i["candidate_span"][1] for r in per_read for ps in r["candidate_placements"].values() for p in ps for i in p["contiguous_q20_exact_intervals"])],
            "gains_with_full_candidate_phase_established":0,
            "observed_full_span_strings":groups, "candidate_block_combinations":combinations,
            "candidates":per_candidate, "reads":per_read,
            "limitations":[
                "The 25 reads were saved because they matched one of eight original Q20 windows. Missing combinations are unobserved in this selected ledger, not proven absent from the million-read prefix or sample.",
                "Every linked span must occur exactly on a single read with Q20 throughout. Separate or overlapping reads are never stitched into a linkage witness.",
                "Overlapping original windows are not independent alleles. Blocks describe literal strings and do not claim unique indel placement.",
                "All exact positions and both orientations are retained. Candidate context has 32 gains and one historical ITS assembly, which is not independent biological truth.",
                "A matching local span can be shared by several candidates. It does not phase the distant differences beginning at position 518 or the full 424–426 bp flanked variable extent.",
                "Variable-core co-occurrence requires all A/C/B blocks and intervening sequence on a single Q20 read but can have fewer than 21 flanking bases. Flank lengths are explicit; this is reported separately from the full flanked assay.",
                "Read counts are records, not independent molecules. Duplicated sequences, PCR copies, sequencing errors, other loci, and unsampled alternatives remain possible.",
                "Full-span candidate mismatches outside a supported interval remain recorded. An incompatible read does not rule out that candidate elsewhere in a mixed sample.",
                "No new FASTQ scan, read recruitment, gap alignment, whole-repeat reconstruction, or product filtering was performed."],
            "input_subset":raw["input_subset"],
            "source_receipts": {name:receipt(PILOT/name) for name in ["manifest.json","derivation.json","run-v1/raw.json","run-v1/execution.json","summary.json"]}}


OUTPUT_NAMES=["report.json","summary.json","candidate_disposition.tsv","block_combinations.tsv","checksums.json"]


def check_output_paths(output):
    existing=[name for name in OUTPUT_NAMES if (output/name).exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite saved evidence: {existing}")
    output.mkdir(parents=True,exist_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True,help="Fresh output directory; saved artifacts cannot be overwritten")
    output=parser.parse_args().output.resolve()
    check_output_paths(output)
    result = analyze()
    result["producer"] = receipt(Path(__file__))
    result["tests_source"] = receipt(BASE/"test_linkage.py")
    with (output/"report.json").open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True)+"\n")
    with (output/"candidate_disposition.tsv").open("x") as stream:
        fields=["candidate_id","role","length","A","C","B","A_C_reads","C_B_reads","A_C_B_reads","A_C_B_variable_core_reads","disposition","full_candidate_phase"]
        writer=csv.DictWriter(stream,fields,delimiter="\t");writer.writeheader()
        for c in result["candidates"]:
            writer.writerow({**{k:c[k] for k in ["candidate_id","role","length","disposition","full_candidate_phase"]},
                             **c["blocks"], **{k+"_reads":','.join(map(str,c["segment_read_ordinals"][k])) for k in ["A_C","C_B","A_C_B","A_C_B_variable_core"]}})
    with (output/"block_combinations.tsv").open("x") as stream:
        fields=["A","C","B","candidate_count","candidate_ids","q20_full_flanked_read_ordinals","q20_variable_core_read_ordinals","full_flanked_disposition","variable_core_disposition"]
        writer=csv.DictWriter(stream,fields,delimiter="\t");writer.writeheader()
        for c in result["candidate_block_combinations"]:
            core_ids=ordinal_ids(c["variable_core_placements"])
            writer.writerow({**c["blocks"],"candidate_count":len(c["candidate_ids"]),"candidate_ids":','.join(c["candidate_ids"]),
                             "q20_full_flanked_read_ordinals":','.join(map(str,c["q20_read_ordinals"])),
                             "q20_variable_core_read_ordinals":','.join(map(str,core_ids)),
                             "full_flanked_disposition":"observed" if c["q20_read_ordinals"] else "unobserved_in_selected_saved_ledger",
                             "variable_core_disposition":"observed" if core_ids else "unobserved_in_selected_saved_ledger"})
    summary={k:v for k,v in result.items() if k not in ["reads","candidates","candidate_block_combinations"]}
    with (output/"summary.json").open("x") as stream:
        stream.write(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    outputs={name:receipt(output/name) for name in OUTPUT_NAMES if name!="checksums.json"}
    with (output/"checksums.json").open("x") as stream:
        stream.write(json.dumps(outputs,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ["retained_ITS_reads","candidate_context_count","gains_with_complete_local_span_observed","gains_with_variable_core_observed","gains_with_full_candidate_phase_established"]}))


if __name__ == "__main__":
    main()

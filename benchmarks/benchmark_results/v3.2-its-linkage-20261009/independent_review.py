#!/usr/bin/env python3
"""Independent direct enumeration of the saved ITS linkage evidence."""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
RUN = BASE / "run-v1"
PILOT = BASE.parent / "v3.2-its-gain-pilot-20261009"
COMP = str.maketrans("ACGTN", "TGCAN")


def receipt(path):
    data = path.read_bytes()
    return {"path": str(path.resolve()), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def load(path):
    return json.loads(path.read_text())


def placements(sequence, word):
    return [i for i in range(len(sequence)-len(word)+1)
            if sequence[i:i+len(word)] == word]


def oriented(read):
    return [("+", read["sequence"], read["quality"]),
            ("-", read["sequence"].translate(COMP)[::-1], read["quality"][::-1])]


def matches(reads, word):
    return [(read["ordinal"], strand, i, i+len(word),
             min(ord(q)-33 for q in quality[i:i+len(word)]))
            for read in reads for strand, seq, quality in oriented(read)
            for i in placements(seq, word)
            if all(ord(q) >= 53 for q in quality[i:i+len(word)])]


def main():
    report = load(RUN / "report.json")
    manifest = load(PILOT / "manifest.json")
    raw = load(PILOT / "run-v1/raw.json")
    derivation = load(PILOT / "derivation.json")
    verification = load(RUN / "verification.json")
    for name, expected in load(RUN / "checksums.json").items():
        actual = receipt(RUN / name)
        assert all(actual[k] == expected[k] for k in ("sha256", "bytes")), name
    for expected in report["source_receipts"].values():
        actual = receipt(Path(expected["path"]))
        assert all(actual[k] == expected[k] for k in ("sha256", "bytes"))
    assert verification["status"] == "pass"
    assert verification["test_exit_code"] == verification["deterministic_replay_exit_code"] == 0
    assert all(verification["deterministic_replay_identical"].values())
    for path, key in [(BASE/"analyze_linkage.py", "producer_sha256"),
                      (BASE/"test_linkage.py", "test_source_sha256"),
                      (RUN/"tests.log", "test_log_sha256"),
                      (RUN/"checksums.json", "checksums_sha256")]:
        assert receipt(path)["sha256"] == verification[key]
    assert receipt(BASE/"test_linkage.py") == report["tests_source"]
    assert receipt(BASE/"analyze_linkage.py") == report["producer"]
    templates = {c["id"]: c["sequence"] for c in manifest["candidates"]}
    for candidate in manifest["candidates"]:
        assert hashlib.sha256(candidate["sequence"].encode()).hexdigest() == candidate["sequence_sha256"]
    reads = [x for x in raw["read_evidence"]
             if any(o["event_id"].startswith("its_local_") for o in x["observations"])]
    assert len(reads) == len({r["ordinal"] for r in reads}) == 25
    assert len(report["candidates"]) == 33
    for event in raw["events"]:
        ids = sorted({m[0] for m in matches(raw["read_evidence"], event["window_sequence"])})
        assert len(ids) == event["records"]
    anchors = report["anchors"]
    left, middle, inner, right = [anchors[k] for k in
                                 ("left", "between_A_C", "between_C_B", "right")]
    assays = 0
    for candidate in report["candidates"]:
        a,c,b = [candidate["blocks"][k] for k in ("A", "C", "B")]
        core = a+middle+c+inner+b
        words = {"A":left+a+middle[:21], "C":middle[-21:]+c+inner,
                 "B":inner+b+right, "A_C":left+a+middle+c+inner,
                 "C_B":middle[-21:]+c+inner+b+right,
                 "A_C_B":left+core+right, "A_C_B_variable_core":core}
        start,end = candidate["A_C_B_candidate_span"]
        assert templates[candidate["candidate_id"]][start:end] == words["A_C_B"]
        for name, word in words.items():
            ids = sorted({m[0] for m in matches(reads, word)})
            assert ids == candidate["segment_read_ordinals"][name], (candidate["candidate_id"], name)
            assays += 1
        assert candidate["full_candidate_phase"] == "unresolved"
    # Enumerate both anchors without using the producer's regular expression.
    anchored = []
    for read in reads:
        for strand, seq, quality in oriented(read):
            for i in placements(seq,left):
                for j in placements(seq,right):
                    if j >= i+len(left):
                        end = j+len(right)
                        anchored.append((read["ordinal"],strand,i,end,seq[i:end],
                                         min(ord(q)-33 for q in quality[i:end])))
    reported = [(p["ordinal"],p["strand"],*p["oriented_read_span"],p["sequence"],p["minimum_phred"])
                for group in report["observed_full_span_strings"] for p in group["all_placements"]]
    assert sorted(anchored) == sorted(reported)
    lookup = {r["ordinal"]:r for r in reads}
    core_placements = 0
    for combination in report["candidate_block_combinations"]:
        a,c,b = [combination["blocks"][k] for k in ("A", "C", "B")]
        core = a+middle+c+inner+b
        actual = matches(reads,core)
        stored = [(p["ordinal"],p["strand"],*p["oriented_read_span"],p["minimum_phred"])
                  for p in combination["variable_core_placements"]]
        assert sorted(actual) == sorted(stored)
        for p in combination["variable_core_placements"]:
            read = lookup[p["ordinal"]]
            strand, seq, quality = next(z for z in oriented(read) if z[0] == p["strand"])
            start,end = p["oriented_read_span"]
            assert p["original_read_span"] == ([start,end] if strand == "+" else [len(seq)-end,len(seq)-start])
            for side, word in (("left",left),("right",right)):
                possible = []
                for n in range(len(word)+1):
                    q = quality[start-n:start] if side == "left" else quality[end:end+n]
                    s = seq[start-n:start] if side == "left" else seq[end:end+n]
                    expected = word[len(word)-n:] if side == "left" else word[:n]
                    if len(s) == n and s == expected and all(ord(x)>=53 for x in q):
                        possible.append(n)
                assert max(possible) == p[f"consecutive_exact_q20_{side}_flank_bases"]
            core_placements += 1
    # Derive every candidate offset directly from seed occurrence pairs, then
    # derive compatible intervals from whole-overlap runs of equal Q20 bases.
    placement_count = 0
    for saved in report["reads"]:
        read = lookup[saved["ordinal"]]
        assert read["sequence"] == saved["sequence"] and read["quality"] == saved["quality"]
        for candidate in report["candidates"]:
            cid = candidate["candidate_id"]
            template = templates[cid]
            expected_keys = set()
            stored_by_key = {(p["strand"],p["candidate_start_offset"]):p
                             for p in saved["candidate_placements"].get(cid,[])}
            for strand,seq,quality in oriented(read):
                offsets = {}
                for win in derivation["windows"]:
                    word = win["sequence"]
                    for i in placements(seq,word):
                        if not all(ord(q)>=53 for q in quality[i:i+len(word)]):
                            continue
                        for j in placements(template,word):
                            offsets.setdefault(j-i,[]).append((i,i+len(word)))
                for offset,seeds in offsets.items():
                    key = (strand,offset)
                    expected_keys.add(key)
                    p = stored_by_key[key]
                    lo,hi = max(0,-offset),min(len(seq),len(template)-offset)
                    assert p["oriented_read_overlap"] == [lo,hi]
                    assert p["candidate_overlap"] == [lo+offset,hi+offset]
                    assert p["overlap_minimum_phred"] == min(ord(q)-33 for q in quality[lo:hi])
                    mismatches = [(i,i+offset,seq[i],template[i+offset],ord(quality[i])-33)
                                  for i in range(lo,hi) if seq[i] != template[i+offset]]
                    assert mismatches == [(v["oriented_read_position"],v["candidate_position"],
                                           v["observed"],v["candidate"],v["phred"]) for v in p["mismatches"]]
                    assert p["exact_overlap"] == (not mismatches)
                    assert p["full_read_within_candidate"] == (lo==0 and hi==len(seq))
                    runs=[]
                    start=None
                    for i in range(lo,hi+1):
                        good=i<hi and seq[i]==template[i+offset] and ord(quality[i])>=53
                        if good and start is None:
                            start=i
                        if not good and start is not None:
                            if any(start<=a and b<=i for a,b in seeds):
                                runs.append((start,i))
                            start=None
                    actual_runs = [tuple(v["oriented_read_span"]) for v in p["contiguous_q20_exact_intervals"]]
                    assert sorted(runs) == sorted(actual_runs)
                    for interval in p["contiguous_q20_exact_intervals"]:
                        assert interval["candidate_span"] == [v+offset for v in interval["oriented_read_span"]]
                    placement_count += 1
            assert expected_keys == set(stored_by_key)
    gains = [c for c in report["candidates"] if c["role"] == "gained_uncertain"]
    full = sum(bool(c["segment_read_ordinals"]["A_C_B"]) for c in gains)
    core = sum(bool(c["segment_read_ordinals"]["A_C_B_variable_core"]) for c in gains)
    assert len(gains) == 32 and full == 11 and core == 15
    assert report["gains_with_complete_local_span_observed"] == full
    assert report["gains_with_variable_core_observed"] == core
    assert report["gains_with_full_candidate_phase_established"] == 0
    result = {"status":"pass", "scope":"Independent review of frozen run-v1; no FASTQ rescan",
              "candidate_segment_assays_recounted":assays, "candidate_contexts":33,
              "saved_ITS_reads":25, "candidate_placements_reconstructed":placement_count,
              "complete_anchor_placements_recounted":len(anchored),
              "variable_core_placements_recounted":core_placements,
              "gains_with_full_flanked_local_span":full,
              "gains_with_variable_core_including_partial_flanks":core,
              "full_candidate_phase_established":0,
              "controls":"All eight ITS and two CO1 original event counts reproduce; CO1 27/0 preserved",
              "synthetic_verification":"11 reported tests pass; source/log receipts and deterministic replay verified",
              "conclusion":"Nearby exact Q20 co-occurrence is observed for 15 gains (11 with full 21-base flanks); the remaining 17 are unobserved in this selected saved ledger. All full hypotheses remain unresolved.",
              "production_fix_assessment":"No concrete production failure demonstrated: bounded uncertain output does not promise exhaustive candidates or validated phase; a novel local read string is not evidence that emission or admission is broken.",
              "limitations":["The saved reads were selected for original local positive windows; zeros cannot estimate prevalence or establish absence in the prefix.",
                             "Observed literal spans may match multiple hypotheses and do not establish biological target, independent molecules, polymorphism, correct repeats, or full-product phase.",
                             "Direct exact ungapped placement does not exhaust alternate gapped alignments; the report preserves mismatches and ambiguous seed placements."],
              "artifacts":{str(p.relative_to(BASE)):receipt(p) for p in
                           [BASE/"analyze_linkage.py",BASE/"test_linkage.py",RUN/"report.json",RUN/"summary.json",
                            RUN/"candidate_disposition.tsv",RUN/"block_combinations.tsv",RUN/"checksums.json",
                            RUN/"verification.json",RUN/"tests.log",Path(__file__)]}}
    with (BASE/"independent_review.json").open("x") as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","candidate_segment_assays_recounted","candidate_placements_reconstructed","gains_with_full_flanked_local_span","gains_with_variable_core_including_partial_flanks")}))


if __name__ == "__main__":
    main()

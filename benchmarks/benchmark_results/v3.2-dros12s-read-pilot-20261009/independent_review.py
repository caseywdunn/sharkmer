#!/usr/bin/env python3
"""Review frozen assay scope, source identities and saved read evidence."""
import csv
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent
COMP = str.maketrans("ACGTN", "TGCAN")


def sha(value):
    return hashlib.sha256(value).hexdigest()


def receipt(path):
    data = path.read_bytes()
    return {"path":str(path.resolve()),"bytes":len(data),"sha256":sha(data)}


def load(name):
    return json.loads((BASE/name).read_text())


def positions(sequence, word):
    return [i for i in range(len(sequence)-len(word)+1)
            if sequence[i:i+len(word)] == word]


def main():
    manifest,derivation,raw,execution,summary,controls = [load(n) for n in
        ("manifest.json","derivation.json","run-v1/raw.json","run-v1/execution.json","summary.json","prefilter_controls.json")]
    assert execution["status"] == "pass" and execution["sources_unchanged"]
    assert execution["returncode"] == 0 and execution["time_limit_seconds"] == 600
    assert execution["source_hashes_before"] == execution["source_hashes_after"]
    for path,expected in execution["source_hashes_after"].items():
        assert receipt(Path(path))["sha256"] == expected,path
    for name,expected in execution["outputs"].items():
        actual=receipt(BASE/"run-v1"/name)
        assert all(actual[k]==expected[k] for k in ("bytes","sha256")),name
    for name,expected in summary["source_receipts"].items():
        assert receipt(BASE/name)["sha256"]==expected,name
    for path,expected in raw["provenance"].items():
        assert receipt(Path(path))["sha256"]==expected,path
    assert raw["input_subset"]==load("expected_prefix.json")==summary["verified_input_subset"]
    assert raw["input_subset"]["records"]==1_000_000
    assert raw["input_subset"]["bases"]==101_000_000
    assert raw["maximum_read_length"]==101 and raw["minimum_quality"]==20
    assert raw["max_substitutions"]==0
    assert raw["prefilter"]["records_examined"]==1_000_000
    assert controls["status"]=="pass" and controls["exact_observation_equality_cases"]==1835
    for name,expected in controls["source_sha256"].items():
        assert receipt(BASE/name)["sha256"]==expected
    candidates={c["id"]:c["sequence"] for c in manifest["candidates"]}
    assert len(candidates)==len(manifest["candidates"])==28
    for c in manifest["candidates"]:
        assert sha(c["sequence"].encode())==c["sequence_sha256"]
    current={c["candidate_id"]:c for c in derivation["candidates"]}
    assert len(current)==26 and sum(c["gained"] for c in current.values())==23
    # Bind each source FASTA header to its own sequence and emitted hash.
    archive_name,member_name=next((key.split("/",1) for key in derivation["archive_members"]
                                   if key.endswith("insecta_12S.fasta")))
    archive=BASE.parent/archive_name/"evidence.tar.gz"
    with tarfile.open(archive,"r:gz") as tar:
        data=tar.extractfile(member_name).read()
    expected=derivation["archive_members"][archive_name+"/"+member_name]
    assert len(data)==expected["bytes"] and sha(data)==expected["sha256"]
    fasta=[]
    for block in data.decode().split(">"):
        if not block.strip():
            continue
        header,*lines=block.splitlines()
        sequence="".join(lines)
        assert "sequence_sha256="+sha(sequence.encode()) in header
        fasta.append(sequence)
    assert sorted(fasta)==sorted(candidates[c] for c in current)
    old=load("old_control_manifest.json")
    historical=[c for c in old["candidates"] if c["id"] in derivation["historical_12s_aliases"]]
    assert sorted(len(c["sequence"]) for c in historical)==[475,487,499]
    for c in historical:
        assert candidates[derivation["historical_12s_aliases"][c["id"]]]==c["sequence"]
    events={e["id"]:e for e in raw["events"]}
    assert len(events)==len(manifest["events"])==139
    patterns={}
    for e in manifest["events"]:
        word=candidates[e["candidate_id"]][e["start"]-21:e["end"]+21]
        saved=events[e["id"]]
        assert word==saved["window_sequence"]
        assert saved["window_start"]==e["start"]-21 and saved["window_end"]==e["end"]+21
        patterns[e["id"]]=[("+",word),("-",word.translate(COMP)[::-1])]
        context=[]
        for cid,sequence in candidates.items():
            for strand,oriented in (("+",sequence),("-",sequence.translate(COMP)[::-1])):
                context.extend({"candidate_id":cid,"strand":strand,"start":i,"substitutions":0}
                               for i in positions(oriented,word))
        assert context==saved["template_placements"]
        assert saved["unique_template_placement"]==(len(context)==1)
        assert saved["assay_callability"]==("longer_than_every_read" if len(word)>101 else "assayed")
    recounted=Counter()
    ledger_matches=0
    for read in raw["read_evidence"]:
        observed={o["event_id"]:o for o in read["observations"]}
        expected_ids=set()
        for eid,oriented_patterns in patterns.items():
            wanted=[]
            for strand,word in oriented_patterns:
                for start in positions(read["sequence"],word):
                    end=start+len(word)
                    if all(ord(q)>=53 for q in read["quality"][start:end]):
                        wanted.append((strand,start,end))
            if wanted:
                expected_ids.add(eid)
                assert wanted==[(m["strand"],m["read_start"],m["read_end"]) for m in observed[eid]["matches"]]
                recounted[eid]+=1
                ledger_matches+=len(wanted)
        assert expected_ids==set(observed)
    assert len(raw["read_evidence"])==135
    for eid,e in events.items():
        assert recounted[eid]==e["records"]
        assert e["unfiltered_exact_records"]>=e["records"]
    for eid,expected in derivation["expected_co1_controls"].items():
        assert all(events[eid][k]==v for k,v in expected.items())
    assert len(derivation["expected_co1_controls"])==16
    assert len(derivation["windows"])==117
    tiles=Counter()
    nearest_memberships=0
    for cid,c in current.items():
        sequence=candidates[cid]
        assert len(sequence)==c["length"] and sha(sequence.encode())==c["sha256"]
        seen=[]
        for window in derivation["windows"]:
            assert window["sequence"]==events[window["event_id"]]["window_sequence"]
            assert events[window["event_id"]]["records"]==events[window["event_id"]]["unfiltered_exact_records"]==0
            for member in window["memberships"]:
                if member["candidate_id"]!=cid:
                    continue
                start,end=member["window_start"],member["window_end"]
                assert sequence[start:end]==window["sequence"]
                if member["role"]=="backbone_tile":
                    assert end-start==61
                    seen.append(start)
                    tiles[cid]+=1
                elif member["role"]=="nearest_local_difference":
                    nearest_memberships+=1
                    assert 42<=end-start<=101
                else:
                    raise AssertionError(member["role"])
        assert sorted(seen)==sorted(set(range(0,len(sequence)-60,21))|{len(sequence)-61})
    assert sum(tiles.values())==568 and nearest_memberships==52
    assert summary["current_12s_candidates"]==26
    assert summary["q20_current_12s_strings_observed"]==summary["unfiltered_current_12s_strings_observed"]==0
    assert sorted(r["candidate_id"] for r in summary["per_candidate"])==sorted(current)
    for name,length in (("per_candidate.tsv",26),("per_event.tsv",139)):
        with (BASE/name).open() as stream:
            assert len(list(csv.DictReader(stream,delimiter="\t")))==length
    structural=[e for e in raw["events"] if e["role"]=="primary_frozen_repeat_interval"]
    assert len(structural)==6
    assert sorted(len(e["window_sequence"]) for e in structural if e["assay_callability"]=="longer_than_every_read")==[123,123,146]
    assert summary["verification"]["supplementary_raw_counter_equivalence_cases"]==835
    result={"status":"pass","scope":"Frozen exact occurrence pilot; source/control review and independent saved-ledger recount, no full-prefix rescan",
            "current_12s_candidates_verified":26,"gained_candidates":23,
            "historical_full_sequence_lengths_verified":[475,487,499],
            "events_verified":139,"local_strings_verified":117,
            "all_candidate_tile_memberships_verified":568,"nearest_endpoint_memberships_verified":52,
            "positive_CO1_controls_reproduced":16,"CO1_Q20_record_range":[9,24],
            "saved_Q20_reads_independently_recounted":135,"exact_match_placements_recounted":ledger_matches,
            "Q20_and_unfiltered_positive_12s_strings":0,
            "uncallable_structural_window_lengths":[123,123,146],
            "prefilter_completeness_review":"Every exact oriented pattern contains its first19 bases. Literal equal-length seed alternatives cannot shadow distinct same-position seeds; overlapping lookahead retains every position. Event ordering is preserved when calling the unchanged scanner.",
            "verification_evidence":{"direct_scanner_prefilter_equality_cases":1835,"supplemental_raw_count_equality_cases":835},
            "interpretation":"The recovered 487/499 and historical475 sequences, together with all26current12S hypotheses, lack corroboration in the assayed local exact strings. Quality exclusion does not explain those zeros; mismatches, indels, sampling, target identity and whole-product validity remain unresolved.",
            "limits":["An independent saved-ledger recount verifies positive observations but cannot independently establish which records were omitted from the prefix. Source completeness proof, frozen execution, consumed-prefix hash and positive controls support the full-prefix negative result.",
                      "Supplemental quality-unfiltered positive reads were not all saved, so those counts cannot all be reconstructed from the Q20 ledger. Their counter was source-reviewed and tested independently on835synthetic cases.",
                      "The CO1 controls establish successful scanning of matching sample reads; no standard12S positive sequence exists in this comparison.",
                      "No genotype, locus, correct repeat length, biological absence or filtering policy follows from these bounded exact-string observations."],
            "production_fix_assessment":"No sequence-spelling defect demonstrated by this assay. A separately recorded bounded19mer invariant check can investigate observed graph words versus unsupported longer combinations.",
            "artifacts":{name:receipt(BASE/name) for name in
                         ("prepare_pilot.py","manifest.json","derivation.json","expected_prefix.json","scan_prefiltered.py",
                          "check_prefilter.py","prefilter_controls.json","run_pilot.py","summarize_pilot.py","summary.json",
                          "per_candidate.tsv","per_event.tsv","run-v1/raw.json","run-v1/execution.json","run-v1/time.json","independent_review.py")}}
    with (BASE/"independent_review.json").open("x") as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","current_12s_candidates_verified","events_verified","saved_Q20_reads_independently_recounted","Q20_and_unfiltered_positive_12s_strings")}))


if __name__=="__main__":
    main()

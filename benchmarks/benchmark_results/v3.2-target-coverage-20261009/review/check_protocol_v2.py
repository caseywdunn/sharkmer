#!/usr/bin/env python3
"""Independently reconstruct all frozen tile and contrast memberships."""
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
COMP=str.maketrans("ACGT","TGCA")
FOCUS={"SRR25378184":{"angiospermae_psbA-trnH":49,"angiospermae_trnV-atpE":41,"angiospermae_atpB-rbcL":1},
       "ERR14009273":{"angiospermae_trnL-F":1},"SRR1057608":{"insecta_ND4":44},"SRR27962769":{"insecta_18S_1":18}}
CONTROLS={"SRR25378184":["angiospermae_psbM-trnD","angiospermae_rpl36-infA-rps8"],
          "ERR14009273":["angiospermae_ITS","angiospermae_atpB-rbcL"],
          "SRR1057608":["insecta_CO1_1","insecta_CO2_1"],"SRR27962769":["insecta_CO1_1","insecta_CO2_1"]}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def receipt(path):
    data=path.read_bytes()
    return {"sha256":sha(data),"bytes":len(data)}


def load(path):
    return json.loads(path.read_text())


def identifier(row):
    return "|".join(row[k] for k in ("accession","gene","sha256"))


def diff(a,b):
    prefix=next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),min(len(a),len(b)))
    suffix=0
    for i in range(1,min(len(a),len(b))-prefix+1):
        if a[-i]!=b[-i]:
            break
        suffix=i
    return prefix,len(a)-suffix,len(b)-suffix


def main():
    inventory=load(BASE/"matrix/sequence_inventory.json")
    matrix_validation=load(BASE/"matrix/validation.json")
    prefixes={row["accession"]:row["prefix"] for row in matrix_validation["current_runs"]}
    synthetic=load(HERE/"exact_independent_controls_v2.json")
    assert synthetic["status"]=="pass" and synthetic["sources_unchanged"]
    for path,expected in synthetic["source_hashes_after"].items():
        assert sha(Path(path).read_bytes())==expected
    results=[]
    artifacts={}
    for accession,expected_targets in FOCUS.items():
        path=BASE/"reads/manifests-v2"/(accession+".json")
        manifest=load(path)
        artifacts[str(path.relative_to(BASE))]=receipt(path)
        for source,expected in manifest["source_files"].items():
            assert receipt(Path(source))==expected,source
        assert manifest["expected_prefix"]==prefixes[accession]
        assert manifest["accession"]==accession and manifest["focus_targets"]==expected_targets
        context={identifier(row):row for row in inventory if row["accession"]==accession}
        candidates={row["id"]:row for row in manifest["candidates"]}
        assert set(candidates)==set(context)
        for cid,row in candidates.items():
            assert all(row[k]==context[cid][k] for k in ("sequence","sha256","gene","assembly_status"))
            assert sha(row["sequence"].encode())==row["sha256"]
        focus={cid:row for cid,row in context.items() if row["gene"] in expected_targets}
        assert Counter(row["gene"] for row in focus.values())==expected_targets
        assert all(row["new_vs_release"] and row["newly_represented_target"] for row in focus.values())
        assert set(manifest["focused_candidate_ids"])==set(focus)
        controls={}
        for gene in CONTROLS[accession]:
            row=min((r for r in context.values() if r["gene"]==gene and r["assembly_status"]=="standard"),key=lambda r:r["sha256"])
            assert row["historical_retained"]
            controls[identifier(row)]=row
        assert set(manifest["positive_control_ids"])==set(controls)
        expected_members=defaultdict(list)
        expected_uncallable=[]
        contrasts=[]

        def add(cid,start,end,role,extra=None):
            row=context[cid]
            member={"candidate_id":cid,"gene":row["gene"],"start":start,"end":end,"role":role,**(extra or {})}
            if start<0 or end>len(row["sequence"]):
                expected_uncallable.append({**member,"reason":"missing_complete_21base_flanks"})
                return
            word=row["sequence"][start:end]
            reverse=word.translate(COMP)[::-1]
            canonical=min(word,reverse)
            member["pattern_strand_in_candidate"]="+" if word==canonical else "-"
            expected_members[canonical].append(member)

        for cid,row in {**focus,**controls}.items():
            starts=sorted(set(range(0,len(row["sequence"])-60,21))|{len(row["sequence"])-61})
            assert starts[0]==0 and starts[-1]+61==len(row["sequence"])
            for start in starts:
                add(cid,start,start+61,"candidate_backbone_tile" if cid in focus else "positive_control_tile")
        for cid,row in focus.items():
            others=[(oid,r) for oid,r in focus.items() if r["gene"]==row["gene"] and oid!=cid]
            if not others:
                continue
            options=[]
            for oid,other in others:
                start,a_end,b_end=diff(row["sequence"],other["sequence"])
                options.append(((max(a_end-start,b_end-start),a_end+b_end-2*start,other["sha256"]),oid,start,a_end,b_end))
            _,oid,start,a_end,b_end=min(options)
            assert row["sequence"][:start]+context[oid]["sequence"][start:b_end]+row["sequence"][a_end:]==context[oid]["sequence"]
            contrast_id=cid+"->"+oid
            contrasts.append((contrast_id,cid,oid,start,a_end,b_end))
            add(cid,start-21,a_end+21,"nearest_difference_endpoint",{"contrast_id":contrast_id,"endpoint":"candidate"})
            add(oid,start-21,b_end+21,"nearest_difference_endpoint",{"contrast_id":contrast_id,"endpoint":"partner"})
        observed_contrasts=[tuple(c[k] for k in ("id","candidate_id","partner_id","start","candidate_end","partner_end")) for c in manifest["nearest_contrasts"]]
        assert sorted(contrasts)==sorted(observed_contrasts)
        sort_members=lambda members:sorted(json.dumps(m,sort_keys=True) for m in members)
        assert sort_members(expected_uncallable)==sort_members(manifest["uncallable_endpoints"])
        actual_events={e["sequence"]:e for e in manifest["events"]}
        assert len(actual_events)==len(manifest["events"]) and set(actual_events)==set(expected_members)
        for word,members in expected_members.items():
            event=actual_events[word]
            assert event["sha256"]==sha(word.encode()) and event["id"]=="window_"+event["sha256"][:16]
            assert sort_members(event["memberships"])==sort_members(members)
        assert len(set(e["id"] for e in manifest["events"]))==len(manifest["events"])
        assert manifest["protocol"]=={"tile_length":61,"tile_step":21,"last_tile_end_aligned":True,"contrast_flank":21,
                                     "minimum_phred":20,"supplementary_quality_unfiltered_counts":True,"exact_only":True,
                                     "both_orientations":True,"all_overlapping_occurrences":True,
                                     "maximum_seconds":600,"maximum_address_space_bytes":2147483648,"cpus":1}
        results.append({"accession":accession,"focused_candidates":len(focus),"full_sample_context":len(context),
                        "unchanged_standard_controls":len(controls),"events":len(actual_events),
                        "tile_memberships":sum(m["role"]!="nearest_difference_endpoint" for members in expected_members.values() for m in members),
                        "nearest_contrasts":len(contrasts),"incomplete_flank_endpoints":len(expected_uncallable),
                        "minimum_window_length":min(map(len,actual_events)),"maximum_window_length":max(map(len,actual_events))})
    for name in ("reads/scan_exact_v2.py","reads/prepare_assays_v2.py","reads/run_assay_v2.py","matrix/sequence_inventory.json","matrix/validation.json",
                 "review/check_exact_independent_v2.py","review/exact_independent_controls_v2.json","review/check_protocol_v2.py"):
        artifacts[name]=receipt(BASE/name)
    result={"status":"approved_for_bounded_execution","samples":results,"focused_candidates":sum(r["focused_candidates"] for r in results),
            "review":"All four manifests independently reconstructed from the frozen complete current inventory. All tiles, contrasts, canonical orientations, endpoint memberships, explicit incomplete-flank exclusions, full sample contexts, control membership and assembly-prefix identities match.",
            "source_review":"Exact matcher uses complete overlapping prefix-seed selection; counts continue after every representative cap. Wrapper enforces one CPU,600 seconds and2GiB and freezes source identities before/after. Both strand matches use the physical read quality span; palindromes count once per position.",
            "conditions":["Each predefined control must have at least one Q20 exact tile before interpreting that sample's negative findings; these new schedules do not have preregistered numerical counts.",
                          "Callability must use actual length histogram; overlength windows and missing complete flanks are uncallable.",
                          "Report complete counts separately from capped examples, including cap reasons and per-event completeness. Saved records are not a complete or representative read population.",
                          "Observed tiles and contrast endpoints do not certify target identity, biological validity or whole-product linkage."],
            "artifacts":artifacts}
    assert result["focused_candidates"]==154
    with (HERE/"protocol_review_v2.json").open("x") as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"focused_candidates":154,"samples":results}))


if __name__=="__main__":
    main()

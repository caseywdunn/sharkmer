#!/usr/bin/env python3
"""Independent adversarial controls for the bounded exact-read aggregator."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import random
import tempfile
import time
import unittest

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parents[2]
SCANNER=BASE/"reads/scan_exact_v2.py"
ORIGINAL=ROOT/"scripts/audit_read_bridges.py"
COMPLEMENT=str.maketrans("ACGT","TGCA")


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


audit=module("new_exact_scanner",SCANNER)
original=module("original_exact_scanner",ORIGINAL)


def sha(value):
    return hashlib.sha256(value).hexdigest()


def rc(word):
    return word.translate(COMPLEMENT)[::-1]


def event(word,identifier="e"):
    return {"id":identifier,"sequence":word,"sha256":sha(word.encode())}


def read(word,ordinal=1,quality=None):
    if quality is None:
        quality="I"*len(word)
    return ordinal,f"@read.{ordinal}".encode(),word.encode(),quality.encode()


def brute_observations(sequence,quality,events):
    result=[]
    for index,e in enumerate(events):
        word=e["sequence"]
        reverse=rc(word)
        patterns=[("both",word)] if word==reverse else [("+",word),("-",reverse)]
        found=[]
        for strand,pattern in patterns:
            for start in range(len(sequence)-len(pattern)+1):
                if sequence[start:start+len(pattern)]==pattern:
                    minimum=min(ord(q)-33 for q in quality[start:start+len(pattern)])
                    found.append({"strand":strand,"start":start,"end":start+len(pattern),
                                  "minimum_phred":minimum,"q20":minimum>=20})
        if found:
            result.append({"event_index":index,"matches":found})
    return result


def brute_counts(records,events):
    result=[{"raw_records":0,"q20_records":0,"raw_only_records":0,
             "raw_occurrences":0,"q20_occurrences":0} for _ in events]
    for _,_,s,q in records:
        for observation in brute_observations(s.decode(),q.decode(),events):
            row=result[observation["event_index"]]
            matches=observation["matches"]
            passing=sum(m["q20"] for m in matches)
            row["raw_records"]+=1
            row["q20_records"]+=bool(passing)
            row["raw_only_records"]+=not bool(passing)
            row["raw_occurrences"]+=len(matches)
            row["q20_occurrences"]+=passing
    return result


class IndependentControls(unittest.TestCase):
    word="ACGTTGAGCATCCGATGTTACGCTTAGCGATTCGAGGCTCAAT"

    def assert_counts(self,records,events,limits):
        result=audit.summarize_records(iter(records),events,limits)
        expected=brute_counts(records,events)
        for observed,correct in zip(result["events"],expected):
            for key,value in correct.items():
                self.assertEqual(observed[key],value,key)
        ledger=result["read_evidence"]
        self.assertEqual(len({r["ordinal"] for r in ledger}),len(ledger))
        self.assertLessEqual(len(ledger),limits["reads"])
        self.assertEqual(result["ledger"]["stored_bases"],sum(len(r["sequence"]) for r in ledger))
        matches=sum(len(o["matches"]) for r in ledger for o in r["observations"])
        self.assertEqual(result["ledger"]["stored_match_positions"],matches)
        self.assertLessEqual(matches,limits["matches"])
        self.assertLessEqual(result["ledger"]["stored_bases"],limits["bases"])
        retained_q=[0]*len(events)
        retained_raw=[0]*len(events)
        for row in ledger:
            self.assertEqual(row["observations"],brute_observations(row["sequence"],row["quality"],events))
            for observation in row["observations"]:
                i=observation["event_index"]
                if any(m["q20"] for m in observation["matches"]):
                    retained_q[i]+=1
                else:
                    retained_raw[i]+=1
        for i,e in enumerate(result["events"]):
            self.assertEqual(e["retained_q20_records"],retained_q[i])
            self.assertEqual(e["retained_raw_only_records"],retained_raw[i])
            self.assertEqual(e["example_ledger_complete_for_event"],e["raw_records"]==retained_q[i]+retained_raw[i])
        return result

    def test_randomized_matcher_equals_independent_brute_force(self):
        randomizer=random.Random(71815)
        words=["".join(randomizer.choices("ACGT",k=n)) for n in (19,20,42,43,61,74)]
        words += ["A"*42,"AT"*21,self.word]
        events=[event(word,f"e{i}") for i,word in enumerate(words)]
        matcher=audit.Matcher(events)
        for iteration in range(600):
            word=randomizer.choice(words)
            if iteration%2:
                word=rc(word)
            sequence="".join(randomizer.choices("ACGTN",k=iteration%11))+word
            if iteration%3==0:
                sequence+=word
            sequence+="".join(randomizer.choices("ACGTN",k=iteration%17))
            quality="".join(randomizer.choices("!45I",weights=[1,1,3,30],k=len(sequence)))
            expected=brute_observations(sequence,quality,events)
            self.assertEqual(matcher.observations(sequence.encode(),quality.encode()),expected)
            self.assertEqual(matcher.observations(sequence.encode(),quality.encode(),False),expected)

    def test_palindromic_overlaps_are_physical_occurrences(self):
        records=[read("AT"*23)]
        result=self.assert_counts(records,[event("AT"*21)],{"reads":10,"bases":1000,"matches":1000})
        self.assertEqual(result["events"][0]["raw_occurrences"],3)
        self.assertEqual(result["events"][0]["raw_records"],1)
        self.assertTrue(all(m["strand"]=="both" for m in result["read_evidence"][0]["observations"][0]["matches"]))

    def test_quality_applies_to_each_complete_physical_span(self):
        sequence=self.word+"NN"+rc(self.word)
        quality="5"*len(self.word)+"II"+"4"+"I"*(len(self.word)-1)
        result=self.assert_counts([read(sequence,quality=quality)],[event(self.word)],{"reads":10,"bases":1000,"matches":1000})
        self.assertEqual((result["events"][0]["raw_records"],result["events"][0]["q20_records"]),(1,1))
        self.assertEqual((result["events"][0]["raw_occurrences"],result["events"][0]["q20_occurrences"]),(2,1))
        self.assertEqual(result["events"][0]["raw_only_records"],0)

    def test_example_quotas_do_not_truncate_counts(self):
        records=[read(self.word,i,"I"*len(self.word) if i<=30 else "!"*len(self.word)) for i in range(1,41)]
        result=self.assert_counts(records,[event(self.word)],{"reads":100,"bases":100000,"matches":100000})
        e=result["events"][0]
        self.assertEqual((e["raw_records"],e["q20_records"],e["raw_only_records"]),(40,30,10))
        self.assertEqual((e["retained_q20_records"],e["retained_raw_only_records"]),(8,2))
        self.assertEqual([r["ordinal"] for r in result["read_evidence"]],[*range(1,9),31,32])
        self.assertFalse(e["example_ledger_complete_for_event"])

    def test_global_read_cap_zero_keeps_all_counts(self):
        records=[read(self.word,i) for i in range(1,21)]
        result=self.assert_counts(records,[event(self.word)],{"reads":0,"bases":100000,"matches":100000})
        self.assertEqual(result["events"][0]["ledger_cap_omitted_requests"],20)
        self.assertEqual(result["ledger"]["cap_omitted_requests_by_reason"],{"reads":20})

    def test_global_base_cap_keeps_all_counts(self):
        records=[read(self.word,i) for i in range(1,21)]
        result=self.assert_counts(records,[event(self.word)],{"reads":100,"bases":len(self.word),"matches":100000})
        self.assertEqual(result["ledger"]["stored_reads"],1)
        self.assertEqual(result["events"][0]["ledger_cap_omitted_requests"],19)

    def test_global_match_cap_counts_full_observation_payload(self):
        events=[event("A"*42),event("A"*43,"longer")]
        records=[read("A"*50,1),read("A"*42,2)]
        result=self.assert_counts(records,events,{"reads":100,"bases":100000,"matches":2})
        self.assertEqual([r["ordinal"] for r in result["read_evidence"]],[2])
        self.assertEqual([e["ledger_cap_omitted_requests"] for e in result["events"]],[1,1])

    def test_shared_read_bonus_preserves_all_event_observations(self):
        longer=self.word+"AACCGGTTACGATCG"
        events=[event(self.word),event(longer,"longer")]
        records=[read(self.word,i) for i in range(1,9)]+[read(longer,i) for i in range(9,17)]
        result=self.assert_counts(records,events,{"reads":100,"bases":100000,"matches":100000})
        self.assertEqual([e["retained_q20_records"] for e in result["events"]],[16,8])
        self.assertEqual(result["ledger"]["stored_reads"],16)

    def test_length_callability_counts_short_records_honestly(self):
        records=[read("A"*20,1),read("C"*50,2),read("G"*80,3)]
        events=[event("T"*n,str(n)) for n in (19,42,61,90)]
        result=self.assert_counts(records,events,{"reads":100,"bases":100000,"matches":100000})
        self.assertEqual([e["eligible_length_records"] for e in result["events"]],[3,2,1,0])
        self.assertEqual([e["callability"] for e in result["events"]],["assayed"]*3+["longer_than_every_read"])

    def test_candidate_placements_preserve_orientations_and_repeats(self):
        sequences=["NN",self.word+"AAA"+self.word,rc(self.word),"AT"*25]
        candidates=[{"id":str(i),"sequence":s,"sha256":sha(s.encode())} for i,s in enumerate(sequences[1:])]
        events=[event(self.word),event("AT"*21,"palindrome")]
        actual=audit.prepare_placements({"candidates":candidates,"events":events})
        expected=[]
        for e in events:
            places=[]
            for c in candidates:
                seq=c["sequence"]
                orientations=[("both",seq)] if seq==rc(seq) else [("+",seq),("-",rc(seq))]
                for strand,oriented in orientations:
                    for start in range(len(oriented)-len(e["sequence"])+1):
                        if oriented[start:start+len(e["sequence"])]==e["sequence"]:
                            places.append({"candidate_id":c["id"],"strand":strand,"start":start,"end":start+len(e["sequence"])})
            expected.append(places)
        self.assertEqual(actual,expected)

    def test_q20_positions_equal_original_scanner_after_palindrome_normalization(self):
        words=[self.word,"AT"*21,"ACGT"*15+"A"]
        candidates=[{"id":str(i),"sequence":w,"sequence_sha256":sha(w.encode())} for i,w in enumerate(words)]
        old_events=[{"id":str(i),"candidate_id":str(i),"start":21,"end":len(w)-21} for i,w in enumerate(words)]
        prepared=original.prepare_events({"schema_version":1,"candidates":candidates,"events":old_events},0)
        events=[event(w,str(i)) for i,w in enumerate(words)]
        for word in words:
            for sequence in ("AA"+word+"CCC","AAA"+rc(word)+"CC",word+word):
                for quality in ("I"*len(sequence),"5"*len(sequence),"!"*len(sequence)):
                    old=original.observations_for_read(sequence.encode(),quality.encode(),prepared,0)
                    new=audit.Matcher(events).observations(sequence.encode(),quality.encode())
                    normalized_old={o["event_id"]:{(m["read_start"],m["read_end"]) for m in o["matches"]} for o in old}
                    normalized_new={str(o["event_index"]):{(m["start"],m["end"]) for m in o["matches"] if m["q20"]} for o in new}
                    normalized_new={k:v for k,v in normalized_new.items() if v}
                    self.assertEqual(normalized_old,normalized_new)

    def test_scan_attests_prefix_and_rejects_changed_expected_hash(self):
        sequence=self.word.lower()
        data=(f"@r\n{sequence}\n+\n{'I'*len(sequence)}\n").encode()
        expected={"records":1,"bases":len(sequence),"size_bytes":len(data),"sha256":sha(data)}
        manifest={"accession":"synthetic","expected_prefix":expected,"candidates":[{"id":"a","sequence":self.word,"sha256":sha(self.word.encode())}],"events":[event(self.word)]}
        with tempfile.TemporaryDirectory(prefix="sharkmer-independent-exact-") as directory:
            path=Path(directory)/"reads.fastq"
            path.write_bytes(data)
            result=audit.scan(path,manifest)
            self.assertEqual(result["events"][0]["q20_records"],1)
            self.assertEqual(result["verified_prefix"],expected)
            manifest["expected_prefix"]={**expected,"sha256":"0"*64}
            with self.assertRaisesRegex(ValueError,"identity"):
                audit.scan(path,manifest)
            manifest["expected_prefix"]={**expected,"records":2}
            with self.assertRaisesRegex(ValueError,"Truncated"):
                audit.scan(path,manifest)


def main():
    paths=[SCANNER,ORIGINAL,Path(__file__)]
    before={str(p.resolve()):sha(p.read_bytes()) for p in paths}
    log=io.StringIO()
    start=time.perf_counter()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(IndependentControls)
    outcome=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
    after={str(p.resolve()):sha(p.read_bytes()) for p in paths}
    result={"status":"pass" if outcome.wasSuccessful() and before==after else "fail",
            "tests_run":outcome.testsRun,"errors":len(outcome.errors),"failures":len(outcome.failures),
            "wall_seconds":time.perf_counter()-start,"randomized_brute_force_matcher_cases":600,
            "original_scanner_position_comparison_cases":27,
            "source_hashes_before":before,"source_hashes_after":after,"sources_unchanged":before==after,
            "scope":"Synthetic independent exact matching, raw/Q20 counts, orientation/overlap, quality, callability, capped ledger and prefix-attestation controls; no experimental FASTQ scan",
            "test_log":log.getvalue()}
    output=HERE/"exact_independent_controls_v2.json"
    with output.open("x") as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(log.getvalue(),end="")
    print(json.dumps({k:result[k] for k in ("status","tests_run","wall_seconds")}))
    if result["status"]!="pass":
        raise SystemExit(1)


if __name__=="__main__":
    main()

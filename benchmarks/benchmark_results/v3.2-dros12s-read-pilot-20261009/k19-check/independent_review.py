#!/usr/bin/env python3
"""Verify frozen canonical-word inventory and all per-candidate summaries."""
import csv
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
RUN=BASE/"run-v1"
COMP=str.maketrans("ACGT","TGCA")


def receipt(path):
    data=path.read_bytes()
    return {"path":str(path.resolve()),"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}


def load(path):
    return json.loads(path.read_text())


def words(sequence):
    return [min(sequence[i:i+19],sequence[i:i+19].translate(COMP)[::-1])
            for i in range(len(sequence)-18)]


def main():
    execution=load(RUN/"execution.json")
    report=load(RUN/"results/report.json")
    candidates=load(BASE.parent/"manifest.json")["candidates"]
    assert execution["status"]=="pass" and execution["returncode"]==execution["test_exit_code"]==0
    assert execution["sources_unchanged"] and execution["source_hashes_before"]==execution["source_hashes_after"]
    for name,expected in execution["source_hashes_after"].items():
        assert receipt(Path(name))["sha256"]==expected
    for name,expected in execution["outputs"].items():
        actual=receipt(RUN/name)
        assert all(actual[k]==expected[k] for k in ("bytes","sha256"))
    assert report["input_prefix"]==load(BASE.parent/"expected_prefix.json")
    assert report["manifest_sha256"]==receipt(BASE.parent/"manifest.json")["sha256"]
    assert report["k"]==19 and report["canonical"] and report["quality_filter"] is None
    with (RUN/"results/raw_counts.tsv").open() as stream:
        table=list(csv.DictReader(stream,delimiter="\t"))
    counts={row["canonical_19mer"]:int(row["raw_occurrences"]) for row in table}
    universe={word for c in candidates for word in words(c["sequence"])}
    assert len(table)==len(counts)==len(universe)==1600 and set(counts)==universe
    assert all(count>0 for count in counts.values())
    summaries={c["candidate_id"]:c for c in report["candidates"]}
    assert len(summaries)==len(candidates)==28
    emitted_words=set()
    for c in candidates:
        seq=c["sequence"]
        keys=words(seq)
        row=summaries[c["id"]]
        assert hashlib.sha256(seq.encode()).hexdigest()==row["sequence_sha256"]==c["sequence_sha256"]
        assert len(seq)==row["length"]
        assert row["valid_19mer_starts"]==row["total_possible_starts"]==row["observed_19mer_starts"]==len(keys)
        assert row["distinct_canonical_19mers"]==len(set(keys))
        assert row["minimum_raw_occurrences"]==min(counts[k] for k in keys)
        assert row["maximum_raw_occurrences"]==max(counts[k] for k in keys)
        assert row["absent_19mers"]==[]
        if row["role"]=="emitted_12S":
            emitted_words.update(keys)
            assert row["minimum_raw_occurrences"]==2
    assert report["emitted_candidates_with_every_19mer_observed"]==26
    assert report["target_canonical_19mers_absent"]==0
    assert report["target_distinct_canonical_19mers"]==1600
    result={"status":"pass","scope":"Independent string-based reconstruction of all candidate canonical-word inventories and reported count summaries; frozen source and execution review",
            "current_12s_candidates":26,"historical_context_candidates":2,
            "distinct_canonical_19mers_all_contexts":len(universe),
            "distinct_canonical_19mers_current_12s":len(emitted_words),
            "current_12s_raw_count_range":[min(counts[k] for k in emitted_words),max(counts[k] for k in emitted_words)],
            "current_12s_candidates_with_every_19mer_observed":26,
            "missing_target_words":0,
            "source_review":"Rolling forward and reverse complement encodings represent the same19-base window; invalid bases reset both states, and each read begins with empty state. Counts include each valid occurrence without quality filtering. The consumed input hash and record/base/byte totals bind the full first-million prefix.",
            "test_review":"Five passing controls include400independent string-based comparisons at k1/3/7/19, reverse complements, ambiguity reset, repeated occurrence counts and no cross-read joining.",
            "conclusion":"Every emitted Drosophila12S 19-base word is observed at least twice in the frozen prefix. This satisfies the tested observed-word invariant and does not expose a sequence-spelling failure. Combined with unobserved longer assay windows, the evidence supports short-word availability while leaving longer combinations uncorroborated.",
            "limits":["This review verifies the frozen scan algorithm, input receipt and reported count-table aggregation; it does not perform another full-prefix scan.",
                      "Raw canonical-word counts cannot validate longer adjacency, target identity, repeat copy number, independent molecules or whole products.",
                      "A spelling error that recombines observed words would not necessarily fail this invariant; no claim is made that all possible assembly bugs are excluded."],
            "artifacts":{str(p.relative_to(BASE)):receipt(p) for p in
                         (BASE/"count_k19.py",BASE/"test_count_k19.py",BASE/"run_check.py",RUN/"execution.json",
                          RUN/"results/report.json",RUN/"results/raw_counts.tsv",RUN/"results/per_candidate.tsv",
                          RUN/"tests.log",RUN/"time.json",Path(__file__))}}
    with (BASE/"independent_review.json").open("x") as stream:
        stream.write(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","distinct_canonical_19mers_current_12s","current_12s_raw_count_range","current_12s_candidates_with_every_19mer_observed")}))


if __name__=="__main__":
    main()

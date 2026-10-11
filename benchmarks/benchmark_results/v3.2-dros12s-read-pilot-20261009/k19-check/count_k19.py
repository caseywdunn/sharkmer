#!/usr/bin/env python3
"""Count selected canonical 19-mers in one frozen R1 prefix, without Q filtering."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PILOT=HERE.parent
K=19
EXPECTED_MANIFEST_SHA256="b95c3ba4e6503ce991c40f94b4e191e7da2b8c859a41fddffe722bbcc96b53e7"
VALUES={65:0,67:1,71:2,84:3}


def kmers(sequence,k=K):
    """Emit (start, min(forward,reverse-complement)) for each valid start."""
    forward=reverse=valid=0
    mask=(1<<(2*k))-1
    high=2*(k-1)
    for index,base in enumerate(sequence):
        value=VALUES.get(base)
        if value is None:
            forward=reverse=valid=0
            continue
        forward=((forward<<2)|value)&mask
        reverse=(reverse>>2)|((3-value)<<high)
        valid+=1
        if valid>=k:
            yield index-k+1,min(forward,reverse)


def decode(value,k=K):
    return ''.join("ACGT"[(value>>(2*i))&3] for i in range(k-1,-1,-1))


def load_candidates():
    data=(PILOT/"manifest.json").read_bytes()
    if hashlib.sha256(data).hexdigest()!=EXPECTED_MANIFEST_SHA256:
        raise ValueError("Frozen Drosophila candidate manifest changed")
    candidates=json.loads(data)["candidates"]
    for candidate in candidates:
        if hashlib.sha256(candidate["sequence"].encode()).hexdigest()!=candidate["sequence_sha256"]:
            raise ValueError("Candidate sequence digest mismatch")
    if sum(not c["id"].startswith("SRR") for c in candidates)!=26:
        raise ValueError("Expected 26 emitted Drosophila 12S strings")
    return candidates


def scan(path,targets,expected):
    counts={key:0 for key in targets}
    digest=hashlib.sha256()
    records=bases=size=valid_starts=0
    with path.open('rb') as stream:
        for record_index in range(expected["records"]):
            lines=[stream.readline() for _ in range(4)]
            if any(not line for line in lines):
                raise ValueError("FASTQ prefix truncated")
            for line in lines:
                digest.update(line);size+=len(line)
            header,sequence,plus,quality=[line.rstrip(b'\r\n') for line in lines]
            if not header.startswith(b'@') or not plus.startswith(b'+') or len(sequence)!=len(quality):
                raise ValueError("Malformed FASTQ")
            records+=1;bases+=len(sequence)
            for _,key in kmers(sequence):
                valid_starts+=1
                if key in counts:
                    counts[key]+=1
        if stream.read(1):
            raise ValueError("Input extends beyond expected frozen prefix")
    observed={"records":records,"bases":bases,"size_bytes":size,"sha256":digest.hexdigest()}
    if observed!=expected:
        raise ValueError(f"Input prefix attestation failed: {observed}")
    return counts,observed,valid_starts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastq',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir()
    candidates=load_candidates()
    expected=json.loads((PILOT/'expected_prefix.json').read_text())
    targets={key for c in candidates for _,key in kmers(c['sequence'].encode())}
    counts,observed,valid_starts=scan(args.fastq,targets,expected)
    rows=[]
    for c in candidates:
        starts=list(kmers(c['sequence'].encode()))
        absent=[{"start":start,"end":start+K,"sequence":c['sequence'][start:start+K],
                 "canonical_19mer":decode(key)} for start,key in starts if counts[key]==0]
        rows.append({"candidate_id":c['id'],"sequence_sha256":c['sequence_sha256'],
                     "role":"emitted_12S" if not c['id'].startswith('SRR') else "historical_control",
                     "length":len(c['sequence']),"total_possible_starts":len(c['sequence'])-K+1,
                     "valid_19mer_starts":len(starts),"observed_19mer_starts":sum(counts[key]>0 for _,key in starts),
                     "distinct_canonical_19mers":len({key for _,key in starts}),
                     "minimum_raw_occurrences":min(counts[key] for _,key in starts),
                     "maximum_raw_occurrences":max(counts[key] for _,key in starts),
                     "absent_19mers":absent})
    report={"status":"completed","k":K,"canonical":True,"quality_filter":None,
            "count_unit":"all occurrences across raw R1 records; repeated positions within one read each count",
            "ambiguous_bases":"any non-uppercase-ACGT byte resets both orientations",
            "coordinates":"zero_based_half_open_in_emitted_sequence",
            "input_prefix":observed,"valid_read_kmer_starts":valid_starts,
            "target_distinct_canonical_19mers":len(targets),
            "target_canonical_19mers_absent":sum(n==0 for n in counts.values()),
            "emitted_candidates":26,"emitted_candidates_with_every_19mer_observed":sum(
                not r['absent_19mers'] for r in rows if r['role']=='emitted_12S'),
            "candidates":rows,"manifest_sha256":EXPECTED_MANIFEST_SHA256,
            "limitations":["Only the frozen first-million R1 prefix was checked; counts are not full-library graph counts.",
                           "A canonical 19-mer occurrence supports a short string in either orientation, not adjacency beyond 19 bases or whole-product phase.",
                           "No quality filter is used; these are not validated molecules or biological sequence truth.",
                           "The 26 emitted 12S strings and two historical controls are labelled separately."]}
    with (args.output/'report.json').open('x') as stream:stream.write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    with (args.output/'raw_counts.tsv').open('x') as stream:
        writer=csv.writer(stream,delimiter='\t');writer.writerow(['canonical_19mer','raw_occurrences'])
        writer.writerows((decode(key),counts[key]) for key in sorted(counts))
    with (args.output/'per_candidate.tsv').open('x') as stream:
        fields=[k for k in rows[0] if k!='absent_19mers']+['absent_19mer_start_count']
        writer=csv.DictWriter(stream,fields,delimiter='\t');writer.writeheader()
        for row in rows:writer.writerow({**{k:v for k,v in row.items() if k!='absent_19mers'},'absent_19mer_start_count':len(row['absent_19mers'])})
    print(json.dumps({key:report[key] for key in ['emitted_candidates_with_every_19mer_observed','target_distinct_canonical_19mers','target_canonical_19mers_absent']}))


if __name__=='__main__':main()

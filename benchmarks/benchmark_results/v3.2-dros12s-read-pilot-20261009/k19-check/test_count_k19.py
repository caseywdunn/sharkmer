import collections
import random
import unittest
import count_k19 as audit


def brute(sequence,k):
    result=[]
    for start in range(len(sequence)-k+1):
        word=sequence[start:start+k]
        if set(word)<=set('ACGT'):
            rc=word.translate(str.maketrans('ACGT','TGCA'))[::-1]
            result.append((start,min(word,rc)))
    return result


class CountingControls(unittest.TestCase):
    def test_randomized_brute_force(self):
        generator=random.Random(29403)
        for k in [1,3,7,19]:
            for _ in range(100):
                sequence=''.join(generator.choices('ACGTN',weights=[10,10,10,10,1],k=150))
                self.assertEqual([(start,audit.decode(value,k)) for start,value in audit.kmers(sequence.encode(),k)],brute(sequence,k))

    def test_reverse_complements_share_counts(self):
        sequence='ACAGATCTGCGAGTGTATTGGGCATAGATTA'
        rc=sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]
        self.assertEqual(collections.Counter(v for _,v in audit.kmers(sequence.encode())),collections.Counter(v for _,v in audit.kmers(rc.encode())))

    def test_ambiguous_base_resets(self):
        sequence='A'*18+'N'+'C'*19
        self.assertEqual([(p,audit.decode(v)) for p,v in audit.kmers(sequence.encode())],[(19,'C'*19)])

    def test_repeated_occurrences_are_counted(self):
        counts=collections.Counter(v for _,v in audit.kmers(b'A'*21))
        self.assertEqual(list(counts.values()),[3])

    def test_no_cross_read_bridge(self):
        self.assertEqual(list(audit.kmers(b'A'*18))+list(audit.kmers(b'A'*18)),[])


if __name__=='__main__':unittest.main()

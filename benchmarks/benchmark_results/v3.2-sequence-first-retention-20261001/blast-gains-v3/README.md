# Local BLAST review of sequence-retention gains

This post-comparison check uses the 200 extra uncertain identities in
`../report-v3/uncertain_gains.tsv`. They came from the corrected pair-1
candidate FASTAs. `../check_gains_blast_sensitive.py` reconstructs the queries
from those FASTAs and checks their hashes, lengths, and uncertainty headers.
The queries and the 121 public-source sequences are saved here. The source
catalog is the repository's pinned `panels/reference_sources.json.gz`, not a
live or exhaustive search of public nucleotide records. The two panel
reference sets contain 19 angiosperm and 42 insect regions.

The producer used BLASTN 2.17.0+ with `-task blastn -dust no -evalue 1e-10
-max_target_seqs 121 -max_hsps 20`. `summary.json` records source, query,
and hit-file checksums. Raw tabular HSPs are in `public_sources_hits.tsv` and
the two `*_reference_hits.tsv` files. `per_sequence.json` maps query IDs to
the gains and selected alignments. The `*_db.*` files are generated BLAST
indexes; the FASTAs and saved tabular HSPs are the portable evidence.

`../review_gains_blast.py` independently checks gain/query sequence hashes,
catalog/source sequences, recorded checksums, every saved HSP's IDs and
coordinates, hit counts, direct same-taxon exact-substring matches, and the
90% query-coverage/90% identity single-HSP criterion. Its output is
`review.json` and passes. This audits the saved search output; it does not
repeat the BLAST execution. `../BLAST_SHA256SUMS` lists SHA-256 digests for
the producer, review, and preserved result files.

| Search set | Sources | HSPs | Queries with a hit | HSPs meeting both 90% cutoffs |
| --- | ---: | ---: | ---: | ---: |
| Pinned public sources | 121 | 5,289 | 168/200 | 0 |
| Panel references | 61 | 1,103 | 88/200 | 0 |

Five public-source HSPs cover at least 90% of a query. All five are Gryllus
18S queries aligned to a Drosophila source at 64.0–64.5% identity. Seventy-eight
public-source HSPs reach at least 90% identity only over shorter segments.
No query is an exact full-length forward or reverse-complement substring of a
same-taxon catalog source. The source catalog contains 2 Liriodendron, 10
Drosophila, 5 Heliconius, and 6 Gryllus records; not every target has an
independently annotated region or same-sample truth.

**Interpretation:** none of these 200 identities has qualifying full-query
support in this bounded local reference set under this search. A miss does
not prove an incorrect sequence or absence from a larger public database.
Short, split, or cross-taxon hits do not establish a correct repeat structure,
haplotype, or locus. Biological and full-product read support remain
unevaluated. The preceding `blast-gains-v1`/`v2` panel-only classifications
are exploratory and are not the conclusion of this review.

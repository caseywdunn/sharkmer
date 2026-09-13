# Public reference provenance and population variation

Date: 2026-09-13. Implementation/evaluation issue: [#129](https://github.com/caseywdunn/sharkmer/issues/129).
Release review remains open under [#153](https://github.com/caseywdunn/sharkmer/issues/153).

## Two different questions

1. **Where did a reference come from?** Its stored sequence must be an exact,
   reproducible extraction from a pinned public nucleotide record. A BLAST-hit
   accession attached to a Sharkmer assembly is not that provenance.
2. **Does a product have support for the target gene?** A different individual
   can carry SNPs, indels, repeat-length differences, or another haplotype.
   Exact equality to an archive sequence is not required. Alignment support
   does not establish the sample's full-length haplotype or repeat-copy number.

Strictness in question 1 must not become a zero-variation requirement in
question 2. Conversely, allowing variation cannot authenticate the origin of
an alleged reference. Unverified source provenance is not evidence that a
historical product was biologically wrong.

## Complete panel inventory

The audit covers all nine built-ins and the documentation example: **159 old
entries**, 111 requested accession strings resolving to **109 versioned public
records** containing 2,123,261 unique-record bases. The raw FASTA transfer has
111 records and 2,143,980 bases because two records occur twice; identical
duplicates are consolidated, not counted as independent evidence. Records were retrieved independently
from NCBI on 2026-09-13; the raw FASTA and metadata are preserved in the
[audit evidence archive](../benchmarks/benchmark_results/v3.2-reference-provenance-20260913/README.md).

| Panel | Old entries | Exact source candidates | No exact source match |
| --- | ---: | ---: | ---: |
| angiospermae | 21 | 12 | 9 |
| bacteria | 0 | 0 | 0 |
| c_elegans | 16 | 13 | 3 |
| cnidaria | 31 | 10 | 21 |
| human | 9 | 9 | 0 |
| hydrozoa | 29 | 28 | 1 |
| insecta | 41 | 17 | 24 |
| metazoa | 0 | 0 | 0 |
| teleostei | 11 | 10 | 1 |
| documentation example | 1 | 1 | 0 |
| **Total** | **159** | **100** | **59** |

Every retained sequence is **rebuilt from the downloaded source**, not copied
from a Sharkmer output. The 100 source-exact candidates comprise 20 whole-record matches,
78 contiguous subregions, and two human mitochondrial origin-spanning regions.
Literal IUPAC symbols are preserved, never wildcard-expanded or imputed.
The two c_elegans 18S references each match two locations in Z92784; both
placements are recorded, one is selected deterministically, and no unique
locus/copy identity is claimed.

The 59 source-quarantined entries have no exact literal extraction from their named
public records. This category can include assembled alleles, trimming or
annotation errors, and other unknown origins; it does **not** mean 59 known
incorrect sequences or 59 proven Sharkmer products. All 159 original entries,
including their sequences and old labels, remain in a regression-only archive.
The Gryllus 978 bp ITS_2 entry is a known Sharkmer bootstrap example, not the
256 bp [AK281180.1 record](https://www.ncbi.nlm.nih.gov/nuccore/AK281180.1).

Thirty-one of the 100 source-exact candidates had organism labels differing
from the public record; approved migrations use the source organism instead.
Validation-sample taxa are unchanged. Exact text taxon matching
does not resolve taxonomy synonyms or establish biological species identity.
Panel patch versions advance only for the eight files with reference changes;
primers, validation samples, assembly settings, and default **k=19** do not change.

**Additional annotation quarantine:** the two 16S-labeled XM_047005986.1 entries
(cnidaria and documentation example) match the public source exactly, but its
title is a predicted Xenia mutS mRNA. Both are excluded pending annotation
review, rather than being certified as 16S merely because their bytes match.
The catalog records this reviewed conflict, so reintroducing the same accession
under 16S cannot silently restore positive gene support. The documentation
example instead uses the existing independently sourced Hydra 16S region from
hydrozoa. Final active references: **99 = 98 retained old regions + one new
example**. Final historical quarantine: **61 = 59 source mismatches + two
source-exact annotation conflicts**. All decisions are archived separately.

Source provenance does not automatically establish the remaining gene
assignments; these remain curator-supplied annotations, with broader annotation
review still open. Public assemblies/transcripts can also contain errors.
The references are independent evidence, not absolute sample truth or a newly
held-out evaluation cohort.

## Reproducible reference contract

`panels/reference_sources.json.gz` is a pinned, checksummed local snapshot.
Each record contains the versioned accession, literal sequence and SHA-256,
organism/taxid, topology, source URL, and retrieval time. Checksums detect
inconsistency; they do not authenticate a maliciously replaced catalog.
Catalog updates are reviewed source acquisitions, not automatic network
refreshes during validation.
Optional record-level `annotation_conflicts` maps logical gene names to reviewed
exclusion reasons; its contents are curated decisions, not raw NCBI fields.

Each active reference records `provenance` with schema version 1, kind
`public_record_region`, `accession_version`, `source_sequence_sha256`,
`source_length`, `start`, `end`, `strand`, `wraps_origin`, and `sequence_sha256`.
Coordinates are **zero-based, half-open** on the source orientation:

- Normally extract `source[start:end]`.
- For circular records, `wraps_origin: true` extracts
  `source[start:] + source[:end]`, spanning at most one source traversal.
- Reverse-complement the extracted region for `strand: "-"`.
- Require literal equality, matching hashes/lengths, and the source organism.

The updated Rust loader accepts optional metadata so legacy panels still assemble.
Older releases that reject unknown fields cannot load new provenance-bearing
panel YAML; use their bundled/pinned older panels or upgrade the loader.
Rust schema validation is **not** public-record verification. Python biological
validation excludes missing, malformed, or inconsistent provenance and never
falls back to an unverified embedded sequence. Reports/results retain the
catalog receipt and excluded-reference reasons, including zero-reference runs.

Offline verification (PyYAML required; no network or BLAST needed):

```bash
python scripts/audit_panel_references.py --verify-existing
```

For external trusted catalogs, pass `--reference-catalog path/catalog.json.gz`
to the audit, `scripts/validate_panel.py`, or `benchmarks/run_benchmark.py`.
The default is the repository catalog, not an automatically trusted file
downloaded alongside an external panel. Both compressed/decompressed catalogs
are capped at 25 MiB; individual source sequences at 5 MiB in this first
bounded implementation. Larger source requirements need an explicit design
change rather than silent truncation.

To propose a source migration, independently acquire NCBI ESummary metadata
as `metadata.json` and EFetch nucleotide FASTA as `records.fasta` in an archive
directory, pin their receipts, then run:

```bash
python scripts/audit_panel_references.py \
  --panels-root panels --archive-dir /path/to/source-archive \
  --output /tmp/reference-audit.json --gzip-output /tmp/reference-catalog.json.gz \
  --migration-output /tmp/reference-migration.json \
  --source-url 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text' \
  --retrieved-at 'YYYY-MM-DDTHH:MM:SSZ'
```

This produces a proposal, not panel edits. Review versions, coordinates,
organisms, gene annotation, nonunique placements, and exclusions; rebuild
approved regions from the catalog, bump panel versions, and run the offline
check before committing both panels and source evidence. To reproduce this
migration rather than audit the new panels, use the archive's frozen panels.
Apply the separately archived annotation decisions after the source-only audit;
an automatically generated catalog must not discard reviewed exclusions.

`bootstrap_from_runs.py` labels its products `regression_only`,
`source_kind: sharkmer_assembly`, and `independent_reference: false`.
`bootstrap_references.py` emits `reference_candidates`, not a publishable
`references` block; direct `--write` publication is disabled. Candidate
discovery and a partial BLAST match do not bypass provenance review.

## Variation-aware product evidence

The validator retains default 90% alignment identity and 90% query coverage
as recorded analysis criteria, **not universal species/gene cutoffs**. It
reports alignment extent on both query and reference, differences, gaps,
unmatched regions, competing hits, and the selected reference provenance.
Compatible collinear HSPs can contribute support without being called a
chimera solely because BLAST splits the alignment. Stronger conflicting gene
evidence, competing ties, and contradictory arrangements remain distinct.

Evidence is separated into:

- `target_support`: supported, conflicting, ambiguous, insufficient, or unavailable.
- `sequence_relationship`: reference-identical, aligned differences, partial
  unresolved, structural conflict, or unavailable.
- `haplotype_truth: not_established`: not supplied by cross-individual alignment.
- `read_support: not_evaluated`: this validator does not independently realign reads.

New same-gene results use `gene_supported_expected_taxon` or
`gene_supported_other_taxon`, not `confirmed_product`. Gene support can
coexist with SNPs/indels; even exact query equality to a reference subregion
does not mean the product covers the whole public record. Missing/inadequate
reference coverage leaves uncertainty, not a false-negative correctness verdict.
Legacy result files are historical artifacts; their old `confirmed_product`
labels must not be promoted into new biological truth claims.

## Reclassification and unchanged assembly checks

The saved 2026-09-11 high-copy comparison was reclassified without rerunning
Sharkmer or altering its original results. All 60 result/FASTA receipts verify;
the three replicates in each of 20 role/cell groups have identical sequence
fingerprints, permitting 150 representative product assessments. Initial BLAST+
2.17.0 executable receipts, the final catalog, and stable database metadata
receipts are recorded. Temporary database artifacts are checksum-verified
during queries, but the indexes and their build-specific hashes are not
preserved by this reclassification; it does not attest the executable again
afterwards.
This uses the new panel **regions**, not a database of whole source genomes.

There are **30 gene-supported product identities in each version, and those
identity sets are equal**: 21 expected-taxon-supported and nine other-taxon-
supported per version. This is reassuring evidence within the covered targets,
not proof that every missing product was wrong. Inadequate references can hide
real recovery differences. The seven previously identified high-copy losses
remain baseline-only sequences:

| Saved product | New reference assessment | What this establishes |
| --- | --- | --- |
| Gryllus ITS_2, 978 bp | `no_verified_reference` | The old bootstrap match no longer supplies circular confirmation; correctness remains unresolved |
| Gryllus CO1_1, 373 bp | 100% identity over 81.501% of query | Partial support, not full-product rejection or confirmation |
| Retained Gryllus CO1_1, 352 bp | 100% identity over 86.364% of query | Also below the 90% query-coverage gate with this short reference |
| Drosophila 12S, 475/487/499 bp | `no_verified_reference` | No independent eligible panel reference; not evidence of incorrect sequence |
| Drosophila 16S_2, 590 bp; Heliconius ND1, 263 bp | `no_significant_hit` | No qualifying local hit in this limited database; not proof of absence/error |

The CO1 comparison currently uses a **304 bp public region** of PP230540.1.
Both the lost 373 bp and retained 352 bp sequences therefore receive
`insufficient_alignment`; this is limited reference span, not an exact-match
requirement or evidence against population variation. The earlier
[whole-public-record assessment](BENCHMARK_reference_gates.md#existing-seven-high-copy-losses)
contains additional evidence beyond that short region. Extend independently
annotated gene coverage and assess read support before deciding whether the
21 bp difference is a genuine allele or assembly error.

A separate before/after implementation smoke test uses the same 100,000-read
Porites fixture at k19 and k31. All five ingestion-count fields and all product
sequence hashes match: five products at k19 and two at k31. Parsed primers and
validation blocks are unchanged across the migrated panels. This protects
against accidental assembly changes in this patch; it is **not** a new timing,
RAM, or held-out biological benchmark. Validation passes 223 Rust unit tests,
22 integration tests, 77 Python tests, format checks, and Clippy on both hash
backends. Astra independently reviewed the classifier, public-source audit,
known annotation exclusions, and migration invariance.

## Remaining release decisions

This change improves the **evaluation contract**, not assembly recovery or
counting throughput. It neither relaxes nor strengthens repeat gates. Historical
product losses remain real output regressions; their biological correctness
requires evidence beyond exact matching to bootstrap products.

Next work under #129/#153 is to review gene annotations, add independent
same-sample truth where available, and evaluate read/pair support for differing
and repeat-spanning regions. Reads from the same dataset can support a
candidate but are not independent truth; reads unable to span a repeat or
phase variants leave the full haplotype unresolved. Register held-out inputs,
callable regions, and locus-appropriate SNP/indel tolerances before using them
to tune assembly gates. No release is authorized by this audit.

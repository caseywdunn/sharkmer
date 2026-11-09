# Accessioned-reference assessment of repeat gates

Date: 2026-09-13. Release gate: [#153](https://github.com/caseywdunn/sharkmer/issues/153).
This is a sequence-correctness assessment, not a performance benchmark or a
release approval. Production code, panel defaults, and default k remain unchanged.

## Question and decision criterion

Does repeat withholding remove sequences that the preceding code reconstructed
correctly? There are two distinct oracles: exact recovery of a known generating
template in a synthetic experiment, and agreement with an independently
deposited biological sequence. A match to a previous Sharkmer assembly stored
in a panel is a regression oracle, but not independent biological confirmation.

**Important correction:** the earlier interpretation of Gryllus ITS_2 as a
full-length accession-confirmed product was too strong. The panel's 978 bp
sequence labeled AK281180 is a bootstrapped Sharkmer amplicon. The actual
public AK281180.1 record is only 256 bp. This is a reference-provenance problem,
not merely a hypothetical read-phasing limitation.

The causal comparison is the parent of #132,
`928491610822ffdb236d0326b32c8dd19fecf66a`, against its immediate successor,
`b728f14219945567b3d78f3602590864a256c809`. The current production implementation
is `0c38d6a051082cc2a1eb961b4206837f82948fda`; later development commits only
record evidence. Exact known-template recovery before #132 and loss immediately
afterward on identical synthetic inputs is a recovery regression. A real-data
loss requires independent sequence evidence before calling it biologically
incorrect withholding. Missing output in both versions is not evidence that
the new gate caused a loss.

The initial oracles are sequences actually embedded in the nine built-in
panels, with accession labels, taxon, panel version, and sequence/file hashes
preserved. A separate public-record audit checks those labels rather than
silently replacing frozen inputs. Existing panel references and real samples
remain calibration data, not independent biological held-out validation.

## Reference provenance discovery

The [public AK281180.1 record](https://www.ncbi.nlm.nih.gov/nuccore/AK281180.1)
retrieved on 2026-09-13 is a 256 bp Gryllus bimaculatus mRNA record,
GBcontig28041. Its reverse complement exactly matches positions 573–828
(1-based, inclusive) of the panel's 978 bp sequence. It does not validate the
other 722 bases or the full reconstructed product.

Git commit `aafe77aef0a87ec4232c5d8463622a3f7fd55cb2` explicitly populated
panel references from NCBI-verified **bootstrap amplicons**. The bootstrap code
collects Sharkmer output sequences, submits them to BLAST, and writes the
amplicon sequence with hit metadata. A hit accession is therefore not necessarily
the source of the entire stored sequence. Exact agreement of a later Sharkmer
product with that stored amplicon can be circular evidence. The earlier
benchmark's frozen `confirmed_product` classification is reproducible, but
does not acquire independent full-length biological truth through repetition.

This provenance distinction must be addressed under #129 as well as the
high-copy release review. Preserve existing regression fixtures, but label
bootstrap-derived sequences separately from exact deposited regions, including
source version, coordinates, strand, and alignment extent where applicable.

A separate public NCBI batch retrieves all 29 hydrozoa accessions (53,994
deposited bases), before the synthetic experiment produces outcomes. Of the
23 ACGT-only embedded templates, **22 are exact contiguous deposited regions**
in the stored or reverse-complement orientation. The exception is
[EU293971.1](https://www.ncbi.nlm.nih.gov/nuccore/EU293971.1): the panel has
538 bp versus the public record's 539 bp, with one missing base within public
positions 279–281 (1-based; repeated bases make the deletion position ambiguous).
This template remains in the frozen diagnostic cohort,
but is not counted as exact deposited-sequence truth. Six additional embedded
references contain ambiguity and remain outside the unambiguous fixture cohort.
The public records, coordinates, hashes, and eligibility distinctions are
preserved separately; frozen panel content is not silently repaired.

## Existing seven high-copy losses

An offline audit inventories 158 accession-labeled references across nine panels and
compares the seven previously lost raw sequence identities with all references
for the same logical gene, including other panels and primer indices. It verifies
the old FASTAs against the three preserved baseline replicas.

| Old product | Reference evidence | Interpretation |
| --- | --- | --- |
| Gryllus ITS_2, 978 bp | Exact equality to the panel's bootstrapped sequence labeled AK281180; public AK281180.1 supports a 256 bp reverse-complemented internal segment | Embedded-sequence regression confirmed; full biological correctness remains unverified |
| Gryllus CO1_1, 373 bp | Contains the 304 bp panel reference; public PP230540.1 exactly matches a 328 bp prefix | Partial external support, not full correctness for the remaining 45 bp or whole product |
| Drosophila 12S, 475/487/499 bp; 16S_2, 590 bp; Heliconius ND1, 263 bp | No full exact or exact-subsequence match in the same-gene embedded reference collection | Unresolved by this oracle, not proven wrong |

Logical gene matching here does not conflate primer index with gene identity.
The comparison is exact sequence matching, not a new BLAST analysis. Lack of
exact agreement cannot distinguish a genuine divergent sequence from an error.
The [public PP230540.1 mitogenome](https://www.ncbi.nlm.nih.gov/nuccore/PP230540.1)
is 15,955 bp; neither the complete lost CO1_1 product nor its reverse complement
is an exact contiguous substring. Its 328 bp matching prefix begins at position
1669 (1-based) of that record. This does not prove the unmatched tail erroneous.

A subsequent exact-block comparison finds that the retained post-#132/current
352 bp CO1_1 product matches 351 contiguous public bases at that same position.
The lost 373 bp product equals 328 shared bases, a 21 bp insertion, then 24
shared bases relative to the retained product. The inserted sequence occurs
nowhere in the retained product or public mitogenome. The old product resumes
with 23 exact public bases after that insertion. This favors the retained
sequence as a closer reference match, but does not rule out a genuine insertion
in the sequenced sample or independently certify either entire product.

## Real Gryllus causal contrast

Three fresh isolated invocations use the same first 1,000,000 records from
SRR27962769, the same whole insecta panel, k=19, two threads, `--chunks 0`, and
threading disabled. Only the binary and output directory change. This is one
correctness run per version, not a timed performance comparison.

| Version | Exact 978 bp panel-embedded ITS_2 sequence |
| --- | --- |
| Immediately before #132 | Recovered |
| Immediately after #132 | Absent |
| Current development implementation | Absent |

All five aggregate ingestion/count fields agree: 1,000,000 reads,
150,000,000 bases read, 1,000,000 subreads, 149,999,387 bases ingested, and
131,993,167 counted k-mer occurrences. Thus this loss is not explained by
different data ingestion or a different k.

Before #132, ITS_2 returns the exact embedded sequence at threshold 4. Immediately after
#132, the log explicitly withholds all 20 enumerated candidates at threshold 4;
the later threshold-2 node-budget failure obscures this earlier cause in the
historical final summary. Current diagnostics retain the distinction: at
threshold 4 all 180 completed candidates touch pre-pruning SCC markers and are
withheld, with no DFS or eligible-path quota exhaustion; threshold 2 then fails
to establish connectivity within the node budget.

This isolates removal of that real-data sequence by #132, not just an
association across an entire release's changes. It does **not** independently
establish that the whole 978 bp sequence was biologically correct, because its
panel reference came from an earlier assembly. Nor does it establish correctness
of the other six missing sequence identities.

## Reference-derived fixture protocol

The primary inventory contains all 29 hydrozoa references plus Gryllus AK281180.
It checks the actual retained 3-prime primer suffix, configured mismatches,
orientation, and length bounds. Full-primer matching alone is not an eligibility
test because production trims primers. Ambiguous reference bases are never
silently resolved.

The panel sequence labeled AK281180 is eligible for a known-template synthetic
oracle at k=19 and k=31 under unchanged panel parameters. None of the hydrozoa references supplies an eligible
unambiguous oracle under those parameters: six contain IUPAC ambiguity, and
the remaining 23 lack a permitted primer pair or permitted product length.
Two ambiguous 28S references do contain suitable primer pairs; their exclusion
does not mean their primer sites are absent. These are fixture-eligibility
limitations, not observed sPCR failures on biological samples.

Therefore a separately specified **native-endpoint diagnostic** uses all 23
ACGT-only hydrozoa embedded references without changing any reference bases. Forward
primers are the first 15 bases; reverse primers are the reverse complement of
the last 15 bases. Trim is 15, mismatches are zero, minimum length is 100 bp,
and maximum length is reference length plus 500 bp. Other search defaults stay
unchanged. This deliberately changes primers/bounds to isolate reconstruction
of the embedded sequence; it is **not recovery with the production hydrozoa
panel**. Internal matches to those endpoint primers are recorded rather than
silently excluded. Broad bounds allow many shortened erroneous products to
remain observable; they cannot admit every possible erroneous length.

For each eligible reference, deterministic 149 bp ACGT flanks are added outside
the unchanged template. Every 150 bp window is emitted forward and
reverse-complemented, with constant quality. Every reference k-mer occurrence
therefore receives 264 observations at k=19 or 240 at k=31; repeated identities
legitimately aggregate. The same FASTQ is used across versions and k values.
Reference-containing reads have no sequencing errors or coverage gaps. Flanks
must not manufacture accepted primer-bounded truth. All emitted products are
checked, separating full-reference exactness, exact internal segments, and
sequences not present contiguously in the generating template. Clean-positive
and known collapsed-repeat controls are included.

## Controlled results

All **150 synthetic invocations** complete successfully. All 50 fixture/k
groups have identical values for the five ingestion/count fields across the
three versions. The 29 unchanged-panel hydrozoa eligibility assessments are
not CLI runs and are not counted as failed assemblies. The 13 harness tests
pass; source, binary, fixture, and protocol receipts remain unchanged during
execution.

| Native-endpoint hydrozoa diagnostic | Pre-#132 exact full templates | Post-#132 exact full templates | Current exact full templates |
| --- | --- | --- | --- |
| k=19 | 23/23 | 23/23 | 23/23 |
| k=31 | 23/23 | 23/23 | 23/23 |
| Public-record-corroborated subset, either k | 22/22 | 22/22 | 22/22 |

Every diagnostic invocation emits exactly one product, identical to its
generating template. Thus the 22 independently corroborated deposited regions
show **no over-conservative rejection in these clean, uniformly covered
fixtures**. The one-base-deleted EU293971 panel derivative also reconstructs
exactly, but that is not recovery of the unchanged deposited record.

The 978 bp Gryllus bootstrap template likewise reconstructs exactly in all
three versions at both k values when used to generate clean reads. Its failure
on the real input is therefore not an unavoidable property of this template
alone; actual coverage, graph context, and read composition must be considered.
This experiment does not identify which of those differences is decisive.

Both repeat controls run at k=19. A clean 138 bp template containing A18 is
recovered exactly in all versions. A 160 bp template containing A40 produces
a **wrong, collapsed 138 bp product before #132**, while post-#132 and current
withhold it. No other emitted synthetic product is noncontiguous with its
generating template, and there are no internal-substring or flank-dependent
outputs. This preserves positive evidence that the gates prevent a genuine
error while retaining all tested correct hydrozoan reconstructions.

## What the mechanism does and does not establish

The current gate marks every oriented node in a cyclic strongly connected
component before graph pruning. It rejects a completed candidate if that path
touches even one marked node. The candidate need not traverse the cycle or
revisit a node. An off-route cycle can therefore veto an otherwise correct path,
and pruning the cycle does not clear its marker.

The 978 bp panel sequence labeled AK281180 has 961 distinct oriented 18-mers at k=19, with no
repeats. It also has no 18-mer shared with its reverse complement; both properties
hold for 30-mers at k=31. This is a reference-only motif check, not a trace of
the actual read-derived graph. The real before/after contrast establishes the
loss; it does not identify its precise SCC intersection, every enumerated path,
or independent correctness of the entire sequence.

Replacing the gate with “reject only paths that revisit nodes” is not a safe
fix: a shortcut through a collapsed tandem repeat can itself be a nonrevisiting
path. Exact known-template recovery and prevention of known collapsed products
must both remain acceptance criteria. Neither circular panel agreement nor
the existence of an SCC alone settles the biological correctness of the missing
real-data product. Any policy revision needs candidate-specific evidence and
valid independent or constructed truth; reference matching must not become an
assembly filter that supplies the desired answer.

## Conclusion and next action

**These tests do not demonstrate that the repeat checks are too conservative.**
They retain every tested correct deposited-region reconstruction, while
rejecting a known collapsed product. They also do not establish that all lost
real-data products were wrong. The strongest earlier apparent biological
counterexample, Gryllus ITS_2, used a circular panel reference: the causal
output loss is real, but full biological correctness is still unresolved.

Before changing policy:

1. Audit reference provenance and label bootstrap amplicons versus independently
   deposited regions; do not treat a partial BLAST-hit accession as provenance
   for every base of the query. Record the hydrozoa one-base discrepancy for
   review rather than silently altering a frozen benchmark.
2. Preserve these exact-reference positive controls and the collapsed-repeat
   negative control as acceptance criteria for any repeat-policy revision.
3. Trace the lost real candidates through their actual threshold graphs and
   seek independent full-region sequence or candidate-specific spanning-read
   evidence. Distinguish a correct route touching an off-route cycle from an
   unsupported repeat-copy shortcut.
4. Only then assess a narrower gate. Do not disable all SCC safeguards, claim
   all seven losses are accuracy improvements, or count agreement with a prior
   assembly as independent validation.

#129 and #153 remain open. No production change, new default k, release, tag,
master merge, or package publication is authorized by this assessment. Stop
for user review before revising assembly policy or releasing.

## Evidence

The [reference-gate archive](../benchmarks/benchmark_results/v3.2-reference-gates-20260913)
preserves the 150 synthetic and three real-data invocations, all raw outputs
and inputs for synthetic cases, public-record responses and provenance audits,
frozen protocols/helpers/build receipts, and per-file checksums. Original
large real FASTQs and compiled binaries are not duplicated in Git. Sol
implements the synthetic harness; Terra independently audits references and
real-run receipts; Astra reviews the protocol and results. No runtime/RAM
comparison is inferred from these correctness runs.

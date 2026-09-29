# Product behavior

The user prefers **simple, inclusive behavior by default**. Out-of-the-box
operation should deliver as much useful sequence information as the bounded
analysis produces, with clear annotations that let users choose what to keep.

- Include standard products and repeat-uncertain hypotheses together in the
  normal per-target output FASTA. Do not require an opt-in flag or a separate
  uncertainty FASTA to obtain those hypotheses.
- Use explicit, stable FASTA header fields to describe assembly status and
  uncertainty. Make filtering with simple shell commands straightforward,
  and provide an explicit opt-out for users who want only standard products.
- Keep delivery simple while preserving provenance and the distinction between
  an emitted hypothesis and a biologically validated sequence. More emitted
  records do not establish more validated loci, haplotypes, or correct repeats.
- Keep search and retention bounded, and report when those limits are reached.
  Inclusive delivery does not imply exhaustive enumeration.

Apply this preference when making future interface and output decisions.

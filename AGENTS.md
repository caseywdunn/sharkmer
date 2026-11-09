# Product behavior

The user prefers **simple, inclusive behavior by default**. Out-of-the-box
operation should deliver as much useful sequence information as the bounded
analysis produces, with clear annotations that let users choose what to keep.

- Include standard products and repeat-uncertain hypotheses together in the
  normal per-target output FASTA. Do not require an opt-in flag or a separate
  uncertainty FASTA to obtain those hypotheses.
- Use stable, meaningful FASTA annotations with defined omission semantics.
  Keep ordinary headers simple and related fields recognizable by shared
  prefixes. Make filtering with simple shell commands straightforward,
  and provide an explicit opt-out for users who want only standard products.
- Keep delivery simple while preserving provenance and the distinction between
  an emitted hypothesis and a biologically validated sequence. More emitted
  records do not establish more validated loci, haplotypes, or correct repeats.
- Keep search and retention bounded, and report when those limits are reached.
  Inclusive delivery does not imply exhaustive enumeration.

Apply this preference when making future interface and output decisions.

## Development documentation

Start with [the developer index](dev_docs/README.md),
[the active plan](dev_docs/PLAN.md), and
[the architecture overview](dev_docs/OVERVIEW.md).
Use [CONTRIBUTING.md](CONTRIBUTING.md) for builds, checks, and releases.
Historical protocols and results preserve their dated evidence scope.

## Sequence annotations

Follow [semantic conventions for sequence annotations](dev_docs/SEQUENCE_ANNOTATIONS.md)
when designing or changing FASTA annotations and their associated evidence.
It is the canonical home for field meanings, shared prefixes, readable
values, and omission rules. It distinguishes agreed design from current
output and scopes breaking changes to released fields for v4.0.

## Commits

End agent-written commit messages with a `Co-Authored-By` trailer naming the
actual model and version. If the version is unavailable, say `version unknown`
rather than guessing.

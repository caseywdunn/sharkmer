# YAML memory ownership

`src/stats/yaml.rs` reserves stats, manifest and generated recovery-journal workspace under
`RetainedResults` when the development parent ledger is enabled. For the same immutable input,
this owner preserves YAML fields, ordering and formatting. This is one allocation owner under
[#139](https://github.com/caseywdunn/sharkmer/issues/139), not a whole-process
memory guarantee or a change to default backend selection.

## Audited implementation and input

The model is tied to exactly `serde_yaml_ng 0.10.0` and
`unsafe-libyaml 0.2.11`, pinned in Cargo.toml. The latter is already a transitive
dependency; the direct dependency supplies public type sizes. Sharkmer makes
no direct unsafe libyaml calls.

Production entry points accept immutable stats, manifest views and generated
recovery-journal views. The sealed manifest family also supports the historical
owned manifest model for compatibility checks. These trees contain derived structs,
maps, sequences, options, unit enum labels, primitive numbers and UTF-8 text.
The only custom container serializers serialize borrowed strings or slices.
Serde flattening emits entries into the enclosing map without collecting a
second tree. The preflight visits the same input through Serde, tracking maximum
scalar byte length and live nesting with fixed stack fields; it creates no
scalar copies, event vectors, YAML values or output strings.

Non-unit tagged variants, byte strings and `collect_str` are outside this model
and are rejected during preflight. Adding such forms requires extending this
audit, rather than using an uncharged fallback. Arbitrary custom serializers
with interior mutation are not supported by the two-pass model. Production
construction is limited to those explicit immutable models; the generic
manifest entry is sealed against unrelated serializers.

## Source-derived bound

The following paths in the pinned dependency sources establish the bound:

| Source | Allocation or lifetime fact |
| --- | --- |
| `serde_yaml_ng/src/ser.rs`, `Serializer::new`, `to_writer` | One boxed borrowed writer; one emitter owner. `to_writer` drops the emitter before returning, including on serialization failure. |
| `serde_yaml_ng/src/libyaml/emitter.rs`, `Emitter::new` | One pinned allocation containing `yaml_emitter_t`, a boxed writer trait object and optional I/O error. |
| `unsafe-libyaml/src/api.rs`, `yaml_emitter_initialize` | Buffers of 16,384 and 32,770 bytes; initial event, state, indentation and tag-directive arrays of 16 elements each. |
| `unsafe-libyaml/src/emitter.rs`, `yaml_emitter_need_more_events` | Mapping-start lookahead is three events, sequence-start two, document-start one; other head events drain immediately. There are at most four events after enqueue, including the current event. |
| `unsafe-libyaml/src/api.rs`, `yaml_scalar_event_initialize` | Every scalar event owns a copy of the UTF-8 bytes plus a NUL terminator. No scalar output-sized buffer is created for quoting or escaping. |
| `unsafe-libyaml/src/emitter.rs`, `yaml_emitter_append_tag_directive` | Untagged documents copy exactly four default directive strings: `!`, `!`, `!!` and `tag:yaml.org,2002:` including their NUL terminators. The directive array does not grow. |
| `unsafe-libyaml/src/macros.rs`, `STACK_INIT`, `PUSH`; `api.rs`, `yaml_stack_extend` | State and indentation arrays start at 16 and double. The bound includes both old and replacement allocations. |
| `unsafe-libyaml/src/emitter.rs`, node/collection/scalar handlers | A document return state and ancestor states bound the live state depth. Indentation includes ancestor collections and a transient scalar entry. Both fit within the scanned node depth plus one. |
| `unsafe-libyaml/src/lib.rs`, allocator shim | Each malloc/realloc allocation includes a stored length header, rounded to the greater alignment of u64 and usize. This is additional to allocator implementation overhead. |
| `serde_yaml_ng/src/ser.rs`, numeric serializers | Integer and float formatting uses stack-based itoa/ryu buffers. A 40-byte scalar floor covers i128::MIN, u128::MAX, floats, booleans and null. |
| `serde_yaml_ng/src/de.rs`, `visit_untagged_scalar`, `visit_int`, `parse_negative_int` | String style inference can format a negative radix-prefixed string. At most one temporary string is live at a time; repeated i64/i128 attempts do not overlap. |

Let S be the maximum UTF-8 scalar byte length, at least 40; D the maximum
serialized node depth; H the allocator shim's header size; and
C = next_power_of_two(max(16, D + 1)). Reserve the sum of:

- The pinned emitter allocation and boxed borrowed writer. The emitter bound
  sums the three field sizes plus padding after every field, without depending
  on private Rust field order.
- Both fixed output buffers, including H for each.
- Sixteen event entries and sixteen tag-directive entries, including H for each.
  Four-event lookahead cannot grow the sixteen-entry queue; tail exhaustion
  compacts existing entries.
- Twice each C-element state/indentation allocation, including headers, to cover
  realloc overlap. Earlier capacities are smaller than these reserved arrays.
- Four default tag strings, each including its terminator and H.
- Four (S + 1 + H) scalar copies, covering queued events and the event being
  initialized before enqueue.
- max(8, 2S) bytes for the one temporary negative-radix formatting string. Its
  writes are the minus sign and one borrowed substring; the allowance covers
  minimum String capacity and geometric growth. The temporary is dropped
  before the next style inference attempt.

Workspace sizing arithmetic is checked. Scalar byte lengths above i32::MAX fail
preflight before the dependency's usize-to-i32 cast. The reservation depends on
maximum scalar length and depth, not the number of targets or records. Original
stats data keeps its separate existing owner; workspace is charged in addition.

## Reservation lifetime and verification

`write_stats` constructs the prepared owner before `stage_file` creates the
staged stats file. Manifest and generated recovery-journal writes prepare the same owner before
creating their YAML temporary file. Views borrow receipt owners instead of
cloning metadata vectors and strings. A single numeric counter plan constructs the
shared ledger before output startup; the counter later consumes that same
ledger. Preparation creates no tables, workers or temporary directories.
Counter-input peak measurement starts with actual counter construction, so
completed startup workspace does not replace the count-phase peak; lifetime
peaks still include both phases. A configured refusal therefore happens before serializer
construction or stats-file creation. The prepared owner holds both the same
borrowed input and the reservation. Its write method calls the original
`serde_yaml_ng::to_writer`; emitter buffers, scalar events and the writer box
are destroyed while the reservation is still live. Success, writer error and
panic unwinding all release the reservation after emitter destruction.

Tests compare exact YAML bytes, exercise UTF-8/multiline/numeric-looking strings
and flattened fields, check typed refusal alongside another live owner, and
check cleanup on writer errors and panics. An isolated child process measures
allocation-free scanning and serializer heap peaks for large Unicode, multiline
and negative-radix strings, wide maps and depths across several doubling
boundaries. These measurements check the source-derived model; they do not
replace the source audit. Updating either pinned dependency requires repeating
both the audit and the allocation/compatibility checks.

The reservation covers logical serializer workspace, including libyaml's own
allocation headers. Allocator fragmentation/rounding, thread/call stacks,
error-reporting allocations, caller-owned writer buffers and upstream
configuration remain separate owners. Historical loader workspace now has
its own reservation described below. New
transaction labels, paths and staged receipts, prior-output cleanup paths and
duplicate tables, recovery receipt indexes and admitted restore plans have
separate parent reservations; they are not counted as serializer workspace.
Historical input and retained typed models now have their own reservations as
described below, in addition to loader workspace.
Tempfile-library path metadata, physical storage and complete process descriptor
ownership still require accounting.
The current startup ledger begins before panel/ENA collection, after raw
CLI/clap parsing. Standalone previews use their separate inherited-cap-aware
metadata ledger; raw CLI allocations remain unowned. The production writer is
a staged file, not an unbounded in-memory YAML buffer. See the
[remaining owner contracts](COUNT_STORE_CONTRACTS.md) before claiming a whole-run
RAM guarantee.

## Historical input and typed model ownership

`output::parsing` admits historical manifest and recovery-journal input before
copying bytes into a parent-owned buffer. It closes the file before decoding
and passes a borrowed slice to the pinned YAML deserializer, avoiding the
loader's separate `read_to_end` input copy. Input growth still reserves old and
replacement storage before allocation.

Concrete Serde seeds build only the existing manifest, journal and receipt
models. A seed reserves each String before copying scalar text and reserves
receipt-vector growth while its old allocation remains live. Scalar aliases
visit the same admitted path for every copied value. Field identifiers resolve
to static field names without another owned string. Optional and missing
fields, mapping and sequence forms, duplicate fields, ignored manifest fields
and rejected unknown journal fields retain their existing parsing behavior.
The original configured-limit error survives the Serde error boundary.

A parsed model keeps its original String/Vec data before an aggregate allocation
owner in field destruction order. The owner holds the actual admitted backing
capacities until all model data drops, including partial-parse failure. Previous
completed models keep that owner through the output transaction. An incomplete
prior model stays owned through its cleanup, then releases before ingestion.

## Historical loader workspace

`output::parsing::workspace` scans the admitted input without allocating and
reserves the pinned loader's workspace before creating its deserializer. The
same guard remains live until the deserialization call has destroyed the loader,
documents, scalar boxes and anchor dictionaries, including errors and unwinding.
Input and retained typed payload keep their distinct owners. Returned error
reporting, allocator rounding/fragmentation and call/thread stacks remain part
of the separate process-workspace requirements.

The audit covers `serde_yaml_ng 0.10.0` loader/parser/de code and
`unsafe-libyaml 0.2.11` reader/scanner/parser/api/macros. Rust container growth
and BTree layouts were checked in the primary
[1.88 Vec](https://raw.githubusercontent.com/rust-lang/rust/1.88.0/library/alloc/src/raw_vec/mod.rs),
[1.98.1 Vec](https://raw.githubusercontent.com/rust-lang/rust/1.98.1/library/alloc/src/raw_vec/mod.rs),
[1.88 BTree](https://raw.githubusercontent.com/rust-lang/rust/1.88.0/library/alloc/src/collections/btree/node.rs)
and [1.98.1 BTree](https://raw.githubusercontent.com/rust-lang/rust/1.98.1/library/alloc/src/collections/btree/node.rs)
sources. Updating dependencies or these standard-library assumptions requires
repeating the source audit and allocator checks. Native heap measurements here
use the recorded Linux/Rust toolchain; no other platform execution is implied.

### Boundary and token bound

Let N be the original byte count, A the raw `&` count and B the raw `!` count.
U is the sum of:

- starts of nonempty ASCII-whitespace-delimited runs;
- every raw occurrence of `[]{}:,?&*!|>%-`, either quote, CR or LF;
- every UTF-8 NEL, line-separator or paragraph-separator sequence.

The scan deliberately counts quoted, commented and block-scalar content too.
This avoids needing to interpret YAML before admission. The reader fixes UTF-8
encoding; invalid encodings terminate parsing instead of introducing hidden
syntactic boundaries.

In the pinned scanner, explicit indicators, quoted/block scalar starts, tags,
anchors, aliases and directives each have a counted indicator. Plain scalar
starts either begin a counted run or follow a counted flow/line boundary.
The implicit key/value and block-start/end tokens can be charged to the key,
colon or entry that introduces them. An implicit block contributes its matched
start and end to that introducing boundary. These charges fit four tokens per
counted unit; eight additional entries cover stream/document framing and an
empty document. Empty scalar events use key/value/entry indicators which have
no other scalar payload. Consequently token and stored event counts are bounded
by T = 4U + 8, even for malformed input before its error.

C = next_power_of_two(max(16, T)) bounds each growable array capacity. The
libyaml queues/stacks start with 16 elements and double. Rust event/skip vectors
use the audited doubling growth with minima below 16. Queue compaction reuses
storage. All multiplication, growth, input and component layout calculations
are checked before the loader is invoked.

### Reserved components

H is libyaml's allocation header, rounded to the greater alignment of u64 and
usize. Public libyaml type sizes and conservative field-padding sums determine
all element sizes; private struct field order is not assumed.

| Source and component | Reservation |
| --- | --- |
| `libyaml/parser.rs::ParserPinned` and `api.rs::yaml_parser_initialize` | One pinned parser/Cow owner, fixed 16,384-byte raw and 49,152-byte decoded buffers, plus headers. Reader refill compacts these buffers. |
| `scanner.rs` token queue, indentation and simple-key stacks; `parser.rs` state and mark stacks | Two C-element allocations per array, including headers, covering old and replacement storage. |
| `parser.rs::yaml_parser_process_directives` | Six C-element tag-directive arrays including headers: parser, temporary document and returned event arrays, with replacement overlap. One version-directive allocation is included. |
| `loader.rs::Document.events` | Four C-element event/mark arrays. Two documents may coexist while a second document is rejected; each permits replacement overlap. The largest Rust event contains scalar boxes, optional borrowed representation and style; its field/padding sum also bounds the smaller variants. |
| `de.rs::ignore_any` | Two C-element usize-sized arrays bound its smaller nesting-enum vector, including growth overlap. Ignored aliases do not recursively expand this vector. |
| `scanner.rs` and `api.rs` string helpers | Four simultaneous scanner strings use at most 8D capacity and one reallocation can overlap by D, where D=2N. The queued scalar strings need 2D, retained scalar boxes D, and anchor/directive copy lifetimes add 8N. Per-token minimum string capacity adds 64T and fixed minima add 512 bytes. |
| Explicit tag prefix expansion in `parser.rs::yaml_parser_parse_node` and Rust tag boxes | 4(N+32)B bytes. Every explicit tag has a raw `!`; its expanded prefix/suffix fits N plus the default-prefix allowance. Retained, C-event, concatenation and directive-copy lifetimes are covered separately. |
| libyaml payload allocation headers | (8T+32)H additional bytes for token payloads, directive copies, scratch and overlap; control-array headers are already included above. |
| `loader.rs` anchor and alias BTreeMaps | Two node allowances per raw source anchor for each dictionary. Each node bound includes 11 keys/values, 12 child pointers, a parent pointer, two u16 fields and padding; it bounds internal and leaf nodes. Splitting does not require copying the existing tree. The source-wide anchor count covers overlapping documents. |

D=2N includes YAML `\L` and `\P`: two original bytes can produce three UTF-8
bytes. Other escapes, raw Unicode and line folding also fit that bound. Scalar
aliases are stored as event indices; their repeated typed text is admitted by
the separate value seeds rather than multiplied into this document bound.
Expanded tag text is different: prefix reuse can repeat a long directive before
value seeds run, so it receives the explicit product term above.

The reservation is conservative and can substantially exceed measured heap for
long quoted/comment/block-scalar text. For example, the retained large block
fixture reserves about 4.9 GB numerically while its loader heap peak is about
6.6 MB. A reservation does not allocate those bytes; it admits the proven upper
bound against the configured ledger. A tighter scan needs its own proof and
compatibility checks. No existing finite-memory or native fault protocol ceiling
is increased to accommodate this change.

### Verification and remaining scope

An isolated process checks allocation-free scanning, loader heap peaks and zero
residual heap over 28 shapes, including deep/wide structures, anchors, Unicode,
block scalars, escaped line separators, many directives, prefix reuse, malformed
aliases and multiple documents. The same process checks refusal without parser
buffer allocation and reservation lifetime through a visitor panic. These
measurements check the source-derived bound; they do not replace the derivation.
Input/typed-model, ordinary output, finite-descriptor and durable recovery gates
remain required. Original process and filesystem ownership gaps still require
separate work before a whole-run RAM guarantee or default-backend promotion.

## Panel configuration models

The development panel loader uses the same audited loader workspace before
calling the existing Serde derives. `src/configuration/tracking.rs` wraps value
visitors, retaining field/default/unknown-field rules while admitting owned
strings, sequences and maps. The supported panel models contain inline
structs/options, scalar or empty defaults, `String`, `Vec` and `BTreeMap`.
Adding boxes, shared pointers, other containers or allocating defaults requires
extending this audit and its allocation checks.

`serde_core 1.0.228/src/de/impls.rs` establishes that its `StringVisitor` copies
borrowed text and accepts owned strings without another copy. Both routes
receive admission before model ownership; parser-owned text remains covered by
the loader guard during transfer. Its `VecVisitor` starts with the supplied
size hint. The wrapper supplies no hint, then admits growing element storage
before each seed. Cumulative previous-capacity charges cover replacement overlap
and are reduced only after an explicit inventory of surviving allocations.

The current compiler is Rust 1.98.1, commit
`48a229ceaefd4985c50990b14116b6d856af0985`. Its
[vector allocation implementation](https://github.com/rust-lang/rust/blob/48a229ceaefd4985c50990b14116b6d856af0985/library/alloc/src/raw_vec/mod.rs)
doubles capacity and uses a minimum initial capacity no greater than eight.
The sequence guard covers those cases using the actual decoded element size.
The compiler's
[B-tree node implementation](https://github.com/rust-lang/rust/blob/48a229ceaefd4985c50990b14116b6d856af0985/library/alloc/src/collections/btree/node.rs)
has eleven key/value cells and twelve child pointers. Two maximum node
allowances per admitted map entry conservatively cover accumulated leaf/branch
allocations and splits for these models. Compiler and dependency sources must
remain bound to the execution evidence; measurements alone do not establish
these growth rules for a different implementation.

Panel input is read into a fallible, charged byte buffer with the previous strict
UTF-8 requirement. URL loading retains its timeouts and decoded-body size limit
and uses a separate automatic HTTP-gzip reservation. Requested built-in panels
are decoded individually. Derived target names and source labels are admitted
before formatting; uniqueness tables borrow their keys. After unused root
metadata and temporary arrays drop, an explicit capacity inventory reduces the
owner to the surviving primer strings and arrays. Collected records retain
those owners through downstream use. Inline configuration and validation
reports use the same run ledger, with errors retaining their text guards until
disposal.

Process-envelope preflight, standalone-preview ceilings, transport/TLS and
complete descriptor/physical filesystem ownership remain distinct requirements.
The active plan records their verification status; this model does not waive
those gates or permit larger benchmark ceilings.

### Retained decoder errors

A decoder failure can retain formatted scalar/path text after its document and
partial typed model have dropped. The shared loader workspace now transfers
with such errors in configuration and historical manifest/journal readers.
The wrapper drops the library error before its reservation and preserves the
original diagnostic text and underlying source chain. It conservatively keeps
the original loader allowance until error disposal rather than guessing a
smaller library-error layout.

A captured typed admission failure follows a separate path: the temporary Serde
error drops while the workspace still lives, then the original resource error
returns with its type and downcast behavior intact. No decoded model survives
that path. The active plan records verification of the lifetime follow-up;
other validation/library/transport owners remain separate work.

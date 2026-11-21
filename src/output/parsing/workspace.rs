//! Precharged workspace for serde_yaml_ng 0.10.0 / unsafe-libyaml 0.2.11.
//! See the source derivation in dev_docs/YAML_WORKSPACE.md. The scan never
//! interprets YAML: it overcounts lexical boundaries, including quoted text.
use super::*;
use std::borrow::Cow;
use std::mem::{align_of, size_of};

/// The deserializer and its loader drop inside the call while this guard lives.
pub(crate) struct Workspace {
    _reservation: Option<MemoryReservation>,
}
impl Workspace {
    pub(crate) fn new(input: &[u8], budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let reservation = budget
            .map(|b| {
                b.reserve(MemoryClass::RetainedResults, bytes(input)?)
                    .context("Unable to admit historical YAML loader workspace")
            })
            .transpose()?;
        Ok(Self {
            _reservation: reservation,
        })
    }

    /// Keep decoder-owned error text covered until the error is disposed.
    /// A typed admission failure is returned separately without this wrapper.
    pub(crate) fn retain_error(self, error: serde_yaml_ng::Error) -> anyhow::Error {
        anyhow::Error::new(Failed {
            inner: error,
            _workspace: self,
        })
    }
}

// Decoder failures can own formatted scalar/path text after the loader returns.
// Drop that payload before releasing its source-derived workspace allowance.
struct Failed {
    inner: serde_yaml_ng::Error,
    _workspace: Workspace,
}
impl std::fmt::Debug for Failed {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(&self.inner, formatter)
    }
}
impl std::fmt::Display for Failed {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.inner, formatter)
    }
}
impl std::error::Error for Failed {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        std::error::Error::source(&self.inner)
    }
}

/// Count potential boundaries conservatively without collecting tokens/scalars.
#[derive(Debug, Default, PartialEq, Eq)]
struct Shape {
    input: usize,
    units: usize,
    anchors: usize,
    bangs: usize,
}
impl Shape {
    fn scan(input: &[u8]) -> Result<Self> {
        let mut shape = Self {
            input: input.len(),
            ..Self::default()
        };
        let mut word = false;
        for (offset, byte) in input.iter().enumerate() {
            let nonwhite = !byte.is_ascii_whitespace();
            let run = usize::from(nonwhite && !word);
            let syntax = usize::from(b"[]{}:,?&*!|>%-\"'\r\n".contains(byte));
            let rest = &input[offset..];
            let unicode_break = usize::from(
                rest.starts_with(b"\xc2\x85")
                    || rest.starts_with(b"\xe2\x80\xa8")
                    || rest.starts_with(b"\xe2\x80\xa9"),
            );
            shape.units = shape
                .units
                .checked_add(run + syntax + unicode_break)
                .context("Historical YAML boundary count overflow")?;
            shape.anchors = shape
                .anchors
                .checked_add(usize::from(*byte == b'&'))
                .context("Historical YAML anchor count overflow")?;
            shape.bangs = shape
                .bangs
                .checked_add(usize::from(*byte == b'!'))
                .context("Historical YAML tag count overflow")?;
            word = nonwhite;
        }
        Ok(shape)
    }

    fn bytes(&self) -> Result<u64> {
        ensure!(
            self.input <= isize::MAX as usize,
            "Historical YAML input layout overflow"
        );
        ensure!(
            self.anchors <= self.input && self.bangs <= self.input,
            "Invalid historical YAML shape"
        );
        let tokens = self
            .units
            .checked_mul(4)
            .and_then(|n| n.checked_add(8))
            .context("Historical YAML token bound overflow")?;
        let capacity = tokens
            .max(16)
            .checked_next_power_of_two()
            .context("Historical YAML workspace capacity overflow")?;
        let n = u64::try_from(self.input)?;
        let t = u64::try_from(tokens)?;
        let c = u64::try_from(capacity)?;
        let alignment = align_of::<u64>().max(align_of::<usize>());
        let header = u64::try_from(size_of::<usize>().next_multiple_of(alignment))?;
        // A libyaml array may retain its old and new allocations during growth.
        let array = |element: usize, copies: u64| -> Result<u64> {
            let one = mul(c, element as u64)?;
            ensure!(
                one <= isize::MAX as u64,
                "Historical YAML array layout overflow"
            );
            mul(add(one, header)?, copies)
        };
        let parser_alignment =
            align_of::<unsafe_libyaml::yaml_parser_t>().max(align_of::<Cow<'_, [u8]>>());
        let parser = size_of::<unsafe_libyaml::yaml_parser_t>()
            + size_of::<Cow<'_, [u8]>>()
            + 2 * parser_alignment;
        // Rust Event's largest variant has three slice boxes, an optional
        // borrowed representation and a style. Include enum/field/tuple padding
        // and the accompanying libyaml Mark without private layout assumptions.
        let event_cell = 3 * size_of::<Option<Box<[u8]>>>()
            + size_of::<Option<&[u8]>>()
            + size_of::<usize>()
            + size_of::<unsafe_libyaml::yaml_mark_t>()
            + 8 * align_of::<usize>().max(align_of::<unsafe_libyaml::yaml_mark_t>());
        let mut total = parser as u64;
        for part in [
            add(16_384, header)?,
            add(49_152, header)?,
            array(size_of::<unsafe_libyaml::yaml_token_t>(), 2)?,
            array(size_of::<i32>(), 2)?,
            array(size_of::<unsafe_libyaml::yaml_simple_key_t>(), 2)?,
            array(size_of::<unsafe_libyaml::yaml_parser_state_t>(), 2)?,
            array(size_of::<unsafe_libyaml::yaml_mark_t>(), 2)?,
            // Parser directives, temporary document directives and their
            // returned event array; each permits replacement overlap.
            array(size_of::<unsafe_libyaml::yaml_tag_directive_t>(), 6)?,
            // At most two loader Documents overlap, including rejected second
            // documents; each array permits old/replacement overlap.
            mul(mul(c, event_cell as u64)?, 4)?,
            // IgnoredAny has one growing Nest vector; usize exceeds its enum.
            mul(mul(c, size_of::<usize>() as u64)?, 2)?,
            add(
                size_of::<unsafe_libyaml::yaml_version_directive_t>() as u64,
                header,
            )?,
        ] {
            total = add(total, part)?;
        }
        // D=2N covers decoded UTF-8, including two-byte \L/\P escapes producing
        // three-byte line separators. Four scanner strings can each approach
        // 2D capacity; one realloc overlaps by at most D (plus initial minima).
        // Queued scalar capacity is at most 2D plus per-token minima; retained
        // Rust scalar boxes need D. Anchor/directive copies contribute 8N.
        let decoded = mul(n, 2)?;
        total = add(total, mul(decoded, 12)?)?; // 9D scratch + 2D queue + D boxes
        total = add(total, mul(n, 8)?)?;
        total = add(total, add(mul(t, 64)?, 512)?)?;
        // Prefix expansion is not linear in input bytes. Each explicit tag has
        // a raw '!' and its fully expanded text is <= N+32 (including defaults).
        // Cover retained, C-event, concatenation and directive-copy lifetimes.
        total = add(total, mul(mul(add(n, 32)?, self.bangs as u64)?, 4)?)?;
        // Allocation headers for token payloads, directive copies and scratch,
        // including overlap. Control-array headers were included above.
        total = add(total, mul(add(mul(t, 8)?, 32)?, header)?)?;
        // Audited Rust BTree nodes: 11 keys/values, 12 child pointers, parent
        // pointer, two u16s and padding. Anchor keys are slice boxes; aliases
        // use usize keys. Two allowances per source anchor cover split overlap
        // and both documents; raw '&' also overcounts quoted/unused anchors.
        let node = |key: usize| -> Result<u64> {
            let entries = mul(11, (key + size_of::<usize>()) as u64)?;
            add(
                entries,
                (13 * size_of::<usize>() + 2 * size_of::<u16>() + 7 * align_of::<usize>()) as u64,
            )
        };
        let tree_pair = add(node(size_of::<Box<[u8]>>())?, node(size_of::<usize>())?)?;
        total = add(total, mul(mul(self.anchors as u64, 2)?, tree_pair)?)?;
        Ok(total)
    }
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b)
        .context("Historical YAML workspace size overflow")
}
fn mul(a: u64, b: u64) -> Result<u64> {
    a.checked_mul(b)
        .context("Historical YAML workspace size overflow")
}

pub(super) fn bytes(input: &[u8]) -> Result<u64> {
    Shape::scan(input)?.bytes()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};
    use serde::de::IgnoredAny;

    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: [bytes; MEMORY_CLASSES],
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }
    fn fixtures() -> Vec<(String, Vec<u8>)> {
        let mut cases = vec![
            ("empty".to_owned(), Vec::new()),
            ("plain".to_owned(), b"key: value\n".to_vec()),
        ];
        for n in [1, 16, 128, 1024] {
            cases.push((
                format!("flow-{n}"),
                format!("[{}]", "[null, null],".repeat(n)).into_bytes(),
            ));
            cases.push((
                format!("anchors-{n}"),
                (0..n)
                    .map(|i| format!("- &a{i} [one, two]\n"))
                    .collect::<String>()
                    .into_bytes(),
            ));
            cases.push((
                format!("deep-{n}"),
                format!("{}null{}", "[".repeat(n), "]".repeat(n)).into_bytes(),
            ));
        }
        for n in [31, 16384, 262144] {
            cases.push((
                format!("unicode-{n}"),
                format!("key: \"{}\"\n", "λ🧬".repeat(n)).into_bytes(),
            ));
            cases.push((
                format!("block-{n}"),
                format!("key: |\n{}", "  row text\n".repeat(n)).into_bytes(),
            ));
        }
        for (length, count) in [(32, 8), (32768, 128)] {
            cases.push((
                format!("tag-{length}-{count}"),
                format!(
                    "%TAG !e! tag:example.org,2026:{}\n---\n{}",
                    "a".repeat(length),
                    "- !e!field value\n".repeat(count)
                )
                .into_bytes(),
            ));
        }
        cases.push((
            "two-documents".to_owned(),
            b"---\nkey: value\n---\nsecond: [one, two]\n".to_vec(),
        ));
        cases.push(("missing-alias".to_owned(), b"key: *missing\n".to_vec()));
        cases.push((
            "escaped-breaks".to_owned(),
            format!("key: \"{}\"\n", "\\L\\P".repeat(4096)).into_bytes(),
        ));
        cases.push((
            "unicode-breaks".to_owned(),
            "key:\u{85}  one\u{2028}next:\u{2029}  two\n"
                .as_bytes()
                .to_vec(),
        ));
        cases.push((
            "quoted-bangs".to_owned(),
            format!("key: \"{}\"\n", "!".repeat(4096)).into_bytes(),
        ));
        cases.push((
            "many-directives".to_owned(),
            format!(
                "{}---\nkey: value\n",
                (0..128)
                    .map(|i| format!("%TAG !e{i}! tag:example.org,{i}:\n"))
                    .collect::<String>()
            )
            .into_bytes(),
        ));
        cases
    }

    #[test]
    fn loader_shape_counts_boundaries_and_checks_allocation_arithmetic() {
        assert_eq!(Shape::scan(b"word").unwrap().units, 1);
        assert_eq!(Shape::scan(b"one two").unwrap().units, 2);
        assert_eq!(Shape::scan(b"[]").unwrap().units, 3);
        assert_eq!(
            Shape::scan("a\u{85}b\u{2028}c\u{2029}d".as_bytes())
                .unwrap()
                .units,
            4
        );
        assert_eq!(Shape::scan(b"\"&not_anchor !not_tag\"").unwrap().anchors, 1);
        assert_eq!(Shape::scan(b"\"&not_anchor !not_tag\"").unwrap().bangs, 1);
        for shape in [
            Shape {
                units: usize::MAX,
                ..Shape::default()
            },
            Shape {
                input: usize::MAX,
                bangs: usize::MAX,
                ..Shape::default()
            },
            Shape {
                anchors: usize::MAX,
                ..Shape::default()
            },
        ] {
            assert!(shape.bytes().is_err());
        }
        let plain = bytes(b"value: plain").unwrap();
        let tagged = bytes(b"value: !custom plain").unwrap();
        assert!(tagged > plain);
        let b = budget(plain);
        let other = b.reserve(MemoryClass::RetainedResults, 1).unwrap();
        let error = Workspace::new(b"value: plain", Some(&b)).err().unwrap();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(b.snapshot().reserved_memory_bytes, 1);
        drop(other);
        let workspace = Workspace::new(b"value: plain", Some(&b)).unwrap();
        assert_eq!(b.snapshot().reserved_memory_bytes, plain);
        drop(workspace);
        assert_eq!(b.snapshot().reserved_memory_bytes, 0);
        assert!(b.snapshot().valid);
    }

    #[test]
    fn loader_workspace_covers_isolated_allocator_peaks_and_cleanup_paths() {
        const ENV: &str = "SHARKMER_LOADER_WORKSPACE_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status=std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact","output::parsing::workspace::tests::loader_workspace_covers_isolated_allocator_peaks_and_cleanup_paths","--nocapture","--test-threads=1"])
                .env(ENV,"1").status().unwrap();
            assert!(status.success());
            return;
        }
        let cases = fixtures();
        assert_eq!(cases.len(), 28);
        for (name, input) in cases {
            let bound = bytes(&input).unwrap();
            let budget = budget(bound + 4096);
            let _ = budget.snapshot();
            let baseline = crate::PEAK_ALLOC.current_usage();
            crate::PEAK_ALLOC.reset_peak_usage();
            let shape = Shape::scan(&input).unwrap();
            assert_eq!(shape.bytes().unwrap(), bound);
            assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
            assert_eq!(crate::PEAK_ALLOC.peak_usage(), baseline);
            let workspace = Workspace::new(&input, Some(&budget)).unwrap();
            let result = serde_yaml_ng::from_slice::<IgnoredAny>(&input);
            let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
            assert!(
                peak as u64 <= bound,
                "loader {name} heap {peak} exceeded {bound}"
            );
            let ok = result.is_ok();
            drop(result);
            assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
            assert_eq!(budget.snapshot().reserved_memory_bytes, bound);
            drop(workspace);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
            println!(
                "loader_allocator case={name} input_bytes={} units={} peak_bytes={peak} reserved_bytes={bound} ok={ok}",
                input.len(),
                shape.units
            );
        }
        let input = b"value: plain";
        let bound = bytes(input).unwrap();
        let small = budget(bound - 1);
        let _ = small.snapshot();
        let baseline = crate::PEAK_ALLOC.current_usage();
        crate::PEAK_ALLOC.reset_peak_usage();
        let error = Workspace::new(input, Some(&small)).err().unwrap();
        let refusal_peak = crate::PEAK_ALLOC.peak_usage() - baseline;
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert!(
            refusal_peak < 65536,
            "refusal allocated parser buffers: {refusal_peak}"
        );
        drop(error);
        assert_eq!(crate::PEAK_ALLOC.current_usage(), baseline);
        assert_eq!(small.snapshot().reserved_memory_bytes, 0);
        struct Panicking;
        impl<'de> Deserialize<'de> for Panicking {
            fn deserialize<D: de::Deserializer<'de>>(de: D) -> std::result::Result<Self, D::Error> {
                struct Panic;
                impl<'de> Visitor<'de> for Panic {
                    type Value = Panicking;
                    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                        f.write_str("injected panic")
                    }
                    fn visit_map<A: MapAccess<'de>>(
                        self,
                        _: A,
                    ) -> std::result::Result<Self::Value, A::Error> {
                        panic!("injected loader visitor panic")
                    }
                }
                de.deserialize_map(Panic)
            }
        }
        let budget = budget(bound);
        let workspace = Workspace::new(input, Some(&budget)).unwrap();
        assert!(
            std::panic::catch_unwind(|| serde_yaml_ng::from_slice::<Panicking>(input)).is_err()
        );
        assert_eq!(budget.snapshot().reserved_memory_bytes, bound);
        drop(workspace);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }
}

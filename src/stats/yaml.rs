//! Precharged workspace for the pinned YAML emitter used by immutable run stats.
//!
//! The scan borrows scalar text and tracks nesting without collecting values.
//! The source audit and allocation model live in dev_docs/YAML_WORKSPACE.md.
use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result};
use serde::Serialize;
use serde::ser::{self, Impossible};
use std::fmt;
use std::io::{self, Write};
use std::mem::{align_of, size_of};
use std::sync::Arc;

/// Reserve before the output transaction creates the staged stats file.
pub(super) fn prepare(stats: &super::RunStats) -> Result<Prepared<'_, super::RunStats>> {
    Prepared::new(stats, stats.sample.budget())
}

/// Reserve the same audited forms for an immutable publication manifest.
pub(crate) fn prepare_manifest<'a, M: crate::output::ManifestModel>(
    manifest: &'a M,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Prepared<'a, M>> {
    Prepared::new(manifest, budget)
}

/// Precharge the immutable generated recovery-journal model before staging it.
pub(crate) fn prepare_recovery<'a, 'b>(
    journal: &'a crate::output::RecoveryJournalView<'b>,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Prepared<'a, crate::output::RecoveryJournalView<'b>>> {
    Prepared::new(journal, budget)
}

/// Reserve the audited immutable emitter for cached-source metadata.
pub(crate) fn prepare_cache<'a>(
    meta: &'a crate::cache::CacheMeta,
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Prepared<'a, crate::cache::CacheMeta>> {
    Prepared::new(meta, budget)
}

/// Keep the exact immutable input and its guard together through emitter drop.
pub(crate) struct Prepared<'a, T: ?Sized> {
    value: &'a T,
    _reservation: Option<MemoryReservation>,
}

impl<'a, T: ?Sized + Serialize> Prepared<'a, T> {
    // Production entry points accept only immutable derived stats/manifests.
    // Custom serializers with interior mutation would invalidate a two-pass scan.
    pub(super) fn new(value: &'a T, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let reservation = budget
            .map(|budget| {
                let bytes = Shape::scan(value)?.workspace_bytes()?;
                budget.reserve(MemoryClass::RetainedResults, bytes)
            })
            .transpose()?;
        Ok(Self {
            value,
            _reservation: reservation,
        })
    }

    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        // to_writer destroys its serializer (including queued events) before
        // returning on either success or error. The reservation still lives.
        serde_yaml_ng::to_writer(writer, self.value).context("Failed to write YAML")
    }
}

#[derive(Debug, Default, PartialEq, Eq)]
struct Shape {
    max_scalar_bytes: usize,
    max_depth: usize,
    depth: usize,
}

impl Shape {
    fn scan(value: &(impl ?Sized + Serialize)) -> Result<Self, ScanError> {
        // 40 bytes covers signed i128; ryu's longest finite float fits in 24.
        let mut shape = Self {
            max_scalar_bytes: 40,
            ..Self::default()
        };
        value.serialize(&mut shape)?;
        Ok(shape)
    }

    fn scalar(&mut self, bytes: usize) -> Result<(), ScanError> {
        // libyaml's wrapper casts scalar lengths to i32. Reject before it can
        // reinterpret a too-long string as a NUL-terminated scalar.
        if bytes > i32::MAX as usize {
            return Err(ScanError(
                "Stats YAML scalar exceeds libyaml's length limit",
            ));
        }
        self.max_scalar_bytes = self.max_scalar_bytes.max(bytes);
        self.max_depth = self.max_depth.max(self.depth + 1);
        Ok(())
    }

    fn enter(&mut self) -> Result<Compound<'_>, ScanError> {
        self.depth = self
            .depth
            .checked_add(1)
            .ok_or(ScanError("Stats YAML nesting overflow"))?;
        self.max_depth = self.max_depth.max(self.depth);
        Ok(Compound(self))
    }

    fn workspace_bytes(&self) -> Result<u64> {
        // unsafe-libyaml's malloc adds a length header rounded to this alignment.
        let alignment = align_of::<u64>().max(align_of::<usize>());
        let header = size_of::<usize>().next_multiple_of(alignment);
        let allocation = |bytes: usize| -> Result<usize> {
            bytes
                .checked_add(header)
                .context("YAML allocation size overflow")
        };
        // One document return state plus at most one per nested node; indentation
        // similarly includes a transient scalar entry. Stack capacities double.
        let entries = self
            .max_depth
            .checked_add(1)
            .and_then(|n| n.max(16).checked_next_power_of_two())
            .context("YAML stack capacity overflow")?;
        let stack = |element: usize| -> Result<usize> {
            let bytes = entries
                .checked_mul(element)
                .context("YAML stack size overflow")?;
            // Both old and replacement allocations can coexist in realloc.
            allocation(bytes)?
                .checked_mul(2)
                .context("YAML stack overlap overflow")
        };
        let scalar = allocation(
            self.max_scalar_bytes
                .checked_add(1)
                .context("YAML scalar size overflow")?,
        )?
        .checked_mul(4)
        .context("YAML queued scalar size overflow")?;
        // Style inference creates one temporary negative radix string at a time.
        // This covers String's minimum capacity and geometric growth as well.
        let style = self
            .max_scalar_bytes
            .checked_mul(2)
            .context("YAML scalar formatting overflow")?
            .max(8);
        // EmitterPinned has these three fields. Allow padding after each field
        // rather than relying on the private Rust struct's field order/layout.
        let emitter_alignment = align_of::<unsafe_libyaml::yaml_emitter_t>()
            .max(align_of::<Box<dyn Write>>())
            .max(align_of::<Option<io::Error>>());
        let emitter = size_of::<unsafe_libyaml::yaml_emitter_t>()
            + size_of::<Box<dyn Write>>()
            + size_of::<Option<io::Error>>()
            + 3 * emitter_alignment;
        let pieces = [
            emitter,
            size_of::<&mut ()>(), // Box<&mut W> in Serializer::new
            allocation(16_384)?,
            allocation(32_770)?,
            allocation(16 * size_of::<unsafe_libyaml::yaml_event_t>())?,
            allocation(16 * size_of::<unsafe_libyaml::yaml_tag_directive_t>())?,
            stack(size_of::<unsafe_libyaml::yaml_emitter_state_t>())?,
            stack(size_of::<i32>())?,
            // The only four strdup calls for an untagged document.
            allocation(b"!\0".len())?,
            allocation(b"!\0".len())?,
            allocation(b"!!\0".len())?,
            allocation(b"tag:yaml.org,2002:\0".len())?,
            scalar,
            style,
        ];
        let bytes = pieces.into_iter().try_fold(0usize, |total, n| {
            total.checked_add(n).context("YAML workspace size overflow")
        })?;
        u64::try_from(bytes).context("YAML workspace exceeds u64")
    }
}

#[derive(Debug)]
struct ScanError(&'static str);
impl fmt::Display for ScanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for ScanError {}
impl ser::Error for ScanError {
    fn custom<T: fmt::Display>(_: T) -> Self {
        Self("Custom stats serializer is outside the audited YAML model")
    }
}

struct Compound<'a>(&'a mut Shape);
impl Drop for Compound<'_> {
    fn drop(&mut self) {
        self.0.depth -= 1;
    }
}

macro_rules! primitive {
    ($($method:ident($ty:ty)),* $(,)?) => {$ (
        fn $method(self, _: $ty) -> Result<(), ScanError> { self.scalar(40) }
    )*};
}

impl<'a> ser::Serializer for &'a mut Shape {
    type Ok = ();
    type Error = ScanError;
    type SerializeSeq = Compound<'a>;
    type SerializeTuple = Compound<'a>;
    type SerializeTupleStruct = Compound<'a>;
    type SerializeMap = Compound<'a>;
    type SerializeStruct = Compound<'a>;
    type SerializeTupleVariant = Impossible<(), ScanError>;
    type SerializeStructVariant = Impossible<(), ScanError>;

    primitive! {
        serialize_bool(bool), serialize_i8(i8), serialize_i16(i16), serialize_i32(i32),
        serialize_i64(i64), serialize_i128(i128), serialize_u8(u8), serialize_u16(u16),
        serialize_u32(u32), serialize_u64(u64), serialize_u128(u128),
        serialize_f32(f32), serialize_f64(f64), serialize_char(char),
    }
    fn serialize_str(self, value: &str) -> Result<(), ScanError> {
        self.scalar(value.len())
    }
    fn serialize_unit(self) -> Result<(), ScanError> {
        self.scalar(4)
    }
    fn serialize_none(self) -> Result<(), ScanError> {
        self.serialize_unit()
    }
    fn serialize_unit_struct(self, _: &'static str) -> Result<(), ScanError> {
        self.serialize_unit()
    }
    fn serialize_unit_variant(
        self,
        _: &'static str,
        _: u32,
        variant: &'static str,
    ) -> Result<(), ScanError> {
        self.serialize_str(variant)
    }
    fn serialize_some<T: ?Sized + Serialize>(self, value: &T) -> Result<(), ScanError> {
        value.serialize(self)
    }
    fn serialize_newtype_struct<T: ?Sized + Serialize>(
        self,
        _: &'static str,
        value: &T,
    ) -> Result<(), ScanError> {
        value.serialize(self)
    }
    fn serialize_seq(self, _: Option<usize>) -> Result<Compound<'a>, ScanError> {
        self.enter()
    }
    fn serialize_tuple(self, _: usize) -> Result<Compound<'a>, ScanError> {
        self.enter()
    }
    fn serialize_tuple_struct(self, _: &'static str, _: usize) -> Result<Compound<'a>, ScanError> {
        self.enter()
    }
    fn serialize_map(self, _: Option<usize>) -> Result<Compound<'a>, ScanError> {
        self.enter()
    }
    fn serialize_struct(self, _: &'static str, _: usize) -> Result<Compound<'a>, ScanError> {
        self.enter()
    }
    // These paths introduce tag/formatting allocations not used by RunStats.
    // A future schema change must extend the audit instead of escaping the cap.
    fn serialize_bytes(self, _: &[u8]) -> Result<(), ScanError> {
        Err(ScanError(
            "Byte strings are outside the audited stats YAML model",
        ))
    }
    fn collect_str<T: ?Sized + fmt::Display>(self, _: &T) -> Result<(), ScanError> {
        Err(ScanError(
            "Formatted strings are outside the audited stats YAML model",
        ))
    }
    fn serialize_newtype_variant<T: ?Sized + Serialize>(
        self,
        _: &'static str,
        _: u32,
        _: &'static str,
        _: &T,
    ) -> Result<(), ScanError> {
        Err(ScanError(
            "Tagged variants are outside the audited stats YAML model",
        ))
    }
    fn serialize_tuple_variant(
        self,
        _: &'static str,
        _: u32,
        _: &'static str,
        _: usize,
    ) -> Result<Self::SerializeTupleVariant, ScanError> {
        Err(ScanError(
            "Tagged variants are outside the audited stats YAML model",
        ))
    }
    fn serialize_struct_variant(
        self,
        _: &'static str,
        _: u32,
        _: &'static str,
        _: usize,
    ) -> Result<Self::SerializeStructVariant, ScanError> {
        Err(ScanError(
            "Tagged variants are outside the audited stats YAML model",
        ))
    }
}

macro_rules! sequence {
    ($trait:ident, $method:ident) => {
        impl ser::$trait for Compound<'_> {
            type Ok = ();
            type Error = ScanError;
            fn $method<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), ScanError> {
                value.serialize(&mut *self.0)
            }
            fn end(self) -> Result<(), ScanError> {
                Ok(())
            }
        }
    };
}
sequence!(SerializeSeq, serialize_element);
sequence!(SerializeTuple, serialize_element);
sequence!(SerializeTupleStruct, serialize_field);
impl ser::SerializeMap for Compound<'_> {
    type Ok = ();
    type Error = ScanError;
    fn serialize_key<T: ?Sized + Serialize>(&mut self, key: &T) -> Result<(), ScanError> {
        key.serialize(&mut *self.0)
    }
    fn serialize_value<T: ?Sized + Serialize>(&mut self, value: &T) -> Result<(), ScanError> {
        value.serialize(&mut *self.0)
    }
    fn end(self) -> Result<(), ScanError> {
        Ok(())
    }
}
impl ser::SerializeStruct for Compound<'_> {
    type Ok = ();
    type Error = ScanError;
    fn serialize_field<T: ?Sized + Serialize>(
        &mut self,
        key: &'static str,
        value: &T,
    ) -> Result<(), ScanError> {
        self.0.scalar(key.len())?;
        value.serialize(&mut *self.0)
    }
    fn end(self) -> Result<(), ScanError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};
    use serde::ser::{SerializeMap, SerializeSeq};
    use std::sync::atomic::{AtomicBool, Ordering};

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

    struct Nested<'a> {
        depth: usize,
        text: &'a str,
        width: usize,
    }
    impl Serialize for Nested<'_> {
        fn serialize<S: ser::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
            if self.depth > 0 {
                let mut seq = serializer.serialize_seq(Some(1))?;
                seq.serialize_element(&Nested {
                    depth: self.depth - 1,
                    ..*self
                })?;
                seq.end()
            } else {
                let mut map = serializer.serialize_map(Some(self.width))?;
                for _ in 0..self.width {
                    map.serialize_entry("text", self.text)?;
                }
                map.end()
            }
        }
    }

    #[test]
    fn yaml_scan_tracks_utf8_nesting_and_ignores_collection_width() {
        let text = "β🧬".repeat(100);
        let narrow = Nested {
            depth: 17,
            text: &text,
            width: 1,
        };
        let wide = Nested {
            width: 10_000,
            ..narrow
        };
        let a = Shape::scan(&narrow).unwrap();
        let b = Shape::scan(&wide).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.max_scalar_bytes, 600);
        assert_eq!(a.max_depth, 19);
        assert_eq!(a.depth, 0);
        assert!(
            a.workspace_bytes().unwrap()
                > Shape::scan(&"short").unwrap().workspace_bytes().unwrap()
        );
    }

    #[derive(Serialize)]
    struct Flat<'a> {
        scalar: &'a str,
        max_integer: u128,
        min_integer: i128,
        float: f64,
        empty: Option<u32>,
    }
    #[derive(Serialize)]
    struct Fixture<'a> {
        rows: Vec<&'a str>,
        #[serde(flatten)]
        flat: Flat<'a>,
    }
    fn fixture() -> Fixture<'static> {
        Fixture {
            rows: vec![
                "true",
                "-0xff",
                "-0o123",
                "-0b101",
                ".nan",
                "λ🧬",
                "line\nline",
                "!tag",
                "001",
                "",
            ],
            flat: Flat {
                scalar: "null",
                max_integer: u128::MAX,
                min_integer: i128::MIN,
                float: f64::MIN,
                empty: None,
            },
        }
    }

    struct CheckedWriter<'a> {
        budget: &'a ResourceBudget,
        required: u64,
        bytes: Vec<u8>,
        fail: bool,
    }
    impl Write for CheckedWriter<'_> {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            assert_eq!(self.budget.snapshot().reserved_memory_bytes, self.required);
            if self.fail {
                return Err(io::Error::from(io::ErrorKind::StorageFull));
            }
            self.bytes.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn yaml_workspace_preserves_bytes_and_stays_reserved_through_writer_errors() {
        let value = fixture();
        let required = Shape::scan(&value).unwrap().workspace_bytes().unwrap();
        let budget = budget(required);
        for fail in [false, true] {
            let prepared = Prepared::new(&value, Some(&budget)).unwrap();
            let mut writer = CheckedWriter {
                budget: &budget,
                required,
                bytes: Vec::new(),
                fail,
            };
            let result = prepared.write(&mut writer);
            assert_eq!(result.is_err(), fail);
            if !fail {
                assert_eq!(
                    writer.bytes,
                    serde_yaml_ng::to_string(&value).unwrap().as_bytes()
                );
            }
            assert_eq!(budget.snapshot().reserved_memory_bytes, required);
            drop(result);
            drop(prepared);
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert!(budget.snapshot().valid);
        }
    }

    #[test]
    fn yaml_workspace_refusal_precedes_serializer_and_preserves_other_owners() {
        let value = fixture();
        let required = Shape::scan(&value).unwrap().workspace_bytes().unwrap();
        let budget = budget(required + 16);
        let other = budget.reserve(MemoryClass::RetainedResults, 17).unwrap();
        let err = Prepared::new(&value, Some(&budget)).err().unwrap();
        assert!(err.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(budget.snapshot().peak_reserved_memory_bytes, 17);
        drop(other);
        assert!(Prepared::new(&value, Some(&budget)).is_ok());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }

    #[test]
    fn yaml_workspace_releases_reservation_on_writer_panic() {
        struct PanicWriter;
        impl Write for PanicWriter {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                panic!("test writer panic")
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let budget = budget(1 << 20);
        let result = std::panic::catch_unwind(|| {
            Prepared::new(&fixture(), Some(&budget))
                .unwrap()
                .write(&mut PanicWriter)
                .unwrap();
        });
        assert!(result.is_err());
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn yaml_preflight_rejects_unaudited_forms_before_running_formatters() {
        #[derive(Serialize)]
        enum Tagged {
            Newtype(u32),
            Tuple(u32, u32),
            Struct { value: u32 },
        }
        for value in [
            Tagged::Newtype(1),
            Tagged::Tuple(1, 2),
            Tagged::Struct { value: 3 },
        ] {
            assert!(Shape::scan(&value).is_err());
        }
        struct Formatted<'a>(&'a AtomicBool);
        impl fmt::Display for Formatted<'_> {
            fn fmt(&self, _: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.store(true, Ordering::SeqCst);
                panic!("formatter must not run during preflight")
            }
        }
        impl Serialize for Formatted<'_> {
            fn serialize<S: ser::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                serializer.collect_str(self)
            }
        }
        let called = AtomicBool::new(false);
        assert!(Shape::scan(&Formatted(&called)).is_err());
        assert!(!called.load(Ordering::SeqCst));
        struct Bytes;
        impl Serialize for Bytes {
            fn serialize<S: ser::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                serializer.serialize_bytes(b"abc")
            }
        }
        assert!(Shape::scan(&Bytes).is_err());
    }

    #[test]
    fn yaml_preflight_checks_scalar_cast_and_workspace_arithmetic() {
        assert_eq!(Shape::scan(&i128::MIN).unwrap().max_scalar_bytes, 40);
        let mut shape = Shape::default();
        shape.scalar(i32::MAX as usize).unwrap();
        assert!(shape.scalar(i32::MAX as usize + 1).is_err());
        assert!(
            Shape {
                max_scalar_bytes: usize::MAX,
                ..Shape::default()
            }
            .workspace_bytes()
            .is_err()
        );
        assert!(
            Shape {
                max_depth: usize::MAX,
                ..Shape::default()
            }
            .workspace_bytes()
            .is_err()
        );
    }

    #[test]
    fn yaml_source_bound_covers_isolated_allocator_peaks_without_scan_allocations() {
        const ENV: &str = "SHARKMER_YAML_ALLOCATION_AUDIT_CHILD";
        if std::env::var_os(ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "stats::yaml::tests::yaml_source_bound_covers_isolated_allocator_peaks_without_scan_allocations", "--nocapture", "--test-threads=1"])
                .env(ENV, "1")
                .status()
                .unwrap();
            assert!(status.success());
            return;
        }
        // A dedicated child runs only this test, so other tests cannot reset the
        // process-global peak or add unrelated allocations to these measurements.
        let texts = [
            "short".to_owned(),
            "β🧬".repeat(65_536),
            format!("-0x{}", "f".repeat(262_144)),
            "line\n".repeat(65_536),
        ];
        for text in &texts {
            for depth in [0, 14, 15, 16, 17, 30, 31, 32, 33, 65] {
                let value = Nested {
                    depth,
                    text,
                    width: 32,
                };
                let baseline = crate::PEAK_ALLOC.current_usage();
                crate::PEAK_ALLOC.reset_peak_usage();
                let shape = Shape::scan(&value).unwrap();
                assert_eq!(
                    crate::PEAK_ALLOC.peak_usage(),
                    baseline,
                    "preflight allocated"
                );
                let required = shape.workspace_bytes().unwrap();
                let budget = budget(required);
                let prepared = Prepared::new(&value, Some(&budget)).unwrap();
                let baseline = crate::PEAK_ALLOC.current_usage();
                crate::PEAK_ALLOC.reset_peak_usage();
                prepared.write(&mut io::sink()).unwrap();
                let peak = crate::PEAK_ALLOC.peak_usage() - baseline;
                let remaining = crate::PEAK_ALLOC.current_usage();
                assert_eq!(remaining, baseline, "serializer leaked heap storage");
                assert!(
                    peak as u64 <= required,
                    "depth={depth}, scalar={}, peak={peak}, reserved={required}",
                    text.len()
                );
                println!(
                    "yaml_allocator depth={depth} scalar_bytes={} peak_bytes={peak} reserved_bytes={required}",
                    text.len()
                );
                drop(prepared);
                assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            }
        }
    }
}

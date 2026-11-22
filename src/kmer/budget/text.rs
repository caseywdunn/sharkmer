//! UTF-8 text with a reservation that outlives its backing allocation.
use super::collections::Charge;
use super::{MemoryClass, ResourceBudget};
use anyhow::{Context, Result};
use std::sync::Arc;

/// Only UTF-8-preserving writes are exposed; every growth is precharged.
pub(crate) struct BudgetString {
    data: String,
    charge: Charge,
}

impl BudgetString {
    /// Start empty in the requested memory class.
    pub(crate) fn new(budget: Option<&Arc<ResourceBudget>>, class: MemoryClass) -> Self {
        Self {
            data: String::new(),
            charge: Charge::in_class(budget, class),
        }
    }

    /// Reserve text backing before a later operation with a known output length.
    pub(crate) fn with_capacity(
        capacity: usize,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let mut result = Self::new(budget, class);
        result.reserve_capacity(capacity)?;
        Ok(result)
    }

    /// Borrow text without exposing growing mutable String operations.
    pub(crate) fn as_str(&self) -> &str {
        &self.data
    }

    /// Borrow the ledger for metadata derived from this text.
    pub(crate) fn budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.charge.budget()
    }

    /// Encode digest bytes without per-byte temporary strings.
    pub(crate) fn hex(
        bytes: &[u8],
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let mut text = Self::new(budget, class);
        text.push_hex(bytes)?;
        Ok(text)
    }

    /// Append digest text, reusing previously admitted backing when available.
    pub(crate) fn push_hex(&mut self, bytes: &[u8]) -> Result<()> {
        for byte in bytes {
            for digit in [byte >> 4, byte & 15] {
                self.push_char(char::from(b"0123456789abcdef"[usize::from(digit)]))?;
            }
        }
        Ok(())
    }

    /// Independently copy text under its original ledger and memory class.
    pub(crate) fn try_clone(&self) -> Result<Self> {
        Self::from_str(&self.data, self.charge.budget(), self.charge.class())
    }

    /// Copy text into an independently charged allocation.
    pub(crate) fn from_str(
        value: &str,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        let mut result = Self::new(budget, class);
        result.push_str(value)?;
        Ok(result)
    }

    /// Append only after reserving old/replacement overlap.
    pub(crate) fn push_str(&mut self, value: &str) -> Result<()> {
        let needed = self
            .data
            .len()
            .checked_add(value.len())
            .context("Text length overflow")?;
        self.reserve_capacity(needed)?;
        self.data.push_str(value);
        Ok(())
    }

    fn reserve_capacity(&mut self, needed: usize) -> Result<()> {
        if needed > self.data.capacity() {
            let bound = needed.max(
                self.data
                    .capacity()
                    .checked_mul(2)
                    .context("Text capacity overflow")?,
            );
            let reservation = self.charge.reserve(u64::try_from(bound)?)?;
            self.data
                .try_reserve_exact(bound - self.data.len())
                .context("Unable to allocate text")?;
            self.charge.replace(
                reservation,
                u64::try_from(self.data.capacity())?,
                u64::try_from(bound)?,
            )?;
        }
        Ok(())
    }

    /// Append UTF-8 using fixed stack scratch.
    pub(crate) fn push_char(&mut self, value: char) -> Result<()> {
        let mut encoded = [0; 4];
        self.push_str(value.encode_utf8(&mut encoded))
    }

    /// Stream formatting into owned text, preserving the actual reservation
    /// error through fmt's smaller error type. No intermediate String is made.
    pub(crate) fn format(
        args: std::fmt::Arguments<'_>,
        budget: Option<&Arc<ResourceBudget>>,
        class: MemoryClass,
    ) -> Result<Self> {
        struct Writer<'a> {
            text: &'a mut BudgetString,
            error: Option<anyhow::Error>,
        }
        impl std::fmt::Write for Writer<'_> {
            fn write_str(&mut self, value: &str) -> std::fmt::Result {
                match self.text.push_str(value) {
                    Ok(()) => Ok(()),
                    Err(error) => {
                        self.error = Some(error);
                        Err(std::fmt::Error)
                    }
                }
            }
        }
        let mut result = Self::new(budget, class);
        let mut writer = Writer {
            text: &mut result,
            error: None,
        };
        if std::fmt::write(&mut writer, args).is_err() {
            return Err(writer
                .error
                .unwrap_or_else(|| anyhow::anyhow!("Unable to format text")));
        }
        Ok(result)
    }
}

impl std::ops::Deref for BudgetString {
    type Target = str;
    fn deref(&self) -> &str {
        &self.data
    }
}
impl AsRef<str> for BudgetString {
    fn as_ref(&self) -> &str {
        &self.data
    }
}
impl std::fmt::Display for BudgetString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}
impl std::fmt::Debug for BudgetString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}
impl PartialEq for BudgetString {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}
impl serde::Serialize for BudgetString {
    fn serialize<S: serde::Serializer>(
        &self,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.data)
    }
}

impl Eq for BudgetString {}
impl std::hash::Hash for BudgetString {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        std::hash::Hash::hash(self.as_str(), state);
    }
}
impl std::borrow::Borrow<str> for BudgetString {
    fn borrow(&self) -> &str {
        self.as_str()
    }
}
impl Ord for BudgetString {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_str().cmp(other.as_str())
    }
}
impl PartialOrd for BudgetString {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl super::collections::TryClone for BudgetString {
    fn try_clone(&self) -> Result<Self> {
        BudgetString::try_clone(self)
    }
}
#[cfg(test)]
impl Clone for BudgetString {
    fn clone(&self) -> Self {
        self.try_clone().unwrap()
    }
}
#[cfg(test)]
impl From<String> for BudgetString {
    fn from(value: String) -> Self {
        Self::from_str(&value, None, MemoryClass::RetainedResults).unwrap()
    }
}
#[cfg(test)]
impl From<&str> for BudgetString {
    fn from(value: &str) -> Self {
        Self::from_str(value, None, MemoryClass::RetainedResults).unwrap()
    }
}
impl PartialEq<str> for BudgetString {
    fn eq(&self, other: &str) -> bool {
        self.as_str() == other
    }
}
impl PartialEq<&str> for BudgetString {
    fn eq(&self, other: &&str) -> bool {
        self.as_str() == *other
    }
}

impl PartialEq<BudgetString> for &str {
    fn eq(&self, other: &BudgetString) -> bool {
        *self == other.as_str()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceKind, ResourceLimits};
    fn budget(bytes: u64) -> Arc<ResourceBudget> {
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = bytes;
        ResourceBudget::new(ResourceLimits {
            memory_bytes: bytes,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap()
    }
    #[test]
    fn text_growth_refusal_preserves_utf8_and_retries_after_competing_owner_drops() {
        let resources = budget(10);
        let mut text =
            BudgetString::from_str("é", Some(&resources), MemoryClass::RetainedResults).unwrap();
        let held = resources.reserve(MemoryClass::RetainedResults, 6).unwrap();
        let before = resources.snapshot();
        let error = text.push_str("abcdef").unwrap_err();
        assert_eq!(
            error.downcast_ref::<ResourceExhausted>().unwrap().resource,
            ResourceKind::StageMemory(MemoryClass::RetainedResults)
        );
        assert_eq!(&*text, "é");
        assert_eq!(
            resources.snapshot().reserved_memory_bytes,
            before.reserved_memory_bytes
        );
        assert_eq!(
            resources.snapshot().peak_reserved_memory_bytes,
            before.peak_reserved_memory_bytes
        );
        drop(held);
        text.push_str("abcdef").unwrap();
        assert_eq!(&*text, "éabcdef");
        assert_eq!(
            resources.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        drop(text);
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
    }
    #[test]
    fn streamed_formatting_matches_string_serialization_and_preserves_typed_failures() {
        let resources = budget(4096);
        let text = BudgetString::format(
            format_args!(
                "sample={} mean={:.2} count={}",
                "échantillon",
                1.234,
                u64::MAX
            ),
            Some(&resources),
            MemoryClass::RetainedResults,
        )
        .unwrap();
        let expected = format!(
            "sample={} mean={:.2} count={}",
            "échantillon",
            1.234,
            u64::MAX
        );
        assert_eq!(&*text, expected);
        assert_eq!(
            serde_yaml_ng::to_string(&text).unwrap(),
            serde_yaml_ng::to_string(&expected).unwrap()
        );
        drop(text);
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
        let tiny = budget(1);
        let error = BudgetString::format(
            format_args!("{} {}", "too long", u64::MAX),
            Some(&tiny),
            MemoryClass::RetainedResults,
        )
        .unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
    }
    #[test]
    fn owned_hex_matches_all_byte_values_and_releases_failed_encoding() {
        let resources = budget(4096);
        let values: Vec<u8> = (0..=255).collect();
        let text =
            BudgetString::hex(&values, Some(&resources), MemoryClass::RetainedResults).unwrap();
        let expected: String = values.iter().map(|byte| format!("{byte:02x}")).collect();
        assert_eq!(text.as_str(), expected);
        assert_eq!(text.len(), 512);
        assert_eq!(
            resources.snapshot().reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        drop(text);
        assert_eq!(resources.snapshot().reserved_memory_bytes, 0);
        let tiny = budget(3);
        let error =
            BudgetString::hex(&values, Some(&tiny), MemoryClass::RetainedResults).unwrap_err();
        assert!(error.downcast_ref::<ResourceExhausted>().is_some());
        assert_eq!(tiny.snapshot().reserved_memory_bytes, 0);
    }
}

//! Admission for owned configuration models, preserving their Serde rules.
//!
//! The audited models contain Strings, Vecs and BTreeMaps, with scalar defaults
//! and inline structs/options. New heap-owning field types or allocating
//! defaults require a corresponding admission rule and allocator check.
use crate::kmer::budget::{MemoryClass, MemoryReservation, ResourceBudget};
use anyhow::{Context, Result, ensure};
use serde::de::{self, DeserializeOwned, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use std::fmt;
use std::sync::Arc;

mod tracking;

/// The data drops before the aggregate backing-allocation reservation.
pub(crate) struct Owned<T> {
    data: T,
    owner: Owner,
}

impl<T> std::ops::Deref for Owned<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.data
    }
}

impl<T> Owned<T> {
    /// Transfer data and its still-live owner together into a collected model.
    pub(crate) fn into_parts(self) -> (T, Owner) {
        (self.data, self.owner)
    }
}

/// Admission shared by a parsed model and subsequent deterministic transforms.
pub(crate) struct Owner {
    reservation: Option<MemoryReservation>,
    bytes: u64,
    failure: Option<anyhow::Error>,
}

impl Owner {
    /// Begin without any backing allocations.
    pub(crate) fn new(budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        Ok(Self {
            reservation: budget
                .map(|b| b.reserve(MemoryClass::RetainedResults, 0))
                .transpose()?,
            bytes: 0,
            failure: None,
        })
    }

    fn resize(&mut self, bytes: u64) -> Result<()> {
        if let Some(reservation) = &mut self.reservation {
            reservation.resize(bytes)?;
        }
        self.bytes = bytes;
        Ok(())
    }

    fn admit(&mut self, bytes: u64) -> Result<()> {
        self.resize(
            self.bytes
                .checked_add(bytes)
                .context("Configuration allocation overflow")?,
        )
    }

    fn record<E: de::Error>(&mut self, error: anyhow::Error) -> E {
        if self.failure.is_none() {
            self.failure = Some(error);
        }
        E::custom("Unable to admit configuration data")
    }

    fn serde_admit<E: de::Error>(&mut self, bytes: usize) -> std::result::Result<(), E> {
        u64::try_from(bytes)
            .context("Configuration allocation exceeds supported range")
            .and_then(|bytes| self.admit(bytes))
            .map_err(|error| self.record(error))
    }

    /// Format deterministic strings/numbers without an unadmitted temporary.
    pub(crate) fn text(&mut self, args: fmt::Arguments<'_>) -> Result<String> {
        struct Count(usize);
        impl fmt::Write for Count {
            fn write_str(&mut self, text: &str) -> fmt::Result {
                self.0 = self.0.checked_add(text.len()).ok_or(fmt::Error)?;
                Ok(())
            }
        }
        struct Write<'a> {
            text: &'a mut String,
            limit: usize,
        }
        impl fmt::Write for Write<'_> {
            fn write_str(&mut self, text: &str) -> fmt::Result {
                if self
                    .text
                    .len()
                    .checked_add(text.len())
                    .is_none_or(|n| n > self.limit)
                {
                    return Err(fmt::Error);
                }
                self.text.push_str(text);
                Ok(())
            }
        }
        let mut length = Count(0);
        fmt::write(&mut length, args).context("Configuration text length overflow")?;
        self.admit(u64::try_from(length.0)?)?;
        let mut text = String::new();
        text.try_reserve_exact(length.0)
            .context("Unable to allocate configuration text")?;
        ensure!(
            text.capacity() <= length.0,
            "Configuration text exceeded admitted capacity"
        );
        fmt::write(
            &mut Write {
                text: &mut text,
                limit: length.0,
            },
            args,
        )
        .context("Configuration formatting changed between passes")?;
        Ok(text)
    }

    /// Preserve Unicode case conversion without intermediate owned text.
    pub(crate) fn case(&mut self, value: &str, uppercase: bool) -> Result<String> {
        struct Case<'a> {
            value: &'a str,
            uppercase: bool,
        }
        impl fmt::Display for Case<'_> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                use fmt::Write;
                for character in self.value.chars() {
                    if self.uppercase {
                        for character in character.to_uppercase() {
                            f.write_char(character)?;
                        }
                    } else {
                        for character in character.to_lowercase() {
                            f.write_char(character)?;
                        }
                    }
                }
                Ok(())
            }
        }
        self.text(format_args!("{}", Case { value, uppercase }))
    }

    /// Release discarded parser/model allocations after their data has dropped.
    /// The caller inventories every surviving owned String/Vec capacity.
    pub(crate) fn retain_bytes(&mut self, bytes: u64) -> Result<()> {
        ensure!(
            bytes <= self.bytes,
            "Configuration retained data exceeded admission"
        );
        self.resize(bytes)
    }

    /// Retain a model produced using this owner.
    pub(crate) fn finish<T>(self, data: T) -> Owned<T> {
        Owned { data, owner: self }
    }
}

/// Deserialize borrowed input while precharging loader and retained model data.
pub(crate) fn yaml<T: DeserializeOwned>(
    input: &[u8],
    budget: Option<&Arc<ResourceBudget>>,
) -> Result<Owned<T>> {
    let workspace = crate::output::parsing::workspace::Workspace::new(input, budget)?;
    let mut owner = Owner::new(budget)?;
    let result = T::deserialize(tracking::Deserializer {
        inner: serde_yaml_ng::Deserializer::from_slice(input),
        owner: &mut owner,
    });
    let data = match result {
        Ok(data) => data,
        Err(error) => {
            if let Some(failure) = owner.failure.take() {
                drop(error);
                return Err(failure);
            }
            return Err(workspace.retain_error(error));
        }
    };
    Ok(owner.finish(data))
}

#[cfg(test)]
mod tests;

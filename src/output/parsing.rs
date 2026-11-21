//! Owned historical YAML input, typed data and precharged pinned-loader workspace.
//! A model's aggregate reservation outlives all of its original String/Vec data;
//! admission happens in Serde seeds before each allocation, including aliases.
use super::recovery::RecoveryJournal;
use super::*;
use crate::kmer::budget::MemoryReservation;
use serde::de::{self, DeserializeSeed, Error as _, MapAccess, SeqAccess, Visitor};
use std::fmt;
use std::marker::PhantomData;

pub(crate) mod workspace;

/// Keep typed data ahead of its backing-allocation owner in destruction order.
pub(super) struct Parsed<T> {
    data: T,
    _owner: Owner,
}
impl<T> std::ops::Deref for Parsed<T> {
    type Target = T;
    fn deref(&self) -> &T {
        &self.data
    }
}

struct Owner {
    reservation: Option<MemoryReservation>,
    bytes: u64,
    failure: Option<anyhow::Error>,
}
impl Owner {
    fn new(budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
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
                .context("Parsed YAML allocation overflow")?,
        )
    }
    fn observe(&mut self, old: u64, bound: u64, actual: u64) -> Result<()> {
        ensure!(
            actual <= bound,
            "Parsed YAML allocation exceeded its admitted bound"
        );
        let bytes = self
            .bytes
            .checked_sub(old)
            .and_then(|n| n.checked_sub(bound))
            .and_then(|n| n.checked_add(actual))
            .context("Parsed YAML accounting overflow")?;
        self.resize(bytes)
    }
    fn text(&mut self, value: &str) -> Result<String> {
        let bound = u64::try_from(value.len())?;
        self.admit(bound)
            .context("Unable to admit parsed YAML text")?;
        let mut text = String::new();
        text.try_reserve_exact(value.len())
            .context("Unable to allocate parsed YAML text")?;
        text.push_str(value);
        self.observe(0, bound, u64::try_from(text.capacity())?)?;
        Ok(text)
    }
    fn push<T>(&mut self, values: &mut Vec<T>, value: T) -> Result<()> {
        if values.len() == values.capacity() {
            let capacity = values
                .capacity()
                .checked_mul(2)
                .context("Parsed YAML array capacity overflow")?
                .max(4);
            let bound = crate::kmer::budget::collections::bytes::<T>(capacity)?;
            let old = crate::kmer::budget::collections::bytes::<T>(values.capacity())?;
            self.admit(bound)?;
            values
                .try_reserve_exact(
                    capacity
                        .checked_sub(values.len())
                        .context("Invalid parsed YAML array")?,
                )
                .context("Unable to allocate parsed YAML receipt array")?;
            self.observe(
                old,
                bound,
                crate::kmer::budget::collections::bytes::<T>(values.capacity())?,
            )?;
        }
        values.push(value);
        Ok(())
    }
    fn record<E: de::Error>(&mut self, error: anyhow::Error) -> E {
        // Serde errors cannot carry the typed budget error. Keep the original
        // beside the owner and recover it at the concrete parser boundary.
        let message = E::custom("Unable to admit parsed YAML data");
        if self.failure.is_none() {
            self.failure = Some(error);
        }
        message
    }
}

struct Seed<'a, T> {
    owner: &'a mut Owner,
    marker: PhantomData<T>,
}
impl<'a, T> Seed<'a, T> {
    fn new(owner: &'a mut Owner) -> Self {
        Self {
            owner,
            marker: PhantomData,
        }
    }
}
trait OwnedDeserialize: Sized {
    fn read<'de, D: de::Deserializer<'de>>(
        owner: &mut Owner,
        de: D,
    ) -> std::result::Result<Self, D::Error>;
    fn missing<E: de::Error>(name: &'static str) -> std::result::Result<Self, E> {
        Err(E::missing_field(name))
    }
}
impl<'de, T: OwnedDeserialize> DeserializeSeed<'de> for Seed<'_, T> {
    type Value = T;
    fn deserialize<D: de::Deserializer<'de>>(self, de: D) -> std::result::Result<T, D::Error> {
        T::read(self.owner, de)
    }
}
impl OwnedDeserialize for u32 {
    fn read<'de, D: de::Deserializer<'de>>(
        _: &mut Owner,
        de: D,
    ) -> std::result::Result<Self, D::Error> {
        Self::deserialize(de)
    }
}
impl OwnedDeserialize for String {
    fn read<'de, D: de::Deserializer<'de>>(
        owner: &mut Owner,
        de: D,
    ) -> std::result::Result<Self, D::Error> {
        struct Text<'a>(&'a mut Owner);
        impl Visitor<'_> for Text<'_> {
            type Value = String;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("a string")
            }
            fn visit_str<E: de::Error>(self, value: &str) -> std::result::Result<String, E> {
                self.0.text(value).map_err(|error| self.0.record(error))
            }
        }
        de.deserialize_string(Text(owner))
    }
}
impl<T: OwnedDeserialize> OwnedDeserialize for Option<T> {
    fn missing<E: de::Error>(_: &'static str) -> std::result::Result<Self, E> {
        Ok(None)
    }
    fn read<'de, D: de::Deserializer<'de>>(
        owner: &mut Owner,
        de: D,
    ) -> std::result::Result<Self, D::Error> {
        struct Optional<'a, T>(Seed<'a, T>);
        impl<'de, T: OwnedDeserialize> Visitor<'de> for Optional<'_, T> {
            type Value = Option<T>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("an optional value")
            }
            fn visit_none<E: de::Error>(self) -> std::result::Result<Self::Value, E> {
                Ok(None)
            }
            fn visit_unit<E: de::Error>(self) -> std::result::Result<Self::Value, E> {
                Ok(None)
            }
            fn visit_some<D: de::Deserializer<'de>>(
                self,
                de: D,
            ) -> std::result::Result<Self::Value, D::Error> {
                self.0.deserialize(de).map(Some)
            }
        }
        de.deserialize_option(Optional(Seed::<T>::new(owner)))
    }
}
impl<T: OwnedDeserialize> OwnedDeserialize for Vec<T> {
    fn read<'de, D: de::Deserializer<'de>>(
        owner: &mut Owner,
        de: D,
    ) -> std::result::Result<Self, D::Error> {
        struct Array<'a, T>(Seed<'a, T>);
        impl<'de, T: OwnedDeserialize> Visitor<'de> for Array<'_, T> {
            type Value = Vec<T>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("a sequence")
            }
            fn visit_seq<A: SeqAccess<'de>>(
                self,
                mut access: A,
            ) -> std::result::Result<Self::Value, A::Error> {
                let mut values = Vec::new();
                while let Some(value) = access.next_element_seed(Seed::<T>::new(self.0.owner))? {
                    self.0
                        .owner
                        .push(&mut values, value)
                        .map_err(|error| self.0.owner.record(error))?;
                }
                Ok(values)
            }
        }
        de.deserialize_seq(Array(Seed::<T>::new(owner)))
    }
}

// Field identifiers resolve to static names. Escaped or aliased YAML keys do
// not create another owned String before the value seed is invoked.
struct Key {
    names: &'static [&'static str],
    deny_unknown: bool,
}
impl Visitor<'_> for Key {
    type Value = Option<&'static str>;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("a field identifier")
    }
    fn visit_str<E: de::Error>(self, name: &str) -> std::result::Result<Self::Value, E> {
        match self
            .names
            .iter()
            .copied()
            .find(|expected| *expected == name)
        {
            Some(name) => Ok(Some(name)),
            None if self.deny_unknown => Err(E::unknown_field(name, self.names)),
            None => Ok(None),
        }
    }
}

impl<'de> DeserializeSeed<'de> for Key {
    type Value = Option<&'static str>;
    fn deserialize<D: de::Deserializer<'de>>(
        self,
        de: D,
    ) -> std::result::Result<Self::Value, D::Error> {
        de.deserialize_identifier(self)
    }
}

macro_rules! owned_model {
    ($model:ident, $deny:literal, { $($field:ident : $ty:ty),+ $(,)? }) => {
        impl OwnedDeserialize for $model {
            fn read<'de, D: de::Deserializer<'de>>(owner: &mut Owner, de: D) -> std::result::Result<Self, D::Error> {
                const FIELDS: &[&str] = &[$(stringify!($field)),+];
                struct Fields<'a>(&'a mut Owner);
                impl<'de> Visitor<'de> for Fields<'_> {
                    type Value = $model;
                    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result { f.write_str(stringify!($model)) }
                    fn visit_map<A: MapAccess<'de>>(self, mut access: A) -> std::result::Result<Self::Value, A::Error> {
                        $(let mut $field: Option<$ty> = None;)+
                        while let Some(key) = access.next_key_seed(Key { names: FIELDS, deny_unknown: $deny })? {
                            match key {
                                $(Some(stringify!($field)) => {
                                    if $field.is_some() { return Err(A::Error::duplicate_field(stringify!($field))); }
                                    $field = Some(access.next_value_seed(Seed::<$ty>::new(self.0))?);
                                },)+
                                _ => { access.next_value::<de::IgnoredAny>()?; }
                            }
                        }
                        Ok($model { $($field: match $field {
                            Some(value) => value,
                            None => <$ty as OwnedDeserialize>::missing(stringify!($field))?,
                        }),+ })
                    }
                    fn visit_seq<A: SeqAccess<'de>>(self, mut access: A) -> std::result::Result<Self::Value, A::Error> {
                        Ok($model { $($field: match access.next_element_seed(Seed::<$ty>::new(self.0))? {
                            Some(value) => value,
                            None => <$ty as OwnedDeserialize>::missing(stringify!($field))?,
                        }),+ })
                    }
                }
                de.deserialize_struct(stringify!($model), FIELDS, Fields(owner))
            }
        }
    };
}
owned_model!(OutputReceipt, false, { path: String, sha256: String });
owned_model!(OutputManifest, false, {
    schema_version: u32, producer: String, sample: String, run_id: String,
    status: String, files: Vec<OutputReceipt>, failure_reason: Option<String>, staging_directory: Option<String>,
});
owned_model!(RecoveryJournal, true, {
    schema_version: u32, producer: String, sample: String, run_id: String,
    staging_directory: String, previous_manifest_sha256: String,
    backup_files: Vec<OutputReceipt>, replacement_files: Vec<OutputReceipt>,
});

fn read<T: OwnedDeserialize>(path: &BudgetPath) -> Result<Parsed<T>> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open historical YAML {}", path.display()))?;
    let capacity = usize::try_from(file.metadata()?.len())
        .context("Historical YAML size exceeds this platform")?;
    let mut input =
        BudgetVec::with_capacity_in_class(capacity, path.budget(), MemoryClass::RetainedResults)?;
    let mut buffer = [0u8; 16_384];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        for byte in &buffer[..n] {
            input.push(*byte)?;
        }
    }
    drop(file);
    let workspace = workspace::Workspace::new(&input, path.budget())?;
    let mut owner = Owner::new(path.budget())?;
    // Borrow the admitted input. The loader drops its events, scanner and
    // anchor dictionaries inside this call while the workspace guard lives.
    let result =
        Seed::<T>::new(&mut owner).deserialize(serde_yaml_ng::Deserializer::from_slice(&input));
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
    Ok(Parsed {
        data,
        _owner: owner,
    })
}

pub(super) fn manifest(path: &BudgetPath) -> Result<Parsed<OutputManifest>> {
    read(path)
}
pub(super) fn recovery(path: &BudgetPath) -> Result<Parsed<RecoveryJournal>> {
    read(path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, ResourceExhausted, ResourceLimits};

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
    fn path(root: &Path, budget: &Arc<ResourceBudget>) -> BudgetPath {
        BudgetPath::from_path(root, Some(budget), MemoryClass::RetainedResults).unwrap()
    }
    fn compare<T: OwnedDeserialize + serde::de::DeserializeOwned + Serialize>(text: &str) {
        let directory = tempfile::tempdir().unwrap();
        let file = directory.path().join("metadata.yaml");
        fs::write(&file, text).unwrap();
        let budget = budget(8 << 20);
        let file = path(&file, &budget);
        let baseline = budget.snapshot().reserved_memory_bytes;
        let old = serde_yaml_ng::from_str::<T>(text);
        let new = read::<T>(&file);
        match (old, new) {
            (Ok(old), Ok(new)) => {
                assert_eq!(
                    serde_yaml_ng::to_string(&old).unwrap(),
                    serde_yaml_ng::to_string(&*new).unwrap(),
                    "{text}"
                );
            }
            (Err(_), Err(_)) => {}
            (old, new) => panic!(
                "Parser compatibility differs: old_ok={} new_ok={} input={text}",
                old.is_ok(),
                new.is_ok()
            ),
        }
        assert_eq!(budget.snapshot().reserved_memory_bytes, baseline);
        drop(file);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn owned_model_seeds_preserve_historical_maps_sequences_aliases_and_validation() {
        let digest = "a".repeat(64);
        for schema in 1..=4 {
            let simple = format!(
                "schema_version: {schema}\nproducer: sharkmer\nsample: sample\nrun_id: abc-1\nstatus: complete\nfiles: []\n"
            );
            compare::<OutputManifest>(&format!("{simple}---\n{simple}"));
            for suffix in [
                "",
                "failure_reason: null\nstaging_directory: null\n",
                "extra: {ignored: [one, two, {three: four}]}\n",
            ] {
                compare::<OutputManifest>(&format!("{simple}{suffix}"));
            }
            compare::<OutputManifest>(&format!(
                "[{schema}, sharkmer, sample, abc-1, complete, [], null, null]"
            ));
            compare::<OutputManifest>(&format!(
                "[{schema}, sharkmer, sample, abc-1, complete, []]"
            ));
            compare::<OutputManifest>(&format!(
                "[{schema}, sharkmer, sample, abc-1, complete, [], null, null, extra]"
            ));
        }
        compare::<OutputManifest>(&format!(
            "\"schema_versi\\u006fn\": 0x4\nproducer: sharkmer\nsample: &sample sample\nrun_id: abc\nstatus: complete\nfiles:\n- &receipt {{path: sample_gene.fasta, sha256: &digest {digest}, ignored: null}}\n- *receipt\nextra: {{a: *digest, b: *sample}}\n"
        ));
        compare::<OutputManifest>(
            "schema_version: 4\nproducer: sharkmer\nsample: λ🧬\nrun_id: abc\nstatus: failed\nfiles: []\nfailure_reason: |\n  first line\n  second line\nstaging_directory: .stage\n",
        );
        for invalid in [
            "{}",
            "[]",
            "null",
            "schema_version: 4294967296",
            "schema_version: -1",
            "schema_version: 4\nschema_version: 3",
            "files: &cycle [*cycle]",
            "schema_version: 4\nproducer: []",
            "schema_version: 4\nfiles: {}",
            "schema_version: 4\nfailure_reason: null\nfailure_reason: null",
            "---\n{}\n---\n{}\n",
            "files: *missing",
            "files: [{path: one, path: two}]",
        ] {
            compare::<OutputManifest>(invalid);
        }
        for receipt in [
            "{path: sample_gene.fasta, sha256: aaa}",
            "[sample_gene.fasta, aaa]",
            "{path: λ, sha256: aaa, extra: true}",
            "{path: one, path: two, sha256: aaa}",
            "{path: []}",
        ] {
            compare::<OutputReceipt>(receipt);
        }
        let journal = format!(
            "schema_version: 1\nproducer: sharkmer\nsample: sample\nrun_id: abc\nstaging_directory: .stage\nprevious_manifest_sha256: &hash {digest}\nbackup_files: [{{path: sample.stats.yaml, sha256: *hash}}]\nreplacement_files: []\n"
        );
        compare::<RecoveryJournal>(&journal);
        compare::<RecoveryJournal>(&format!("{journal}unknown: true\n"));
        compare::<RecoveryJournal>(&format!(
            "[1, sharkmer, sample, abc, .stage, {digest}, [], []]"
        ));
    }

    #[test]
    fn alias_payload_admission_preserves_typed_refusal_and_retries_after_release() {
        let directory = tempfile::tempdir().unwrap();
        let file = directory.path().join("aliases.yaml");
        let text = format!(
            "schema_version: 4\nproducer: sharkmer\nsample: sample\nrun_id: abc\nstatus: complete\nfiles:\n- &receipt {{path: sample_gene.fasta, sha256: {}}}\n{}",
            "f".repeat(65536),
            "- *receipt\n".repeat(63)
        );
        fs::write(&file, &text).unwrap();
        let budget = budget(16 << 20);
        let file = path(&file, &budget);
        let baseline = budget.snapshot().reserved_memory_bytes;
        let workspace = workspace::bytes(text.as_bytes()).unwrap();
        let typed_window = 1 << 20;
        let competing = budget
            .reserve(
                MemoryClass::RetainedResults,
                (16 << 20) - baseline - text.len() as u64 - workspace - typed_window,
            )
            .unwrap();
        let before = budget.snapshot().reserved_memory_bytes;
        let error = read::<OutputManifest>(&file).err().unwrap();
        assert!(
            error.downcast_ref::<ResourceExhausted>().is_some(),
            "{error:#}"
        );
        assert!(format!("{error:#}").contains("Unable to admit parsed YAML text"));
        assert_eq!(budget.snapshot().reserved_memory_bytes, before);
        assert_eq!(fs::read(&file).unwrap(), text.as_bytes());
        drop(competing);
        let parsed = read::<OutputManifest>(&file).unwrap();
        assert_eq!(parsed.files.len(), 64);
        for receipt in &parsed.files {
            assert_eq!(receipt.sha256.len(), 65536);
        }
        assert!(budget.snapshot().reserved_memory_bytes > baseline + (4 << 20));
        assert_eq!(
            budget.snapshot().peak_reserved_by_class[MemoryClass::TargetGraphs as usize],
            0
        );
        drop(parsed);
        assert_eq!(budget.snapshot().reserved_memory_bytes, baseline);
        drop(file);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }

    #[test]
    fn parsed_owner_matches_live_model_capacities_and_outlives_destruction() {
        let directory = tempfile::tempdir().unwrap();
        let file = directory.path().join("manifest.yaml");
        let text = format!(
            "schema_version: 4\nproducer: sharkmer\nsample: sample\nrun_id: abc\nstatus: failed\nfiles:\n{}failure_reason: first line\nstaging_directory: .stage\n",
            (0..19)
                .map(|n| format!("- path: sample_{n}.fasta\n  sha256: {}\n", "f".repeat(64)))
                .collect::<String>()
        );
        fs::write(&file, &text).unwrap();
        let budget = budget(2 << 20);
        let file = path(&file, &budget);
        let baseline = budget.snapshot().reserved_memory_bytes;
        let parsed = read::<OutputManifest>(&file).unwrap();
        let expected = [
            &parsed.producer,
            &parsed.sample,
            &parsed.run_id,
            &parsed.status,
        ]
        .into_iter()
        .map(|s| s.capacity())
        .sum::<usize>()
            + parsed.failure_reason.as_ref().map_or(0, String::capacity)
            + parsed
                .staging_directory
                .as_ref()
                .map_or(0, String::capacity)
            + parsed.files.capacity() * size_of::<OutputReceipt>()
            + parsed
                .files
                .iter()
                .map(|r| r.path.capacity() + r.sha256.capacity())
                .sum::<usize>();
        assert_eq!(parsed._owner.bytes, expected as u64);
        assert_eq!(
            budget.snapshot().reserved_memory_bytes,
            baseline + expected as u64
        );
        struct Dropping {
            data: OutputManifest,
            budget: Arc<ResourceBudget>,
            required: u64,
        }
        impl Drop for Dropping {
            fn drop(&mut self) {
                assert!(self.budget.snapshot().reserved_memory_bytes >= self.required);
                assert_eq!(self.data.files.len(), 19);
            }
        }
        let Parsed { data, _owner } = parsed;
        let wrapped = Parsed {
            data: Dropping {
                data,
                budget: Arc::clone(&budget),
                required: baseline + expected as u64,
            },
            _owner,
        };
        drop(wrapped);
        assert_eq!(budget.snapshot().reserved_memory_bytes, baseline);
        drop(file);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
        assert!(budget.snapshot().valid);
    }
}

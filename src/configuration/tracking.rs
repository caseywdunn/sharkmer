//! Wrap Serde visitors; field interpretation and validation stay in the derives.
use super::*;

pub(super) struct Deserializer<'a, D> {
    pub(super) inner: D,
    pub(super) owner: &'a mut Owner,
}

#[derive(Clone, Copy)]
enum Storage {
    Inline,
    Text,
    Sequence,
    Map,
    Dynamic,
}

struct TrackedVisitor<'a, V> {
    inner: V,
    owner: &'a mut Owner,
    storage: Storage,
}

macro_rules! dispatch {
    ($method:ident, $storage:ident $(, $arg:ident : $ty:ty)*) => {
        fn $method<V: Visitor<'de>>(self, $($arg: $ty,)* visitor: V) -> std::result::Result<V::Value, D::Error> {
            self.inner.$method($($arg,)* TrackedVisitor { inner: visitor, owner: self.owner, storage: Storage::$storage })
        }
    };
}
impl<'de, D: de::Deserializer<'de>> de::Deserializer<'de> for Deserializer<'_, D> {
    type Error = D::Error;
    dispatch!(deserialize_any, Dynamic);
    dispatch!(deserialize_bool, Inline);
    dispatch!(deserialize_i8, Inline);
    dispatch!(deserialize_i16, Inline);
    dispatch!(deserialize_i32, Inline);
    dispatch!(deserialize_i64, Inline);
    dispatch!(deserialize_i128, Inline);
    dispatch!(deserialize_u8, Inline);
    dispatch!(deserialize_u16, Inline);
    dispatch!(deserialize_u32, Inline);
    dispatch!(deserialize_u64, Inline);
    dispatch!(deserialize_u128, Inline);
    dispatch!(deserialize_f32, Inline);
    dispatch!(deserialize_f64, Inline);
    dispatch!(deserialize_char, Inline);
    dispatch!(deserialize_str, Text);
    dispatch!(deserialize_string, Text);
    dispatch!(deserialize_bytes, Text);
    dispatch!(deserialize_byte_buf, Text);
    dispatch!(deserialize_option, Inline);
    dispatch!(deserialize_unit, Inline);
    dispatch!(deserialize_unit_struct, Inline, name: &'static str);
    dispatch!(deserialize_newtype_struct, Inline, name: &'static str);
    dispatch!(deserialize_seq, Sequence);
    dispatch!(deserialize_tuple, Inline, len: usize);
    dispatch!(deserialize_tuple_struct, Inline, name: &'static str, len: usize);
    dispatch!(deserialize_map, Map);
    dispatch!(deserialize_struct, Inline, name: &'static str, fields: &'static [&'static str]);
    dispatch!(deserialize_enum, Inline, name: &'static str, variants: &'static [&'static str]);
    dispatch!(deserialize_identifier, Inline);
    dispatch!(deserialize_ignored_any, Inline);
    fn is_human_readable(&self) -> bool {
        self.inner.is_human_readable()
    }
}

struct Seed<'a, T> {
    inner: T,
    owner: &'a mut Owner,
}
impl<'de, T: DeserializeSeed<'de>> DeserializeSeed<'de> for Seed<'_, T> {
    type Value = T::Value;
    fn deserialize<D: de::Deserializer<'de>>(
        self,
        inner: D,
    ) -> std::result::Result<Self::Value, D::Error> {
        self.inner.deserialize(Deserializer {
            inner,
            owner: self.owner,
        })
    }
}

macro_rules! primitive {
    ($method:ident, $ty:ty) => {
        fn $method<E: de::Error>(self, value: $ty) -> std::result::Result<V::Value, E> {
            self.inner.$method(value)
        }
    };
}
impl<'de, V: Visitor<'de>> Visitor<'de> for TrackedVisitor<'_, V> {
    type Value = V::Value;
    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.inner.expecting(f)
    }
    primitive!(visit_bool, bool);
    primitive!(visit_i8, i8);
    primitive!(visit_i16, i16);
    primitive!(visit_i32, i32);
    primitive!(visit_i64, i64);
    primitive!(visit_i128, i128);
    primitive!(visit_u8, u8);
    primitive!(visit_u16, u16);
    primitive!(visit_u32, u32);
    primitive!(visit_u64, u64);
    primitive!(visit_u128, u128);
    primitive!(visit_f32, f32);
    primitive!(visit_f64, f64);
    primitive!(visit_char, char);
    fn visit_str<E: de::Error>(self, value: &str) -> std::result::Result<V::Value, E> {
        if matches!(self.storage, Storage::Text | Storage::Dynamic) {
            self.owner.serde_admit(value.len())?;
        }
        self.inner.visit_str(value)
    }
    fn visit_borrowed_str<E: de::Error>(self, value: &'de str) -> std::result::Result<V::Value, E> {
        if matches!(self.storage, Storage::Text | Storage::Dynamic) {
            self.owner.serde_admit(value.len())?;
        }
        self.inner.visit_borrowed_str(value)
    }
    fn visit_string<E: de::Error>(self, value: String) -> std::result::Result<V::Value, E> {
        self.owner.serde_admit(value.capacity())?;
        self.inner.visit_string(value)
    }
    fn visit_bytes<E: de::Error>(self, value: &[u8]) -> std::result::Result<V::Value, E> {
        if matches!(self.storage, Storage::Text | Storage::Dynamic) {
            self.owner.serde_admit(value.len())?;
        }
        self.inner.visit_bytes(value)
    }
    fn visit_borrowed_bytes<E: de::Error>(
        self,
        value: &'de [u8],
    ) -> std::result::Result<V::Value, E> {
        if matches!(self.storage, Storage::Text | Storage::Dynamic) {
            self.owner.serde_admit(value.len())?;
        }
        self.inner.visit_borrowed_bytes(value)
    }
    fn visit_byte_buf<E: de::Error>(self, value: Vec<u8>) -> std::result::Result<V::Value, E> {
        self.owner.serde_admit(value.capacity())?;
        self.inner.visit_byte_buf(value)
    }
    fn visit_none<E: de::Error>(self) -> std::result::Result<V::Value, E> {
        self.inner.visit_none()
    }
    fn visit_unit<E: de::Error>(self) -> std::result::Result<V::Value, E> {
        self.inner.visit_unit()
    }
    fn visit_some<D: de::Deserializer<'de>>(
        self,
        inner: D,
    ) -> std::result::Result<V::Value, D::Error> {
        self.inner.visit_some(Deserializer {
            inner,
            owner: self.owner,
        })
    }
    fn visit_newtype_struct<D: de::Deserializer<'de>>(
        self,
        inner: D,
    ) -> std::result::Result<V::Value, D::Error> {
        self.inner.visit_newtype_struct(Deserializer {
            inner,
            owner: self.owner,
        })
    }
    fn visit_seq<A: SeqAccess<'de>>(self, inner: A) -> std::result::Result<V::Value, A::Error> {
        self.inner.visit_seq(Sequence {
            inner,
            owner: self.owner,
            length: 0,
            capacity: 0,
            growing: matches!(self.storage, Storage::Sequence | Storage::Dynamic),
        })
    }
    fn visit_map<A: MapAccess<'de>>(self, inner: A) -> std::result::Result<V::Value, A::Error> {
        self.inner.visit_map(Map {
            inner,
            owner: self.owner,
            key_bytes: 0,
            growing: matches!(self.storage, Storage::Map | Storage::Dynamic),
        })
    }
    fn visit_enum<A: de::EnumAccess<'de>>(
        self,
        inner: A,
    ) -> std::result::Result<V::Value, A::Error> {
        self.inner.visit_enum(Enum {
            inner,
            owner: self.owner,
        })
    }
}

struct Sequence<'a, A> {
    inner: A,
    owner: &'a mut Owner,
    length: usize,
    capacity: usize,
    growing: bool,
}
impl<'de, A: SeqAccess<'de>> SeqAccess<'de> for Sequence<'_, A> {
    type Error = A::Error;
    fn next_element_seed<T: DeserializeSeed<'de>>(
        &mut self,
        seed: T,
    ) -> std::result::Result<Option<T::Value>, A::Error> {
        if self.growing && self.length == self.capacity {
            let capacity = self
                .capacity
                .checked_mul(2)
                .ok_or_else(|| {
                    self.owner
                        .record(anyhow::anyhow!("Configuration sequence capacity overflow"))
                })?
                .max(8);
            let bytes = capacity.checked_mul(size_of::<T::Value>()).ok_or_else(|| {
                self.owner
                    .record(anyhow::anyhow!("Configuration sequence layout overflow"))
            })?;
            // Serde Vec starts empty because size_hint is None. Its doubling
            // growth and minimum capacity <=8 are covered before each element.
            // Keep old capacity charges as conservative overlap/retention.
            self.owner.serde_admit(bytes)?;
            self.capacity = capacity;
        }
        let result = self.inner.next_element_seed(Seed {
            inner: seed,
            owner: self.owner,
        })?;
        if result.is_some() {
            self.length += 1;
        }
        Ok(result)
    }
    fn size_hint(&self) -> Option<usize> {
        None
    }
}

struct Map<'a, A> {
    inner: A,
    owner: &'a mut Owner,
    key_bytes: usize,
    growing: bool,
}
impl<'de, A: MapAccess<'de>> MapAccess<'de> for Map<'_, A> {
    type Error = A::Error;
    fn next_key_seed<T: DeserializeSeed<'de>>(
        &mut self,
        seed: T,
    ) -> std::result::Result<Option<T::Value>, A::Error> {
        self.key_bytes = size_of::<T::Value>();
        self.inner.next_key_seed(Seed {
            inner: seed,
            owner: self.owner,
        })
    }
    fn next_value_seed<T: DeserializeSeed<'de>>(
        &mut self,
        seed: T,
    ) -> std::result::Result<T::Value, A::Error> {
        if self.growing {
            // Two BTree nodes per new entry cover leaf/branch split overlap:
            // 11 key/value cells, 13 pointers, u16 lengths and alignment.
            let bytes = self
                .key_bytes
                .checked_add(size_of::<T::Value>())
                .and_then(|n| n.checked_mul(11))
                .and_then(|n| n.checked_add(16 * size_of::<usize>()))
                .and_then(|n| n.checked_mul(2))
                .ok_or_else(|| {
                    self.owner
                        .record(anyhow::anyhow!("Configuration map layout overflow"))
                })?;
            self.owner.serde_admit(bytes)?;
        }
        self.inner.next_value_seed(Seed {
            inner: seed,
            owner: self.owner,
        })
    }
    fn size_hint(&self) -> Option<usize> {
        None
    }
}

struct Enum<'a, A> {
    inner: A,
    owner: &'a mut Owner,
}
impl<'a, 'de, A: de::EnumAccess<'de>> de::EnumAccess<'de> for Enum<'a, A> {
    type Error = A::Error;
    type Variant = Variant<'a, A::Variant>;
    fn variant_seed<T: DeserializeSeed<'de>>(
        self,
        seed: T,
    ) -> std::result::Result<(T::Value, Self::Variant), A::Error> {
        let (value, inner) = self.inner.variant_seed(Seed {
            inner: seed,
            owner: self.owner,
        })?;
        Ok((
            value,
            Variant {
                inner,
                owner: self.owner,
            },
        ))
    }
}
struct Variant<'a, A> {
    inner: A,
    owner: &'a mut Owner,
}
impl<'de, A: de::VariantAccess<'de>> de::VariantAccess<'de> for Variant<'_, A> {
    type Error = A::Error;
    fn unit_variant(self) -> std::result::Result<(), A::Error> {
        self.inner.unit_variant()
    }
    fn newtype_variant_seed<T: DeserializeSeed<'de>>(
        self,
        seed: T,
    ) -> std::result::Result<T::Value, A::Error> {
        self.inner.newtype_variant_seed(Seed {
            inner: seed,
            owner: self.owner,
        })
    }
    fn tuple_variant<V: Visitor<'de>>(
        self,
        len: usize,
        visitor: V,
    ) -> std::result::Result<V::Value, A::Error> {
        self.inner.tuple_variant(
            len,
            TrackedVisitor {
                inner: visitor,
                owner: self.owner,
                storage: Storage::Inline,
            },
        )
    }
    fn struct_variant<V: Visitor<'de>>(
        self,
        fields: &'static [&'static str],
        visitor: V,
    ) -> std::result::Result<V::Value, A::Error> {
        self.inner.struct_variant(
            fields,
            TrackedVisitor {
                inner: visitor,
                owner: self.owner,
                storage: Storage::Inline,
            },
        )
    }
}

//! Project-owned sequence records avoid hidden FASTA-library copies.
use crate::kmer::budget::text::BudgetString;

/// All header and sequence backing buffers carry their own run-level charge.
#[derive(Debug, PartialEq)]
pub(crate) struct SequenceRecord {
    id: BudgetString,
    description: Option<BudgetString>,
    sequence: BudgetString,
}

impl SequenceRecord {
    pub(super) fn from_parts(
        id: BudgetString,
        description: Option<BudgetString>,
        sequence: BudgetString,
    ) -> Self {
        Self {
            id,
            description,
            sequence,
        }
    }
    pub(crate) fn id(&self) -> &str {
        self.id.as_str()
    }
    pub(crate) fn desc(&self) -> Option<&str> {
        self.description.as_deref()
    }
    pub(crate) fn seq(&self) -> &[u8] {
        self.sequence.as_bytes()
    }
    /// Move a new header into place; the sequence allocation is retained.
    pub(super) fn replace_header(&mut self, id: BudgetString, description: Option<BudgetString>) {
        self.id = id;
        self.description = description;
    }
}

#[cfg(test)]
impl From<bio::io::fasta::Record> for SequenceRecord {
    fn from(record: bio::io::fasta::Record) -> Self {
        use crate::kmer::budget::MemoryClass;
        let text = |value: &str| {
            BudgetString::from_str(value, None, MemoryClass::RetainedResults).unwrap()
        };
        Self::from_parts(
            text(record.id()),
            record.desc().map(text),
            text(std::str::from_utf8(record.seq()).unwrap()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::budget::{MEMORY_CLASSES, MemoryClass, ResourceBudget, ResourceLimits};
    #[test]
    fn header_replacement_keeps_sequence_allocation_and_borrowed_fasta_bytes() {
        let mut classes = [0; MEMORY_CLASSES];
        classes[MemoryClass::RetainedResults as usize] = 8192;
        let budget = ResourceBudget::new(ResourceLimits {
            memory_bytes: 8192,
            process_overhead_bytes: 0,
            class_bytes: classes,
            temporary_disk_bytes: 0,
            temporary_open_files: 0,
        })
        .unwrap();
        let text = |s: &str| {
            BudgetString::from_str(s, Some(&budget), MemoryClass::RetainedResults).unwrap()
        };
        let sequence = "ACGT".repeat(23);
        let mut record = SequenceRecord::from_parts(
            text("old"),
            Some(text("sequence_index=7")),
            text(&sequence),
        );
        let pointer = record.seq().as_ptr();
        record.replace_header(text("new"), Some(text("sequence_index=0")));
        assert_eq!(record.seq().as_ptr(), pointer);
        let mut bytes = Vec::new();
        crate::io::write_fasta_fields(&mut bytes, record.id(), record.desc(), record.seq())
            .unwrap();
        assert_eq!(
            bytes,
            format!(
                ">new sequence_index=0\n{}\n{}\n",
                &sequence[..80],
                &sequence[80..]
            )
            .into_bytes()
        );
        assert!(budget.snapshot().reserved_by_class[MemoryClass::RetainedResults as usize] > 0);
        drop(record);
        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
    }
}

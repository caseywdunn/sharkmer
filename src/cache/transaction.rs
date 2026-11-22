//! Same-filesystem cache publication with receipted, repeatable rollback.
//! The cooperating cache lock is held throughout every operation in this module.

use super::*;
use crate::kmer::budget::MemoryClass;
use crate::kmer::budget::collections::BudgetVec;
use crate::kmer::budget::path::BudgetPath;
use crate::kmer::budget::text::BudgetString;
use crate::output::{RecoveryCharge, StorageInventory, sync_directory};
use anyhow::ensure;

mod record;
use record::{CONTEXT_BYTES, JOURNAL_BYTES, Record, hash};

const ROOT: &str = ".sharkmer-cache-transaction-v1";
const OWNER: &str = ".sharkmer-cache-owner-v1";
const CONTEXT: &str = "context";
const JOURNAL: &str = "journal";
const ACCEPTED: &str = "accepted";
const MAX_FILES: usize = 16;
#[cfg(test)]
pub(super) const CONTEXT_DISK_BYTES: u64 = CONTEXT_BYTES as u64;
#[cfg(test)]
pub(super) const CONTROL_BYTES: u64 = (CONTEXT_BYTES + JOURNAL_BYTES) as u64;
#[cfg(test)]
pub(super) const CONTROL_PEAK_BYTES: u64 = if cfg!(unix) {
    CONTROL_BYTES
} else {
    2 * CONTROL_BYTES
};

struct Locations {
    cache: BudgetPath,
    owner: BudgetPath,
    root: BudgetPath,
    context: BudgetPath,
    journal: BudgetPath,
    accepted: BudgetPath,
    old: [BudgetPath; 2],
    link: [BudgetPath; 2],
    restore: [BudgetPath; 2],
}

impl Locations {
    fn new(cache: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<Self> {
        let cache = BudgetPath::from_path(cache, budget, MemoryClass::RetainedResults)?;
        let root = cache.try_join(ROOT)?;
        Ok(Self {
            owner: cache.try_join(OWNER)?,
            context: root.try_join(CONTEXT)?,
            journal: root.try_join(JOURNAL)?,
            accepted: root.try_join(ACCEPTED)?,
            old: [root.try_join("old.data")?, root.try_join("old.meta")?],
            link: [
                root.try_join("publish.data")?,
                root.try_join("publish.meta")?,
            ],
            restore: [
                root.try_join("restore.data")?,
                root.try_join("restore.meta")?,
            ],
            cache,
            root,
        })
    }

    fn public(&self, key: &[u8]) -> Result<[BudgetPath; 2]> {
        ensure!(record::valid_key(key), "Invalid cache recovery key");
        let key = std::str::from_utf8(key)?;
        let name = |suffix| {
            BudgetString::format(
                format_args!("{key}.{suffix}"),
                self.cache.budget(),
                MemoryClass::RetainedResults,
            )
        };
        Ok([
            self.cache.try_join(name("fastq.gz")?.as_str())?,
            self.cache.try_join(name("meta.yaml")?.as_str())?,
        ])
    }
}

/// Registered private storage persists after failed cleanup or interruption.
pub(super) struct Transaction {
    // All files are closed before public links exist. Their charges survive
    // failed recovery; only confirmed cleanup releases them.
    data: Option<CacheTemporary>,
    metadata: Option<CacheTemporary>,
    context_file: Option<CacheTemporary>,
    journal_file: Option<CacheTemporary>,
    backup_charge: RecoveryCharge,
    locations: Locations,
    public: [BudgetPath; 2],
    record: Record,
    registered: bool,
    finished: bool,
}

impl Transaction {
    /// Register an exact private child before creating any download payload.
    pub(super) fn begin(
        cache: &Path,
        key: &str,
        previous: bool,
        budget: Option<&Arc<ResourceBudget>>,
    ) -> Result<Self> {
        recover(cache, budget)?;
        let locations = Locations::new(cache, budget)?;
        let public = locations.public(key.as_bytes())?;
        let record = Record::initial(key, &public, previous)?;
        let mut transaction = Self {
            data: None,
            metadata: None,
            context_file: None,
            journal_file: None,
            backup_charge: RecoveryCharge::default(),
            locations,
            public,
            record,
            registered: false,
            finished: false,
        };
        let owner = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&transaction.locations.owner)?;
        transaction.registered = true;
        owner.sync_all()?;
        drop(owner);
        sync_directory(&transaction.locations.cache)?;
        fs::create_dir(&transaction.locations.root)?;
        sync_directory(&transaction.locations.cache)?;
        let bytes = transaction.record.encode::<CONTEXT_BYTES>()?;
        transaction.context_file = Some(write_control(
            &transaction.locations.root,
            ".context-",
            &bytes,
        )?);
        let context = transaction
            .context_file
            .as_mut()
            .context("Cache context is unavailable")?;
        context.begin_link_publication()?;
        fs::hard_link(context.path(), &transaction.locations.context)?;
        sync_directory(&transaction.locations.root)?;
        Ok(transaction)
    }

    /// Borrow the registered directory on the cache filesystem.
    pub(super) fn directory(&self) -> &Path {
        &self.locations.root
    }

    /// Bind both new files and old backups before the first public mutation.
    pub(super) fn publish(mut self, data: CacheTemporary, metadata: CacheTemporary) -> Result<()> {
        self.data = Some(data);
        self.metadata = Some(metadata);
        self.prepare()?;
        for index in 0..2 {
            let temporary = if index == 0 {
                self.data.as_mut()
            } else {
                self.metadata.as_mut()
            }
            .context("Cache payload is unavailable")?;
            temporary.begin_link_publication()?;
            fs::hard_link(temporary.path(), &self.locations.link[index])?;
            fs::rename(&self.locations.link[index], &self.public[index])?;
        }
        sync_directory(&self.locations.cache)?;
        let accepted = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&self.locations.accepted)?;
        accepted.sync_all()?;
        drop(accepted);
        sync_directory(&self.locations.root)?;
        self.finish_recovery()
    }

    fn prepare(&mut self) -> Result<()> {
        ensure!(
            public_matches(&self.public, self.record.previous.as_ref())?,
            "Cache changed before publication; cache preserved"
        );
        let budget = self.locations.cache.budget();
        if self.record.previous.is_some() {
            let bytes = fs::metadata(&self.public[0])?
                .len()
                .checked_add(fs::metadata(&self.public[1])?.len())
                .context("Cache backup extent overflow")?;
            self.backup_charge = RecoveryCharge::reserve(bytes, budget)?;
            for index in 0..2 {
                fs::hard_link(&self.public[index], &self.locations.old[index])?;
            }
            sync_directory(&self.locations.root)?;
        }
        let next = [
            hash(self.data.as_ref().context("Missing cache data")?.path())?,
            hash(
                self.metadata
                    .as_ref()
                    .context("Missing cache metadata")?
                    .path(),
            )?,
        ];
        self.record.next = Some(next);
        let bytes = self.record.encode::<JOURNAL_BYTES>()?;
        self.journal_file = Some(write_control(&self.locations.root, ".journal-", &bytes)?);
        let journal = self
            .journal_file
            .as_mut()
            .context("Cache journal is unavailable")?;
        journal.begin_link_publication()?;
        fs::hard_link(journal.path(), &self.locations.journal)?;
        sync_directory(&self.locations.root)?;
        Ok(())
    }

    fn credit(&self) -> Result<u64> {
        [
            &self.data,
            &self.metadata,
            &self.context_file,
            &self.journal_file,
        ]
        .into_iter()
        .try_fold(self.backup_charge.bytes(), |sum, item| {
            sum.checked_add(item.as_ref().map_or(0, CacheTemporary::reserved_bytes))
                .context("Cache live credit overflow")
        })
    }

    fn finish_recovery(&mut self) -> Result<()> {
        recover_locations(&self.locations, self.credit()?)?;
        for temporary in [
            &mut self.data,
            &mut self.metadata,
            &mut self.context_file,
            &mut self.journal_file,
        ]
        .into_iter()
        .flatten()
        {
            temporary.confirm_removed_or_accepted();
        }
        self.backup_charge.release();
        self.finished = true;
        Ok(())
    }
}

impl Drop for Transaction {
    fn drop(&mut self) {
        if self.registered
            && !self.finished
            && let Err(error) = self.finish_recovery()
        {
            log::error!(
                "Cache recovery remains pending in {}: {error:#}",
                self.locations.root.display()
            );
        }
    }
}

fn write_control(root: &BudgetPath, prefix: &str, bytes: &[u8]) -> Result<CacheTemporary> {
    let mut file = CacheTemporary::create(root, prefix, root.budget())?;
    let written = file.writer().write_all(bytes).map_err(anyhow::Error::from);
    file.finish(written)?;
    Ok(file)
}

fn regular(path: &Path) -> Result<Option<fs::Metadata>> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            ensure!(
                metadata.file_type().is_file(),
                "Non-regular cache recovery path {}; cache preserved",
                path.display()
            );
            Ok(Some(metadata))
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(error).context("Unable to inspect cache recovery path"),
    }
}

fn zero_marker(path: &Path) -> Result<bool> {
    let Some(metadata) = regular(path)? else {
        return Ok(false);
    };
    ensure!(
        metadata.len() == 0,
        "Modified cache ownership/acceptance marker; cache preserved"
    );
    Ok(true)
}

fn public_matches(public: &[BudgetPath; 2], expected: Option<&[[u8; 32]; 2]>) -> Result<bool> {
    for index in 0..2 {
        match (regular(&public[index])?, expected) {
            (None, None) => {}
            (Some(_), Some(expected)) if hash(&public[index])? == expected[index] => {}
            _ => return Ok(false),
        }
    }
    Ok(true)
}

/// Recover only the exact registered child under the cooperating cache lock.
pub(super) fn recover(cache: &Path, budget: Option<&Arc<ResourceBudget>>) -> Result<()> {
    recover_locations(&Locations::new(cache, budget)?, 0)
}

fn recover_locations(locations: &Locations, credit: u64) -> Result<()> {
    let owned = zero_marker(&locations.owner)?;
    match fs::symlink_metadata(&locations.root) {
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            if owned {
                fs::remove_file(&locations.owner)?;
                sync_directory(&locations.cache)?;
            }
            return Ok(());
        }
        Ok(metadata) => ensure!(
            owned && metadata.file_type().is_dir(),
            "Unregistered/non-directory cache transaction; cache preserved"
        ),
        Err(error) => return Err(error).context("Unable to inspect cache transaction"),
    }
    let budget = locations.cache.budget();
    let mut files =
        BudgetVec::with_capacity_in_class(MAX_FILES, budget, MemoryClass::RetainedResults)?;
    let mut inventory = StorageInventory::private(budget);
    for entry in fs::read_dir(&locations.root)? {
        let name = entry?.file_name();
        let text = name
            .to_str()
            .context("Unrecognized cache transaction path; cache preserved")?;
        ensure!(
            files.len() < MAX_FILES && allowed_name(text),
            "Unexpected cache transaction contents; cache preserved"
        );
        let path = locations.root.try_join(&name)?;
        inventory.include(&path)?;
        files.push(path)?;
    }
    // Before any unlink/rename, admit every remaining private inode. Live
    // owners supply only credit belonging to this transaction and ledger.
    let context = Record::read::<CONTEXT_BYTES>(&locations.context)?;
    let journal = Record::read::<JOURNAL_BYTES>(&locations.journal)?;
    let accepted = zero_marker(&locations.accepted)?;
    if let Some(journal) = &journal {
        validate_private_files(&files, journal)?;
    }
    #[cfg(unix)]
    let restoration_bytes = 0;
    #[cfg(not(unix))]
    let restoration_bytes = if !accepted && let Some(previous) = journal.and_then(|j| j.previous) {
        let public = locations.public(&journal.context("Missing cache journal")?.key)?;
        if public_matches(&public, Some(&previous))? {
            0
        } else {
            fs::metadata(&locations.old[0])?
                .len()
                .checked_add(fs::metadata(&locations.old[1])?.len())
                .context("Cache restoration size overflow")?
        }
    } else {
        0
    };
    let mut charge = inventory.admit(restoration_bytes, credit, 0, budget)?;
    match (context, journal) {
        (Some(context), Some(journal)) => {
            ensure!(
                context.key == journal.key && context.previous == journal.previous,
                "Conflicting cache receipts; cache preserved"
            );
            let public = locations.public(&context.key)?;
            let next = journal.next.context("Cache journal lacks the new pair")?;
            if accepted {
                ensure!(
                    public_matches(&public, Some(&next))?,
                    "Accepted cache changed; recovery evidence preserved"
                );
            } else {
                rollback(locations, &public, &journal)?;
            }
        }
        (Some(context), None) => {
            let public = locations.public(&context.key)?;
            let payload = files
                .iter()
                .any(|p| !control_name(p.file_name().and_then(|n| n.to_str()).unwrap_or("")));
            if payload {
                ensure!(
                    !accepted && public_matches(&public, context.previous.as_ref())?,
                    "Missing cache recovery journal after publication; evidence preserved"
                );
            } else if accepted {
                // Journal removal follows all payload removal. Validate the
                // surviving pair before finishing this control-only cleanup.
                let key = std::str::from_utf8(&context.key)?;
                validate_entry_with_budget(key, &public[0], &public[1], budget)?
                    .context("Accepted cache pair disappeared")?;
            } else {
                ensure!(
                    public_matches(&public, context.previous.as_ref())?,
                    "Cache baseline changed without a journal; evidence preserved"
                );
            }
        }
        (None, None) => ensure!(
            files.iter().all(|p| p
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with(".context-") || n == ACCEPTED)),
            "Missing cache preparation context; evidence preserved"
        ),
        (None, Some(_)) => anyhow::bail!("Missing cache preparation context; journal preserved"),
    }
    cleanup(locations, &files)?;
    charge.release();
    Ok(())
}

fn rollback(locations: &Locations, public: &[BudgetPath; 2], record: &Record) -> Result<()> {
    let next = record
        .next
        .context("Cache recovery journal lacks new hashes")?;
    for index in 0..2 {
        if regular(&public[index])?.is_some() {
            let actual = hash(&public[index])?;
            ensure!(
                actual == next[index]
                    || record
                        .previous
                        .as_ref()
                        .is_some_and(|old| actual == old[index]),
                "Modified public cache file; recovery evidence preserved"
            );
        }
    }
    if let Some(previous) = &record.previous {
        if !public_matches(public, Some(previous))? {
            // Neither backup is consumed. A failure at any restoration step
            // leaves both original receipts available to the next attempt.
            for (index, expected) in previous.iter().enumerate() {
                ensure!(
                    hash(&locations.old[index])? == *expected,
                    "Modified cache backup; cache preserved"
                );
                if regular(&locations.restore[index])?.is_some() {
                    ensure!(
                        hash(&locations.restore[index])? == *expected,
                        "Modified cache restoration link; cache preserved"
                    );
                }
            }
            for (index, destination) in public.iter().enumerate() {
                if regular(destination)?.is_some() && hash(destination)? == previous[index] {
                    continue;
                }
                if regular(&locations.restore[index])?.is_none() {
                    fs::hard_link(&locations.old[index], &locations.restore[index])?;
                }
                fs::rename(&locations.restore[index], destination)?;
            }
        }
    } else {
        for path in public {
            if regular(path)?.is_some() {
                fs::remove_file(path)?;
            }
        }
    }
    sync_directory(&locations.cache)
}

fn validate_private_files(files: &[BudgetPath], journal: &Record) -> Result<()> {
    let next = journal.next.context("Cache journal lacks payload hashes")?;
    for path in files {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .context("Unrecognized cache payload")?;
        let expected = if name == "old.data" || name == "restore.data" {
            Some(journal.previous.context("Unexpected previous cache data")?[0])
        } else if name == "old.meta" || name == "restore.meta" {
            Some(
                journal
                    .previous
                    .context("Unexpected previous cache metadata")?[1],
            )
        } else if name == "publish.data" || name.starts_with(".new-data-") {
            Some(next[0])
        } else if name == "publish.meta" || name.starts_with(".sharkmer-cache-meta-") {
            Some(next[1])
        } else {
            None
        };
        if let Some(expected) = expected {
            ensure!(
                hash(path)? == expected,
                "Modified private cache payload; recovery evidence preserved"
            );
        }
    }
    Ok(())
}

fn control_name(name: &str) -> bool {
    matches!(name, CONTEXT | JOURNAL | ACCEPTED)
        || name.starts_with(".context-")
        || name.starts_with(".journal-")
}

fn allowed_name(name: &str) -> bool {
    control_name(name)
        || matches!(
            name,
            "old.data"
                | "old.meta"
                | "publish.data"
                | "publish.meta"
                | "restore.data"
                | "restore.meta"
        )
        || name.starts_with(".new-data-")
        || name.starts_with(".sharkmer-cache-meta-")
}

fn cleanup(locations: &Locations, files: &[BudgetPath]) -> Result<()> {
    for path in files {
        if path.as_ref() != locations.context.as_ref()
            && path.as_ref() != locations.journal.as_ref()
            && path.as_ref() != locations.accepted.as_ref()
        {
            match fs::remove_file(path) {
                Ok(()) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => return Err(error).context("Unable to remove owned cache storage"),
            }
        }
    }
    sync_directory(&locations.root)?;
    for path in [&locations.journal, &locations.context, &locations.accepted] {
        if regular(path)?.is_some() {
            fs::remove_file(path)?;
            sync_directory(&locations.root)?;
        }
    }
    fs::remove_dir(&locations.root)?;
    sync_directory(&locations.cache)?;
    fs::remove_file(&locations.owner)?;
    sync_directory(&locations.cache)
}

#[cfg(test)]
mod tests;

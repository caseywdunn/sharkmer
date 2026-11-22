//! Verify this transaction's immutable receipts without parsing its own YAML.

use super::*;

/// Clean an accepted generation using hashes captured before publication.
/// Historical restart recovery still validates and parses its on-disk journal.
pub(in crate::output) fn accept(
    transaction: &OutputTransaction,
    manifest_digest: &[u8; 32],
) -> Result<()> {
    ensure!(
        transaction.committed && transaction.recovery_pending,
        "Output transaction has not accepted a replacement generation"
    );
    let budget = transaction.budget.as_ref();
    let journal = journal_path(&transaction.output_directory, &transaction.sample)?;
    require_regular(&journal)?;
    ensure!(
        Some(sha256_bytes(&journal)?) == transaction.recovery_journal_sha256,
        "Modified output recovery receipt; backups preserved"
    );
    require_regular(&transaction.manifest_path)?;
    ensure!(
        sha256_bytes(&transaction.manifest_path)? == *manifest_digest,
        "Modified accepted output manifest; backups preserved"
    );
    require_directory(&transaction.staging_directory)?;
    let backup = transaction.staging_directory.try_join(BACKUP_DIRECTORY)?;
    require_directory(&backup)?;
    let manifest = backup.try_join(BACKUP_MANIFEST)?;
    require_regular(&manifest)?;
    ensure!(
        transaction.previous_manifest_sha256.as_deref()
            == Some(sha256_file_owned(&manifest, budget)?.as_str()),
        "Modified output recovery manifest; backups preserved"
    );
    let previous = transaction
        .previous_complete
        .as_ref()
        .context("Missing previous generation during accepted cleanup")?;
    let mut inventory = backup_inventory();
    for receipt in &previous.files {
        let path = backup.try_join(receipt.path.as_str())?;
        match fs::symlink_metadata(&path) {
            Ok(_) => {
                verify_file(&backup, receipt)?;
                inventory_entry(&mut inventory, receipt)?;
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error).context("Failed to inspect accepted recovery backup"),
        }
    }
    ensure!(
        Some(<[u8; 32]>::from(inventory.finalize())) == transaction.backup_inventory_sha256,
        "Modified output recovery inventory; backups preserved"
    );
    for staged in &transaction.staged_outputs {
        let path = transaction
            .output_directory
            .try_join(staged.receipt.path.as_str())?;
        require_regular(&path)?;
        ensure!(
            sha256_file_owned(&path, budget)? == staged.receipt.sha256,
            "Modified accepted output {}; backups preserved",
            path.display()
        );
    }
    let old = receipt_map(&transaction.sample, &previous.files, budget)?;
    let new = receipt_index(
        &transaction.sample,
        transaction
            .staged_outputs
            .iter()
            .map(|staged| (staged.receipt.path.as_str(), staged.receipt.sha256.as_str())),
        transaction.staged_outputs.len(),
        budget,
    )?;
    let cleanup = storage::Cleanup::read(
        &transaction.staging_directory,
        &backup,
        &old,
        &new,
        &transaction.run_id,
        manifest_digest,
    )?;
    sync_directory(&transaction.output_directory)?;
    cleanup.discard(
        &transaction.output_directory,
        &journal,
        &transaction.staging_directory,
        &backup,
    )
}

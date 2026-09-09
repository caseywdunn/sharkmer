use anyhow::{Context, Result, bail, ensure};
use fs2::FileExt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Component, Path, PathBuf};

const MANIFEST_SCHEMA_VERSION: u32 = 1;
const MANIFEST_PRODUCER: &str = "sharkmer";

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct OutputReceipt {
    pub(crate) path: String,
    pub(crate) sha256: String,
}

#[derive(Debug, Deserialize, Serialize)]
struct OutputManifest {
    schema_version: u32,
    producer: String,
    sample: String,
    run_id: String,
    status: String,
    files: Vec<OutputReceipt>,
    #[serde(skip_serializing_if = "Option::is_none")]
    failure_reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    staging_directory: Option<String>,
}

struct StagedOutput {
    receipt: OutputReceipt,
    staging_path: PathBuf,
}

pub(crate) struct OutputTransaction {
    output_directory: PathBuf,
    sample: String,
    run_id: String,
    manifest_basename: String,
    manifest_path: PathBuf,
    staging_basename: String,
    staging_directory: PathBuf,
    staged_outputs: Vec<StagedOutput>,
    failure_receipts: Vec<OutputReceipt>,
    _lock_file: File,
    committed: bool,
}

impl OutputTransaction {
    pub(crate) fn begin(output_directory: &Path, sample: &str) -> Result<Self> {
        fs::create_dir_all(output_directory).with_context(|| {
            format!(
                "Failed to create output directory {}",
                output_directory.display()
            )
        })?;
        let lock_path = output_directory.join(format!("{sample}.lock"));
        reject_non_regular_existing_path(&lock_path, "output lock")?;
        let lock_file = OpenOptions::new()
            .create(true)
            .truncate(false)
            .read(true)
            .write(true)
            .open(&lock_path)
            .with_context(|| format!("Failed to open output lock {}", lock_path.display()))?;
        lock_file
            .lock_exclusive()
            .with_context(|| format!("Failed to lock output sample {sample}"))?;

        let manifest_basename = format!("{sample}.manifest.yaml");
        let manifest_path = output_directory.join(&manifest_basename);
        let previous_manifest = read_previous_manifest(&manifest_path, sample)?;
        if previous_manifest.is_none() {
            reject_unowned_stats_output(output_directory, sample)?;
        }
        let previous_receipts = previous_manifest
            .as_ref()
            .map(|manifest| manifest.files.clone())
            .unwrap_or_default();
        let cleanup = previous_manifest
            .as_ref()
            .map(|manifest| verify_previous_outputs(output_directory, manifest))
            .transpose()?;

        let run_id = new_run_id();
        let staging_basename = format!(".sharkmer-stage-{sample}-{run_id}");
        let staging_directory = output_directory.join(&staging_basename);
        fs::create_dir(&staging_directory).with_context(|| {
            format!(
                "Failed to create staging directory {}",
                staging_directory.display()
            )
        })?;

        let mut transaction = Self {
            output_directory: output_directory.to_path_buf(),
            sample: sample.to_string(),
            run_id,
            manifest_basename,
            manifest_path,
            staging_basename,
            staging_directory,
            staged_outputs: Vec::new(),
            failure_receipts: previous_receipts,
            _lock_file: lock_file,
            committed: false,
        };
        transaction.write_manifest("in_progress", transaction.failure_receipts.clone(), None)?;
        if let Some(cleanup) = cleanup {
            cleanup.remove()?;
        }
        transaction.failure_receipts.clear();
        transaction.write_manifest("in_progress", Vec::new(), None)?;
        Ok(transaction)
    }

    pub(crate) fn run_id(&self) -> &str {
        &self.run_id
    }

    pub(crate) fn manifest_basename(&self) -> &str {
        &self.manifest_basename
    }

    pub(crate) fn stage_file(
        &mut self,
        basename: &str,
        write_file: impl FnOnce(&mut File) -> Result<()>,
    ) -> Result<()> {
        ensure_owned_output_basename(&self.sample, basename)?;
        ensure!(
            !self
                .staged_outputs
                .iter()
                .any(|staged| staged.receipt.path == basename),
            "Output file {} was staged more than once",
            basename
        );

        let staging_path = self.staging_directory.join(basename);
        let mut file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&staging_path)
            .with_context(|| format!("Failed to create staged output {basename}"))?;
        write_file(&mut file)
            .with_context(|| format!("Failed to write staged output {basename}"))?;
        file.flush()
            .with_context(|| format!("Failed to flush staged output {basename}"))?;
        file.sync_all()
            .with_context(|| format!("Failed to sync staged output {basename}"))?;
        drop(file);
        let sha256 = sha256_file(&staging_path)?;
        self.staged_outputs.push(StagedOutput {
            receipt: OutputReceipt {
                path: basename.to_string(),
                sha256,
            },
            staging_path,
        });
        Ok(())
    }

    pub(crate) fn commit(mut self) -> Result<()> {
        self.staged_outputs
            .sort_by(|left, right| left.receipt.path.cmp(&right.receipt.path));
        for staged in &self.staged_outputs {
            reject_existing_destination(
                &self.output_directory.join(&staged.receipt.path),
                &staged.receipt.path,
            )?;
        }

        let receipts = self.receipts();
        self.failure_receipts = receipts.clone();
        self.write_manifest("publishing", receipts.clone(), None)?;
        for staged in &self.staged_outputs {
            publish_no_clobber(
                &staged.staging_path,
                &self.output_directory.join(&staged.receipt.path),
            )?;
        }
        fs::remove_dir(&self.staging_directory).with_context(|| {
            format!(
                "Failed to remove staging directory {}",
                self.staging_directory.display()
            )
        })?;
        sync_directory(&self.output_directory)?;
        self.write_manifest("complete", receipts, None)?;
        self.committed = true;
        Ok(())
    }

    fn receipts(&self) -> Vec<OutputReceipt> {
        self.staged_outputs
            .iter()
            .map(|staged| staged.receipt.clone())
            .collect()
    }

    fn write_manifest(
        &mut self,
        status: &str,
        files: Vec<OutputReceipt>,
        failure_reason: Option<String>,
    ) -> Result<()> {
        let manifest = OutputManifest {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER.to_string(),
            sample: self.sample.clone(),
            run_id: self.run_id.clone(),
            status: status.to_string(),
            files,
            failure_reason,
            staging_directory: (status != "complete").then(|| self.staging_basename.clone()),
        };
        write_manifest_atomic(&self.output_directory, &self.manifest_path, &manifest)
    }
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        if self.committed {
            return;
        }
        let _ = self.write_manifest(
            "failed",
            self.failure_receipts.clone(),
            Some("run ended before output publication completed".to_string()),
        );
        let _ = fs::remove_dir_all(&self.staging_directory);
    }
}

struct VerifiedCleanup {
    files: Vec<PathBuf>,
    staging_directory: Option<PathBuf>,
}

impl VerifiedCleanup {
    fn remove(self) -> Result<()> {
        for path in self.files {
            fs::remove_file(&path).with_context(|| {
                format!("Failed to remove prior owned output {}", path.display())
            })?;
        }
        if let Some(path) = self.staging_directory {
            fs::remove_dir_all(&path).with_context(|| {
                format!(
                    "Failed to remove prior staging directory {}",
                    path.display()
                )
            })?;
        }
        Ok(())
    }
}

fn read_previous_manifest(path: &Path, sample: &str) -> Result<Option<OutputManifest>> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => {
            return Err(error)
                .with_context(|| format!("Failed to inspect output manifest {}", path.display()));
        }
    };
    ensure!(
        metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
        "Refusing non-regular output manifest {}",
        path.display()
    );
    let file = File::open(path)
        .with_context(|| format!("Failed to open output manifest {}", path.display()))?;
    let manifest: OutputManifest = serde_yaml_ng::from_reader(file)
        .with_context(|| format!("Failed to parse output manifest {}", path.display()))?;
    ensure!(
        manifest.schema_version == MANIFEST_SCHEMA_VERSION
            && manifest.producer == MANIFEST_PRODUCER
            && manifest.sample == sample,
        "Refusing foreign or unsupported output manifest {}",
        path.display()
    );
    ensure!(
        matches!(
            manifest.status.as_str(),
            "in_progress" | "publishing" | "complete" | "failed"
        ),
        "Output manifest {} has unknown status {}",
        path.display(),
        manifest.status
    );
    ensure!(
        !manifest.run_id.is_empty()
            && manifest
                .run_id
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() || byte == b'-'),
        "Output manifest {} has an invalid run ID",
        path.display()
    );
    if manifest.status == "complete" {
        ensure!(
            manifest.staging_directory.is_none() && manifest.failure_reason.is_none(),
            "Complete output manifest {} contains incomplete-run metadata",
            path.display()
        );
    } else {
        ensure!(
            manifest.staging_directory.is_some(),
            "Noncomplete output manifest {} is missing its staging directory",
            path.display()
        );
    }
    Ok(Some(manifest))
}

fn verify_previous_outputs(
    output_directory: &Path,
    manifest: &OutputManifest,
) -> Result<VerifiedCleanup> {
    let mut files = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let stats_basename = format!("{}.stats.yaml", manifest.sample);
    if matches!(manifest.status.as_str(), "publishing" | "complete") {
        ensure!(
            manifest
                .files
                .iter()
                .filter(|receipt| receipt.path == stats_basename)
                .count()
                == 1,
            "{} output manifest must contain exactly one stats receipt",
            manifest.status
        );
    }
    for receipt in &manifest.files {
        ensure_owned_output_basename(&manifest.sample, &receipt.path)?;
        ensure!(
            receipt.sha256.len() == 64
                && receipt
                    .sha256
                    .bytes()
                    .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()),
            "Output manifest receipt for {} has an invalid SHA-256",
            receipt.path
        );
        ensure!(
            seen.insert(receipt.path.as_str()),
            "Output manifest contains duplicate receipt {}",
            receipt.path
        );
        let path = output_directory.join(&receipt.path);
        let metadata = match fs::symlink_metadata(&path) {
            Ok(metadata) => metadata,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(error) => {
                return Err(error)
                    .with_context(|| format!("Failed to inspect prior output {}", path.display()));
            }
        };
        ensure!(
            metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
            "Refusing to replace non-regular prior output {}",
            path.display()
        );
        let actual_sha256 = sha256_file(&path)?;
        ensure!(
            actual_sha256 == receipt.sha256,
            "Refusing to replace modified prior output {}",
            path.display()
        );
        files.push(path);
    }

    let staging_directory = match &manifest.staging_directory {
        Some(basename) => {
            ensure_safe_basename(basename)?;
            ensure!(
                basename == &format!(".sharkmer-stage-{}-{}", manifest.sample, manifest.run_id),
                "Output manifest has an invalid staging directory {}",
                basename
            );
            let path = output_directory.join(basename);
            match fs::symlink_metadata(&path) {
                Ok(metadata) => {
                    ensure!(
                        metadata.file_type().is_dir() && !metadata.file_type().is_symlink(),
                        "Refusing non-directory prior staging path {}",
                        path.display()
                    );
                    Some(path)
                }
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!("Failed to inspect prior staging path {}", path.display())
                    });
                }
            }
        }
        None => None,
    };
    Ok(VerifiedCleanup {
        files,
        staging_directory,
    })
}

fn reject_non_regular_existing_path(path: &Path, label: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) => {
            ensure!(
                metadata.file_type().is_file() && !metadata.file_type().is_symlink(),
                "Refusing non-regular {} {}",
                label,
                path.display()
            );
            Ok(())
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect {}", path.display())),
    }
}

fn reject_existing_destination(path: &Path, basename: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(_) => bail!(
            "Refusing to overwrite unowned or concurrently created output {}",
            basename
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect output {basename}")),
    }
}

fn reject_unowned_stats_output(output_directory: &Path, sample: &str) -> Result<()> {
    let basename = format!("{sample}.stats.yaml");
    let path = output_directory.join(&basename);
    match fs::symlink_metadata(&path) {
        Ok(_) => bail!(
            "Refusing existing unowned stats output {}. Use a fresh output directory or move the legacy file.",
            path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("Failed to inspect output {basename}")),
    }
}

fn publish_no_clobber(staging_path: &Path, destination: &Path) -> Result<()> {
    fs::hard_link(staging_path, destination).with_context(|| {
        format!(
            "Failed to publish output {} without overwriting an existing file",
            destination.display()
        )
    })?;
    fs::remove_file(staging_path).with_context(|| {
        format!(
            "Failed to remove staged output after publishing {}",
            destination.display()
        )
    })?;
    Ok(())
}

fn write_manifest_atomic(
    output_directory: &Path,
    manifest_path: &Path,
    manifest: &OutputManifest,
) -> Result<()> {
    let mut temporary = tempfile::Builder::new()
        .prefix(".sharkmer-manifest-")
        .tempfile_in(output_directory)
        .context("Failed to create temporary output manifest")?;
    serde_yaml_ng::to_writer(temporary.as_file_mut(), manifest)
        .context("Failed to write output manifest")?;
    temporary
        .as_file_mut()
        .flush()
        .context("Failed to flush output manifest")?;
    temporary
        .as_file_mut()
        .sync_all()
        .context("Failed to sync output manifest")?;
    temporary
        .persist(manifest_path)
        .map_err(|error| error.error)
        .with_context(|| {
            format!(
                "Failed to atomically publish output manifest {}",
                manifest_path.display()
            )
        })?;
    sync_directory(output_directory)?;
    Ok(())
}

fn ensure_owned_output_basename(sample: &str, basename: &str) -> Result<()> {
    ensure_safe_basename(basename)?;
    let stats_basename = format!("{sample}.stats.yaml");
    let fasta_prefix = format!("{sample}_");
    let gene_name = basename
        .strip_prefix(&fasta_prefix)
        .and_then(|rest| rest.strip_suffix(".fasta"));
    ensure!(
        basename == stats_basename || gene_name.is_some_and(|gene| !gene.is_empty()),
        "Output receipt path is outside the owned sample namespace: {}",
        basename
    );
    Ok(())
}

fn ensure_safe_basename(basename: &str) -> Result<()> {
    let path = Path::new(basename);
    let mut components = path.components();
    ensure!(
        matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none(),
        "Output manifest path must be a basename: {}",
        basename
    );
    Ok(())
}

fn sha256_file(path: &Path) -> Result<String> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open {} for SHA-256", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let bytes_read = file
            .read(&mut buffer)
            .with_context(|| format!("Failed to read {} for SHA-256", path.display()))?;
        if bytes_read == 0 {
            break;
        }
        hasher.update(&buffer[..bytes_read]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

#[cfg(unix)]
fn sync_directory(path: &Path) -> Result<()> {
    File::open(path)
        .with_context(|| format!("Failed to open directory {} for sync", path.display()))?
        .sync_all()
        .with_context(|| format!("Failed to sync directory {}", path.display()))
}

#[cfg(not(unix))]
fn sync_directory(_path: &Path) -> Result<()> {
    Ok(())
}

fn new_run_id() -> String {
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    format!(
        "{timestamp:x}-{:x}-{:016x}",
        std::process::id(),
        rand::random::<u64>()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stage_text(transaction: &mut OutputTransaction, basename: &str, value: &str) {
        transaction
            .stage_file(basename, |file| {
                file.write_all(value.as_bytes())?;
                Ok(())
            })
            .unwrap();
    }

    fn read_manifest(directory: &Path, sample: &str) -> OutputManifest {
        serde_yaml_ng::from_reader(
            File::open(directory.join(format!("{sample}.manifest.yaml"))).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn changed_output_set_removes_only_receipted_prior_files() {
        let directory = tempfile::tempdir().unwrap();
        let sentinel = directory.path().join("sentinel.txt");
        fs::write(&sentinel, "keep").unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_a.fasta", "a");
        stage_text(&mut first, "sample_b.fasta", "b");
        stage_text(&mut first, "sample.stats.yaml", "first");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "sample").unwrap();
        assert!(!directory.path().join("sample_a.fasta").exists());
        assert!(!directory.path().join("sample_b.fasta").exists());
        stage_text(&mut second, "sample_a.fasta", "new-a");
        stage_text(&mut second, "sample.stats.yaml", "second");
        second.commit().unwrap();

        assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        assert!(!directory.path().join("sample_b.fasta").exists());
        let manifest = read_manifest(directory.path(), "sample");
        assert_eq!(manifest.status, "complete");
        assert_eq!(manifest.files.len(), 2);
    }

    #[test]
    fn success_followed_by_no_products_removes_prior_fasta() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "product");
        stage_text(&mut first, "sample.stats.yaml", "success");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut second, "sample.stats.yaml", "no product");
        second.commit().unwrap();

        assert!(!directory.path().join("sample_gene.fasta").exists());
        assert!(directory.path().join("sample.stats.yaml").exists());
    }

    #[test]
    fn modified_owned_file_is_preserved_and_refused() {
        let directory = tempfile::tempdir().unwrap();
        let mut first = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut first, "sample_gene.fasta", "owned");
        stage_text(&mut first, "sample.stats.yaml", "stats");
        first.commit().unwrap();
        fs::write(directory.path().join("sample_gene.fasta"), "modified").unwrap();

        let error = OutputTransaction::begin(directory.path(), "sample")
            .err()
            .expect("modified file must be refused");
        assert!(error.to_string().contains("modified prior output"));
        assert_eq!(
            fs::read_to_string(directory.path().join("sample_gene.fasta")).unwrap(),
            "modified"
        );
    }

    #[test]
    fn legacy_collision_is_not_overwritten() {
        let directory = tempfile::tempdir().unwrap();
        let stats_path = directory.path().join("sample.stats.yaml");
        fs::write(&stats_path, "legacy").unwrap();

        let error = OutputTransaction::begin(directory.path(), "sample")
            .err()
            .expect("legacy stats must be refused during begin");
        assert!(error.to_string().contains("existing unowned stats"));
        assert_eq!(fs::read_to_string(stats_path).unwrap(), "legacy");
        assert!(!directory.path().join("sample.manifest.yaml").exists());
    }

    #[test]
    fn publishing_transaction_with_partial_files_is_invalidated_on_restart() {
        let directory = tempfile::tempdir().unwrap();
        let mut interrupted = OutputTransaction::begin(directory.path(), "sample").unwrap();
        stage_text(&mut interrupted, "sample_gene.fasta", "partial");
        stage_text(&mut interrupted, "sample.stats.yaml", "partial stats");
        interrupted
            .staged_outputs
            .sort_by(|left, right| left.receipt.path.cmp(&right.receipt.path));
        let receipts = interrupted.receipts();
        interrupted.failure_receipts = receipts.clone();
        interrupted
            .write_manifest("publishing", receipts, None)
            .unwrap();
        let published_basename = interrupted.staged_outputs[0].receipt.path.clone();
        publish_no_clobber(
            &interrupted.staged_outputs[0].staging_path,
            &directory.path().join(&published_basename),
        )
        .unwrap();
        interrupted.committed = true;
        drop(interrupted);
        assert_eq!(
            read_manifest(directory.path(), "sample").status,
            "publishing"
        );
        assert!(directory.path().join(&published_basename).exists());

        let mut restarted = OutputTransaction::begin(directory.path(), "sample").unwrap();
        assert!(!directory.path().join("sample_gene.fasta").exists());
        assert!(!directory.path().join("sample.stats.yaml").exists());
        stage_text(&mut restarted, "sample.stats.yaml", "restarted");
        restarted.commit().unwrap();
        assert_eq!(read_manifest(directory.path(), "sample").status, "complete");
    }

    #[test]
    fn overlapping_sample_and_gene_names_cannot_clobber_outputs() {
        let directory = tempfile::tempdir().unwrap();
        let collision_basename = "a_b_c.fasta";
        let mut first = OutputTransaction::begin(directory.path(), "a").unwrap();
        stage_text(&mut first, collision_basename, "first product");
        stage_text(&mut first, "a.stats.yaml", "first stats");
        first.commit().unwrap();

        let mut second = OutputTransaction::begin(directory.path(), "a_b").unwrap();
        stage_text(&mut second, collision_basename, "second product");
        stage_text(&mut second, "a_b.stats.yaml", "second stats");
        assert!(second.commit().is_err());

        assert_eq!(
            fs::read_to_string(directory.path().join(collision_basename)).unwrap(),
            "first product"
        );
        assert_eq!(read_manifest(directory.path(), "a").status, "complete");
    }

    #[test]
    fn sample_lock_is_held_for_transaction_lifetime() {
        let directory = tempfile::tempdir().unwrap();
        let transaction = OutputTransaction::begin(directory.path(), "sample").unwrap();
        let lock_path = directory.path().join("sample.lock");
        let second_handle = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&lock_path)
            .unwrap();
        assert!(second_handle.try_lock_exclusive().is_err());
        drop(transaction);
        second_handle.try_lock_exclusive().unwrap();
        second_handle.unlock().unwrap();
        assert!(lock_path.exists());
    }

    #[test]
    fn tampered_receipt_cannot_claim_unrelated_or_lock_files() {
        assert!(ensure_owned_output_basename("sample", "sample.manifest.yaml").is_err());
        for basename in ["sentinel.txt", "sample.lock"] {
            let directory = tempfile::tempdir().unwrap();
            let claimed_path = directory.path().join(basename);
            fs::write(&claimed_path, "keep").unwrap();
            let manifest = OutputManifest {
                schema_version: MANIFEST_SCHEMA_VERSION,
                producer: MANIFEST_PRODUCER.to_string(),
                sample: "sample".to_string(),
                run_id: "abc-123".to_string(),
                status: "complete".to_string(),
                files: vec![OutputReceipt {
                    path: basename.to_string(),
                    sha256: sha256_file(&claimed_path).unwrap(),
                }],
                failure_reason: None,
                staging_directory: None,
            };
            let manifest_path = directory.path().join("sample.manifest.yaml");
            serde_yaml_ng::to_writer(File::create(&manifest_path).unwrap(), &manifest).unwrap();

            assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
            assert_eq!(fs::read_to_string(claimed_path).unwrap(), "keep");
        }
    }

    #[test]
    fn malformed_complete_manifest_is_refused() {
        let directory = tempfile::tempdir().unwrap();
        let manifest = OutputManifest {
            schema_version: MANIFEST_SCHEMA_VERSION,
            producer: MANIFEST_PRODUCER.to_string(),
            sample: "sample".to_string(),
            run_id: "abc-123".to_string(),
            status: "complete".to_string(),
            files: Vec::new(),
            failure_reason: None,
            staging_directory: None,
        };
        let manifest_path = directory.path().join("sample.manifest.yaml");
        serde_yaml_ng::to_writer(File::create(&manifest_path).unwrap(), &manifest).unwrap();

        assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
        assert_eq!(read_manifest(directory.path(), "sample").status, "complete");
    }

    #[cfg(unix)]
    #[test]
    fn symlink_collision_is_preserved_and_refused() {
        use std::os::unix::fs::symlink;

        let directory = tempfile::tempdir().unwrap();
        let sentinel = directory.path().join("sentinel");
        fs::write(&sentinel, "keep").unwrap();
        symlink(&sentinel, directory.path().join("sample.stats.yaml")).unwrap();

        assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
        assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        assert!(!directory.path().join("sample.manifest.yaml").exists());
    }

    #[cfg(unix)]
    #[test]
    fn symlink_lock_and_manifest_are_refused() {
        use std::os::unix::fs::symlink;

        for basename in ["sample.lock", "sample.manifest.yaml"] {
            let directory = tempfile::tempdir().unwrap();
            let sentinel = directory.path().join("sentinel");
            fs::write(&sentinel, "keep").unwrap();
            symlink(&sentinel, directory.path().join(basename)).unwrap();

            assert!(OutputTransaction::begin(directory.path(), "sample").is_err());
            assert_eq!(fs::read_to_string(sentinel).unwrap(), "keep");
        }
    }
}

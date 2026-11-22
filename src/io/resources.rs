//! One numeric counter plan and one shared ledger, prepared before publication.
//! Preparing does not allocate tables, start workers or create temporary paths.
use super::*;

#[derive(Clone, Copy, Debug)]
pub(super) enum CounterLayout {
    Table {
        table_bytes: u64,
        packed: bool,
    },
    Parallel {
        table_bytes: u64,
        workers: usize,
    },
    Spill {
        table_bytes: u64,
        sort_bytes: u64,
        workers: usize,
        identity_worker: bool,
        owner_transfer_threshold: Option<usize>,
    },
}

#[derive(Clone, Copy, Debug)]
pub(super) struct CounterPlan {
    pub(super) k: usize,
    pub(super) record_bytes: usize,
    pub(super) disk_bytes: u64,
    pub(super) replay_bytes: u64,
    pub(super) graph_bytes: u64,
    pub(super) layout: CounterLayout,
}

/// Keep the validated counter plan, exact original args and parent ledger together.
pub(crate) struct PreparedIngestion<'a> {
    args: &'a Args,
    count: Option<(CounterPlan, Arc<ResourceBudget>)>,
    replay_requested: bool,
    workspace_root: Option<crate::kmer::budget::path::BudgetPath>,
}

fn requests_replay(args: &Args) -> bool {
    args.read_threading
        && (!args.pcr_panel.is_empty()
            || !args.pcr_panel_file.is_empty()
            || !args.pcr_primers.is_empty())
}

fn counter_plan(args: &Args, replay_requested: bool) -> Result<Option<CounterPlan>> {
    let mib = |value: u64, message: &'static str| value.checked_mul(1 << 20).context(message);
    let graph_bytes = mib(
        args.graph_buffer_budget_mib(),
        "Graph-buffer budget overflow",
    )?;
    let disk_bytes = mib(
        args.count_disk_budget_mib,
        "Count/replay disk budget overflow",
    )?;
    let replay_bytes = if replay_requested { 64 << 20 } else { 0 };
    let table = || mib(args.count_table_budget_mib, "Count table budget overflow");
    let sort = || mib(args.count_sort_buffer_mib, "External sort budget overflow");
    let layout = match args.count_backend() {
        "legacy" => None,
        "memory" | "packed" => Some(CounterLayout::Table {
            table_bytes: table()?,
            packed: args.count_backend() == "packed",
        }),
        "parallel" if args.threads <= 1 => Some(CounterLayout::Table {
            table_bytes: table()?,
            packed: false,
        }),
        "parallel" => Some(CounterLayout::Parallel {
            table_bytes: table()?,
            workers: args.threads - 1,
        }),
        "external" => Some(CounterLayout::Spill {
            table_bytes: 0,
            sort_bytes: sort()?,
            workers: 0,
            identity_worker: false,
            owner_transfer_threshold: None,
        }),
        "adaptive" | "adaptive-parallel" | "pipeline" | "hybrid" => {
            let pipelined = matches!(args.count_backend(), "pipeline" | "hybrid");
            Some(CounterLayout::Spill {
                table_bytes: table()?,
                sort_bytes: sort()?,
                workers: if pipelined {
                    args.threads.saturating_sub(2)
                } else if args.count_backend() == "adaptive-parallel" {
                    args.threads.saturating_sub(1)
                } else {
                    0
                },
                identity_worker: pipelined && args.threads > 1,
                owner_transfer_threshold: if args.count_backend() == "hybrid" && args.threads == 2 {
                    Some(usize::try_from(mib(
                        args.count_owner_threshold_mib,
                        "Count owner transfer threshold overflow",
                    )?)?)
                } else {
                    None
                },
            })
        }
        _ => anyhow::bail!("Unknown count backend: {}", args.count_backend()),
    };
    Ok(layout.map(|layout| CounterPlan {
        k: args.k,
        record_bytes: args.fastq_record_bytes,
        disk_bytes,
        replay_bytes,
        graph_bytes,
        layout,
    }))
}

/// Fit the existing count and assembly allowances inside a requested
/// process envelope after reserving runtime headroom.
pub(crate) fn fit_memory_budget(args: &mut Args, ceiling: u64) -> Result<()> {
    use crate::process_memory::{NORMAL_ADDRESS_SPACE_BYTES, NORMAL_PROCESS_HEADROOM_BYTES};

    let plan = counter_plan(args, requests_replay(args))?.context("Exact count plan is missing")?;
    let limits = plan.limits()?;
    let usable = ceiling.saturating_sub(NORMAL_PROCESS_HEADROOM_BYTES);
    if ceiling == NORMAL_ADDRESS_SPACE_BYTES && limits.memory_bytes <= usable {
        return Ok(());
    }
    let CounterLayout::Spill {
        table_bytes,
        sort_bytes,
        ..
    } = plan.layout
    else {
        anyhow::bail!("--memory-budget requires normal counting");
    };
    // Derive the fixed worker, record, lookup and replay allocations from
    // the same plan the counter uses. Split only the remaining memory.
    let fixed = limits.memory_bytes - table_bytes - sort_bytes - plan.graph_bytes;
    let available_mib = usable.saturating_sub(fixed) >> 20;
    ensure!(
        available_mib >= 10,
        "Memory budget is too small for counting and assembly with runtime headroom; raise --memory-budget or request more job memory"
    );
    // Use the default sort share up to 64 MiB. Split the rest at the
    // default 8:1 table-to-graph ratio, with up to one quarter per target.
    let sort_mib = (available_mib / 144).clamp(1, 64);
    let graph_mib = (available_mib - sort_mib) / 9;
    args.count_sort_buffer_mib = sort_mib;
    args.graph_buffer_budget_mib = Some(graph_mib);
    args.graph_buffer_target_mib = (graph_mib / 4).max(1);
    args.count_table_budget_mib = available_mib - sort_mib - graph_mib;
    Ok(())
}

/// Validate the numeric plan and observe process limits before output creation.
pub(crate) fn prepare_ingestion(args: &Args) -> Result<PreparedIngestion<'_>> {
    let replay_requested = requests_replay(args);
    let count = counter_plan(args, replay_requested)?
        .map(|plan| -> Result<_> {
            let budget = ResourceBudget::for_counting(plan.limits()?, 1)?;
            Ok((plan, budget))
        })
        .transpose()?;
    Ok(PreparedIngestion {
        args,
        count,
        replay_requested,
        workspace_root: None,
    })
}

impl<'a> PreparedIngestion<'a> {
    /// Borrow the original exact-run ledger before metadata resolution starts.
    /// This does not construct tables, workers or temporary directories.
    pub(crate) fn input_budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.count.as_ref().map(|(_, budget)| budget)
    }

    /// Only storage-capable counters or requested replay need a private root.
    pub(crate) fn needs_workspace(&self) -> bool {
        self.count.as_ref().is_some_and(|(plan, _)| {
            self.replay_requested || matches!(plan.layout, CounterLayout::Spill { .. })
        })
    }

    /// Preserve original arguments while directing temporary owners to a grant.
    pub(crate) fn with_workspace(mut self, root: crate::kmer::budget::path::BudgetPath) -> Self {
        self.workspace_root = Some(root);
        self
    }

    /// Refuse an oversized normal-run plan before output creation. Headroom is
    /// additional to every class in the shared ledger, not a second table limit.
    pub(crate) fn check_memory_envelope(&self, ceiling: u64, headroom: u64) -> Result<()> {
        let (plan, _) = self.count.as_ref().context("Exact count plan is missing")?;
        let required = plan
            .limits()?
            .memory_bytes
            .checked_add(headroom)
            .context("Production memory envelope overflow")?;
        ensure!(
            required <= ceiling,
            "Normal sPCR resource plan requires {required} bytes including process headroom, exceeding the {ceiling}-byte process envelope; use --memory-budget to plan a smaller run or request more job memory"
        );
        Ok(())
    }

    /// Borrow result ownership under the existing opt-in assembly allowance.
    pub(crate) fn result_budget(&self) -> Option<&Arc<ResourceBudget>> {
        self.count
            .as_ref()
            .and_then(|(plan, budget)| (plan.graph_bytes > 0).then_some(budget))
    }

    pub(super) fn start(
        self,
    ) -> Result<(
        &'a Args,
        Option<DirectCounter>,
        bool,
        Option<crate::kmer::budget::path::BudgetPath>,
    )> {
        let counter = self
            .count
            .map(|(plan, budget)| {
                plan.start(
                    budget,
                    self.workspace_root
                        .as_deref()
                        .or(self.args.count_temp_dir.as_deref()),
                )
            })
            .transpose()?;
        Ok((
            self.args,
            counter,
            self.replay_requested,
            self.workspace_root,
        ))
    }
}

/// Preserve the original combined ceiling; result limits add no allowance.
pub(super) fn resource_limits(
    mut classes: [u64; crate::kmer::budget::MEMORY_CLASSES],
    disk_bytes: u64,
    graph_bytes: u64,
) -> Result<ResourceLimits> {
    classes[MemoryClass::RetainedResults as usize] = 0;
    classes[MemoryClass::TargetGraphs as usize] = graph_bytes;
    let memory_bytes = classes.iter().try_fold(0u64, |sum, bytes| {
        sum.checked_add(*bytes)
            .context("Count/assembly resource budget overflow")
    })?;
    classes[MemoryClass::RetainedResults as usize] = memory_bytes;
    Ok(ResourceLimits {
        memory_bytes,
        process_overhead_bytes: 0,
        class_bytes: classes,
        temporary_disk_bytes: disk_bytes,
        temporary_open_files: 128,
    })
}

impl CounterPlan {
    /// Pure layout-specific class formulas, shared by startup and all constructors.
    pub(super) fn limits(self) -> Result<ResourceLimits> {
        let mut classes = [0; crate::kmer::budget::MEMORY_CLASSES];
        let key_bytes = u64::try_from(KMER_BATCH_CAPACITY * size_of::<CanonicalKey>())?;
        classes[MemoryClass::InputBuffer as usize] = u64::try_from(self.record_bytes)?
            .checked_mul(4)
            .context("FASTQ buffer budget overflow")?;
        classes[MemoryClass::WorkerQueues as usize] = key_bytes;
        classes[MemoryClass::ReplayBuffer as usize] = self.replay_bytes;
        match self.layout {
            CounterLayout::Table { table_bytes, .. } => {
                classes[MemoryClass::CountSlots as usize] = table_bytes
            }
            CounterLayout::Parallel {
                table_bytes,
                workers,
            } => {
                classes[MemoryClass::CountSlots as usize] = table_bytes;
                classes[MemoryClass::WorkerQueues as usize] = key_bytes
                    .checked_add(ParallelCountBuilder::reservation_bytes(workers)?.min(64 << 20))
                    .context("Counting worker reservation overflow")?;
            }
            CounterLayout::Spill {
                table_bytes,
                sort_bytes,
                workers,
                identity_worker,
                owner_transfer_threshold,
            } => {
                ensure!(
                    owner_transfer_threshold
                        .is_none_or(|bytes| bytes > 0 && workers == 0 && identity_worker),
                    "Owner transfer requires a positive threshold and one identity worker"
                );
                classes[MemoryClass::CountSlots as usize] = table_bytes;
                let owners = u64::try_from(workers.max(1))?;
                let future_workers = workers.max(usize::from(owner_transfer_threshold.is_some()));
                let pool_bytes = if future_workers == 0 {
                    0
                } else {
                    AdaptiveParallelCountBuilder::reservation_bytes(future_workers)?.min(64 << 20)
                };
                let identity_bytes = if identity_worker {
                    IdentityWriter::reservation_bytes()?
                } else {
                    0
                };
                classes[MemoryClass::WorkerQueues as usize] =
                    ExternalCountBuilder::run_slot_bytes()?
                        .checked_mul(owners)
                        .and_then(|n| n.checked_add(key_bytes))
                        .and_then(|n| n.checked_add(pool_bytes))
                        .and_then(|n| n.checked_add(identity_bytes))
                        .context("Count worker reservation overflow")?;
                classes[MemoryClass::SortBuffer as usize] = sort_bytes;
                classes[MemoryClass::LookupCache as usize] =
                    ExternalCountBuilder::cache_limit_bytes()
                        .checked_mul(owners)
                        .context("Count page reservation overflow")?;
            }
        }
        resource_limits(classes, self.disk_bytes, self.graph_bytes)
    }

    /// Start actual owners only after startup/publication can use the same ledger.
    pub(super) fn start(
        self,
        budget: Arc<ResourceBudget>,
        directory: Option<&std::path::Path>,
    ) -> Result<DirectCounter> {
        let memory_peak = budget.begin_memory_peak()?;
        let disk_peak = budget.begin_disk_peak()?;
        let mut owned_directory = None;
        let mut identity_worker = false;
        let mut transfer = None;
        let builder = match self.layout {
            CounterLayout::Table {
                table_bytes: _,
                packed,
            } => {
                if packed {
                    DirectBuilder::Packed(PackedCountBuilder::with_budget(
                        KmerSize::new(self.k)?,
                        Arc::clone(&budget),
                        1024,
                    )?)
                } else {
                    DirectBuilder::Memory(MemoryCountBuilder::with_budget(
                        KmerSize::new(self.k)?,
                        Arc::clone(&budget),
                        1024,
                    )?)
                }
            }
            CounterLayout::Parallel { workers, .. } => {
                DirectBuilder::Parallel(ParallelCountBuilder::new(
                    KmerSize::new(self.k)?,
                    Arc::clone(&budget),
                    workers,
                    1024,
                )?)
            }
            CounterLayout::Spill {
                table_bytes,
                sort_bytes,
                workers,
                identity_worker: identity,
                owner_transfer_threshold,
            } => {
                let fallback;
                let directory = match directory {
                    Some(directory) => directory,
                    None => {
                        fallback = std::env::temp_dir();
                        &fallback
                    }
                };
                let owned = tempfile::Builder::new()
                    .prefix(".sharkmer-count-run-")
                    .tempdir_in(directory)
                    .context("Unable to create external-count run directory")?;
                let builder = if workers > 0 {
                    DirectBuilder::AdaptiveParallel(AdaptiveParallelCountBuilder::adaptive(
                        owned.path(),
                        KmerSize::new(self.k)?,
                        Arc::clone(&budget),
                        workers,
                        usize::try_from(sort_bytes)?,
                        1024,
                    )?)
                } else if table_bytes == 0 {
                    DirectBuilder::External(ExternalCountBuilder::new(
                        owned.path(),
                        KmerSize::new(self.k)?,
                        Arc::clone(&budget),
                        usize::try_from(sort_bytes)?,
                        1024,
                    )?)
                } else {
                    DirectBuilder::Adaptive(AdaptiveCountBuilder::new(
                        owned.path(),
                        KmerSize::new(self.k)?,
                        Arc::clone(&budget),
                        usize::try_from(sort_bytes)?,
                        1024,
                    )?)
                };
                identity_worker = identity;
                transfer = owner_transfer_threshold;
                owned_directory = Some(owned);
                builder
            }
        };
        let mut counter = DirectCounter::with_builder(
            self.k,
            self.record_bytes,
            budget,
            builder,
            memory_peak,
            disk_peak,
        )?;
        counter.owned_count_directory = owned_directory;
        if identity_worker {
            counter.identity_writer = Some(IdentityWriter::new(&counter.budget)?);
        }
        counter.owner_transfer_threshold = transfer;
        Ok(counter)
    }
}

#[cfg(test)]
mod reference;

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    fn args(backend: &str, threads: usize, graph: u64, replay: bool) -> Args {
        let mut values = vec![
            "sharkmer".to_owned(),
            "--count-backend".to_owned(),
            backend.to_owned(),
            "--threads".to_owned(),
            threads.to_string(),
            "--chunks".to_owned(),
            "0".to_owned(),
            "-k".to_owned(),
            "19".to_owned(),
            "--count-table-budget-mib".to_owned(),
            "1".to_owned(),
            "--count-sort-buffer-mib".to_owned(),
            "1".to_owned(),
            "--count-disk-budget-mib".to_owned(),
            "32".to_owned(),
            "--fastq-record-bytes".to_owned(),
            "4096".to_owned(),
            "--graph-buffer-budget-mib".to_owned(),
            graph.to_string(),
            "--pcr-primers".to_owned(),
            "name=target,forward=AAAA,reverse=TTTT".to_owned(),
        ];
        if replay {
            values.push("--read-threading".to_owned());
        }
        Args::parse_from(values)
    }

    #[test]
    fn normal_plan_fits_twelve_gib_with_headroom_and_refuses_a_smaller_envelope() {
        let args = Args::parse_from(["sharkmer", "-s", "sample"]);
        let prepared = prepare_ingestion(&args).unwrap();
        assert!(prepared.result_budget().is_some());
        prepared.check_memory_envelope(12 << 30, 2 << 30).unwrap();
        assert!(prepared.check_memory_envelope(8 << 30, 2 << 30).is_err());
        assert!(prepared.check_memory_envelope(u64::MAX, u64::MAX).is_err());
    }

    #[test]
    fn explicit_memory_budget_preserves_the_default_plan_and_fits_other_envelopes() {
        let args = || {
            Args::parse_from([
                "sharkmer",
                "-s",
                "sample",
                "--threads",
                "8",
                "--pcr-panel",
                "cnidaria",
                "--read-threading",
            ])
        };
        let standard = counter_plan(&args(), true)
            .unwrap()
            .unwrap()
            .limits()
            .unwrap();
        for mib in [2304, 8192, 12288, 16384] {
            let mut args = args();
            fit_memory_budget(&mut args, mib << 20).unwrap();
            let prepared = prepare_ingestion(&args).unwrap();
            prepared.check_memory_envelope(mib << 20, 2 << 30).unwrap();
            let budget = prepared.result_budget().unwrap();
            assert!(args.graph_buffer_target_mib > 0);
            assert!(args.graph_buffer_target_mib <= args.graph_buffer_budget_mib());
            if mib == 12288 {
                let effective = prepared.count.as_ref().unwrap().0.limits().unwrap();
                assert_eq!(effective.class_bytes, standard.class_bytes);
                assert_eq!(args.graph_buffer_target_mib, 256);
            } else {
                let table_bytes = budget.class_limit(MemoryClass::CountSlots);
                assert_eq!(
                    table_bytes < standard.class_bytes[MemoryClass::CountSlots as usize],
                    mib < 12288
                );
            }
        }
        assert!(fit_memory_budget(&mut args(), 2048 << 20).is_err());
    }

    #[test]
    fn startup_resource_limits_match_frozen_formulas_for_every_backend_and_role_layout() {
        for backend in [
            "memory",
            "packed",
            "parallel",
            "external",
            "adaptive",
            "adaptive-parallel",
            "pipeline",
            "hybrid",
        ] {
            for threads in [1, 2, 4, 8] {
                for graph in [0, 8] {
                    for replay in [false, true] {
                        let args = args(backend, threads, graph, replay);
                        let prepared = prepare_ingestion(&args).unwrap();
                        let (plan, budget) = prepared.count.as_ref().unwrap();
                        let actual = plan.limits().unwrap();
                        let replay_bytes = if replay { 64 << 20 } else { 0 };
                        let graph_bytes = graph << 20;
                        let expected = match backend {
                            "memory" | "packed" => reference::table(
                                19,
                                1 << 20,
                                4096,
                                32 << 20,
                                replay_bytes,
                                backend == "packed",
                                graph_bytes,
                            ),
                            "parallel" => reference::parallel(
                                19,
                                1 << 20,
                                4096,
                                threads,
                                32 << 20,
                                replay_bytes,
                                graph_bytes,
                            ),
                            _ => reference::spill(
                                19,
                                if backend == "external" { 0 } else { 1 << 20 },
                                1 << 20,
                                4096,
                                32 << 20,
                                std::path::Path::new("unused"),
                                replay_bytes,
                                if matches!(backend, "pipeline" | "hybrid") {
                                    threads.saturating_sub(2)
                                } else if backend == "adaptive-parallel" {
                                    threads.saturating_sub(1)
                                } else {
                                    0
                                },
                                matches!(backend, "pipeline" | "hybrid") && threads > 1,
                                (backend == "hybrid" && threads == 2).then_some(16 << 20),
                                graph_bytes,
                            ),
                        }
                        .unwrap();
                        assert_eq!(
                            actual.class_bytes, expected.class_bytes,
                            "{backend} t{threads} graph{graph} replay{replay}"
                        );
                        assert_eq!(actual.memory_bytes, expected.memory_bytes);
                        assert_eq!(actual.temporary_disk_bytes, expected.temporary_disk_bytes);
                        assert_eq!(
                            actual.process_overhead_bytes,
                            expected.process_overhead_bytes
                        );
                        assert_eq!(actual.temporary_open_files, expected.temporary_open_files);
                        assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
                        assert_eq!(budget.snapshot().reserved_temporary_open_files, 0);
                        assert_eq!(prepared.result_budget().is_some(), graph > 0);
                    }
                }
            }
        }
        let args = args("legacy", 2, 0, true);
        let prepared = prepare_ingestion(&args).unwrap();
        assert!(prepared.count.is_none());
        assert!(prepared.result_budget().is_none());
    }

    #[test]
    fn prepared_counter_reuses_startup_ledger_without_early_tables_workers_or_directories() {
        let directory = tempfile::tempdir().unwrap();
        for backend in [
            "memory",
            "packed",
            "parallel",
            "external",
            "adaptive",
            "adaptive-parallel",
            "pipeline",
            "hybrid",
        ] {
            let mut args = args(backend, 2, 8, false);
            let missing = directory.path().join(backend);
            args.count_temp_dir = Some(missing.clone());
            let prepared = prepare_ingestion(&args).unwrap();
            let budget = Arc::clone(prepared.result_budget().unwrap());
            assert!(!missing.exists());
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            // A finished startup owner must not dominate the count phase metric.
            let earlier = budget.class_limit(MemoryClass::RetainedResults) - 1;
            drop(
                budget
                    .reserve(MemoryClass::RetainedResults, earlier)
                    .unwrap(),
            );
            let startup = directory.path().join(format!("{backend}-output"));
            drop(
                crate::output::OutputTransaction::begin_with_budget(
                    &startup,
                    "sample",
                    Some(&budget),
                )
                .unwrap(),
            );
            let earlier_disk = budget.snapshot().peak_reserved_temporary_disk_bytes;
            assert!(earlier_disk > 0);
            assert_eq!(budget.snapshot().reserved_temporary_disk_bytes, 0);
            std::fs::create_dir(&missing).unwrap();
            let (_, counter, _, _) = prepared.start().unwrap();
            let mut counter = counter.unwrap();
            assert!(Arc::ptr_eq(&budget, &counter.budget));
            assert!(budget.snapshot().reserved_memory_bytes > 0);
            let measured = counter.memory_peak.take().unwrap().finish().unwrap();
            assert!(measured > 0 && measured < earlier);
            let measured_disk = counter.disk_peak.take().unwrap().finish().unwrap();
            assert!(measured_disk < earlier_disk);
            drop(counter);
            assert!(std::fs::read_dir(&missing).unwrap().next().is_none());
            assert_eq!(budget.snapshot().reserved_memory_bytes, 0);
            assert_eq!(budget.begin_memory_peak().unwrap().finish().unwrap(), 0);
            assert_eq!(budget.begin_disk_peak().unwrap().finish().unwrap(), 0);
            assert!(budget.snapshot().peak_reserved_temporary_disk_bytes >= earlier_disk);
            assert!(budget.snapshot().peak_reserved_memory_bytes >= earlier);
        }
    }

    #[test]
    fn invalid_startup_plans_refuse_before_creating_output_or_temporary_paths() {
        let directory = tempfile::tempdir().unwrap();
        let mut args = args("hybrid", 2, 8, false);
        args.outdir = directory.path().join("output");
        args.count_temp_dir = Some(directory.path().join("count"));
        args.count_table_budget_mib = u64::MAX;
        assert!(prepare_ingestion(&args).is_err());
        assert!(!args.outdir.exists());
        assert!(!args.count_temp_dir.as_ref().unwrap().exists());
        args.count_table_budget_mib = 1;
        args.count_owner_threshold_mib = 0;
        assert!(prepare_ingestion(&args).is_err());
        assert!(!args.outdir.exists());
    }
}

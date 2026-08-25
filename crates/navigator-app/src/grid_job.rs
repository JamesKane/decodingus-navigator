//! The unit driver: the work that a volunteer node does for one Grid work unit.
//!
//! ```text
//! claim ─► preflight ─► fetch from ENA ─► import ─► analyze ─► digest ─► submit ─► clean
//! ```
//!
//! Each step is a method that already exists somewhere in this crate. This module puts them in
//! order, reports progress, and makes sure that a lease always ends. Design
//! `documents/design/distributed-compute-grid.md` §7.2.
//!
//! # A lease must always end
//!
//! A node that stops in the middle keeps a unit out of the catalogue until the lease time ends. So
//! each path out of [`App::run_grid_unit`] gives the lease back. Success closes the lease through
//! the submit call, and each failure sends a release call. The AppView also has a reaper for the node that
//! disappears. But a node that is still alive must not need it.
//!
//! # What this node can do
//!
//! [`supported_data_kinds`] gives the list that the node advertises, and it is the only place that
//! decides. A node never receives work that this module can not do, because the AppView selects
//! work with that same list. Today the list holds `CRAM` only. See [`map_reads`] for what a FASTQ
//! unit needs.

use super::*;
use crate::ena::{self, ManifestFile};
use crate::grid::ClaimedUnit;
use du_domain::fed::Provenance;
use navigator_analysis::CancelToken;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// The stage that a node is on. The node sends this with each heartbeat, and the fleet view of the
/// AppView shows it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridStage {
    Fetch,
    Import,
    Map,
    Analyze,
    Ancestry,
    Publish,
    Submit,
}

impl GridStage {
    pub fn as_str(self) -> &'static str {
        match self {
            GridStage::Fetch => "fetch",
            GridStage::Import => "import",
            GridStage::Map => "map",
            GridStage::Analyze => "analyze",
            GridStage::Ancestry => "ancestry",
            GridStage::Publish => "publish",
            GridStage::Submit => "submit",
        }
    }
}

/// What one unit gave.
#[derive(Debug, Clone)]
pub struct UnitOutcome {
    pub sample_accession: String,
    /// The submission id, when the node sent a result.
    pub submission_id: Option<i64>,
    /// Why the unit did not finish. The node then gives the lease back.
    pub error: Option<String>,
}

/// How often a node tells the AppView that it is alive, while it works on a unit.
///
/// The value is far below the shortest lease. A node that stops between two beats is still inside
/// its lease, so a lost beat costs nothing.
const HEARTBEAT_EVERY: std::time::Duration = std::time::Duration::from_secs(60);

/// How a node contributes.
#[derive(Debug, Clone)]
pub struct GridJobParams {
    /// How many units to take in one claim call.
    pub max_units: i32,
    /// How long to hold each lease, in seconds.
    pub lease_secs: i64,
    /// Where the files of a unit go. Each unit gets its own directory below this one.
    pub scratch_root: PathBuf,
    /// The build that the analysis is against, such as `chm13v2.0`.
    pub reference_build: String,
}

/// The data kinds that this node can process. These are the kinds that it advertises.
///
/// A `CRAM` unit arrives with an alignment that a laboratory already made. The node imports that
/// file and analyzes it.
///
/// A `FASTQ` unit holds reads only, so the node maps them first. See [`map_unit_reads`].
///
/// This list is the one place that decides. The AppView selects work with it, so a node never
/// receives work that this module can not do.
pub fn supported_data_kinds() -> Vec<String> {
    vec!["CRAM".to_string(), "FASTQ".to_string()]
}

/// Sort the read files of a FASTQ unit into the first mate file, the second mate file, and the
/// reads with no mate.
///
/// ENA gives the mate number in the file name, as `_1` and `_2` before the extension. A run with
/// one file only is a set of reads with no mate, and a long-read run is always such a set.
fn split_mates(files: &[PathBuf]) -> (Option<PathBuf>, Option<PathBuf>, Vec<PathBuf>) {
    let (mut r1, mut r2, mut singles) = (None, None, Vec::new());
    for f in files {
        let name = f.file_name().and_then(|n| n.to_str()).unwrap_or_default();
        // `_1.` and not `_1`: a run accession such as `ERR1_1.fastq.gz` must not match on the
        // accession itself.
        if name.contains("_1.") && r1.is_none() {
            r1 = Some(f.clone());
        } else if name.contains("_2.") && r2.is_none() {
            r2 = Some(f.clone());
        } else {
            singles.push(f.clone());
        }
    }
    (r1, r2, singles)
}

/// Map the reads of a FASTQ unit to the target build, and give back a finished CRAM file.
///
/// # Why this does not call the realignment job
///
/// The realignment job (`realign_job`) does the same four operations, and the first plan was to
/// call it here. A reading of that module changed the plan. That module holds its stages together
/// with the machinery that continues a job which stopped: `Resumed`, `ScratchState`, and the rules
/// about which file each stage may remove. A comment in that file records a fault in exactly those
/// rules. That fault destroyed a 59 GB file and about four hours of work.
///
/// A Grid unit wants none of that machinery. It has no source alignment, so there is no revert
/// stage. It does not continue a job that stopped, because a unit that fails gives its lease back
/// and another node takes it from the start. And it registers no alignment against a source row.
/// The only common part is the four operations below, and each one is already public.
///
/// So this function calls those four operations directly. That leaves the realignment module
/// exactly as `v0.1.0-alpha.17` validated it on a full genome. The other method was to divide that
/// module along its most dangerous line, with no way to run that validation again here.
async fn map_unit_reads(
    app: &App,
    files: &[PathBuf],
    dir: &Path,
    target_build: &str,
    cancel: &CancelToken,
    report: &mut (dyn FnMut(GridStage, &str) + Send),
) -> Result<PathBuf, AppError> {
    use navigator_analysis::postprocess::{self, MarkDupParams, SortParams};

    let (r1, r2, singles) = split_mates(files);
    let paired = r1.is_some() && r2.is_some();
    // A short-read preset for a set of reads with a mate, and a long-read preset for a set with no
    // mate. A map of long reads under a short-read preset does not fail. It gives alignments that
    // look correct and are wrong.
    let preset = if paired {
        navigator_align::Preset::ShortRead
    } else {
        navigator_align::Preset::MapHifi
    };

    report(GridStage::Map, "reference");
    let reference = app.resolve_reference(target_build, &mut |_, _| {}).await?;

    report(GridStage::Map, "index");
    let index = {
        let (build, reference) = (target_build.to_string(), reference.clone());
        let batch = navigator_align::batch::BatchSize::for_this_machine();
        tokio::task::spawn_blocking(move || {
            navigator_align::index::ensure_index(
                &navigator_align::index::cache_root(),
                &build,
                &reference,
                preset,
                batch,
                &mut |_, _| {},
            )
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??
    };

    let mapped = dir.join("mapped.bam");
    let sorted = dir.join("sorted.bam");
    let marked = dir.join("marked.bam");
    let output = dir.join("aligned.cram");

    report(GridStage::Map, "map");
    {
        let (out, work) = (mapped.clone(), dir.join("map"));
        let token = cancel.clone();
        let map_params = navigator_align::MapParams {
            preset,
            threads: 0,
            read_group: None,
            format: navigator_align::OutputFormat::Bam,
            reference: None,
        };
        let (r1c, r2c, singlesc) = (r1.clone(), r2.clone(), singles.clone());
        tokio::task::spawn_blocking(move || -> Result<(), AppError> {
            let cancelled = move || token.is_cancelled();
            if let (Some(a), Some(b)) = (&r1c, &r2c) {
                navigator_align::map_pairs(&index, a, b, &out, &work, &map_params, &cancelled, &mut |_, _, _| {})?;
            } else {
                let single = singlesc
                    .first()
                    .or(r1c.as_ref())
                    .ok_or_else(|| AppError::Import("the unit holds no read file".into()))?;
                navigator_align::map_reads(&index, single, &out, &work, &map_params, &cancelled, &mut |_, _, _| {})?;
            }
            Ok(())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
    }
    // The reads have no more use, and a set of read files for a whole genome is tens of GB.
    for f in files {
        let _ = std::fs::remove_file(f);
    }

    report(GridStage::Map, "sort");
    {
        let (input, out, work) = (mapped.clone(), sorted.clone(), dir.join("sort"));
        let token = cancel.clone();
        tokio::task::spawn_blocking(move || {
            postprocess::sort_alignment(&input, &out, &work, &SortParams::default(), &token, &mut |_| {})
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
    }
    let _ = std::fs::remove_file(&mapped);

    report(GridStage::Map, "duplicates");
    {
        let (input, out) = (sorted.clone(), marked.clone());
        let token = cancel.clone();
        // A long-read library usually needs no PCR step, and two long reads rarely have the same
        // end points. So a mark on those reads removes real coverage.
        let md_params = MarkDupParams {
            enabled: paired,
            ..Default::default()
        };
        tokio::task::spawn_blocking(move || {
            postprocess::mark_duplicates(&input, &out, &md_params, &token, &mut |_| {})
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
    }
    let _ = std::fs::remove_file(&sorted);

    report(GridStage::Map, "compress");
    let finalized = {
        let (input, out) = (marked.clone(), output.clone());
        tokio::task::spawn_blocking(move || postprocess::finalize_bam(&input, &out))
            .await
            .map_err(|e| AppError::Join(e.to_string()))??
    };
    Ok(finalized.bam)
}

/// The values that go into the digest of a result.
///
/// Each one is optional. A sample can have no Y chromosome. The autosomal consensus of a fresh
/// sample can be absent. An absent value is not a failure of the unit, and the AppView compares two
/// absent values as equal.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UnitResults {
    pub sex: Option<String>,
    pub y_terminal: Option<String>,
    pub ancestry_superpop_argmax: Option<String>,
    pub coverage_mean: Option<f64>,
    pub callable_fraction: Option<f64>,
}

/// Build the digest that the node signs and sends.
///
/// The values are **raw**. The AppView puts the continuous values into groups when it compares two
/// results. A node that made the groups itself would put the group rule into two repositories.
///
/// There is no `mt_terminal` field. `App::analyze_biosample` does not give an mtDNA value, because
/// that value is not final on CHM13, and the Grid analyzes against CHM13. A digest can not ask for
/// a value that the analysis does not make. See design §12.3.
pub fn build_digest(
    sample_accession: &str,
    reference_build: &str,
    stack_version: &str,
    aligner: Option<&str>,
    r: &UnitResults,
) -> serde_json::Value {
    let mut calls = serde_json::Map::new();
    if let Some(v) = &r.sex {
        calls.insert("sex".into(), serde_json::json!(v));
    }
    if let Some(v) = &r.y_terminal {
        calls.insert("y_terminal".into(), serde_json::json!(v));
    }
    if let Some(v) = &r.ancestry_superpop_argmax {
        calls.insert("ancestry_superpop_argmax".into(), serde_json::json!(v));
    }
    if let Some(v) = r.coverage_mean {
        calls.insert("coverage_mean".into(), serde_json::json!(v));
    }
    if let Some(v) = r.callable_fraction {
        calls.insert("callable_fraction".into(), serde_json::json!(v));
    }
    serde_json::json!({
        "unit": sample_accession,
        "reference_build": reference_build,
        "stack_version": stack_version,
        "aligner": aligner,
        "calls": serde_json::Value::Object(calls),
    })
}

/// The key that the `provenance` block takes in a published record.
///
/// It is the serde name of the field on the record types of `du_domain::fed`. A test below makes a
/// record with the typed method and then reads the key back. So a test checks this value against
/// the type, and this value is not an assumption.
const PROVENANCE_KEY: &str = "provenance";

/// Put the Grid provenance block into a record that a builder already made.
///
/// The record builders of `publish.rs` serve the ordinary path, where a user publishes a record
/// about their own genome. Such a record carries no provenance, and those builders have thirteen
/// call sites. A new argument on each builder would put a `None` at each of those call sites, for a
/// value that only the Grid supplies. That `None` would mean nothing to any of them.
///
/// So the Grid adds the block after the builder finishes. The value comes from the typed
/// [`Provenance`] of `du_domain::fed`, so the shape and the field names come from the shared
/// contract. Only the key is a string here, and a test checks that string against the type.
fn attach_provenance(mut value: serde_json::Value, p: &Provenance) -> Result<serde_json::Value, AppError> {
    let block = serde_json::to_value(p).map_err(|e| AppError::Import(e.to_string()))?;
    match &mut value {
        serde_json::Value::Object(map) => {
            map.insert(PROVENANCE_KEY.to_string(), block);
            Ok(value)
        }
        _ => Err(AppError::Import("a published record must be a JSON object".into())),
    }
}

/// The provenance of a result that this node computed for the Grid.
fn grid_provenance(did: &str, reference_build: &str, aligner: Option<&str>) -> Provenance {
    Provenance::new(
        did,
        "navigator",
        env!("CARGO_PKG_VERSION"),
        reference_build,
        // How the input arrived. It names the public origin, so a later reader can honour any term
        // that the study of that sample sets.
        "ena:read_run",
    )
    .with_aligner(aligner.map(str::to_string))
}

/// The accession as a directory name, or `None` when it is not a safe name.
///
/// The accession arrives from the AppView, and the code makes a path from it and later **removes
/// that path and everything below it**. A value with `..` in it would leave the scratch directory,
/// and the remove would then delete a directory of the user.
///
/// The AppView is not an attacker. But a value from a server is still a value from outside, and a
/// recursive delete is on the other side of it. Each archive that this code reads gives an
/// accession of the form `[A-Za-z0-9_.-]+`, so this check refuses nothing real.
fn safe_dir_name(accession: &str) -> Option<&str> {
    let a = accession.trim();
    let ok = !a.is_empty()
        && a.len() <= 64
        && a != "."
        && a != ".."
        && !a.contains("..")
        && a.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'));
    ok.then_some(a)
}

/// The value that the digest carries for the sex of a sample.
///
/// The mapping is explicit, and it does not use the `Debug` form of the enum. The AppView compares
/// this value between two nodes. A new name for a variant would change the value, and no other
/// thing would change. Two versions of Navigator would then disagree about one sample, and no
/// reader could find the cause.
///
/// An uncertain result gives `None`, and the digest then holds no sex value. That is honest: two
/// nodes that both could not tell agree, and a node that could tell does not agree with one that
/// could not.
fn sex_for_digest(sex: navigator_analysis::sex::InferredSex) -> Option<String> {
    use navigator_analysis::sex::InferredSex;
    match sex {
        InferredSex::Male => Some("XY".to_string()),
        InferredSex::Female => Some("XX".to_string()),
        InferredSex::Unknown => None,
    }
}

/// How many sequencing runs a manifest holds.
///
/// The count is of the **run accessions** and not of the files. A single run with two mates gives
/// two files. ENA frequently gives a third file, for the reads of that run that lost their mate. A
/// count of files would refuse such a run as though it held three runs.
///
/// An entry with no run accession counts as one run. So a manifest with no accession at all is one
/// run. That is the safe reading: the node does the work, and it does not refuse a unit because a
/// field was empty.
fn runs_in(manifest: &[ManifestFile]) -> usize {
    let named: std::collections::BTreeSet<&str> = manifest
        .iter()
        .map(|m| m.run_accession.trim())
        .filter(|a| !a.is_empty())
        .collect();
    named.len().max(1)
}

/// A short class for the release call, from the error of a unit.
///
/// The message of an error can hold a local file path. The server keeps this value and the node
/// signs it, so it must hold no name from the machine of the volunteer.
fn release_reason(e: &AppError) -> &'static str {
    let text = e.to_string().to_ascii_lowercase();
    if text.contains("cancel") {
        "stopped"
    } else if text.contains("checksum") {
        "checksum"
    } else if text.contains("room") || text.contains("space") {
        "disk"
    } else if text.contains("merge runs") {
        "unsupported"
    } else if text.contains("analysis") {
        "analysis"
    } else {
        "error"
    }
}

/// The primary data file of a unit: the alignment for a CRAM unit, or the first read file for a
/// FASTQ unit. An index file is never the primary file.
fn primary_file<'a>(manifest: &'a [ManifestFile], files: &'a [PathBuf]) -> Option<&'a PathBuf> {
    manifest
        .iter()
        .zip(files)
        .find(|(m, _)| matches!(m.format.as_str(), "CRAM" | "BAM" | "FASTQ"))
        .map(|(_, p)| p)
}

/// The free space of the volume that holds `dir`, in bytes. Zero means that the code can not
/// measure it.
pub fn free_space_for(dir: &Path) -> u64 {
    crate::realign_job::free_space(dir)
}

/// The physical memory of this machine, in bytes. Zero means that the code can not measure it.
pub fn machine_memory_bytes() -> u64 {
    navigator_align::batch::detect_memory().map(|m| m.total).unwrap_or(0)
}

/// Remove each unit directory that is older than `max_age`.
///
/// A unit that the user stopped keeps its files, so that the transfer can continue. This
/// removes the directories that no run continued. Call it when a node starts.
pub async fn sweep_old_scratch(scratch_root: &Path, max_age: std::time::Duration) -> usize {
    let Ok(mut entries) = tokio::fs::read_dir(scratch_root).await else {
        return 0;
    };
    let mut removed = 0;
    while let Ok(Some(entry)) = entries.next_entry().await {
        let old = entry
            .metadata()
            .await
            .ok()
            .and_then(|m| m.modified().ok())
            .and_then(|t| t.elapsed().ok())
            .is_some_and(|age| age > max_age);
        if old && tokio::fs::remove_dir_all(entry.path()).await.is_ok() {
            removed += 1;
        }
    }
    removed
}

impl App {
    /// Do the work of one unit, and always give the lease back.
    ///
    /// The node reports each stage through `report`, and the caller sends those to the AppView as a
    /// heartbeat. A heartbeat that answers `false` means that this node lost the lease. The node
    /// then stops, because it receives no credit for more work on that unit.
    pub async fn run_grid_unit(
        &self,
        unit: &ClaimedUnit,
        params: &GridJobParams,
        cancel: &CancelToken,
        report: &mut (dyn FnMut(GridStage, &str) + Send),
    ) -> UnitOutcome {
        let mut outcome = UnitOutcome {
            sample_accession: unit.sample_accession.clone(),
            submission_id: None,
            error: None,
        };
        // The subject that the unit makes, so the code can remove it at the end. See
        // [`discard_unit_subject`].
        let subject = Arc::new(Mutex::new(None));
        let Some(safe_name) = safe_dir_name(&unit.sample_accession) else {
            let e = AppError::Import(format!(
                "the accession \"{}\" is not a name that this node will make a directory from",
                unit.sample_accession
            ));
            outcome.error = Some(e.to_string());
            let _ = self.grid_release(unit.lease_id, release_reason(&e)).await;
            return outcome;
        };
        let dir = params.scratch_root.join(safe_name);

        // The stage that the beat reports. The work writes this value and the beat reads it. So
        // the AppView shows the stage that the node is on now, and not the first stage.
        let stage = Arc::new(Mutex::new(GridStage::Fetch));
        // Set by the beat when the AppView says that another node holds this lease. It separates
        // that event from a stop by the user, and the two want different treatment of the files.
        let lease_lost = Arc::new(std::sync::atomic::AtomicBool::new(false));

        // The work and the beat run together. Neither one can be a separate task, because both use
        // `&self`, and a task needs a value that lives for the whole program. `select!` needs no
        // such value: it drives two futures that borrow the same data.
        let result = {
            let stage_for_work = Arc::clone(&stage);
            let mut record = |s: GridStage, detail: &str| {
                if let Ok(mut cur) = stage_for_work.lock() {
                    *cur = s;
                }
                report(s, detail);
            };
            let work = self.grid_unit_inner(unit, params, &dir, cancel, &mut record, &subject);
            let beat = self.beat_while_working(unit.lease_id, &stage, cancel, &lease_lost);
            tokio::pin!(work);
            tokio::pin!(beat);
            tokio::select! {
                r = &mut work => r,
                lost = &mut beat => {
                    // The beat already set the token. Wait for the work to see it and return.
                    //
                    // A `select!` that ends here would drop the work in the middle of an `await`.
                    // The heavy stages run in `spawn_blocking`. A dropped handle does not stop such
                    // a task. So the sort or the duplicate mark would continue, while the code
                    // below removes the directory that it writes into. The code waits instead, and
                    // the work stops at its next test of the token.
                    // Take the result of the work when it has one.
                    //
                    // The work can reach the submit call inside this window and complete. A node
                    // that reported a failure then would lose the credit for work that it
                    // finished. It would also call release on a lease that the submit call had
                    // already closed.
                    match (&mut work).await {
                        Ok(id) => Ok(id),
                        Err(_) => Err(lost),
                    }
                }
            }
        };

        match result {
            Ok(id) => outcome.submission_id = Some(id),
            Err(e) => {
                outcome.error = Some(e.to_string());
                // The unit goes back to the catalogue at once. Without this call, it waits for the
                // full lease time, and no other node can take it.
                //
                // The reason that goes to the server is a short class and not the full message.
                // The full message can hold a local path, because `ena` puts the path of a file in
                // the text of an I/O error. The server keeps the reason, and the node signs it. So
                // a path in that text would send the directory names of a volunteer to a public
                // service. The full message stays here, in `outcome.error`.
                let _ = self.grid_release(unit.lease_id, release_reason(&e)).await;
            }
        }
        // A unit is work, and it is not the data of the user. Remove the subject before the files,
        // because the subject names those files.
        let made = subject.lock().ok().and_then(|g| *g);
        if let Some(guid) = made {
            self.discard_unit_subject(guid).await;
        }

        // Keep the files of a unit that **the user** stopped. Remove them in each other case.
        //
        // `ena` can continue a transfer that stopped. It keeps each `.part` file, and it reads the
        // md5 state of that prefix again. A remove of the directory here would make all of that
        // work impossible: a stop at 25 GB of a 30 GB file would lose those 25 GB.
        //
        // A unit that failed is different. Another node takes it, and this node may never see it
        // again. So those files stay on the disk with no purpose, and a whole genome is tens of GB.
        //
        // [`sweep_old_scratch`] removes a directory that a stop left, after some days. Without that
        // step, a user who stops a run and never continues it keeps those files for ever.
        //
        // The beat also cancels this token, when another node takes the lease. That is not a stop
        // by the user: this node never sees that unit again, so its files have no purpose. The
        // caller says which of the two occurred.
        let user_stopped = cancel.is_cancelled() && !lease_lost.load(std::sync::atomic::Ordering::Relaxed);
        if !user_stopped {
            let _ = tokio::fs::remove_dir_all(&dir).await;
        }
        outcome
    }

    /// Remove the subject that a unit made, with each row and each cached result below it.
    ///
    /// **A Grid unit must leave no subject in the workspace.** The subject exists only because the
    /// analysis works on a subject. Its alignment names a file in the temporary directory of the
    /// unit, and that directory goes away at the end of the unit. A subject that stayed would name
    /// a file that is not there.
    ///
    /// Without this step, a node that contributes for one week puts some thousands of such
    /// subjects among the true subjects of its owner. Each one holds no data that a person can use,
    /// and each one is difficult to tell from a real subject. The result of the unit is already
    /// safe: the digest went to the AppView, and the records went to the publish queue.
    ///
    /// The delete of a subject refuses while the subject holds data, so this removes each sequence
    /// run first. A failure gives no message to the user, because the unit is already complete. A
    /// subject that stays is a fault for a later version to correct. It is not a reason to report
    /// a unit as failed.
    async fn discard_unit_subject(&self, guid: SampleGuid) {
        if let Ok(runs) = self.list_sequence_runs(guid).await {
            for run in runs {
                let _ = self.delete_sequence_run(run.id).await;
            }
        }
        let _ = self.delete_biosample(guid).await;
    }

    /// Tell the AppView that this node is alive, until the unit ends.
    ///
    /// This future never finishes on its own. It ends when the work beside it finishes, and
    /// `select!` then drops it. It returns only when the node **loses** the lease, which is a
    /// reason to stop the work at once.
    ///
    /// A node that lost its lease receives no credit for more work on that unit. Without this
    /// check, such a node can spend hours on a unit that another node already finished. The value
    /// that the AppView sends back is the only way for the node to learn that.
    async fn beat_while_working(
        &self,
        lease_id: i64,
        stage: &Arc<Mutex<GridStage>>,
        cancel: &CancelToken,
        lease_lost: &Arc<std::sync::atomic::AtomicBool>,
    ) -> AppError {
        loop {
            tokio::time::sleep(HEARTBEAT_EVERY).await;
            if cancel.is_cancelled() {
                // The work stops by itself. This future must not end the unit with an error that
                // hides the true reason.
                continue;
            }
            let now = stage.lock().map(|s| *s).unwrap_or(GridStage::Analyze);
            match self.grid_heartbeat(lease_id, now.as_str(), None).await {
                Ok(true) => {}
                Ok(false) => {
                    lease_lost.store(true, std::sync::atomic::Ordering::Relaxed);
                    cancel.cancel();
                    return AppError::Import("another node now holds this unit".into());
                }
                // A beat that did not arrive is not proof that the lease is gone. The network of a
                // volunteer is not always available, and the work continues. The lease has its own
                // time limit, and the AppView reclaims it if this node truly stopped.
                Err(_) => {}
            }
        }
    }

    async fn grid_unit_inner(
        &self,
        unit: &ClaimedUnit,
        params: &GridJobParams,
        dir: &Path,
        cancel: &CancelToken,
        report: &mut (dyn FnMut(GridStage, &str) + Send),
        subject: &Arc<Mutex<Option<SampleGuid>>>,
    ) -> Result<i64, AppError> {
        // A sample with more than one run needs each run mapped and then all of them merged into
        // one alignment. There is no merge stage here yet.
        //
        // The check occurs **before** the fetch, on purpose. The manifest arrives with the claim,
        // so the node knows the count at no cost.
        //
        // An earlier version took the first file and ignored the others. It pulled every one of
        // them first, from a public archive that gives us its bandwidth at no charge. It then gave
        // a coverage value from one part of the sample, as a value for the whole sample.
        if runs_in(&unit.manifest) > 1 {
            return Err(AppError::Import(format!(
                "{} holds more than one sequencing run, and this node can not merge runs yet",
                unit.sample_accession
            )));
        }

        // ---- fetch ----
        report(GridStage::Fetch, &unit.sample_accession);
        let client = self.auth.http.clone();
        let mut on_bytes = |name: &str, recv: u64, total: Option<u64>| {
            let pct = total.filter(|t| *t > 0).map(|t| recv * 100 / t).unwrap_or(0);
            report(GridStage::Fetch, &format!("{name} {pct}%"));
        };
        let files = ena::fetch_unit(&client, dir, &unit.manifest, cancel, &mut on_bytes).await?;
        let primary = primary_file(&unit.manifest, &files)
            .ok_or_else(|| AppError::Import(format!("unit {} has no data file", unit.sample_accession)))?;

        // A FASTQ unit needs a map stage that does not exist yet. The node must never reach this
        // point, because `supported_data_kinds` does not advertise FASTQ. The check stays, because
        // a wrong advertisement must give a clear message and not a strange failure much later.
        let aligned = if unit.data_kind == "FASTQ" {
            map_unit_reads(self, &files, dir, &params.reference_build, cancel, report).await?
        } else {
            primary.clone()
        };

        // ---- import ----
        report(GridStage::Import, &unit.sample_accession);
        // The ENA accession is the identity of the subject. It is a public catalogue id and not
        // personal data, so it is safe as the name that a user sees.
        let biosample = self
            .add_biosample(None, &unit.sample_accession, Some(unit.sample_accession.clone()), None)
            .await?;
        // Record the subject at once, and before the import. A failure in any step after this
        // point must still remove it.
        if let Ok(mut g) = subject.lock() {
            *g = Some(biosample.guid);
        }
        // The ENA accession as an external id. This is what makes the published record *about*
        // the public sample.
        //
        // `BiosampleRecord` carries no accession field. That rule keeps personal data out of a
        // published record. The record reads `external_ids` from this table instead.
        //
        // Without this call, the record goes out with an empty `externalIds`. No reader can then
        // connect it to the sample, or join it with the record of a second contributor.
        self.add_external_id(biosample.guid, "ENA", &unit.sample_accession)
            .await?;
        self.add_data(biosample.guid, &aligned).await?;

        // Make the coordinate index when the alignment has none.
        //
        // The fetch step takes the index of ENA when ENA has one, and that path costs least. ENA
        // does not always publish one, and that fetch can fail with no result.
        //
        // Without an index, each step that asks for a region fails. `analyze_biosample` puts those
        // failures in `errors`, and the check below then fails the unit. That occurs **after** the
        // walk over the whole file already succeeded. So the index comes first.
        report(GridStage::Import, "index");
        for aln in self.list_alignments_for_biosample(biosample.guid).await? {
            if let Err(e) = self.ensure_alignment_index(aln.id, |_, _| {}).await {
                return Err(AppError::Import(format!(
                    "{}: no coordinate index, and this node could not make one: {e}",
                    unit.sample_accession
                )));
            }
        }

        // ---- analyze ----
        report(GridStage::Analyze, &unit.sample_accession);
        let analyzed = self.analyze_biosample(&biosample, cancel.clone()).await?;
        if !analyzed.had_alignment {
            return Err(AppError::Import(format!(
                "no alignment for {} after import",
                unit.sample_accession
            )));
        }
        // `analyze_biosample` puts the failure of one step in `errors` and still gives `Ok`. That
        // is correct for a batch over the subjects of a user, where the other steps still give a
        // result that a person can use. It is **not** correct here.
        //
        // A step that failed leaves its value out of the digest. The agreement test compares an
        // absent value with an absent value as equal. So two nodes that both failed would agree,
        // reach a quorum on a result with no content, and receive credit for it. A unit must fail
        // instead, and another node then does the work.
        if !analyzed.errors.is_empty() {
            return Err(AppError::Import(format!(
                "analysis of {} did not complete: {}",
                unit.sample_accession,
                analyzed.errors.join("; ")
            )));
        }

        let mut results = UnitResults::default();
        let alignments = self.list_alignments_for_biosample(biosample.guid).await?;

        // The build that the calls are truly against.
        //
        // A `CRAM` unit is a passthrough. The submitter of that file chose its build, and that
        // build is GRCh37 or GRCh38 for most of the archive. Nothing here maps it again. The header
        // probe reads the true build during the import, and the row keeps it.
        //
        // An earlier version reported `params.reference_build` for each unit. That value is the
        // default of the command line. So a result on GRCh38 went to the AppView as a result on
        // CHM13.
        //
        // The AppView compares two results only when the build agrees. So such a result
        // joined a group of true CHM13 results. It then compared the coverage and the Y value of
        // two different references as one measurement.
        let reference_build = alignments
            .first()
            .map(|a| a.reference_build.clone())
            .unwrap_or_else(|| params.reference_build.clone());

        if let Some(aln) = alignments.first() {
            if let Some(cov) = self.cached_coverage(aln.id).await? {
                results.coverage_mean = Some(cov.mean_coverage);
                if cov.genome_territory > 0 {
                    results.callable_fraction = Some(cov.callable_bases as f64 / cov.genome_territory as f64);
                }
            }
            if let Some(sex) = self.cached_sex(aln.id).await? {
                results.sex = sex_for_digest(sex.inferred_sex);
            }
        }
        let y_calls = self.haplogroup_calls(biosample.guid, DnaType::Y).await?;
        results.y_terminal = y_calls.first().map(|c| c.haplogroup.clone());

        // ---- ancestry ----
        //
        // The autosomal consensus of a new sample can be absent, and the estimate then fails. That
        // is not a failure of the unit. The digest holds no ancestry value, and two results with no
        // ancestry value still agree. A unit that failed here would waste the hours of analysis
        // that are already complete.
        report(GridStage::Ancestry, &unit.sample_accession);
        // An error here gives no message to the user. An absent estimate is a normal result, and
        // the digest then holds no ancestry value.
        if let Ok(a) = self.estimate_ancestry_from_consensus(biosample.guid).await {
            results.ancestry_superpop_argmax = a
                .super_population_summary
                .iter()
                .max_by(|x, y| x.percentage.total_cmp(&y.percentage))
                .map(|s| s.super_population.clone());
        }

        let aligner = (unit.data_kind == "FASTQ").then_some("minimap2-pure-rs");

        // ---- publish ----
        //
        // The records go to the repository of the contributor, and each one carries the provenance
        // block. That block is what makes a record *about* a public sample that nobody owns while
        // it is *made by* this node. See design §5.1.
        //
        // A failure here does not fail the unit. The analysis is complete and its digest is the
        // thing that the quorum reads. The records are the full result behind that digest, and the
        // outbox sends them again later. A unit that failed here would discard hours of work
        // because a network call did not answer.
        report(GridStage::Publish, &unit.sample_accession);
        let record_refs = self
            .publish_grid_records(&biosample, &results, &reference_build, aligner)
            .await
            .unwrap_or_default();

        // ---- submit ----
        report(GridStage::Submit, &unit.sample_accession);
        let stack_version = env!("CARGO_PKG_VERSION");
        let digest = build_digest(
            &unit.sample_accession,
            &reference_build,
            stack_version,
            aligner,
            &results,
        );
        self.grid_submit(
            unit.work_unit_id,
            Some(unit.lease_id),
            &digest,
            stack_version,
            &reference_build,
            aligner,
            &record_refs,
        )
        .await
    }

    /// Put the records of a finished unit in the publish queue, and give back the `at://` address
    /// of each one.
    ///
    /// The queue is the durable path that the rest of the application uses. It repeats a call that
    /// failed. A second publish of the same record replaces the first record, and adds no second
    /// record. A volunteer machine goes offline, and a direct write would then lose records that a
    /// queue keeps.
    ///
    /// Each record here uses a **fixed** record key. That key gives the address of the record
    /// before the write occurs. So this method can give those addresses to the submit call in the
    /// same run, and it does not wait for the queue to empty.
    async fn publish_grid_records(
        &self,
        biosample: &Biosample,
        results: &UnitResults,
        reference_build: &str,
        aligner: Option<&str>,
    ) -> Result<Vec<String>, AppError> {
        let did = self.require_account()?;
        let prov = grid_provenance(&did, reference_build, aligner);
        let mut refs = Vec::new();

        // The biosample record is the anchor. It carries the ENA accession as an external id.
        // That id makes the record about the public sample, and not about this contributor.
        let anchor = attach_provenance(self.biosample_record(&did, biosample.guid).await?, &prov)?;
        self.enqueue_publish(
            "biosample",
            &format!("biosample:{}", biosample.guid),
            NS_BIOSAMPLE,
            Some(&biosample_rkey(biosample.guid)),
            anchor,
        )
        .await?;
        refs.push(biosample_at_uri(&did, biosample.guid));

        // The coverage record holds the measurements behind the digest. A digest says that two
        // nodes agree; this record says what they agree about.
        //
        // A failure on one record must not discard the records that already went in the queue.
        //
        // `coverage_record` gives an error in two cases. The first is an alignment with no cached
        // coverage. The second is a file that names a whole genome while its reads cover chrY only.
        // Both occur on real ENA samples.
        //
        // An earlier version used `?` here. One such error then gave an empty list, and the
        // submission named **no** record at all. It did not even name the anchor, which was
        // already in the queue and which the AppView was going to publish.
        for aln in self.list_alignments_for_biosample(biosample.guid).await? {
            if results.coverage_mean.is_none() {
                break;
            }
            let built = match self.coverage_record(&did, aln.id).await {
                Ok(v) => attach_provenance(v, &prov),
                Err(e) => Err(e),
            };
            let Ok(value) = built else { continue };
            if self
                .enqueue_publish(
                    "coverage",
                    &format!("alignment:{}", aln.id),
                    NS_ALIGNMENT,
                    Some(&alignment_rkey(aln.id)),
                    value,
                )
                .await
                .is_ok()
            {
                refs.push(format!("at://{did}/{NS_ALIGNMENT}/{}", alignment_rkey(aln.id)));
            }
        }
        Ok(refs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn full() -> UnitResults {
        UnitResults {
            sex: Some("Male".into()),
            y_terminal: Some("R-FGC29071".into()),
            ancestry_superpop_argmax: Some("EUR".into()),
            coverage_mean: Some(30.4),
            callable_fraction: Some(0.9412),
        }
    }

    #[test]
    fn the_digest_holds_the_raw_values() {
        let d = build_digest("SAMEA1", "chm13v2.0", "1.7.0", None, &full());
        assert_eq!(d["unit"], "SAMEA1");
        assert_eq!(d["reference_build"], "chm13v2.0");
        assert_eq!(d["calls"]["y_terminal"], "R-FGC29071");
        // Raw, and not in a group. The AppView makes the groups, so the rule has one home.
        assert_eq!(d["calls"]["coverage_mean"], 30.4);
        assert_eq!(d["calls"]["callable_fraction"], 0.9412);
    }

    /// The analysis path does not give an mtDNA value on CHM13, so the digest must not ask for one.
    #[test]
    fn the_digest_has_no_mtdna_field() {
        let d = build_digest("SAMEA1", "chm13v2.0", "1.7.0", None, &full());
        assert!(d["calls"].get("mt_terminal").is_none());
        assert!(!d.to_string().contains("mt_terminal"));
    }

    /// An absent value is absent from the digest. It is not `null`. Two results that both have no
    /// Y value then agree, and a result with a Y value does not agree with one that has none.
    #[test]
    fn an_absent_value_is_left_out_and_not_sent_as_null() {
        let d = build_digest("SAMEA1", "chm13v2.0", "1.7.0", None, &UnitResults::default());
        assert!(d["calls"].as_object().unwrap().is_empty());
        assert!(!d.to_string().contains("null") || d["aligner"].is_null());
    }

    /// A unit with no new alignment names no mapper. That value records how the node made the
    /// result, and a later check reads it.
    #[test]
    fn a_passthrough_unit_names_no_mapper() {
        let d = build_digest("SAMEA1", "chm13v2.0", "1.7.0", None, &full());
        assert!(d["aligner"].is_null());
        let d = build_digest("SAMEA1", "chm13v2.0", "1.7.0", Some("minimap2-pure-rs"), &full());
        assert_eq!(d["aligner"], "minimap2-pure-rs");
    }

    /// This node advertises each kind that it can process, and no other. The AppView selects work
    /// with this list, so a wrong entry here gives a node work that it can not do.
    #[test]
    fn the_node_advertises_each_kind_that_it_can_do() {
        let kinds = supported_data_kinds();
        assert!(kinds.contains(&"CRAM".to_string()), "a unit with an alignment");
        assert!(kinds.contains(&"FASTQ".to_string()), "a unit with reads only");
    }

    /// ENA names the two mate files with `_1` and `_2` before the extension.
    #[test]
    fn the_two_mate_files_are_found_by_name() {
        let f = |n: &str| PathBuf::from(format!("/x/{n}"));
        let (r1, r2, singles) = split_mates(&[f("ERR1_1.fastq.gz"), f("ERR1_2.fastq.gz")]);
        assert_eq!(r1, Some(f("ERR1_1.fastq.gz")));
        assert_eq!(r2, Some(f("ERR1_2.fastq.gz")));
        assert!(singles.is_empty());
    }

    /// A run with one file has reads with no mate, and a long-read run is always such a run.
    #[test]
    fn one_file_gives_reads_with_no_mate() {
        let f = PathBuf::from("/x/ERR1.fastq.gz");
        let (r1, r2, singles) = split_mates(std::slice::from_ref(&f));
        assert!(r1.is_none());
        assert!(r2.is_none());
        assert_eq!(singles, vec![f]);
    }

    /// The match is on `_1.` and not on `_1`. A run accession can hold those two characters, and a
    /// file that matched on the accession would go to the wrong mate.
    #[test]
    fn the_mate_match_does_not_read_the_accession() {
        let f = PathBuf::from("/x/ERR1_1_1.fastq.gz");
        let (r1, _, singles) = split_mates(std::slice::from_ref(&f));
        assert_eq!(r1, Some(f), "the mate marker is the one before the extension");
        assert!(singles.is_empty());
    }

    /// `PROVENANCE_KEY` must be the serde name of the field on the record types. This test makes a
    /// record with the typed method and then reads the key back, so it checks the string against
    /// the type. A new name in `du-domain` then fails here. Without this test, it would give a
    /// record that the AppView reads and does not understand, with no message.
    #[test]
    fn the_provenance_key_matches_the_shared_type() {
        let rec = du_domain::fed::BiosampleRecord::new(None, None, None, None, "2026-08-25T00:00:00Z")
            .with_provenance(Some(grid_provenance("did:plc:x", "chm13v2.0", None)));
        let value = serde_json::to_value(&rec).expect("serialize");
        assert!(
            value.get(PROVENANCE_KEY).is_some(),
            "the typed record wrote its provenance under a different key: {value}"
        );
    }

    /// The block that this module adds must equal the block that the typed method writes. If the
    /// two differ, a Grid record and an ordinary record carry different shapes for one idea.
    #[test]
    fn the_added_block_equals_the_block_that_the_type_writes() {
        let prov = grid_provenance("did:plc:x", "chm13v2.0", Some("minimap2-pure-rs"));
        let typed = du_domain::fed::BiosampleRecord::new(None, None, None, None, "2026-08-25T00:00:00Z")
            .with_provenance(Some(prov.clone()));
        let from_type = serde_json::to_value(&typed).unwrap()[PROVENANCE_KEY].clone();

        let plain = serde_json::to_value(du_domain::fed::BiosampleRecord::new(
            None,
            None,
            None,
            None,
            "2026-08-25T00:00:00Z",
        ))
        .unwrap();
        let added = attach_provenance(plain, &prov).unwrap()[PROVENANCE_KEY].clone();

        assert_eq!(added, from_type);
    }

    /// A passthrough unit names no mapper, and the block then has no `aligner` key at all. That is
    /// a fact about how the node made the result, and not a value that is missing.
    #[test]
    fn provenance_from_a_passthrough_unit_names_no_mapper() {
        let p = grid_provenance("did:plc:x", "chm13v2.0", None);
        assert!(p.aligner.is_none());
        let v = serde_json::to_value(&p).unwrap();
        assert!(v.get("aligner").is_none(), "an absent mapper is left out: {v}");
        assert_eq!(v["computedBy"], "did:plc:x");
        assert_eq!(v["source"], "ena:read_run");
    }

    /// A record that is not an object can not take a provenance block. That is a fault in the
    /// builder, and it must give an error and not a record with no provenance.
    #[test]
    fn a_record_that_is_not_an_object_is_refused() {
        let p = grid_provenance("did:plc:x", "chm13v2.0", None);
        assert!(attach_provenance(serde_json::json!("not a record"), &p).is_err());
    }

    /// A run with two mates and a file of reads that lost their mate is **one** run. ENA gives
    /// three files for such a run, and a count of files would refuse it.
    #[test]
    fn three_files_of_one_run_are_one_run() {
        let f = |run: &str, name: &str| ManifestFile {
            run_accession: run.into(),
            url: format!("ftp/{name}"),
            index_url: None,
            md5: None,
            bytes: None,
            format: "FASTQ".into(),
        };
        let one_run = vec![
            f("ERR1", "ERR1_1.fastq.gz"),
            f("ERR1", "ERR1_2.fastq.gz"),
            f("ERR1", "ERR1.fastq.gz"),
        ];
        assert_eq!(runs_in(&one_run), 1);

        let two_runs = vec![f("ERR1", "ERR1_1.fastq.gz"), f("ERR2", "ERR2_1.fastq.gz")];
        assert_eq!(runs_in(&two_runs), 2);
    }

    /// An empty accession must not refuse the unit. The safe reading is one run.
    #[test]
    fn a_manifest_with_no_accession_counts_as_one_run() {
        let f = ManifestFile {
            run_accession: String::new(),
            url: "ftp/x.cram".into(),
            index_url: None,
            md5: None,
            bytes: None,
            format: "CRAM".into(),
        };
        assert_eq!(runs_in(std::slice::from_ref(&f)), 1);
        assert_eq!(runs_in(&[]), 1);
    }

    /// The release reason that goes to the server must carry no local path. The message of an
    /// error can hold one, because the fetch module puts the path of a file in its error text.
    #[test]
    fn the_release_reason_carries_no_local_path() {
        let leaky = AppError::Import(
            "/Users/someone/Library/navigator-grid/SAMEA1/x.cram.part: No space left on device".into(),
        );
        let reason = release_reason(&leaky);
        assert_eq!(reason, "disk");
        assert!(!reason.contains('/'), "a path must never reach the server");
        assert!(!reason.contains("Users"));
    }

    /// Each class is short, and each one tells the operator of the AppView something different.
    #[test]
    fn each_failure_gives_its_own_short_class() {
        let r = |m: &str| release_reason(&AppError::Import(m.into()));
        assert_eq!(r("cancelled"), "stopped");
        assert_eq!(r("checksum mismatch for x.cram"), "checksum");
        assert_eq!(r("not enough room for this work unit"), "disk");
        assert_eq!(r("this node can not merge runs yet"), "unsupported");
        assert_eq!(r("analysis of SAMEA1 did not complete"), "analysis");
        assert_eq!(r("something else"), "error");
    }

    /// An index file is never the primary file of a unit.
    #[test]
    fn the_index_file_is_not_the_primary_file() {
        let m = |fmt: &str, url: &str| ManifestFile {
            run_accession: "ERR1".into(),
            url: url.into(),
            index_url: None,
            md5: None,
            bytes: None,
            format: fmt.into(),
        };
        let manifest = vec![m("CRAI", "a.cram.crai"), m("CRAM", "a.cram")];
        let files = vec![PathBuf::from("/x/a.cram.crai"), PathBuf::from("/x/a.cram")];
        assert_eq!(primary_file(&manifest, &files).unwrap(), &PathBuf::from("/x/a.cram"));
    }
}

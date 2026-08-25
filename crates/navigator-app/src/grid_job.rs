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

/// The primary data file of a unit: the alignment for a CRAM unit, or the first read file for a
/// FASTQ unit. An index file is never the primary file.
fn primary_file<'a>(manifest: &'a [ManifestFile], files: &'a [PathBuf]) -> Option<&'a PathBuf> {
    manifest
        .iter()
        .zip(files)
        .find(|(m, _)| matches!(m.format.as_str(), "CRAM" | "BAM" | "FASTQ"))
        .map(|(_, p)| p)
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
        let dir = params.scratch_root.join(&unit.sample_accession);

        match self.grid_unit_inner(unit, params, &dir, cancel, report).await {
            Ok(id) => outcome.submission_id = Some(id),
            Err(e) => {
                outcome.error = Some(e.to_string());
                // The unit goes back to the catalogue at once. Without this call, it waits for the
                // full lease time, and no other node can take it.
                let _ = self.grid_release(unit.lease_id, &e.to_string()).await;
            }
        }
        // The files of a unit are large. Remove them whatever the result, or a node that runs for a
        // week fills the disk of its owner.
        let _ = tokio::fs::remove_dir_all(&dir).await;
        outcome
    }

    async fn grid_unit_inner(
        &self,
        unit: &ClaimedUnit,
        params: &GridJobParams,
        dir: &Path,
        cancel: &CancelToken,
        report: &mut (dyn FnMut(GridStage, &str) + Send),
    ) -> Result<i64, AppError> {
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
        self.add_data(biosample.guid, &aligned).await?;

        // ---- analyze ----
        report(GridStage::Analyze, &unit.sample_accession);
        let analyzed = self.analyze_biosample(&biosample, cancel.clone()).await?;
        if !analyzed.had_alignment {
            return Err(AppError::Import(format!(
                "no alignment for {} after import",
                unit.sample_accession
            )));
        }

        let mut results = UnitResults::default();
        let alignments = self.list_alignments_for_biosample(biosample.guid).await?;
        if let Some(aln) = alignments.first() {
            if let Some(cov) = self.cached_coverage(aln.id).await? {
                results.coverage_mean = Some(cov.mean_coverage);
                if cov.genome_territory > 0 {
                    results.callable_fraction = Some(cov.callable_bases as f64 / cov.genome_territory as f64);
                }
            }
            if let Some(sex) = self.cached_sex(aln.id).await? {
                results.sex = Some(format!("{:?}", sex.inferred_sex));
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
            .publish_grid_records(&biosample, &results, params, aligner)
            .await
            .unwrap_or_default();

        // ---- submit ----
        report(GridStage::Submit, &unit.sample_accession);
        let stack_version = env!("CARGO_PKG_VERSION");
        let digest = build_digest(
            &unit.sample_accession,
            &params.reference_build,
            stack_version,
            aligner,
            &results,
        );
        self.grid_submit(
            unit.work_unit_id,
            Some(unit.lease_id),
            &digest,
            stack_version,
            &params.reference_build,
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
        params: &GridJobParams,
        aligner: Option<&str>,
    ) -> Result<Vec<String>, AppError> {
        let did = self.require_account()?;
        let prov = grid_provenance(&did, &params.reference_build, aligner);
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
        for aln in self.list_alignments_for_biosample(biosample.guid).await? {
            if results.coverage_mean.is_none() {
                break;
            }
            let value = attach_provenance(self.coverage_record(&did, aln.id).await?, &prov)?;
            self.enqueue_publish(
                "coverage",
                &format!("alignment:{}", aln.id),
                NS_ALIGNMENT,
                Some(&alignment_rkey(aln.id)),
                value,
            )
            .await?;
            refs.push(format!("at://{did}/{NS_ALIGNMENT}/{}", alignment_rkey(aln.id)));
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

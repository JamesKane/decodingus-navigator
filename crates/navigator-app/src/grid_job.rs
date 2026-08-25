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
use navigator_analysis::CancelToken;
use std::path::{Path, PathBuf};

/// The stage that a node is on. The node sends this with each heartbeat, and the fleet view of the
/// AppView shows it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridStage {
    Fetch,
    Import,
    Analyze,
    Ancestry,
    Submit,
}

impl GridStage {
    pub fn as_str(self) -> &'static str {
        match self {
            GridStage::Fetch => "fetch",
            GridStage::Import => "import",
            GridStage::Analyze => "analyze",
            GridStage::Ancestry => "ancestry",
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
/// A `FASTQ` unit needs the node to map the reads first. That path is not here yet. So this list
/// does not hold `FASTQ`, and the AppView never offers such a unit to this node. The capability
/// filter is not a suggestion. It is what stops a node when it can not do the work. See
/// [`map_reads`].
pub fn supported_data_kinds() -> Vec<String> {
    vec!["CRAM".to_string()]
}

/// Map the reads of a FASTQ unit to the target build. **This is not written yet.**
///
/// The work is small but it is not zero, and it touches a module that is already in use. The
/// realignment job (`realign_job`) has the stages that a FASTQ unit needs: index, map, sort, mark
/// duplicates, and finalize. Its stage A recovers reads from an alignment and writes them as FASTQ
/// files. Its stage B then maps those files. So a FASTQ unit from ENA is the same pipeline with a
/// different source for stage A.
///
/// To make that possible, stage A of `realign_job` must accept read files from outside. That is a
/// change to a module that shipped in `v0.1.0-alpha.17` and that phase 5 validated on a full
/// genome. Such a change belongs in its own commit, where a reviewer can compare it against that
/// validated behaviour. It does not belong inside a first version of this driver.
///
/// Until then [`supported_data_kinds`] does not hold `FASTQ`, so no node claims such a unit.
fn map_reads(_reads: &[PathBuf], _target_build: &str) -> Result<PathBuf, AppError> {
    Err(AppError::Import(
        "this node can not map FASTQ reads yet; it must not have claimed a FASTQ unit".into(),
    ))
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
            map_reads(&files, &params.reference_build)?
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

        // ---- submit ----
        report(GridStage::Submit, &unit.sample_accession);
        let stack_version = env!("CARGO_PKG_VERSION");
        let aligner = (unit.data_kind == "FASTQ").then_some("minimap2-pure-rs");
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
            &[],
        )
        .await
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

    /// This node advertises only the kinds that it can process. A FASTQ unit needs a map stage that
    /// does not exist yet, so the node must not advertise FASTQ.
    #[test]
    fn the_node_advertises_only_what_it_can_do() {
        let kinds = supported_data_kinds();
        assert!(kinds.contains(&"CRAM".to_string()));
        assert!(
            !kinds.contains(&"FASTQ".to_string()),
            "the map stage is not written, so a FASTQ unit must never reach this node"
        );
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

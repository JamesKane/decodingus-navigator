//! The FTDNA project import. This module holds the engine that matches a kit to a subject and
//! finds a duplicate. It also holds the two steps of the import, which are the plan and the commit.
//! See design §5 and §6.
//!
//! Phase 1 covers the roster and the ancestry, which are the base of the import. The module does
//! four steps. It parses the batch CSV files, joins them by the kit number, matches each kit
//! against the workspace, and makes a **plan** for the administrator. This phase writes nothing.
//!
//! A separate commit step applies the plan. That step uses the decisions of the administrator for
//! each candidate that the engine is not sure about.
//!
//! A later change adds the deep data of each member, which is Big Y, mtDNA, and Family Finder. It
//! also adds the wide Y-STR chart. This module only connects the identity, the MDKA rows, and the
//! membership.

use std::collections::BTreeMap;
use std::path::PathBuf;

use chrono::Utc;
use du_domain::ids::SampleGuid;
use navigator_domain::ftdna::{self, AncestryRow, MemberRow};
use navigator_domain::identity::{FtdnaMember, IdSource, Lineage, NewMdka};
use navigator_domain::reconciliation::DnaType;
use navigator_domain::strprofile::{NewStrProfile, StrMarker};
use navigator_store::{biosample, biosample_project, external_id, ftdna_member, mdka, str_profile, sync_state};

use crate::{decodingus_appview_url, AccessionBackfill, App, AppError, CatalogBackfill};

/// The subset of the AppView `/api/v1/samples/{alias}` response the accession backfill reads. The
/// full record carries haplogroups/publications/etc.; we only need the authoritative accession.
#[derive(serde::Deserialize)]
struct CatalogSample {
    #[serde(default)]
    accession: Option<String>,
}

/// The values that control the engine that matches a kit to a subject.
#[derive(Debug, Clone)]
pub struct FtdnaImportOptions {
    /// The lowest score, from 0 to 1, that lets the engine offer a subject as a merge candidate.
    /// This score is not exact, and the engine calculates it from the name.
    pub fuzzy_threshold: f32,
}

impl Default for FtdnaImportOptions {
    fn default() -> Self {
        // Conservative: a Y-terminal match alone qualifies; a weak name-only hint does not.
        Self { fuzzy_threshold: 0.5 }
    }
}

/// The parsed + cross-file-joined data for one kit (the payload a plan row commits).
#[derive(Debug, Clone)]
pub struct FtdnaSubjectInput {
    pub kit_number: String,
    pub member: Option<MemberRow>,
    pub paternal: Option<AncestryRow>,
    pub maternal: Option<AncestryRow>,
    /// Y-STR markers from the wide overview (empty if no Y-STR file / no row for this kit).
    pub ystr_markers: Vec<StrMarker>,
}

/// A workspace Subject offered as a fuzzy merge candidate, with why.
#[derive(Debug, Clone)]
pub struct FuzzyCandidate {
    pub guid: SampleGuid,
    pub donor_identifier: String,
    pub score: f32,
    pub reasons: Vec<String>,
}

/// The action that the engine proposes for a kit. The engine merges without a question only for an
/// exact match on the vendor id. For a match that is not exact, it adds the kit to a list, and the
/// administrator decides. The engine never merges such a kit on its own.
#[derive(Debug, Clone)]
pub enum MatchKind {
    /// No workspace match → create a new Subject.
    New,
    /// Exact `external_id(FTDNA, kit)` hit → reuse that Subject (design §5.1, locked).
    AutoMerge { guid: SampleGuid, donor_identifier: String },
    /// Fuzzy candidates above threshold → the admin confirms/rejects each (design §5.2).
    NeedsConfirm { candidates: Vec<FuzzyCandidate> },
}

/// One row of the dry-run plan.
#[derive(Debug, Clone)]
pub struct FtdnaPlanRow {
    pub kit_number: String,
    /// Best display label (kit + ancestor/member name).
    pub label: String,
    /// FTDNA-reported Y terminal SNP from the paternal clade (provisional label until the YDNA
    /// overview supplies the full `R-…` haplogroup).
    pub y_terminal: Option<String>,
    /// `false` = ancestry data for a kit absent from the roster (orphan; still importable, flagged).
    pub in_roster: bool,
    /// Number of Y-STR markers that will attach (from the wide overview).
    pub ystr_count: usize,
    pub kind: MatchKind,
    pub input: FtdnaSubjectInput,
}

/// The counts of the input files that the code recognized, and of the rows that it read. The review
/// header shows these counts. So the administrator sees an absent file, or a file with the wrong
/// class, at once. One example is an import with no roster. Without these counts, such an import
/// gives rows with no subject and no message.
#[derive(Debug, Clone, Default)]
pub struct FtdnaPlanStats {
    /// Roster rows parsed from `Member_Information`.
    pub roster: usize,
    /// Rows parsed from the paternal ancestry file.
    pub paternal: usize,
    /// Rows parsed from the maternal ancestry file.
    pub maternal: usize,
    /// Kits with Y-STR markers from the wide overview.
    pub ystr: usize,
    /// Workspace Subjects scanned for matches.
    pub scanned_subjects: usize,
}

/// The reviewable plan: every kit with its proposed disposition. No writes happen until commit.
#[derive(Debug, Clone)]
pub struct FtdnaImportPlan {
    /// Target project, or `None` to create one named [`Self::project_name`] at commit (so a cancelled
    /// dry-run leaves no empty project behind).
    pub project_id: Option<i64>,
    /// The target/derived project name (shown in the review header).
    pub project_name: String,
    /// Recognized-input counts (header diagnostics).
    pub stats: FtdnaPlanStats,
    pub rows: Vec<FtdnaPlanRow>,
}

impl FtdnaImportPlan {
    /// `(new, auto_merge, needs_confirm)` counts for the review header.
    pub fn counts(&self) -> (usize, usize, usize) {
        let mut c = (0, 0, 0);
        for r in &self.rows {
            match r.kind {
                MatchKind::New => c.0 += 1,
                MatchKind::AutoMerge { .. } => c.1 += 1,
                MatchKind::NeedsConfirm { .. } => c.2 += 1,
            }
        }
        c
    }
}

/// The admin's decision for a fuzzy row, keyed by kit number at commit.
#[derive(Debug, Clone)]
pub enum FtdnaResolution {
    /// Merge this kit into an existing Subject.
    Merge(SampleGuid),
    /// Treat as a new Subject.
    New,
    /// Do not import this kit at all.
    Skip,
}

/// What the commit did.
#[derive(Debug, Clone, Default)]
pub struct FtdnaImportSummary {
    /// The project that received the kits. The commit step finds this project or makes it.
    pub project_id: i64,
    pub created: usize,
    pub merged: usize,
    pub memberships_added: usize,
    pub mdka_written: usize,
    /// Kits that had ancestry data but no roster row.
    pub orphans: usize,
    /// Kits the admin chose to skip.
    pub skipped: usize,
    /// Y-STR profiles attached (from the wide overview).
    pub str_profiles: usize,
    pub errors: Vec<String>,
}

/// The genealogy data that the app imported for one subject. It holds the vendor ids, the FTDNA
/// member labels, and the MDKA rows.
///
/// This data is personal. The app shows it on this machine only, and it never sends it to the
/// network. The value is empty when the app imported nothing for the subject.
#[derive(Debug, Clone, Default)]
pub struct FtdnaGenealogy {
    pub external_ids: Vec<navigator_domain::identity::ExternalId>,
    pub member: Option<navigator_domain::identity::FtdnaMember>,
    pub mdka: Vec<navigator_domain::identity::Mdka>,
}

impl FtdnaGenealogy {
    /// Shows that the app imported nothing. The UI can then skip the detail card.
    pub fn is_empty(&self) -> bool {
        self.external_ids.is_empty() && self.member.is_none() && self.mdka.is_empty()
    }
}

impl App {
    /// One-shot read of a Subject's imported genealogy (vendor ids + FTDNA member + MDKA) for the
    /// subject-detail card.
    pub async fn subject_genealogy(&self, guid: SampleGuid) -> Result<FtdnaGenealogy, AppError> {
        Ok(FtdnaGenealogy {
            external_ids: self.external_ids(guid).await?,
            member: self.ftdna_member(guid).await?,
            mdka: self.mdka_for(guid).await?,
        })
    }

    /// Add a vendor id, which is a kit number, to a subject from the subject editor.
    ///
    /// The method refuses a blank source and a blank id. It also refuses a `(source, external_id)`
    /// pair that belongs to a *different* subject.
    ///
    /// That pair is unique, and the app uses it to find a duplicate donor. The method must never
    /// move the pair to another subject without a message. The caller resolves such a conflict.
    ///
    /// A second call for the same subject is safe.
    pub async fn add_external_id(
        &self,
        guid: SampleGuid,
        source: &str,
        external_id: &str,
    ) -> Result<navigator_domain::identity::ExternalId, AppError> {
        let (source, external_id) = (source.trim(), external_id.trim());
        if source.is_empty() || external_id.is_empty() {
            return Err(AppError::Import("vendor source and id are both required".into()));
        }
        let row = external_id::add(self.store.pool(), guid, source, external_id).await?;
        if row.biosample_guid != guid {
            return Err(AppError::Import(format!(
                "{source} id \"{external_id}\" is already linked to another subject"
            )));
        }
        // The published dedup anchor changed → refresh the biosample record (best-effort).
        let _ = self.republish_biosample_ids(guid).await;
        Ok(row)
    }

    /// Detach a vendor id (by row id) from a Subject.
    pub async fn delete_external_id(&self, id: i64) -> Result<(), AppError> {
        // Read the subject of this row before the code deletes the row. The app then refreshes
        // the published record of that subject.
        let guid = external_id::get(self.store.pool(), id).await?.map(|e| e.biosample_guid);
        external_id::delete(self.store.pool(), id).await?;
        if let Some(guid) = guid {
            let _ = self.republish_biosample_ids(guid).await;
        }
        Ok(())
    }

    /// Add the public-catalog external ids that the code can derive from the local provenance of
    /// each subject. The namespaces are `IGSR`, `HGDP`, and INSDC, and
    /// [`navigator_domain::identity::catalog_ids_from_provenance`] derives them.
    ///
    /// A public dataset that a user imported in bulk then publishes ids that match its rows in the
    /// catalog of the AppView.
    ///
    /// The method is deterministic and uses no network. A sample with only a friendly name gives no
    /// id. A second call is safe, because the method skips an id that already exists.
    ///
    /// The method counts a `(namespace, value)` pair that belongs to a *different* subject as a
    /// conflict, and it changes nothing. It never moves such a pair without a message.
    /// `apply == false` makes the method report the changes and write nothing.
    ///
    /// The method writes to the store directly, and it does not publish a record for each id.
    /// Publish the subjects that changed after the method completes.
    pub async fn backfill_catalog_ids(
        &self,
        project_id: Option<i64>,
        apply: bool,
    ) -> Result<CatalogBackfill, AppError> {
        let mut out = CatalogBackfill {
            applied: apply,
            ..CatalogBackfill::default()
        };
        for b in self.list_all_biosamples().await? {
            if let Some(pid) = project_id {
                if b.project_id != Some(pid) {
                    continue;
                }
            }
            out.subjects_examined += 1;
            let derived = navigator_domain::identity::catalog_ids_from_provenance(
                &b.donor_identifier,
                b.sample_accession.as_deref(),
            );
            if derived.is_empty() {
                continue;
            }
            out.subjects_matched += 1;
            let existing: std::collections::HashSet<(String, String)> =
                external_id::list_for(self.store.pool(), b.guid)
                    .await?
                    .into_iter()
                    .map(|e| (e.source, e.external_id))
                    .collect();
            for (ns, val) in derived {
                if existing.contains(&(ns.clone(), val.clone())) {
                    continue;
                }
                out.ids_to_add += 1;
                if apply {
                    let row = external_id::add(self.store.pool(), b.guid, &ns, &val).await?;
                    if row.biosample_guid == b.guid {
                        out.ids_added += 1;
                    } else {
                        // This (namespace, value) pair belongs to another subject. The user
                        // imported the same data twice. Change nothing.
                        out.conflicts += 1;
                    }
                }
            }
        }
        Ok(out)
    }

    /// Read one public-catalog sample record from the samples API of the AppView. The path is
    /// `/api/v1/samples/{alias}`, the read is public, and the alias is our `donor_identifier`.
    ///
    /// A 404 response gives `Ok(None)`, which means that the catalog does not know the alias. That
    /// result is normal while a correction on the server is not complete.
    ///
    /// The `accession` value in the response has authority. Our local `sample_accession` field does
    /// not hold it.
    async fn fetch_catalog_sample(&self, base: &str, alias: &str) -> Result<Option<CatalogSample>, AppError> {
        let url = format!("{}/api/v1/samples/{alias}", base.trim_end_matches('/'));
        let resp = self
            .auth
            .http
            .get(&url)
            .send()
            .await
            .map_err(|e| AppError::Import(format!("catalog API {alias}: {e}")))?;
        if resp.status() == reqwest::StatusCode::NOT_FOUND {
            return Ok(None);
        }
        if !resp.status().is_success() {
            return Err(AppError::Import(format!("catalog API {alias}: HTTP {}", resp.status())));
        }
        let s = resp
            .json::<CatalogSample>()
            .await
            .map_err(|e| AppError::Import(format!("catalog API {alias}: {e}")))?;
        Ok(Some(s))
    }

    /// Look up each subject in the samples API of the AppView and add the full set of
    /// public-catalog ids **in one pass**.
    ///
    /// The set holds two kinds of id. The first is the catalog *name* id, which is an IGSR id or an
    /// HGDP id, and the code derives it from the donor id. The second is the INSDC *accession* that
    /// the API returns, and that value has authority. A `SAMN` prefix gives BIOSAMPLE, an `ERS`
    /// prefix gives ENA, and an `SRS` prefix gives SRA. The method also corrects the local
    /// `sample_accession` field, which holds a temporary value.
    ///
    /// This method does more than [`backfill_catalog_ids`](Self::backfill_catalog_ids), which uses
    /// the name only and needs no network. Use that method when the API is not available.
    ///
    /// By default the method queries only a subject whose `donor_identifier` looks like a catalog
    /// alias, which is an IGSR alias or an HGDP alias. The `all` option removes that limit. The
    /// default stops many 404 responses for a friendly name.
    ///
    /// `apply == false` makes the method write nothing. `limit` sets the maximum count of
    /// queries.
    pub async fn backfill_accessions(
        &self,
        project_id: Option<i64>,
        apply: bool,
        all: bool,
        limit: Option<usize>,
    ) -> Result<AccessionBackfill, AppError> {
        let base = decodingus_appview_url();
        let mut out = AccessionBackfill {
            applied: apply,
            ..AccessionBackfill::default()
        };
        for b in self.list_all_biosamples().await? {
            if let Some(pid) = project_id {
                if b.project_id != Some(pid) {
                    continue;
                }
            }
            // Skip samples whose name is not a recognizable catalog alias unless `--all`.
            if !all && navigator_domain::identity::catalog_ids_from_provenance(&b.donor_identifier, None).is_empty() {
                continue;
            }
            if limit.is_some_and(|n| out.examined >= n) {
                break;
            }
            out.examined += 1;
            let sample = match self.fetch_catalog_sample(&base, &b.donor_identifier).await {
                Ok(Some(s)) => s,
                Ok(None) => {
                    out.not_found += 1;
                    continue;
                }
                Err(_) => {
                    out.errors += 1;
                    continue;
                }
            };
            out.resolved += 1;
            let fetched_acc = sample.accession.as_deref().map(str::trim).filter(|a| !a.is_empty());
            // One pass gives both ids. The catalog *name* id comes from the donor id. The INSDC
            // *accession* comes from the API, when the API holds a real one. The shared helper
            // joins the two sources.
            let ids = navigator_domain::identity::catalog_ids_from_provenance(&b.donor_identifier, fetched_acc);
            if ids.is_empty() {
                continue;
            }
            // Surface the accession resolution (the API's contribution) in the examples.
            if let Some(acc) = fetched_acc {
                if navigator_domain::identity::insdc_sample_namespace(acc).is_some() && out.examples.len() < 10 {
                    out.examples.push(format!("{} → {acc}", b.donor_identifier));
                }
            }
            let existing: std::collections::HashSet<(String, String)> =
                external_id::list_for(self.store.pool(), b.guid)
                    .await?
                    .into_iter()
                    .map(|e| (e.source, e.external_id))
                    .collect();
            for (ns, val) in &ids {
                if existing.contains(&(ns.clone(), val.clone())) {
                    continue;
                }
                out.ids_to_add += 1;
                if apply {
                    let row = external_id::add(self.store.pool(), b.guid, ns, val).await?;
                    if row.biosample_guid == b.guid {
                        out.ids_added += 1;
                    } else {
                        out.conflicts += 1;
                    }
                }
            }
            // Correct the local placeholder accession to the authoritative INSDC one.
            if apply {
                if let Some(acc) = fetched_acc {
                    if navigator_domain::identity::insdc_sample_namespace(acc).is_some()
                        && b.sample_accession.as_deref() != Some(acc)
                    {
                        biosample::set_accession(self.store.pool(), b.guid, acc).await?;
                        out.accession_updated += 1;
                    }
                }
            }
        }
        Ok(out)
    }

    /// Publish the biosample anchor of a subject again, after the set of identifiers of that
    /// subject changed. The mirror of the AppView replaces the full `external_ids` field, so it
    /// then holds each addition and each removal.
    ///
    /// The method uses a fixed rkey, so the second publish replaces the record.
    ///
    /// The method acts **only for a subject that the app already published**, and only while an
    /// account is active. With no active account, or for a subject that the app never published, it
    /// does nothing. A new local id must not put a donor on the network for the first time.
    async fn republish_biosample_ids(&self, guid: SampleGuid) -> Result<(), AppError> {
        let Some(did) = self.current_account() else {
            return Ok(());
        };
        if sync_state::get(self.store.pool(), &did, &format!("biosample:{guid}"))
            .await?
            .is_none()
        {
            return Ok(());
        }
        self.publish_biosample(guid).await
    }

    /// Insert or change the MDKA of a subject for one lineage, from the subject editor. There is
    /// one row for each lineage, and the method sets `updated_at`. Give a `source` of `MANUAL` for
    /// a row that the user typed.
    pub async fn upsert_mdka(&self, guid: SampleGuid, mdka: NewMdka) -> Result<(), AppError> {
        let now = Utc::now().to_rfc3339();
        mdka::upsert(self.store.pool(), guid, &mdka, &now).await?;
        Ok(())
    }

    /// Remove a Subject's MDKA for one lineage.
    pub async fn delete_mdka(&self, guid: SampleGuid, lineage: &str) -> Result<(), AppError> {
        mdka::delete(self.store.pool(), guid, lineage).await?;
        Ok(())
    }

    /// Parse the FTDNA batch files, join by kit, and match against the workspace → a dry-run plan.
    /// Any of the three files may be absent (a roster-only or ancestry-only import is valid).
    ///
    /// `project_id` targets an existing project; pass `None` to import into a new project (created at
    /// commit, named `project_name` or a default). Matching is workspace-global, so no project need
    /// exist yet for the plan.
    #[allow(clippy::too_many_arguments)] // distinct optional file paths + target + options
    pub async fn plan_ftdna_import(
        &self,
        project_id: Option<i64>,
        project_name: Option<String>,
        member_path: Option<PathBuf>,
        paternal_path: Option<PathBuf>,
        maternal_path: Option<PathBuf>,
        ystr_path: Option<PathBuf>,
        options: FtdnaImportOptions,
    ) -> Result<FtdnaImportPlan, AppError> {
        // Resolve a display name: the existing project's name, else the caller's, else a default.
        let resolved_name = match project_id {
            Some(id) => navigator_store::project::get(self.store.pool(), id)
                .await?
                .map(|p| p.name)
                .unwrap_or_else(|| "FTDNA Project".to_string()),
            None => project_name
                .filter(|n| !n.trim().is_empty())
                .unwrap_or_else(|| "FTDNA Project".to_string()),
        };
        let members = match member_path {
            Some(p) => ftdna::parse_member_information(&std::fs::read_to_string(p)?).map_err(AppError::Import)?,
            None => Vec::new(),
        };
        let paternal = match paternal_path {
            Some(p) => ftdna::parse_ancestry(&std::fs::read_to_string(p)?).map_err(AppError::Import)?,
            None => Vec::new(),
        };
        let maternal = match maternal_path {
            Some(p) => ftdna::parse_ancestry(&std::fs::read_to_string(p)?).map_err(AppError::Import)?,
            None => Vec::new(),
        };
        let ystr = match ystr_path {
            Some(p) => ftdna::parse_ydna_overview(&std::fs::read_to_string(p)?).map_err(AppError::Import)?,
            None => Vec::new(),
        };
        let mut stats = FtdnaPlanStats {
            roster: members.len(),
            paternal: paternal.len(),
            maternal: maternal.len(),
            ystr: ystr.len(),
            scanned_subjects: 0,
        };
        // The import holds a roster only when it holds member rows. The "orphan" mark applies
        // only in that case. An orphan is data with no roster row.
        let roster_provided = !members.is_empty();

        // Join by kit number (BTreeMap → stable, kit-sorted plan).
        let mut inputs: BTreeMap<String, FtdnaSubjectInput> = BTreeMap::new();
        let mut roster: std::collections::HashSet<String> = std::collections::HashSet::new();
        for m in members {
            let kit = m.kit_number.clone();
            roster.insert(kit.clone());
            inputs.entry(kit.clone()).or_insert_with(|| empty_input(&kit)).member = Some(m);
        }
        for a in paternal {
            let kit = a.kit_number.clone();
            inputs.entry(kit.clone()).or_insert_with(|| empty_input(&kit)).paternal = Some(a);
        }
        for a in maternal {
            let kit = a.kit_number.clone();
            inputs.entry(kit.clone()).or_insert_with(|| empty_input(&kit)).maternal = Some(a);
        }
        for (kit, markers) in ystr {
            inputs
                .entry(kit.clone())
                .or_insert_with(|| empty_input(&kit))
                .ystr_markers = markers;
        }

        // Precompute each workspace Subject's Y terminal once (avoids O(kits × subjects) consensus reads).
        let existing = self.existing_subject_index().await?;
        stats.scanned_subjects = existing.len();

        let mut rows = Vec::with_capacity(inputs.len());
        for (kit, input) in inputs {
            let y_terminal = input
                .paternal
                .as_ref()
                .and_then(|a| a.sub_group.as_deref())
                .and_then(terminal_snp);
            let kind = self
                .match_kit(&kit, &input, y_terminal.as_deref(), &existing, options.fuzzy_threshold)
                .await?;
            rows.push(FtdnaPlanRow {
                label: display_label(&kit, &input),
                kit_number: kit,
                y_terminal,
                // Mark the kit as an orphan only when the import holds a roster and that roster
                // does not name the kit.
                in_roster: !roster_provided || roster.contains(&input.kit_number),
                ystr_count: input.ystr_markers.len(),
                kind,
                input,
            });
        }
        Ok(FtdnaImportPlan {
            project_id,
            project_name: resolved_name,
            stats,
            rows,
        })
    }

    /// Apply a plan. `resolutions` holds the decision of the administrator for each kit with the
    /// `NeedsConfirm` mark. A kit with that mark and no decision becomes a **New** subject. This
    /// default is the safe one, because the method must never merge a kit without a decision.
    pub async fn commit_ftdna_import(
        &self,
        plan: &FtdnaImportPlan,
        resolutions: &BTreeMap<String, FtdnaResolution>,
    ) -> Result<FtdnaImportSummary, AppError> {
        let mut summary = FtdnaImportSummary::default();
        let now = Utc::now().to_rfc3339();

        // Find the target project. Make the project now when the plan names a new one.
        let project_id = match plan.project_id {
            Some(id) => id,
            None => {
                self.create_project(navigator_domain::workspace::NewProject {
                    name: plan.project_name.clone(),
                    description: None,
                    administrator: "unknown".to_string(),
                })
                .await?
                .id
            }
        };
        summary.project_id = project_id;

        for row in &plan.rows {
            // An explicit Skip (only meaningful for a fuzzy row) drops the kit entirely.
            if matches!(resolutions.get(&row.kit_number), Some(FtdnaResolution::Skip)) {
                summary.skipped += 1;
                continue;
            }
            // Resolve the effective target: existing guid (merge) or None (create new).
            let target = match &row.kind {
                MatchKind::AutoMerge { guid, .. } => Some(*guid),
                MatchKind::New => None,
                MatchKind::NeedsConfirm { .. } => match resolutions.get(&row.kit_number) {
                    Some(FtdnaResolution::Merge(g)) => Some(*g),
                    _ => None,
                },
            };

            let result = self.commit_one(project_id, row, target, &now).await;
            match result {
                Ok((wrote_mdka, wrote_str)) => {
                    if target.is_some() {
                        summary.merged += 1;
                    } else {
                        summary.created += 1;
                    }
                    summary.memberships_added += 1;
                    summary.mdka_written += wrote_mdka;
                    summary.str_profiles += wrote_str as usize;
                    if !row.in_roster {
                        summary.orphans += 1;
                    }
                }
                Err(e) => summary.errors.push(format!("{}: {e}", row.kit_number)),
            }
        }
        Ok(summary)
    }

    /// Commit one plan row to `guid` (merge) or a fresh Subject (create). Returns
    /// `(mdka_rows_written, str_profile_created)`.
    async fn commit_one(
        &self,
        project_id: i64,
        row: &FtdnaPlanRow,
        target: Option<SampleGuid>,
        now: &str,
    ) -> Result<(usize, bool), AppError> {
        let pool = self.store.pool();
        let input = &row.input;

        // Resolve the Subject: reuse on merge, else create with the kit as the stable donor id.
        let guid = match target {
            Some(g) => g,
            None => {
                self.add_biosample(Some(project_id), input.kit_number.clone(), None, None)
                    .await?
                    .guid
            }
        };

        // The vendor identity. A second call is safe, and the code never moves an id that
        // belongs to another subject.
        external_id::add(pool, guid, IdSource::FTDNA, &input.kit_number).await?;

        // FTDNA-reported member labels.
        let member_name = input.member.as_ref().and_then(|m| clean_name(m.name.as_deref()));
        ftdna_member::upsert(
            pool,
            &FtdnaMember {
                biosample_guid: guid,
                member_name,
                y_haplogroup_ftdna: row.y_terminal.clone(),
                mt_haplogroup_ftdna: None,
                haplo_status: None,
                access_granted: input.member.as_ref().and_then(|m| m.access_granted.clone()),
                publicly_shares: input.member.as_ref().and_then(|m| m.publicly_shares),
            },
        )
        .await?;

        // The MDKA rows from the paternal (Y) ancestry and the maternal (Mt) ancestry. The code
        // writes a row only when the ancestry holds a value.
        let mut wrote = 0;
        if let Some(m) = input.paternal.as_ref().and_then(|a| mdka_from(a, Lineage::Y)) {
            mdka::upsert(pool, guid, &m, now).await?;
            wrote += 1;
        }
        if let Some(m) = input.maternal.as_ref().and_then(|a| mdka_from(a, Lineage::Mt)) {
            mdka::upsert(pool, guid, &m, now).await?;
            wrote += 1;
        }

        // Project membership (the M:N link; role = the clade subgroup label if present).
        let role = input
            .paternal
            .as_ref()
            .and_then(|a| a.sub_group.as_deref())
            .map(subgroup_role);
        biosample_project::add(pool, guid, project_id, role.as_deref(), now).await?;

        // The Y-STR profile from the wide overview, which is Phase 2. The code adds the profile
        // only when it makes a new subject. On a merge, the subject already holds its own data
        // sources. So the code adds the FTDNA identity, the membership, and the MDKA data above,
        // and it does not add a second Y-STR profile.
        let wrote_str = !input.ystr_markers.is_empty() && target.is_none();
        if wrote_str {
            str_profile::create(
                pool,
                &NewStrProfile {
                    biosample_guid: guid,
                    panel_name: panel_name_for_count(input.ystr_markers.len()),
                    provider: Some(IdSource::FTDNA.to_string()),
                    source: Some("IMPORTED".to_string()),
                    markers: input.ystr_markers.clone(),
                },
            )
            .await?;
        }

        // A Y-STR profile (or any Y-targeted run) means the subject is male.
        self.assign_male_for_y_evidence(guid).await?;

        Ok((wrote, wrote_str))
    }

    /// Vendor identifiers attached to a Subject.
    pub async fn external_ids(
        &self,
        guid: SampleGuid,
    ) -> Result<Vec<navigator_domain::identity::ExternalId>, AppError> {
        Ok(external_id::list_for(self.store.pool(), guid).await?)
    }

    /// The reverse of [`external_ids`]. The method returns the subject of a
    /// `(source, external_id)` vendor id, when one exists.
    ///
    /// This lookup is the exact-match anchor that finds a duplicate donor, in design §5.1. One use
    /// is to find the biosample of an FTDNA kit number. The method returns `None` when the
    /// workspace does not hold the id.
    pub async fn find_biosample_by_external_id(
        &self,
        source: &str,
        external_id: &str,
    ) -> Result<Option<SampleGuid>, AppError> {
        Ok(
            navigator_store::external_id::find(self.store.pool(), source, external_id)
                .await?
                .map(|e| e.biosample_guid),
        )
    }

    /// FTDNA-reported member labels for a Subject, if imported.
    pub async fn ftdna_member(&self, guid: SampleGuid) -> Result<Option<FtdnaMember>, AppError> {
        Ok(ftdna_member::get(self.store.pool(), guid).await?)
    }

    /// MDKA rows (paternal/maternal) for a Subject.
    pub async fn mdka_for(&self, guid: SampleGuid) -> Result<Vec<navigator_domain::identity::Mdka>, AppError> {
        Ok(mdka::list_for(self.store.pool(), guid).await?)
    }

    /// The ids of the projects that hold this subject. The method reads the M:N membership
    /// table.
    pub async fn project_membership_ids(&self, guid: SampleGuid) -> Result<Vec<i64>, AppError> {
        Ok(biosample_project::list_projects_for(self.store.pool(), guid).await?)
    }

    /// Group the members of a project by their Y-STR values, and copy an SNP branch to a member
    /// that has only STR values. The project cluster view shows this result.
    ///
    /// The branch of a member is the terminal SNP that FTDNA reports for it. The markers are the
    /// merged Y-STR profiles. The calculation costs O(n²), so it runs on its own thread.
    pub async fn cluster_project_ystr(
        &self,
        project_id: i64,
    ) -> Result<navigator_domain::ystr_cluster::YstrClustering, AppError> {
        use navigator_domain::ystr_cluster::{cluster_ystr, ClusterMember, ClusterOpts};

        // All project members (M:N membership ∪ legacy home column).
        let subjects = biosample::list_members_for_project(self.store.pool(), project_id).await?;

        let mut members = Vec::with_capacity(subjects.len());
        for b in subjects {
            let guid = b.guid;
            let fm = ftdna_member::get(self.store.pool(), guid).await?;
            let label = fm
                .as_ref()
                .and_then(|m| m.member_name.clone())
                .unwrap_or(b.donor_identifier);
            let branch = fm.and_then(|m| m.y_haplogroup_ftdna);
            // Merge the subject's Y-STR markers (dedup by name).
            let mut markers = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for p in self.list_str_profiles(guid).await? {
                for m in p.markers {
                    if seen.insert(m.marker.to_ascii_uppercase()) {
                        markers.push(m);
                    }
                }
            }
            members.push(ClusterMember {
                guid,
                label,
                branch,
                markers,
            });
        }

        tokio::task::spawn_blocking(move || cluster_ystr(&members, &ClusterOpts::default()))
            .await
            .map_err(|e| AppError::Join(e.to_string()))
    }

    /// Build a one-shot index of workspace Subjects with their Y terminal SNP + merged Y-STR markers
    /// (for fuzzy matching). Computed once to avoid O(kits × subjects) DB reads.
    async fn existing_subject_index(&self) -> Result<Vec<ExistingSubject>, AppError> {
        let mut out = Vec::new();
        for b in biosample::list_all(self.store.pool()).await? {
            let y_terminal = self
                .haplogroup_consensus(b.guid, DnaType::Y)
                .await?
                .map(|c| c.haplogroup)
                .as_deref()
                .and_then(terminal_snp);
            // Merge all of the subject's Y-STR profiles into one marker set (dedup by name).
            let mut ystr: Vec<StrMarker> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for p in self.list_str_profiles(b.guid).await? {
                for m in p.markers {
                    if seen.insert(m.marker.to_ascii_uppercase()) {
                        ystr.push(m);
                    }
                }
            }
            out.push(ExistingSubject {
                guid: b.guid,
                donor_identifier: b.donor_identifier,
                y_terminal,
                ystr,
            });
        }
        Ok(out)
    }

    /// Match one kit: exact vendor-id first (auto-merge), else fuzzy candidates, else new.
    async fn match_kit(
        &self,
        kit: &str,
        input: &FtdnaSubjectInput,
        y_terminal: Option<&str>,
        existing: &[ExistingSubject],
        threshold: f32,
    ) -> Result<MatchKind, AppError> {
        // 1. Exact vendor id → locked auto-merge.
        if let Some(hit) = external_id::find(self.store.pool(), IdSource::FTDNA, kit).await? {
            if let Some(b) = biosample::get(self.store.pool(), hit.biosample_guid).await? {
                return Ok(MatchKind::AutoMerge {
                    guid: b.guid,
                    donor_identifier: b.donor_identifier,
                });
            }
        }

        // 2. Fuzzy candidates.
        let incoming_name = input
            .member
            .as_ref()
            .and_then(|m| clean_name(m.name.as_deref()))
            .or_else(|| input.paternal.as_ref().and_then(|a| a.ancestor_name.clone()));
        let mut candidates: Vec<FuzzyCandidate> = Vec::new();
        for e in existing {
            let mut score = 0.0f32;
            let mut reasons = Vec::new();
            if let (Some(inc), Some(ex)) = (y_terminal, e.y_terminal.as_deref()) {
                if inc.eq_ignore_ascii_case(ex) {
                    score += 0.6;
                    reasons.push(format!("same Y terminal {ex}"));
                }
            }
            // The genetic distance of the Y-STR values. This distance shows the SAME PERSON only
            // when it is zero, or almost zero, across many markers.
            //
            // A high limit gives many false results in a project with one haplogroup. Each member
            // of such a project is a relative, and a distance of 3 to 11 across 100 markers is
            // normal there.
            //
            // Only an exact haplotype, or a haplotype with one difference, names the same person. A
            // larger distance names a cousin in the same clade.
            if !input.ystr_markers.is_empty() && !e.ystr.is_empty() {
                let (diff, compared) = navigator_domain::strprofile::str_distance(&input.ystr_markers, &e.ystr);
                if compared >= 67 && diff <= 1 {
                    score += 0.8;
                    reasons.push(format!("Y-STR GD {diff}/{compared}"));
                } else if compared >= 37 && diff == 0 {
                    score += 0.7;
                    reasons.push(format!("Y-STR exact ({compared} markers)"));
                }
            }
            if let Some(name) = incoming_name.as_deref() {
                let sim = name_similarity(name, &e.donor_identifier);
                if sim > 0.0 {
                    score += 0.3 * sim;
                    reasons.push("name overlap".to_string());
                }
            }
            if score >= threshold {
                candidates.push(FuzzyCandidate {
                    guid: e.guid,
                    donor_identifier: e.donor_identifier.clone(),
                    score,
                    reasons,
                });
            }
        }
        if candidates.is_empty() {
            Ok(MatchKind::New)
        } else {
            candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
            Ok(MatchKind::NeedsConfirm { candidates })
        }
    }
}

/// Workspace Subject summary used for fuzzy matching.
struct ExistingSubject {
    guid: SampleGuid,
    donor_identifier: String,
    /// The terminal SNP of the Y consensus that the app calculated for the subject. The value can
    /// be a long ISOGG label with no SNP inside it. In that case the Y-STR values are the signal
    /// that the code can trust.
    y_terminal: Option<String>,
    /// The subject's merged Y-STR markers (across all imported profiles), for genetic-distance match.
    ystr: Vec<StrMarker>,
}

fn empty_input(kit: &str) -> FtdnaSubjectInput {
    FtdnaSubjectInput {
        kit_number: kit.to_string(),
        member: None,
        paternal: None,
        maternal: None,
        ystr_markers: Vec::new(),
    }
}

/// The terminal SNP token of a haplogroup label or a clade path. The function splits the text on
/// `>` for a clade, or on `-` for a haplogroup prefix, and returns the last part. Both
/// `"R-FGC29071"` and `"CTS4466>S1115>FGC29071"` give `FGC29071`.
fn terminal_snp(label: &str) -> Option<String> {
    let t = label.rsplit(['>', '-']).next()?.trim();
    (!t.is_empty()).then(|| t.to_string())
}

/// A friendly plan label: `kit — ancestor/member name` (or just the kit).
fn display_label(kit: &str, input: &FtdnaSubjectInput) -> String {
    let name = input
        .member
        .as_ref()
        .and_then(|m| clean_name(m.name.as_deref()))
        .or_else(|| input.paternal.as_ref().and_then(|a| a.ancestor_name.clone()));
    match name {
        Some(n) => format!("{kit} — {n}"),
        None => kit.to_string(),
    }
}

/// Drop FTDNA redaction/placeholder names so they do not pollute identifiers or matching.
fn clean_name(name: Option<&str>) -> Option<String> {
    let n = name?.trim();
    if n.is_empty() || n.eq_ignore_ascii_case("REDACTED") {
        None
    } else {
        Some(n.to_string())
    }
}

/// Build an MDKA value from an ancestry row. The function returns `None` when the row holds no
/// data for the store.
fn mdka_from(a: &AncestryRow, lineage: Lineage) -> Option<NewMdka> {
    if a.ancestor_name.is_none() && a.origin_place.is_none() && a.country.is_none() && a.latitude.is_none() {
        return None;
    }
    Some(NewMdka {
        lineage: lineage.as_str().to_string(),
        ancestor_name: a.ancestor_name.clone(),
        birth_year: a.birth_year,
        death_year: a.death_year,
        origin_place: a.origin_place.clone(),
        origin_country: a.country.clone(),
        latitude: a.latitude,
        longitude: a.longitude,
        source: Some(IdSource::FTDNA.to_string()),
        notes: None,
    })
}

/// FTDNA Y-STR panel name from the count of populated markers (the standard tier boundaries).
fn panel_name_for_count(n: usize) -> String {
    let tier = [12, 25, 37, 67, 111].into_iter().find(|&t| n <= t);
    match tier {
        Some(t) => format!("Y-{t}"),
        None => "Y-700".to_string(),
    }
}

/// The `Sub Group` value of a clade, as a membership role. The function keeps the last part only,
/// and it removes the sort number at the start.
fn subgroup_role(sub_group: &str) -> String {
    terminal_snp(sub_group).unwrap_or_else(|| sub_group.trim().to_string())
}

/// The Jaccard overlap of the word tokens, in lower case, with a length of 2 or more. The value is
/// a fast measurement of the similarity of two names, from 0 to 1.
fn name_similarity(a: &str, b: &str) -> f32 {
    let toks = |s: &str| -> std::collections::HashSet<String> {
        s.split(|c: char| !c.is_ascii_alphanumeric())
            .filter(|t| t.len() >= 2)
            .map(|t| t.to_ascii_lowercase())
            .collect()
    };
    let (ta, tb) = (toks(a), toks(b));
    if ta.is_empty() || tb.is_empty() {
        return 0.0;
    }
    let inter = ta.intersection(&tb).count() as f32;
    let union = ta.union(&tb).count() as f32;
    inter / union
}

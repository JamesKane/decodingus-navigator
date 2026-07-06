//! FTDNA project import — the matching/dedup engine + two-phase plan/commit (design §5/§6).
//!
//! Phase 1 scope (roster + ancestry, the spine): parse the batch CSVs, join by kit number, match
//! each kit against the workspace, and produce a reviewable **plan** (dry-run, no writes). A separate
//! commit step applies the plan with the admin's resolutions for fuzzy candidates.
//!
//! Deep per-member data (Big Y / mtDNA / Family Finder) and the wide Y-STR chart are layered on by
//! later slices; this module only wires identity + MDKA + membership.

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

/// Tuning for the matching engine.
#[derive(Debug, Clone)]
pub struct FtdnaImportOptions {
    /// Minimum fuzzy score (0..1) for a workspace Subject to be offered as a merge candidate.
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

/// How the matcher proposes to handle a kit. Auto-merge is locked for an exact vendor-id hit; fuzzy
/// hits are queued for the admin (never auto-merged).
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

/// Recognized-input + scan counts for the review header — so a missing/misclassified file (e.g. no
/// roster) is immediately visible rather than silently producing all-orphan rows.
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
    /// Don't import this kit at all.
    Skip,
}

/// What the commit did.
#[derive(Debug, Clone, Default)]
pub struct FtdnaImportSummary {
    /// The project the kits were imported into (resolved/created at commit).
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

/// A Subject's imported genealogy bundle: vendor ids, FTDNA member labels, and MDKA rows. PII —
/// for local display only (never federated). Empty when nothing was imported for the Subject.
#[derive(Debug, Clone, Default)]
pub struct FtdnaGenealogy {
    pub external_ids: Vec<navigator_domain::identity::ExternalId>,
    pub member: Option<navigator_domain::identity::FtdnaMember>,
    pub mdka: Vec<navigator_domain::identity::Mdka>,
}

impl FtdnaGenealogy {
    /// Nothing imported → the detail card can be skipped.
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

    /// Attach a vendor id (kit number) to a Subject from the subject editor. Rejects a blank
    /// source/id, and refuses to bind a `(source, external_id)` that already belongs to a *different*
    /// Subject (the `(source, external_id)` uniqueness is the dedup anchor — never silently re-point
    /// it; the caller resolves the conflict). Idempotent for the same Subject.
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
        // Recover the owning subject before the row is gone, so we can refresh its published record.
        let guid = external_id::get(self.store.pool(), id).await?.map(|e| e.biosample_guid);
        external_id::delete(self.store.pool(), id).await?;
        if let Some(guid) = guid {
            let _ = self.republish_biosample_ids(guid).await;
        }
        Ok(())
    }

    /// Backfill public-catalog external ids (`IGSR`/`HGDP`/INSDC) derivable from each subject's local
    /// provenance ([`navigator_domain::identity::catalog_ids_from_provenance`]) — so bulk-imported
    /// public datasets publish ids that match their existing AppView catalog rows. Deterministic and
    /// network-free; a friendly-name-only sample contributes nothing. Idempotent (skips ids already
    /// present); a `(namespace, value)` already owned by a *different* subject is counted as a
    /// conflict and left untouched (never silently re-pointed). `apply == false` is a dry run.
    /// Adds via the store directly (no per-id re-publish); re-publish the affected subjects after.
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
            let derived =
                navigator_domain::identity::catalog_ids_from_provenance(&b.donor_identifier, b.sample_accession.as_deref());
            if derived.is_empty() {
                continue;
            }
            out.subjects_matched += 1;
            let existing: std::collections::HashSet<(String, String)> = external_id::list_for(self.store.pool(), b.guid)
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
                        // (namespace,value) already belongs to another subject — a dup import; leave it.
                        out.conflicts += 1;
                    }
                }
            }
        }
        Ok(out)
    }

    /// Fetch one public-catalog sample record from the AppView samples API (`/api/v1/samples/{alias}`,
    /// public read) by its alias (= our `donor_identifier`). `Ok(None)` for a 404 (alias unknown to
    /// the catalog — expected while server-side corrections are pending). The authoritative
    /// `accession` it returns is the datum our local `sample_accession` lacks.
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

    /// Resolve each subject against the AppView samples API and attach, **in one pass**, the full set
    /// of public-catalog ids: the catalog *name* id (IGSR/HGDP, derived from the donor id) **and** the
    /// authoritative INSDC *accession* the API returns (`SAMN…` → BIOSAMPLE, `ERS…` → ENA, `SRS…` →
    /// SRA) — plus correcting the local `sample_accession` placeholder. A superset of
    /// [`backfill_catalog_ids`](Self::backfill_catalog_ids) (the offline name-only path); use that one
    /// when the API is unavailable. By default only subjects whose `donor_identifier` looks like a
    /// catalog alias (IGSR/HGDP) are queried (`all` overrides), to avoid hammering the API with
    /// friendly-name 404s. `apply == false` is a dry run. `limit` caps how many are queried.
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
            // Skip samples whose name isn't a recognizable catalog alias unless `--all`.
            if !all
                && navigator_domain::identity::catalog_ids_from_provenance(&b.donor_identifier, None).is_empty()
            {
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
            // One pass: the catalog *name* id (from the donor id) + the authoritative INSDC *accession*
            // (from the API, when it's a real one) — the union of both sources via the shared helper.
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
            let existing: std::collections::HashSet<(String, String)> = external_id::list_for(self.store.pool(), b.guid)
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

    /// Re-publish a subject's biosample anchor after its identifier set changed, so the AppView's
    /// mirror (which full-replaces `external_ids`) honors the add/remove. Deterministic rkey → the
    /// re-publish overwrites in place. **Only for a subject already federated** and while signed in —
    /// signed out, or a never-published subject, is a no-op (we don't newly federate a donor just
    /// because a local id was attached).
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

    /// Insert or update a Subject's MDKA for one lineage from the subject editor (one row per
    /// lineage; stamps `updated_at`). Pass a `source` of `MANUAL` for hand-entered rows.
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
        // A roster was provided iff there are member rows — only then is "orphan" (data without a
        // roster row) a meaningful flag.
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
                // Orphan only when a roster was provided but this kit isn't in it.
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

    /// Apply a plan. `resolutions` carries the admin's choice for each fuzzy (`NeedsConfirm`) kit;
    /// an unresolved fuzzy row defaults to **New** (conservative — never silently merge).
    pub async fn commit_ftdna_import(
        &self,
        plan: &FtdnaImportPlan,
        resolutions: &BTreeMap<String, FtdnaResolution>,
    ) -> Result<FtdnaImportSummary, AppError> {
        let mut summary = FtdnaImportSummary::default();
        let now = Utc::now().to_rfc3339();

        // Resolve the target project, creating it now if the plan targeted a new one.
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

        // Vendor identity (idempotent; never steals a conflicting id).
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

        // MDKA from paternal (Y) + maternal (Mt) ancestry, when there's anything worth storing.
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

        // Y-STR profile from the wide overview (Phase 2). Attached only when CREATING a new Subject;
        // on a merge the existing Subject already carries its own data sources, so we add the FTDNA
        // identity/membership/MDKA metadata above but skip duplicating the Y-STR profile.
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

    /// Reverse of [`external_ids`]: the Subject bound to a `(source, external_id)` vendor id, if any.
    /// This is the exact-match dedup anchor (design §5.1) — e.g. resolve an FTDNA kit number to the
    /// biosample it was imported under. Returns `None` when the id is unknown to the workspace.
    pub async fn find_biosample_by_external_id(
        &self,
        source: &str,
        external_id: &str,
    ) -> Result<Option<SampleGuid>, AppError> {
        Ok(navigator_store::external_id::find(self.store.pool(), source, external_id)
            .await?
            .map(|e| e.biosample_guid))
    }

    /// FTDNA-reported member labels for a Subject, if imported.
    pub async fn ftdna_member(&self, guid: SampleGuid) -> Result<Option<FtdnaMember>, AppError> {
        Ok(ftdna_member::get(self.store.pool(), guid).await?)
    }

    /// MDKA rows (paternal/maternal) for a Subject.
    pub async fn mdka_for(&self, guid: SampleGuid) -> Result<Vec<navigator_domain::identity::Mdka>, AppError> {
        Ok(mdka::list_for(self.store.pool(), guid).await?)
    }

    /// Project ids a Subject belongs to (via the M:N membership table).
    pub async fn project_membership_ids(&self, guid: SampleGuid) -> Result<Vec<i64>, AppError> {
        Ok(biosample_project::list_projects_for(self.store.pool(), guid).await?)
    }

    /// Autocluster a project's members by Y-STR and propagate SNP branches to STR-only members
    /// (the project clustering view). Branch per member = its FTDNA-reported terminal SNP; markers =
    /// the merged Y-STR profiles. The O(n²) compute runs on a blocking thread.
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
            // Y-STR genetic distance — a SAME-PERSON signal only at (near-)zero GD over many markers.
            // A loose GD threshold floods inside a single-haplogroup project, where every member is
            // related and within-project distances of GD 3–11 over 100 markers are normal. Only an
            // exact (or off-by-one) haplotype uniquely identifies the same person, not a clade cousin.
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
    /// Terminal SNP of the subject's computed Y consensus (may be an ISOGG long-form label that
    /// doesn't reduce to an SNP — then Y-STR is the reliable signal).
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

/// The terminal SNP token of a haplogroup label or clade path: the last segment after splitting on
/// `>` (clade) or `-` (haplogroup prefix). `"R-FGC29071"` and `"CTS4466>S1115>FGC29071"` → `FGC29071`.
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

/// Drop FTDNA redaction/placeholder names so they don't pollute identifiers or matching.
fn clean_name(name: Option<&str>) -> Option<String> {
    let n = name?.trim();
    if n.is_empty() || n.eq_ignore_ascii_case("REDACTED") {
        None
    } else {
        Some(n.to_string())
    }
}

/// Build an MDKA payload from an ancestry row, or `None` if it carries nothing worth storing.
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

/// The clade `Sub Group` value as a membership role: keep it compact (the terminal segment), dropping
/// the leading sort number.
fn subgroup_role(sub_group: &str) -> String {
    terminal_snp(sub_group).unwrap_or_else(|| sub_group.trim().to_string())
}

/// Jaccard overlap of lowercased word tokens (len ≥ 2) — a cheap name-similarity proxy in `0..=1`.
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

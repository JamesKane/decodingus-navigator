//! Navigator application/command layer — the single API the UI dispatches to, and the
//! antidote to the `WorkbenchViewModel` god object. Orchestrates `navigator-store` (and
//! later analysis/sync) behind commands and queries; holds policy the old dialogs
//! embedded (identity assignment, existence checks, result (de)serialization). The UI
//! holds only view-state and dispatch — no DB calls or domain decisions in widgets.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

use chrono::Utc;
use du_domain::ids::SampleGuid;
use navigator_analysis::caller::{self, HaploidCallerParams, SiteGenotype, Site, VariantCall};
use navigator_analysis::coverage::{self, CallableLociParams, CoverageResult};
use navigator_analysis::ibd::{
    ChromosomeGenotypes, GeneticMap, IbdSegment, MatchSummary, PairwiseIbdDetector,
};
use navigator_domain::workspace::{Panel, PanelSite};
use navigator_store::panel;

// Re-export the analysis result types the command API returns, so the UI depends only
// on navigator-app (ui -> app), not directly on navigator-analysis.
pub use navigator_analysis::caller::SiteGenotype as PanelGenotype;
pub use navigator_analysis::caller::VariantCall as DenovoCall;
pub use navigator_analysis::coverage::CoverageResult as Coverage;
pub use navigator_analysis::haplo::{BranchEvidence, CallState, ScoredHaplogroup, SnpEvidence};

/// A haplogroup assignment: the ranked candidates plus, for the reported terminal, the
/// child branches with per-SNP evidence (why descent stopped — unsupported splits show
/// ancestral SNPs, unresolved ones show no-calls).
#[derive(Debug, Clone)]
pub struct HaploAssignment {
    pub ranked: Vec<ScoredHaplogroup>,
    pub branches: Vec<BranchEvidence>,
}

/// How a private (off-backbone) variant relates to the tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrivateClass {
    /// A known tree SNP off the assigned path — supports a finer/sibling branch.
    OffPathKnown(String),
    /// Not in the tree at all — a candidate for proposing a new branch.
    Novel,
}

/// A derived variant the sample carries that the haplogroup placement doesn't explain.
#[derive(Debug, Clone)]
pub struct PrivateVariant {
    pub position: i64,
    pub reference: char,
    pub alternate: char,
    pub depth: u32,
    pub allele_fraction: f64,
    pub class: PrivateClass,
}

/// The private bucket for an alignment: de-novo Y calls not on the assigned backbone,
/// split into off-path-known (finer branches) and novel (new-branch candidates).
#[derive(Debug, Clone)]
pub struct PrivateBucket {
    pub terminal: String,
    pub variants: Vec<PrivateVariant>,
}

impl PrivateBucket {
    pub fn novel(&self) -> usize {
        self.variants.iter().filter(|v| v.class == PrivateClass::Novel).count()
    }
    pub fn off_path(&self) -> usize {
        self.variants.iter().filter(|v| matches!(v.class, PrivateClass::OffPathKnown(_))).count()
    }
}
pub use navigator_analysis::ibd::{
    IbdDetectorConfig, IbdSegment as Segment, MatchSummary as IbdSummary, RelationshipEstimate,
};
// Sync/publish types the command API uses, re-exported so the UI depends only on navigator-app.
pub use navigator_sync::{
    CoverageSummaryRecord, PdsClient, PrivateVariantsRecord, RecordRef, VariantCallEntry,
    COVERAGE_SUMMARY_COLLECTION, PRIVATE_VARIANTS_COLLECTION,
};
use navigator_sync::{dev_http_client, login_default, AsyncSync, OAuthConfig, RetryPolicy, TokenStore};

/// Keychain service namespace for stored sessions (plan §7).
const KEYCHAIN_SERVICE: &str = "decodingus-navigator";

/// IBD comparison result between two genotyped alignments.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdComparison {
    pub summary: MatchSummary,
    pub segments: Vec<IbdSegment>,
}
use navigator_domain::workspace::{
    Alignment, AnalysisArtifact, Biosample, NewAlignment, NewProject, NewSequenceRun, Project,
    SequenceRun,
};
use navigator_domain::chipprofile::{self, ChipProfile, NewChipProfile};
use navigator_domain::filetype;
pub use navigator_domain::filetype::DetectedData;
use navigator_domain::mtdna::{self, MtdnaSequence, NewMtdnaSequence};
use navigator_domain::reconciliation::{self, RunHaplogroupCall};
pub use navigator_domain::reconciliation::{CompatibilityLevel, Consensus, DnaType};
use navigator_domain::strprofile::{self, NewStrProfile, StrProfile};
use navigator_domain::variants::{self, NewVariantSet, VariantSet};
use navigator_store::{
    alignment, artifact, biosample, chip_profile, haplogroup_call, mtdna as mtdna_store, project,
    sequence_run, str_profile, variant_set, Store, StoreError,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use uuid::Uuid;

pub mod error;
pub use error::AppError;

/// Artifact kind for de-novo calls, keyed per contig so different contigs don't
/// overwrite each other in the cache.
fn denovo_kind(contig: &str) -> String {
    format!("denovo_snps:{contig}")
}

/// On-disk cache path for a downloaded haplotree, under `$NAVIGATOR_TREE_DIR` (tests/
/// overrides) or `~/.decodingus/trees`.
fn tree_cache_path(file: &str) -> PathBuf {
    let dir = std::env::var("NAVIGATOR_TREE_DIR").ok().map(PathBuf::from).unwrap_or_else(|| {
        let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
        PathBuf::from(home).join(".decodingus").join("trees")
    });
    dir.join(file)
}

/// Score a tree against the sample calls and attach the terminal's child-branch evidence.
fn assemble_assignment(tree: &navigator_analysis::haplo::HaploTree, calls: &HashMap<i64, char>) -> HaploAssignment {
    let ranked = navigator_analysis::haplo::score(tree, calls);
    let branches = ranked
        .first()
        .map(|t| navigator_analysis::haplo::child_evidence(tree, calls, t.id))
        .unwrap_or_default();
    HaploAssignment { ranked, branches }
}

/// Haploid-caller params adapted to the sample's read tech: long, accurate reads (HiFi,
/// mean read length > 1 kb) make confident haploid calls at much lower depth, so halve
/// `min_depth` (floor 2). Sampled from the BAM head; falls back to defaults on any error.
/// Blocking (reads the BAM) — call inside `spawn_blocking`.
fn adaptive_haploid_params(bam_path: &Path) -> HaploidCallerParams {
    let mut params = HaploidCallerParams::default();
    if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(bam_path) {
        if read_len > 1000.0 {
            params.min_depth = (params.min_depth / 2).max(2);
        }
    }
    params
}

/// Read the first 64 KiB of a file as lossy UTF-8 — enough to fingerprint a text file's
/// type without slurping a multi-MB chip export.
fn read_head(path: &Path) -> Result<String, AppError> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut buf = vec![0u8; 64 * 1024];
    let n = f.read(&mut buf)?;
    buf.truncate(n);
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Artifact kind for panel genotypes, keyed by panel + ploidy.
fn panel_kind(panel_id: i64, ploidy: u8) -> String {
    format!("panel:{panel_id}:p{ploidy}")
}

/// Group per-site genotypes into per-chromosome dosage arrays (sorted by position) for
/// the IBD detector.
fn group_chrom_genotypes(genotypes: &[SiteGenotype]) -> std::collections::HashMap<String, ChromosomeGenotypes> {
    let mut by_contig: BTreeMap<String, Vec<(i64, i32)>> = BTreeMap::new();
    for g in genotypes {
        by_contig.entry(g.contig.clone()).or_default().push((g.position, g.dosage));
    }
    by_contig
        .into_iter()
        .map(|(chrom, mut v)| {
            v.sort_by_key(|(p, _)| *p);
            let positions = v.iter().map(|(p, _)| *p as i32).collect();
            let dosages = v.iter().map(|(_, d)| *d as i8).collect();
            (chrom.clone(), ChromosomeGenotypes { chromosome: chrom, positions, dosages })
        })
        .collect()
}

/// A project plus a rolled-up count for list/dashboard views.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectOverview {
    pub project: Project,
    pub sample_count: i64,
}

/// AT Proto auth state: keychain-backed sessions + the in-memory active account. Shared
/// (cheaply cloned with the `App`); the active DID is the only mutable bit.
#[derive(Clone)]
struct Auth {
    tokens: TokenStore,
    config: OAuthConfig,
    http: reqwest::Client,
    /// The signed-in account's DID, or `None`. `Arc<Mutex>` so clones of `App` share it.
    active: Arc<Mutex<Option<String>>>,
    /// Offline indicator shared with every [`AsyncSync`] this app builds: cleared on a
    /// transient write failure, set on success. Starts optimistic (`true`).
    online: Arc<AtomicBool>,
}

impl Auth {
    fn new() -> Self {
        let tokens = TokenStore::new(KEYCHAIN_SERVICE);
        // Reload whoever was signed in last launch; a keychain error just means "nobody".
        let active = tokens.active().ok().flatten();
        Auth {
            tokens,
            config: OAuthConfig::loopback("atproto"),
            http: dev_http_client(),
            active: Arc::new(Mutex::new(active)),
            online: Arc::new(AtomicBool::new(true)),
        }
    }
}

/// The application. Cheap to clone (the store wraps a connection pool).
#[derive(Clone)]
pub struct App {
    store: Store,
    auth: Auth,
}

impl App {
    pub fn new(store: Store) -> Self {
        App { store, auth: Auth::new() }
    }

    /// Open/create the workspace database and build the app.
    pub async fn open(path: &std::path::Path) -> Result<Self, AppError> {
        Ok(App::new(Store::open(path).await?))
    }

    // ---- commands ----------------------------------------------------------

    pub async fn create_project(&self, new: NewProject) -> Result<Project, AppError> {
        Ok(project::create(self.store.pool(), &new).await?)
    }

    /// Register a biosample, assigning its stable `SampleGuid` here (identity is an
    /// app-layer decision, not the UI's). Verifies the target project exists first so
    /// the caller gets a clear `NotFound` rather than a raw foreign-key error.
    pub async fn add_biosample(
        &self,
        project_id: Option<i64>,
        donor_identifier: impl Into<String>,
        sample_accession: Option<String>,
        sex: Option<String>,
    ) -> Result<Biosample, AppError> {
        if let Some(pid) = project_id {
            if project::get(self.store.pool(), pid).await?.is_none() {
                return Err(AppError::Store(StoreError::NotFound(format!("project {pid}"))));
            }
        }
        let b = Biosample {
            guid: SampleGuid(Uuid::new_v4()),
            sample_accession,
            donor_identifier: donor_identifier.into(),
            description: None,
            center_name: None,
            sex,
            project_id,
        };
        biosample::create(self.store.pool(), &b).await?;
        Ok(b)
    }

    pub async fn record_sequence_run(&self, run: NewSequenceRun) -> Result<SequenceRun, AppError> {
        Ok(sequence_run::create(self.store.pool(), &run).await?)
    }

    pub async fn record_alignment(&self, aln: NewAlignment) -> Result<Alignment, AppError> {
        Ok(alignment::create(self.store.pool(), &aln).await?)
    }

    /// Persist a typed analysis result as a versioned artifact (JSON payload). The
    /// `algorithm_version` is part of the cache key, so a newer version supersedes the
    /// old entry. Pair with [`App::load_analysis`].
    pub async fn save_analysis<T: Serialize>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
        result: &T,
    ) -> Result<AnalysisArtifact, AppError> {
        let payload = serde_json::to_string(result)?;
        Ok(artifact::upsert(self.store.pool(), alignment_id, kind, algorithm_version, Utc::now(), &payload).await?)
    }

    /// Load and deserialize a stored analysis result, if present for this version.
    pub async fn load_analysis<T: DeserializeOwned>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
    ) -> Result<Option<T>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, kind, algorithm_version).await? {
            Some(a) => Ok(Some(serde_json::from_str(&a.payload)?)),
            None => Ok(None),
        }
    }

    // ---- analysis (compute + persist) --------------------------------------

    /// Run the coverage + callable walker on an alignment's BAM and persist the result
    /// as a versioned `coverage` artifact. The blocking noodles I/O runs on a blocking
    /// thread so the async runtime is not stalled.
    pub async fn run_coverage(
        &self,
        alignment_id: i64,
        bam: PathBuf,
        reference: PathBuf,
        contig_allowlist: Option<HashSet<String>>,
        params: CallableLociParams,
    ) -> Result<CoverageResult, AppError> {
        let result = tokio::task::spawn_blocking(move || {
            coverage::collect_coverage_callable(&bam, &reference, &params, contig_allowlist.as_ref())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result).await?;
        Ok(result)
    }

    /// Cached `coverage` result for the current algorithm version, if present.
    pub async fn cached_coverage(&self, alignment_id: i64) -> Result<Option<CoverageResult>, AppError> {
        self.load_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION).await
    }

    /// Run coverage using the alignment's own stored BAM/reference paths, then persist.
    /// Errors if the alignment is unknown or has no paths recorded.
    pub async fn run_coverage_for_alignment(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let (Some(bam), Some(reference)) = (aln.bam_path, aln.reference_path) else {
            return Err(AppError::MissingPaths(alignment_id));
        };
        self.run_coverage(
            alignment_id,
            PathBuf::from(bam),
            PathBuf::from(reference),
            None,
            CallableLociParams::default(),
        )
        .await
    }

    /// Run de-novo haploid calling on a contig and persist the SNP calls as a versioned
    /// `denovo_snps` artifact.
    pub async fn run_denovo_caller(
        &self,
        alignment_id: i64,
        bam: PathBuf,
        reference: PathBuf,
        contig: String,
        params: HaploidCallerParams,
    ) -> Result<Vec<VariantCall>, AppError> {
        let kind = denovo_kind(&contig);
        let calls = tokio::task::spawn_blocking(move || {
            caller::call_denovo(&bam, &reference, &contig, &params)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, &kind, caller::DENOVO_VERSION, &calls).await?;
        Ok(calls)
    }

    /// Cached de-novo calls for `contig` at the current caller version, if present.
    pub async fn cached_denovo(&self, alignment_id: i64, contig: &str) -> Result<Option<Vec<VariantCall>>, AppError> {
        self.load_analysis(alignment_id, &denovo_kind(contig), caller::DENOVO_VERSION).await
    }

    /// Run de-novo calling on `contig` using the alignment's own stored paths.
    pub async fn run_denovo_for_alignment(&self, alignment_id: i64, contig: String) -> Result<Vec<VariantCall>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let (Some(bam), Some(reference)) = (aln.bam_path, aln.reference_path) else {
            return Err(AppError::MissingPaths(alignment_id));
        };
        let bam = PathBuf::from(bam);
        let probe = bam.clone();
        let params = tokio::task::spawn_blocking(move || adaptive_haploid_params(&probe))
            .await
            .map_err(|e| AppError::Join(e.to_string()))?; // HiFi -> lower min_depth
        self.run_denovo_caller(alignment_id, bam, PathBuf::from(reference), contig, params).await
    }

    // ---- publish -----------------------------------------------------------

    /// Build the coverage-summary record JSON for an alignment (floats encoded as strings).
    async fn coverage_record(&self, alignment_id: i64) -> Result<serde_json::Value, AppError> {
        let cov = self
            .cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))?;
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let record = CoverageSummaryRecord::new(
            aln.reference_build,
            cov.mean_coverage,
            cov.median_coverage,
            cov.sd_coverage,
            cov.pct_10x,
            cov.pct_20x,
            cov.pct_30x,
            cov.genome_territory,
            cov.callable_bases,
            Utc::now().to_rfc3339(),
        );
        Ok(serde_json::to_value(&record)?)
    }

    /// Build the private-variants record JSON for an alignment's cached de-novo calls.
    async fn variants_record(&self, alignment_id: i64, contig: &str) -> Result<serde_json::Value, AppError> {
        let calls = self
            .cached_denovo(alignment_id, contig)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("de-novo calls for alignment {alignment_id} {contig}"))))?;
        let variants = calls
            .iter()
            .map(|c| VariantCallEntry::new(c.position, c.reference_allele, c.alternate_allele, c.depth, c.alt_depth, c.allele_fraction))
            .collect();
        let record = PrivateVariantsRecord::new(contig, caller::DENOVO_VERSION, Utc::now().to_rfc3339(), variants);
        Ok(serde_json::to_value(&record)?)
    }

    /// Publish an alignment's cached coverage summary using an explicit `client` (the
    /// testable core; production callers use [`publish_coverage`](Self::publish_coverage)).
    pub async fn publish_coverage_summary(
        &self,
        client: &PdsClient,
        alignment_id: i64,
    ) -> Result<RecordRef, AppError> {
        let value = self.coverage_record(alignment_id).await?;
        Ok(client.create_record(COVERAGE_SUMMARY_COLLECTION, value, None).await?)
    }

    /// Publish an alignment's cached de-novo calls for `contig` using an explicit `client`
    /// (the testable core; production callers use [`publish_variants`](Self::publish_variants)).
    pub async fn publish_private_variants(
        &self,
        client: &PdsClient,
        alignment_id: i64,
        contig: &str,
    ) -> Result<RecordRef, AppError> {
        let value = self.variants_record(alignment_id, contig).await?;
        Ok(client.create_record(PRIVATE_VARIANTS_COLLECTION, value, None).await?)
    }

    // ---- authentication ----------------------------------------------------

    /// Run the public-client OAuth login for `handle` (handle or DID): browser authorize →
    /// loopback callback → token exchange. On success the DPoP-bound session is persisted
    /// to the OS keychain and becomes the active account. Returns the authenticated DID.
    pub async fn login(&self, handle: &str) -> Result<String, AppError> {
        let session = login_default(&self.auth.http, &self.auth.config, handle).await?;
        let did = session.did.clone();
        self.auth.tokens.save(&did, &session)?;
        self.auth.tokens.set_active(&did)?;
        *self.auth.active.lock().unwrap() = Some(did.clone());
        Ok(did)
    }

    /// The signed-in account's DID, or `None`.
    pub fn current_account(&self) -> Option<String> {
        self.auth.active.lock().unwrap().clone()
    }

    /// Sign out: drop the active account and delete its stored session.
    pub async fn logout(&self) -> Result<(), AppError> {
        let did = self.auth.active.lock().unwrap().take();
        if let Some(did) = did {
            self.auth.tokens.delete(&did)?;
        }
        self.auth.tokens.clear_active()?;
        Ok(())
    }

    /// Build the resilient sync engine for the active account, loading its session from the
    /// keychain. Errors with [`AppError::NotAuthenticated`] when no one is signed in. The
    /// engine auto-refreshes on 401 and retries transient failures with backoff.
    fn sync_engine(&self) -> Result<AsyncSync, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let session = self.auth.tokens.load(&did)?.ok_or(AppError::NotAuthenticated)?;
        Ok(AsyncSync::new(
            self.auth.http.clone(),
            self.auth.tokens.clone(),
            session,
            RetryPolicy::default(),
            self.auth.online.clone(),
        ))
    }

    /// Whether the last PDS write reached the server. Drives the UI's offline indicator;
    /// optimistic (`true`) until a transient write failure.
    pub fn is_online(&self) -> bool {
        self.auth.online.load(Ordering::Relaxed)
    }

    /// Publish the alignment's coverage summary to the signed-in account's PDS (with
    /// refresh-on-expiry and retry/backoff via [`AsyncSync`]).
    pub async fn publish_coverage(&self, alignment_id: i64) -> Result<RecordRef, AppError> {
        let mut engine = self.sync_engine()?; // auth check before touching the DB
        let value = self.coverage_record(alignment_id).await?;
        Ok(engine.push_create(COVERAGE_SUMMARY_COLLECTION, value).await?)
    }

    /// Publish the alignment's de-novo calls for `contig` to the signed-in account's PDS.
    pub async fn publish_variants(&self, alignment_id: i64, contig: &str) -> Result<RecordRef, AppError> {
        let mut engine = self.sync_engine()?; // auth check before touching the DB
        let value = self.variants_record(alignment_id, contig).await?;
        Ok(engine.push_create(PRIVATE_VARIANTS_COLLECTION, value).await?)
    }

    // ---- panels + IBD ------------------------------------------------------

    /// Create a genotyping panel from explicit sites.
    pub async fn import_panel(&self, name: &str, sites: &[PanelSite]) -> Result<Panel, AppError> {
        Ok(panel::create(self.store.pool(), name, sites).await?)
    }

    /// Create a panel from a (plain-text) sites VCF — biallelic SNP rows only.
    pub async fn import_panel_from_vcf(&self, name: &str, vcf_path: &Path) -> Result<Panel, AppError> {
        let variants = navigator_analysis::parity::parse_truth_vcf(vcf_path)?;
        let sites: Vec<PanelSite> = variants
            .iter()
            .filter_map(|v| {
                let alt = v.alternate.first()?;
                (v.reference.len() == 1 && alt.len() == 1).then(|| PanelSite {
                    chrom: v.chrom.clone(),
                    position: v.pos,
                    reference_allele: v.reference.clone(),
                    alternate_allele: alt.clone(),
                    name: v.ids.first().cloned().unwrap_or_else(|| format!("{}:{}", v.chrom, v.pos)),
                })
            })
            .collect();
        self.import_panel(name, &sites).await
    }

    pub async fn list_panels(&self) -> Result<Vec<Panel>, AppError> {
        Ok(panel::list(self.store.pool()).await?)
    }

    // ---- STR profiles ------------------------------------------------------

    /// Import a Y-STR profile for a subject from an exported marker table (CSV/TSV).
    pub async fn import_str_profile_from_csv(
        &self,
        biosample_guid: SampleGuid,
        panel_name: &str,
        provider: Option<String>,
        source: Option<String>,
        csv_path: &Path,
    ) -> Result<StrProfile, AppError> {
        let text = std::fs::read_to_string(csv_path)?;
        let markers = strprofile::parse_csv(&text).map_err(AppError::Import)?;
        let new = NewStrProfile { biosample_guid, panel_name: panel_name.to_string(), provider, source, markers };
        Ok(str_profile::create(self.store.pool(), &new).await?)
    }

    /// All STR profiles for a subject.
    pub async fn list_str_profiles(&self, biosample_guid: SampleGuid) -> Result<Vec<StrProfile>, AppError> {
        Ok(str_profile::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- SNP variants ------------------------------------------------------

    /// Import a subject's SNP variant calls from a file. `.vcf` is parsed as a VCF (reusing
    /// the shared column parser); `.csv`/`.tsv` as a `contig,position,ref,alt[,rsid][,gt]`
    /// table. Indels/symbolic alleles are dropped (SNP-only). `source_label` defaults to the
    /// file name.
    pub async fn import_variants_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<VariantSet, AppError> {
        let label = path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "variants".into());
        let is_vcf = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("vcf"))
            .unwrap_or(false);

        let calls = if is_vcf {
            navigator_analysis::parity::parse_truth_vcf(path)?
                .into_iter()
                .filter_map(|v| {
                    let alt = v.alternate.first()?;
                    variants::snp_call(&v.chrom, v.pos, &v.reference, alt, v.ids.first().cloned(), None)
                })
                .collect()
        } else {
            let text = std::fs::read_to_string(path)?;
            variants::parse_csv(&text).map_err(AppError::Import)?
        };
        if calls.is_empty() {
            return Err(AppError::Import("no SNP variants found in file".into()));
        }
        let new = NewVariantSet { biosample_guid, source_label: label, calls };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// All variant sets for a subject.
    pub async fn list_variant_sets(&self, biosample_guid: SampleGuid) -> Result<Vec<VariantSet>, AppError> {
        Ok(variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- chip / array profiles ---------------------------------------------

    /// Import a genotyping-array raw-data export (CSV/TSV) and store its QC summary.
    /// `provider` overrides vendor detection when given; `chip_version` is optional.
    pub async fn import_chip_profile_from_csv(
        &self,
        biosample_guid: SampleGuid,
        provider: Option<String>,
        chip_version: Option<String>,
        path: &Path,
    ) -> Result<ChipProfile, AppError> {
        let text = std::fs::read_to_string(path)?;
        let (summary, detected) = chipprofile::summarize(&text).map_err(AppError::Import)?;
        let provider = provider.or(detected).unwrap_or_else(|| "OTHER".into());
        let source_file_name = path.file_name().map(|s| s.to_string_lossy().into_owned());
        let new = NewChipProfile { biosample_guid, provider, chip_version, summary, source_file_name };
        Ok(chip_profile::create(self.store.pool(), &new).await?)
    }

    /// All chip profiles for a subject.
    pub async fn list_chip_profiles(&self, biosample_guid: SampleGuid) -> Result<Vec<ChipProfile>, AppError> {
        Ok(chip_profile::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- mtDNA sequences ---------------------------------------------------

    /// Import a vendor mtDNA FASTA (~16,569 bp) for a subject. Validates the header,
    /// length, and bases; stores the sequence + N count.
    pub async fn import_mtdna_from_fasta(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<MtdnaSequence, AppError> {
        let text = std::fs::read_to_string(path)?;
        let parsed = mtdna::parse_fasta(&text).map_err(AppError::Import)?;
        let source_file_name = path.file_name().map(|s| s.to_string_lossy().into_owned());
        let new = NewMtdnaSequence {
            biosample_guid,
            defline: parsed.defline,
            sequence: parsed.sequence,
            n_count: parsed.n_count,
            source_file_name,
        };
        Ok(mtdna_store::create(self.store.pool(), &new).await?)
    }

    /// All mtDNA sequences for a subject.
    pub async fn list_mtdna_sequences(&self, biosample_guid: SampleGuid) -> Result<Vec<MtdnaSequence>, AppError> {
        Ok(mtdna_store::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Derive mtDNA variants for a stored sequence by comparing it to an rCRS reference
    /// FASTA, and save them as a variant set (contig `rCRS`) so they appear alongside the
    /// subject's other variants. The reference is validated as an mtDNA FASTA.
    pub async fn derive_mtdna_variants(&self, mtdna_id: i64, rcrs_path: &Path) -> Result<VariantSet, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;
        let rcrs_text = std::fs::read_to_string(rcrs_path)?;
        let rcrs = mtdna::parse_fasta(&rcrs_text).map_err(|e| AppError::Import(format!("rCRS reference: {e}")))?;

        let derived = navigator_analysis::mtvariants::derive(&rcrs.sequence, &seq.sequence);
        let calls = derived
            .iter()
            .map(|v| variants::VariantCall {
                contig: "rCRS".to_string(),
                position: v.position,
                reference: v.reference.to_string(),
                alternate: v.alternate.to_string(),
                rs_id: None,
                genotype: None,
            })
            .collect();
        let label = format!("mtDNA vs rCRS ({} variants)", derived.len());
        let new = NewVariantSet { biosample_guid: seq.biosample_guid, source_label: label, calls };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// Assign an mtDNA haplogroup to a stored sequence: fetch (and cache) the FTDNA mt-DNA
    /// haplotree and rank haplogroups by the Kulczynski measure over the sample's base
    /// calls. RSRS-anchored and reference-free (no rCRS needed). Best first.
    pub async fn assign_mtdna_haplogroup(&self, mtdna_id: i64) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let assignment = self.assign_mtdna_haplogroup_with_tree(mtdna_id, &tree_json).await?;
        if let Some(seq) = mtdna_store::get(self.store.pool(), mtdna_id).await? {
            self.record_call(seq.biosample_guid, DnaType::Mt, &format!("mtseq:{mtdna_id}"), format!("mtDNA seq #{mtdna_id}"), &assignment).await?;
        }
        Ok(assignment)
    }

    /// The biosample a alignment belongs to (alignment → sequencing run → biosample).
    async fn biosample_of_alignment(&self, alignment_id: i64) -> Result<SampleGuid, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let run = sequence_run::get(self.store.pool(), aln.sequence_run_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("sequence run {}", aln.sequence_run_id))))?;
        Ok(run.biosample_guid)
    }

    /// Record (upsert) a source's haplogroup call for donor-level reconciliation.
    pub async fn record_haplogroup_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        call: &RunHaplogroupCall,
    ) -> Result<(), AppError> {
        haplogroup_call::upsert(self.store.pool(), biosample_guid, dna_type, source_key, call).await?;
        Ok(())
    }

    /// Record an assignment's top candidate as a per-source call (no-op if no match).
    async fn record_call(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        source_label: String,
        assignment: &HaploAssignment,
    ) -> Result<(), AppError> {
        if let Some(top) = assignment.ranked.first() {
            let call = RunHaplogroupCall {
                source_label,
                haplogroup: top.name.clone(),
                lineage: top.lineage.clone(),
                score: top.score,
                matched: top.matched as i64,
                expected: top.expected as i64,
            };
            self.record_haplogroup_call(biosample_guid, dna_type, source_key, &call).await?;
        }
        Ok(())
    }

    /// The reconciled donor-level haplogroup consensus across all recorded sources.
    pub async fn haplogroup_consensus(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<Consensus>, AppError> {
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;
        Ok(reconciliation::reconcile(&calls))
    }

    /// All recorded per-source calls for a subject + DNA type (for display / audit).
    pub async fn haplogroup_calls(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Vec<RunHaplogroupCall>, AppError> {
        Ok(haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?)
    }

    /// Like [`assign_mtdna_haplogroup`](Self::assign_mtdna_haplogroup) but with the tree
    /// JSON supplied directly (no network) — the testable core.
    pub async fn assign_mtdna_haplogroup_with_tree(
        &self,
        mtdna_id: i64,
        tree_json: &str,
    ) -> Result<HaploAssignment, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;

        // Sample base at each (rCRS-coordinate) position, straight from the full sequence.
        let mut calls: HashMap<i64, char> = HashMap::new();
        for (i, b) in seq.sequence.bytes().enumerate() {
            let u = b.to_ascii_uppercase();
            if matches!(u, b'A' | b'C' | b'G' | b'T') {
                calls.insert((i + 1) as i64, u as char);
            }
        }

        let tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// FTDNA mt-DNA haplotree JSON, from the on-disk cache or freshly downloaded + cached.
    async fn fetch_ftdna_mt_tree(&self) -> Result<String, AppError> {
        self.fetch_tree("https://www.familytreedna.com/public/mt-dna-haplotree/get", "ftdna-mttree.json")
            .await
    }

    /// FTDNA Y-DNA haplotree JSON, from the on-disk cache or freshly downloaded + cached.
    async fn fetch_ftdna_y_tree(&self) -> Result<String, AppError> {
        self.fetch_tree("https://www.familytreedna.com/public/y-dna-haplotree/get", "ftdna-ytree.json")
            .await
    }

    /// A cached-or-downloaded haplotree JSON (cache hit short-circuits the network).
    async fn fetch_tree(&self, url: &str, cache_file: &str) -> Result<String, AppError> {
        let path = tree_cache_path(cache_file);
        if let Ok(cached) = std::fs::read_to_string(&path) {
            if !cached.trim().is_empty() {
                return Ok(cached);
            }
        }
        let body = self
            .auth
            .http
            .get(url)
            .send()
            .await
            .and_then(|r| r.error_for_status())
            .map_err(|e| AppError::Import(format!("downloading {url}: {e}")))?
            .text()
            .await
            .map_err(|e| AppError::Import(format!("reading {url}: {e}")))?;
        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        let _ = std::fs::write(&path, &body);
        Ok(body)
    }

    /// Assign an mtDNA haplogroup directly from an alignment's chrM reads (FTDNA mt tree),
    /// the BAM-based counterpart to [`assign_mtdna_haplogroup`]. Requires a GRCh38/rCRS
    /// chrM (the tree is in rCRS coordinates).
    pub async fn assign_mtdna_haplogroup_from_alignment(
        &self,
        alignment_id: i64,
    ) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let assignment = self.assign_haplogroup_from_alignment(alignment_id, "chrM", &tree_json).await?;
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            self.record_call(bio, DnaType::Mt, &format!("aln:{alignment_id}:mt"), format!("aln #{alignment_id} mtDNA"), &assignment).await?;
        }
        Ok(assignment)
    }

    /// Assign a Y haplogroup to an alignment: fetch (and cache) the FTDNA Y-DNA haplotree,
    /// call the sample's base at each tree position on chrY, and rank by Kulczynski. Best
    /// first. Requires the alignment to have a recorded BAM/CRAM path.
    pub async fn assign_y_haplogroup(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_y_tree().await?;
        let assignment = self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json).await?;
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            self.record_call(bio, DnaType::Y, &format!("aln:{alignment_id}"), format!("aln #{alignment_id} Y"), &assignment).await?;
        }
        Ok(assignment)
    }

    /// Genotype an alignment at a haplotree's positions on `contig` and rank haplogroups by
    /// the Kulczynski measure. The networkless core shared by [`assign_y_haplogroup`] (also
    /// directly testable with a local tree + contig).
    pub async fn assign_haplogroup_from_alignment(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<HaploAssignment, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?;

        let tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        let targets: HashSet<i64> =
            tree.nodes.values().flat_map(|n| n.loci.iter().map(|l| l.position)).collect();

        let bam = PathBuf::from(bam);
        let contig_owned = contig.to_string();
        let calls = tokio::task::spawn_blocking(move || {
            let params = adaptive_haploid_params(&bam); // HiFi -> lower min_depth
            caller::call_bases_at(&bam, &contig_owned, &targets, &params)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;

        Ok(assemble_assignment(&tree, &calls))
    }

    /// Self-referential callable intervals (BED 0-based half-open) for `contig` from the
    /// alignment's own reads. Parameters adapt to the sample: long reads (HiFi) earn
    /// callability at lower depth, and the CALLABLE-run gate scales with molecule length
    /// (`f`·fragment), so long molecules clear it over far more of chrY. Requires the BAM.
    pub async fn callable_chr_intervals(&self, alignment_id: i64, contig: &str) -> Result<Vec<(i64, i64)>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let contig = contig.to_string();
        tokio::task::spawn_blocking(move || {
            let (read_len, frag_len) = coverage::estimate_molecule_lengths(&bam)?;
            let molecule = frag_len.max(read_len);
            let mut params = CallableLociParams::default();
            // Long, accurate reads (HiFi) are callable at well under half the short-read depth.
            if read_len > 1000.0 {
                params.min_depth = (params.min_depth / 2).max(2);
            }
            let min_run_len = molecule.round().max(1.0) as u32; // f = 1.0
            coverage::callable_intervals(&bam, &contig, &params, min_run_len)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?
        .map_err(Into::into)
    }

    /// The **private bucket**: de-novo SNP calls on chrY that the Y placement doesn't
    /// explain (not on the assigned backbone), classified as off-path-known (a finer/
    /// sibling FTDNA branch) or novel (a new-branch candidate). With `callable_bed` (e.g.
    /// the Poznik/1KG `b38_sites.bed`), calls outside reliable regions are dropped.
    pub async fn private_y_variants(
        &self,
        alignment_id: i64,
        callable_bed: Option<&Path>,
    ) -> Result<PrivateBucket, AppError> {
        let mask = match callable_bed {
            Some(p) => Some(navigator_analysis::mask::RegionMask::from_bed(p, "chrY")?),
            None => None,
        };
        self.private_y_core(alignment_id, mask).await
    }

    /// [`private_y_variants`] using the sample's **own** callable-Y BED as the mask
    /// (self-referential — adapts to the sample's depth and read tech; no external file).
    pub async fn private_y_variants_self_masked(&self, alignment_id: i64) -> Result<PrivateBucket, AppError> {
        let intervals = self.callable_chr_intervals(alignment_id, "chrY").await?;
        let mask = navigator_analysis::mask::RegionMask::from_intervals(intervals);
        self.private_y_core(alignment_id, Some(mask)).await
    }

    /// Shared core: assign Y, de-novo chrY, subtract the backbone, optionally mask, classify.
    async fn private_y_core(
        &self,
        alignment_id: i64,
        mask: Option<navigator_analysis::mask::RegionMask>,
    ) -> Result<PrivateBucket, AppError> {
        let tree_json = self.fetch_ftdna_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;

        let assignment = self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json).await?;
        let terminal = assignment
            .ranked
            .first()
            .ok_or_else(|| AppError::Import("no Y haplogroup match".into()))?;
        let path = navigator_analysis::haplo::path_positions(&tree, terminal.id);
        let known = navigator_analysis::haplo::tree_positions(&tree);

        // De-novo chrY (cached as an artifact), then keep only off-backbone, callable calls.
        let denovo = self.run_denovo_for_alignment(alignment_id, "chrY".to_string()).await?;
        let mut variants: Vec<PrivateVariant> = denovo
            .iter()
            .filter(|c| !path.contains(&c.position))
            .filter(|c| mask.as_ref().map_or(true, |m| m.contains(c.position)))
            .map(|c| PrivateVariant {
                position: c.position,
                reference: c.reference_allele,
                alternate: c.alternate_allele,
                depth: c.depth,
                allele_fraction: c.allele_fraction,
                class: match known.get(&c.position) {
                    Some(name) => PrivateClass::OffPathKnown(name.clone()),
                    None => PrivateClass::Novel,
                },
            })
            .collect();
        variants.sort_by_key(|v| v.position);
        Ok(PrivateBucket { terminal: terminal.name.clone(), variants })
    }

    // ---- unified import ----------------------------------------------------

    /// Detect a file's type and route it to the right subject importer (STR / variants /
    /// chip / mtDNA), using sensible defaults. Returns the detected type. Alignment files
    /// are rejected here — they attach to a sequencing test, not directly to a subject.
    pub async fn add_data(&self, biosample_guid: SampleGuid, path: &Path) -> Result<DetectedData, AppError> {
        let name = path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
        let lower = name.to_ascii_lowercase();
        // Binary/structured formats are detected by extension; only text needs a sniff.
        let by_ext = lower.ends_with(".bam")
            || lower.ends_with(".cram")
            || lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || [".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz"].iter().any(|e| lower.ends_with(e));
        let head = if by_ext { String::new() } else { read_head(path)? };
        let detected = filetype::detect(&name, &head);

        match detected {
            DetectedData::Variants => {
                self.import_variants_from_file(biosample_guid, path).await?;
            }
            DetectedData::StrProfile => {
                self.import_str_profile_from_csv(biosample_guid, "CUSTOM", None, Some("IMPORTED".into()), path).await?;
            }
            DetectedData::ChipData => {
                self.import_chip_profile_from_csv(biosample_guid, None, None, path).await?;
            }
            DetectedData::MtdnaFasta => {
                self.import_mtdna_from_fasta(biosample_guid, path).await?;
            }
            DetectedData::Alignment => {
                return Err(AppError::Import(
                    "alignment files attach to a sequencing test — add it under that test".into(),
                ));
            }
            DetectedData::Unknown => {
                return Err(AppError::Import(format!("could not recognize the data in {name}")));
            }
        }
        Ok(detected)
    }

    pub async fn panel_site_count(&self, panel_id: i64) -> Result<i64, AppError> {
        Ok(panel::site_count(self.store.pool(), panel_id).await?)
    }

    /// Genotype an alignment against a panel at the given ploidy and persist the dosages
    /// (one artifact per alignment+panel+ploidy). Runs the blocking caller off-thread.
    pub async fn genotype_panel(
        &self,
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<Vec<SiteGenotype>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?;
        let sites: Vec<Site> = panel::sites(self.store.pool(), panel_id)
            .await?
            .into_iter()
            .map(|s| Site {
                name: s.name,
                contig: s.chrom,
                position: s.position,
                reference_allele: s.reference_allele,
                alternate_allele: s.alternate_allele,
            })
            .collect();

        let bam_pb = PathBuf::from(bam);
        let params = HaploidCallerParams::default();
        let genotypes = tokio::task::spawn_blocking(move || {
            let contigs: std::collections::BTreeSet<&str> = sites.iter().map(|s| s.contig.as_str()).collect();
            let mut all = Vec::new();
            for contig in contigs {
                all.extend(caller::genotype_sites(&bam_pb, contig, &sites, ploidy, &params)?);
            }
            Ok::<_, navigator_analysis::AnalysisError>(all)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;

        self.save_analysis(alignment_id, &panel_kind(panel_id, ploidy), caller::GENOTYPE_VERSION, &genotypes).await?;
        Ok(genotypes)
    }

    /// Cached panel genotypes for an alignment, if present.
    pub async fn cached_panel_genotypes(
        &self,
        alignment_id: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<Option<Vec<SiteGenotype>>, AppError> {
        self.load_analysis(alignment_id, &panel_kind(panel_id, ploidy), caller::GENOTYPE_VERSION).await
    }

    /// Compare two alignments for IBD, using each one's cached panel genotypes. Both must
    /// have been genotyped against `panel_id` at `ploidy` first.
    pub async fn compare_ibd(
        &self,
        alignment_a: i64,
        alignment_b: i64,
        panel_id: i64,
        ploidy: u8,
        config: IbdDetectorConfig,
    ) -> Result<IbdComparison, AppError> {
        let ga = self
            .cached_panel_genotypes(alignment_a, panel_id, ploidy)
            .await?
            .ok_or_else(|| AppError::NotGenotyped(alignment_a))?;
        let gb = self
            .cached_panel_genotypes(alignment_b, panel_id, ploidy)
            .await?
            .ok_or_else(|| AppError::NotGenotyped(alignment_b))?;

        let sample_a = group_chrom_genotypes(&ga);
        let sample_b = group_chrom_genotypes(&gb);

        // Uniform 1 cM/Mb map over the chromosomes present (max observed position as the
        // length). A real HapMap genetic map can replace this later.
        let mut lengths: BTreeMap<String, i32> = BTreeMap::new();
        for sample in [&sample_a, &sample_b] {
            for (chr, cg) in sample {
                let m = cg.positions.last().copied().unwrap_or(1);
                lengths.entry(chr.clone()).and_modify(|e| *e = (*e).max(m)).or_insert(m);
            }
        }
        let pairs: Vec<(&str, i32)> = lengths.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        let gmap = GeneticMap::uniform(1.0, &pairs);

        let segments = PairwiseIbdDetector::new(config).detect_segments(&sample_a, &sample_b, &gmap);
        let summary = MatchSummary::from_segments(&segments);
        Ok(IbdComparison { summary, segments })
    }

    // ---- queries -----------------------------------------------------------

    /// Biosamples belonging to a project.
    pub async fn list_biosamples(&self, project_id: i64) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_for_project(self.store.pool(), project_id).await?)
    }

    /// Every biosample (subject), regardless of project association.
    pub async fn list_all_biosamples(&self) -> Result<Vec<Biosample>, AppError> {
        Ok(biosample::list_all(self.store.pool()).await?)
    }

    /// Sequence runs for a biosample.
    pub async fn list_sequence_runs(&self, biosample_guid: SampleGuid) -> Result<Vec<SequenceRun>, AppError> {
        Ok(sequence_run::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Alignments for a sequence run.
    pub async fn list_alignments(&self, sequence_run_id: i64) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_for_run(self.store.pool(), sequence_run_id).await?)
    }

    /// Every alignment in the workspace (for cross-sample selection like IBD compare).
    pub async fn list_all_alignments(&self) -> Result<Vec<Alignment>, AppError> {
        Ok(alignment::list_all(self.store.pool()).await?)
    }

    /// Projects with their sample counts, for a dashboard/list view.
    pub async fn project_overview(&self) -> Result<Vec<ProjectOverview>, AppError> {
        let mut out = Vec::new();
        for project in project::list(self.store.pool()).await? {
            let sample_count = biosample::count_for_project(self.store.pool(), project.id).await?;
            out.push(ProjectOverview { project, sample_count });
        }
        Ok(out)
    }
}

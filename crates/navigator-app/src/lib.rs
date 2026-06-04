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
use navigator_analysis::heteroplasmy::{self, HeteroplasmyParams};
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
pub use navigator_analysis::heteroplasmy::HeteroplasmySite;
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
use navigator_refgenome::{cache as refgenome_cache, canonical_build, Build as ReferenceBuild, LiftedPos, ReferenceGateway};
pub use navigator_refgenome::RefStatus;
use navigator_sync::{
    AuditEntryRecord, HaplogroupReconciliationRecord, HeteroplasmyObservationRecord,
    IdentityVerificationRecord, ManualOverrideRecord, ReconciliationStatusRecord,
    RunHaplogroupCallRecord, HAPLOGROUP_RECONCILIATION_COLLECTION,
};

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
pub use navigator_domain::reconciliation::{
    AuditEntry, CompatibilityLevel, Consensus, DnaType, IdentityVerification, ReconciledVariant, VariantStatus,
    VerificationStatus,
};
use navigator_domain::strprofile::{self, NewStrProfile, StrProfile};
use navigator_domain::variants::{self, NewVariantSet, VariantSet};
pub use navigator_domain::variants::SourceType;
use navigator_store::{
    alignment, artifact, biosample, chip_profile, haplogroup_call, mtdna as mtdna_store, project,
    reconciliation as recon_store, sequence_run, str_profile, variant_set, Store, StoreError,
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
fn adaptive_haploid_params(bam_path: &Path, reference: Option<&Path>) -> HaploidCallerParams {
    let mut params = HaploidCallerParams::default();
    if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(bam_path, reference) {
        if read_len > 1000.0 {
            params.min_depth = (params.min_depth / 2).max(2);
        }
    }
    params
}

/// The lexicon's UPPER_SNAKE compatibility level (matches the AppView's knownValues).
fn compat_lexicon(c: CompatibilityLevel) -> &'static str {
    match c {
        CompatibilityLevel::Compatible => "COMPATIBLE",
        CompatibilityLevel::MinorDivergence => "MINOR_DIVERGENCE",
        CompatibilityLevel::MajorDivergence => "MAJOR_DIVERGENCE",
        CompatibilityLevel::Incompatible => "INCOMPATIBLE",
    }
}

/// The lexicon's DNA-type token for the reconciliation record (`Y_DNA`/`MT_DNA`).
fn dna_type_lexicon(d: DnaType) -> &'static str {
    match d {
        DnaType::Y => "Y_DNA",
        DnaType::Mt => "MT_DNA",
    }
}

/// The lexicon's UPPER_SNAKE verification status.
fn verification_lexicon(s: VerificationStatus) -> &'static str {
    match s {
        VerificationStatus::VerifiedSame => "VERIFIED_SAME",
        VerificationStatus::LikelySame => "LIKELY_SAME",
        VerificationStatus::Uncertain => "UNCERTAIN",
        VerificationStatus::LikelyDifferent => "LIKELY_DIFFERENT",
        VerificationStatus::VerifiedDifferent => "VERIFIED_DIFFERENT",
    }
}

/// Reference build inferred from an alignment filename (`*.chm13.*` → CHM13v2.0, else
/// unknown). A best-effort label; the actual decode uses the supplied reference FASTA.
fn reference_build_for(path: &Path) -> String {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("").to_ascii_lowercase();
    if name.contains("chm13") {
        "chm13v2.0".to_string()
    } else {
        "unknown".to_string()
    }
}

/// Watson–Crick complement of a base (for reverse-strand lifts); non-ACGT passes through.
fn complement_base(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// The build a haplotree's positions are in, by contig: the FTDNA Y tree is GRCh38; mtDNA
/// (`chrM`) is rCRS and stays a direct query (no chain), so it returns `None`.
fn tree_build_for_contig(contig: &str) -> Option<&'static str> {
    if contig.eq_ignore_ascii_case("chrY") {
        Some("GRCh38")
    } else {
        None
    }
}

/// Whether `<alignment>.crai`/`.bai` is present among the discovered index files.
fn has_sibling_index(aln_path: &Path, index_files: &[PathBuf]) -> bool {
    let Some(aln_name) = aln_path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    index_files.iter().filter_map(|i| i.file_name().and_then(|n| n.to_str())).any(|n| {
        n == format!("{aln_name}.crai") || n == format!("{aln_name}.bai")
    })
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

/// Autosomal genotype concordance between two genotyped alignments: (matched, compared)
/// over sites both called (dosage within ploidy). ~1.0 ⇒ same individual; relatives lower.
fn genotype_concordance(a: &[SiteGenotype], b: &[SiteGenotype]) -> (i64, i64) {
    let called = |g: &SiteGenotype| (0..=g.ploidy as i32).contains(&g.dosage);
    let idx: HashMap<(&str, i64), i32> =
        b.iter().filter(|g| called(g)).map(|g| ((g.contig.as_str(), g.position), g.dosage)).collect();
    let (mut matched, mut sites) = (0i64, 0i64);
    for g in a.iter().filter(|g| called(g)) {
        if let Some(&db) = idx.get(&(g.contig.as_str(), g.position)) {
            sites += 1;
            if db == g.dosage {
                matched += 1;
            }
        }
    }
    (matched, sites)
}

/// A project plus a rolled-up count for list/dashboard views.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProjectOverview {
    pub project: Project,
    pub sample_count: i64,
}

/// One row of a project's per-sample report: coverage roll-up + haplogroup consensus.
/// Coverage fields are `None` when no coverage has been computed yet; haplogroup fields
/// are `None` until calls are recorded (deferred this slice).
#[derive(Debug, Clone)]
pub struct ProjectSampleReport {
    pub biosample: Biosample,
    /// An alignment to drive "recompute coverage" from (the coverage-bearing one if any,
    /// else the first); `None` if the sample has no alignments.
    pub primary_alignment_id: Option<i64>,
    pub alignment_count: usize,
    pub mean_coverage: Option<f64>,
    pub median_coverage: Option<f64>,
    pub pct_10x: Option<f64>,
    pub pct_20x: Option<f64>,
    pub callable_bases: Option<u64>,
    pub y_haplogroup: Option<String>,
    pub mt_haplogroup: Option<String>,
}

/// A reference build an import needs but doesn't have cached — surfaced so the UI can
/// prompt and download it before retrying.
#[derive(Debug, Clone)]
pub struct BuildNeed {
    pub build: String,
    pub url: String,
    pub est_bytes: u64,
}

/// Outcome of a project-wide analyze pass (coverage + Y haplogroup per sample).
#[derive(Debug, Clone)]
pub struct AnalyzeSummary {
    pub project_id: i64,
    pub samples: usize,
    pub coverage_done: usize,
    pub y_done: usize,
    /// Per-sample failures (best-effort: one sample's error doesn't abort the rest).
    pub errors: Vec<String>,
}

/// Outcome of a batch project-directory import (idempotent — counts only what's new).
#[derive(Debug, Clone)]
pub struct ProjectImportSummary {
    pub project: Project,
    pub samples_total: usize,
    pub samples_created: usize,
    pub alignments_created: usize,
    pub alignments_skipped: usize,
    /// Sample ids whose alignment had no sibling index (.crai/.bai) — coverage needs one.
    pub missing_index: Vec<String>,
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
    gateway: ReferenceGateway,
}

impl App {
    pub fn new(store: Store) -> Self {
        let gateway = ReferenceGateway::new(refgenome_cache::base_dir(), dev_http_client());
        App { store, auth: Auth::new(), gateway }
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
        let reference = PathBuf::from(reference);
        let probe = bam.clone();
        let probe_ref = reference.clone();
        let params = tokio::task::spawn_blocking(move || adaptive_haploid_params(&probe, Some(&probe_ref)))
            .await
            .map_err(|e| AppError::Join(e.to_string()))?; // HiFi -> lower min_depth
        self.run_denovo_caller(alignment_id, bam, reference, contig, params).await
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
    /// table (a YSEQ/Sanger panel export fits this). Indels/symbolic alleles are dropped
    /// (SNP-only). `source_type` sets the concordance weight (Sanger = gold standard).
    pub async fn import_variants_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        source_type: SourceType,
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
        let new = NewVariantSet { biosample_guid, source_label: label, source_type, calls };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// Add a manually-entered variant set — paste `contig,position,ref,alt` rows (e.g.
    /// Sanger/YSEQ confirmations). `source_type` sets the weight (Sanger = 1.0).
    pub async fn add_variants(
        &self,
        biosample_guid: SampleGuid,
        source_label: &str,
        source_type: SourceType,
        text: &str,
    ) -> Result<VariantSet, AppError> {
        let calls = variants::parse_csv(text).map_err(AppError::Import)?;
        let new = NewVariantSet {
            biosample_guid,
            source_label: source_label.to_string(),
            source_type,
            calls,
        };
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
        let new = NewVariantSet {
            biosample_guid: seq.biosample_guid,
            source_label: label,
            source_type: variants::SourceType::Imported,
            calls,
        };
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
        self.audit(biosample_guid, dna_type, "RUN_RECORDED", &format!("{source_key}: {}", call.haplogroup)).await?;
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

    /// The reconciled donor-level haplogroup consensus across all recorded sources. A user
    /// manual override, when set, replaces the computed terminal (flagged `overridden`).
    pub async fn haplogroup_consensus(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<Consensus>, AppError> {
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;
        let mut consensus = reconciliation::reconcile(&calls);

        if let Some((hg, reason)) = recon_store::get_override(self.store.pool(), biosample_guid, dna_type).await? {
            let mut c = consensus.unwrap_or(Consensus {
                haplogroup: hg.clone(),
                lineage: vec![hg.clone()],
                compatibility: CompatibilityLevel::Compatible,
                divergence_point: None,
                confidence: 1.0,
                run_count: 0,
                overridden: true,
                warnings: Vec::new(),
            });
            c.haplogroup = hg;
            c.overridden = true;
            c.confidence = 1.0;
            c.warnings.push(match reason {
                Some(r) => format!("manual override: {r}"),
                None => "manual override".to_string(),
            });
            consensus = Some(c);
        }
        Ok(consensus)
    }

    /// Manually override the consensus haplogroup for a subject + DNA type.
    pub async fn set_manual_override(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        haplogroup: &str,
        reason: Option<&str>,
    ) -> Result<(), AppError> {
        recon_store::set_override(self.store.pool(), biosample_guid, dna_type, haplogroup, reason).await?;
        self.audit(biosample_guid, dna_type, "MANUAL_OVERRIDE", &format!("override to {haplogroup}")).await
    }

    /// Clear a manual override.
    pub async fn clear_manual_override(&self, biosample_guid: SampleGuid, dna_type: DnaType) -> Result<(), AppError> {
        recon_store::clear_override(self.store.pool(), biosample_guid, dna_type).await?;
        self.audit(biosample_guid, dna_type, "OVERRIDE_CLEARED", "cleared manual override").await
    }

    /// The reconciliation audit log for a subject + DNA type.
    pub async fn reconciliation_audit(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Vec<AuditEntry>, AppError> {
        Ok(recon_store::list_audit(self.store.pool(), biosample_guid, dna_type).await?)
    }

    async fn audit(&self, biosample_guid: SampleGuid, dna_type: DnaType, action: &str, note: &str) -> Result<(), AppError> {
        let entry = AuditEntry { timestamp: Utc::now().to_rfc3339(), action: action.to_string(), note: note.to_string() };
        recon_store::append_audit(self.store.pool(), biosample_guid, dna_type, &entry).await?;
        Ok(())
    }

    /// Reconcile the subject's variant sets at the variant level — which positions are
    /// confirmed across sources, in conflict, or single-source (Sanger-confirmation
    /// candidates).
    pub async fn reconcile_variants(&self, biosample_guid: SampleGuid) -> Result<Vec<ReconciledVariant>, AppError> {
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let sources: Vec<(String, f64, &[variants::VariantCall])> = sets
            .iter()
            .map(|s| (s.source_label.clone(), s.source_type.snp_weight(), s.calls.as_slice()))
            .collect();
        Ok(reconciliation::reconcile_variants(&sources))
    }

    /// Build the `com.decodingus.atmosphere.haplogroupReconciliation` record JSON for a
    /// subject + DNA type from the stored consensus, per-run calls, manual override, and
    /// audit log. mtDNA heteroplasmy observations and an optional identity-verification
    /// result are passed in (the caller computes them from the relevant alignments).
    async fn reconciliation_record(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<serde_json::Value, AppError> {
        let consensus = self.haplogroup_consensus(biosample_guid, dna_type).await?.ok_or_else(|| {
            AppError::Store(StoreError::NotFound(format!(
                "no {} haplogroup calls for {}",
                dna_type.as_str(),
                biosample_guid.0
            )))
        })?;
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;

        let run_calls = calls
            .iter()
            .map(|c| RunHaplogroupCallRecord {
                source_ref: c.source_label.clone(),
                haplogroup: c.haplogroup.clone(),
                confidence: c.score.to_string(),
                call_method: "SNP_PHYLOGENETIC".into(),
                score: Some(c.score.to_string()),
                supporting_snps: Some(c.matched),
                conflicting_snps: Some((c.expected - c.matched).max(0)),
            })
            .collect();

        let status = ReconciliationStatusRecord {
            compatibility_level: compat_lexicon(consensus.compatibility).into(),
            consensus_haplogroup: consensus.haplogroup.clone(),
            confidence: Some(consensus.confidence.to_string()),
            divergence_point: consensus.divergence_point.clone(),
            branch_compatibility_score: None,
            snp_concordance: identity.and_then(|i| i.snp_concordance).map(|f| f.to_string()),
            run_count: consensus.run_count as i64,
            warnings: consensus.warnings.clone(),
        };

        // Heteroplasmy is mtDNA-only; major frequency is 1 − minor fraction.
        let heteroplasmy_observations = if dna_type == DnaType::Mt {
            heteroplasmy
                .iter()
                .map(|h| HeteroplasmyObservationRecord {
                    position: h.position,
                    major_allele: h.major_base.to_string(),
                    minor_allele: h.minor_base.to_string(),
                    major_allele_frequency: (1.0 - h.minor_fraction).to_string(),
                    depth: Some(h.depth as i64),
                    is_defining_snp: None,
                    affected_haplogroup: None,
                })
                .collect()
        } else {
            Vec::new()
        };

        let identity_verification = identity.map(|i| IdentityVerificationRecord {
            kinship_coefficient: None,
            fingerprint_snp_concordance: i.snp_concordance.map(|f| f.to_string()),
            y_str_distance: i.y_str_distance,
            verification_status: Some(verification_lexicon(i.status).into()),
            verification_method: Some(i.method.clone()),
        });

        let manual_override = recon_store::get_override(self.store.pool(), biosample_guid, dna_type)
            .await?
            .map(|(hg, reason)| ManualOverrideRecord {
                overridden_haplogroup: hg,
                reason,
                overridden_at: Utc::now().to_rfc3339(),
                overridden_by: self.current_account(),
            });

        let audit_log = self
            .reconciliation_audit(biosample_guid, dna_type)
            .await?
            .into_iter()
            .map(|a| AuditEntryRecord {
                timestamp: a.timestamp,
                action: a.action,
                previous_consensus: None,
                new_consensus: None,
                run_ref: None,
                notes: Some(a.note),
            })
            .collect();

        let record = HaplogroupReconciliationRecord::new(
            biosample_guid.0.to_string(),
            dna_type_lexicon(dna_type),
            Utc::now().to_rfc3339(),
            status,
            run_calls,
            heteroplasmy_observations,
            identity_verification,
            manual_override,
            audit_log,
        );
        Ok(serde_json::to_value(&record)?)
    }

    /// Publish a subject's haplogroup reconciliation using an explicit `client` (the
    /// testable core; production callers use [`publish_reconciliation`](Self::publish_reconciliation)).
    pub async fn publish_reconciliation_with(
        &self,
        client: &PdsClient,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<RecordRef, AppError> {
        let value = self.reconciliation_record(biosample_guid, dna_type, heteroplasmy, identity).await?;
        Ok(client.create_record(HAPLOGROUP_RECONCILIATION_COLLECTION, value, None).await?)
    }

    /// Publish a subject's haplogroup reconciliation record to the signed-in account's PDS
    /// (with refresh-on-expiry and retry/backoff via [`AsyncSync`]).
    pub async fn publish_reconciliation(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        heteroplasmy: &[HeteroplasmySite],
        identity: Option<&IdentityVerification>,
    ) -> Result<RecordRef, AppError> {
        let mut engine = self.sync_engine()?; // auth check before touching the DB
        let value = self.reconciliation_record(biosample_guid, dna_type, heteroplasmy, identity).await?;
        Ok(engine.push_create(HAPLOGROUP_RECONCILIATION_COLLECTION, value).await?)
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

    /// Scan an alignment's chrM pileup for heteroplasmic positions — sites where a second
    /// mitochondrial allele coexists above the noise floor. A screening pass for the
    /// reconciliation view (a curator judges real heteroplasmy vs. artefacts); ascending
    /// by position. Requires a chrM-bearing BAM.
    pub async fn mtdna_heteroplasmy(&self, alignment_id: i64) -> Result<Vec<HeteroplasmySite>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = aln.reference_path.map(PathBuf::from);
        tokio::task::spawn_blocking(move || {
            heteroplasmy::detect_heteroplasmy(&bam, "chrM", &HeteroplasmyParams::default(), reference.as_deref())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?
        .map_err(Into::into)
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
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = aln.reference_path.map(PathBuf::from);

        let tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        let targets: HashSet<i64> =
            tree.nodes.values().flat_map(|n| n.loci.iter().map(|l| l.position)).collect();

        // The tree's positions are in its own build (Y → GRCh38, mt → rCRS). If the alignment
        // is a different build, lift the positions onto it, query there, and map back;
        // otherwise query directly (Y-on-GRCh38, mt-on-GRCh38/rCRS).
        let lifted = self
            .lifted_targets(&aln.reference_build, reference.as_deref(), contig, &targets)
            .await?;

        let calls = match lifted {
            Some(lifted) => self.build_calls_from_lifted(&bam, reference.as_deref(), lifted).await?,
            None => {
                let bam = bam.clone();
                let reference = reference.clone();
                let targets = targets.clone();
                let contig_owned = contig.to_string();
                tokio::task::spawn_blocking(move || {
                    let params = adaptive_haploid_params(&bam, reference.as_deref()); // HiFi -> lower min_depth
                    caller::call_bases_at(&bam, &contig_owned, &targets, &params, reference.as_deref())
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))??
            }
        };

        Ok(assemble_assignment(&tree, &calls))
    }

    /// Lift the haplotree's positions onto the alignment's build, or `None` to query the tree
    /// positions directly. **chrY**: uses the (auto-downloaded) GRCh38→build liftover chain.
    /// **chrM**: a self-generated rCRS↔`chrM` map — bundled rCRS aligned to *this* reference's
    /// `chrM` (CHM13 builds only; GRCh38/rCRS `chrM` is already rCRS → direct).
    async fn lifted_targets(
        &self,
        reference_build: &str,
        reference: Option<&Path>,
        contig: &str,
        targets: &HashSet<i64>,
    ) -> Result<Option<Vec<LiftedPos>>, AppError> {
        if targets.is_empty() {
            return Ok(None);
        }

        // chrY: downloaded nuclear chain.
        if let Some(src) = tree_build_for_contig(contig) {
            let differ = matches!((canonical_build(src), canonical_build(reference_build)), (Some(s), Some(t)) if s != t);
            if differ && self.gateway.chain_available(src, reference_build) {
                self.gateway.resolve_chain(src, reference_build, &mut |_, _| {}).await?;
                let targets_vec: Vec<i64> = targets.iter().copied().collect();
                return Ok(Some(self.gateway.lift_positions(src, reference_build, contig, &targets_vec)?));
            }
            return Ok(None);
        }

        // chrM on CHM13: self-generated rCRS↔chrM alignment map (no chain exists).
        if contig.eq_ignore_ascii_case("chrM") && canonical_build(reference_build) == Some(ReferenceBuild::Chm13v2) {
            let Some(reference) = reference else { return Ok(None) };
            let reference = reference.to_path_buf();
            // Align bundled rCRS to this reference's chrM (cheap, ~16.5 kb) → (rcrs, chrM) pairs.
            let map = tokio::task::spawn_blocking(move || {
                navigator_analysis::reader::read_contig_sequence(&reference, "chrM").map(|chrm| {
                    let chrm = String::from_utf8_lossy(&chrm).into_owned();
                    navigator_analysis::mtvariants::align_positions(navigator_analysis::mtvariants::rcrs(), &chrm)
                })
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))?;
            let Ok(pairs) = map else { return Ok(None) }; // chrM absent/unreadable → direct fallback
            // rcrs_idx/chrm_idx are 0-based; tree + query positions are 1-based.
            let by_rcrs: HashMap<i64, i64> = pairs.into_iter().map(|(r, c)| (r as i64 + 1, c as i64 + 1)).collect();
            let lifted = targets
                .iter()
                .filter_map(|&t| by_rcrs.get(&t).map(|&c| LiftedPos { tree_pos: t, contig: "chrM".to_string(), pos: c, reverse: false }))
                .collect();
            return Ok(Some(lifted));
        }

        Ok(None)
    }

    /// Query the already-lifted positions and map observed bases back to the original tree
    /// positions so [`assemble_assignment`] (which keys on tree positions) scores unchanged.
    /// Queries each lifted contig present in the BAM header; minus-strand lifts are
    /// reverse-complemented.
    async fn build_calls_from_lifted(
        &self,
        bam: &Path,
        reference: Option<&Path>,
        lifted: Vec<LiftedPos>,
    ) -> Result<HashMap<i64, char>, AppError> {
        // Group lifted positions by their target contig + a back-map (lifted → tree position,
        // plus whether the lift was to the minus strand → the base needs complementing).
        let mut by_contig: HashMap<String, HashSet<i64>> = HashMap::new();
        let mut back: HashMap<(String, i64), (i64, bool)> = HashMap::new();
        for lp in lifted {
            by_contig.entry(lp.contig.clone()).or_default().insert(lp.pos);
            back.insert((lp.contig, lp.pos), (lp.tree_pos, lp.reverse));
        }

        // Only query contigs the alignment actually has (drop off-target lifts).
        let header_contigs: HashSet<String> = {
            let bam = bam.to_path_buf();
            let reference = reference.map(|p| p.to_path_buf());
            tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, reference.as_deref()))
                .await
                .map_err(|e| AppError::Join(e.to_string()))??
                .into_iter()
                .collect()
        };

        let mut calls: HashMap<i64, char> = HashMap::new();
        for (qcontig, set) in by_contig {
            if !header_contigs.contains(&qcontig) {
                continue;
            }
            let bam = bam.to_path_buf();
            let reference = reference.map(|p| p.to_path_buf());
            let qc = qcontig.clone();
            let lifted_calls = tokio::task::spawn_blocking(move || {
                let params = adaptive_haploid_params(&bam, reference.as_deref());
                caller::call_bases_at(&bam, &qc, &set, &params, reference.as_deref())
            })
            .await
            .map_err(|e| AppError::Join(e.to_string()))??;
            for (lpos, base) in lifted_calls {
                if let Some(&(tree_pos, reverse)) = back.get(&(qcontig.clone(), lpos)) {
                    // Inverted tracts (common on the CHM13 Y): the tree allele is GRCh38-forward,
                    // so complement the base read off the minus-strand-lifted CHM13 position.
                    calls.insert(tree_pos, if reverse { complement_base(base) } else { base });
                }
            }
        }
        Ok(calls)
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
        let reference = aln.reference_path.map(PathBuf::from);
        let contig = contig.to_string();
        tokio::task::spawn_blocking(move || {
            let (read_len, frag_len) = coverage::estimate_molecule_lengths(&bam, reference.as_deref())?;
            let molecule = frag_len.max(read_len);
            let mut params = CallableLociParams::default();
            // Long, accurate reads (HiFi) are callable at well under half the short-read depth.
            if read_len > 1000.0 {
                params.min_depth = (params.min_depth / 2).max(2);
            }
            let min_run_len = molecule.round().max(1.0) as u32; // f = 1.0
            coverage::callable_intervals(&bam, &contig, &params, min_run_len, reference.as_deref())
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
                self.import_variants_from_file(biosample_guid, path, variants::SourceType::Imported).await?;
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

    /// Batch-import a NAS project directory: scan `{dir}/{sample}/…` and create the Project
    /// plus its Biosample → SequenceRun → Alignment rows. The reference is resolved per
    /// alignment: pass `Some(fasta)` to use a specific FASTA (validated with its `.fai`) for
    /// every alignment, or `None` to let the gateway resolve each file's inferred build from
    /// the cache. If a needed build isn't cached, returns [`AppError::ReferenceNeeded`]
    /// **before any DB writes** so the UI can prompt + download, then retry. Idempotent: an
    /// existing project (by name), biosample (by donor id), or alignment (by path) is reused.
    /// Coverage is NOT computed here — run it per alignment or via the project report.
    pub async fn import_project_dir(
        &self,
        dir: &Path,
        reference: Option<PathBuf>,
        administrator: String,
    ) -> Result<ProjectImportSummary, AppError> {
        // An explicit FASTA must exist and be indexed; it applies to every alignment.
        if let Some(path) = &reference {
            if !path.exists() {
                return Err(AppError::Import(format!("reference FASTA not found: {}", path.display())));
            }
            let fai = PathBuf::from(format!("{}.fai", path.display()));
            if !fai.exists() {
                return Err(AppError::Import(format!("reference FASTA index (.fai) not found: {}", fai.display())));
            }
        }

        let scan_dir = dir.to_path_buf();
        let discovered = tokio::task::spawn_blocking(move || navigator_analysis::scan::scan(&scan_dir))
            .await
            .map_err(|e| AppError::Join(e.to_string()))??;

        // Resolve each alignment's reference build to a path (explicit FASTA, else the cache).
        // Collect any builds that need downloading and bail before writing anything.
        let explicit = reference.as_ref().map(|p| p.to_string_lossy().into_owned());
        let mut resolved: HashMap<String, String> = HashMap::new();
        let mut needs: Vec<BuildNeed> = Vec::new();
        for sample in &discovered.samples {
            for aln_path in &sample.alignment_files {
                let build = reference_build_for(aln_path);
                if resolved.contains_key(&build) || needs.iter().any(|n| n.build == build) {
                    continue;
                }
                if let Some(ref path) = explicit {
                    resolved.insert(build, path.clone());
                } else if let Some(p) = self.gateway.cached_reference(&build) {
                    resolved.insert(build, p.to_string_lossy().into_owned());
                } else {
                    match self.gateway.reference_status(&build) {
                        RefStatus::NeedsDownload { url, est_bytes } => {
                            needs.push(BuildNeed { build, url, est_bytes })
                        }
                        RefStatus::Unknown => {
                            return Err(AppError::Import(format!(
                                "unknown reference build '{build}' — supply a reference FASTA explicitly"
                            )))
                        }
                        RefStatus::Cached(p) | RefStatus::LocalOverride(p) => {
                            resolved.insert(build, p.to_string_lossy().into_owned());
                        }
                    }
                }
            }
        }
        if !needs.is_empty() {
            return Err(AppError::ReferenceNeeded(needs));
        }

        // Project: reuse an existing one with the same name.
        let project = match project::list(self.store.pool()).await?.into_iter().find(|p| p.name == discovered.project_id) {
            Some(p) => p,
            None => {
                self.create_project(NewProject {
                    name: discovered.project_id.clone(),
                    description: None,
                    administrator,
                })
                .await?
            }
        };

        let mut summary = ProjectImportSummary {
            project: project.clone(),
            samples_total: discovered.samples.len(),
            samples_created: 0,
            alignments_created: 0,
            alignments_skipped: 0,
            missing_index: Vec::new(),
        };

        for sample in &discovered.samples {
            // Biosample: reuse by donor identifier within the project.
            let biosample = match biosample::list_for_project(self.store.pool(), project.id)
                .await?
                .into_iter()
                .find(|b| b.donor_identifier == sample.sample_id)
            {
                Some(b) => b,
                None => {
                    summary.samples_created += 1;
                    self.add_biosample(Some(project.id), sample.sample_id.clone(), Some(sample.sample_id.clone()), None)
                        .await?
                }
            };

            // SequenceRun: reuse the first existing run, else create one (defaults to WGS).
            let run = match sequence_run::list_for_biosample(self.store.pool(), biosample.guid).await?.into_iter().next() {
                Some(r) => r,
                None => {
                    self.record_sequence_run(NewSequenceRun {
                        biosample_guid: biosample.guid,
                        platform_name: "UNKNOWN".into(),
                        instrument_model: None,
                        test_type: "WGS".into(),
                        library_layout: None,
                        total_reads: None,
                        pf_reads_aligned: None,
                        mean_read_length: None,
                        mean_insert_size: None,
                    })
                    .await?
                }
            };

            let existing = alignment::list_for_run(self.store.pool(), run.id).await?;
            for aln_path in &sample.alignment_files {
                let path_str = aln_path.to_string_lossy().into_owned();
                if existing.iter().any(|a| a.bam_path.as_deref() == Some(path_str.as_str())) {
                    summary.alignments_skipped += 1;
                    continue;
                }
                if !has_sibling_index(aln_path, &sample.index_files) {
                    summary.missing_index.push(sample.sample_id.clone());
                }
                let build = reference_build_for(aln_path);
                let reference_path = resolved.get(&build).cloned();
                self.record_alignment(NewAlignment {
                    sequence_run_id: run.id,
                    reference_build: build,
                    aligner: "unknown".into(),
                    variant_caller: None,
                    bam_path: Some(path_str),
                    reference_path,
                })
                .await?;
                summary.alignments_created += 1;
            }
        }
        Ok(summary)
    }

    /// Cache/override status of a reference build (no network).
    pub fn reference_status(&self, build: &str) -> RefStatus {
        self.gateway.reference_status(build)
    }

    /// Resolve a reference build to a cached, indexed `.fa`, downloading on a miss.
    /// `progress(received, total)` is invoked as bytes arrive.
    pub async fn resolve_reference(
        &self,
        build: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_reference(build, progress).await?)
    }

    /// Resolve (and cache) a liftover chain for a build pair, downloading on a miss. The
    /// cached `.chain` is then available for the haplogroup/liftover path.
    pub async fn resolve_chain(
        &self,
        from: &str,
        to: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, AppError> {
        Ok(self.gateway.resolve_chain(from, to, progress).await?)
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
        let reference = aln.reference_path.map(PathBuf::from);
        let params = HaploidCallerParams::default();
        let genotypes = tokio::task::spawn_blocking(move || {
            let contigs: std::collections::BTreeSet<&str> = sites.iter().map(|s| s.contig.as_str()).collect();
            let mut all = Vec::new();
            for contig in contigs {
                all.extend(caller::genotype_sites(&bam_pb, contig, &sites, ploidy, &params, reference.as_deref())?);
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

    /// Identity verification — are two alignments the same individual? Autosomal genotype
    /// concordance at the panel sites (primary), corroborated by Y-STR distance when both
    /// subjects have an STR profile. Both alignments must be genotyped against the panel.
    pub async fn verify_identity(
        &self,
        alignment_a: i64,
        alignment_b: i64,
        panel_id: i64,
        ploidy: u8,
    ) -> Result<IdentityVerification, AppError> {
        let ga = self
            .cached_panel_genotypes(alignment_a, panel_id, ploidy)
            .await?
            .ok_or(AppError::NotGenotyped(alignment_a))?;
        let gb = self
            .cached_panel_genotypes(alignment_b, panel_id, ploidy)
            .await?
            .ok_or(AppError::NotGenotyped(alignment_b))?;
        let (matched, sites) = genotype_concordance(&ga, &gb);
        let concordance = (sites > 0).then(|| matched as f64 / sites as f64);

        // Optional Y-STR corroboration from each subject's first STR profile.
        let (mut y_dist, mut y_markers) = (None, 0i64);
        if let (Ok(ba), Ok(bb)) = (self.biosample_of_alignment(alignment_a).await, self.biosample_of_alignment(alignment_b).await) {
            let (pa, pb) = (self.list_str_profiles(ba).await?, self.list_str_profiles(bb).await?);
            if let (Some(a), Some(b)) = (pa.first(), pb.first()) {
                let (d, c) = strprofile::str_distance(&a.markers, &b.markers);
                if c > 0 {
                    y_dist = Some(d);
                    y_markers = c;
                }
            }
        }
        Ok(reconciliation::classify_identity(concordance, sites, y_dist, y_markers))
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

    /// Per-sample report for a project: each biosample's alignment count, coverage roll-up
    /// (the first alignment with cached coverage), and Y/mtDNA haplogroup consensus.
    /// Composes existing per-subject queries (no new join) — coverage/haplogroup cells are
    /// `None` until those analyses have run.
    pub async fn project_report(&self, project_id: i64) -> Result<Vec<ProjectSampleReport>, AppError> {
        let mut out = Vec::new();
        for biosample in biosample::list_for_project(self.store.pool(), project_id).await? {
            let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
            let mut coverage = None;
            let mut coverage_aln = None;
            for a in &alignments {
                if let Some(c) = self.cached_coverage(a.id).await? {
                    coverage = Some(c);
                    coverage_aln = Some(a.id);
                    break;
                }
            }
            // Prefer the coverage-bearing alignment; else fall back to the first.
            let primary_alignment_id = coverage_aln.or_else(|| alignments.first().map(|a| a.id));
            let y_haplogroup = self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.map(|c| c.haplogroup);
            let mt_haplogroup = self.haplogroup_consensus(biosample.guid, DnaType::Mt).await?.map(|c| c.haplogroup);
            out.push(ProjectSampleReport {
                primary_alignment_id,
                alignment_count: alignments.len(),
                mean_coverage: coverage.as_ref().map(|c| c.mean_coverage),
                median_coverage: coverage.as_ref().map(|c| c.median_coverage),
                pct_10x: coverage.as_ref().map(|c| c.pct_10x),
                pct_20x: coverage.as_ref().map(|c| c.pct_20x),
                callable_bases: coverage.as_ref().map(|c| c.callable_bases),
                y_haplogroup,
                mt_haplogroup,
                biosample,
            });
        }
        Ok(out)
    }

    /// Analyze every sample in a project: compute coverage and assign the Y haplogroup on each
    /// sample's primary (first BAM-bearing) alignment, so the project report fills in. Coverage
    /// already cached and Y already recorded are skipped (idempotent re-run). Best-effort: one
    /// sample's failure is recorded and the rest continue. mtDNA is intentionally not assigned
    /// here (provisional on CHM13 — see the reconciliation/liftover notes).
    pub async fn analyze_project(&self, project_id: i64) -> Result<AnalyzeSummary, AppError> {
        let mut summary = AnalyzeSummary { project_id, samples: 0, coverage_done: 0, y_done: 0, errors: Vec::new() };
        for biosample in biosample::list_for_project(self.store.pool(), project_id).await? {
            let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
            let Some(aln) = alignments.iter().find(|a| a.bam_path.is_some()) else {
                continue;
            };
            summary.samples += 1;
            let label = &biosample.donor_identifier;

            if self.cached_coverage(aln.id).await?.is_some() {
                summary.coverage_done += 1;
            } else {
                match self.run_coverage_for_alignment(aln.id).await {
                    Ok(_) => summary.coverage_done += 1,
                    Err(e) => summary.errors.push(format!("{label} coverage: {e}")),
                }
            }

            if self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.is_some() {
                summary.y_done += 1;
            } else {
                match self.assign_y_haplogroup(aln.id).await {
                    Ok(_) => summary.y_done += 1,
                    Err(e) => summary.errors.push(format!("{label} Y: {e}")),
                }
            }
        }
        Ok(summary)
    }
}

/// Render a project report as CSV (one header row + one row per sample). Empty cells for
/// not-yet-computed coverage/haplogroup. Hand-formatted to avoid a CSV dependency; values
/// containing a comma or quote are quoted.
pub fn report_csv(rows: &[ProjectSampleReport]) -> String {
    fn field(s: &str) -> String {
        if s.contains([',', '"', '\n']) {
            format!("\"{}\"", s.replace('"', "\"\""))
        } else {
            s.to_string()
        }
    }
    fn num(o: Option<f64>) -> String {
        o.map(|v| format!("{v:.4}")).unwrap_or_default()
    }

    let mut s = String::from(
        "sample_id,alignment_count,mean_coverage,median_coverage,pct_10x,pct_20x,callable_bases,y_haplogroup,mt_haplogroup\n",
    );
    for r in rows {
        s.push_str(&field(&r.biosample.donor_identifier));
        s.push(',');
        s.push_str(&r.alignment_count.to_string());
        s.push(',');
        s.push_str(&num(r.mean_coverage));
        s.push(',');
        s.push_str(&num(r.median_coverage));
        s.push(',');
        s.push_str(&num(r.pct_10x));
        s.push(',');
        s.push_str(&num(r.pct_20x));
        s.push(',');
        s.push_str(&r.callable_bases.map(|v| v.to_string()).unwrap_or_default());
        s.push(',');
        s.push_str(&field(r.y_haplogroup.as_deref().unwrap_or("")));
        s.push(',');
        s.push_str(&field(r.mt_haplogroup.as_deref().unwrap_or("")));
        s.push('\n');
    }
    s
}

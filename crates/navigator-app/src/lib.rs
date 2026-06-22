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
use navigator_analysis::ancestry::{self as ancestry_analysis};
use navigator_analysis::caller::{self, HaploidCallerParams, Site, SiteGenotype, VariantCall};
use navigator_analysis::coverage::{self, CallableLociParams, CoverageResult};
use navigator_analysis::gvcf;
use navigator_analysis::heteroplasmy::{self, HeteroplasmyParams};
use navigator_analysis::ibd::{ChromosomeGenotypes, GeneticMap, MatchSummary, PairwiseIbdDetector};
use navigator_analysis::scan::SampleSidecars;
use navigator_analysis::sidecar;
use navigator_domain::workspace::{Panel, PanelSite};
use navigator_store::panel;

// Re-export the analysis result types the command API returns, so the UI depends only
// on navigator-app (ui -> app), not directly on navigator-analysis.
pub use navigator_analysis::caller::SiteGenotype as PanelGenotype;
pub use navigator_analysis::caller::VariantCall as DenovoCall;
pub use navigator_analysis::coverage::CoverageResult as Coverage;
pub use navigator_analysis::haplo::{BranchEvidence, CallState, ScoredHaplogroup, SnpEvidence};
pub use navigator_analysis::heteroplasmy::HeteroplasmySite;
pub use navigator_analysis::mask::YRegionClass;
pub use navigator_analysis::mtvariants::{MtRegion, MtVariant, MtVariantKind};
pub use navigator_analysis::probe::AlignmentProbe;
pub use navigator_analysis::read_metrics::{PairOrientation, ReadMetrics};
pub use navigator_analysis::sex::{Confidence as SexConfidence, InferredSex, SexInferenceResult};
pub use navigator_analysis::sv::types::{SvAnalysisResult, SvCall, SvType};
pub use navigator_analysis::unified::UnifiedMetricsResult;
pub use navigator_domain::ancestry::{
    AncestryResult, AncestrySegment, ConfidenceInterval, PopulationComponent, SuperPopulationSummary,
};
// The ancestry panel format, re-exported so panel tooling/tests depend only on navigator-app.
pub use navigator_analysis::ancestry::{AncestryPanel, PanelSite as AncestryPanelSite};

/// A haplogroup assignment: the ranked candidates plus, for the reported terminal, the
/// child branches with per-SNP evidence (why descent stopped — unsupported splits show
/// ancestral SNPs, unresolved ones show no-calls).
#[derive(Debug, Clone)]
pub struct HaploAssignment {
    pub ranked: Vec<ScoredHaplogroup>,
    pub branches: Vec<BranchEvidence>,
    /// Per-SNP evidence along the placed lineage (root→terminal): every defining mutation the
    /// sample carries (or doesn't), Derived/Ancestral/NoCall. This is the set the multi-source
    /// variant/mutation **profile** reconciles — distinct from `branches`, which is the *untaken*
    /// child branches (explaining why descent stopped, hence largely ancestral/no-call).
    pub lineage: Vec<SnpEvidence>,
}

/// How a private (off-backbone) variant relates to the tree.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum PrivateClass {
    /// A known tree SNP off the assigned path — supports a finer/sibling branch.
    OffPathKnown(String),
    /// Not in the tree at all — a candidate for proposing a new branch.
    Novel,
}

/// A derived variant the sample carries that the haplogroup placement doesn't explain.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrivateVariant {
    pub position: i64,
    pub reference: char,
    pub alternate: char,
    pub depth: u32,
    pub allele_fraction: f64,
    pub class: PrivateClass,
    /// Curated CHM13 chrY structural class at this position (palindrome / amplicon / AZF-DYZ),
    /// if any — a paralog-prone zone where short-read mapping is unreliable, so the call is
    /// suspect (annotation only; not dropped). `None` = unique sequence, or a non-CHM13 build.
    #[serde(default)]
    pub region: Option<navigator_analysis::mask::YRegionClass>,
}

/// The private bucket for an alignment: de-novo Y calls not on the assigned backbone,
/// split into off-path-known (finer branches) and novel (new-branch candidates).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrivateBucket {
    pub terminal: String,
    pub variants: Vec<PrivateVariant>,
}

impl PrivateBucket {
    pub fn novel(&self) -> usize {
        self.variants.iter().filter(|v| v.class == PrivateClass::Novel).count()
    }
    pub fn off_path(&self) -> usize {
        self.variants
            .iter()
            .filter(|v| matches!(v.class, PrivateClass::OffPathKnown(_)))
            .count()
    }
    /// Calls that fall in a curated chrY structural (paralog-prone) region — suspect, to be
    /// down-weighted in reports rather than treated as confident new variants.
    pub fn in_structural_region(&self) -> usize {
        self.variants.iter().filter(|v| v.region.is_some()).count()
    }
    /// Novel calls in *unique* sequence (no structural-region flag) — the high-confidence
    /// new-branch candidates, separated from the paralog-zone noise.
    pub fn novel_in_unique_sequence(&self) -> usize {
        self.variants
            .iter()
            .filter(|v| v.class == PrivateClass::Novel && v.region.is_none())
            .count()
    }
}
pub use navigator_analysis::ibd::IbdSegment;
pub use navigator_analysis::ibd::{
    IbdDetectorConfig, IbdSegment as Segment, MatchSummary as IbdSummary, RelationshipEstimate,
};
// Sync/publish types the command API uses, re-exported so the UI depends only on navigator-app.
pub use ftdna_import::{
    FtdnaGenealogy, FtdnaImportOptions, FtdnaImportPlan, FtdnaImportSummary, FtdnaPlanRow, FtdnaPlanStats,
    FtdnaResolution, FtdnaSubjectInput, FuzzyCandidate, MatchKind,
};
pub use navigator_domain::identity::{ExternalId, FtdnaMember, Mdka};
pub use navigator_domain::ystr_cluster::{BranchSuggestion, ClusteredMember, YstrCluster, YstrClustering};
pub use navigator_refgenome::vcf_lift::infer_source_build as infer_vcf_source_build;
pub use navigator_refgenome::RefStatus;
use navigator_refgenome::{
    cache as refgenome_cache, canonical_build, Build as ReferenceBuild, LiftedPos, ReferenceGateway,
};
pub use navigator_refgenome::{ChromosomeRegions, Cytoband, GenomeRegions, RegionAnnotation};
pub use navigator_refgenome::{VcfLiftOpts, VcfLiftStats, VerifyOutcome};
use navigator_sync::exchange::{self, ExchangeKey};
use navigator_sync::{
    dev_http_client, login_default, AsyncSync, DeviceKey, OAuthConfig, RetryPolicy, TokenStore, DEVICE_KEY_COLLECTION,
};
pub use navigator_sync::{
    AlignmentRecord, BiosampleRecord, FeedPostRecord, PdsClient, PopulationBreakdownRecord, PrivateVariantsRecord,
    RecordRef, SequenceRunRecord, VariantCallEntry, NS_ALIGNMENT, NS_BIOSAMPLE, NS_FEED_POST, NS_POPULATION_BREAKDOWN,
    NS_SEQUENCERUN, PRIVATE_VARIANTS_COLLECTION,
};
use navigator_sync::{
    AuditEntryRecord, HaplogroupReconciliationRecord, HeteroplasmyObservationRecord, IdentityVerificationRecord,
    ManualOverrideRecord, ReconciliationStatusRecord, RunHaplogroupCallRecord, HAPLOGROUP_RECONCILIATION_COLLECTION,
};
use navigator_sync::{FedPopulationComponent, FedSuperPopulationSummary};
pub use recruitment::RecruitmentInvitation;
pub use social::{
    FederatedItem, FeedItem, FeedView, NotificationList, SocialMessage, SocialNotification, SocialThreadSummary,
};

/// Keychain service namespace for stored sessions (plan §7).
const KEYCHAIN_SERVICE: &str = "decodingus-navigator";

/// IBD comparison result between two samples.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdComparison {
    pub summary: MatchSummary,
    pub segments: Vec<IbdSegment>,
    /// Sites called in **both** samples — the effective comparison size. Sparse overlap (a
    /// chip↔chip pair, or chip↔WGS limited to the chip's sites) weakens short-segment calls, so
    /// it's surfaced rather than hidden.
    pub overlapping_sites: usize,
}

/// A sample for an IBD comparison — either a WGS/CRAM **alignment** (genotyped at the IBD-panel
/// sites) or an imported **chip** profile (resolved to the same CHM13 sites). Both yield dosages
/// over the canonical IBD panel, so the comparison is data-type-agnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IbdSource {
    Alignment(i64),
    Chip(i64),
}

/// The outcome of a federated IBD exchange over the encrypted channel (gap §4): the locally computed
/// match plus both signed [`IbdAttestation`]s. `agreed` ⇒ the partner's signature verified AND both
/// peers' summary hashes match (they computed the same result).
#[derive(Debug, Clone, PartialEq)]
pub struct IbdExchangeResult {
    pub summary: MatchSummary,
    pub segments: Vec<IbdSegment>,
    pub overlapping_sites: usize,
    pub my_attestation: IbdAttestation,
    pub partner_attestation: IbdAttestation,
    pub agreed: bool,
}

/// Presence + integrity of one ancestry/IBD reference asset, for the "data sources" transparency
/// affordance. `verified` is true only when a manifest lists the file and its SHA-256 matches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssetStatus {
    pub name: String,
    pub present: bool,
    pub verified: bool,
}

/// A pseudonymous federated-IBD candidate from the AppView's match engine. The
/// `suggested_sample_guid` is the AppView's opaque handle for the counterpart (not a DID,
/// not PII) — used to request an introduction. `signals` names the sources that contributed
/// (e.g. `POPULATION_OVERLAP`, `HAPLOGROUP`, `SHARED_MATCH`) behind the composite `score`.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdSuggestion {
    pub suggested_sample_guid: String,
    pub suggestion_type: String,
    pub score: f64,
    pub signals: Vec<String>,
}

/// Result of requesting an introduction to a candidate: the AppView's request URI and its
/// status (initially `PENDING`, awaiting the consent round-trip).
#[derive(Debug, Clone, PartialEq)]
pub struct IbdIntroResult {
    pub request_uri: String,
    pub status: String,
}

/// An inbound, **symmetric-blind** exchange request awaiting this account's consent (the initiator
/// is hidden until both parties consent). From `GET /api/v1/exchange/incoming`.
#[derive(Debug, Clone, PartialEq)]
pub struct IncomingRequest {
    pub request_uri: String,
    pub purpose: String,
    pub created_at: String,
}

/// A consent-ready exchange session (both parties consented): the partner's DID and their published
/// X25519 key URI are now revealed. From `GET /api/v1/exchange/pending`.
#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeSessionInfo {
    pub session_id: String,
    pub request_uri: String,
    pub purpose: String,
    pub partner_did: String,
    pub partner_key_uri: Option<String>,
}

/// Outcome of `POST /api/v1/exchange/consent`: `CONSENTED` (with the opened `session_id`),
/// `DECLINED`, or `PENDING` (recorded, awaiting the counterpart).
#[derive(Debug, Clone, PartialEq)]
pub struct ConsentOutcome {
    pub status: String,
    pub session_id: Option<String>,
}

/// A pulled relay envelope: the opaque ciphertext `blob` plus its routing (`from_did`/`seq`) and the
/// broker `id` to ack. From `GET /api/v1/exchange/relay/pull`.
#[derive(Debug, Clone, PartialEq)]
pub struct RelayEnvelope {
    pub id: i64,
    pub from_did: String,
    pub seq: i32,
    pub blob: String,
}

/// A live exchange session with a derived shared key, ready to seal/open payloads. Holds key
/// material, so it is deliberately not `Debug`/`Serialize` and should be kept in memory only.
#[derive(Clone)]
pub struct EstablishedSession {
    pub session_id: String,
    pub partner_did: String,
    key: [u8; 32],
}

/// Jetstream-ingest retry budget for a freshly-published device key: a 403 right after
/// publishing means the AppView hasn't ingested our `deviceKey` record yet. Exponential
/// backoff 1+2+4+8 s ≈ 15 s total before giving up.
const DEVICE_KEY_INGEST_RETRIES: u32 = 4;

/// Poll rounds (≈1s each) an IBD exchange waits for the partner's dosages / attestation.
const EXCHANGE_POLL_ROUNDS: u32 = 30;

/// The fed-record collections this client publishes — a PULL reconcile scans each (mirrors the
/// `publish_*` NSIDs). Derived-summary collections are tracked but not overwritten locally.
const PUBLISHED_COLLECTIONS: &[&str] = &[
    NS_BIOSAMPLE,
    NS_ALIGNMENT,
    NS_POPULATION_BREAKDOWN,
    NS_SEQUENCERUN,
    PRIVATE_VARIANTS_COLLECTION,
    HAPLOGROUP_RECONCILIATION_COLLECTION,
];

/// PDS collection NSID for a published IBD match attestation (the AppView indexes these via Jetstream).
const IBD_ATTESTATION_COLLECTION: &str = "com.decodingus.atmosphere.ibdAttestation";

/// Above this many sites, the exchanged dosage vector is decimated to fit the relay's 1 MiB envelope.
const EXCHANGE_SITE_BUDGET: usize = 100_000;
/// Decimation stride when over budget: keep sites at `position % N == 0`. A **position-based** rule
/// (not index) so both peers keep the *same physical sites* — preserving the IBD intersection — even
/// when their panels differ in size (WGS vs chip). Yields ~1/N of the canonical panel.
const EXCHANGE_DECIMATE: i64 = 16;

/// Downsample a dosage vector to fit the relay envelope, deterministically + cross-peer-aligned.
/// Small sets (synthetic tests, sparse chips) pass through untouched.
fn decimate_for_exchange(sites: Vec<IbdSite>) -> Vec<IbdSite> {
    if sites.len() <= EXCHANGE_SITE_BUDGET {
        return sites;
    }
    sites
        .into_iter()
        .filter(|s| s.position % EXCHANGE_DECIMATE == 0)
        .collect()
}

/// Parse the AppView's `/api/v1/ibd/suggestions` body into [`IbdSuggestion`]s. Lenient on
/// field casing (camel/snake) and on the `signals` shape (object map or array) so a minor
/// contract drift degrades gracefully rather than dropping every candidate.
fn parse_ibd_suggestions(body: &serde_json::Value) -> Vec<IbdSuggestion> {
    let Some(items) = body
        .get("items")
        .or_else(|| body.get("suggestions"))
        .and_then(|v| v.as_array())
    else {
        return Vec::new();
    };
    items
        .iter()
        .filter_map(|it| {
            let suggested_sample_guid = it
                .get("suggestedSampleGuid")
                .or_else(|| it.get("suggested_sample_guid"))
                .or_else(|| it.get("sampleGuid"))
                .and_then(|v| v.as_str())?
                .to_string();
            let suggestion_type = it
                .get("suggestionType")
                .or_else(|| it.get("suggestion_type"))
                .or_else(|| it.get("type"))
                .and_then(|v| v.as_str())
                .unwrap_or("unknown")
                .to_string();
            let score = it.get("score").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let signals = it
                .get("metadata")
                .and_then(|m| m.get("signals"))
                .or_else(|| it.get("signals"))
                .map(parse_ibd_signals)
                .unwrap_or_default();
            Some(IbdSuggestion {
                suggested_sample_guid,
                suggestion_type,
                score,
                signals,
            })
        })
        .collect()
}

/// Signal names. The AppView emits an array of plain strings
/// (`["POPULATION_OVERLAP", "HAPLOGROUP"]`); also tolerate an array of `{name|source}`
/// objects or an object map (keys) so a contract tweak degrades gracefully.
fn parse_ibd_signals(v: &serde_json::Value) -> Vec<String> {
    if let Some(arr) = v.as_array() {
        arr.iter()
            .filter_map(|s| {
                s.as_str()
                    .or_else(|| s.get("name").and_then(|x| x.as_str()))
                    .or_else(|| s.get("source").and_then(|x| x.as_str()))
                    .map(str::to_string)
            })
            .collect()
    } else if let Some(obj) = v.as_object() {
        obj.keys().cloned().collect()
    } else {
        Vec::new()
    }
}

/// Classify a non-2xx AppView response into a user-facing [`AppError::AppView`]. Consumes
/// `resp` to read the body (so capture the status first at the call site if also needed).
async fn appview_status_error(api: &str, resp: reqwest::Response) -> AppError {
    let status = resp.status();
    let body = resp.text().await.unwrap_or_default();
    match status.as_u16() {
        403 => AppError::AppView(format!(
            "{api}: device key not yet registered or verified by the AppView (403)"
        )),
        422 => AppError::AppView(format!(
            "{api}: request rejected, likely clock skew (422) — check the system clock"
        )),
        _ => AppError::AppView(format!("{api}: {status}: {body}")),
    }
}
pub use navigator_analysis::ibd_attest::{IbdAttestation, IbdExchangeMsg, IbdSite};
use navigator_domain::bisdna;
pub use navigator_domain::brief::{Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief};
use navigator_domain::chipprofile::{self, ChipProfile, NewChipProfile};
pub use navigator_domain::consensus::{DiploidSourceObs, DiploidVariant};
use navigator_domain::filetype;
pub use navigator_domain::filetype::DetectedData;
use navigator_domain::mtdna::{self, MtdnaSequence, NewMtdnaSequence};
use navigator_domain::reconciliation::{self, RunHaplogroupCall};
pub use navigator_domain::reconciliation::{
    AuditEntry, CompatibilityLevel, Consensus, DnaType, IdentityVerification, VerificationStatus,
};
use navigator_domain::strprofile::{self, NewStrProfile, StrProfile};
pub use navigator_domain::variants::SourceType;
use navigator_domain::variants::{self, NewVariantSet, VariantSet};
use navigator_domain::workspace::{
    Alignment, AnalysisArtifact, Biosample, NewAlignment, NewProject, NewSequenceRun, Project, SequenceRun,
};
pub use navigator_domain::ymatch::{Tmrca, YMatch, YSignal};
use navigator_domain::yprofile::{self, YObsInput};
pub use navigator_domain::yprofile::{YProfileSummary, YProfileVariant, YSourceObs, YState, YVariantStatus};
use navigator_domain::ysnp_dict::{self, YsnpDictionary};
pub use navigator_store::dm::{DmConversationSummary, DmMessage};
pub use navigator_store::ibd_exchange::StoredIbdExchange;
pub use navigator_store::source_file::SourceFile;
use navigator_store::{
    alignment, ancestry_result, artifact, biosample, chip_profile, consensus_painting, consensus_profile,
    haplogroup_call, mtdna as mtdna_store, project, reconciliation as recon_store, sequence_run, source_file,
    str_profile, sync_history, sync_outbox, sync_state, variant_set, Store, StoreError,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use uuid::Uuid;

pub mod error;
pub use error::AppError;
pub mod export;
pub mod settings;
pub mod sync_reconcile;
pub use settings::AppSettings;

/// Artifact kind for de-novo calls, keyed per contig so different contigs don't
/// overwrite each other in the cache.
fn denovo_kind(contig: &str) -> String {
    format!("denovo_snps:{contig}")
}

/// On-disk cache path for a downloaded haplotree, under `$NAVIGATOR_TREE_DIR` (tests/
/// overrides) or `~/.decodingus/trees`.
fn tree_cache_path(file: &str) -> PathBuf {
    let dir = std::env::var("NAVIGATOR_TREE_DIR")
        .ok()
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
            PathBuf::from(home).join(".decodingus").join("trees")
        });
    dir.join(file)
}

/// How long a cached haplotree is trusted before [`App::fetch_tree`] re-downloads it. The
/// AppView's curated tree changes slowly (curator review, periodic builds), so a weekly
/// refresh keeps placements current without hitting the network on every run. Override with
/// `NAVIGATOR_TREE_TTL_DAYS` (0 = always refetch).
const TREE_CACHE_TTL_DAYS_DEFAULT: u64 = 7;

/// Is the cached tree at `path` still within its TTL (default 7 days; `NAVIGATOR_TREE_TTL_DAYS`
/// overrides)? Unknown mtime / unreadable metadata → not fresh (forces a refresh attempt).
fn tree_cache_is_fresh(path: &Path) -> bool {
    let days = std::env::var("NAVIGATOR_TREE_TTL_DAYS")
        .ok()
        .and_then(|v| v.trim().parse::<u64>().ok())
        .or_else(|| AppSettings::load().tree_ttl_days)
        .unwrap_or(TREE_CACHE_TTL_DAYS_DEFAULT);
    let ttl = std::time::Duration::from_secs(days * 24 * 3600);
    std::fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|mtime| std::time::SystemTime::now().duration_since(mtime).ok())
        .map(|age| age < ttl)
        .unwrap_or(false)
}

/// Score a tree against the sample calls and attach the terminal's child-branch evidence.
///
/// The Kulczynski `score` ranks the candidates by proportional similarity (and supplies the
/// alternatives list), but the *reported terminal* is chosen in two steps: (1) the best-ranked
/// candidate the path-supported parsimony guard admits — i.e. whose lineage doesn't tunnel
/// through a branch the sample contradicts (the distal-Y paralog artifact); then (2)
/// [`haplo::deepen_terminal`] descends further into any child the sample clearly entered,
/// correcting under-calls at **unsplit tree nodes** (a half-ancestral SNP block scores below
/// its parent). The chosen node is moved to the front so every `ranked.first()` consumer
/// transparently gets it. See `documents/design/PangenomeExpansion.md`.
/// Pool every source's vote into one consensus map by a `SourceType`-weighted majority, keyed by
/// `K` (SNP **name** for Y — build-portable; rCRS **position** for mt) over value `V` (a per-SNP
/// **state** for Y — strand-/build-independent, since CHM13 vs GRCh38 can flip a base but not
/// "carries the derived allele"; a **base** for mt, which has one coordinate system). The weight
/// matches the variant reconcile's [`navigator_domain::consensus::obs_weight`] `SourceType` term;
/// the highest-weight value wins per key. The pooled set is placed on the tree **once** (genome-
/// level placement) instead of voting among per-run terminal labels.
/// Is a variant set's stored build GRCh38 (the FTDNA tree's native Y coordinate space)? `None`
/// (unknown build) is treated as GRCh38, matching the import default for a vendor Y VCF. Used to
/// gate which variant sets may pool into the GRCh38 genome consensus without liftover.
fn is_grch38_build(build: &Option<String>) -> bool {
    match build {
        None => true,
        Some(b) => {
            let b = b.to_ascii_lowercase();
            b.contains("grch38") || b.contains("hg38") || b == "38" || b == "b38"
        }
    }
}

fn pool_votes<K, V>(sources: &[(SourceType, HashMap<K, V>)]) -> HashMap<K, V>
where
    K: std::hash::Hash + Eq + Clone,
    V: std::hash::Hash + Eq + Clone,
{
    let mut tally: HashMap<K, HashMap<V, f64>> = HashMap::new();
    for (st, calls) in sources {
        let w = st.snp_weight();
        for (k, v) in calls {
            *tally.entry(k.clone()).or_default().entry(v.clone()).or_insert(0.0) += w;
        }
    }
    tally
        .into_iter()
        .filter_map(|(k, votes)| {
            votes
                .into_iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(v, _)| (k, v))
        })
        .collect()
}

fn assemble_assignment(tree: &navigator_analysis::haplo::HaploTree, calls: &HashMap<i64, char>) -> HaploAssignment {
    use navigator_analysis::haplo;
    let mut ranked = haplo::score(tree, calls);
    let terminal_id = ranked
        .iter()
        .find(|r| haplo::path_admissible(tree, calls, r.id))
        .map(|r| haplo::deepen_terminal(tree, calls, r.id));
    if let Some(tid) = terminal_id {
        if let Some(idx) = ranked.iter().position(|r| r.id == tid) {
            if idx != 0 {
                let chosen = ranked.remove(idx);
                ranked.insert(0, chosen);
            }
        }
    }
    let top = ranked.first().map(|t| t.id);
    let branches = top.map(|id| haplo::child_evidence(tree, calls, id)).unwrap_or_default();
    let lineage = top
        .map(|id| haplo::lineage_evidence(tree, calls, id))
        .unwrap_or_default();
    HaploAssignment {
        ranked,
        branches,
        lineage,
    }
}

/// Terminal selection for **named Y-SNP panel** data (BISDNA chip), as opposed to the
/// alignment-tuned [`assemble_assignment`]. Such panels give confident but sparse genotype
/// calls: a handful of recurrent or mis-probed ancestral calls on backbone nodes can make the
/// strict `path_admissible` guard (designed to kill distal tunnel artifacts in *coverage-
/// limited* alignment data) veto the genuine deep lineage, dropping the call to a shallow node
/// (e.g. A1). With confident chip calls that failure mode dominates, so here we trust the
/// proportional Kulczynski top — robust to a few stray calls — then [`deepen_terminal`] into
/// clearly-entered children. (Validated: this kit's chromo2 export → R-S1121 on both the
/// DecodingUs/hs1 and FTDNA/GRCh38 trees, on the lineage to its WGS-confirmed R-FGC29071.)
fn assemble_assignment_robust(
    tree: &navigator_analysis::haplo::HaploTree,
    calls: &HashMap<i64, char>,
) -> HaploAssignment {
    use navigator_analysis::haplo;
    let mut ranked = haplo::score(tree, calls);
    if let Some(top_id) = ranked.first().map(|r| r.id) {
        let terminal_id = haplo::deepen_terminal(tree, calls, top_id);
        // Parsimony back-off: don't report a deeper terminal than the evidence supports. Trim any
        // net-contradicted tail of the lineage (sparse-panel / damaged-aDNA over-deepening) while
        // a lone contradiction outweighed by deeper derived support still reaches the deep terminal.
        let chosen_id = support_backoff_terminal(tree, calls, terminal_id);
        if let Some(idx) = ranked.iter().position(|r| r.id == chosen_id) {
            if idx != 0 {
                let chosen = ranked.remove(idx);
                ranked.insert(0, chosen);
            }
        }
    }
    let top = ranked.first().map(|t| t.id);
    let branches = top.map(|id| haplo::child_evidence(tree, calls, id)).unwrap_or_default();
    let lineage = top
        .map(|id| haplo::lineage_evidence(tree, calls, id))
        .unwrap_or_default();
    HaploAssignment {
        ranked,
        branches,
        lineage,
    }
}

/// The root→`target` path of node ids (inclusive), or empty if `target` isn't reachable.
fn lineage_ids(tree: &navigator_analysis::haplo::HaploTree, target: i64) -> Vec<i64> {
    fn dfs(tree: &navigator_analysis::haplo::HaploTree, id: i64, target: i64, acc: &mut Vec<i64>) -> bool {
        let Some(node) = tree.nodes.get(&id) else { return false };
        acc.push(id);
        if id == target {
            return true;
        }
        for &c in &node.children {
            if dfs(tree, c, target, acc) {
                return true;
            }
        }
        acc.pop();
        false
    }
    let mut roots: Vec<i64> = tree.nodes.values().filter(|n| n.is_root).map(|n| n.id).collect();
    roots.sort_unstable();
    for r in roots {
        let mut acc = Vec::new();
        if dfs(tree, r, target, &mut acc) {
            return acc;
        }
    }
    Vec::new()
}

/// Root→`name` lineage of haplogroup names from the tree (empty if the name isn't found). Used to
/// derive a placed terminal's lineage path for cross-subject divergence/LCA without re-genotyping.
fn lineage_names(tree: &navigator_analysis::haplo::HaploTree, name: &str) -> Vec<String> {
    let Some(id) = tree.nodes.values().find(|n| n.name == name).map(|n| n.id) else {
        return Vec::new();
    };
    lineage_ids(tree, id)
        .into_iter()
        .filter_map(|i| tree.nodes.get(&i).map(|n| n.name.clone()))
        .collect()
}

/// Back off an over-deepened terminal to the node that maximizes running support along its
/// lineage. Walking root→terminal, each node contributes `(covered derived − covered ancestral)`
/// over its defining SNPs the sample has a call for; the chosen terminal is the deepest node at
/// which that running balance peaks. A net-contradicted tail (more ancestral than derived calls —
/// a sparse chip or degraded aDNA sample tunnelling into a wrong sub-clade) is trimmed, but a tail
/// whose deeper derived calls outweigh a shallow contradiction is kept (ties favour the deeper
/// node, preserving the robust "survive a lone backbone contradiction" behaviour). Returns
/// `terminal_id` unchanged when its lineage can't be traced.
fn support_backoff_terminal(
    tree: &navigator_analysis::haplo::HaploTree,
    calls: &HashMap<i64, char>,
    terminal_id: i64,
) -> i64 {
    let path = lineage_ids(tree, terminal_id);
    if path.is_empty() {
        return terminal_id;
    }
    let (mut balance, mut best_balance, mut best_id) = (0i32, i32::MIN, terminal_id);
    for &id in &path {
        let mut node_derived = false;
        if let Some(node) = tree.nodes.get(&id) {
            for l in &node.loci {
                let (Some(der), Some(anc)) = (l.derived.chars().next(), l.ancestral.chars().next()) else {
                    continue;
                };
                match calls.get(&l.position).map(|c| c.to_ascii_uppercase()) {
                    Some(b) if b == der.to_ascii_uppercase() => {
                        balance += 1;
                        node_derived = true;
                    }
                    Some(b) if b == anc.to_ascii_uppercase() => balance -= 1,
                    Some(_) => balance -= 1, // a third allele contradicts this branch
                    None => {}
                }
            }
        }
        // Deepen on strictly more support, or on a tie *only* when this node is itself
        // derived-supported. So a contradiction recovered by a deeper derived call still reaches
        // the deep terminal, while a net-negative tail or a flat run of marker-less nodes (the
        // sparse-panel / aDNA tunnel) is trimmed back to the last positively-supported node.
        if balance > best_balance || (balance == best_balance && node_derived) {
            best_balance = balance;
            best_id = id;
        }
    }
    best_id
}

/// Reconcile chip genotype calls to a haplotree's strand. Consumer arrays report alleles on the
/// reference plus strand, but a subset of sites sit on the opposite strand from the tree's
/// ancestral/derived convention. For each call at a tree position: keep the observed base if it
/// already equals the ancestral or derived allele; else substitute its complement when *that*
/// matches; else keep it (a genuine no-match the scorer will count against the branch). Positions
/// absent from the tree pass through unchanged (they don't affect scoring). This is a no-op for
/// dictionary-reconciled BISDNA calls (their base is always the derived allele), so it's safe to
/// apply on the shared chip-placement path.
fn strand_reconcile_to_tree(
    tree: &navigator_analysis::haplo::HaploTree,
    calls: HashMap<i64, char>,
) -> HashMap<i64, char> {
    let mut allowed: HashMap<i64, (char, char)> = HashMap::new();
    for node in tree.nodes.values() {
        for l in &node.loci {
            if let (Some(a), Some(d)) = (l.ancestral.chars().next(), l.derived.chars().next()) {
                allowed
                    .entry(l.position)
                    .or_insert((a.to_ascii_uppercase(), d.to_ascii_uppercase()));
            }
        }
    }
    calls
        .into_iter()
        .map(|(pos, base)| match allowed.get(&pos) {
            Some(&(a, d)) if base != a && base != d => {
                let c = complement_base(base);
                if c == a || c == d {
                    (pos, c)
                } else {
                    (pos, base)
                }
            }
            _ => (pos, base),
        })
        .collect()
}

/// Map GVCF-decoded bases at *lifted* positions back to tree positions (the GVCF-sourced
/// analogue of [`App::build_calls_from_lifted`]). A variant base wins; otherwise a callable
/// hom-ref lifted site takes the **reference base** at that lifted position — both reverse-
/// complemented for a minus-strand lift; otherwise the position is a no-call. `ref_base` is
/// keyed by lifted position (the GVCF/reference coordinate), not the tree position.
fn assemble_calls_lifted(
    called: &gvcf::CalledBases,
    lifted: &[LiftedPos],
    ref_base: &HashMap<i64, char>,
) -> HashMap<i64, char> {
    let mut calls = HashMap::new();
    for lp in lifted {
        let base = called.variant_bases.get(&lp.pos).copied().or_else(|| {
            called
                .callable
                .contains(&lp.pos)
                .then(|| ref_base.get(&lp.pos).copied())
                .flatten()
        });
        if let Some(b) = base {
            calls.insert(lp.tree_pos, if lp.reverse { complement_base(b) } else { b });
        }
    }
    calls
}

/// Minimum callable/calling depth adapted to read technology. The default (4) is a
/// short-read assumption — ~4 reads to call a base confidently. Long, accurate reads (HiFi,
/// mean read length > 1 kb) make a confident haploid observation from a *single* read, so a
/// ~4× HiFi sample is callable at 1×; clamping the floor at 2 needlessly threw away half its
/// already-shallow coverage. (ONT long reads are less accurate — revisit if we ever adapt by
/// platform rather than read length.)
fn adaptive_min_depth(base: u32, read_len: f64) -> u32 {
    if read_len > 1000.0 {
        1
    } else {
        base
    }
}

/// Haploid-caller params adapted to the sample's read tech (see [`adaptive_min_depth`]).
/// Sampled from the BAM head; falls back to defaults on any error. Blocking (reads the BAM)
/// — call inside `spawn_blocking`.
fn adaptive_haploid_params(bam_path: &Path, reference: Option<&Path>) -> HaploidCallerParams {
    let mut params = HaploidCallerParams::default();
    if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(bam_path, reference) {
        params.min_depth = adaptive_min_depth(params.min_depth, read_len);
    }
    params
}

/// Minimum genotyped sites for a reliable AIMs ancestry estimate (Scala `minSnpsAims`).
/// Overridable via `$NAVIGATOR_ANCESTRY_MIN_SNPS` (tests use a small panel).
fn ancestry_min_snps() -> usize {
    std::env::var("NAVIGATOR_ANCESTRY_MIN_SNPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3000)
}

/// Resolve an ancestry/IBD asset path under `<refgenome base>/ancestry/`: an `$<env_var>` override
/// (when non-empty) wins, else `<base>/ancestry/<stem>_<build>.<ext>`. The per-asset wrappers below
/// delegate here so the override+join+format pattern lives in one place.
fn ancestry_asset_path(env_var: &str, stem: &str, build: ReferenceBuild, ext: &str) -> PathBuf {
    if !env_var.is_empty() {
        if let Ok(p) = std::env::var(env_var) {
            return PathBuf::from(p);
        }
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("{stem}_{}.{ext}", build.as_str()))
}

/// Where the ancestry panel for `build` lives: `$NAVIGATOR_ANCESTRY_PANEL` (override), else
/// `<refgenome base>/ancestry/ancestry_panel_<build>.bin`. The offline `navigator-panelbuild`
/// tool writes it; install/ship copies it into the cache dir.
fn ancestry_panel_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_PANEL", "ancestry_panel", build, "bin")
}

/// Where the PCA loadings for `build` live: `$NAVIGATOR_ANCESTRY_PCA` (override), else
/// `<refgenome base>/ancestry/ancestry_pca_<build>.bin`. Optional — absent means the
/// AF-likelihood estimate runs without PCA coordinates.
fn ancestry_pca_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_PCA", "ancestry_pca", build, "bin")
}

/// Where the **ancient** PCA loadings for `build` live: `$NAVIGATOR_ANCESTRY_PCA_ANCIENT`
/// (override), else `<refgenome base>/ancestry/ancestry_pca_ancient_<build>.bin`. Optional —
/// present means the PCA-projection GMM runs against ancient reference components
/// (Steppe/EEF/WHG) instead of the modern super-populations. Must be built over the same panel
/// sites the AF panel genotypes (so the single genotyping pass covers it).
fn ancestry_pca_ancient_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_PCA_ANCIENT", "ancestry_pca_ancient", build, "bin")
}

/// The fine-population frequency asset path (`$NAVIGATOR_ANCESTRY_FREQ` override, else
/// `<base>/ancestry/ancestry_freq_global_<build>.bin`). Optional — fine admixture is skipped if absent.
fn ancestry_freq_global_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_FREQ", "ancestry_freq_global", build, "bin")
}

/// The chip-compatible IBD panel asset path (`$NAVIGATOR_IBD_PANEL` override, else
/// `<base>/ancestry/ibd_panel_<build>.bin`).
fn ibd_panel_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_IBD_PANEL", "ibd_panel", build, "bin")
}

/// The ancestry/IBD reference assets for the analysis build (CHM13), each with presence + manifest
/// verification — the "data sources" transparency line. Pure filesystem inspection (no analysis).
pub fn ancestry_asset_status() -> Vec<AssetStatus> {
    let build = ReferenceBuild::Chm13v2;
    let manifest = load_asset_manifest(build);
    [
        ("super-pop panel", ancestry_panel_path(build)),
        ("PCA (modern)", ancestry_pca_path(build)),
        ("PCA (ancient)", ancestry_pca_ancient_path(build)),
        ("fine frequencies", ancestry_freq_global_path(build)),
        ("genetic map", genetic_map_path(build)),
        ("IBD panel", ibd_panel_path(build)),
    ]
    .into_iter()
    .map(|(name, path)| {
        let bytes = std::fs::read(&path).ok();
        let verified = match (&manifest, &bytes, path.file_name().and_then(|n| n.to_str())) {
            (Some(m), Some(b), Some(fname)) => m.assets.contains_key(fname) && m.verify(fname, b).is_ok(),
            _ => false,
        };
        AssetStatus {
            name: name.to_string(),
            present: bytes.is_some(),
            verified,
        }
    })
    .collect()
}

/// What [`seed_bundled_assets`] copied into the cache on first run.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SeedSummary {
    pub copied: usize,
    pub skipped: usize,
}

/// Copy every regular file in `src_dir` into `dest_dir` that isn't already present there. Never
/// overwrites an existing file — a CDN-refreshed asset must win over the bundled one. Creates
/// `dest_dir`. A missing/unreadable `src_dir` is a no-op (returns the empty summary). Pure over the
/// two directories (no globals) so it's unit-testable.
pub fn seed_assets_from(src_dir: &Path, dest_dir: &Path) -> std::io::Result<SeedSummary> {
    let mut summary = SeedSummary::default();
    let Ok(entries) = std::fs::read_dir(src_dir) else {
        return Ok(summary); // no bundle present
    };
    std::fs::create_dir_all(dest_dir)?;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(name) = path.file_name() else { continue };
        // Skip hidden files (e.g. the staging `.staged` marker) — only real assets are seeded.
        if name.to_string_lossy().starts_with('.') {
            continue;
        }
        let dest = dest_dir.join(name);
        if dest.exists() {
            summary.skipped += 1;
        } else {
            std::fs::copy(&path, &dest)?;
            summary.copied += 1;
        }
    }
    Ok(summary)
}

/// Locate the bundled ancestry-asset resource directory shipped inside the installed image:
/// `$NAVIGATOR_BUNDLED_ASSETS` (override), else candidates relative to the running executable for
/// each packaged layout (macOS `.app` Resources, Linux `usr/lib|share/<app>`, Windows alongside).
/// `None` when running from a dev `target/` build with no bundle.
fn bundled_assets_dir() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_BUNDLED_ASSETS") {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
    }
    let exe = std::env::current_exe().ok()?;
    let dir = exe.parent()?;
    [
        dir.join("../Resources/ancestry"),         // macOS .app/Contents/MacOS → ../Resources
        dir.join("ancestry"),                      // Windows (alongside) / portable
        dir.join("../lib/DUNavigator/ancestry"),   // Linux .deb/AppImage usr/bin → usr/lib/<app>
        dir.join("../share/DUNavigator/ancestry"), // Linux usr/share/<app>
        dir.join("resources/ancestry"),            // generic
    ]
    .into_iter()
    .find(|c| c.is_dir())
}

/// Seed the bundled ancestry/IBD assets into `<cache base>/ancestry/` on first run (the offline
/// installer ships them as image resources; the runtime read path stays `~/.decodingus/...`). Copies
/// only the files missing from the cache, so a later manifest-verified CDN download transparently
/// overrides a bundled asset. Best-effort + non-fatal: no bundle (dev build) ⇒ empty summary.
pub fn seed_bundled_assets() -> SeedSummary {
    let Some(src) = bundled_assets_dir() else {
        return SeedSummary::default();
    };
    let dest = refgenome_cache::base_dir().join("ancestry");
    seed_assets_from(&src, &dest).unwrap_or_default()
}

/// The asset integrity manifest path for a build (`<base>/ancestry/ancestry_manifest_<build>.json`).
fn ancestry_manifest_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("", "ancestry_manifest", build, "json")
}

/// Load the build's asset manifest, if one is published. `None` (absent / unparseable) ⇒ integrity
/// checks are skipped (advisory).
fn load_asset_manifest(build: ReferenceBuild) -> Option<navigator_analysis::manifest::AssetManifest> {
    std::fs::read_to_string(ancestry_manifest_path(build))
        .ok()
        .and_then(|s| navigator_analysis::manifest::AssetManifest::from_json(&s).ok())
}

/// Read an asset file (`None` if absent), verifying its SHA-256 against the build manifest when one
/// is present. A **checksum mismatch is a hard error** — refuse a corrupt / truncated asset rather
/// than analyze against it. A missing manifest (or an unlisted file) passes through unverified.
fn read_verified_asset(build: ReferenceBuild, path: &Path) -> Result<Option<Vec<u8>>, AppError> {
    let Ok(bytes) = std::fs::read(path) else {
        return Ok(None);
    };
    if let Some(manifest) = load_asset_manifest(build) {
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if let Err((expected, got)) = manifest.verify(name, &bytes) {
                return Err(AppError::Import(format!(
                    "asset {name} failed its integrity check (manifest sha256 {expected}, file {got}) — re-download it"
                )));
            }
        }
    }
    Ok(Some(bytes))
}

/// The genetic-map asset path for a build (`$NAVIGATOR_GENETIC_MAP` override, else
/// `<base>/ancestry/genetic_map_<build>.bin`). Optional — IBD falls back to a uniform map if absent.
fn genetic_map_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_GENETIC_MAP", "genetic_map", build, "bin")
}

/// Load the real recombination map for IBD if the asset is present, else fall back to a uniform
/// 1 cM/Mb map over `lengths` (logged). `lengths` is the observed `(chromosome, max_bp)` per the
/// compared samples — used only for the uniform fallback.
fn load_genetic_map(build: ReferenceBuild, lengths: &[(&str, i32)]) -> GeneticMap {
    let path = genetic_map_path(build);
    let bytes = read_verified_asset(build, &path).unwrap_or_else(|e| {
        eprintln!("{e}"); // integrity mismatch on an optional asset → fall through to uniform
        None
    });
    match bytes.and_then(|b| GeneticMap::from_bytes(&b).ok()) {
        Some(m) => m,
        None => {
            eprintln!(
                "genetic map {} not found — IBD using uniform 1 cM/Mb (segment cM + relationship bands are approximate)",
                path.display()
            );
            GeneticMap::uniform(1.0, lengths)
        }
    }
}

/// Map a computed [`AncestryResult`] onto the shared federated wire record. The analysis
/// method is carried verbatim from the estimator that produced the result (never inferred),
/// so the published `analysisMethod` always matches the composition shown.
/// How many outbox rows a single [`App::drain_outbox`] pass attempts.
const OUTBOX_BATCH: i64 = 16;

/// Exponential backoff for a transient publish failure: `2^attempt` minutes, capped at 1 hour
/// (mirrors the legacy Scala sync queue). `attempt` is the 1-based retry count.
fn backoff_secs(attempt: i64) -> i64 {
    let minutes = 1i64.checked_shl(attempt.clamp(0, 16) as u32).unwrap_or(i64::MAX);
    (minutes.saturating_mul(60)).min(3600)
}

/// A request to export a cached result as a file body (gap §6). The id is the alignment id, except
/// [`Self::MtdnaTsv`] whose id is the mtDNA-sequence id. Carries enough for the UI to suggest a
/// filename + dialog filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportRequest {
    CoverageTsv(i64),
    CoverageHtml(i64),
    ReadMetricsTsv(i64),
    AncestryTsv(i64),
    AncestryHtml(i64),
    CallableBed(i64),
    MtdnaTsv(i64),
    /// Whole-genome diploid variant calls (SNV + indel) for an alignment, as a VCF. Heavy — re-walks
    /// the BAM per primary chromosome (cached).
    DiploidVcf(i64),
    /// The subject-level **consensus** diploid VCF — the joint genotype across the subject's
    /// same-build alignments. Heavy (call + force-call per alignment).
    ConsensusDiploidVcf(SampleGuid),
}

impl ExportRequest {
    /// File extension (no dot) for the save dialog + filter.
    pub fn extension(&self) -> &'static str {
        match self {
            ExportRequest::CoverageHtml(_) | ExportRequest::AncestryHtml(_) => "html",
            ExportRequest::CallableBed(_) => "bed",
            ExportRequest::DiploidVcf(_) | ExportRequest::ConsensusDiploidVcf(_) => "vcf",
            _ => "tsv",
        }
    }

    /// A short human label for the kind of export (status messages).
    pub fn label(&self) -> &'static str {
        match self {
            ExportRequest::CoverageTsv(_) => "coverage (TSV)",
            ExportRequest::CoverageHtml(_) => "coverage (HTML)",
            ExportRequest::ReadMetricsTsv(_) => "read metrics (TSV)",
            ExportRequest::AncestryTsv(_) => "ancestry (TSV)",
            ExportRequest::AncestryHtml(_) => "ancestry (HTML)",
            ExportRequest::CallableBed(_) => "callable loci (BED)",
            ExportRequest::MtdnaTsv(_) => "mtDNA variants (TSV)",
            ExportRequest::DiploidVcf(_) => "diploid variants (VCF)",
            ExportRequest::ConsensusDiploidVcf(_) => "consensus diploid (VCF)",
        }
    }

    /// A suggested default filename (`<stem>_<id>.<ext>`) for the save dialog.
    pub fn default_filename(&self) -> String {
        let (stem, id) = match self {
            ExportRequest::CoverageTsv(id) | ExportRequest::CoverageHtml(id) => ("coverage", id),
            ExportRequest::ReadMetricsTsv(id) => ("read_metrics", id),
            ExportRequest::AncestryTsv(id) | ExportRequest::AncestryHtml(id) => ("ancestry", id),
            ExportRequest::CallableBed(id) => ("callable", id),
            ExportRequest::MtdnaTsv(id) => ("mtdna_variants", id),
            ExportRequest::DiploidVcf(id) => ("diploid_variants", id),
            ExportRequest::ConsensusDiploidVcf(_) => return format!("consensus_diploid.{}", self.extension()),
        };
        format!("{stem}_{id}.{}", self.extension())
    }
}

/// One row of the Y-STR concordance view: a marker called from sequence (FTDNA-convention value +
/// calibration status) alongside the subject's imported vendor value, and whether they agree.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StrConcordanceRow {
    pub marker: String,
    /// Called FTDNA-convention value, or `None` if the marker wasn't called from sequence.
    pub called: Option<i32>,
    /// Calibration status: `Reliable` | `ConventionOffset` | `Excluded` | `Uncalibrated` | `NotCalled`.
    pub status: String,
    /// Whether the called value is corpus-calibrated (Reliable/ConventionOffset) — i.e. comparable.
    pub calibrated: bool,
    /// Imported vendor value (e.g. `"13"`, or a multi-copy `"11-15"`), or `None` if not in the profile.
    pub imported: Option<String>,
    pub depth: u32,
    /// Calibrated call whose value matches the imported single value.
    pub agree: bool,
}

/// Outcome of a batch [`App::add_data_batch`] run over multiple files / folders — for the import
/// summary the GUI shows after a multi-file Add Data or a drag-and-drop.
#[derive(Debug, Clone, Default)]
pub struct BatchImportSummary {
    /// `(filename, detected-type description)` for each successfully imported file.
    pub imported: Vec<(String, String)>,
    /// `(filename, reason)` for each file skipped or errored (unrecognized / import failure).
    pub skipped: Vec<(String, String)>,
}

/// A file's cheap signature (`mtime_secs:size`) for analysis-cache staleness — no content read.
/// `None` if the file is missing / unstattable.
fn file_signature(path: &Path) -> Option<String> {
    let meta = std::fs::metadata(path).ok()?;
    let mtime = meta
        .modified()
        .ok()?
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?
        .as_secs();
    Some(format!("{mtime}:{}", meta.len()))
}

/// Whether a cached artifact is still fresh for the current source signature: a stale entry (both
/// signatures known and differing) is rejected; an unknown stored sig (legacy / non-file source) or
/// an unknown current sig (file gone → nothing to recompute against) is trusted.
fn artifact_is_fresh(stored: Option<&str>, current: Option<&str>) -> bool {
    match (stored, current) {
        (Some(s), Some(c)) => s == c,
        _ => true,
    }
}

/// A recognized data-file extension — the pre-filter for directory expansion (a dropped folder is
/// walked for these). `add_data` re-sniffs text files (csv/tsv/txt) to route chip / STR / variants.
fn is_recognized_data_file(path: &Path) -> bool {
    let n = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    [
        ".bam",
        ".cram",
        ".vcf",
        ".vcf.gz",
        ".fasta",
        ".fa",
        ".fna",
        ".fas",
        ".fasta.gz",
        ".fa.gz",
        ".fna.gz",
        ".csv",
        ".tsv",
        ".txt",
    ]
    .iter()
    .any(|e| n.ends_with(e))
}

/// Collect recognized data files from `path`: a file yields itself (if recognized); a directory is
/// walked recursively. Bounded depth + file count so dropping a large tree can't run away.
fn collect_data_files(path: &Path, out: &mut Vec<PathBuf>, depth: usize) {
    const MAX_DEPTH: usize = 4;
    const MAX_FILES: usize = 2000;
    if out.len() >= MAX_FILES {
        return;
    }
    if path.is_dir() {
        if depth > MAX_DEPTH {
            return;
        }
        if let Ok(rd) = std::fs::read_dir(path) {
            let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).collect();
            entries.sort();
            for e in entries {
                collect_data_files(&e, out, depth + 1);
                if out.len() >= MAX_FILES {
                    break;
                }
            }
        }
    } else if is_recognized_data_file(path) {
        out.push(path.to_path_buf());
    }
}

/// Immediate child directories of `root` that contributed at least one collected file — the
/// "this folder holds several samples" signal. A single sample's folder fans data into at most one
/// subdirectory (e.g. FTDNA `<sample>/<kit>/<uuid>.bam` plus a top-level results CSV → just the kit
/// dir); a *parent* of many per-sample folders spreads it across several. Files sitting directly in
/// `root` aren't counted (they belong to the picked folder itself).
fn contributing_subdirs(root: &std::path::Path, files: &[PathBuf]) -> std::collections::BTreeSet<String> {
    use std::path::Component;
    let mut set = std::collections::BTreeSet::new();
    for f in files {
        if let Ok(rel) = f.strip_prefix(root) {
            let mut comps = rel.components();
            if let Some(Component::Normal(first)) = comps.next() {
                // Count it only when the file lives *inside* this child dir (a further component
                // follows), not when it sits directly at the root level.
                if comps.next().is_some() {
                    set.insert(first.to_string_lossy().into_owned());
                }
            }
        }
    }
    set
}

/// Whether `name` is a primary chromosome (1–22, X, Y, M/MT), with or without a `chr` prefix —
/// the contigs the whole-genome diploid caller runs over (skipping alts / decoys / unplaced).
fn is_primary_contig(name: &str) -> bool {
    let s = name.strip_prefix("chr").unwrap_or(name).to_ascii_uppercase();
    matches!(s.as_str(), "X" | "Y" | "M" | "MT") || s.parse::<u32>().map(|n| (1..=22).contains(&n)).unwrap_or(false)
}

/// The result of one [`App::drain_outbox`] pass — what the UI reports / shows in its indicator.
#[derive(Debug, Clone, Default)]
pub struct DrainOutcome {
    /// `(kind, at-uri)` of each row published this pass.
    pub published: Vec<(String, String)>,
    /// Rows that hit a non-transient error and were marked FAILED.
    pub failed: usize,
    /// Whether a transient failure rescheduled a row (i.e. we're likely offline).
    pub retry_scheduled: usize,
    /// Rows still awaiting a successful push after this pass.
    pub pending: i64,
}

/// The result of a [`App::pull_sync`] reconcile pass over the account's PDS records.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PullOutcome {
    /// Records unchanged since our last push.
    pub in_sync: usize,
    /// Records changed on the PDS and applied locally (where applicable).
    pub applied: usize,
    /// Remote records with no local mapping (PII-free summaries — tracked, not reconstructed).
    pub adopted: usize,
    /// Locally-published records missing on the PDS — flagged for re-publish.
    pub repushed: usize,
    /// Records that diverged on both sides (remote won, logged).
    pub conflicts: usize,
}

fn population_breakdown_record(result: &AncestryResult) -> PopulationBreakdownRecord {
    let components = result
        .components
        .iter()
        .map(|c| FedPopulationComponent {
            population: c.population_code.clone(),
            population_name: Some(c.population_name.clone()),
            percentage: c.percentage.into(),
            rank: Some(c.rank as i64),
        })
        .collect();
    let super_population_summary = result
        .super_population_summary
        .iter()
        .map(|s| FedSuperPopulationSummary {
            super_population: s.super_population.clone(),
            percentage: s.percentage.into(),
            populations: s.populations.clone(),
        })
        .collect();
    PopulationBreakdownRecord::new(
        result.method.clone(),
        result.panel_type.clone(),
        Some(result.reference_version.clone()),
        result.snps_analyzed as i64,
        result.snps_with_genotype as i64,
        result.snps_missing as i64,
        result.confidence_level,
        components,
        super_population_summary,
        result.pca_coordinates.clone(),
        Utc::now().to_rfc3339(),
    )
    .with_fit_distance(result.fit_distance)
}

/// Build a community feed-post record — the shared `com.decodingus.atmosphere.feed.post`
/// contract the AppView mirrors into its feed (top-level `createdAt`, optional topic /
/// reply pointers). PII-free beyond the text the user chose to publish. `reply` is the
/// `(root_uri, parent_uri)` pair on a threaded reply (`None` for a top-level post).
fn feed_post_record(content: &str, topic: Option<&str>, reply: Option<(&str, &str)>) -> serde_json::Value {
    let mut rec = FeedPostRecord::new(content, Utc::now().to_rfc3339()).with_topic(topic.map(str::to_string));
    if let Some((root, parent)) = reply {
        rec = rec.with_reply(root, parent);
    }
    // A struct of plain strings always serializes; surface a build bug loudly rather than
    // silently dropping the post.
    serde_json::to_value(&rec).expect("feed-post record serializes")
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
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if name.contains("chm13") {
        "chm13v2.0".to_string()
    } else {
        "unknown".to_string()
    }
}

/// Cheap VCF header peek: the `##` meta block (joined) + the contig names from `##contig=<ID=…>`.
/// Reads only the header (stops at the first data line). Plain text — matches the import parser,
/// which doesn't decompress either; a gzipped VCF simply yields an empty peek (→ generic).
fn peek_vcf_header(path: &Path) -> (String, Vec<String>) {
    use std::io::BufRead;
    let Ok(file) = std::fs::File::open(path) else {
        return (String::new(), Vec::new());
    };
    let mut meta = String::new();
    let mut contigs = Vec::new();
    for line in std::io::BufReader::new(file).lines().map_while(Result::ok) {
        if let Some(rest) = line.strip_prefix("##") {
            meta.push_str(&line);
            meta.push('\n');
            // ##contig=<ID=chrY,length=…>
            if let Some(after) = rest.strip_prefix("contig=<ID=") {
                let id: String = after.chars().take_while(|&c| c != ',' && c != '>').collect();
                if !id.is_empty() {
                    contigs.push(id);
                }
            }
        } else if line.starts_with('#') {
            continue; // the #CHROM column line — header still, no useful meta
        } else {
            break; // first data record → header done
        }
    }
    (meta, contigs)
}

/// Parse a VCF into the **subject's** SNP calls, honoring the genotype.
///
/// A vendor VCF (FTDNA Big Y, YSEQ, …) lists a site's REF/ALT even where the sample is
/// homozygous-reference (`GT 0/0`) — e.g. `chrY 2781955 C T … 0/0`, where the sample is C, not T.
/// Taking `ALT[0]` blindly (as a sites-only VCF parser does) records that T as a derived call, and
/// a Big Y export carries thousands of such reference sites → the placement deepens into branches
/// the sample doesn't actually carry. So when a genotyped sample column is present we read its `GT`
/// and keep a single-base ALT only when the genotype selects it (the first non-zero allele,
/// multi-allelic-aware); `0/0` and `./.` rows are dropped. A VCF with no FORMAT/sample column
/// (a sites-only list) keeps its old meaning: every listed ALT is one of the subject's variants.
fn parse_vcf_subject_snps(path: &Path) -> Result<Vec<variants::VariantCall>, AppError> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let Ok(pos) = f[1].parse::<i64>() else { continue };
        let (chrom, id, reference, alt_field) = (f[0], f[2], f[3], f[4]);
        if alt_field == "." {
            continue; // no ALT listed → nothing to call
        }
        let alts: Vec<&str> = alt_field.split(',').collect();

        // Genotyped (FORMAT + ≥1 sample, with a GT key) → honor the call; else sites-only.
        let gt = (f.len() >= 10)
            .then(|| {
                f[8].split(':')
                    .position(|k| k == "GT")
                    .and_then(|i| f[9].split(':').nth(i))
            })
            .flatten();
        let (alt, genotype) = match gt {
            Some(gt) => {
                // First non-zero allele index selects the carried ALT; all-zero (0/0) or no-call
                // (./.) means the subject is reference here — skip it.
                match gt
                    .split(['/', '|'])
                    .filter_map(|a| a.parse::<usize>().ok())
                    .find(|&a| a > 0)
                {
                    Some(idx) => match alts.get(idx - 1) {
                        Some(&a) => (a, Some(gt.to_string())),
                        None => continue,
                    },
                    None => continue,
                }
            }
            None => (alts[0], None), // sites-only VCF: the listed variant is the subject's
        };
        if let Some(call) = variants::snp_call(chrom, pos, reference, alt, Some(id.to_string()), genotype) {
            out.push(call);
        }
    }
    Ok(out)
}

/// Detect the reference build from VCF meta lines (`##reference=…`, `##contig assembly=…`).
fn detect_vcf_build(meta: &str) -> Option<String> {
    let l = meta.to_lowercase();
    if l.contains("chm13") || l.contains("t2t") || l.contains("hs1") {
        Some("chm13v2.0".into())
    } else if l.contains("hg38") || l.contains("grch38") {
        Some("GRCh38".into())
    } else if l.contains("hg19") || l.contains("grch37") {
        Some("GRCh37".into())
    } else {
        None
    }
}

/// Read a sibling `readme.txt` (FTDNA Big Y bundles one beside `variants.vcf`), if present.
fn sibling_readme(path: &Path) -> Option<String> {
    let dir = path.parent()?;
    for name in ["readme.txt", "README.txt", "README"] {
        if let Ok(text) = std::fs::read_to_string(dir.join(name)) {
            return Some(text);
        }
    }
    None
}

/// Best-effort vendor label for an mtDNA FASTA, from the file name + defline (FTDNA mtFull / YSEQ).
fn mt_vendor_label(filename: Option<&str>, defline: Option<&str>) -> &'static str {
    let hay = format!("{} {}", filename.unwrap_or(""), defline.unwrap_or("")).to_lowercase();
    if hay.contains("ftdna") || hay.contains("familytreedna") || hay.contains("mtfull") {
        "FTDNA mtFull Sequence"
    } else if hay.contains("yseq") {
        "YSEQ mtDNA"
    } else {
        "mtDNA FASTA"
    }
}

/// A disambiguating label context for a vendor VCF: the parent directory when the file name is the
/// generic vendor name (`variants.vcf`), else the file name itself.
fn vcf_label_context(path: &Path, filename: &str) -> String {
    let generic = matches!(
        filename.to_ascii_lowercase().as_str(),
        "variants.vcf" | "variants.vcf.gz"
    );
    if generic {
        if let Some(parent) = path.parent().and_then(|p| p.file_name()).and_then(|n| n.to_str()) {
            if !parent.is_empty() {
                return parent.to_string();
            }
        }
    }
    filename.to_string()
}

/// SHA-256 of a file's content (hex), computed off the async runtime (streamed — safe for large
/// alignments).
async fn sha256_file_async(path: PathBuf) -> Result<String, AppError> {
    let hash = tokio::task::spawn_blocking(move || du_bio::hash::sha256_file(&path)).await??;
    Ok(hash)
}

/// SHA-256 of an in-memory string (hex) — for hashing tree JSON / small content.
fn sha256_str(s: &str) -> String {
    du_bio::hash::sha256_hex(s.as_bytes())
}

/// Reconstruct a minimal [`HaploAssignment`] from a recorded call — the terminal + lineage,
/// without the full ranked list or branch evidence. Returned on a scoring cache hit (the
/// recorded call is the source of truth; the detail is only needed on a fresh score).
fn assignment_from_call(call: &navigator_domain::reconciliation::RunHaplogroupCall) -> HaploAssignment {
    HaploAssignment {
        ranked: vec![navigator_analysis::haplo::ScoredHaplogroup {
            id: 0,
            name: call.haplogroup.clone(),
            score: call.score,
            depth: call.lineage.len(),
            lineage: call.lineage.clone(),
            matched: call.matched.max(0) as usize,
            expected: call.expected.max(0) as usize,
            found: 0,
        }],
        branches: Vec::new(),
        lineage: Vec::new(),
    }
}

// Watson–Crick complement of a base (for reverse-strand lifts) — the shared helper in
// navigator-domain (also used by the chip dosage + BISDNA QC paths).
use navigator_domain::seq::complement_base;

/// The build a haplotree's positions are in, by contig: the FTDNA Y tree is GRCh38; mtDNA
/// (`chrM`) is rCRS and stays a direct query (no chain), so it returns `None`.
fn tree_build_for_contig(contig: &str) -> Option<&'static str> {
    if contig.eq_ignore_ascii_case("chrY") {
        Some("GRCh38")
    } else {
        None
    }
}

/// Which Y-DNA haplogroup tree to place against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum YTreeProvider {
    /// The DecodingUs tree served by our AppView (native multi-build coords incl. CHM13/`hs1`).
    DecodingUs,
    /// FTDNA's public Y-DNA haplotree (GRCh38; lifted onto the alignment build).
    Ftdna,
}

/// Selected Y-tree provider. Defaults to **DecodingUs** (our tree; native CHM13 coordinates →
/// no liftover). Override with `NAVIGATOR_Y_TREE_PROVIDER=ftdna|decodingus`.
/// Resolve the Y-tree provider given the env override and the settings value (pure; env wins →
/// settings → default DecodingUs).
fn resolve_y_provider(env: Option<&str>, settings: Option<&str>) -> YTreeProvider {
    match env.or(settings).map(str::trim) {
        Some(v) if v.eq_ignore_ascii_case("ftdna") => YTreeProvider::Ftdna,
        _ => YTreeProvider::DecodingUs,
    }
}

fn y_tree_provider() -> YTreeProvider {
    let env = std::env::var("NAVIGATOR_Y_TREE_PROVIDER").ok();
    let settings = AppSettings::load().y_tree_provider;
    resolve_y_provider(env.as_deref(), settings.as_deref())
}

/// Application interface mode: a casual, single-person experience with plain-language briefs
/// (`Simple`) vs. the full power-user UI with projects and per-source analysis (`Advanced`). This is
/// app-level UI state, persisted in [`AppSettings`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UiMode {
    Simple,
    Advanced,
}

impl UiMode {
    /// The settings-file token (`"simple"` / `"advanced"`).
    pub fn as_str(self) -> &'static str {
        match self {
            UiMode::Simple => "simple",
            UiMode::Advanced => "advanced",
        }
    }

    /// Parse a settings/env token; unrecognized values yield `None`.
    pub fn parse(s: &str) -> Option<UiMode> {
        match s.trim().to_ascii_lowercase().as_str() {
            "simple" => Some(UiMode::Simple),
            "advanced" => Some(UiMode::Advanced),
            _ => None,
        }
    }
}

/// Resolve an explicitly-configured UI mode (pure; env wins → settings → `None`). `None` means the
/// user has never pinned a mode, so the UI applies its first-run heuristic.
fn resolve_ui_mode(env: Option<&str>, settings: Option<&str>) -> Option<UiMode> {
    // A recognized env value wins; an unrecognized one is ignored and falls through to settings.
    env.and_then(UiMode::parse).or_else(|| settings.and_then(UiMode::parse))
}

/// The configured UI mode, honoring `NAVIGATOR_UI_MODE` over the persisted setting. `None` when
/// neither is set (first run → the UI defaults from a workspace heuristic).
pub fn configured_ui_mode() -> Option<UiMode> {
    let env = std::env::var("NAVIGATOR_UI_MODE").ok();
    let settings = AppSettings::load().ui_mode;
    resolve_ui_mode(env.as_deref(), settings.as_deref())
}

/// Persist the chosen UI mode, preserving all other settings.
pub fn persist_ui_mode(mode: UiMode) -> std::io::Result<()> {
    let mut s = AppSettings::load();
    s.ui_mode = Some(mode.as_str().to_string());
    s.save()
}

/// Base URL of the DecodingUs AppView serving the tree API. Local by default for testing;
/// switch with `DECODINGUS_APPVIEW_URL` (e.g. the production host at cutover).
/// Resolve the AppView base URL (pure; env wins → settings → default localhost; trailing slash
/// trimmed; blank values ignored).
fn resolve_appview_url(env: Option<String>, settings: Option<String>) -> String {
    env.or(settings)
        .map(|s| s.trim_end_matches('/').to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "http://localhost:9000".to_string())
}

fn decodingus_appview_url() -> String {
    resolve_appview_url(
        std::env::var("DECODINGUS_APPVIEW_URL").ok(),
        AppSettings::load().appview_url,
    )
}

/// A subject's multi-source variant **consensus profile** for one DNA type (Y today; mtDNA /
/// autosomal adapters reuse this aggregate + the generic engine). Persisted as a snapshot (serialized
/// to the `consensus_profile` table's payload, keyed by `(biosample, dna_type)`) so
/// [`App::cached_consensus_profile`] can reload it without re-genotyping.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ConsensusProfile {
    pub variants: Vec<YProfileVariant>,
    pub summary: YProfileSummary,
    /// Consensus lineage label (terminal Y/mt haplogroup) across sources, if any. `None` for DNA
    /// types without a lineage label (e.g. autosomal).
    pub terminal: Option<String>,
    /// Per-source provenance (which tests contributed, and how many variants each).
    #[serde(default)]
    pub sources: Vec<YSourceSummary>,
}

/// The Y-DNA view of a [`ConsensusProfile`] — the Y adapter is the first consumer of the generic
/// consensus aggregate; the name is kept for the Y-DNA tab + worker contract.
pub type YProfile = ConsensusProfile;

/// One contributing source in a [`ConsensusProfile`] (provenance display).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct YSourceSummary {
    pub label: String,
    pub source_type: SourceType,
    pub variant_count: usize,
}

/// A subject's **autosomal** multi-source consensus profile — the diploid (0/1/2) sibling of
/// [`ConsensusProfile`], over the canonical CHM13 IBD-panel sites. Persisted in the same
/// `consensus_profile` table under `dna_type='Auto'`. No lineage label (autosomes have no haplogroup).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DiploidProfile {
    pub variants: Vec<DiploidVariant>,
    pub summary: YProfileSummary,
    /// Per-source provenance (which tests contributed, and how many sites each).
    #[serde(default)]
    pub sources: Vec<YSourceSummary>,
}

/// `ancestry_result.alignment_id` sentinel marking a result derived from the subject's autosomal
/// **consensus** (pooled across all sources) rather than a single sequencing alignment.
pub const CONSENSUS_SOURCE_ID: i64 = 0;

/// Bridge the autosomal consensus to the genotype carrier the ancestry estimators + IBD detector
/// consume. Each reconciled site becomes a [`SiteGenotype`] with the consensus dosage (count of the
/// CHM13 ALT — the canonical orientation the AIM freq / PCA assets are keyed against). No-calls
/// (dosage -1) are carried through and ignored downstream like any missing genotype.
pub fn consensus_genotypes(profile: &DiploidProfile) -> Vec<SiteGenotype> {
    profile
        .variants
        .iter()
        .map(|v| SiteGenotype {
            name: v.name.clone(),
            contig: v.contig.clone(),
            position: v.position,
            reference_allele: v.reference.clone(),
            alternate_allele: v.alternate.clone(),
            ploidy: 2,
            dosage: v.consensus_dosage as i32,
            gq: 0,
            depth: 0,
            ref_depth: 0,
            alt_depth: 0,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        })
        .collect()
}

/// Flatten a placement's branch SNP evidence into per-SNP observations (deduped by name; a SNP
/// defines one branch, but guard against duplicates). `in_tree` is true for tree-defining SNPs.
/// Build per-SNP observations for the multi-source consensus profile from a placement's **lineage**
/// (root→terminal defining mutations the sample carries) — not its child branches, which are the
/// *untaken* deeper splits and are by construction ancestral/no-call. (Using `branches` made a
/// single-source profile read as all-no-call even when the terminal placed cleanly.)
fn snp_obs_from_assignment(assignment: &HaploAssignment, in_tree: bool) -> Vec<YObsInput> {
    let mut by_name: std::collections::HashMap<String, YObsInput> = std::collections::HashMap::new();
    for snp in &assignment.lineage {
        let state = match snp.state {
            CallState::Derived => YState::Derived,
            CallState::Ancestral => YState::Ancestral,
            CallState::NoCall => YState::NoCall,
        };
        by_name.entry(snp.name.clone()).or_insert_with(|| {
            YObsInput::snp(
                snp.name.clone(),
                snp.position,
                snp.ancestral.clone(),
                snp.derived.clone(),
                state,
                in_tree,
            )
        });
    }
    by_name.into_values().collect()
}

/// Per-build reference-genome status + override for the Settings UI.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefBuildStatus {
    /// Canonical build label (e.g. "GRCh38").
    pub build: String,
    /// Human-readable cache/override status.
    pub status: String,
    /// User-pinned local FASTA, if any.
    pub local_path: Option<String>,
    /// Whether a missing reference may be auto-downloaded.
    pub auto_download: bool,
}

mod analysis;
mod auth;
mod brief;
mod commands;
mod dm;
mod fastpath;
mod ftdna_import;
mod haplogroup;
mod ibd_exchange;
mod import_profiles;
mod import_unified;
mod publish;
mod queries;
mod recruitment;
mod social;
mod sync;

impl App {
    /// Reference-genome settings + cache status, one row per supported build.
    pub fn reference_settings(&self) -> Vec<RefBuildStatus> {
        let cfg = navigator_refgenome::UserConfig::load(&self.gateway.config_path());
        ReferenceBuild::all()
            .iter()
            .map(|&b| {
                let name = b.as_str();
                let ov = cfg.references.get(name);
                let status = match self.gateway.reference_status(name) {
                    RefStatus::LocalOverride(p) => format!("local file: {}", p.display()),
                    RefStatus::Cached(_) => "in cache".to_string(),
                    RefStatus::NeedsDownload { est_bytes, .. } => {
                        format!("not downloaded (~{} MB)", est_bytes / 1_000_000)
                    }
                    RefStatus::Unknown => "unknown".to_string(),
                };
                RefBuildStatus {
                    build: name.to_string(),
                    status,
                    local_path: ov.and_then(|o| o.local_path.clone()),
                    auto_download: ov.map(|o| o.auto_download).unwrap_or(true),
                }
            })
            .collect()
    }

    /// Set the local-FASTA override + auto-download flag for a build, persisting
    /// `reference_sources.json`. Applies on the next reference resolve (no restart needed — the
    /// gateway re-reads the file).
    pub fn set_reference_override(
        &self,
        build: &str,
        local_path: Option<String>,
        auto_download: bool,
    ) -> Result<(), AppError> {
        let path = self.gateway.config_path();
        let mut cfg = navigator_refgenome::UserConfig::load(&path);
        let key = canonical_build(build)
            .map(|b| b.as_str().to_string())
            .unwrap_or_else(|| build.to_string());
        let entry = cfg.references.entry(key).or_default();
        entry.local_path = local_path.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        entry.auto_download = auto_download;
        cfg.save(&path)?;
        Ok(())
    }
}

/// One instrument→lab association from the AppView `sequencer` endpoints (D8). Mirrors the
/// `SequencerLabDto` shape; extra fields are tolerated.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct SequencerLabInfo {
    pub instrument_id: String,
    pub lab_name: String,
    #[serde(default)]
    pub is_d2c: bool,
    #[serde(default)]
    pub manufacturer: Option<String>,
    #[serde(default)]
    pub model_name: Option<String>,
    #[serde(default)]
    pub website_url: Option<String>,
}

/// Map an alignment's reference build to the DecodingUs coordinate key (`"hs1"` for CHM13,
/// `"GRCh38"`, `"GRCh37"`). `None` for builds the tree has no coordinates for. Drives the
/// native-build (no-liftover) placement in `assign_y_decodingus`.
fn decodingus_build_key(reference_build: &str) -> Option<&'static str> {
    match canonical_build(reference_build) {
        Some(ReferenceBuild::Grch38) => Some("GRCh38"),
        Some(ReferenceBuild::Grch37) => Some("GRCh37"),
        Some(ReferenceBuild::Chm13v2) | Some(ReferenceBuild::Chm13v2MaskedRcrs) => Some("hs1"),
        None => None,
    }
}

/// Whether an alignment's reference build matches a GVCF name's build token (e.g. `chm13`),
/// compared on the canonical build so `chm13`/`chm13v2`/`hs1` all agree. A token that doesn't
/// resolve to a known build is treated as a non-match (fall back to the first alignment).
fn build_hint_matches(reference_build: &str, hint: &str) -> bool {
    match (canonical_build(reference_build), canonical_build(hint)) {
        (Some(a), Some(b)) => a == b,
        _ => false,
    }
}

/// Whether `<alignment>.crai`/`.bai` is present among the discovered index files.
fn has_sibling_index(aln_path: &Path, index_files: &[PathBuf]) -> bool {
    let Some(aln_name) = aln_path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    index_files
        .iter()
        .filter_map(|i| i.file_name().and_then(|n| n.to_str()))
        .any(|n| n == format!("{aln_name}.crai") || n == format!("{aln_name}.bai"))
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
/// Artifact kind for an alignment's cached IBD-panel genotypes (distinct from the store-panel
/// genotype cache, which is keyed by `panel_kind(panel_id, ploidy)`).
const IBD_PANEL_KIND: &str = "ibd_panel_genotypes";

/// The IBD-panel genotype cache kind, salted with the CHM13 panel asset's manifest sha256 (first 16
/// hex chars) so regenerating the panel auto-invalidates stale per-alignment genotypes. Falls back to
/// the bare kind when no manifest is published (genotypes then keyed only by `GENOTYPE_VERSION`).
fn ibd_panel_cache_kind() -> String {
    let build = ReferenceBuild::Chm13v2;
    let name = ibd_panel_path(build)
        .file_name()
        .and_then(|n| n.to_str())
        .map(String::from);
    match (load_asset_manifest(build), name) {
        (Some(m), Some(n)) => match m.assets.get(&n) {
            Some(e) => format!("{IBD_PANEL_KIND}:{}", &e.sha256[..16.min(e.sha256.len())]),
            None => IBD_PANEL_KIND.to_string(),
        },
        _ => IBD_PANEL_KIND.to_string(),
    }
}

/// Count of sites called (dosage within ploidy) in **both** samples — the effective IBD comparison
/// size, surfaced so a sparse chip↔chip / chip↔WGS overlap isn't mistaken for a confident result.
fn overlapping_called_sites(a: &[SiteGenotype], b: &[SiteGenotype]) -> usize {
    let called = |g: &SiteGenotype| (0..=g.ploidy as i32).contains(&g.dosage);
    let set: std::collections::HashSet<(&str, i64)> = a
        .iter()
        .filter(|g| called(g))
        .map(|g| (g.contig.as_str(), g.position))
        .collect();
    b.iter()
        .filter(|g| called(g))
        .filter(|g| set.contains(&(g.contig.as_str(), g.position)))
        .count()
}

/// Group two samples' dosages, load the genetic map for `build`, detect IBD segments, and record
/// the overlapping-site count. Shared by the alignment-pair and chip-or-WGS compare paths.
fn detect_ibd(
    ga: &[SiteGenotype],
    gb: &[SiteGenotype],
    build: ReferenceBuild,
    config: IbdDetectorConfig,
) -> IbdComparison {
    let overlapping_sites = overlapping_called_sites(ga, gb);
    let sample_a = group_chrom_genotypes(ga);
    let sample_b = group_chrom_genotypes(gb);
    let mut lengths: BTreeMap<String, i32> = BTreeMap::new();
    for sample in [&sample_a, &sample_b] {
        for (chr, cg) in sample {
            let m = cg.positions.last().copied().unwrap_or(1);
            lengths.entry(chr.clone()).and_modify(|e| *e = (*e).max(m)).or_insert(m);
        }
    }
    let pairs: Vec<(&str, i32)> = lengths.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    let gmap = load_genetic_map(build, &pairs);
    let segments = PairwiseIbdDetector::new(config).detect_segments(&sample_a, &sample_b, &gmap);
    let summary = MatchSummary::from_segments(&segments);
    IbdComparison {
        summary,
        segments,
        overlapping_sites,
    }
}

fn group_chrom_genotypes(genotypes: &[SiteGenotype]) -> std::collections::HashMap<String, ChromosomeGenotypes> {
    let mut by_contig: BTreeMap<String, Vec<(i64, i32)>> = BTreeMap::new();
    for g in genotypes {
        by_contig
            .entry(g.contig.clone())
            .or_default()
            .push((g.position, g.dosage));
    }
    by_contig
        .into_iter()
        .map(|(chrom, mut v)| {
            v.sort_by_key(|(p, _)| *p);
            let positions = v.iter().map(|(p, _)| *p as i32).collect();
            let dosages = v.iter().map(|(_, d)| *d as i8).collect();
            (
                chrom.clone(),
                ChromosomeGenotypes {
                    chromosome: chrom,
                    positions,
                    dosages,
                },
            )
        })
        .collect()
}

/// IBD detection over two [`IbdSite`] dosage vectors (the federated-exchange path — the partner's
/// dosages arrive as `IbdSite`, not [`SiteGenotype`]). Mirrors [`detect_ibd`] but groups directly
/// from the compact wire type.
fn detect_ibd_sites(
    my: &[IbdSite],
    partner: &[IbdSite],
    build: ReferenceBuild,
    config: IbdDetectorConfig,
) -> IbdComparison {
    let group = |sites: &[IbdSite]| -> std::collections::HashMap<String, ChromosomeGenotypes> {
        let mut by: BTreeMap<String, Vec<(i64, i32)>> = BTreeMap::new();
        for s in sites {
            by.entry(s.contig.clone()).or_default().push((s.position, s.dosage));
        }
        by.into_iter()
            .map(|(chrom, mut v)| {
                v.sort_by_key(|(p, _)| *p);
                let positions = v.iter().map(|(p, _)| *p as i32).collect();
                let dosages = v.iter().map(|(_, d)| *d as i8).collect();
                (
                    chrom.clone(),
                    ChromosomeGenotypes {
                        chromosome: chrom,
                        positions,
                        dosages,
                    },
                )
            })
            .collect()
    };
    // Overlapping called sites (dosage 0..=2 in both).
    let partner_called: HashMap<(&str, i64), ()> = partner
        .iter()
        .filter(|s| (0..=2).contains(&s.dosage))
        .map(|s| ((s.contig.as_str(), s.position), ()))
        .collect();
    let overlapping_sites = my
        .iter()
        .filter(|s| (0..=2).contains(&s.dosage) && partner_called.contains_key(&(s.contig.as_str(), s.position)))
        .count();

    let sample_a = group(my);
    let sample_b = group(partner);
    let mut lengths: BTreeMap<String, i32> = BTreeMap::new();
    for sample in [&sample_a, &sample_b] {
        for (chr, cg) in sample {
            let m = cg.positions.last().copied().unwrap_or(1);
            lengths.entry(chr.clone()).and_modify(|e| *e = (*e).max(m)).or_insert(m);
        }
    }
    let pairs: Vec<(&str, i32)> = lengths.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    let gmap = load_genetic_map(build, &pairs);
    let segments = PairwiseIbdDetector::new(config).detect_segments(&sample_a, &sample_b, &gmap);
    let summary = MatchSummary::from_segments(&segments);
    IbdComparison {
        summary,
        segments,
        overlapping_sites,
    }
}

/// Autosomal genotype concordance between two genotyped alignments: (matched, compared)
/// over sites both called (dosage within ploidy). ~1.0 ⇒ same individual; relatives lower.
fn genotype_concordance(a: &[SiteGenotype], b: &[SiteGenotype]) -> (i64, i64) {
    let called = |g: &SiteGenotype| (0..=g.ploidy as i32).contains(&g.dosage);
    let idx: HashMap<(&str, i64), i32> = b
        .iter()
        .filter(|g| called(g))
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();
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
    /// Inferred sex (M/F/U) from the `sex` artifact, if computed.
    pub sex: Option<String>,
    /// Mean read length (read-metrics artifact).
    pub mean_read_length: Option<f64>,
    /// % PF reads aligned (read-metrics artifact).
    pub pct_aligned: Option<f64>,
    /// Median insert size (read-metrics artifact).
    pub median_insert_size: Option<f64>,
    /// Number of structural variants called (`sv` artifact); `None` if not run.
    pub sv_count: Option<usize>,
    /// The coverage shown is a `partial` (lite sidecar) result — upgradeable by a deep walk.
    /// `false` when full (or no coverage yet).
    pub coverage_partial: bool,
}

/// One member row of a project's FTDNA-style Y-DNA STR overview: identity columns + the subject's
/// consensus STR marker values (normalized marker name → value) + terminal Y haplogroup. Members
/// without any STR profile are omitted by the query.
#[derive(Debug, Clone)]
pub struct ProjectStrMember {
    pub guid: SampleGuid,
    /// Display name (donor identifier).
    pub name: String,
    /// Kit / accession identifier, if recorded.
    pub kit: Option<String>,
    /// Loose origin (center name), if recorded.
    pub origin: Option<String>,
    /// Paternal-ancestor / free-text note (biosample description), if recorded.
    pub ancestor: Option<String>,
    /// Terminal Y haplogroup from the genome-level consensus, if placed.
    pub y_haplogroup: Option<String>,
    /// True when the haplogroup is SNP-backed (we only place from SNP evidence, so any placed
    /// haplogroup is confirmed); `false` leaves room for future STR-predicted labels.
    pub y_confirmed: bool,
    /// Highest STR panel/tier reached (e.g. "Y-111", "Alpha") — the "Test" column.
    pub test: Option<String>,
    /// Consensus STR values keyed by uppercase marker name (DYS393 → "13", DYS385 → "11-15").
    pub markers: std::collections::HashMap<String, String>,
}

/// One rendered cell of the project Y-STR chart: the marker value text plus its deviation from the
/// subgroup's modal value (the colour coding), precomputed so the UI does no per-frame work.
#[derive(Debug, Clone)]
pub struct StrChartCell {
    pub text: String,
    pub dev: navigator_domain::strchart::Deviation,
}

/// What a [`StrChartRow`] represents — drives how the UI styles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrRowKind {
    /// A subgroup banner (haplogroup heading); `cells` are empty.
    Group,
    /// The subgroup's per-marker MIN row.
    Min,
    /// The subgroup's per-marker MAX row.
    Max,
    /// The subgroup's per-marker MODE row.
    Mode,
    /// One member's marker values.
    Member,
}

/// One fully-prepared row of the project Y-STR overview, in display order. The whole table is
/// computed once off the UI thread; the renderer only iterates these.
#[derive(Debug, Clone)]
pub struct StrChartRow {
    pub kind: StrRowKind,
    /// Nesting depth for indentation (subgroups nest under their tree ancestors).
    pub depth: usize,
    /// Group banner text, member name, or the MIN/MAX/MODE label.
    pub label: String,
    /// Member kit/accession (member rows only).
    pub kit: String,
    /// Member terminal haplogroup (member rows only).
    pub haplogroup: String,
    /// True when the member's haplogroup is SNP-backed (drives green/red).
    pub confirmed: bool,
    /// Member's reached STR panel/tier (member rows only).
    pub test: String,
    /// Per-marker cells, aligned to [`ProjectStrChart::markers`]; empty for non-member rows that
    /// don't fill every column.
    pub cells: Vec<StrChartCell>,
}

/// The precomputed FTDNA-style project Y-DNA STR overview: marker column order + a flat, ordered
/// list of render-ready rows (subgroup banners, MIN/MAX/MODE, members), grouped by assigned Y
/// haplogroup and ordered by tree topology (basal → derived, children nested under parents).
#[derive(Debug, Clone, Default)]
pub struct ProjectStrChart {
    pub markers: Vec<String>,
    pub rows: Vec<StrChartRow>,
    pub member_count: usize,
    pub group_count: usize,
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
    pub sex_done: usize,
    pub metrics_done: usize,
    pub sv_done: usize,
    /// Per-sample failures (best-effort: one sample's error doesn't abort the rest).
    pub errors: Vec<String>,
}

/// What the deep analyze pass filled (or skipped as already-present) for one biosample.
/// `had_alignment` is false when the sample has no BAM-bearing alignment to walk — the
/// caller skips it without counting. Each `*_done` is true when that artifact is now present
/// (freshly computed or already cached); failures land in `errors`.
#[derive(Debug, Clone, Default)]
pub struct SampleAnalyzeOutcome {
    pub had_alignment: bool,
    pub coverage_done: bool,
    pub y_done: bool,
    pub sex_done: bool,
    pub metrics_done: bool,
    pub sv_done: bool,
    pub errors: Vec<String>,
}

/// Outcome of a BISDNA chromo2 Y-SNP import: the variant set created plus a per-category
/// tally so the UI/CLI can surface coverage and any names the dictionary couldn't place.
#[derive(Debug, Clone)]
pub struct BisdnaImportSummary {
    pub variant_set: VariantSet,
    /// Reference build the calls were emitted on (e.g. `"hs1"`).
    pub build: String,
    /// Total marker rows parsed from the file.
    pub total_markers: usize,
    /// Positive (derived) calls resolved to a locus and emitted as variant calls.
    pub derived_calls: usize,
    /// Negative (ancestral) markers — not variants, so not emitted (still counted).
    pub ancestral: usize,
    /// `no_call` markers (genotype `00`).
    pub no_call: usize,
    /// Back-mutated markers — flagged and excluded from placement.
    pub back_mutated: usize,
    /// Markers whose name was absent from the dictionary on this build (cannot be placed).
    pub unresolved: usize,
    /// A sample of unresolved names for diagnostics (capped).
    pub unresolved_names: Vec<String>,
    /// Positive calls whose genotype disagreed with the dictionary alleles on either strand
    /// (a QC signal — the call is still emitted, trusting the file's verdict).
    pub strand_mismatches: usize,
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
    /// Roll-up of the fast-path sidecar ingest across the imported samples.
    pub fast_path: FastPathSummary,
}

/// What the fast-path sidecar ingest filled across a project import (one tally per result
/// kind), so the import returns immediately with the report already populated.
#[derive(Debug, Clone, Default)]
pub struct FastPathSummary {
    /// Samples that had pipeline sidecars to ingest.
    pub samples_with_sidecars: usize,
    pub y_placed: usize,
    pub mt_placed: usize,
    pub sex_filled: usize,
    pub metrics_filled: usize,
    pub coverage_filled: usize,
    /// Per-sample ingest errors (`"<sample>: <detail>"`), non-fatal.
    pub errors: Vec<String>,
}

/// What [`App::ingest_sidecars`] managed to fill for one alignment.
#[derive(Debug, Clone, Default)]
pub struct SidecarIngest {
    pub y_haplogroup: Option<String>,
    pub mt_haplogroup: Option<String>,
    pub sex: Option<String>,
    pub read_metrics: bool,
    pub lite_coverage: bool,
    pub errors: Vec<String>,
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
        App {
            store,
            auth: Auth::new(),
            gateway,
        }
    }

    /// Open/create the workspace database and build the app.
    pub async fn open(path: &std::path::Path) -> Result<Self, AppError> {
        Ok(App::new(Store::open(path).await?))
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
        "sample_id,alignment_count,mean_coverage,median_coverage,pct_10x,pct_20x,callable_bases,\
         y_haplogroup,mt_haplogroup,sex,mean_read_length,pct_aligned,median_insert_size,sv_count\n",
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
        s.push(',');
        s.push_str(&field(r.sex.as_deref().unwrap_or("")));
        s.push(',');
        s.push_str(&num(r.mean_read_length));
        s.push(',');
        s.push_str(&num(r.pct_aligned));
        s.push(',');
        s.push_str(&num(r.median_insert_size));
        s.push(',');
        s.push_str(&r.sv_count.map(|v| v.to_string()).unwrap_or_default());
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod placement_tests {
    use super::SourceType;
    use super::{
        assemble_assignment, assemble_assignment_robust, contributing_subdirs, pool_votes, snp_obs_from_assignment,
        strand_reconcile_to_tree, support_backoff_terminal,
    };
    use navigator_analysis::haplo::parse_ftdna_json;
    use std::collections::HashMap;
    use std::path::PathBuf;

    #[test]
    fn contributing_subdirs_flags_multi_sample_folders() {
        let root = std::path::Path::new("/data/FTDNA");
        // A single sample's folder: a top-level results CSV + one kit subdir with the BAM.
        let one_sample = [
            PathBuf::from("/data/FTDNA/42048/42048_YDNA_DYS_Results.csv"),
            PathBuf::from("/data/FTDNA/42048/2691/abc.bam"),
        ];
        // Relative to the *sample* folder, only the kit dir contributes → not multi-sample.
        let sample_root = std::path::Path::new("/data/FTDNA/42048");
        assert_eq!(contributing_subdirs(sample_root, &one_sample).len(), 1);

        // Relative to the parent download root, each sample folder contributes → multi-sample.
        let many = [
            PathBuf::from("/data/FTDNA/42048/42048_YDNA.csv"),
            PathBuf::from("/data/FTDNA/166433/166433_YDNA.csv"),
            PathBuf::from("/data/FTDNA/166433/9369/x.bam"),
        ];
        let dirs = contributing_subdirs(root, &many);
        assert_eq!(dirs.len(), 2);
        assert!(dirs.contains("42048") && dirs.contains("166433"));

        // Files sitting directly in the picked folder don't count as a subdir.
        let flat = [PathBuf::from("/data/FTDNA/a.csv"), PathBuf::from("/data/FTDNA/b.bam")];
        assert!(contributing_subdirs(root, &flat).is_empty());
    }

    // A six-node spine for the back-off tests: root(1) → A(2,@146) → B(3,@263) → C(4,@750)
    // → D(5,@1000) → F(6,@1100), one defining SNP per node.
    const SPINE6: &str = r#"{ "allNodes": {
      "1": {"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
      "2": {"haplogroupId":2,"name":"A","isRoot":false,"variants":[{"variant":"a","position":146,"ancestral":"A","derived":"G"}],"children":[3]},
      "3": {"haplogroupId":3,"name":"B","isRoot":false,"variants":[{"variant":"b","position":263,"ancestral":"A","derived":"G"}],"children":[4]},
      "4": {"haplogroupId":4,"name":"C","isRoot":false,"variants":[{"variant":"c","position":750,"ancestral":"C","derived":"T"}],"children":[5]},
      "5": {"haplogroupId":5,"name":"D","isRoot":false,"variants":[{"variant":"d","position":1000,"ancestral":"G","derived":"A"}],"children":[6]},
      "6": {"haplogroupId":6,"name":"F","isRoot":false,"variants":[{"variant":"f","position":1100,"ancestral":"C","derived":"T"}],"children":[]}
    }}"#;

    /// The consensus-profile observations come from the placed **lineage** (the derived mutations
    /// root→terminal), not the untaken child branches. Regression for the all-no-call profile: a
    /// clean single-source placement at C must yield Derived obs for a/b/c, and must NOT carry the
    /// child branch D's defining SNP `d` (which is no-call because descent stopped).
    #[test]
    fn profile_obs_follow_the_lineage_not_the_untaken_children() {
        let tree = parse_ftdna_json(SPINE6).unwrap();
        // Derived A+B+C; D/F not covered → terminal C(4), child D(5,@1000) is no-call.
        let calls: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'T')].into_iter().collect();
        let assignment = assemble_assignment(&tree, &calls);
        // The lineage carries a,b,c (derived); the child branch D carries d.
        assert!(
            assignment.branches.iter().any(|b| b.snps.iter().any(|s| s.name == "d")),
            "D is a child branch"
        );
        let obs = snp_obs_from_assignment(&assignment, true);
        let mut names: Vec<&str> = obs.iter().map(|o| o.name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(
            names,
            ["a", "b", "c"],
            "obs are the lineage mutations, not the untaken child d"
        );
    }

    /// The parsimony back-off trims a net-contradicted deep tail (the sparse-panel / aDNA
    /// over-deepening) to the node where running (derived − ancestral) support peaks, but keeps a
    /// clean deep path and tolerates a lone contradiction outweighed by deeper support.
    #[test]
    fn support_backoff_trims_net_negative_tail_but_keeps_supported_depth() {
        let tree = parse_ftdna_json(SPINE6).unwrap();
        // Derived A+B (peak at B), then below B: ancestral@750, contradiction@1000 (G≠der A),
        // a lone derived@1100 — tail net −1. Should back off F(6) → B(3).
        let sparse: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'C'), (1000, 'G'), (1100, 'T')]
            .into_iter()
            .collect();
        assert_eq!(
            support_backoff_terminal(&tree, &sparse, 6),
            3,
            "net-negative tail trimmed to B"
        );

        // A clean fully-derived path keeps the deepest terminal F.
        let clean: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A'), (1100, 'T')]
            .into_iter()
            .collect();
        assert_eq!(
            support_backoff_terminal(&tree, &clean, 6),
            6,
            "clean path keeps the terminal"
        );

        // A lone contradiction (@750) outweighed by deeper derived calls still reaches F.
        let recovered: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'C'), (1000, 'A'), (1100, 'T')]
            .into_iter()
            .collect();
        assert_eq!(
            support_backoff_terminal(&tree, &recovered, 6),
            6,
            "deeper support recovers depth"
        );
    }

    /// Genome-level pooling: a sparse source that alone stops shallow, combined with a dense source
    /// that confirms the deep branches, places the *pooled* call set (by position on one tree) at
    /// the deep terminal — a per-run sparse call no longer drags the consensus shallow.
    #[test]
    fn pooled_consensus_places_deeper_than_a_sparse_source_alone() {
        let tree = parse_ftdna_json(SPINE6).unwrap();
        // Sparse chip: only the shallow A SNP derived (146). Alone → A(2).
        let sparse: HashMap<i64, char> = [(146, 'G')].into_iter().collect();
        let sparse_only = assemble_assignment_robust(&tree, &sparse);
        assert_eq!(sparse_only.ranked.first().unwrap().name, "A");

        // Dense WGS: every spine SNP derived. Pool the two by position → the deep terminal F(6).
        let dense: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A'), (1100, 'T')]
            .into_iter()
            .collect();
        let pooled = pool_votes(&[(SourceType::Chip, sparse), (SourceType::WgsShortRead, dense)]);
        let placed = assemble_assignment(&tree, &pooled);
        assert_eq!(
            placed.ranked.first().unwrap().name,
            "F",
            "pooled evidence reaches the deep terminal"
        );
    }

    /// A higher-weight source wins the per-position vote: WGS (0.85) derived outvotes a Chip
    /// ancestral call at the same SNP.
    #[test]
    fn pool_vote_prefers_the_higher_weight_source() {
        let wgs: HashMap<i64, char> = [(750, 'T')].into_iter().collect(); // derived
        let chip: HashMap<i64, char> = [(750, 'C')].into_iter().collect(); // ancestral
        let pooled = pool_votes(&[(SourceType::WgsShortRead, wgs), (SourceType::Chip, chip)]);
        assert_eq!(pooled.get(&750), Some(&'T'), "WGS derived outweighs chip ancestral");
    }

    /// Chip alleles on the tree's opposite strand are flipped to the matching ancestral/derived
    /// allele; in-tree matches and out-of-tree positions are untouched. A flipped derived call
    /// then places as deep as the forward one would.
    #[test]
    fn strand_reconcile_flips_only_opposite_strand_calls() {
        let tree = parse_ftdna_json(TREE).unwrap();
        // 146 der=G observed as C (= complement of G) → flips to G; 263 der=G observed forward;
        // 999 absent from the tree → passthrough unchanged.
        let calls: HashMap<i64, char> = [(146, 'C'), (263, 'G'), (999, 'C')].into_iter().collect();
        let fixed = strand_reconcile_to_tree(&tree, calls);
        assert_eq!(fixed[&146], 'G', "complement matched the derived allele");
        assert_eq!(fixed[&263], 'G', "already matched → unchanged");
        assert_eq!(fixed[&999], 'C', "not in the tree → passthrough");

        // The reconciled calls place to B (derived at 146 + 263), same as forward-strand input.
        assert_eq!(
            assemble_assignment_robust(
                &tree,
                &strand_reconcile_to_tree(&tree, [(146, 'C'), (263, 'G')].into_iter().collect())
            )
            .ranked
            .first()
            .unwrap()
            .name,
            "B"
        );
    }

    // root → A(146) → B(263) → C(750) → D(1000). A single defining SNP per node.
    const TREE: &str = r#"{ "allNodes": {
      "1": {"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
      "2": {"haplogroupId":2,"name":"A","isRoot":false,"variants":[{"variant":"a","position":146,"ancestral":"A","derived":"G"}],"children":[3]},
      "3": {"haplogroupId":3,"name":"B","isRoot":false,"variants":[{"variant":"b","position":263,"ancestral":"A","derived":"G"}],"children":[4]},
      "4": {"haplogroupId":4,"name":"C","isRoot":false,"variants":[{"variant":"c","position":750,"ancestral":"C","derived":"T"}],"children":[5]},
      "5": {"haplogroupId":5,"name":"D","isRoot":false,"variants":[{"variant":"d","position":1000,"ancestral":"G","derived":"A"}],"children":[]}
    }}"#;

    /// A deep lineage with a single stray ancestral call on a backbone node (C) — the sparse-
    /// chip failure mode. Strict selection vetoes the whole lineage and stops shallow (B);
    /// robust selection trusts the proportional top and reaches the deep terminal (D).
    #[test]
    fn robust_selection_survives_a_backbone_contradiction() {
        let tree = parse_ftdna_json(TREE).unwrap();
        // Derived at 146, 263, 1000; but ANCESTRAL (C) at 750 — a lone contradiction on node C.
        let calls: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'C'), (1000, 'A')].into_iter().collect();

        let strict = assemble_assignment(&tree, &calls);
        let robust = assemble_assignment_robust(&tree, &calls);

        // Strict stops above the contradicted node C → terminal B (shallow).
        assert_eq!(strict.ranked.first().unwrap().name, "B");
        // Robust reaches the genuine deep terminal D despite the stray ancestral.
        assert_eq!(robust.ranked.first().unwrap().name, "D");
    }

    /// With a clean lineage (no contradiction) both selectors agree on the deep terminal.
    #[test]
    fn robust_and_strict_agree_when_path_is_clean() {
        let tree = parse_ftdna_json(TREE).unwrap();
        let calls: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A')].into_iter().collect();
        assert_eq!(assemble_assignment(&tree, &calls).ranked.first().unwrap().name, "D");
        assert_eq!(
            assemble_assignment_robust(&tree, &calls).ranked.first().unwrap().name,
            "D"
        );
    }

    /// The GVCF fast path reconstructs exactly the `calls` a pileup would yield. A fully
    /// derived path (every defining SNP a variant) places to the deep terminal D.
    #[test]
    fn gvcf_derived_path_places_deep() {
        use navigator_analysis::gvcf;
        let tree = parse_ftdna_json(TREE).unwrap();
        let mut called = gvcf::CalledBases::default();
        called
            .variant_bases
            .extend([(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A')]);
        called.callable.extend([146, 263, 750, 1000]);
        // Reference bases are irrelevant here (every site is a variant).
        let calls = gvcf::assemble_calls(&called, &HashMap::new());
        let expected: HashMap<i64, char> = [(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A')].into_iter().collect();
        assert_eq!(calls, expected);
        assert_eq!(assemble_assignment(&tree, &calls).ranked.first().unwrap().name, "D");
    }

    /// A hom-ref (callable, no variant) tree SNP reconstructs as the **reference base** — which
    /// on a real reference can be the *derived* allele (CHM13 Y = J1). Here position 750's
    /// reference base is the derived T, so node C is supported and placement reaches D — the
    /// exact case the old "assume ancestral" logic got wrong (it stopped at B).
    #[test]
    fn gvcf_homref_site_takes_reference_base_not_ancestral() {
        use navigator_analysis::gvcf;
        let tree = parse_ftdna_json(TREE).unwrap();
        let mut called = gvcf::CalledBases::default();
        called.variant_bases.extend([(146, 'G'), (263, 'G'), (1000, 'A')]);
        called.callable.extend([146, 263, 750, 1000]); // 750 hom-ref → its reference base
                                                       // The reference carries the *derived* T at 750 (shared backbone the sample also has).
        let ref_base: HashMap<i64, char> = [(750, 'T')].into_iter().collect();
        let calls = gvcf::assemble_calls(&called, &ref_base);
        assert_eq!(
            calls.get(&750),
            Some(&'T'),
            "hom-ref site takes the reference base (derived here)"
        );
        assert_eq!(assemble_assignment(&tree, &calls).ranked.first().unwrap().name, "D");
    }

    /// Lifted assembly maps GVCF observations back to tree positions, reverse-complementing a
    /// minus-strand lift; hom-ref lifted sites take the reference base at the lifted position.
    #[test]
    fn lifted_assembly_maps_back_and_revcomps() {
        use navigator_analysis::gvcf;
        use navigator_refgenome::LiftedPos;
        let mut called = gvcf::CalledBases::default();
        called.variant_bases.insert(500, 'G'); // tree 146 → derived G (forward)
        called.callable.extend([500, 900]); // 900 hom-ref → reference base, minus strand
        let ref_base: HashMap<i64, char> = [(900, 'C')].into_iter().collect();
        let lifted = vec![
            LiftedPos {
                tree_pos: 146,
                contig: "chrM".into(),
                pos: 500,
                reverse: false,
            },
            LiftedPos {
                tree_pos: 263,
                contig: "chrM".into(),
                pos: 900,
                reverse: true,
            },
        ];
        let calls = super::assemble_calls_lifted(&called, &lifted, &ref_base);
        assert_eq!(calls.get(&146), Some(&'G'));
        assert_eq!(calls.get(&263), Some(&'G'), "minus-strand reference C → complement G");
    }
}

#[cfg(test)]
mod publish_tests {
    use super::*;
    use navigator_domain::workspace::NewSequenceRun;
    use navigator_store::Store;

    /// The published sequence-run record carries the inferred `instrumentId` (camelCase) so the
    /// AppView can crowd-source the instrument→lab map. Regression guard: this field was hardcoded
    /// to `None` while the lab inference was being restored.
    #[tokio::test]
    async fn sequence_run_record_publishes_instrument_id() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: Some("NovaSeq".into()),
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await
            .unwrap();
        sequence_run::set_library_stats(
            app.store.pool(),
            run.id,
            Some("A00182"),
            None,
            None,
            None,
            Some("H5WLTDMXX"),
        )
        .await
        .unwrap();
        let reloaded = sequence_run::get(app.store.pool(), run.id).await.unwrap().unwrap();

        let value = app.sequence_run_record(&reloaded).await.unwrap();
        assert_eq!(value.get("instrumentId").and_then(|v| v.as_str()), Some("A00182"));
        assert_eq!(value.get("$type").and_then(|v| v.as_str()), Some(NS_SEQUENCERUN));
    }
}

#[cfg(test)]
mod ymatch_tests {
    use super::*;
    use navigator_domain::strprofile::StrMarker;
    use navigator_store::Store;

    async fn seed_str(app: &App, guid: SampleGuid, markers: &[(&str, &str)]) {
        let new = NewStrProfile {
            biosample_guid: guid,
            panel_name: "Y-37".into(),
            provider: Some("FTDNA".into()),
            source: None,
            markers: markers
                .iter()
                .map(|(m, v)| StrMarker {
                    marker: (*m).into(),
                    value: (*v).into(),
                })
                .collect(),
        };
        str_profile::create(app.store.pool(), &new).await.unwrap();
    }

    /// End-to-end app orchestration for STR-only subjects: enumerate the workspace, assemble each
    /// subject's match profile from cached data, and rank by ascending Y-STR genetic distance. (No
    /// Y tree / SNP profiles → exercises the offline, graceful-degradation path.)
    #[tokio::test]
    async fn y_matches_ranks_str_only_by_distance() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let q = app.add_biosample(None, "Query", None, None).await.unwrap();
        let near = app.add_biosample(None, "Near", None, None).await.unwrap();
        let mid = app.add_biosample(None, "Mid", None, None).await.unwrap();
        let far = app.add_biosample(None, "Far", None, None).await.unwrap();

        let base = [("DYS393", "13"), ("DYS390", "24"), ("DYS19", "14")];
        seed_str(&app, q.guid, &base).await;
        seed_str(&app, near.guid, &base).await; // GD 0
        seed_str(&app, mid.guid, &[("DYS393", "13"), ("DYS390", "25"), ("DYS19", "14")]).await; // GD 1
        seed_str(&app, far.guid, &[("DYS393", "12"), ("DYS390", "26"), ("DYS19", "14")]).await; // GD 2

        let matches = app.y_matches(q.guid, None).await.unwrap();
        assert_eq!(matches.len(), 3, "query is excluded; three candidates ranked");
        assert_eq!(
            matches.iter().map(|m| m.donor.as_str()).collect::<Vec<_>>(),
            ["Near", "Mid", "Far"]
        );
        assert_eq!(matches[0].str_gd, Some(0));
        assert_eq!(matches[2].str_gd, Some(2));
        assert!(matches.iter().all(|m| m.signal == YSignal::Str));
        // STR TMRCA present and monotonic with distance.
        assert!(matches[0].str_tmrca.is_some());
        assert!(matches[2].str_tmrca.unwrap().generations > matches[0].str_tmrca.unwrap().generations);
    }

    /// A subject with no comparable Y data is dropped; an empty workspace yields no matches.
    #[tokio::test]
    async fn y_matches_drops_incomparable_and_self() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let q = app.add_biosample(None, "Query", None, None).await.unwrap();
        seed_str(&app, q.guid, &[("DYS393", "13")]).await;
        // A subject with no STR / Y data at all.
        let _empty = app.add_biosample(None, "Empty", None, None).await.unwrap();
        // A subject whose markers don't overlap the query's.
        let other = app.add_biosample(None, "Other", None, None).await.unwrap();
        seed_str(&app, other.guid, &[("DYS999", "10")]).await;

        let matches = app.y_matches(q.guid, None).await.unwrap();
        assert!(matches.is_empty(), "no comparable candidates");
    }
}

#[cfg(test)]
mod sync_pull_tests {
    use super::*;
    use navigator_store::Store;

    /// PULL applies a remote biosample record's editable summary fields (sex / center) onto the local
    /// subject. (The fed record is PII-free, so only these fields are present to apply.)
    #[tokio::test]
    async fn apply_remote_updates_biosample_fields() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, Some("F".into())).await.unwrap();
        let value = serde_json::json!({ "sex": "M", "center_name": "LabX" });
        app.apply_remote(NS_BIOSAMPLE, &format!("biosample:{}", b.guid), &value)
            .await
            .unwrap();
        let updated = app
            .list_all_biosamples()
            .await
            .unwrap()
            .into_iter()
            .find(|x| x.guid == b.guid)
            .unwrap();
        assert_eq!(updated.sex.as_deref(), Some("M"));
        assert_eq!(updated.center_name.as_deref(), Some("LabX"));
        assert_eq!(
            updated.donor_identifier, "S1",
            "identity (donor_identifier) is preserved — not in the PII-free record"
        );
    }

    /// A derived-summary collection is a no-op on apply (recomputed locally, never overwritten).
    #[tokio::test]
    async fn apply_remote_derived_is_noop() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        // No panic / no error for a collection we only track.
        app.apply_remote(NS_ALIGNMENT, "alignment:1", &serde_json::json!({}))
            .await
            .unwrap();
    }
}

#[cfg(test)]
mod ibd_attest_tests {
    use super::*;
    use navigator_analysis::ibd::{IbdSegment, MatchSummary};
    use navigator_sync::DeviceKey;

    fn summary(cm: f64) -> MatchSummary {
        MatchSummary::from_segments(&[IbdSegment {
            chromosome: "chr1".into(),
            start_position: 1,
            end_position: 10_000_000,
            length_cm: cm,
            snp_count: Some(500),
            is_half_identical: None,
        }])
    }

    /// A signed attestation verifies with the same code path the AppView runs; tampering breaks it.
    #[test]
    fn attestation_signs_and_verifies() {
        let key = DeviceKey::generate();
        let mut att = IbdAttestation::unsigned(
            "exchange:r1",
            "sess-1",
            key.did_key(),
            Some("bio-a".into()),
            Some("bio-b".into()),
            &summary(42.0),
            "2026-06-17T00:00:00Z",
        );
        att.signature = key.sign(&att.canonical());
        att.signing_public_key = key.did_key();

        assert!(
            du_atproto::verify_did_key(&att.signing_public_key, att.canonical().as_bytes(), &att.signature).is_ok()
        );
        // Tamper a signed field → canonical changes → verification fails.
        let mut bad = att.clone();
        bad.total_shared_cm = 999.0;
        assert!(
            du_atproto::verify_did_key(&bad.signing_public_key, bad.canonical().as_bytes(), &att.signature).is_err()
        );
    }

    /// Two peers computing the same summary produce the same agreement hash; different summaries don't.
    #[test]
    fn summary_hash_drives_agreement() {
        use navigator_analysis::ibd_attest::summary_hash;
        assert_eq!(summary_hash(&summary(42.0)), summary_hash(&summary(42.0)));
        assert_ne!(summary_hash(&summary(42.0)), summary_hash(&summary(7.0)));
    }

    /// An exchange result persists and reads back per subject (the UI's saved-results path).
    #[tokio::test]
    async fn exchange_result_persists_and_lists() {
        use navigator_store::Store;
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, None).await.unwrap();
        let session = EstablishedSession {
            session_id: "sess-x".into(),
            partner_did: "did:key:zB".into(),
            key: [0u8; 32],
        };
        let result = IbdExchangeResult {
            summary: summary(75.0),
            segments: vec![],
            overlapping_sites: 100,
            my_attestation: IbdAttestation::unsigned(
                "exchange:r",
                "sess-x",
                "did:key:zA",
                Some(b.guid.to_string()),
                Some("bio-B".into()),
                &summary(75.0),
                "t",
            ),
            partner_attestation: IbdAttestation::unsigned(
                "exchange:r",
                "sess-x",
                "did:key:zB",
                Some("bio-B".into()),
                Some(b.guid.to_string()),
                &summary(75.0),
                "t",
            ),
            agreed: true,
        };
        app.record_ibd_exchange(b.guid, &session, "exchange:r", &result)
            .await
            .unwrap();
        let rows = app.list_ibd_exchanges_for_subject(b.guid).await.unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].total_shared_cm, 75.0);
        assert!(rows[0].agreed);
        assert_eq!(rows[0].partner_did, "did:key:zB");
    }
}

#[cfg(test)]
mod ibd_federated_tests {
    use super::*;
    use navigator_sync::DeviceKey;

    #[test]
    fn ibd_poll_message_signs_and_verifies() {
        // The exact canonical bytes a device signs for the suggestions poll, end-to-end
        // verifiable by the AppView's own verifier (proves the wire contract).
        let key = DeviceKey::generate();
        let msg = format!("ibd-poll\n{}\n{}", "did:plc:abc123", "1718000000");
        assert_eq!(msg, "ibd-poll\ndid:plc:abc123\n1718000000");
        let sig = key.sign(&msg);
        assert!(du_atproto::verify_did_key(&key.did_key(), msg.as_bytes(), &sig).is_ok());
    }

    #[test]
    fn ibd_introduce_message_shape() {
        let msg = format!("ibd-introduce\n{}\n{}", "did:plc:abc123", "sample-xyz");
        assert_eq!(msg, "ibd-introduce\ndid:plc:abc123\nsample-xyz");
    }

    #[test]
    fn query_sig_is_url_encoded() {
        // STANDARD base64 (`+` `/` `=`) must be percent-escaped in the GET query string.
        let req = reqwest::Client::new()
            .get("http://x/api")
            .query(&[("sig", "a+b/c=")])
            .build()
            .unwrap();
        let q = req.url().query().unwrap();
        assert!(q.contains("a%2Bb%2Fc%3D"), "sig not URL-encoded: {q}");
    }

    #[test]
    fn parse_suggestions_appview_snake_case_string_signals() {
        // The exact shape the AppView emits: snake_case keys, metadata.signals as a
        // plain string array.
        let body = serde_json::json!({
            "items": [{
                "suggested_sample_guid": "g1",
                "suggestion_type": "POPULATION_OVERLAP",
                "score": 0.82,
                "metadata": { "signals": ["POPULATION_OVERLAP", "HAPLOGROUP"] }
            }]
        });
        let out = parse_ibd_suggestions(&body);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].suggested_sample_guid, "g1");
        assert_eq!(out[0].suggestion_type, "POPULATION_OVERLAP");
        assert!((out[0].score - 0.82).abs() < 1e-9);
        assert_eq!(
            out[0].signals,
            vec!["POPULATION_OVERLAP".to_string(), "HAPLOGROUP".to_string()]
        );
    }

    #[test]
    fn parse_suggestions_camel_case_and_object_signals_tolerated() {
        let body = serde_json::json!({
            "suggestions": [{
                "suggestedSampleGuid": "g2",
                "type": "SHARED_MATCH",
                "score": 1.0,
                "signals": { "sharedMatches": 3.0 }
            }]
        });
        let out = parse_ibd_suggestions(&body);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].suggested_sample_guid, "g2");
        assert_eq!(out[0].suggestion_type, "SHARED_MATCH");
        assert_eq!(out[0].signals, vec!["sharedMatches".to_string()]);
    }

    #[test]
    fn parse_suggestions_empty_or_malformed_is_empty() {
        assert!(parse_ibd_suggestions(&serde_json::json!({})).is_empty());
        assert!(parse_ibd_suggestions(&serde_json::json!({ "items": "nope" })).is_empty());
    }
}

#[cfg(test)]
mod export_tests {
    use super::*;

    #[test]
    fn primary_contig_filter() {
        for ok in ["chr1", "1", "chr22", "22", "chrX", "X", "chrY", "Y", "chrM", "MT"] {
            assert!(is_primary_contig(ok), "{ok} should be primary");
        }
        for no in [
            "chr23",
            "chrUn_KI270302v1",
            "chr1_KI270706v1_random",
            "HLA-A",
            "GL000220.1",
            "",
        ] {
            assert!(!is_primary_contig(no), "{no} should not be primary");
        }
    }

    #[test]
    fn diploid_vcf_export_metadata() {
        let r = ExportRequest::DiploidVcf(7);
        assert_eq!(r.extension(), "vcf");
        assert_eq!(r.default_filename(), "diploid_variants_7.vcf");
        assert!(r.label().contains("VCF"));
    }
}

#[cfg(test)]
mod ibd_tests {
    use super::*;

    fn sg(contig: &str, pos: i64, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: contig.into(),
            position: pos,
            reference_allele: "A".into(),
            alternate_allele: "G".into(),
            ploidy: 2,
            dosage,
            gq: 0,
            depth: 0,
            ref_depth: 0,
            alt_depth: 0,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        }
    }

    #[test]
    fn overlapping_sites_counts_both_called_intersection() {
        let a = vec![sg("chr1", 100, 0), sg("chr1", 200, 1), sg("chr1", 300, -1)]; // 300 no-call
        let b = vec![
            sg("chr1", 100, 2),
            sg("chr1", 200, 1),
            sg("chr1", 300, 0),
            sg("chr1", 400, 0),
        ];
        // Shared & called in both: 100, 200 (300 is a no-call in a; 400 absent in a).
        assert_eq!(overlapping_called_sites(&a, &b), 2);
        assert_eq!(overlapping_called_sites(&a, &[]), 0);
    }
}

#[cfg(test)]
mod outbox_tests {
    use super::*;

    #[test]
    fn backoff_doubles_per_attempt_and_caps_at_one_hour() {
        assert_eq!(backoff_secs(1), 120); // 2 min
        assert_eq!(backoff_secs(2), 240); // 4 min
        assert_eq!(backoff_secs(3), 480); // 8 min
        assert_eq!(backoff_secs(5), 1920); // 32 min
        assert_eq!(backoff_secs(6), 3600); // 64 min → capped at 1 h
        assert_eq!(backoff_secs(40), 3600); // huge attempt → still capped, no overflow
        assert_eq!(backoff_secs(0), 60); // defensive: 1 min
    }

    #[tokio::test]
    async fn publish_while_signed_out_is_not_authenticated_and_queues_nothing() {
        // Use the in-memory keychain so this never touches the OS keychain (no prompts) and is
        // hermetic regardless of an ambient session (e.g. a GUI signed in on this machine).
        navigator_sync::TokenStore::use_in_memory_for_tests();
        let app = App::new(Store::open_in_memory().await.unwrap());
        assert!(matches!(app.publish_coverage(1).await, Err(AppError::NotAuthenticated)));
        // No account → nothing enqueued, and the accessors degrade gracefully.
        assert_eq!(app.outbox_pending_count().await.unwrap(), 0);
        assert!(app.outbox_entries().await.unwrap().is_empty());
        assert!(app.sync_history(10).await.unwrap().is_empty());
        // Draining without an account is a harmless no-op.
        let outcome = app.drain_outbox().await.unwrap();
        assert_eq!(outcome.pending, 0);
        assert!(outcome.published.is_empty());
    }
}

#[cfg(test)]
mod settings_tests {
    use super::*;

    #[test]
    fn y_provider_precedence_env_then_settings_then_default() {
        // env wins even when settings disagree
        assert!(matches!(
            resolve_y_provider(Some("ftdna"), Some("decodingus")),
            YTreeProvider::Ftdna
        ));
        assert!(matches!(
            resolve_y_provider(Some("decodingus"), Some("ftdna")),
            YTreeProvider::DecodingUs
        ));
        // settings used when env absent
        assert!(matches!(resolve_y_provider(None, Some("ftdna")), YTreeProvider::Ftdna));
        // default when neither
        assert!(matches!(resolve_y_provider(None, None), YTreeProvider::DecodingUs));
        // unrecognized value falls back to default
        assert!(matches!(
            resolve_y_provider(Some("bogus"), None),
            YTreeProvider::DecodingUs
        ));
    }

    #[test]
    fn appview_url_precedence_and_normalization() {
        assert_eq!(
            resolve_appview_url(Some("https://av.example/".into()), Some("http://x".into())),
            "https://av.example"
        );
        assert_eq!(
            resolve_appview_url(None, Some("http://host:9000".into())),
            "http://host:9000"
        );
        assert_eq!(resolve_appview_url(None, None), "http://localhost:9000");
        // blank values are ignored (fall through to default)
        assert_eq!(resolve_appview_url(Some("".into()), None), "http://localhost:9000");
    }

    #[test]
    fn app_settings_serde_round_trip_and_defaults() {
        let s = AppSettings {
            y_tree_provider: Some("ftdna".into()),
            appview_url: Some("https://av.example".into()),
            tree_ttl_days: Some(3),
            theme: Some("light".into()),
            prompt_before_download: Some(false),
            ui_scale: Some(1.5),
            ui_mode: Some("simple".into()),
        };
        let json = serde_json::to_string(&s).unwrap();
        assert_eq!(serde_json::from_str::<AppSettings>(&json).unwrap(), s);

        // Missing/partial fields default to None (forward/backward compatible).
        let partial: AppSettings = serde_json::from_str(r#"{"appview_url":"http://h"}"#).unwrap();
        assert_eq!(partial.appview_url.as_deref(), Some("http://h"));
        assert_eq!(partial.y_tree_provider, None);
        assert_eq!(
            AppSettings::default(),
            serde_json::from_str::<AppSettings>("{}").unwrap()
        );
    }

    #[test]
    fn ui_mode_resolution_env_over_settings() {
        use super::{resolve_ui_mode, UiMode};
        // env wins
        assert_eq!(
            resolve_ui_mode(Some("advanced"), Some("simple")),
            Some(UiMode::Advanced)
        );
        // settings used when no env
        assert_eq!(resolve_ui_mode(None, Some("Simple")), Some(UiMode::Simple));
        // neither set → None (UI applies its heuristic)
        assert_eq!(resolve_ui_mode(None, None), None);
        // unrecognized tokens are ignored
        assert_eq!(resolve_ui_mode(Some("expert"), Some("simple")), Some(UiMode::Simple));
        assert_eq!(resolve_ui_mode(Some("expert"), None), None);
    }
}

#[cfg(test)]
mod import_tests {
    use super::{artifact_is_fresh, collect_data_files, file_signature, is_recognized_data_file};
    use std::path::Path;

    #[test]
    fn artifact_freshness_only_rejects_a_known_mismatch() {
        assert!(artifact_is_fresh(Some("100:5"), Some("100:5")), "matching sig → fresh");
        assert!(
            !artifact_is_fresh(Some("100:5"), Some("200:5")),
            "changed mtime → stale"
        );
        assert!(!artifact_is_fresh(Some("100:5"), Some("100:9")), "changed size → stale");
        assert!(
            artifact_is_fresh(None, Some("100:5")),
            "legacy row (no stored sig) → trusted"
        );
        assert!(
            artifact_is_fresh(Some("100:5"), None),
            "source gone (no current sig) → trusted"
        );
    }

    #[test]
    fn file_signature_changes_when_content_grows() {
        let p = std::env::temp_dir().join(format!("dun-sig-{}.bin", std::process::id()));
        std::fs::write(&p, b"abc").unwrap();
        let s1 = file_signature(&p).expect("sig");
        std::fs::write(&p, b"abcdef").unwrap(); // size changes → signature changes
        let s2 = file_signature(&p).expect("sig");
        assert_ne!(s1, s2);
        assert!(file_signature(Path::new("/no/such/file/xyz")).is_none());
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn recognizes_data_extensions() {
        for ok in [
            "x.bam", "x.cram", "x.vcf", "x.vcf.gz", "x.fasta", "x.fa", "x.csv", "x.tsv", "x.txt",
        ] {
            assert!(is_recognized_data_file(Path::new(ok)), "{ok} should be recognized");
        }
        for no in ["x.png", "x.pdf", "x", "x.bai", "x.crai"] {
            assert!(!is_recognized_data_file(Path::new(no)), "{no} should not be recognized");
        }
    }

    #[test]
    fn collect_walks_dir_and_filters() {
        // A temp tree: top-level a.bam + ignore.png, nested sub/b.vcf + sub/c.txt.
        let dir = std::env::temp_dir().join(format!("dun-import-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        for f in ["a.bam", "ignore.png"] {
            std::fs::write(dir.join(f), b"x").unwrap();
        }
        std::fs::write(dir.join("sub/b.vcf"), b"x").unwrap();
        std::fs::write(dir.join("sub/c.txt"), b"x").unwrap();

        let mut out = Vec::new();
        collect_data_files(&dir, &mut out, 0);
        let names: std::collections::BTreeSet<String> = out
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
            .collect();
        assert_eq!(
            names,
            ["a.bam", "b.vcf", "c.txt"].iter().map(|s| s.to_string()).collect()
        );

        // A single recognized file yields itself; an unrecognized file yields nothing.
        let mut one = Vec::new();
        collect_data_files(&dir.join("a.bam"), &mut one, 0);
        assert_eq!(one.len(), 1);
        let mut none = Vec::new();
        collect_data_files(&dir.join("ignore.png"), &mut none, 0);
        assert!(none.is_empty());

        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
mod seed_tests {
    use super::{seed_assets_from, SeedSummary};

    #[test]
    fn seeds_missing_files_and_never_overwrites() {
        let base = std::env::temp_dir().join(format!("dun-seed-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let src = base.join("bundle");
        let dest = base.join("cache");
        std::fs::create_dir_all(&src).unwrap();
        std::fs::write(src.join("a.bin"), b"alpha").unwrap();
        std::fs::write(src.join("b.bin"), b"bravo").unwrap();
        std::fs::write(src.join("manifest.json"), b"{}").unwrap();
        // `b.bin` already exists in the cache (e.g. a CDN-refreshed copy) — must be preserved.
        std::fs::create_dir_all(&dest).unwrap();
        std::fs::write(dest.join("b.bin"), b"NEWER").unwrap();

        let s = seed_assets_from(&src, &dest).unwrap();
        assert_eq!(s, SeedSummary { copied: 2, skipped: 1 }); // a.bin + manifest copied; b.bin skipped
        assert_eq!(std::fs::read(dest.join("a.bin")).unwrap(), b"alpha");
        assert_eq!(std::fs::read(dest.join("b.bin")).unwrap(), b"NEWER"); // not overwritten

        // Idempotent: a second run copies nothing.
        let again = seed_assets_from(&src, &dest).unwrap();
        assert_eq!(again, SeedSummary { copied: 0, skipped: 3 });

        // A missing bundle dir is a harmless no-op.
        let absent = seed_assets_from(&base.join("nope"), &dest).unwrap();
        assert_eq!(absent, SeedSummary::default());

        let _ = std::fs::remove_dir_all(&base);
    }
}

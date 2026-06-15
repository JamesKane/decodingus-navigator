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
use navigator_analysis::caller::{self, HaploidCallerParams, SiteGenotype, Site, VariantCall};
use navigator_analysis::coverage::{self, CallableLociParams, CoverageResult};
use navigator_analysis::gvcf;
use navigator_analysis::heteroplasmy::{self, HeteroplasmyParams};
use navigator_analysis::scan::SampleSidecars;
use navigator_analysis::sidecar;
use navigator_analysis::ibd::{
    ChromosomeGenotypes, GeneticMap, IbdSegment, MatchSummary, PairwiseIbdDetector,
};
use navigator_domain::workspace::{Panel, PanelSite};
use navigator_store::panel;

// Re-export the analysis result types the command API returns, so the UI depends only
// on navigator-app (ui -> app), not directly on navigator-analysis.
pub use navigator_analysis::probe::AlignmentProbe;
pub use navigator_analysis::caller::SiteGenotype as PanelGenotype;
pub use navigator_analysis::caller::VariantCall as DenovoCall;
pub use navigator_analysis::coverage::CoverageResult as Coverage;
pub use navigator_analysis::read_metrics::{PairOrientation, ReadMetrics};
pub use navigator_analysis::unified::UnifiedMetricsResult;
pub use navigator_analysis::sex::{Confidence as SexConfidence, InferredSex, SexInferenceResult};
pub use navigator_analysis::sv::types::{SvAnalysisResult, SvCall, SvType};
pub use navigator_analysis::heteroplasmy::HeteroplasmySite;
pub use navigator_analysis::mtvariants::{MtRegion, MtVariant, MtVariantKind};
pub use navigator_analysis::haplo::{BranchEvidence, CallState, ScoredHaplogroup, SnpEvidence};
pub use navigator_analysis::mask::YRegionClass;
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
        self.variants.iter().filter(|v| matches!(v.class, PrivateClass::OffPathKnown(_))).count()
    }
    /// Calls that fall in a curated chrY structural (paralog-prone) region — suspect, to be
    /// down-weighted in reports rather than treated as confident new variants.
    pub fn in_structural_region(&self) -> usize {
        self.variants.iter().filter(|v| v.region.is_some()).count()
    }
    /// Novel calls in *unique* sequence (no structural-region flag) — the high-confidence
    /// new-branch candidates, separated from the paralog-zone noise.
    pub fn novel_in_unique_sequence(&self) -> usize {
        self.variants.iter().filter(|v| v.class == PrivateClass::Novel && v.region.is_none()).count()
    }
}
pub use navigator_analysis::ibd::{
    IbdDetectorConfig, IbdSegment as Segment, MatchSummary as IbdSummary, RelationshipEstimate,
};
// Sync/publish types the command API uses, re-exported so the UI depends only on navigator-app.
pub use navigator_sync::{
    AlignmentRecord, BiosampleRecord, PdsClient, PopulationBreakdownRecord, PrivateVariantsRecord,
    RecordRef, SequenceRunRecord, VariantCallEntry, NS_ALIGNMENT, NS_BIOSAMPLE,
    NS_POPULATION_BREAKDOWN, NS_SEQUENCERUN, PRIVATE_VARIANTS_COLLECTION,
};
use navigator_sync::{FedPopulationComponent, FedSuperPopulationSummary};
use navigator_sync::{dev_http_client, login_default, AsyncSync, DeviceKey, OAuthConfig, RetryPolicy, TokenStore, DEVICE_KEY_COLLECTION};
use navigator_sync::exchange::{self, ExchangeKey};
use navigator_refgenome::{cache as refgenome_cache, canonical_build, Build as ReferenceBuild, LiftedPos, ReferenceGateway};
pub use navigator_refgenome::RefStatus;
pub use navigator_refgenome::{ChromosomeRegions, Cytoband, GenomeRegions, RegionAnnotation};
use navigator_sync::{
    AuditEntryRecord, HaplogroupReconciliationRecord, HeteroplasmyObservationRecord,
    IdentityVerificationRecord, ManualOverrideRecord, ReconciliationStatusRecord,
    RunHaplogroupCallRecord, HAPLOGROUP_RECONCILIATION_COLLECTION,
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
            Some(IbdSuggestion { suggested_sample_guid, suggestion_type, score, signals })
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
use navigator_domain::yprofile::{self, YObsInput};
pub use navigator_domain::yprofile::{YProfileSummary, YProfileVariant, YSourceObs, YState, YVariantStatus};
pub use navigator_domain::consensus::{DiploidSourceObs, DiploidVariant};
use navigator_domain::bisdna;
use navigator_domain::ysnp_dict::{self, YsnpDictionary};
use navigator_store::{
    alignment, ancestry_result, artifact, biosample, chip_profile, haplogroup_call,
    mtdna as mtdna_store, project, reconciliation as recon_store, sequence_run, str_profile,
    sync_history, sync_outbox, variant_set, Store, StoreError,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use uuid::Uuid;

pub mod error;
pub use error::AppError;
pub mod export;
pub mod settings;
pub use settings::AppSettings;

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
    let branches = ranked
        .first()
        .map(|t| haplo::child_evidence(tree, calls, t.id))
        .unwrap_or_default();
    HaploAssignment { ranked, branches }
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
    let branches = ranked
        .first()
        .map(|t| haplo::child_evidence(tree, calls, t.id))
        .unwrap_or_default();
    HaploAssignment { ranked, branches }
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
                allowed.entry(l.position).or_insert((a.to_ascii_uppercase(), d.to_ascii_uppercase()));
            }
        }
    }
    calls
        .into_iter()
        .map(|(pos, base)| match allowed.get(&pos) {
            Some(&(a, d)) if base != a && base != d => {
                let c = complement_base(base);
                if c == a || c == d { (pos, c) } else { (pos, base) }
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
        let base = called
            .variant_bases
            .get(&lp.pos)
            .copied()
            .or_else(|| called.callable.contains(&lp.pos).then(|| ref_base.get(&lp.pos).copied()).flatten());
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

/// Where the ancestry panel for `build` lives: `$NAVIGATOR_ANCESTRY_PANEL` (override), else
/// `<refgenome base>/ancestry/ancestry_panel_<build>.bin`. The offline `navigator-panelbuild`
/// tool writes it; install/ship copies it into the cache dir.
fn ancestry_panel_path(build: ReferenceBuild) -> PathBuf {
    if let Ok(p) = std::env::var("NAVIGATOR_ANCESTRY_PANEL") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("ancestry_panel_{}.bin", build.as_str()))
}

/// Where the PCA loadings for `build` live: `$NAVIGATOR_ANCESTRY_PCA` (override), else
/// `<refgenome base>/ancestry/ancestry_pca_<build>.bin`. Optional — absent means the
/// AF-likelihood estimate runs without PCA coordinates.
fn ancestry_pca_path(build: ReferenceBuild) -> PathBuf {
    if let Ok(p) = std::env::var("NAVIGATOR_ANCESTRY_PCA") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("ancestry_pca_{}.bin", build.as_str()))
}

/// Where the **ancient** PCA loadings for `build` live: `$NAVIGATOR_ANCESTRY_PCA_ANCIENT`
/// (override), else `<refgenome base>/ancestry/ancestry_pca_ancient_<build>.bin`. Optional —
/// present means the PCA-projection GMM runs against ancient reference components
/// (Steppe/EEF/WHG) instead of the modern super-populations. Must be built over the same panel
/// sites the AF panel genotypes (so the single genotyping pass covers it).
fn ancestry_pca_ancient_path(build: ReferenceBuild) -> PathBuf {
    if let Ok(p) = std::env::var("NAVIGATOR_ANCESTRY_PCA_ANCIENT") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("ancestry_pca_ancient_{}.bin", build.as_str()))
}

/// The fine-population frequency asset path (`$NAVIGATOR_ANCESTRY_FREQ` override, else
/// `<base>/ancestry/ancestry_freq_global_<build>.bin`). Optional — fine admixture is skipped if absent.
fn ancestry_freq_global_path(build: ReferenceBuild) -> PathBuf {
    if let Ok(p) = std::env::var("NAVIGATOR_ANCESTRY_FREQ") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("ancestry_freq_global_{}.bin", build.as_str()))
}

/// The chip-compatible IBD panel asset path (`$NAVIGATOR_IBD_PANEL` override, else
/// `<base>/ancestry/ibd_panel_<build>.bin`).
fn ibd_panel_path(build: ReferenceBuild) -> PathBuf {
    if let Ok(p) = std::env::var("NAVIGATOR_IBD_PANEL") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("ibd_panel_{}.bin", build.as_str()))
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
        AssetStatus { name: name.to_string(), present: bytes.is_some(), verified }
    })
    .collect()
}

/// The asset integrity manifest path for a build (`<base>/ancestry/ancestry_manifest_<build>.json`).
fn ancestry_manifest_path(build: ReferenceBuild) -> PathBuf {
    refgenome_cache::base_dir().join("ancestry").join(format!("ancestry_manifest_{}.json", build.as_str()))
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
    let Ok(bytes) = std::fs::read(path) else { return Ok(None) };
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
    if let Ok(p) = std::env::var("NAVIGATOR_GENETIC_MAP") {
        return PathBuf::from(p);
    }
    refgenome_cache::base_dir()
        .join("ancestry")
        .join(format!("genetic_map_{}.bin", build.as_str()))
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
}

impl ExportRequest {
    /// File extension (no dot) for the save dialog + filter.
    pub fn extension(&self) -> &'static str {
        match self {
            ExportRequest::CoverageHtml(_) | ExportRequest::AncestryHtml(_) => "html",
            ExportRequest::CallableBed(_) => "bed",
            ExportRequest::DiploidVcf(_) => "vcf",
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
        };
        format!("{stem}_{id}.{}", self.extension())
    }
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

/// Cheap VCF header peek: the `##` meta block (joined) + the contig names from `##contig=<ID=…>`.
/// Reads only the header (stops at the first data line). Plain text — matches the import parser,
/// which doesn't decompress either; a gzipped VCF simply yields an empty peek (→ generic).
fn peek_vcf_header(path: &Path) -> (String, Vec<String>) {
    use std::io::BufRead;
    let Ok(file) = std::fs::File::open(path) else { return (String::new(), Vec::new()) };
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
    let generic = matches!(filename.to_ascii_lowercase().as_str(), "variants.vcf" | "variants.vcf.gz");
    if generic {
        if let Some(parent) = path.parent().and_then(|p| p.file_name()).and_then(|n| n.to_str()) {
            if !parent.is_empty() {
                return parent.to_string();
            }
        }
    }
    filename.to_string()
}

/// Stream a file through SHA-256 and return the lowercase hex digest. Blocking (reads the
/// whole file in 1 MiB chunks) — call via [`sha256_file_async`] for large alignments.
fn sha256_file(path: &Path) -> std::io::Result<String> {
    use sha2::{Digest, Sha256};
    use std::io::Read;
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let digest = hasher.finalize();
    let mut hex = String::with_capacity(digest.len() * 2);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(hex, "{b:02x}");
    }
    Ok(hex)
}

/// SHA-256 of a file's content (hex), computed off the async runtime.
async fn sha256_file_async(path: PathBuf) -> Result<String, AppError> {
    let hash = tokio::task::spawn_blocking(move || sha256_file(&path))
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
    Ok(hash)
}

/// SHA-256 of an in-memory string (hex) — for hashing tree JSON / small content.
fn sha256_str(s: &str) -> String {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(s.as_bytes());
    let mut hex = String::with_capacity(digest.len() * 2);
    for b in digest {
        use std::fmt::Write;
        let _ = write!(hex, "{b:02x}");
    }
    hex
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
    resolve_appview_url(std::env::var("DECODINGUS_APPVIEW_URL").ok(), AppSettings::load().appview_url)
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
fn snp_obs_from_assignment(assignment: &HaploAssignment, in_tree: bool) -> Vec<YObsInput> {
    let mut by_name: std::collections::HashMap<String, YObsInput> = std::collections::HashMap::new();
    for branch in &assignment.branches {
        for snp in &branch.snps {
            let state = match snp.state {
                CallState::Derived => YState::Derived,
                CallState::Ancestral => YState::Ancestral,
                CallState::NoCall => YState::NoCall,
            };
            by_name.entry(snp.name.clone()).or_insert_with(|| {
                YObsInput::snp(snp.name.clone(), snp.position, snp.ancestral.clone(), snp.derived.clone(), state, in_tree)
            });
        }
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
/// Artifact kind for an alignment's cached IBD-panel genotypes (distinct from the store-panel
/// genotype cache, which is keyed by `panel_kind(panel_id, ploidy)`).
const IBD_PANEL_KIND: &str = "ibd_panel_genotypes";

/// Count of sites called (dosage within ploidy) in **both** samples — the effective IBD comparison
/// size, surfaced so a sparse chip↔chip / chip↔WGS overlap isn't mistaken for a confident result.
fn overlapping_called_sites(a: &[SiteGenotype], b: &[SiteGenotype]) -> usize {
    let called = |g: &SiteGenotype| (0..=g.ploidy as i32).contains(&g.dosage);
    let set: std::collections::HashSet<(&str, i64)> =
        a.iter().filter(|g| called(g)).map(|g| (g.contig.as_str(), g.position)).collect();
    b.iter()
        .filter(|g| called(g))
        .filter(|g| set.contains(&(g.contig.as_str(), g.position)))
        .count()
}

/// Group two samples' dosages, load the genetic map for `build`, detect IBD segments, and record
/// the overlapping-site count. Shared by the alignment-pair and chip-or-WGS compare paths.
fn detect_ibd(ga: &[SiteGenotype], gb: &[SiteGenotype], build: ReferenceBuild, config: IbdDetectorConfig) -> IbdComparison {
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
    IbdComparison { summary, segments, overlapping_sites }
}

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

    /// Update a project's editable fields (name required; description optional; administrator
    /// defaults to "unknown" when blank). Returns the updated record.
    pub async fn update_project(
        &self,
        id: i64,
        name: String,
        description: Option<String>,
        administrator: String,
    ) -> Result<Project, AppError> {
        let name = name.trim();
        if name.is_empty() {
            return Err(AppError::Conflict("project name cannot be empty".into()));
        }
        let desc = description.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let admin = administrator.trim();
        let admin = if admin.is_empty() { "unknown" } else { admin };
        let updated = project::update(self.store.pool(), id, name, desc.as_deref(), admin).await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("project {id}"))));
        }
        project::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("project {id}"))))
    }

    /// Delete a project. Refused (with a clear message) while subjects still belong to it, so
    /// the user reassigns them first rather than orphaning the rows.
    pub async fn delete_project(&self, id: i64) -> Result<(), AppError> {
        let members = biosample::count_for_project(self.store.pool(), id).await?;
        if members > 0 {
            return Err(AppError::Conflict(format!(
                "cannot delete project: {members} subject(s) still belong to it — reassign them first"
            )));
        }
        if !project::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("project {id}"))));
        }
        Ok(())
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

    /// Update a subject's editable fields (identity, accession, description, center, sex).
    /// Empty strings are normalized to NULL. Returns the updated record.
    pub async fn update_biosample(
        &self,
        guid: SampleGuid,
        donor_identifier: String,
        sample_accession: Option<String>,
        description: Option<String>,
        center_name: Option<String>,
        sex: Option<String>,
    ) -> Result<Biosample, AppError> {
        let donor = donor_identifier.trim();
        if donor.is_empty() {
            return Err(AppError::Conflict("subject identifier cannot be empty".into()));
        }
        let norm = |o: Option<String>| o.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let (acc, desc, center, sex) = (norm(sample_accession), norm(description), norm(center_name), norm(sex));
        let updated = biosample::update(
            self.store.pool(),
            guid,
            donor,
            acc.as_deref(),
            desc.as_deref(),
            center.as_deref(),
            sex.as_deref(),
        )
        .await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        biosample::get(self.store.pool(), guid)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))))
    }

    /// Assign a subject to a project (validating the project exists). `None` clears it.
    pub async fn add_biosample_to_project(&self, guid: SampleGuid, project_id: Option<i64>) -> Result<(), AppError> {
        if let Some(pid) = project_id {
            if project::get(self.store.pool(), pid).await?.is_none() {
                return Err(AppError::Store(StoreError::NotFound(format!("project {pid}"))));
            }
        }
        if !biosample::set_project(self.store.pool(), guid, project_id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        Ok(())
    }

    /// Delete a subject. Refused (with a clear message) when it still has dependent data —
    /// sequencing runs or any imported profile — so the user removes data first rather than
    /// silently orphaning rows.
    pub async fn delete_biosample(&self, guid: SampleGuid) -> Result<(), AppError> {
        let runs = self.list_sequence_runs(guid).await?.len();
        let strs = self.list_str_profiles(guid).await?.len();
        let variants = self.list_variant_sets(guid).await?.len();
        let chips = self.list_chip_profiles(guid).await?.len();
        let mt = self.list_mtdna_sequences(guid).await?.len();
        let total = runs + strs + variants + chips + mt;
        if total > 0 {
            return Err(AppError::Conflict(format!(
                "cannot delete subject: it still has {runs} sequencing run(s), {strs} STR, \
                 {variants} variant-set, {chips} chip, {mt} mtDNA record(s) — remove its data first"
            )));
        }
        if !biosample::delete(self.store.pool(), guid).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("biosample {}", guid.0))));
        }
        Ok(())
    }

    pub async fn record_sequence_run(&self, run: NewSequenceRun) -> Result<SequenceRun, AppError> {
        Ok(sequence_run::create(self.store.pool(), &run).await?)
    }

    pub async fn record_alignment(&self, aln: NewAlignment) -> Result<Alignment, AppError> {
        Ok(alignment::create(self.store.pool(), &aln).await?)
    }

    /// Update a sequence run's descriptive fields (test type required; platform defaults to
    /// "UNKNOWN" when blank; instrument/layout optional). Read metrics are preserved. Returns
    /// the updated record.
    pub async fn update_sequence_run(
        &self,
        id: i64,
        platform_name: String,
        instrument_model: Option<String>,
        test_type: String,
        library_layout: Option<String>,
        sequencing_facility: Option<String>,
    ) -> Result<SequenceRun, AppError> {
        let test_type = test_type.trim();
        if test_type.is_empty() {
            return Err(AppError::Conflict("test type cannot be empty".into()));
        }
        let norm = |o: Option<String>| o.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let platform = platform_name.trim();
        let platform = if platform.is_empty() { "UNKNOWN" } else { platform };
        let updated = sequence_run::update(
            self.store.pool(),
            id,
            platform,
            norm(instrument_model).as_deref(),
            test_type,
            norm(library_layout).as_deref(),
            norm(sequencing_facility).as_deref(),
        )
        .await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("sequence run {id}"))));
        }
        sequence_run::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("sequence run {id}"))))
    }

    /// Update an alignment's descriptive fields (reference build + aligner required; variant
    /// caller optional). File paths are managed by import/probe. Returns the updated record.
    pub async fn update_alignment(
        &self,
        id: i64,
        reference_build: String,
        aligner: String,
        variant_caller: Option<String>,
    ) -> Result<Alignment, AppError> {
        let build = reference_build.trim();
        let aligner = aligner.trim();
        if build.is_empty() || aligner.is_empty() {
            return Err(AppError::Conflict("reference build and aligner are required".into()));
        }
        let caller = variant_caller.map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        let updated = alignment::update(self.store.pool(), id, build, aligner, caller.as_deref()).await?;
        if !updated {
            return Err(AppError::Store(StoreError::NotFound(format!("alignment {id}"))));
        }
        alignment::get(self.store.pool(), id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {id}"))))
    }

    /// Delete a sequence run and everything beneath it (its alignments + cached analysis
    /// artifacts). This is how a mistaken BAM/CRAM import is undone.
    pub async fn delete_sequence_run(&self, id: i64) -> Result<(), AppError> {
        if !sequence_run::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("sequence run {id}"))));
        }
        Ok(())
    }

    /// Delete a single alignment and its cached analysis artifacts (the parent run is kept).
    pub async fn delete_alignment(&self, id: i64) -> Result<(), AppError> {
        if !alignment::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("alignment {id}"))));
        }
        Ok(())
    }

    /// Delete an imported STR profile (and its markers).
    pub async fn delete_str_profile(&self, id: i64) -> Result<(), AppError> {
        if !str_profile::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("STR profile {id}"))));
        }
        Ok(())
    }

    /// Delete an imported variant set (and its calls).
    pub async fn delete_variant_set(&self, id: i64) -> Result<(), AppError> {
        if !variant_set::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("variant set {id}"))));
        }
        Ok(())
    }

    /// Delete an imported chip/array profile.
    pub async fn delete_chip_profile(&self, id: i64) -> Result<(), AppError> {
        if !chip_profile::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("chip profile {id}"))));
        }
        Ok(())
    }

    /// Delete an imported mtDNA sequence.
    pub async fn delete_mtdna_sequence(&self, id: i64) -> Result<(), AppError> {
        if !mtdna_store::delete(self.store.pool(), id).await? {
            return Err(AppError::Store(StoreError::NotFound(format!("mtDNA sequence {id}"))));
        }
        Ok(())
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
        // Default provenance: a full result from a Navigator CRAM walk.
        self.save_analysis_with_provenance(alignment_id, kind, algorithm_version, result, "navigator-walk", "full")
            .await
    }

    /// Like [`save_analysis`] but stamps provenance: `source` (`navigator-walk` |
    /// `pipeline-sidecar`) and `completeness` (`full` | `partial`). The fast-path sidecar
    /// ingest uses this so the manual deep pass can tell a sidecar/partial result apart from a
    /// full walk and upgrade it rather than skip it.
    pub async fn save_analysis_with_provenance<T: Serialize>(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
        result: &T,
        source: &str,
        completeness: &str,
    ) -> Result<AnalysisArtifact, AppError> {
        let payload = serde_json::to_string(result)?;
        Ok(artifact::upsert(self.store.pool(), alignment_id, kind, algorithm_version, Utc::now(), &payload, source, completeness).await?)
    }

    /// `(source, completeness)` of a cached artifact, defaulting `None` columns to
    /// `("navigator-walk", "full")` (pre-provenance rows). `None` when no artifact exists.
    pub async fn analysis_provenance(
        &self,
        alignment_id: i64,
        kind: &str,
        algorithm_version: &str,
    ) -> Result<Option<(String, String)>, AppError> {
        Ok(artifact::get(self.store.pool(), alignment_id, kind, algorithm_version).await?.map(|a| {
            (
                a.source.unwrap_or_else(|| "navigator-walk".into()),
                a.completeness.unwrap_or_else(|| "full".into()),
            )
        }))
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
        self.run_coverage_for_alignment_with_progress(alignment_id, |_, _| {}).await
    }

    /// Like [`run_coverage_for_alignment`], reporting `progress(contigs_done, contigs_total)` as
    /// the whole-genome pass walks each contig (the slow step — minutes on a real WGS BAM — so a
    /// progress bar can advance instead of sitting frozen). The callback runs on the blocking
    /// thread.
    pub async fn run_coverage_for_alignment_with_progress(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(usize, usize) + Send + 'static,
    ) -> Result<CoverageResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        // The reference isn't asked for at import — resolve the alignment's build via the gateway
        // (cached, else download) when no FASTA was stored.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        let mut params = CallableLociParams::default();
        let result = tokio::task::spawn_blocking(move || {
            // Adapt the callable threshold to read tech (HiFi → 1×; see adaptive_min_depth).
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, Some(&reference)) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            coverage::collect_coverage_callable_with_progress(&bam, &reference, &params, None, &mut progress)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result).await?;
        Ok(result)
    }

    /// Infer biological sex from the alignment's chrX:autosome read-density ratio, persisting
    /// the result as a `sex` artifact. Cheap (BAI fast-path for BAM). `reference` is used only
    /// for CRAM decode.
    pub async fn run_sex(&self, alignment_id: i64) -> Result<navigator_analysis::sex::SexInferenceResult, AppError> {
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let result = tokio::task::spawn_blocking(move || {
            navigator_analysis::sex::infer_from_bam(&bam, reference.as_deref())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "sex", "1", &result).await?;
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(result)
    }

    /// Write the inferred sex back to the biosample when the user didn't provide one, so it
    /// shows in the subjects table + header instead of "Unknown". No-op for Unknown sex or
    /// when the biosample already carries a sex.
    async fn write_back_inferred_sex(
        &self,
        alignment_id: i64,
        result: &navigator_analysis::sex::SexInferenceResult,
    ) -> Result<(), AppError> {
        let label = match result.inferred_sex {
            InferredSex::Male => Some("Male"),
            InferredSex::Female => Some("Female"),
            InferredSex::Unknown => None,
        };
        if let (Some(label), Ok(guid)) = (label, self.biosample_of_alignment(alignment_id).await) {
            if let Ok(Some(bio)) = biosample::get(self.store.pool(), guid).await {
                if bio.sex.as_deref().map(str::trim).unwrap_or("").is_empty() {
                    biosample::set_sex(self.store.pool(), guid, label).await?;
                }
            }
        }
        Ok(())
    }

    /// Cached `sex` inference, if present.
    pub async fn cached_sex(&self, alignment_id: i64) -> Result<Option<navigator_analysis::sex::SexInferenceResult>, AppError> {
        self.load_analysis(alignment_id, "sex", "1").await
    }

    /// Collect read-level QC metrics (alignment summary + read-length/insert-size distributions,
    /// pair orientation, mean MAPQ) and persist as a `read_metrics` artifact.
    pub async fn run_read_metrics(&self, alignment_id: i64) -> Result<navigator_analysis::read_metrics::ReadMetrics, AppError> {
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let result = tokio::task::spawn_blocking(move || {
            navigator_analysis::read_metrics::collect_read_metrics(&bam, reference.as_deref())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "read_metrics", "1", &result).await?;
        self.write_back_read_stats(alignment_id, &result).await?;
        Ok(result)
    }

    /// Mirror an alignment's library-level read stats onto its owning sequence run (`total_reads`,
    /// `mean_read_length`, `mean_insert_size`) so the Data Sources run card shows them without
    /// re-walking. Best-effort: a missing alignment/run is ignored. When a run has several
    /// alignments the last write wins — these are per-library properties, so any pass is
    /// representative.
    async fn write_back_read_stats(
        &self,
        alignment_id: i64,
        m: &navigator_analysis::read_metrics::ReadMetrics,
    ) -> Result<(), AppError> {
        if let Some(aln) = alignment::get(self.store.pool(), alignment_id).await? {
            // Paired-end evidence: any reads aligned in pairs ⇒ PAIRED. Only overrides the stored
            // layout when we have aligned reads to judge (else leave the import-time flag value).
            let layout = (m.pf_reads_aligned > 0)
                .then_some(if m.reads_aligned_in_pairs > 0 { "PAIRED" } else { "SINGLE" });
            sequence_run::set_read_stats(
                self.store.pool(),
                aln.sequence_run_id,
                Some(m.total_reads as i64),
                (m.mean_read_length > 0.0).then_some(m.mean_read_length),
                (m.mean_insert_size > 0.0).then_some(m.mean_insert_size),
                layout,
            )
            .await?;
        }
        Ok(())
    }

    /// Cached `read_metrics`, if present.
    pub async fn cached_read_metrics(&self, alignment_id: i64) -> Result<Option<navigator_analysis::read_metrics::ReadMetrics>, AppError> {
        self.load_analysis(alignment_id, "read_metrics", "1").await
    }

    /// Run the unified quality-metrics walker — coverage + callable, read-level QC metrics, and
    /// sex inference in **one pass** over the alignment's BAM/CRAM (vs. the separate passes
    /// `run_coverage` + `run_read_metrics` + `run_sex` cost: 2 reads for BAM, 3 for CRAM). All
    /// three sub-results are persisted under their existing artifact keys (`coverage`/
    /// `COVERAGE_VERSION`, `read_metrics`/`"1"`, `sex`/`"1"`), so `cached_coverage`/
    /// `cached_read_metrics`/`cached_sex` and the SV step's reuse logic keep working unchanged.
    pub async fn run_unified_metrics(&self, alignment_id: i64) -> Result<UnifiedMetricsResult, AppError> {
        self.run_unified_metrics_with_progress(alignment_id, |_, _| {}).await
    }

    /// Like [`run_unified_metrics`], reporting `progress(contigs_done, contigs_total)` as the
    /// (slow) whole-genome coverage portion finalizes each contig. Uses the per-contig parallel
    /// walker (falling back to a sequential pass for CRAM / unindexed BAM); the callback is
    /// `Fn + Sync` because it's invoked concurrently from the fan-out's worker threads.
    pub async fn run_unified_metrics_with_progress(
        &self,
        alignment_id: i64,
        progress: impl Fn(usize, usize) + Send + Sync + 'static,
    ) -> Result<UnifiedMetricsResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        // The walker requires a reference (CRAM decode + reference-N detection); resolve the
        // build via the gateway when no FASTA was stored at import.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        let mut params = CallableLociParams::default();
        let result = tokio::task::spawn_blocking(move || {
            // Adapt the callable threshold to read tech (HiFi → 1×; see adaptive_min_depth).
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, Some(&reference)) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            navigator_analysis::unified::collect_unified_metrics_parallel_with_progress(
                &bam,
                &reference,
                &params,
                None,
                &progress,
            )
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;

        // Persist each sub-result under its own existing cache key.
        self.save_analysis(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result.coverage).await?;
        self.save_analysis(alignment_id, "read_metrics", "1", &result.read_metrics).await?;
        self.write_back_read_stats(alignment_id, &result.read_metrics).await?;
        if let Some(sex) = &result.sex {
            self.save_analysis(alignment_id, "sex", "1", sex).await?;
            self.write_back_inferred_sex(alignment_id, sex).await?;
        }
        Ok(result)
    }

    /// Call structural variants (depth-segmentation + paired-end/split-read evidence) and
    /// persist as an `sv` artifact. Needs coverage + insert-size inputs (computed/loaded here)
    /// and **≥10× mean coverage** (the caller errors below that).
    pub async fn run_sv(&self, alignment_id: i64) -> Result<navigator_analysis::sv::types::SvAnalysisResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = aln.reference_path.clone().map(PathBuf::from);
        let reference_build = aln.reference_build.clone();

        let cov = match self.cached_coverage(alignment_id).await? {
            Some(c) => c,
            None => self.run_coverage_for_alignment(alignment_id).await?,
        };
        let rm = match self.cached_read_metrics(alignment_id).await? {
            Some(m) => m,
            None => self.run_read_metrics(alignment_id).await?,
        };
        let (mean_cov, mean_ins, sd_ins, mean_rl) =
            (cov.mean_coverage, rm.mean_insert_size, rm.std_insert_size, rm.mean_read_length);

        let result = tokio::task::spawn_blocking(move || {
            let lengths = caller::header_contig_lengths(&bam, reference.as_deref())?;
            navigator_analysis::sv::caller::call_structural_variants(
                &bam,
                &lengths,
                &reference_build,
                mean_cov,
                mean_ins,
                sd_ins,
                mean_rl,
                &navigator_analysis::sv::types::SvCallerConfig::default(),
            )
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, "sv", "1", &result).await?;
        Ok(result)
    }

    /// Cached `sv` result, if present.
    pub async fn cached_sv(&self, alignment_id: i64) -> Result<Option<navigator_analysis::sv::types::SvAnalysisResult>, AppError> {
        self.load_analysis(alignment_id, "sv", "1").await
    }

    /// The alignment's BAM (required) + reference (optional; required only for CRAM).
    async fn alignment_paths(&self, alignment_id: i64) -> Result<(PathBuf, Option<PathBuf>), AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        Ok((bam, aln.reference_path.map(PathBuf::from)))
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

    /// Whole-contig **de-novo diploid** SNV calling (het 0/1 + hom-alt 1/1) on `contig`, cached per
    /// alignment+contig. Reuses the alignment's BAM + reference (resolved from the build). Returns
    /// [`SiteGenotype`]s in position order — feed to [`Self::diploid_vcf`].
    pub async fn run_diploid_calls(&self, alignment_id: i64, contig: String) -> Result<Vec<SiteGenotype>, AppError> {
        let kind = format!("diploid_denovo:{contig}");
        if let Some(c) = self.load_analysis(alignment_id, &kind, caller::GENOTYPE_VERSION).await? {
            return Ok(c);
        }
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let params = adaptive_haploid_params(&bam, Some(&reference));
        let calls = tokio::task::spawn_blocking(move || caller::call_denovo_diploid(&bam, &reference, &contig, &params))
            .await
            .map_err(|e| AppError::Join(e.to_string()))??;
        self.save_analysis(alignment_id, &kind, caller::GENOTYPE_VERSION, &calls).await?;
        Ok(calls)
    }

    /// A diploid VCF (VCFv4.2, `GT:AD:DP:GQ:PL`) of the de-novo diploid SNV calls for `contig`
    /// (computing + caching them if needed). The sample column is `aln<id>`.
    pub async fn diploid_vcf(&self, alignment_id: i64, contig: String) -> Result<String, AppError> {
        let calls = self.run_diploid_calls(alignment_id, contig).await?;
        Ok(navigator_analysis::vcf::write_diploid_vcf(&format!("aln{alignment_id}"), &calls))
    }

    /// A **whole-genome** diploid VCF: de-novo SNV + indel calls over every primary chromosome
    /// (1–22, X, Y, M) of the alignment, per-contig cached. Heavy (a real WGS calling pass); the
    /// caller runs it off the UI thread (the export path).
    pub async fn diploid_vcf_genome(&self, alignment_id: i64) -> Result<String, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let contigs = tokio::task::spawn_blocking(move || caller::header_contig_names(&bam, Some(&reference)))
            .await
            .map_err(|e| AppError::Join(e.to_string()))??;
        let mut all = Vec::new();
        for contig in contigs.into_iter().filter(|c| is_primary_contig(c)) {
            all.extend(self.run_diploid_calls(alignment_id, contig).await?);
        }
        Ok(navigator_analysis::vcf::write_diploid_vcf(&format!("aln{alignment_id}"), &all))
    }

    /// Run de-novo calling on `contig` using the alignment's own stored paths.
    /// The alignment's BAM + a usable reference FASTA: the stored path, else resolved from the
    /// alignment's build via the gateway (cached, else downloaded). Errors only if no BAM is
    /// recorded. Use this in steps that *require* the reference, so the user never has to supply
    /// one (it follows from the header-detected build).
    async fn alignment_bam_reference(&self, alignment_id: i64) -> Result<(PathBuf, PathBuf), AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        Ok((bam, reference))
    }

    pub async fn run_denovo_for_alignment(&self, alignment_id: i64, contig: String) -> Result<Vec<VariantCall>, AppError> {
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let probe = bam.clone();
        let probe_ref = reference.clone();
        let params = tokio::task::spawn_blocking(move || adaptive_haploid_params(&probe, Some(&probe_ref)))
            .await
            .map_err(|e| AppError::Join(e.to_string()))?; // HiFi -> lower min_depth
        self.run_denovo_caller(alignment_id, bam, reference, contig, params).await
    }

    // ---- publish -----------------------------------------------------------

    /// Build the alignment (coverage) record JSON for an alignment — the shared
    /// `com.decodingus.atmosphere.alignment` contract the AppView ingests (floats as strings).
    async fn coverage_record(&self, alignment_id: i64) -> Result<serde_json::Value, AppError> {
        let cov = self
            .cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))?;
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let record = AlignmentRecord::new(
            aln.reference_build,
            Some(aln.aligner),
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

    /// Build the population-breakdown (ancestry) record JSON for an alignment from its
    /// persisted estimate — the shared `com.decodingus.atmosphere.populationBreakdown`
    /// contract the AppView ingests (floats as strings).
    /// The populationBreakdown record JSON for each persisted estimate of an alignment (one
    /// per method — e.g. ADMIXTURE + PCA_PROJECTION_GMM). Empty if none computed.
    async fn ancestry_records(&self, alignment_id: i64) -> Result<Vec<serde_json::Value>, AppError> {
        let results = ancestry_result::list_for_alignment(self.store.pool(), alignment_id).await?;
        results
            .iter()
            .map(|r| serde_json::to_value(population_breakdown_record(r)).map_err(AppError::from))
            .collect()
    }

    /// Build the anonymized biosample record JSON — sex, center, and best-effort Y/mt
    /// haplogroup calls. Donor identifiers / accession / description are never carried.
    async fn biosample_record(&self, biosample_guid: SampleGuid) -> Result<serde_json::Value, AppError> {
        let bio = biosample::get(self.store.pool(), biosample_guid)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("biosample {biosample_guid:?}"))))?;
        let y = self.consensus_haplogroup(biosample_guid, DnaType::Y).await?;
        let mt = self.consensus_haplogroup(biosample_guid, DnaType::Mt).await?;
        let runs = self.list_sequence_runs(biosample_guid).await?;
        let record = BiosampleRecord::new(bio.sex, y, mt, bio.center_name, Utc::now().to_rfc3339())
            .with_refs(runs.iter().map(|r| r.id.to_string()).collect(), None, None);
        Ok(serde_json::to_value(&record)?)
    }

    /// Build a sequence-run characterization record JSON (platform/instrument/test — no files).
    /// `instrument_id` (the sequencer serial inferred from read names) is published so the AppView
    /// can grow its crowd-sourced instrument→lab map (`fed.sequencerun.instrument_id` → the
    /// `instrument_observation`→proposal→accept consensus). It identifies the physical sequencer,
    /// not the donor — no PII, consistent with the anonymized fed-record posture.
    async fn sequence_run_record(&self, run: &SequenceRun) -> Result<serde_json::Value, AppError> {
        let record = SequenceRunRecord::new(
            None,
            Some(run.platform_name.clone()),
            run.instrument_model.clone(),
            run.instrument_id.clone(),
            Some(run.test_type.clone()),
            run.library_layout.clone(),
            run.total_reads,
            run.mean_read_length.map(|l| l.round() as i32),
            run.mean_insert_size,
            Utc::now().to_rfc3339(),
        );
        Ok(serde_json::to_value(&record)?)
    }

    /// Best-effort consensus haplogroup for a subject arm: the manual override if set,
    /// else the first recorded per-source call. `None` when nothing has been called.
    async fn consensus_haplogroup(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<String>, AppError> {
        if let Some((hg, _)) = recon_store::get_override(self.store.pool(), biosample_guid, dna_type).await? {
            return Ok(Some(hg));
        }
        let calls = haplogroup_call::list_for(self.store.pool(), biosample_guid, dna_type).await?;
        Ok(reconciliation::reconcile(&calls).map(|c| c.haplogroup))
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
        Ok(client.create_record(NS_ALIGNMENT, value, None).await?)
    }

    /// Publish every persisted ancestry estimate for an alignment (one populationBreakdown per
    /// method — admixture + PCA-GMM) using an explicit `client` (the testable core; production
    /// callers use [`publish_ancestry`](Self::publish_ancestry)). Returns a ref per record.
    pub async fn publish_ancestry_with(
        &self,
        client: &PdsClient,
        alignment_id: i64,
    ) -> Result<Vec<RecordRef>, AppError> {
        let mut refs = Vec::new();
        for value in self.ancestry_records(alignment_id).await? {
            refs.push(client.create_record(NS_POPULATION_BREAKDOWN, value, None).await?);
        }
        Ok(refs)
    }

    /// Publish the anonymized biosample summary using an explicit `client`.
    pub async fn publish_biosample_with(
        &self,
        client: &PdsClient,
        biosample_guid: SampleGuid,
    ) -> Result<RecordRef, AppError> {
        let value = self.biosample_record(biosample_guid).await?;
        Ok(client.create_record(NS_BIOSAMPLE, value, None).await?)
    }

    /// Publish a sequence-run characterization using an explicit `client`.
    pub async fn publish_sequence_run_with(
        &self,
        client: &PdsClient,
        run: &SequenceRun,
    ) -> Result<RecordRef, AppError> {
        let value = self.sequence_run_record(run).await?;
        Ok(client.create_record(NS_SEQUENCERUN, value, None).await?)
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

    /// The signed-in account's DID, or [`AppError::NotAuthenticated`] — the cheap auth guard publish
    /// methods run before building a record / touching the DB.
    fn require_account(&self) -> Result<String, AppError> {
        self.current_account().ok_or(AppError::NotAuthenticated)
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

    // ---- sync durability: outbox enqueue + drain (gap §5) -------------------

    /// Enqueue a built record for publishing to the signed-in account's PDS. The publish becomes
    /// durable: it survives restart and retries automatically (with backoff) on a transient/offline
    /// failure instead of being lost. Re-enqueuing the same `entity_ref` coalesces (newest wins).
    /// Errors [`AppError::NotAuthenticated`] when signed out (we need the destination DID).
    async fn enqueue_publish(
        &self,
        kind: &str,
        entity_ref: &str,
        collection: &str,
        rkey: Option<&str>,
        value: serde_json::Value,
    ) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let entry = sync_outbox::NewOutboxEntry {
            account_did: did,
            kind: kind.to_string(),
            entity_ref: entity_ref.to_string(),
            collection: collection.to_string(),
            rkey: rkey.map(str::to_string),
            payload: serde_json::to_string(&value)?,
        };
        sync_outbox::enqueue(self.store.pool(), &entry, &Utc::now().to_rfc3339()).await?;
        Ok(())
    }

    /// Pending (not-yet-published) outbox rows for the signed-in account — drives the UI's
    /// "N pending" indicator. `0` when signed out.
    pub async fn outbox_pending_count(&self) -> Result<i64, AppError> {
        let Some(did) = self.current_account() else { return Ok(0) };
        Ok(sync_outbox::pending_count(self.store.pool(), &did).await?)
    }

    /// All non-completed outbox rows (PENDING + FAILED) for the signed-in account — a sync detail view.
    pub async fn outbox_entries(&self) -> Result<Vec<sync_outbox::OutboxEntry>, AppError> {
        let Some(did) = self.current_account() else { return Ok(Vec::new()) };
        Ok(sync_outbox::list(self.store.pool(), &did).await?)
    }

    /// Recent publish outcomes (success/failure) for the signed-in account — the audit trail.
    pub async fn sync_history(&self, limit: i64) -> Result<Vec<sync_history::HistoryEntry>, AppError> {
        let Some(did) = self.current_account() else { return Ok(Vec::new()) };
        Ok(sync_history::recent(self.store.pool(), &did, limit).await?)
    }

    /// Attempt to publish the ready outbox rows for the signed-in account. Each success is logged to
    /// history and its row removed; a transient failure reschedules the row with exponential backoff
    /// and stops the batch (we're likely offline); a non-transient failure marks the row `FAILED`.
    /// A no-op (and `Ok`) when signed out. Safe to call repeatedly (periodically + after a publish).
    pub async fn drain_outbox(&self) -> Result<DrainOutcome, AppError> {
        let Some(did) = self.current_account() else { return Ok(DrainOutcome::default()) };
        let mut outcome = DrainOutcome::default();
        // Build the resilient engine once (loads the session). Signed-out / no session → nothing to do.
        let mut engine = match self.sync_engine() {
            Ok(e) => e,
            Err(_) => return Ok(outcome),
        };
        let now = Utc::now();
        let batch = sync_outbox::ready(self.store.pool(), &did, &now.to_rfc3339(), OUTBOX_BATCH).await?;
        for entry in batch {
            let value: serde_json::Value = serde_json::from_str(&entry.payload)?;
            let result = match &entry.rkey {
                Some(rk) => engine.push_create_rkey(&entry.collection, value, rk).await,
                None => engine.push_create(&entry.collection, value).await,
            };
            let attempt = entry.attempt_count + 1;
            match result {
                Ok(rref) => {
                    self.log_history(&entry, "SUCCESS", Some(&rref), attempt, None).await?;
                    sync_outbox::complete(self.store.pool(), entry.id).await?;
                    outcome.published.push((entry.kind.clone(), rref.uri));
                }
                Err(e) if e.is_transient() => {
                    // Offline / 5xx / timeout: back off and stop — the rest of the batch will wait too.
                    let next = now + chrono::Duration::seconds(backoff_secs(attempt));
                    sync_outbox::reschedule(
                        self.store.pool(),
                        entry.id,
                        attempt,
                        &next.to_rfc3339(),
                        &e.to_string(),
                        &now.to_rfc3339(),
                    )
                    .await?;
                    outcome.retry_scheduled += 1;
                    break;
                }
                Err(e) => {
                    // Validation / auth / other terminal error: give up on this row (visible as FAILED).
                    self.log_history(&entry, "FAILED", None, attempt, Some(&e.to_string())).await?;
                    sync_outbox::mark_failed(self.store.pool(), entry.id, attempt, &e.to_string(), &now.to_rfc3339()).await?;
                    outcome.failed += 1;
                }
            }
        }
        outcome.pending = sync_outbox::pending_count(self.store.pool(), &did).await?;
        Ok(outcome)
    }

    /// Append a sync-history row for a finished push attempt.
    async fn log_history(
        &self,
        entry: &sync_outbox::OutboxEntry,
        status: &str,
        rref: Option<&RecordRef>,
        attempt_count: i64,
        error: Option<&str>,
    ) -> Result<(), AppError> {
        let h = sync_history::NewHistoryEntry {
            account_did: entry.account_did.clone(),
            kind: entry.kind.clone(),
            entity_ref: entry.entity_ref.clone(),
            collection: entry.collection.clone(),
            status: status.to_string(),
            at_uri: rref.map(|r| r.uri.clone()),
            at_cid: rref.map(|r| r.cid.clone()),
            attempt_count,
            error: error.map(str::to_string),
        };
        sync_history::record(self.store.pool(), &h, &Utc::now().to_rfc3339()).await?;
        Ok(())
    }

    /// Load (or, on first use, generate + publish) this installation's Ed25519 **device key**
    /// — the signing key that authenticates Edge→AppView calls (federated IBD and, later, the
    /// whole signed surface). The key seed lives in the OS keychain scoped to the signed-in
    /// DID; its public half is published once to the user's PDS as a
    /// [`DEVICE_KEY_COLLECTION`] record so the AppView (which ingests it via Jetstream) can
    /// verify our signatures. Idempotent: the record is keyed by its own `did:key`, so a
    /// re-publish overwrites rather than duplicates, and an already-present record is left
    /// alone. Errors [`AppError::NotAuthenticated`] when signed out.
    ///
    /// This does *not* wait for ingest — the signed AppView calls absorb the 403→200 lag with
    /// bounded retries (see the IBD client).
    pub async fn ensure_device_key(&self) -> Result<DeviceKey, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = DeviceKey::load_or_generate(KEYCHAIN_SERVICE, &did)?;

        // Publish the public key once. A public getRecord on the deterministic rkey tells us
        // whether it already exists; only create it when absent (keeps re-launches quiet).
        let rkey = key.record_rkey();
        let session = self.auth.tokens.load(&did)?.ok_or(AppError::NotAuthenticated)?;
        let client = PdsClient::from_session(self.auth.http.clone(), &session)?;
        let already_published = client.get_record(DEVICE_KEY_COLLECTION, &rkey).await.is_ok();
        if !already_published {
            let record = serde_json::json!({
                "publicKey": key.did_key(),
                "createdAt": Utc::now().to_rfc3339(),
            });
            let mut engine = self.sync_engine()?;
            engine.push_create_rkey(DEVICE_KEY_COLLECTION, record, &rkey).await?;
        }
        Ok(key)
    }

    /// Federated IBD — **Step 1**: fetch this account's pseudonymous match suggestions from
    /// the AppView (`GET /api/v1/ibd/suggestions`).
    ///
    /// The AppView mines our already-published `fed.*` records into a top-K candidate list;
    /// no genotypes leave the device here. The call is authenticated by signing
    /// `"ibd-poll\n<DID>\n<ts>"` with the device key (registered on first use). A 403 right
    /// after first-time registration means the AppView hasn't ingested the device-key record
    /// yet, so it's retried with exponential backoff.
    pub async fn ibd_suggestions(&self) -> Result<Vec<IbdSuggestion>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/ibd/suggestions", decodingus_appview_url());

        let mut attempt = 0u32;
        loop {
            let ts = Utc::now().timestamp().to_string();
            let sig = key.sign(&format!("ibd-poll\n{did}\n{ts}"));
            // reqwest URL-encodes query values, so the STANDARD-base64 sig (`+` `/` `=`) is
            // safely escaped.
            let resp = self
                .auth
                .http
                .get(&url)
                .query(&[("did", did.as_str()), ("ts", ts.as_str()), ("sig", sig.as_str())])
                .send()
                .await
                .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
            let status = resp.status();
            if status.is_success() {
                let body: serde_json::Value =
                    resp.json().await.map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
                return Ok(parse_ibd_suggestions(&body));
            }
            if status.as_u16() == 403 && attempt < DEVICE_KEY_INGEST_RETRIES {
                tokio::time::sleep(std::time::Duration::from_secs(1u64 << attempt)).await;
                attempt += 1;
                continue;
            }
            return Err(appview_status_error("ibd/suggestions", resp).await);
        }
    }

    /// Federated IBD — **Step 2**: request an introduction to a suggested candidate
    /// (`POST /api/v1/ibd/introduce`).
    ///
    /// Signs `"ibd-introduce\n<DID>\n<suggested_sample_guid>"` and posts
    /// `{ did, suggestedSampleGuid, signature }`. Returns the AppView's `request_uri` and
    /// status (`PENDING`). The downstream consent round-trip + key exchange are deferred
    /// (gated on the AppView's symmetric-blind counterpart discovery), so this only opens the
    /// request — it does not exchange any genetic data.
    pub async fn ibd_introduce(&self, suggested_sample_guid: &str) -> Result<IbdIntroResult, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let key = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/ibd/introduce", decodingus_appview_url());
        let sig = key.sign(&format!("ibd-introduce\n{did}\n{suggested_sample_guid}"));
        // The AppView's IntroduceBody deserializes plain snake_case (no serde rename), and
        // parses the guid as a UUID — send it verbatim from the suggestion.
        let body = serde_json::json!({
            "did": did,
            "suggested_sample_guid": suggested_sample_guid,
            "signature": sig,
        });
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&body)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error("ibd/introduce", resp).await);
        }
        let v: serde_json::Value =
            resp.json().await.map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        let request_uri = v
            .get("requestUri")
            .or_else(|| v.get("request_uri"))
            .and_then(|x| x.as_str())
            .unwrap_or_default()
            .to_string();
        let status = v.get("status").and_then(|x| x.as_str()).unwrap_or("PENDING").to_string();
        Ok(IbdIntroResult { request_uri, status })
    }

    // ---- IBD Phase 2: encrypted edge-to-edge exchange (D1 substrate) -------
    //
    // The AppView brokers discovery/consent + relays opaque ciphertext (never decrypts). These
    // wrap the `/api/v1/exchange/*` endpoints; the crypto (X25519/X3DH-lite/AES-GCM) lives in
    // `navigator_sync::exchange`. All calls are device-key-signed (no per-call OAuth).

    /// The signed-in account's X25519 identity key (load-or-generate), with its public half
    /// published to the AppView (`POST /exchange/key`, idempotent upsert) so partners can fetch it.
    pub async fn ensure_exchange_key(&self) -> Result<ExchangeKey, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ik = ExchangeKey::load_or_generate(KEYCHAIN_SERVICE, &did)?;
        let pub_b64 = ik.public_b64();
        let sig = dev.sign(&exchange::messages::publickey(&did, &pub_b64, None));
        let body = serde_json::json!({ "did": did, "x25519_pub": pub_b64, "signature": sig });
        let v = self.exchange_post("exchange/key", body).await?;
        let _ = v; // { did, status: "published" }
        Ok(ik)
    }

    /// Fetch a peer's published X25519 public key (STANDARD base64), or `None` if they haven't
    /// published one. Public read — no signature.
    pub async fn fetch_exchange_key(&self, did: &str) -> Result<Option<String>, AppError> {
        let url = format!("{}/api/v1/exchange/key", decodingus_appview_url());
        let resp = self
            .auth
            .http
            .get(&url)
            .query(&[("did", did)])
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if resp.status().as_u16() == 404 {
            return Ok(None);
        }
        if !resp.status().is_success() {
            return Err(appview_status_error("exchange/key", resp).await);
        }
        let v: serde_json::Value = resp.json().await.map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        Ok(v.get("x25519_pub").and_then(|x| x.as_str()).map(str::to_string))
    }

    /// Consent to (or decline) an exchange request. On mutual consent the AppView opens a session
    /// and returns its id.
    pub async fn exchange_consent(&self, request_uri: &str, given: bool) -> Result<ConsentOutcome, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let sig = dev.sign(&exchange::messages::consent(request_uri, &did, given));
        let body = serde_json::json!({
            "request_uri": request_uri,
            "consenting_did": did,
            "consent_given": given,
            "signature": sig,
        });
        let v = self.exchange_post("exchange/consent", body).await?;
        Ok(ConsentOutcome {
            status: v.get("status").and_then(|x| x.as_str()).unwrap_or("PENDING").to_string(),
            session_id: v.get("session_id").and_then(|x| x.as_str()).map(str::to_string),
        })
    }

    /// Poll for inbound (symmetric-blind) exchange requests awaiting this account's consent.
    pub async fn exchange_incoming(&self) -> Result<Vec<IncomingRequest>, AppError> {
        let v = self.exchange_get_poll("exchange/incoming", &[]).await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| IncomingRequest {
                        request_uri: i.get("request_uri").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        purpose: i.get("purpose").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        created_at: i.get("created_at").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Poll for consent-ready sessions (both parties consented; partner identity now revealed).
    pub async fn exchange_pending(&self) -> Result<Vec<ExchangeSessionInfo>, AppError> {
        let v = self.exchange_get_poll("exchange/pending", &[]).await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| ExchangeSessionInfo {
                        session_id: i.get("session_id").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        request_uri: i.get("request_uri").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        purpose: i.get("purpose").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        partner_did: i.get("partner_did").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        partner_key_uri: i.get("partner_key_uri").and_then(|x| x.as_str()).map(str::to_string),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Relay an opaque ciphertext `blob` to `to_did` in a session. The signed hash binds the blob to
    /// its routing (the broker stores ciphertext only). Returns the broker envelope id.
    pub async fn exchange_relay(&self, session_id: &str, to_did: &str, seq: i32, blob: &str) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let hash = exchange::blob_sha256_b64(blob).map_err(AppError::Sync)?;
        let sig = dev.sign(&exchange::messages::relay(session_id, &did, to_did, seq, &hash));
        let body = serde_json::json!({
            "session_id": session_id,
            "from_did": did,
            "to_did": to_did,
            "seq": seq,
            "blob": blob,
            "signature": sig,
        });
        let v = self.exchange_post("exchange/relay", body).await?;
        Ok(v.get("id").and_then(|x| x.as_i64()).unwrap_or_default())
    }

    /// Pull undelivered relay envelopes for a session (ordered by seq).
    pub async fn exchange_relay_pull(&self, session_id: &str) -> Result<Vec<RelayEnvelope>, AppError> {
        let v = self.exchange_get_poll("exchange/relay/pull", &[("session_id", session_id)]).await?;
        Ok(v.get("items")
            .and_then(|x| x.as_array())
            .map(|items| {
                items
                    .iter()
                    .map(|i| RelayEnvelope {
                        id: i.get("id").and_then(|x| x.as_i64()).unwrap_or_default(),
                        from_did: i.get("from_did").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                        seq: i.get("seq").and_then(|x| x.as_i64()).unwrap_or_default() as i32,
                        blob: i.get("blob").and_then(|x| x.as_str()).unwrap_or_default().to_string(),
                    })
                    .collect()
            })
            .unwrap_or_default())
    }

    /// Acknowledge a delivered relay envelope (the broker drops it).
    pub async fn exchange_relay_ack(&self, envelope_id: i64) -> Result<(), AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let sig = dev.sign(&exchange::messages::ack(&did, envelope_id));
        let body = serde_json::json!({ "envelope_id": envelope_id, "did": did, "signature": sig });
        self.exchange_post("exchange/ack", body).await.map(|_| ())
    }

    /// Establish a shared session key for a consent-ready session: publish/load our identity key,
    /// fetch the partner's, exchange ephemeral keys via the relay (handshake, seq 0), and derive the
    /// X3DH-lite session key. Polls the relay up to ~15s for the partner's handshake. The returned
    /// [`EstablishedSession`] then seals/opens payloads. (Live-only — needs a running AppView + the
    /// partner edge online to complete the handshake.)
    pub async fn open_exchange_session(&self, info: &ExchangeSessionInfo) -> Result<EstablishedSession, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let ik = self.ensure_exchange_key().await?;
        let partner_ik = self
            .fetch_exchange_key(&info.partner_did)
            .await?
            .ok_or_else(|| AppError::AppView(format!("partner {} has not published an X25519 key", info.partner_did)))?;

        let ek = exchange::EphemeralKey::generate();
        let hs = exchange::Envelope::handshake(&ek).to_blob().map_err(AppError::Sync)?;
        self.exchange_relay(&info.session_id, &info.partner_did, 0, &hs).await?;

        // Wait for the partner's handshake (seq 0 / a Handshake envelope), acking just it.
        let mut their_ek: Option<String> = None;
        for _ in 0..15 {
            for env in self.exchange_relay_pull(&info.session_id).await? {
                if let Ok(exchange::Envelope::Handshake { ek, .. }) = exchange::Envelope::from_blob(&env.blob) {
                    their_ek = Some(ek);
                    let _ = self.exchange_relay_ack(env.id).await;
                    break;
                }
            }
            if their_ek.is_some() {
                break;
            }
            tokio::time::sleep(std::time::Duration::from_secs(1)).await;
        }
        let their_ek = their_ek.ok_or_else(|| AppError::AppView("partner handshake not received (peer offline?)".into()))?;

        let key = exchange::derive_session_key(&ik, &ek, &partner_ik, &their_ek, exchange::role_is_a(&did, &info.partner_did))
            .map_err(AppError::Sync)?;
        Ok(EstablishedSession { session_id: info.session_id.clone(), partner_did: info.partner_did.clone(), key })
    }

    /// Seal `plaintext` and relay it on an established session (data starts at seq 1).
    pub async fn exchange_send(&self, session: &EstablishedSession, seq: i32, plaintext: &[u8]) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let aad = exchange::relay_aad(&session.session_id, &did, &session.partner_did, seq);
        let blob = exchange::seal(&session.key, &aad, plaintext).and_then(|e| e.to_blob()).map_err(AppError::Sync)?;
        self.exchange_relay(&session.session_id, &session.partner_did, seq, &blob).await
    }

    /// Pull + decrypt + ack the data payloads waiting on an established session (returns plaintexts
    /// in pull order). Non-data / undecryptable envelopes are left un-acked.
    pub async fn exchange_receive(&self, session: &EstablishedSession) -> Result<Vec<Vec<u8>>, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let mut out = Vec::new();
        for env in self.exchange_relay_pull(&session.session_id).await? {
            let Ok(parsed) = exchange::Envelope::from_blob(&env.blob) else { continue };
            // AAD binds the sender's routing: from = the partner (sender), to = us.
            let aad = exchange::relay_aad(&session.session_id, &env.from_did, &did, env.seq);
            if let Ok(pt) = exchange::open(&session.key, &aad, &parsed) {
                out.push(pt);
                let _ = self.exchange_relay_ack(env.id).await;
            }
        }
        Ok(out)
    }

    /// POST a JSON body to an `/api/v1/<path>` exchange endpoint, mapping non-2xx to an AppView error.
    async fn exchange_post(&self, path: &str, body: serde_json::Value) -> Result<serde_json::Value, AppError> {
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&body)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json().await.map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }

    /// Issue a device-key-signed `exchange-poll` GET to an `/api/v1/<path>` endpoint, with `extra`
    /// query params appended. Shared by incoming / pending / relay-pull.
    async fn exchange_get_poll(&self, path: &str, extra: &[(&str, &str)]) -> Result<serde_json::Value, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&exchange::messages::poll(&did, ts));
        let ts_s = ts.to_string();
        let mut query: Vec<(&str, &str)> = vec![("did", did.as_str()), ("ts", ts_s.as_str()), ("sig", sig.as_str())];
        query.extend_from_slice(extra);
        let resp = self
            .auth
            .http
            .get(&url)
            .query(&query)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json().await.map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }

    /// Publish the alignment's coverage summary to the signed-in account's PDS (with
    /// refresh-on-expiry and retry/backoff via [`AsyncSync`]).
    pub async fn publish_coverage(&self, alignment_id: i64) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.coverage_record(alignment_id).await?;
        self.enqueue_publish("coverage", &format!("alignment:{alignment_id}"), NS_ALIGNMENT, None, value).await
    }

    /// Publish every persisted ancestry estimate (admixture + PCA-GMM) for the alignment to the
    /// signed-in account's PDS — one populationBreakdown record per method. This is the researcher
    /// opt-in act for the ancestry section — anonymized population proportions only.
    pub async fn publish_ancestry(&self, alignment_id: i64) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        // One outbox row per method, keyed by method so re-publishing coalesces per estimate.
        let results = ancestry_result::list_for_alignment(self.store.pool(), alignment_id).await?;
        for r in &results {
            let value = serde_json::to_value(population_breakdown_record(r))?;
            let entity_ref = format!("ancestry:{alignment_id}:{}", r.method);
            self.enqueue_publish("ancestry", &entity_ref, NS_POPULATION_BREAKDOWN, None, value).await?;
        }
        Ok(())
    }

    /// Publish the anonymized biosample summary to the signed-in account's PDS.
    pub async fn publish_biosample(&self, biosample_guid: SampleGuid) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.biosample_record(biosample_guid).await?;
        self.enqueue_publish("biosample", &format!("biosample:{biosample_guid}"), NS_BIOSAMPLE, None, value).await
    }

    /// Publish a sequence-run characterization to the signed-in account's PDS.
    pub async fn publish_sequence_run(&self, run: &SequenceRun) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.sequence_run_record(run).await?;
        self.enqueue_publish("seqrun", &format!("seqrun:{}", run.id), NS_SEQUENCERUN, None, value).await
    }

    /// Publish the alignment's de-novo calls for `contig` to the signed-in account's PDS.
    pub async fn publish_variants(&self, alignment_id: i64, contig: &str) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.variants_record(alignment_id, contig).await?;
        let entity_ref = format!("variants:{alignment_id}:{contig}");
        self.enqueue_publish("variants", &entity_ref, PRIVATE_VARIANTS_COLLECTION, None, value).await
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

        // Vendor-aware tagging for VCFs: recognize FTDNA Big Y / Y Elite / YSEQ / mtFull from the
        // header + filename + sibling readme, and record the vendor label, a meaningful SourceType,
        // and the reference build (feeds Y/mt placement liftover). A generic VCF keeps the caller's
        // label/source_type. CSV imports are unchanged.
        let (source_label, source_type, reference_build) = if is_vcf {
            let (meta, contigs) = peek_vcf_header(path);
            let vendor = navigator_domain::vendorvcf::classify(&meta, &contigs, &label, sibling_readme(path).as_deref());
            let build = detect_vcf_build(&meta);
            if vendor.is_recognized() {
                (format!("{} ({})", vendor.display(), vcf_label_context(path, &label)), vendor.source_type(), build)
            } else {
                (label, source_type, build)
            }
        } else {
            (label, source_type, None)
        };

        let new = NewVariantSet { biosample_guid, source_label, source_type, reference_build, calls };
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
            reference_build: None,
            calls,
        };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// The build to emit a subject's BISDNA calls on: the first of its alignments whose
    /// reference build maps to a dictionary key, else `"hs1"` (the project default).
    async fn bisdna_target_build(&self, biosample_guid: SampleGuid) -> String {
        if let Ok(aligns) = alignment::list_for_biosample(self.store.pool(), biosample_guid).await {
            for a in &aligns {
                if let Some(key) = decodingus_build_key(&a.reference_build) {
                    return key.to_string();
                }
            }
        }
        "hs1".to_string()
    }

    /// Import a BISDNA chromo2 Y-SNP export. Each named marker is resolved to a locus via the
    /// Y-SNP dictionary on `build` (when `None`, the subject's alignment build, else `"hs1"`).
    /// Only **positive** (derived) calls become variant calls: a negative is not a variant, and
    /// [`reconciliation::reconcile_variants`] weights every stored call as a carried allele.
    /// `no_call`, back-mutated, and dictionary-unresolved markers are tallied but not emitted.
    /// The genotype is a QC cross-check only — the file's verdict (independent of the Illumina
    /// TOP strand) decides derived/ancestral. Stored as a `Chip`-weighted [`VariantSet`].
    pub async fn import_bisdna_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        build: Option<&str>,
    ) -> Result<BisdnaImportSummary, AppError> {
        let text = std::fs::read_to_string(path)?;
        let calls = bisdna::parse(&text).map_err(AppError::Import)?;
        let build = match build {
            Some(b) => b.to_string(),
            None => self.bisdna_target_build(biosample_guid).await,
        };

        let dict_dir = ysnp_dict::asset_dir();
        let dict = YsnpDictionary::load(&dict_dir).map_err(|e| {
            AppError::Import(format!(
                "{e}. Build the Y-SNP dictionary with scripts/ysnp-dictionary (expected under {})",
                dict_dir.display()
            ))
        })?;

        const UNRESOLVED_SAMPLE_CAP: usize = 25;
        let outcome = bisdna::resolve_calls(&calls, &dict, &build, UNRESOLVED_SAMPLE_CAP);

        let derived_calls = outcome.calls.len();
        let label =
            path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "BISDNA".into());

        // Also record an array QC summary so the chromo2 chip appears under Data Sources →
        // Chip / Array Profiles (the placeable per-SNP calls live in the variant set below; a
        // genotyping array legitimately has both a QC/provenance summary and its calls). BISDNA
        // is a Y-only haploid panel: every called marker is a Y marker, heterozygosity is n/a.
        let total = calls.len() as i64;
        let called = total - outcome.no_call as i64;
        let chip = NewChipProfile {
            biosample_guid,
            provider: "BISDNA".into(),
            chip_version: Some("chromo2".into()),
            summary: chipprofile::ChipSummary {
                total_markers_possible: total,
                total_markers_called: called,
                no_call_rate: if total > 0 { outcome.no_call as f64 / total as f64 } else { 0.0 },
                het_rate: None,
                y_markers_called: called,
                mt_markers_called: 0,
                autosomal_markers_called: 0,
            },
            source_file_name: Some(label.clone()),
            source_path: None, // BISDNA is a Y-only panel — no autosomal genotypes for ancestry
        };
        chip_profile::create(self.store.pool(), &chip).await?;

        let new = NewVariantSet {
            biosample_guid,
            source_label: label,
            source_type: SourceType::Chip,
            reference_build: Some(build.clone()),
            calls: outcome.calls,
        };
        let variant_set = variant_set::create(self.store.pool(), &new).await?;

        Ok(BisdnaImportSummary {
            variant_set,
            build,
            total_markers: calls.len(),
            derived_calls,
            ancestral: outcome.ancestral,
            no_call: outcome.no_call,
            back_mutated: outcome.back_mutated,
            unresolved: outcome.unresolved,
            unresolved_names: outcome.unresolved_names,
            strand_mismatches: outcome.strand_mismatches,
        })
    }

    /// All variant sets for a subject.
    pub async fn list_variant_sets(&self, biosample_guid: SampleGuid) -> Result<Vec<VariantSet>, AppError> {
        Ok(variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- chip / array profiles ---------------------------------------------

    /// Import a genotyping-array raw-data export (CSV/TSV) and store its QC summary.
    /// `provider` overrides vendor detection when given; `chip_version` is optional.
    /// Import a genotyping-array raw-data export and (1) store its QC summary as a [`ChipProfile`],
    /// (2) store the haploid Y/MT genotype rows as a `Chip`-source [`VariantSet`], and (3)
    /// best-effort place the Y (and, where present, mtDNA) haplogroup on import — the consumer-array
    /// counterpart to BISDNA's chromo2 path. 23andMe carries both Y and MT rows; AncestryDNA carries
    /// Y but no usable mtDNA. The stored observed bases flow through the same
    /// [`assign_y_bisdna`](Self::assign_y_bisdna) / [`assign_mt_chip`](Self::assign_mt_chip) +
    /// `assemble_assignment_robust` placement as BISDNA, with plus-strand reconciliation to the tree.
    /// Placement is best-effort: an unreachable tree (offline) leaves the calls stored for a later
    /// manual "Assign … (panel)" — it does not fail the import.
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
        let label = source_file_name.clone().unwrap_or_else(|| provider.clone());
        // Record the absolute path so ancestry-from-chip can re-read the autosomal genotypes later
        // (like alignments re-read bam_path). Canonicalize best-effort; fall back to the given path.
        let source_path = Some(
            std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf()).to_string_lossy().into_owned(),
        );
        let new = NewChipProfile {
            biosample_guid,
            provider: provider.clone(),
            chip_version,
            summary,
            source_file_name,
            source_path,
        };
        let profile = chip_profile::create(self.store.pool(), &new).await?;

        // Pull the haploid Y/MT genotype rows and store them as Chip-source variant calls so the
        // haplogroup placement (and later re-placement) has them without re-reading the file. The
        // observed allele goes in both `reference` and `alternate` (we don't know the ancestral);
        // the placement reads `alternate`.
        let haplo = chipprofile::haplo_calls(&text);
        if !haplo.is_empty() {
            let build = chipprofile::detect_build(&text);
            let (mut y_count, mut mt_count) = (0usize, 0usize);
            let mut variant_calls = Vec::with_capacity(haplo.len());
            for c in &haplo {
                let (contig, is_y) = match c.dna {
                    chipprofile::ChipDna::Y => ("chrY", true),
                    chipprofile::ChipDna::Mt => ("chrM", false),
                };
                let b = c.base.to_string();
                if let Some(call) = variants::snp_call(contig, c.position, &b, &b, Some(c.rsid.clone()), Some("1".into())) {
                    if is_y {
                        y_count += 1;
                    } else {
                        mt_count += 1;
                    }
                    variant_calls.push(call);
                }
            }
            let set = NewVariantSet {
                biosample_guid,
                source_label: format!("{label} Y/MT calls"),
                source_type: SourceType::Chip,
                reference_build: Some(build.clone()),
                calls: variant_calls,
            };
            variant_set::create(self.store.pool(), &set).await?;

            // Compute the haplogroups on import (best-effort; an offline tree just leaves the calls).
            if y_count > 0 {
                if let Err(e) = self.assign_y_bisdna(biosample_guid, Some(&build)).await {
                    eprintln!("chip Y placement deferred ({e})");
                }
            }
            // AncestryDNA's stray MT rows aren't a usable mtDNA panel — only place mtDNA when the
            // array carries a real MT marker set (23andMe has thousands; the threshold filters noise).
            const MIN_MT_CALLS: usize = 20;
            if mt_count >= MIN_MT_CALLS {
                if let Err(e) = self.assign_mt_chip(biosample_guid).await {
                    eprintln!("chip mtDNA placement deferred ({e})");
                }
            }
        }

        Ok(profile)
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
        let seq = mtdna_store::create(self.store.pool(), &new).await?;

        // Derive rCRS-relative variants and persist them, so an mtDNA FASTA yields a variant set on
        // import (not only on the on-demand "show mutations" view) — like a chip/VCF import does.
        let derived = navigator_analysis::mtvariants::derive(navigator_analysis::mtvariants::rcrs(), &seq.sequence);
        if !derived.is_empty() {
            let label = mt_vendor_label(seq.source_file_name.as_deref(), seq.defline.as_deref());
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
            let set = NewVariantSet {
                biosample_guid,
                source_label: format!("{label} ({} variants vs rCRS)", derived.len()),
                // A full-mtDNA consensus is authoritative for its calls (gold-standard weight).
                source_type: variants::SourceType::Sanger,
                reference_build: None, // calls are rCRS-relative (contig "rCRS"), not a nuclear build
                calls,
            };
            // Best-effort: a variant-set hiccup must not lose the stored sequence.
            let _ = variant_set::create(self.store.pool(), &set).await;
        }

        // Haplogroup placement is intentionally NOT run here: it needs the mt haplotree (network),
        // and coupling a deterministic import to a network fetch is what the alignment import
        // deliberately avoids too. The mtDNA tab's "Assign mtDNA haplogroup" places it on demand.
        Ok(seq)
    }

    /// All mtDNA sequences for a subject.
    pub async fn list_mtdna_sequences(&self, biosample_guid: SampleGuid) -> Result<Vec<MtdnaSequence>, AppError> {
        Ok(mtdna_store::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Derive mtDNA variants for a stored sequence by comparing it to an rCRS reference
    /// FASTA, and save them as a variant set (contig `rCRS`) so they appear alongside the
    /// subject's other variants. The reference is validated as an mtDNA FASTA.
    /// The mtDNA mutation list for a stored sequence: variants relative to the **bundled** rCRS
    /// (NC_012920.1), via banded alignment — substitutions, insertions, and deletions in standard
    /// mtDNA notation. On-demand (one ~16.5 kb alignment), not stored. The classic mtDNA result.
    pub async fn mtdna_variants(&self, mtdna_id: i64) -> Result<Vec<MtVariant>, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;
        Ok(navigator_analysis::mtvariants::derive(navigator_analysis::mtvariants::rcrs(), &seq.sequence))
    }

    // ---- result exports (gap §6) -------------------------------------------

    /// Format a cached result as a shareable file body (TSV / HTML / BED). The UI writes the
    /// returned string to the user-chosen path. Errors when the source result hasn't been computed
    /// yet (`NotFound`). [`ExportRequest::CallableBed`] re-walks the BAM (no cached intervals).
    pub async fn export_content(&self, req: &ExportRequest) -> Result<String, AppError> {
        match req {
            ExportRequest::CoverageTsv(id) => Ok(export::coverage_tsv(&self.require_coverage(*id).await?)),
            ExportRequest::CoverageHtml(id) => {
                Ok(export::coverage_html(&self.require_coverage(*id).await?, &format!("alignment {id}")))
            }
            ExportRequest::ReadMetricsTsv(id) => {
                let m = self
                    .cached_read_metrics(*id)
                    .await?
                    .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("read metrics for alignment {id}"))))?;
                Ok(export::read_metrics_tsv(&m))
            }
            ExportRequest::AncestryTsv(id) => Ok(export::ancestry_tsv(&self.require_ancestry(*id).await?)),
            ExportRequest::AncestryHtml(id) => Ok(export::ancestry_html(&self.require_ancestry(*id).await?)),
            ExportRequest::MtdnaTsv(id) => Ok(export::mtdna_variants_tsv(&self.mtdna_variants(*id).await?)),
            ExportRequest::CallableBed(id) => {
                let per_contig = self.callable_intervals_all(*id).await?;
                Ok(export::callable_bed(&per_contig))
            }
            ExportRequest::DiploidVcf(id) => self.diploid_vcf_genome(*id).await,
        }
    }

    async fn require_coverage(&self, alignment_id: i64) -> Result<CoverageResult, AppError> {
        self.cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))
    }

    async fn require_ancestry(&self, alignment_id: i64) -> Result<AncestryResult, AppError> {
        self.ancestry_for_alignment(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("ancestry for alignment {alignment_id}"))))
    }

    /// Walk each analyzed contig for its CALLABLE intervals (BED export). Re-reads the BAM — the
    /// coverage artifact stores only per-contig callable *counts*, not the intervals. Uses the
    /// contig list from the cached coverage result, so coverage must have been run first.
    async fn callable_intervals_all(&self, alignment_id: i64) -> Result<Vec<(String, Vec<(i64, i64)>)>, AppError> {
        let cov = self.require_coverage(alignment_id).await?;
        let contigs: Vec<String> = cov.contig_coverage_stats.iter().map(|s| s.contig.clone()).collect();
        let (bam, reference) = self.alignment_paths(alignment_id).await?;
        let out = tokio::task::spawn_blocking(move || {
            let mut params = CallableLociParams::default();
            if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(&bam, reference.as_deref()) {
                params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            }
            let mut per_contig = Vec::new();
            for contig in contigs {
                // A contig with no aligned reads / bad region just contributes no intervals.
                let intervals =
                    coverage::callable_intervals(&bam, &contig, &params, 1, reference.as_deref()).unwrap_or_default();
                if !intervals.is_empty() {
                    per_contig.push((contig, intervals));
                }
            }
            per_contig
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?;
        Ok(out)
    }

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
            reference_build: None,
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
        self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, call, None).await
    }

    /// Like [`record_haplogroup_call`](Self::record_haplogroup_call) but stamps the input
    /// fingerprint (file + tree content hashes) so a later run can skip re-scoring.
    async fn record_haplogroup_call_fp(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        call: &RunHaplogroupCall,
        fingerprint: Option<&str>,
    ) -> Result<(), AppError> {
        haplogroup_call::upsert(self.store.pool(), biosample_guid, dna_type, source_key, call, fingerprint).await?;
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
        self.record_call_fp(biosample_guid, dna_type, source_key, source_label, assignment, None).await
    }

    /// Like [`record_call`](Self::record_call) but stamps the input fingerprint.
    async fn record_call_fp(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        source_key: &str,
        source_label: String,
        assignment: &HaploAssignment,
        fingerprint: Option<&str>,
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
            self.record_haplogroup_call_fp(biosample_guid, dna_type, source_key, &call, fingerprint).await?;
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

    /// Donor-level Y and mtDNA terminal haplogroups for **every** subject, for the subjects
    /// list. Reconciles each subject's recorded calls (and applies any manual override) in
    /// memory from two bulk queries. `(guid → (Y terminal, mt terminal))`; either is `None`
    /// when nothing is recorded.
    pub async fn haplogroup_terminals(
        &self,
    ) -> Result<HashMap<SampleGuid, (Option<String>, Option<String>)>, AppError> {
        let mut groups: HashMap<(SampleGuid, DnaType), Vec<RunHaplogroupCall>> = HashMap::new();
        for (guid, dna_type, call) in haplogroup_call::list_all(self.store.pool()).await? {
            groups.entry((guid, dna_type)).or_default().push(call);
        }
        let mut out: HashMap<SampleGuid, (Option<String>, Option<String>)> = HashMap::new();
        for ((guid, dna_type), calls) in groups {
            if let Some(c) = reconciliation::reconcile(&calls) {
                let entry = out.entry(guid).or_default();
                match dna_type {
                    DnaType::Y => entry.0 = Some(c.haplogroup),
                    DnaType::Mt => entry.1 = Some(c.haplogroup),
                }
            }
        }
        // Manual overrides win over the reconciled terminal.
        for (guid, dna_type, hg) in recon_store::list_all_overrides(self.store.pool()).await? {
            let entry = out.entry(guid).or_default();
            match dna_type {
                DnaType::Y => entry.0 = Some(hg),
                DnaType::Mt => entry.1 = Some(hg),
            }
        }
        Ok(out)
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

    /// The persisted consensus-profile snapshot for a subject + DNA type, if built — cheap (no
    /// genotyping). The shared loader behind [`cached_y_profile`](Self::cached_y_profile); the future
    /// mtDNA / autosomal tabs reuse it with a different [`DnaType`]. `None` until a build runs.
    pub async fn cached_consensus_profile(&self, biosample_guid: SampleGuid, dna_type: DnaType) -> Result<Option<ConsensusProfile>, AppError> {
        match navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, dna_type.as_str()).await? {
            Some(row) => Ok(Some(serde_json::from_str(&row.payload)?)),
            None => Ok(None),
        }
    }

    /// Persist a reconciled consensus snapshot — the low-level row writer shared by every DNA type
    /// (Y / mt key on [`DnaType`], autosomal keys on `"Auto"`; the payload is whatever profile shape
    /// that type uses). The scalar columns mirror the summary header for quick listing.
    #[allow(clippy::too_many_arguments)]
    async fn persist_consensus_row(
        &self,
        biosample_guid: SampleGuid,
        dna_type: &str,
        consensus_label: Option<String>,
        summary: &navigator_domain::consensus::ConsensusSummary,
        source_count: usize,
        tree_provider: Option<String>,
        payload: String,
    ) -> Result<(), AppError> {
        let stored = navigator_store::consensus_profile::StoredConsensusProfile {
            biosample_guid: biosample_guid.0.to_string(),
            dna_type: dna_type.to_string(),
            consensus_label,
            overall_confidence: summary.overall_confidence,
            source_count: source_count as i64,
            total: summary.total as i64,
            confirmed: summary.confirmed as i64,
            novel: summary.novel as i64,
            conflict: summary.conflict as i64,
            single_source: summary.single_source as i64,
            tree_provider,
            payload,
            last_reconciled_at: Utc::now().to_rfc3339(),
        };
        navigator_store::consensus_profile::upsert(self.store.pool(), &stored).await?;
        Ok(())
    }

    /// Persist a Y/mt [`ConsensusProfile`] snapshot via [`persist_consensus_row`].
    async fn persist_consensus_profile(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
        profile: &ConsensusProfile,
        tree_provider: Option<String>,
    ) -> Result<(), AppError> {
        self.persist_consensus_row(
            biosample_guid,
            dna_type.as_str(),
            profile.terminal.clone(),
            &profile.summary,
            profile.sources.len(),
            tree_provider,
            serde_json::to_string(profile)?,
        )
        .await
    }

    /// The persisted Y-profile snapshot for a subject, if one has been built — cheap (no
    /// genotyping). `None` until [`build_y_profile`](Self::build_y_profile) runs.
    pub async fn cached_y_profile(&self, biosample_guid: SampleGuid) -> Result<Option<YProfile>, AppError> {
        self.cached_consensus_profile(biosample_guid, DnaType::Y).await
    }

    /// Build (and persist) the multi-source Y-variant profile: reconcile each Y-bearing source's
    /// per-SNP calls — every alignment's haplogroup placement, the combined chip/BISDNA placement,
    /// and the private-Y bucket — into one concordance view (confirmed / novel / conflict /
    /// single-source per SNP, with per-source provenance + per-observation quality weighting).
    /// Expensive (re-genotypes each alignment), so it's an explicit action; the result is persisted
    /// so [`cached_y_profile`](Self::cached_y_profile) reloads it instantly. Sources without Y data
    /// are skipped.
    pub async fn build_y_profile(&self, biosample_guid: SampleGuid) -> Result<YProfile, AppError> {
        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();

        // One source per alignment — a *fresh* placement (the cached terminal-only path lacks the
        // per-SNP branch evidence we reconcile here). Expensive; this is why the profile is built
        // on explicit request. Alignments that error / lack chrY are skipped.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok(assignment) = self.y_assignment_full(a.id).await else { continue };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
            }
        }

        // One source *per chip/BISDNA panel* (a distinct VariantSet per import — 23andMe,
        // AncestryDNA, BISDNA chromo2, …), so the profile shows which test confirmed each SNP and a
        // single mistyped panel surfaces as a conflict rather than being averaged into "consumer tests".
        let vsets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let chip_sets: Vec<&VariantSet> = vsets.iter().filter(|s| s.source_type == SourceType::Chip).collect();
        if !chip_sets.is_empty() {
            // Resolve the placement build once (a chip set's stored build, else the alignment's), and
            // fetch the tree once for all panels.
            let build = chip_sets
                .iter()
                .find_map(|s| s.reference_build.clone())
                .unwrap_or(self.bisdna_target_build(biosample_guid).await);
            if let Ok(tree) = self.chip_y_tree(&build).await {
                for set in &chip_sets {
                    let calls: HashMap<i64, char> = set
                        .calls
                        .iter()
                        .filter(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"))
                        .filter_map(|c| c.alternate.chars().next().map(|b| (c.position, b.to_ascii_uppercase())))
                        .collect();
                    if calls.is_empty() {
                        continue;
                    }
                    let assignment = Self::place_chip_panel(&tree, calls);
                    let obs = snp_obs_from_assignment(&assignment, true);
                    if !obs.is_empty() {
                        sources.push((set.source_label.clone(), SourceType::Chip, obs));
                    }
                }
            }
        }

        // Private-Y union: off-path / novel calls (not in the tree).
        if let Some(bucket) = self.donor_private_y(biosample_guid).await? {
            let obs: Vec<YObsInput> = bucket
                .variants
                .iter()
                .map(|v| {
                    let name = match &v.class {
                        PrivateClass::OffPathKnown(n) => n.clone(),
                        PrivateClass::Novel => String::new(), // keyed by position
                    };
                    let mut o = YObsInput::snp(name, v.position, v.reference.to_string(), v.alternate.to_string(), YState::Derived, false);
                    // De-novo calls carry read depth; a structural-region (palindrome/amplicon) call
                    // is paralog-suspect → down-weight via the region modifier.
                    o.depth = Some(v.depth);
                    // Down-weight by the structural-region quality modifier (palindrome 0.4,
                    // ampliconic 0.3, heterochromatin/centromere 0.1…); unique sequence = 1.0.
                    o.region_modifier = v.region.map_or(1.0, |c| c.modifier());
                    o
                })
                .collect();
            if !obs.is_empty() {
                sources.push(("private".to_string(), SourceType::WgsShortRead, obs));
            }
        }

        // Provenance: one entry per contributing source (label, type, SNP count).
        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary { label: label.clone(), source_type: *st, variant_count: obs.len() })
            .collect();

        let variants = yprofile::reconcile_y(&sources);
        let summary = yprofile::summarize(&variants);
        let terminal = self
            .haplogroup_consensus(biosample_guid, DnaType::Y)
            .await?
            .map(|c| c.haplogroup);
        let profile = ConsensusProfile { variants, summary, terminal, sources: source_summaries };

        // Persist the snapshot (keyed dna_type='Y') so the tab reloads it without re-genotyping.
        let provider = Some(match y_tree_provider() {
            YTreeProvider::DecodingUs => "decodingus".to_string(),
            YTreeProvider::Ftdna => "ftdna".to_string(),
        });
        self.persist_consensus_profile(biosample_guid, DnaType::Y, &profile, provider).await?;
        Ok(profile)
    }

    /// Fresh mtDNA placement against the FTDNA mt tree (chrM) with full branch evidence — the mt
    /// counterpart to [`y_assignment_full`](Self::y_assignment_full). Bypasses
    /// [`assign_mtdna_haplogroup_from_alignment`](Self::assign_mtdna_haplogroup_from_alignment)'s
    /// cached terminal-only path (which has `branches: []`) so the consensus has per-mutation evidence.
    async fn mt_assignment_full(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        self.assign_haplogroup_from_alignment(alignment_id, "chrM", &tree_json).await
    }

    /// The persisted mtDNA consensus-profile snapshot for a subject, if one has been built — cheap
    /// (no genotyping). `None` until [`build_mt_profile`](Self::build_mt_profile) runs.
    pub async fn cached_mt_profile(&self, biosample_guid: SampleGuid) -> Result<Option<ConsensusProfile>, AppError> {
        self.cached_consensus_profile(biosample_guid, DnaType::Mt).await
    }

    /// Build (and persist) the multi-source mtDNA consensus profile — the mtDNA adapter over the
    /// generic [`navigator_domain::consensus`] engine. Reconciles each mt-bearing source's
    /// defining-mutation calls (every alignment's chrM placement, each imported mtDNA FASTA
    /// sequence's placement, and the combined chip mtDNA placement) into one concordance view,
    /// keyed by phylotree **mutation name** (rCRS-coordinate, build-independent). Persisted with
    /// `dna_type='Mt'` so [`cached_mt_profile`](Self::cached_mt_profile) reloads it instantly.
    /// Expensive (re-places each alignment's chrM), so it's an explicit action; mt-less sources skip.
    pub async fn build_mt_profile(&self, biosample_guid: SampleGuid) -> Result<ConsensusProfile, AppError> {
        let mut sources: Vec<(String, SourceType, Vec<YObsInput>)> = Vec::new();

        // One source per alignment with chrM — a fresh placement (branch evidence, not the cached
        // terminal). Alignments that error / lack chrM are skipped.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            let Ok(assignment) = self.mt_assignment_full(a.id).await else { continue };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
            }
        }

        // One source per imported mtDNA FASTA sequence (FTDNA mtFull / YSEQ) — a finished consensus
        // sequence we ingested, so weight it as `Imported` (provenance/method not ours to vouch for).
        let seqs = self.list_mtdna_sequences(biosample_guid).await?;
        for s in &seqs {
            let Ok(assignment) = self.assign_mtdna_haplogroup(s.id).await else { continue };
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                let vendor = mt_vendor_label(s.source_file_name.as_deref(), s.defline.as_deref());
                sources.push((format!("{vendor} (mt seq #{})", s.id), SourceType::Imported, obs));
            }
        }

        // The combined chip mtDNA panel (23andMe carries a sparse mt panel). One source — the
        // per-panel split is a deferred follow-on (same as the Y profile).
        if let Ok(assignment) = self.assign_mt_chip(biosample_guid).await {
            let obs = snp_obs_from_assignment(&assignment, true);
            if !obs.is_empty() {
                sources.push(("Chip mtDNA panel".to_string(), SourceType::Chip, obs));
            }
        }

        // Provenance: one entry per contributing source (label, type, mutation count).
        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary { label: label.clone(), source_type: *st, variant_count: obs.len() })
            .collect();

        let variants = yprofile::reconcile_y(&sources);
        let summary = yprofile::summarize(&variants);
        let terminal = self.haplogroup_consensus(biosample_guid, DnaType::Mt).await?.map(|c| c.haplogroup);
        let profile = ConsensusProfile { variants, summary, terminal, sources: source_summaries };

        // Persist (keyed dna_type='Mt'); the mt tree is FTDNA-sourced.
        self.persist_consensus_profile(biosample_guid, DnaType::Mt, &profile, Some("ftdna".to_string())).await?;
        Ok(profile)
    }

    /// The persisted autosomal consensus-profile snapshot for a subject, if built — cheap (no
    /// genotyping). `None` until [`build_autosomal_profile`](Self::build_autosomal_profile) runs.
    pub async fn cached_autosomal_profile(&self, biosample_guid: SampleGuid) -> Result<Option<DiploidProfile>, AppError> {
        match navigator_store::consensus_profile::get(self.store.pool(), biosample_guid, "Auto").await? {
            Some(row) => Ok(Some(serde_json::from_str(&row.payload)?)),
            None => Ok(None),
        }
    }

    /// Build (and persist) the multi-source **autosomal** consensus profile — the diploid (0/1/2)
    /// adapter over the generic [`navigator_domain::consensus`] engine. Genotypes every WGS alignment
    /// and imported chip over the canonical CHM13 **IBD panel** ([`ibd_panel_dosages`](Self::ibd_panel_dosages))
    /// and reconciles the per-site dosages into a voted genotype (confirmed where sources agree,
    /// conflict where they don't), keyed by rsID. Persisted with `dna_type='Auto'`. Requires the IBD
    /// panel asset (built with `panelbuild ibd-panel`); errors if it's missing.
    pub async fn build_autosomal_profile(&self, biosample_guid: SampleGuid) -> Result<DiploidProfile, AppError> {
        use navigator_domain::consensus::{reconcile_diploid, summarize_diploid, DiploidObs};

        let to_obs = |gts: Vec<SiteGenotype>| -> Vec<DiploidObs> {
            gts.into_iter()
                .map(|g| DiploidObs {
                    name: g.name,
                    contig: g.contig,
                    position: g.position,
                    reference: g.reference_allele,
                    alternate: g.alternate_allele,
                    dosage: g.dosage as i8,
                    depth: (g.depth > 0).then_some(g.depth),
                })
                .collect()
        };

        let mut sources: Vec<(String, SourceType, Vec<DiploidObs>)> = Vec::new();
        // Remember the last source error: if *every* source fails (e.g. the panel asset is missing),
        // surface it rather than silently returning an empty profile; a one-off per-source failure
        // (a chip with no stored raw file, an alignment lacking a BAM) is just skipped.
        let mut last_err: Option<AppError> = None;

        // One source per WGS alignment (panel-genotyped, cached per alignment). The IBD panel is
        // CHM13-coordinate and the alignment path has no liftover, so only CHM13 alignments can be
        // genotyped directly — non-CHM13 builds reach the panel via the chip path (multi-build
        // coordinates) or a future lift, and are skipped here rather than yielding wrong-locus calls.
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        for a in &alignments {
            if !matches!(canonical_build(&a.reference_build), Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs)) {
                continue;
            }
            match self.ibd_panel_dosages(IbdSource::Alignment(a.id)).await {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((format!("aln #{} · {}", a.id, a.aligner), SourceType::WgsShortRead, obs));
                    }
                }
                Err(e) => last_err = Some(e),
            }
        }

        // One source per imported chip (resolved to canonical panel dosages, no alignment needed).
        let chips = self.list_chip_profiles(biosample_guid).await?;
        for c in &chips {
            match self.ibd_panel_dosages(IbdSource::Chip(c.id)).await {
                Ok(gts) => {
                    let obs = to_obs(gts);
                    if !obs.is_empty() {
                        sources.push((format!("{} (chip #{})", c.provider, c.id), SourceType::Chip, obs));
                    }
                }
                Err(e) => last_err = Some(e),
            }
        }

        if sources.is_empty() {
            if let Some(e) = last_err {
                return Err(e); // e.g. the IBD panel asset isn't built yet
            }
        }

        let source_summaries: Vec<YSourceSummary> = sources
            .iter()
            .map(|(label, st, obs)| YSourceSummary { label: label.clone(), source_type: *st, variant_count: obs.len() })
            .collect();
        let variants = reconcile_diploid(&sources);
        let summary = summarize_diploid(&variants);
        let profile = DiploidProfile { variants, summary, sources: source_summaries };

        // Persist (keyed dna_type='Auto'; no lineage label, no tree provider).
        self.persist_consensus_row(biosample_guid, "Auto", None, &profile.summary, profile.sources.len(), None, serde_json::to_string(&profile)?)
            .await?;
        Ok(profile)
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
    ) -> Result<(), AppError> {
        self.require_account()?; // auth check before touching the DB
        let value = self.reconciliation_record(biosample_guid, dna_type, heteroplasmy, identity).await?;
        let entity_ref = format!("reconciliation:{biosample_guid}:{dna_type:?}");
        self.enqueue_publish("reconciliation", &entity_ref, HAPLOGROUP_RECONCILIATION_COLLECTION, None, value).await
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

    /// DecodingUs Y-DNA tree-with-variants JSON from our AppView (`/api/v1/y-tree/full`),
    /// host from [`decodingus_appview_url`]. On-disk cached like the FTDNA tree.
    async fn fetch_decodingus_y_tree(&self) -> Result<String, AppError> {
        let url = format!("{}/api/v1/y-tree/full", decodingus_appview_url());
        self.fetch_tree(&url, "decodingus-ytree.json").await
    }

    /// FTDNA Y-DNA haplotree JSON, from the on-disk cache or freshly downloaded + cached.
    async fn fetch_ftdna_y_tree(&self) -> Result<String, AppError> {
        self.fetch_tree("https://www.familytreedna.com/public/y-dna-haplotree/get", "ftdna-ytree.json")
            .await
    }

    /// The AppView's full instrument→lab map (`GET /api/v1/sequencer/lab-instruments`), on-disk
    /// cached like the trees (7-day TTL + offline fallback). Looked up locally so a batch import
    /// makes one network call, not one per sample.
    async fn fetch_lab_instruments(&self) -> Result<Vec<SequencerLabInfo>, AppError> {
        let url = format!("{}/api/v1/sequencer/lab-instruments", decodingus_appview_url());
        let json = self.fetch_tree(&url, "sequencer-lab-instruments.json").await?;
        serde_json::from_str(&json).map_err(|e| AppError::Import(format!("parsing lab-instruments: {e}")))
    }

    /// Resolve an instrument id to a lab display name via the AppView (cached). Normalizes the
    /// returned name to the local [`labs`] catalog's canonical display name when it matches.
    /// `None` if the instrument has no association or the AppView is unreachable (best-effort).
    pub async fn lookup_lab_by_instrument(&self, instrument_id: &str) -> Option<String> {
        let id = instrument_id.trim();
        if id.is_empty() {
            return None;
        }
        let list = self.fetch_lab_instruments().await.ok()?;
        let raw = list.into_iter().find(|l| l.instrument_id == id)?.lab_name;
        Some(navigator_domain::labs::find(&raw).map(|l| l.display_name.to_string()).unwrap_or(raw))
    }

    /// Resolve the sequencing lab for every run that has an inferred `instrument_id` but no facility
    /// yet, via the AppView (one cached fetch). Best-effort; returns how many were filled. Run after
    /// import and on startup so pre-existing runs pick up newly-seeded associations.
    pub async fn backfill_run_labs(&self) -> Result<usize, AppError> {
        // One network/cache fetch, then resolve locally.
        let Ok(list) = self.fetch_lab_instruments().await else { return Ok(0) };
        let by_instrument: HashMap<&str, &str> =
            list.iter().map(|l| (l.instrument_id.as_str(), l.lab_name.as_str())).collect();
        let mut filled = 0usize;
        for biosample in biosample::list_all(self.store.pool()).await? {
            for run in sequence_run::list_for_biosample(self.store.pool(), biosample.guid).await? {
                if run.sequencing_facility.is_some() {
                    continue;
                }
                let Some(inst) = run.instrument_id.as_deref() else { continue };
                if let Some(raw) = by_instrument.get(inst.trim()) {
                    let lab = navigator_domain::labs::find(raw).map(|l| l.display_name).unwrap_or(raw);
                    if sequence_run::set_facility(self.store.pool(), run.id, lab).await.unwrap_or(false) {
                        filled += 1;
                    }
                }
            }
        }
        Ok(filled)
    }

    /// A cached-or-downloaded haplotree JSON. The on-disk cache has a **7-day life** (see
    /// [`TREE_CACHE_TTL`]): a fresh cache short-circuits the network; a stale or missing cache
    /// triggers a re-download (and refresh). If the re-download fails (e.g. the AppView is
    /// unreachable) but a stale copy exists, the stale copy is used rather than failing — so the
    /// app keeps working offline, just on an older tree. (A server-side ETag/version would let us
    /// revalidate without a full re-download; tracked as an AppView backlog item.)
    async fn fetch_tree(&self, url: &str, cache_file: &str) -> Result<String, AppError> {
        let path = tree_cache_path(cache_file);
        let cached = std::fs::read_to_string(&path).ok().filter(|c| !c.trim().is_empty());
        if let Some(cached) = &cached {
            if tree_cache_is_fresh(&path) {
                return Ok(cached.clone());
            }
        }
        // Stale or absent → (re)download, falling back to a stale copy on network failure.
        let downloaded = self
            .auth
            .http
            .get(url)
            .send()
            .await
            .and_then(|r| r.error_for_status())
            .map_err(|e| AppError::Import(format!("downloading {url}: {e}")));
        match downloaded {
            Ok(resp) => {
                let body = resp.text().await.map_err(|e| AppError::Import(format!("reading {url}: {e}")))?;
                if let Some(parent) = path.parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                let _ = std::fs::write(&path, &body);
                Ok(body)
            }
            Err(e) => match cached {
                Some(stale) => {
                    eprintln!("tree refresh failed ({e}); using the cached copy at {}", path.display());
                    Ok(stale)
                }
                None => Err(e),
            },
        }
    }

    /// Assign an mtDNA haplogroup directly from an alignment's chrM reads (FTDNA mt tree),
    /// the BAM-based counterpart to [`assign_mtdna_haplogroup`]. Requires a GRCh38/rCRS
    /// chrM (the tree is in rCRS coordinates).
    pub async fn assign_mtdna_haplogroup_from_alignment(
        &self,
        alignment_id: i64,
    ) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();
        let source_key = format!("aln:{alignment_id}:mt");
        let tree_json = self.fetch_ftdna_mt_tree().await?;

        // Cache: skip re-scoring when the file and the mt tree are unchanged.
        let fingerprint = self
            .alignment_content_hash(alignment_id)
            .await
            .ok()
            .map(|file_hash| format!("f:{}|mt:{}", &file_hash[..16], &sha256_str(&tree_json)[..16]));
        if let (Some(bio), Some(fp)) = (bio, fingerprint.as_deref()) {
            if haplogroup_call::stored_fingerprint(self.store.pool(), bio, DnaType::Mt, &source_key).await?.as_deref()
                == Some(fp)
            {
                if let Some(call) = haplogroup_call::get_one(self.store.pool(), bio, DnaType::Mt, &source_key).await? {
                    return Ok(assignment_from_call(&call));
                }
            }
        }

        let assignment = self.assign_haplogroup_from_alignment(alignment_id, "chrM", &tree_json).await?;
        if let Some(bio) = bio {
            self.record_call_fp(bio, DnaType::Mt, &source_key, format!("aln #{alignment_id} mtDNA"), &assignment, fingerprint.as_deref()).await?;
        }
        Ok(assignment)
    }

    /// mtDNA assignment + per-SNP lineage evidence (for exact GRCh38-vs-CHM13 comparison).
    pub async fn assign_mtdna_haplogroup_detail(
        &self,
        alignment_id: i64,
    ) -> Result<(HaploAssignment, Vec<SnpEvidence>, HashMap<i64, char>), AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        self.assign_haplogroup_detail(alignment_id, "chrM", &tree_json).await
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

    /// Estimate the donor's ancestry for an alignment by the allele-frequency likelihood: load
    /// the (build-matched) AIMs panel, genotype the sample at its sites with the GL caller, and
    /// score each super-population's binomial likelihood. Persists the result; returns it for
    /// display. Requires a recorded BAM/CRAM and a resolvable reference (CRAM/genotyping).
    /// Estimate autosomal ancestry from the subject's **consensus** — no BAM genotyping. Reads the
    /// cached autosomal [`DiploidProfile`] (reconciled 0/1/2 dosages over the probe panel, pooled
    /// across all WGS + chip sources), bridges it to genotypes, and runs the same estimators as the
    /// per-alignment path used to. Persisted under the consensus pseudo-source
    /// ([`CONSENSUS_SOURCE_ID`]). Errors if the autosomal consensus hasn't been built yet.
    pub async fn estimate_ancestry_from_consensus(&self, biosample_guid: SampleGuid) -> Result<AncestryResult, AppError> {
        let profile = self
            .cached_autosomal_profile(biosample_guid)
            .await?
            .ok_or_else(|| AppError::Import("build the autosomal consensus first (Autosomal tab) before estimating ancestry".into()))?;
        let genotypes = consensus_genotypes(&profile);

        // The consensus is canonical CHM13; the AIM freq / PCA assets are keyed by (contig,pos) there.
        let build = ReferenceBuild::Chm13v2;
        let reference_version = "chm13v2.0".to_string();
        let panel_path = ancestry_panel_path(build);
        let panel_bytes =
            read_verified_asset(build, &panel_path)?.ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        let optional = |path: PathBuf| read_verified_asset(build, &path).unwrap_or_else(|e| { eprintln!("{e}"); None });
        let pca_bytes = optional(ancestry_pca_path(build));
        let ancient_pca_bytes = optional(ancestry_pca_ancient_path(build));
        let fine_bytes = optional(ancestry_freq_global_path(build));

        let (result, pca_gmm, nmonte, fine) = tokio::task::spawn_blocking(move || {
            let mut result = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let fine = fine_bytes
                .and_then(|b| ancestry_analysis::AncestryPanel::from_bytes(&b).ok())
                .map(|fp| ancestry_analysis::estimate_fine_admixture(&genotypes, &fp, &reference_version));
            let modern_pca = pca_bytes.and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok());
            if let Some(pca) = &modern_pca {
                result.pca_coordinates = Some(ancestry_analysis::project_pca(&genotypes, pca));
            }
            let gmm_pca = ancient_pca_bytes
                .and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok())
                .or(modern_pca);
            let (pca_gmm, nmonte) = match &gmm_pca {
                Some(pca) => (
                    Some(ancestry_analysis::estimate_pca_gmm(&genotypes, pca, &reference_version)),
                    Some(ancestry_analysis::estimate_nmonte(&genotypes, pca, &reference_version)),
                ),
                None => (None, None),
            };
            (result, pca_gmm, nmonte, fine)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?;

        let required = ancestry_min_snps();
        if result.snps_with_genotype < required {
            return Err(AppError::InsufficientAncestryData { genotyped: result.snps_with_genotype, required });
        }
        ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, &result).await?;
        for extra in [pca_gmm.as_ref(), nmonte.as_ref(), fine.as_ref()].into_iter().flatten() {
            ancestry_result::upsert(self.store.pool(), biosample_guid, CONSENSUS_SOURCE_ID, extra).await?;
        }
        Ok(result)
    }

    pub async fn estimate_ancestry(&self, alignment_id: i64) -> Result<AncestryResult, AppError> {
        self.estimate_ancestry_with_progress(alignment_id, |_, _| {}).await
    }

    /// Like [`estimate_ancestry`], reporting `progress(contigs_done, contigs_total)` as the
    /// per-contig genotyping pass advances — the slow step on a whole-genome BAM (minutes), so
    /// the UI shows a bar. The callback runs on the blocking genotyping thread.
    pub async fn estimate_ancestry_with_progress(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(usize, usize) + Send + 'static,
    ) -> Result<AncestryResult, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);

        // Load the panel for the alignment's build and verify they agree.
        let build = canonical_build(&aln.reference_build)
            .ok_or_else(|| AppError::Refgenome(navigator_refgenome::RefgenomeError::UnknownBuild(aln.reference_build.clone())))?;
        let panel_path = ancestry_panel_path(build);
        let panel_bytes = read_verified_asset(build, &panel_path)?
            .ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        if canonical_build(&panel.build) != Some(build) {
            return Err(AppError::AncestryPanelBuildMismatch {
                panel: panel.build.clone(),
                alignment: aln.reference_build.clone(),
            });
        }

        // Resolve the reference (recorded path, else gateway — needed for CRAM + genotyping).
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        let reference_version = aln.reference_build.clone();
        // Optional PCA loadings (same build): the modern asset projects the sample onto PC space
        // for the scatter; the ancient asset (if present) drives a PCA-projection GMM over ancient
        // components (Steppe/EEF/WHG). The GMM runs against the ancient asset when available, else
        // the modern one — so PCA_PROJECTION_GMM is always over the best available reference.
        // Optional assets: an integrity mismatch on any drops it (logged) rather than failing the run.
        let optional = |path: PathBuf| read_verified_asset(build, &path).unwrap_or_else(|e| { eprintln!("{e}"); None });
        let pca_bytes = optional(ancestry_pca_path(build));
        let ancient_pca_bytes = optional(ancestry_pca_ancient_path(build));
        // Optional fine-population frequencies (same panel sites) → a fine modern admixture over a
        // curated subset. Absent ⇒ silently skipped.
        let fine_bytes = optional(ancestry_freq_global_path(build));

        // Returns (admixture, optional PCA-GMM, optional nMonte, optional fine admixture).
        let (result, pca_gmm, nmonte, fine) = tokio::task::spawn_blocking(move || {
            let params = adaptive_haploid_params(&bam, Some(&reference));
            let genotypes =
                ancestry_analysis::genotype_panel(&bam, Some(&reference), &panel, &params, &mut progress)?;
            // Supervised admixture → 100%-summing composition (the consumer-report shape).
            let mut result =
                ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            // Fine modern admixture over the same genotypes (fine asset shares the AIM sites).
            let fine = fine_bytes
                .and_then(|b| ancestry_analysis::AncestryPanel::from_bytes(&b).ok())
                .map(|fp| ancestry_analysis::estimate_fine_admixture(&genotypes, &fp, &reference_version));
            let modern_pca = pca_bytes.and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok());
            if let Some(pca) = &modern_pca {
                result.pca_coordinates = Some(ancestry_analysis::project_pca(&genotypes, pca));
            }
            // PCA-projection models: prefer the ancient reference asset, else the modern one.
            // Both the GMM (cluster assignment) and the nMonte (distance-minimizing mixture fit)
            // run against the same asset, so a richer/global asset widens both.
            let gmm_pca = ancient_pca_bytes
                .and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok())
                .or(modern_pca);
            let (pca_gmm, nmonte) = match &gmm_pca {
                Some(pca) => (
                    Some(ancestry_analysis::estimate_pca_gmm(&genotypes, pca, &reference_version)),
                    Some(ancestry_analysis::estimate_nmonte(&genotypes, pca, &reference_version)),
                ),
                None => (None, None),
            };
            Ok::<_, navigator_analysis::AnalysisError>((result, pca_gmm, nmonte, fine))
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;

        let required = ancestry_min_snps();
        if result.snps_with_genotype < required {
            return Err(AppError::InsufficientAncestryData {
                genotyped: result.snps_with_genotype,
                required,
            });
        }

        // Persist every estimate (keyed by method) so "publish both" can read each back.
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            ancestry_result::upsert(self.store.pool(), bio, alignment_id, &result).await?;
            for extra in [pca_gmm.as_ref(), nmonte.as_ref(), fine.as_ref()].into_iter().flatten() {
                ancestry_result::upsert(self.store.pool(), bio, alignment_id, extra).await?;
            }
        }
        Ok(result)
    }

    /// Estimate autosomal ancestry from an imported **chip** (23andMe / AncestryDNA), reusing the
    /// same estimators as the alignment path. The chip's GRCh37 SNP genotypes are lifted to the AIMs
    /// panel's CHM13 coordinates (`gateway.lift_positions`, the `hg19-chm13v2` chain) and intersected
    /// with the panel; the per-site alt-allele dosage feeds `estimate_admixture` (+ PCA-GMM / nMonte
    /// when those assets exist). Computed on demand and returned (not persisted — the
    /// `ancestry_result` table is alignment-keyed; chip persistence is a follow-on).
    pub async fn estimate_ancestry_from_chip(&self, chip_profile_id: i64) -> Result<AncestryResult, AppError> {
        let chip = chip_profile::get(self.store.pool(), chip_profile_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("chip profile {chip_profile_id}"))))?;
        let path = chip
            .source_path
            .clone()
            .ok_or_else(|| AppError::Import("this chip has no stored raw-data file — re-import it to enable ancestry".into()))?;
        let text = std::fs::read_to_string(&path).map_err(|e| AppError::Import(format!("chip file {path}: {e}")))?;
        let from_build = chipprofile::detect_build(&text);
        let calls = chipprofile::autosomal_calls(&text);

        // The chip-ancestry target is the CHM13 AIMs panel.
        let build = ReferenceBuild::Chm13v2;
        let panel_path = ancestry_panel_path(build);
        let panel_bytes =
            read_verified_asset(build, &panel_path)?.ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        let reference_version = panel.build.clone();

        // Ensure the GRCh37→CHM13 chain is cached (downloads on first use), then index panel sites.
        self.gateway.resolve_chain(&from_build, &panel.build, &mut |_, _| {}).await?;
        let site_index: HashMap<(&str, i64), &navigator_analysis::ancestry::PanelSite> =
            panel.sites.iter().map(|s| ((s.contig.as_str(), s.position), s)).collect();

        // Group chip calls by contig, lift each contig's positions in one batch, then intersect.
        let mut by_contig: HashMap<String, HashMap<i64, (char, char)>> = HashMap::new();
        for c in &calls {
            by_contig.entry(c.contig.clone()).or_default().insert(c.position, (c.a1, c.a2));
        }
        let mut genotypes: Vec<navigator_analysis::caller::SiteGenotype> = Vec::new();
        for (contig, alleles) in &by_contig {
            let positions: Vec<i64> = alleles.keys().copied().collect();
            let lifted = self.gateway.lift_positions(&from_build, &panel.build, contig, &positions)?;
            for lp in lifted {
                let Some(site) = site_index.get(&(lp.contig.as_str(), lp.pos)) else { continue };
                let Some(&(a1, a2)) = alleles.get(&lp.tree_pos) else { continue };
                // `dosage_from_alleles` rev-comps as a fallback, so an inverted lift resolves too.
                let Some(dosage) =
                    ancestry_analysis::dosage_from_alleles(a1, a2, site.reference_allele, site.alternate_allele)
                else {
                    continue;
                };
                genotypes.push(navigator_analysis::caller::SiteGenotype {
                    name: format!("{}:{}", site.contig, site.position),
                    contig: site.contig.clone(),
                    position: site.position,
                    reference_allele: site.reference_allele.to_string(),
                    alternate_allele: site.alternate_allele.to_string(),
                    ploidy: 2,
                    dosage,
                    gq: 0,
                    depth: 0,
                    ref_depth: 0,
                    alt_depth: 0,
                    pls: Vec::new(),
                    gt: None,
                    allele_depths: None,
                });
            }
        }

        // A chip intersects far fewer AIMs than a WGS BAM, so a lower floor than the WGS minimum.
        const CHIP_ANCESTRY_MIN_SNPS: usize = 100;
        if genotypes.len() < CHIP_ANCESTRY_MIN_SNPS {
            return Err(AppError::InsufficientAncestryData { genotyped: genotypes.len(), required: CHIP_ANCESTRY_MIN_SNPS });
        }

        let mut result = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
        // Project onto PC space for the scatter when the modern PCA asset is present.
        let modern_pca =
            std::fs::read(ancestry_pca_path(build)).ok().and_then(|b| ancestry_analysis::PcaLoadings::from_bytes(&b).ok());
        if let Some(pca) = &modern_pca {
            result.pca_coordinates = Some(ancestry_analysis::project_pca(&genotypes, pca));
        }
        Ok(result)
    }

    /// The persisted ancestry estimate for an alignment, if one has been computed.
    pub async fn ancestry_for_alignment(
        &self,
        alignment_id: i64,
    ) -> Result<Option<AncestryResult>, AppError> {
        Ok(ancestry_result::get_for_alignment(self.store.pool(), alignment_id).await?)
    }

    /// The persisted **fine-population** admixture estimate for an alignment, if one was computed
    /// (the `ancestry_freq_global` asset was present at estimation time). Drives the super→fine
    /// hierarchy rows; the super-pop donut keeps using the primary ([`ancestry_for_alignment`]).
    pub async fn fine_ancestry_for_alignment(
        &self,
        alignment_id: i64,
    ) -> Result<Option<AncestryResult>, AppError> {
        Ok(ancestry_result::get_for_alignment_method(self.store.pool(), alignment_id, "FINE_ADMIXTURE").await?)
    }

    /// Reference population centroids on (PC1, PC2) for the alignment's build — the backdrop
    /// for the PCA scatter. `(population_code, pc1, pc2)`; empty if no PCA loadings are present.
    pub async fn ancestry_pca_reference(
        &self,
        alignment_id: i64,
    ) -> Result<Vec<(String, f64, f64)>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let Some(build) = canonical_build(&aln.reference_build) else { return Ok(Vec::new()) };
        let Ok(bytes) = std::fs::read(ancestry_pca_path(build)) else { return Ok(Vec::new()) };
        let pca = navigator_analysis::ancestry::PcaLoadings::from_bytes(&bytes)?;
        Ok(pca
            .populations
            .iter()
            .enumerate()
            .map(|(p, code)| {
                let c = pca.centroid(p);
                (code.clone(), c.first().copied().unwrap_or(0.0) as f64, c.get(1).copied().unwrap_or(0.0) as f64)
            })
            .collect())
    }

    /// Paint each chromosome with local ancestry from the subject's **consensus** — no BAM walk.
    /// Reads the cached autosomal [`DiploidProfile`], anchors on the genome-wide admixture
    /// composition (estimated from the same consensus dosages), and runs the per-chromosome HMM.
    pub async fn paint_local_ancestry_from_consensus(&self, biosample_guid: SampleGuid) -> Result<Vec<AncestrySegment>, AppError> {
        let profile = self
            .cached_autosomal_profile(biosample_guid)
            .await?
            .ok_or_else(|| AppError::Import("build the autosomal consensus first (Autosomal tab) before painting".into()))?;
        let genotypes = consensus_genotypes(&profile);
        let build = ReferenceBuild::Chm13v2;
        let reference_version = "chm13v2.0".to_string();
        let panel_path = ancestry_panel_path(build);
        let panel_bytes =
            read_verified_asset(build, &panel_path)?.ok_or_else(|| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        let segments = tokio::task::spawn_blocking(move || {
            let composition = ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let prior: Vec<(String, f64)> =
                composition.components.iter().map(|c| (c.population_code.clone(), c.percentage / 100.0)).collect();
            ancestry_analysis::paint_local_ancestry(&genotypes, &panel, &prior, &ancestry_analysis::PaintParams::default())
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?;
        Ok(segments)
    }

    /// Paint each chromosome with local ancestry (the "DNA painting"): genotype the AIMs panel,
    /// anchor on the genome-wide admixture composition, and run the per-chromosome HMM. Returns
    /// the ancestry segments. `progress(contigs_done, total)` reports the genotyping pass.
    pub async fn local_ancestry(&self, alignment_id: i64) -> Result<Vec<AncestrySegment>, AppError> {
        self.local_ancestry_with_progress(alignment_id, |_, _| {}).await
    }

    /// [`local_ancestry`] with a genotyping-progress callback (for the UI bar).
    pub async fn local_ancestry_with_progress(
        &self,
        alignment_id: i64,
        mut progress: impl FnMut(usize, usize) + Send + 'static,
    ) -> Result<Vec<AncestrySegment>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let build = canonical_build(&aln.reference_build)
            .ok_or_else(|| AppError::Refgenome(navigator_refgenome::RefgenomeError::UnknownBuild(aln.reference_build.clone())))?;
        let panel_path = ancestry_panel_path(build);
        let panel_bytes = std::fs::read(&panel_path).map_err(|_| AppError::AncestryPanelMissing(panel_path.clone()))?;
        let panel = AncestryPanel::from_bytes(&panel_bytes)?;
        if canonical_build(&panel.build) != Some(build) {
            return Err(AppError::AncestryPanelBuildMismatch {
                panel: panel.build.clone(),
                alignment: aln.reference_build.clone(),
            });
        }
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        let reference_version = aln.reference_build.clone();

        let segments = tokio::task::spawn_blocking(move || {
            let params = adaptive_haploid_params(&bam, Some(&reference));
            let genotypes =
                ancestry_analysis::genotype_panel(&bam, Some(&reference), &panel, &params, &mut progress)?;
            // Genome-wide composition → the HMM's switch prior.
            let composition =
                ancestry_analysis::estimate_admixture(&genotypes, &panel, &reference_version);
            let prior: Vec<(String, f64)> = composition
                .components
                .iter()
                .map(|c| (c.population_code.clone(), c.percentage / 100.0))
                .collect();
            let segs = ancestry_analysis::paint_local_ancestry(
                &genotypes,
                &panel,
                &prior,
                &ancestry_analysis::PaintParams::default(),
            );
            Ok::<_, navigator_analysis::AnalysisError>(segs)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;

        Ok(segments)
    }

    /// An alignment's content SHA-256, computed once at import. Read from the record if present,
    /// else computed now (hashing the file) and stored — so batch-imported alignments are hashed
    /// lazily on first analysis, then cached on the row.
    async fn alignment_content_hash(&self, alignment_id: i64) -> Result<String, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        if let Some(h) = aln.content_sha256 {
            return Ok(h);
        }
        let bam = aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?;
        let hash = sha256_file_async(PathBuf::from(bam)).await?;
        let _ = alignment::set_content_hash(self.store.pool(), alignment_id, &hash).await;
        Ok(hash)
    }

    /// Fingerprint of the inputs to a Y-haplogroup score: the alignment's content hash + the
    /// active Y-tree's content hash. Unchanged inputs → a re-score is unnecessary. Errors (e.g.
    /// the tree is unreachable and uncached) disable caching for this run rather than failing.
    async fn y_score_fingerprint(&self, alignment_id: i64) -> Result<String, AppError> {
        let file_hash = self.alignment_content_hash(alignment_id).await?;
        let tree_json = match y_tree_provider() {
            YTreeProvider::DecodingUs => self.fetch_decodingus_y_tree().await?,
            YTreeProvider::Ftdna => self.fetch_ftdna_y_tree().await?,
        };
        let tree_hash = sha256_str(&tree_json);
        Ok(format!("f:{}|yt:{}", &file_hash[..16], &tree_hash[..16]))
    }

    /// Assign a Y haplogroup to an alignment: place the sample against the configured Y tree
    /// (DecodingUs by default — our tree, native CHM13 coords, no liftover — falling back to
    /// FTDNA if the AppView is unreachable), call the sample's base at each tree position on
    /// chrY, and rank by Kulczynski. Requires a recorded BAM/CRAM path. Skips re-scoring when
    /// the alignment file and tree are unchanged since the last run (see [`Self::y_score_fingerprint`]).
    pub async fn assign_y_haplogroup(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let bio = self.biosample_of_alignment(alignment_id).await.ok();
        let source_key = format!("aln:{alignment_id}");

        // Input fingerprint = alignment content hash + active Y-tree content hash. If it matches
        // the recorded call's stamp, neither the file nor the tree changed → return the recorded
        // call without re-scoring (the expensive BAM genotyping).
        let fingerprint = self.y_score_fingerprint(alignment_id).await.ok();
        if let (Some(bio), Some(fp)) = (bio, fingerprint.as_deref()) {
            if haplogroup_call::stored_fingerprint(self.store.pool(), bio, DnaType::Y, &source_key).await?.as_deref()
                == Some(fp)
            {
                if let Some(call) = haplogroup_call::get_one(self.store.pool(), bio, DnaType::Y, &source_key).await? {
                    return Ok(assignment_from_call(&call));
                }
            }
        }

        let assignment = self.y_assignment_full(alignment_id).await?;
        if let Some(bio) = bio {
            self.record_call_fp(bio, DnaType::Y, &source_key, format!("aln #{alignment_id} Y"), &assignment, fingerprint.as_deref()).await?;
        }
        Ok(assignment)
    }

    /// Freshly place an alignment against the configured Y tree, returning the **full** assignment
    /// **including per-branch SNP evidence** (the cached [`assign_y_haplogroup`] path returns only
    /// the terminal). Expensive (genotypes chrY tree sites in the BAM) — used by the Y-variant
    /// profile, which the user builds explicitly.
    async fn y_assignment_full(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        match y_tree_provider() {
            YTreeProvider::DecodingUs => match self.assign_y_decodingus(alignment_id).await {
                Ok(a) => Ok(a),
                Err(e) => {
                    // AppView unreachable / build unsupported / parse failure → FTDNA fallback.
                    eprintln!("DecodingUs Y tree unavailable ({e}); falling back to FTDNA");
                    let tree_json = self.fetch_ftdna_y_tree().await?;
                    self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json).await
                }
            },
            YTreeProvider::Ftdna => {
                let tree_json = self.fetch_ftdna_y_tree().await?;
                self.assign_haplogroup_from_alignment(alignment_id, "chrY", &tree_json).await
            }
        }
    }

    /// Place against the DecodingUs Y tree from our AppView, using the alignment's **native**
    /// build coordinates (`hs1` for CHM13, `GRCh38`, `GRCh37`) — queried directly, **no
    /// liftover**. This is the intended architecture (the AppView owns multi-build coordinates;
    /// Navigator stays liftover-free). Today the AppView's `hs1` coords cover the decoding-us
    /// backbone but not the FTDNA-grafted tips, so deep CHM13 placement is limited until the
    /// AppView enriches `hs1` for every variant (lift GRCh38→hs1 at ingest or on the fly).
    async fn assign_y_decodingus(&self, alignment_id: i64) -> Result<HaploAssignment, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
            AppError::Import(format!("no DecodingUs tree coordinates for build {}", aln.reference_build))
        })?;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key).map_err(AppError::Import)?;
        // Native build → no liftover (tree_source_build = None → direct query).
        let calls = self.base_calls(alignment_id, "chrY", &tree, None).await?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// Assign a Y haplogroup from the subject's imported **BISDNA / Y-SNP-panel** calls — no
    /// alignment required. Builds a derived-allele call map from the subject's `Chip`-sourced
    /// variant sets (the panel's positive calls, each `position → derived base`) and scores it
    /// against the Y tree on `build` (the subject's alignment build, else `"hs1"`). Uses the
    /// DecodingUs tree at the native build (FTDNA fallback only on GRCh38, where positions
    /// match), and the chip-robust terminal selection ([`assemble_assignment_robust`]). The
    /// call is recorded as a reconciliation source. Only derived (positive) calls drive the
    /// Kulczynski ranking, so the stored positives-only variant set is sufficient.
    pub async fn assign_y_bisdna(
        &self,
        biosample_guid: SampleGuid,
        build: Option<&str>,
    ) -> Result<HaploAssignment, AppError> {
        // Derived-allele calls from the subject's chip-sourced variant sets (BISDNA positives).
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;

        // Placement build: explicit override, else the build stored on a chip set at import,
        // else (pre-migration sets with no stored build) re-derive from the subject's alignment.
        let build = match build {
            Some(b) => b.to_string(),
            None => match sets.iter().filter(|s| s.source_type == SourceType::Chip).find_map(|s| s.reference_build.clone()) {
                Some(b) => b,
                None => self.bisdna_target_build(biosample_guid).await,
            },
        };

        let mut calls: HashMap<i64, char> = HashMap::new();
        for s in &sets {
            if s.source_type != SourceType::Chip {
                continue;
            }
            for c in &s.calls {
                if !c.contig.eq_ignore_ascii_case("chrY") && !c.contig.eq_ignore_ascii_case("y") {
                    continue;
                }
                if let Some(b) = c.alternate.chars().next() {
                    calls.insert(c.position, b.to_ascii_uppercase());
                }
            }
        }
        if calls.is_empty() {
            return Err(AppError::Import(
                "no Y-SNP panel calls to place — import a BISDNA file for this subject first".into(),
            ));
        }

        let tree = self.chip_y_tree(&build).await?;
        // Chip alleles (BISDNA + consumer arrays) are plus-strand; flip the minority recorded on
        // the tree's opposite strand so they score against the right allele. No-op for BISDNA.
        let calls = strand_reconcile_to_tree(&tree, calls);
        let assignment = assemble_assignment_robust(&tree, &calls);
        self.record_call(biosample_guid, DnaType::Y, "bisdna", "Chip Y-SNP panel".into(), &assignment).await?;
        Ok(assignment)
    }

    /// Fetch + parse the Y haplotree for a chip placement on `build`. DecodingUs is native multi-build
    /// (no liftover); the FTDNA tree is GRCh38-only, so it's a fallback only when the calls are GRCh38.
    /// Shared by the combined [`assign_y_bisdna`](Self::assign_y_bisdna) placement and the per-panel
    /// Y-profile sources, so the tree is fetched once.
    async fn chip_y_tree(&self, build: &str) -> Result<navigator_analysis::haplo::HaploTree, AppError> {
        match self.fetch_decodingus_y_tree().await {
            Ok(json) => navigator_analysis::haplo::parse_decodingus_json(&json, build).map_err(AppError::Import),
            Err(e) if build == "GRCh38" => {
                eprintln!("DecodingUs Y tree unavailable ({e}); falling back to FTDNA (GRCh38)");
                let json = self.fetch_ftdna_y_tree().await?;
                navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)
            }
            Err(e) => Err(e),
        }
    }

    /// Place one chip/BISDNA panel's chrY calls on `tree` (strand-reconciled), without persisting —
    /// for assembling the per-panel sources of the Y-variant profile.
    fn place_chip_panel(tree: &navigator_analysis::haplo::HaploTree, calls: HashMap<i64, char>) -> HaploAssignment {
        let calls = strand_reconcile_to_tree(tree, calls);
        assemble_assignment_robust(tree, &calls)
    }

    /// Place an mtDNA haplogroup from the subject's chip-sourced MT genotype calls (e.g. 23andMe
    /// `MT` rows) against the FTDNA mt tree. Consumer-array MT positions are rCRS coordinates,
    /// which the tree uses directly (no liftover). Reads every `Chip`-source variant set's chrM
    /// calls, reconciles strand, and uses the robust (sparse-chip) terminal selection. Records a
    /// donor call. The counterpart to [`assign_y_bisdna`](Self::assign_y_bisdna) for mtDNA.
    pub async fn assign_mt_chip(&self, biosample_guid: SampleGuid) -> Result<HaploAssignment, AppError> {
        let sets = variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut calls: HashMap<i64, char> = HashMap::new();
        for s in &sets {
            if s.source_type != SourceType::Chip {
                continue;
            }
            for c in &s.calls {
                let mt = c.contig.eq_ignore_ascii_case("chrM")
                    || c.contig.eq_ignore_ascii_case("chrMT")
                    || c.contig.eq_ignore_ascii_case("mt")
                    || c.contig.eq_ignore_ascii_case("m");
                if !mt {
                    continue;
                }
                if let Some(b) = c.alternate.chars().next() {
                    calls.insert(c.position, b.to_ascii_uppercase());
                }
            }
        }
        if calls.is_empty() {
            return Err(AppError::Import(
                "no chip mtDNA calls to place — import a 23andMe raw-data file for this subject first".into(),
            ));
        }
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let calls = strand_reconcile_to_tree(&tree, calls);
        let assignment = assemble_assignment_robust(&tree, &calls);
        self.record_call(biosample_guid, DnaType::Mt, "chip-mt", "Chip mtDNA panel".into(), &assignment).await?;
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
        let (tree, calls) = self.tree_base_calls(alignment_id, contig, tree_json).await?;
        Ok(assemble_assignment(&tree, &calls))
    }

    /// Like [`assign_haplogroup_from_alignment`], but also returns the per-SNP evidence along
    /// the called terminal's lineage (each defining mutation's Derived/Ancestral/NoCall state).
    /// For exact comparisons (e.g. GRCh38 vs a lifted CHM13 call).
    pub async fn assign_haplogroup_detail(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<(HaploAssignment, Vec<SnpEvidence>, HashMap<i64, char>), AppError> {
        let (tree, calls) = self.tree_base_calls(alignment_id, contig, tree_json).await?;
        let assignment = assemble_assignment(&tree, &calls);
        let lineage = match assignment.ranked.first() {
            Some(top) => navigator_analysis::haplo::lineage_evidence(&tree, &calls, top.id),
            None => Vec::new(),
        };
        Ok((assignment, lineage, calls))
    }

    /// Parse the tree, build the per-position base calls (lifting onto the alignment's build
    /// when needed), and return both. Shared by the assignment + detail entry points.
    async fn tree_base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        tree_json: &str,
    ) -> Result<(navigator_analysis::haplo::HaploTree, HashMap<i64, char>), AppError> {
        let tree = navigator_analysis::haplo::parse_ftdna_json(tree_json).map_err(AppError::Import)?;
        // FTDNA tree positions are in the tree's own build (Y → GRCh38, mt → rCRS/direct).
        let source_build = tree_build_for_contig(contig);
        let calls = self.base_calls(alignment_id, contig, &tree, source_build).await?;
        Ok((tree, calls))
    }

    /// Base-call an alignment at a parsed tree's positions on `contig`. `tree_source_build` is
    /// the build the tree's positions are in: when it differs from the alignment build the
    /// positions are lifted (chrY chain), queried there, and mapped back; `None` (e.g. a
    /// DecodingUs tree already in the alignment's build, or mt/rCRS-direct) queries directly.
    async fn base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(alignment_id))?);
        let reference = aln.reference_path.map(PathBuf::from);

        let targets: HashSet<i64> =
            tree.nodes.values().flat_map(|n| n.loci.iter().map(|l| l.position)).collect();

        let lifted = self
            .lifted_targets(&aln.reference_build, reference.as_deref(), contig, &targets, tree_source_build)
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
        Ok(calls)
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
        tree_source_build: Option<&str>,
    ) -> Result<Option<Vec<LiftedPos>>, AppError> {
        if targets.is_empty() {
            return Ok(None);
        }

        // chrY: downloaded nuclear chain (when the tree build differs from the alignment).
        if let Some(src) = tree_source_build {
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
                    // Rotation-aware: CHM13's chrM is a circular permutation of rCRS.
                    navigator_analysis::mtvariants::mt_position_map(navigator_analysis::mtvariants::rcrs(), &chrm)
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

    // ---- fast path: place haplogroups from precomputed pipeline GVCFs ---------

    /// Build a tree's per-position base calls for an alignment from a **precomputed GVCF**
    /// (the fast path — no CRAM pileup). Lifts tree positions onto the GVCF's build when the
    /// tree's coordinates differ (mt rCRS-tree vs CHM13 `chrM`), exactly as the CRAM path does,
    /// then reads the GVCF instead of walking reads.
    async fn gvcf_base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        gvcf: &Path,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        // The reference is required: a GVCF hom-ref site means "the sample's base == the
        // reference base" — and the reference (e.g. CHM13 = HG002/J1 Y) is itself deep in the
        // tree, so its base there is often the *derived* allele, not the ancestral. We read the
        // reference base at every callable tree position (exactly what call_bases_at observes).
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => self.gateway.resolve_reference(&aln.reference_build, &mut |_, _| {}).await?,
        };
        let targets: HashSet<i64> =
            tree.nodes.values().flat_map(|n| n.loci.iter().map(|l| l.position)).collect();
        if targets.is_empty() {
            return Ok(HashMap::new());
        }
        let params = gvcf::GvcfReadParams::default();

        let lifted = self
            .lifted_targets(&aln.reference_build, Some(&reference), contig, &targets, tree_source_build)
            .await?;

        match lifted {
            // Native: tree positions are already in the GVCF's coordinates → direct read, then
            // resolve hom-ref bases from the reference at the same positions.
            None => {
                let gvcf = gvcf.to_path_buf();
                let contig_s = contig.to_string();
                let targets2 = targets.clone();
                let called = tokio::task::spawn_blocking(move || {
                    gvcf::read_called_bases(&gvcf, &contig_s, &targets2, &params)
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))??;
                let ref_base = self.reference_bases(&reference, contig, &called.callable).await?;
                Ok(gvcf::assemble_calls(&called, &ref_base))
            }
            // Lifted: read the GVCF at each lifted contig + the reference bases there, then map
            // observations back to tree positions (reverse-complementing minus-strand lifts).
            Some(lifted) => {
                let mut by_contig: HashMap<String, HashSet<i64>> = HashMap::new();
                for lp in &lifted {
                    by_contig.entry(lp.contig.clone()).or_default().insert(lp.pos);
                }
                let mut all = gvcf::CalledBases::default();
                let mut ref_base: HashMap<i64, char> = HashMap::new();
                for (qcontig, set) in by_contig {
                    let gvcf = gvcf.to_path_buf();
                    let qc = qcontig.clone();
                    let set2 = set.clone();
                    let called = tokio::task::spawn_blocking(move || {
                        gvcf::read_called_bases(&gvcf, &qc, &set2, &params)
                    })
                    .await
                    .map_err(|e| AppError::Join(e.to_string()))??;
                    ref_base.extend(self.reference_bases(&reference, &qcontig, &called.callable).await?);
                    all.variant_bases.extend(called.variant_bases);
                    all.callable.extend(called.callable);
                }
                Ok(assemble_calls_lifted(&all, &lifted, &ref_base))
            }
        }
    }

    /// Reference genome bases (uppercase A/C/G/T) at `positions` on `contig`. Reads the contig
    /// sequence once off-thread; positions are 1-based. Non-ACGT / out-of-range positions are
    /// omitted. Used by the GVCF fast path to resolve hom-ref tree sites to the actual base.
    async fn reference_bases(
        &self,
        reference: &Path,
        contig: &str,
        positions: &HashSet<i64>,
    ) -> Result<HashMap<i64, char>, AppError> {
        if positions.is_empty() {
            return Ok(HashMap::new());
        }
        let reference = reference.to_path_buf();
        let contig = contig.to_string();
        let positions: Vec<i64> = positions.iter().copied().collect();
        let map = tokio::task::spawn_blocking(move || -> Result<HashMap<i64, char>, navigator_analysis::AnalysisError> {
            let seq = navigator_analysis::reader::read_contig_sequence(&reference, &contig)?;
            let mut m = HashMap::with_capacity(positions.len());
            for p in positions {
                if p >= 1 && (p as usize) <= seq.len() {
                    let b = seq[p as usize - 1].to_ascii_uppercase();
                    if matches!(b, b'A' | b'C' | b'G' | b'T') {
                        m.insert(p, b as char);
                    }
                }
            }
            Ok(m)
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))??;
        Ok(map)
    }

    /// Fingerprint of a GVCF-sourced placement: the GVCF's content hash ⊕ the tree's hash.
    /// Distinct from the CRAM-based [`Self::y_score_fingerprint`] (`gv:` vs `f:` prefix) so a
    /// later deep analyze can tell the call came from a sidecar (phase: deep-pass skip logic).
    async fn gvcf_fingerprint(&self, gvcf: &Path, tree_json: &str, tag: &str) -> Result<String, AppError> {
        let h = sha256_file_async(gvcf.to_path_buf()).await?;
        Ok(format!("gv:{}|{}:{}", &h[..16], tag, &sha256_str(tree_json)[..16]))
    }

    /// Assign a Y haplogroup from a precomputed chrY GVCF — no CRAM walk. Places against the
    /// DecodingUs tree at the alignment's native build (liftover-free), records the call under
    /// the same source key as the CRAM path (`aln:{id}`) with a `gv:`-prefixed fingerprint.
    /// Errors if the build has no DecodingUs coordinates or the tree is unreachable; the caller
    /// (`ingest_sidecars`) treats that as "leave Y for the deep pass".
    pub async fn assign_y_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
            AppError::Import(format!("no DecodingUs tree coordinates for build {}", aln.reference_build))
        })?;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key).map_err(AppError::Import)?;
        let calls = self.gvcf_base_calls(alignment_id, "chrY", gvcf, &tree, None).await?;
        // Robust (proportional-top) selection, not the strict alignment-tuned guard. A
        // joint-genotyped GVCF gives confident calls that include a few stray ancestral
        // contradictions on the deep backbone (recurrent sites, the CHM13=J1 reference, joint
        // hard-filters); strict `path_admissible` then vetoes the genuine deep lineage and
        // drops to a shallow node (HG00096 → A1b instead of its true R1b1a1b1a1a, which `score`
        // ranks top at 344/364). This is the same confident-but-sparse-contradiction regime as
        // BISDNA chip data — see [`assemble_assignment_robust`].
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "yt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Y,
                &format!("aln:{alignment_id}"),
                format!("aln #{alignment_id} Y (pipeline GVCF)"),
                &assignment,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Assign an mtDNA haplogroup from a precomputed chrM GVCF — no CRAM walk. Places against
    /// the FTDNA mt tree; on CHM13 the tree's rCRS positions are lifted onto `chrM` (the cheap
    /// self-generated rCRS↔chrM map), on GRCh38 they're read directly. Recorded under the CRAM
    /// path's mt source key (`aln:{id}:mt`) with a `gv:`-prefixed fingerprint.
    pub async fn assign_mt_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let source_build = tree_build_for_contig("chrM"); // None → rCRS-direct / chrM lift
        let calls = self.gvcf_base_calls(alignment_id, "chrM", gvcf, &tree, source_build).await?;
        // Robust selection, as for Y (see assign_y_from_gvcf) — the GVCF's confident calls fit
        // the proportional-top regime better than the strict alignment guard.
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "mt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Mt,
                &format!("aln:{alignment_id}:mt"),
                format!("aln #{alignment_id} mtDNA (pipeline GVCF)"),
                &assignment,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Fast-path ingest of a sample's pipeline sidecars onto one alignment: place Y + mt from
    /// the GVCFs, and fill sex / read-metrics / lite-coverage from the text sidecars — all
    /// without touching the CRAM. Each step is independent and best-effort: a failure is
    /// recorded in the returned report and the rest proceed (a missing/!matching sidecar just
    /// leaves that result for the deep pass). Returns what it managed to fill.
    pub async fn ingest_sidecars(
        &self,
        alignment_id: i64,
        sidecars: &SampleSidecars,
    ) -> Result<SidecarIngest, AppError> {
        let mut out = SidecarIngest::default();

        if let Some(gvcf) = &sidecars.chr_y_gvcf {
            match self.assign_y_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.y_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("Y from GVCF: {e}")),
            }
        }
        if let Some(gvcf) = &sidecars.chr_m_gvcf {
            match self.assign_mt_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.mt_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("mt from GVCF: {e}")),
            }
        }
        if let Some(path) = &sidecars.sex {
            match self.ingest_sex_sidecar(alignment_id, path).await {
                Ok(s) => out.sex = Some(s),
                Err(e) => out.errors.push(format!("sex: {e}")),
            }
        }
        // Read metrics: richest source wins — samtools `stats` (full, with histograms) > Picard
        // AlignmentSummaryMetrics > samtools `flagstat` (counts only).
        match self.ingest_read_metrics(alignment_id, sidecars).await {
            Ok(true) => out.read_metrics = true,
            Ok(false) => {}
            Err(e) => out.errors.push(format!("read metrics: {e}")),
        }
        // Coverage: samtools `coverage` gives per-contig stats; Picard CollectWgsMetrics gives the
        // genome-wide depth distribution (median/sd/MAD, exclusion fractions, pct_Nx). Use whichever
        // are present, overlaying the distribution onto the per-contig breakdown.
        if sidecars.coverage.is_some() || sidecars.wgs_metrics.is_some() {
            match self.ingest_coverage_sidecar(alignment_id, sidecars).await {
                Ok(()) => out.lite_coverage = true,
                Err(e) => out.errors.push(format!("coverage: {e}")),
            }
        }
        Ok(out)
    }

    async fn ingest_sex_sidecar(&self, alignment_id: i64, path: &Path) -> Result<String, AppError> {
        let text = tokio::fs::read_to_string(path).await.map_err(|e| AppError::Import(format!("{}: {e}", path.display())))?;
        let result = sidecar::parse_sex(&text);
        self.save_analysis_with_provenance(alignment_id, "sex", "1", &result, "pipeline-sidecar", "full").await?;
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(match result.inferred_sex {
            InferredSex::Male => "M",
            InferredSex::Female => "F",
            InferredSex::Unknown => "U",
        }
        .to_string())
    }

    /// Ingest read metrics from the best available sidecar (priority: samtools `stats` →
    /// Picard AlignmentSummaryMetrics → samtools `flagstat`). Returns whether one was found.
    async fn ingest_read_metrics(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<bool, AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move { tokio::fs::read_to_string(&p).await.map_err(|e| AppError::Import(format!("{}: {e}", p.display()))) }
        };
        // (metrics, completeness): samtools stats is full (carries histograms); the others are
        // counts/scalars only, so `partial` lets a deep read-metrics walk upgrade them later.
        let (metrics, completeness) = if let Some(p) = &sidecars.stats {
            (sidecar::parse_samtools_stats(&read(p).await?), "full")
        } else if let Some(p) = &sidecars.alignment_summary {
            match sidecar::parse_alignment_summary(&read(p).await?) {
                Some(m) => (m, "partial"),
                None => return Ok(false),
            }
        } else if let Some(p) = &sidecars.flagstat {
            (sidecar::parse_flagstat(&read(p).await?), "partial")
        } else {
            return Ok(false);
        };
        self.save_analysis_with_provenance(alignment_id, "read_metrics", "1", &metrics, "pipeline-sidecar", completeness).await?;
        Ok(true)
    }

    async fn ingest_coverage_sidecar(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<(), AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move { tokio::fs::read_to_string(&p).await.map_err(|e| AppError::Import(format!("{}: {e}", p.display()))) }
        };
        // Per-contig stats + callable counts from samtools coverage (empty base if absent).
        let lite = match &sidecars.coverage {
            Some(cp) => {
                let cov = read(cp).await?;
                let summary = match &sidecars.callable_summary {
                    Some(p) => Some(read(p).await?),
                    None => None,
                };
                sidecar::lite_coverage(&cov, summary.as_deref())
            }
            None => CoverageResult::default(),
        };
        // Overlay Picard's genome-wide depth distribution onto the per-contig breakdown: start from
        // the Picard result (median/sd/MAD, exclusion fractions, pct_Nx) and graft the contig stats.
        let result = match &sidecars.wgs_metrics {
            Some(wp) => match sidecar::parse_wgs_metrics(&read(wp).await?) {
                Some(mut w) => {
                    w.contig_coverage_stats = lite.contig_coverage_stats;
                    w.contig_callable = lite.contig_callable;
                    w.callable_bases = lite.callable_bases;
                    if w.genome_territory == 0 {
                        w.genome_territory = lite.genome_territory;
                    }
                    if w.mean_coverage == 0.0 {
                        w.mean_coverage = lite.mean_coverage;
                    }
                    w
                }
                None => lite,
            },
            None => lite,
        };
        // Still `partial`: no per-base depth histogram (only the deep walk produces that), so the
        // deep pass still upgrades this. Stored under the standard coverage key.
        self.save_analysis_with_provenance(alignment_id, "coverage", coverage::COVERAGE_VERSION, &result, "pipeline-sidecar", "partial")
            .await?;
        Ok(())
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
            // Long, accurate reads (HiFi) are callable from a single read (see adaptive_min_depth).
            params.min_depth = adaptive_min_depth(params.min_depth, read_len);
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
        let bucket = self.private_y_core(alignment_id, Some(mask)).await?;
        // Persist the self-masked bucket so it reloads instead of recomputing next session.
        self.save_analysis(alignment_id, "private_y", "1", &bucket).await?;
        Ok(bucket)
    }

    /// Cached self-masked private-Y bucket for an alignment, if previously computed.
    pub async fn cached_private_y(&self, alignment_id: i64) -> Result<Option<PrivateBucket>, AppError> {
        self.load_analysis(alignment_id, "private_y", "1").await
    }

    /// Shared core: assign Y, de-novo chrY, subtract the backbone, optionally mask, classify.
    /// The curated CHM13 chrY structural regions (palindrome/amplicon/AZF-DYZ), resolving +
    /// caching the three BEDs on first use. Best-effort: any download/parse failure yields
    /// `None` so the annotation never blocks the analysis.
    /// Genome-region metadata (centromere/telomere/cytoband/PAR) for a build, via the gateway's
    /// 2-layer cache (fetches the UCSC cytoBand table on a cold miss). For QC / display context.
    pub async fn genome_regions(&self, build: &str) -> Result<std::sync::Arc<GenomeRegions>, AppError> {
        Ok(self.gateway.genome_regions(build, &mut |_, _| {}).await?)
    }

    /// Region annotation for a 1-based `position` on `contig` in `build` (centromere/telomere/PAR
    /// membership + cytoband name). Uses the cached regions only — `None` if not yet fetched.
    pub fn region_annotation(&self, build: &str, contig: &str, position: i64) -> Option<RegionAnnotation> {
        self.gateway.cached_genome_regions(build).map(|r| r.annotate(contig, position))
    }

    async fn y_structural_regions(&self) -> Option<navigator_analysis::mask::YStructuralRegions> {
        let amplicon = self.gateway.resolve_mask("chm13v2.0Y_amplicons_v1", &mut |_, _| {}).await.ok()?;
        let palindrome = self.gateway.resolve_mask("chm13v2.0Y_inverted_repeats_v1", &mut |_, _| {}).await.ok()?;
        let azf_dyz = self.gateway.resolve_mask("chm13v2.0Y_AZF_DYZ_v1", &mut |_, _| {}).await.ok()?;
        navigator_analysis::mask::YStructuralRegions::from_beds(&amplicon, &palindrome, &azf_dyz).ok()
    }

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

        // The structural BEDs are in CHM13 chrY coordinates, so they only apply to a CHM13
        // alignment (the de-novo positions are in the alignment's build). Best-effort.
        let aln = alignment::get(self.store.pool(), alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {alignment_id}"))))?;
        let regions = match canonical_build(&aln.reference_build) {
            Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs) => self.y_structural_regions().await,
            _ => None,
        };

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
                region: regions.as_ref().and_then(|r| r.classify(c.position)),
            })
            .collect();
        variants.sort_by_key(|v| v.position);
        Ok(PrivateBucket { terminal: terminal.name.clone(), variants })
    }

    // ---- unified import ----------------------------------------------------

    /// Detect a file's type and route it to the right subject importer (STR / variants /
    /// chip / mtDNA), using sensible defaults. Returns the detected type. Alignment files
    /// are rejected here — they attach to a sequencing test, not directly to a subject.
    /// Probe a BAM/CRAM header for the build/aligner/platform/test-type (best-effort).
    pub async fn probe_alignment(&self, path: PathBuf) -> Result<AlignmentProbe, AppError> {
        tokio::task::spawn_blocking(move || navigator_analysis::probe::probe_alignment(&path))
            .await
            .map_err(|e| AppError::Join(e.to_string()))?
            .map_err(AppError::from)
    }

    /// Scan a bounded prefix of an alignment's reads to infer the instrument/library identity —
    /// the `@RG SM/LB/PU` tags plus the most-frequent instrument/flowcell/platform from read names
    /// (the crowd-source input for resolving the lab). Off-thread (blocking IO + CRAM decode);
    /// `reference` is required for CRAM. Best-effort — callers tolerate an error.
    pub async fn library_stats(
        &self,
        path: PathBuf,
        reference: Option<PathBuf>,
    ) -> Result<navigator_analysis::library_stats::LibraryStats, AppError> {
        tokio::task::spawn_blocking(move || {
            navigator_analysis::library_stats::scan_library_stats(
                &path,
                reference.as_deref(),
                navigator_analysis::library_stats::DEFAULT_MAX_READS,
            )
        })
        .await
        .map_err(|e| AppError::Join(e.to_string()))?
        .map_err(AppError::from)
    }

    /// Auto-import an alignment file by probing its header: create the sequencing run (test type,
    /// platform, instrument) and the alignment (reference build + aligner) with no questions
    /// asked. The reference FASTA is **not** required — it's resolved from the build on demand;
    /// if already cached it's stored so every analysis step has it immediately.
    async fn import_alignment_file(&self, biosample_guid: SampleGuid, path: &Path) -> Result<(), AppError> {
        // Idempotent: skip if this exact BAM/CRAM is already recorded as an alignment.
        let path_str = path.to_string_lossy().into_owned();
        if alignment::list_all(self.store.pool())
            .await?
            .iter()
            .any(|a| a.bam_path.as_deref() == Some(path_str.as_str()))
        {
            return Ok(());
        }
        // Best-effort: a probe failure falls back to filename/defaults rather than aborting.
        let probe = self.probe_alignment(path.to_path_buf()).await.unwrap_or_default();

        // Resolve the reference first — the read-name scan needs it to decode a CRAM.
        let reference_build = probe.reference_build.clone().unwrap_or_else(|| reference_build_for(path));
        // Store the cached reference path if we have it; otherwise leave it unset (resolved on
        // demand) — never block import on a download.
        let reference_path = self
            .gateway
            .cached_reference(&reference_build)
            .map(|p| p.to_string_lossy().into_owned());

        // Read-name scan → instrument/library identity (the lab crowd-source input). Best-effort:
        // it fills the platform/model the header `@RG` left blank, and the instrument/flowcell that
        // never live in the header. Skipped silently if the file can't be read (e.g. CRAM with no
        // resolved reference yet).
        let stats = self
            .library_stats(path.to_path_buf(), reference_path.as_deref().map(PathBuf::from))
            .await
            .ok();

        // Platform/model: prefer the header `@RG` (PL/PM); fall back to the read-name inference.
        let platform_name = probe
            .platform
            .clone()
            .or_else(|| stats.as_ref().and_then(|s| s.platform.clone()).map(|p| p.to_uppercase()))
            .unwrap_or_else(|| "UNKNOWN".into());
        let instrument_model = probe
            .instrument_model
            .clone()
            .or_else(|| stats.as_ref().and_then(|s| s.instrument_model.clone()));

        // Test type: refine the header/platform guess with coverage *shape* from the BAI index —
        // a targeted-Y pile-up (autosomes empty) → Big Y / Y Elite / YSEQ; an mtDNA pile-up →
        // mtFull. Best-effort and cheap (O(contigs), no read scan); CRAM / unindexed BAMs have no
        // profile and keep the platform-based guess.
        let test_type = {
            let p = path.to_path_buf();
            let profile = tokio::task::spawn_blocking(move || {
                navigator_analysis::testtype::coverage_profile_from_bai(&p, None)
            })
            .await
            .ok()
            .flatten();
            navigator_analysis::testtype::infer_test_type(
                profile.as_ref(),
                probe.platform.as_deref(),
                probe.vendor_hint.as_deref(),
                None,
            )
            .or_else(|| probe.test_type.clone())
            .unwrap_or_else(|| "WGS".into())
        };

        let run = self
            .record_sequence_run(NewSequenceRun {
                biosample_guid,
                platform_name,
                instrument_model,
                test_type,
                library_layout: stats.as_ref().and_then(|s| s.library_layout.clone()),
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await?;

        // Persist the inferred lab/instrument identity block (the crowd-source key). The lab
        // (`sequencing_facility`) stays unset — set manually, or resolved from `instrument_id`
        // once the AppView lookup ships (roadmap D8).
        if let Some(s) = &stats {
            let _ = sequence_run::set_library_stats(
                self.store.pool(),
                run.id,
                s.instrument_id.as_deref(),
                s.sample_name.as_deref(),
                s.library_id.as_deref(),
                s.platform_unit.as_deref(),
                s.flowcell_id.as_deref(),
            )
            .await;
            // Resolve the lab from the instrument id via the AppView (best-effort, cached).
            if let Some(inst) = s.instrument_id.as_deref() {
                if let Some(lab) = self.lookup_lab_by_instrument(inst).await {
                    let _ = sequence_run::set_facility(self.store.pool(), run.id, &lab).await;
                }
            }
        }

        // Defer the content hash (the file's identity, used to invalidate cached analyses): a
        // whole-file SHA-256 of a multi-GB alignment would block this import for minutes with no
        // feedback. Like the batch path, leave it `None` — `alignment_content_hash` computes and
        // caches it lazily on the first analysis that needs it.
        self.record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build,
            aligner: probe.aligner.clone().unwrap_or_else(|| "unknown".into()),
            variant_caller: None,
            bam_path: Some(path.to_string_lossy().into_owned()),
            reference_path,
            content_sha256: None,
        })
        .await?;
        Ok(())
    }

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
            DetectedData::YSnpPanel => {
                // Build resolved from the subject's alignment, else "hs1" (project default).
                self.import_bisdna_from_file(biosample_guid, path, None).await?;
            }
            DetectedData::ChipData => {
                self.import_chip_profile_from_csv(biosample_guid, None, None, path).await?;
            }
            DetectedData::MtdnaFasta => {
                self.import_mtdna_from_fasta(biosample_guid, path).await?;
            }
            DetectedData::Alignment => {
                self.import_alignment_file(biosample_guid, path).await?;
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
        fast_path: bool,
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
            fast_path: FastPathSummary::default(),
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
                    // Batch import: hash lazily on first analysis (don't stall a bulk NAS import
                    // hashing every multi-GB file up front).
                    content_sha256: None,
                })
                .await?;
                summary.alignments_created += 1;
            }

            // Fast path: ingest the pipeline sidecars onto the build-matching alignment —
            // places Y + mt from the GVCFs and fills sex/metrics/lite-coverage from the text
            // sidecars, no CRAM walk. Best-effort; a failure is tallied and import continues.
            if fast_path && sample.sidecars.has_haplogroup_gvcf() {
                let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
                let chosen = sample
                    .sidecars
                    .build_hint
                    .as_deref()
                    .and_then(|hint| alns.iter().find(|a| build_hint_matches(&a.reference_build, hint)))
                    .or_else(|| alns.iter().find(|a| a.bam_path.is_some()))
                    .or_else(|| alns.first());
                if let Some(a) = chosen {
                    summary.fast_path.samples_with_sidecars += 1;
                    match self.ingest_sidecars(a.id, &sample.sidecars).await {
                        Ok(ing) => {
                            summary.fast_path.y_placed += ing.y_haplogroup.is_some() as usize;
                            summary.fast_path.mt_placed += ing.mt_haplogroup.is_some() as usize;
                            summary.fast_path.sex_filled += ing.sex.is_some() as usize;
                            summary.fast_path.metrics_filled += ing.read_metrics as usize;
                            summary.fast_path.coverage_filled += ing.lite_coverage as usize;
                            for e in ing.errors {
                                summary.fast_path.errors.push(format!("{}: {e}", sample.sample_id));
                            }
                        }
                        Err(e) => summary.fast_path.errors.push(format!("{}: {e}", sample.sample_id)),
                    }
                }
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
            caller::genotype_sites_all_contigs(&bam_pb, &sites, ploidy, &params, reference.as_deref())
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

    /// Resolve an imported chip's genotypes to canonical CHM13 **IBD-panel** dosages — the chip→IBD
    /// path (no alignment, no runtime liftover: the multi-build panel pre-computes coordinates). The
    /// output [`SiteGenotype`]s are over the same CHM13 sites a WGS caller would hit, so a chip and a
    /// WGS sample compare uniformly. Errors if the IBD panel asset isn't built yet.
    pub async fn chip_ibd_dosages(&self, chip_profile_id: i64) -> Result<Vec<SiteGenotype>, AppError> {
        let chip = chip_profile::get(self.store.pool(), chip_profile_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("chip profile {chip_profile_id}"))))?;
        let path = chip
            .source_path
            .clone()
            .ok_or_else(|| AppError::Import("this chip has no stored raw-data file — re-import it to enable IBD".into()))?;
        let text = std::fs::read_to_string(&path).map_err(|e| AppError::Import(format!("chip file {path}: {e}")))?;
        let from_build = chipprofile::detect_build(&text);
        let calls = chipprofile::autosomal_calls(&text);

        let panel_path = ibd_panel_path(ReferenceBuild::Chm13v2);
        let bytes = read_verified_asset(ReferenceBuild::Chm13v2, &panel_path)?.ok_or_else(|| {
            AppError::Import(format!("IBD panel asset not found at {} — build it with `panelbuild ibd-panel`", panel_path.display()))
        })?;
        let panel = navigator_analysis::ibd_panel::IbdPanel::from_bytes(&bytes)?;

        let tuples: Vec<(String, i64, char, char)> =
            calls.into_iter().map(|c| (c.contig, c.position, c.a1, c.a2)).collect();
        Ok(panel.resolve_chip(&from_build, &tuples))
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

        let build = alignment::get(self.store.pool(), alignment_a)
            .await?
            .and_then(|a| canonical_build(&a.reference_build))
            .unwrap_or(ReferenceBuild::Chm13v2);
        Ok(detect_ibd(&ga, &gb, build, config))
    }

    /// IBD comparison over the **chip-compatible IBD panel** for two samples that may each be a
    /// WGS alignment *or* an imported chip (the volume case). Each source resolves to dosages over
    /// the canonical CHM13 IBD-panel sites ([`Self::ibd_panel_dosages`]); the comparison is then
    /// data-type-agnostic. Requires the IBD panel asset (for the WGS-genotyping / chip-resolve path).
    pub async fn compare_ibd_sources(
        &self,
        a: IbdSource,
        b: IbdSource,
        config: IbdDetectorConfig,
    ) -> Result<IbdComparison, AppError> {
        let ga = self.ibd_panel_dosages(a).await?;
        let gb = self.ibd_panel_dosages(b).await?;
        // The IBD panel is CHM13-coordinate, so the CHM13 genetic map applies to both sources.
        Ok(detect_ibd(&ga, &gb, ReferenceBuild::Chm13v2, config))
    }

    /// Dosages over the canonical CHM13 IBD-panel sites for a comparison source. A chip resolves
    /// directly ([`Self::chip_ibd_dosages`]); an alignment genotypes the panel's CHM13 sites from
    /// its BAM (cached per alignment, ploidy-2 autosomal).
    pub async fn ibd_panel_dosages(&self, source: IbdSource) -> Result<Vec<SiteGenotype>, AppError> {
        match source {
            IbdSource::Chip(id) => self.chip_ibd_dosages(id).await,
            IbdSource::Alignment(id) => {
                if let Some(g) = self.load_analysis(id, IBD_PANEL_KIND, caller::GENOTYPE_VERSION).await? {
                    return Ok(g);
                }
                let aln = alignment::get(self.store.pool(), id)
                    .await?
                    .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("alignment {id}"))))?;
                let bam = PathBuf::from(aln.bam_path.ok_or(AppError::MissingPaths(id))?);
                let reference = aln.reference_path.map(PathBuf::from);
                let panel_path = ibd_panel_path(ReferenceBuild::Chm13v2);
                let bytes = read_verified_asset(ReferenceBuild::Chm13v2, &panel_path)?.ok_or_else(|| {
                    AppError::Import(format!("IBD panel asset not found at {} — build it with `panelbuild ibd-panel`", panel_path.display()))
                })?;
                let panel = navigator_analysis::ibd_panel::IbdPanel::from_bytes(&bytes)?;
                let sites: Vec<Site> = panel
                    .sites
                    .iter()
                    .map(|s| Site {
                        name: s.rsid.clone(),
                        contig: s.chm13.contig.clone(),
                        position: s.chm13.position,
                        reference_allele: s.chm13.reference.to_string(),
                        alternate_allele: s.chm13.alternate.to_string(),
                    })
                    .collect();
                let genotypes = tokio::task::spawn_blocking(move || {
                    let params = HaploidCallerParams::default();
                    let contigs: std::collections::BTreeSet<&str> = sites.iter().map(|s| s.contig.as_str()).collect();
                    let mut all = Vec::new();
                    for contig in contigs {
                        all.extend(caller::genotype_sites(&bam, contig, &sites, 2, &params, reference.as_deref())?);
                    }
                    Ok::<_, navigator_analysis::AnalysisError>(all)
                })
                .await
                .map_err(|e| AppError::Join(e.to_string()))??;
                self.save_analysis(id, IBD_PANEL_KIND, caller::GENOTYPE_VERSION, &genotypes).await?;
                Ok(genotypes)
            }
        }
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
        let mut runs = sequence_run::list_for_biosample(self.store.pool(), biosample_guid).await?;
        // One-time backfill: runs analyzed before read stats were mirrored onto the run carry no
        // `total_reads` (and older imports no `library_layout`). Recover them from a cached
        // `read_metrics` artifact on any of the run's alignments and persist, so the card shows
        // library stats + PE/SE without a re-walk.
        for run in &mut runs {
            if run.total_reads.is_some() && run.library_layout.is_some() {
                continue;
            }
            let alns = alignment::list_for_run(self.store.pool(), run.id).await?;
            for a in &alns {
                if let Some(m) = self.cached_read_metrics(a.id).await? {
                    self.write_back_read_stats(a.id, &m).await?;
                    run.total_reads = Some(m.total_reads as i64);
                    run.mean_read_length = (m.mean_read_length > 0.0).then_some(m.mean_read_length);
                    run.mean_insert_size = (m.mean_insert_size > 0.0).then_some(m.mean_insert_size);
                    if m.pf_reads_aligned > 0 {
                        run.library_layout =
                            Some(if m.reads_aligned_in_pairs > 0 { "PAIRED" } else { "SINGLE" }.into());
                    }
                    break;
                }
            }
        }
        Ok(runs)
    }

    /// Cached coverage for several alignments at once (Data Sources alignment rows). `None` for any
    /// alignment without a persisted coverage artifact. No genotyping/walking — pure cache reads.
    pub async fn cached_coverage_bulk(
        &self,
        alignment_ids: &[i64],
    ) -> Result<Vec<(i64, Option<CoverageResult>)>, AppError> {
        let mut out = Vec::with_capacity(alignment_ids.len());
        for &id in alignment_ids {
            out.push((id, self.cached_coverage(id).await?));
        }
        Ok(out)
    }

    /// Alignments for a sequence run.
    /// The best alignment to drive a subject's analysis tabs (subject-centric default): the
    /// highest mean-coverage alignment with a cached coverage result, else the first with a BAM,
    /// else the first. Returns `(sequence_run_id, alignment_id)` so the UI can select the run then
    /// the alignment without the user navigating Data Sources.
    pub async fn default_alignment_for_subject(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Option<(i64, i64)>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        if alignments.is_empty() {
            return Ok(None);
        }
        let mut best: Option<(f64, &Alignment)> = None;
        for a in &alignments {
            if let Some(c) = self.cached_coverage(a.id).await? {
                if best.as_ref().map_or(true, |(cov, _)| c.mean_coverage > *cov) {
                    best = Some((c.mean_coverage, a));
                }
            }
        }
        let chosen = best
            .map(|(_, a)| a)
            .or_else(|| alignments.iter().find(|a| a.bam_path.is_some()))
            .or_else(|| alignments.first());
        Ok(chosen.map(|a| (a.sequence_run_id, a.id)))
    }

    /// Donor-level ancestry: the **consensus** estimate ([`CONSENSUS_SOURCE_ID`]) when present —
    /// it pools all sources, so it's authoritative — else the best-quality per-alignment estimate
    /// (most genotyped SNPs) for back-compat with results predating the consensus path.
    pub async fn donor_ancestry(&self, biosample_guid: SampleGuid) -> Result<Option<(i64, AncestryResult)>, AppError> {
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        if let Some(c) = all.iter().find(|(id, _)| *id == CONSENSUS_SOURCE_ID) {
            return Ok(Some(c.clone()));
        }
        Ok(all.into_iter().max_by_key(|(_, r)| r.snps_with_genotype))
    }

    /// Donor-level private-Y: the **union** of cached (self-masked) private-Y calls across all of
    /// the subject's alignments, deduped by position (keeping the deepest observation). The
    /// terminal is taken from the deepest-covered source bucket.
    pub async fn donor_private_y(&self, biosample_guid: SampleGuid) -> Result<Option<PrivateBucket>, AppError> {
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;
        let mut by_pos: std::collections::HashMap<i64, PrivateVariant> = std::collections::HashMap::new();
        let mut terminal: Option<String> = None;
        let mut any = false;
        for a in &alignments {
            let Some(bucket) = self.cached_private_y(a.id).await? else { continue };
            any = true;
            terminal.get_or_insert_with(|| bucket.terminal.clone());
            for v in bucket.variants {
                by_pos
                    .entry(v.position)
                    .and_modify(|cur| {
                        if v.depth > cur.depth {
                            *cur = v.clone();
                        }
                    })
                    .or_insert(v);
            }
        }
        if !any {
            return Ok(None);
        }
        let mut variants: Vec<PrivateVariant> = by_pos.into_values().collect();
        variants.sort_by_key(|v| v.position);
        Ok(Some(PrivateBucket { terminal: terminal.unwrap_or_default(), variants }))
    }

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
            // A lite (sidecar) coverage is flagged so the UI can badge it and offer a deep walk.
            let coverage_partial = match coverage_aln {
                Some(id) => matches!(
                    self.analysis_provenance(id, "coverage", coverage::COVERAGE_VERSION).await?,
                    Some((_, ref c)) if c == "partial"
                ),
                None => false,
            };
            // Prefer the coverage-bearing alignment; else fall back to the first.
            let primary_alignment_id = coverage_aln.or_else(|| alignments.first().map(|a| a.id));
            let y_haplogroup = self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.map(|c| c.haplogroup);
            let mt_haplogroup = self.haplogroup_consensus(biosample.guid, DnaType::Mt).await?.map(|c| c.haplogroup);
            // Sex + read-metrics from whichever alignment has them cached.
            let mut sex = None;
            let mut metrics = None;
            let mut sv_count = None;
            for a in &alignments {
                if sex.is_none() {
                    sex = self.cached_sex(a.id).await?;
                }
                if metrics.is_none() {
                    metrics = self.cached_read_metrics(a.id).await?;
                }
                if sv_count.is_none() {
                    sv_count = self.cached_sv(a.id).await?.map(|s| s.sv_calls.len());
                }
            }
            let sex = sex.map(|s| match s.inferred_sex {
                navigator_analysis::sex::InferredSex::Male => "M".to_string(),
                navigator_analysis::sex::InferredSex::Female => "F".to_string(),
                navigator_analysis::sex::InferredSex::Unknown => "U".to_string(),
            });
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
                sex,
                mean_read_length: metrics.as_ref().map(|m| m.mean_read_length),
                pct_aligned: metrics.as_ref().map(|m| m.pct_pf_reads_aligned),
                median_insert_size: metrics.as_ref().map(|m| m.median_insert_size),
                sv_count,
                coverage_partial,
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
        let mut summary = AnalyzeSummary {
            project_id,
            samples: 0,
            coverage_done: 0,
            y_done: 0,
            sex_done: 0,
            metrics_done: 0,
            sv_done: 0,
            errors: Vec::new(),
        };
        for biosample in biosample::list_for_project(self.store.pool(), project_id).await? {
            let o = self.analyze_biosample(&biosample).await?;
            if !o.had_alignment {
                continue;
            }
            summary.samples += 1;
            summary.coverage_done += o.coverage_done as usize;
            summary.y_done += o.y_done as usize;
            summary.sex_done += o.sex_done as usize;
            summary.metrics_done += o.metrics_done as usize;
            summary.sv_done += o.sv_done as usize;
            summary.errors.extend(o.errors);
        }
        Ok(summary)
    }

    /// Deep-analyze one biosample's primary (first BAM-bearing) alignment: coverage, Y
    /// haplogroup, sex, read metrics, and SV (≥10× only). Idempotent — a *full* coverage and a
    /// recorded Y/sex/metrics/SV are skipped; a `partial` (lite sidecar) coverage is upgraded by
    /// the per-base walk, which overwrites it. Best-effort: a per-step failure is recorded in
    /// `errors` (prefixed with the donor id) and the remaining steps still run. This is the
    /// per-sample unit the project pass and the streaming deep-analyze job both drive.
    pub async fn analyze_biosample(&self, biosample: &Biosample) -> Result<SampleAnalyzeOutcome, AppError> {
        let mut o = SampleAnalyzeOutcome::default();
        let alignments = alignment::list_for_biosample(self.store.pool(), biosample.guid).await?;
        let Some(aln) = alignments.iter().find(|a| a.bam_path.is_some()) else {
            return Ok(o); // had_alignment stays false
        };
        o.had_alignment = true;
        let label = &biosample.donor_identifier;

        let coverage_full = matches!(
            self.analysis_provenance(aln.id, "coverage", coverage::COVERAGE_VERSION).await?,
            Some((_, ref c)) if c == "full"
        );
        if coverage_full {
            o.coverage_done = true;
        } else {
            match self.run_coverage_for_alignment(aln.id).await {
                Ok(_) => o.coverage_done = true,
                Err(e) => o.errors.push(format!("{label} coverage: {e}")),
            }
        }

        if self.haplogroup_consensus(biosample.guid, DnaType::Y).await?.is_some() {
            o.y_done = true;
        } else {
            match self.assign_y_haplogroup(aln.id).await {
                Ok(_) => o.y_done = true,
                Err(e) => o.errors.push(format!("{label} Y: {e}")),
            }
        }

        if self.cached_sex(aln.id).await?.is_some() {
            o.sex_done = true;
        } else {
            match self.run_sex(aln.id).await {
                Ok(_) => o.sex_done = true,
                Err(e) => o.errors.push(format!("{label} sex: {e}")),
            }
        }

        if self.cached_read_metrics(aln.id).await?.is_some() {
            o.metrics_done = true;
        } else {
            match self.run_read_metrics(aln.id).await {
                Ok(_) => o.metrics_done = true,
                Err(e) => o.errors.push(format!("{label} metrics: {e}")),
            }
        }

        // SV needs ≥10× — only attempt when coverage clears the threshold (avoids logging a
        // "coverage too low" error for every low-coverage sample).
        if self.cached_sv(aln.id).await?.is_some() {
            o.sv_done = true;
        } else if self.cached_coverage(aln.id).await?.map(|c| c.mean_coverage >= 10.0).unwrap_or(false) {
            match self.run_sv(aln.id).await {
                Ok(_) => o.sv_done = true,
                Err(e) => o.errors.push(format!("{label} SV: {e}")),
            }
        }
        Ok(o)
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
    use super::{assemble_assignment, assemble_assignment_robust, strand_reconcile_to_tree, support_backoff_terminal};
    use navigator_analysis::haplo::parse_ftdna_json;
    use std::collections::HashMap;

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

    /// The parsimony back-off trims a net-contradicted deep tail (the sparse-panel / aDNA
    /// over-deepening) to the node where running (derived − ancestral) support peaks, but keeps a
    /// clean deep path and tolerates a lone contradiction outweighed by deeper support.
    #[test]
    fn support_backoff_trims_net_negative_tail_but_keeps_supported_depth() {
        let tree = parse_ftdna_json(SPINE6).unwrap();
        // Derived A+B (peak at B), then below B: ancestral@750, contradiction@1000 (G≠der A),
        // a lone derived@1100 — tail net −1. Should back off F(6) → B(3).
        let sparse: HashMap<i64, char> =
            [(146, 'G'), (263, 'G'), (750, 'C'), (1000, 'G'), (1100, 'T')].into_iter().collect();
        assert_eq!(support_backoff_terminal(&tree, &sparse, 6), 3, "net-negative tail trimmed to B");

        // A clean fully-derived path keeps the deepest terminal F.
        let clean: HashMap<i64, char> =
            [(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A'), (1100, 'T')].into_iter().collect();
        assert_eq!(support_backoff_terminal(&tree, &clean, 6), 6, "clean path keeps the terminal");

        // A lone contradiction (@750) outweighed by deeper derived calls still reaches F.
        let recovered: HashMap<i64, char> =
            [(146, 'G'), (263, 'G'), (750, 'C'), (1000, 'A'), (1100, 'T')].into_iter().collect();
        assert_eq!(support_backoff_terminal(&tree, &recovered, 6), 6, "deeper support recovers depth");
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
            assemble_assignment_robust(&tree, &strand_reconcile_to_tree(&tree, [(146, 'C'), (263, 'G')].into_iter().collect()))
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
        assert_eq!(assemble_assignment_robust(&tree, &calls).ranked.first().unwrap().name, "D");
    }

    /// The GVCF fast path reconstructs exactly the `calls` a pileup would yield. A fully
    /// derived path (every defining SNP a variant) places to the deep terminal D.
    #[test]
    fn gvcf_derived_path_places_deep() {
        use navigator_analysis::gvcf;
        let tree = parse_ftdna_json(TREE).unwrap();
        let mut called = gvcf::CalledBases::default();
        called.variant_bases.extend([(146, 'G'), (263, 'G'), (750, 'T'), (1000, 'A')]);
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
        assert_eq!(calls.get(&750), Some(&'T'), "hom-ref site takes the reference base (derived here)");
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
            LiftedPos { tree_pos: 146, contig: "chrM".into(), pos: 500, reverse: false },
            LiftedPos { tree_pos: 263, contig: "chrM".into(), pos: 900, reverse: true },
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
        sequence_run::set_library_stats(app.store.pool(), run.id, Some("A00182"), None, None, None, Some("H5WLTDMXX"))
            .await
            .unwrap();
        let reloaded = sequence_run::get(app.store.pool(), run.id).await.unwrap().unwrap();

        let value = app.sequence_run_record(&reloaded).await.unwrap();
        assert_eq!(value.get("instrumentId").and_then(|v| v.as_str()), Some("A00182"));
        assert_eq!(value.get("$type").and_then(|v| v.as_str()), Some(NS_SEQUENCERUN));
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
        assert_eq!(out[0].signals, vec!["POPULATION_OVERLAP".to_string(), "HAPLOGROUP".to_string()]);
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
        for no in ["chr23", "chrUn_KI270302v1", "chr1_KI270706v1_random", "HLA-A", "GL000220.1", ""] {
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
        let b = vec![sg("chr1", 100, 2), sg("chr1", 200, 1), sg("chr1", 300, 0), sg("chr1", 400, 0)];
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
        assert!(matches!(resolve_y_provider(Some("ftdna"), Some("decodingus")), YTreeProvider::Ftdna));
        assert!(matches!(resolve_y_provider(Some("decodingus"), Some("ftdna")), YTreeProvider::DecodingUs));
        // settings used when env absent
        assert!(matches!(resolve_y_provider(None, Some("ftdna")), YTreeProvider::Ftdna));
        // default when neither
        assert!(matches!(resolve_y_provider(None, None), YTreeProvider::DecodingUs));
        // unrecognized value falls back to default
        assert!(matches!(resolve_y_provider(Some("bogus"), None), YTreeProvider::DecodingUs));
    }

    #[test]
    fn appview_url_precedence_and_normalization() {
        assert_eq!(resolve_appview_url(Some("https://av.example/".into()), Some("http://x".into())), "https://av.example");
        assert_eq!(resolve_appview_url(None, Some("http://host:9000".into())), "http://host:9000");
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
        };
        let json = serde_json::to_string(&s).unwrap();
        assert_eq!(serde_json::from_str::<AppSettings>(&json).unwrap(), s);

        // Missing/partial fields default to None (forward/backward compatible).
        let partial: AppSettings = serde_json::from_str(r#"{"appview_url":"http://h"}"#).unwrap();
        assert_eq!(partial.appview_url.as_deref(), Some("http://h"));
        assert_eq!(partial.y_tree_provider, None);
        assert_eq!(AppSettings::default(), serde_json::from_str::<AppSettings>("{}").unwrap());
    }
}

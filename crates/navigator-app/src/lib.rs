//! The application layer of Navigator, and its command layer. This crate is the one API that the UI
//! calls. It takes the place of the `WorkbenchViewModel` type, which held too much.
//!
//! The crate controls `navigator-store`, and later the analysis code and the sync code, behind its
//! commands and queries. It also holds each policy that an old dialog held. Those policies assign an
//! identity, test that a record exists, and read and write a result.
//!
//! The UI holds the state of its views and sends commands. No widget calls the database, and no
//! widget makes a decision about the domain.

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

// Re-export the analysis result types the command API returns, so the UI depends only
// on navigator-app (ui -> app), not directly on navigator-analysis.
pub use navigator_analysis::caller::VariantCall as DenovoCall;
pub use navigator_analysis::coverage::CoverageResult as Coverage;
pub use navigator_analysis::haplo::{BranchEvidence, CallState, Locus, NodeEvidence, ScoredHaplogroup, SnpEvidence};
pub use navigator_analysis::heteroplasmy::HeteroplasmySite;
pub use navigator_analysis::mask::YRegionClass;
pub use navigator_analysis::mtvariants::{MtRegion, MtVariant, MtVariantKind};
pub use navigator_analysis::preflight::{
    Check as PreflightCheck, Report as PreflightReport, Status as PreflightStatus,
};
pub use navigator_analysis::CancelToken;

/// Report the state of a BAM **path** or CRAM **path** that has no record in the workspace.
///
/// This case matters when a user reports a file that the app refuses to read. The team needs the
/// answer before it decides whether an import of that file is possible.
///
/// The function blocks, so call it away from the async runtime.
pub fn diagnose_alignment_file(alignment: &std::path::Path, reference: Option<&std::path::Path>) -> PreflightReport {
    navigator_analysis::preflight::diagnose(alignment, reference)
}
pub use navigator_analysis::archaic::{
    ArchaicCallable, ArchaicClassify, ArchaicCountDistribution, ArchaicMarkerPanel, ArchaicMarkerResult,
    ArchaicOutgroup, DiagnosticClass,
};
pub use navigator_analysis::archaic_segments::{
    ArchaicConfig, ArchaicSegment, ArchaicSegmentResult, ArchaicSource, ArchaicSummary,
};
pub use navigator_analysis::probe::AlignmentProbe;
pub use navigator_analysis::read_metrics::{PairOrientation, ReadMetrics};
pub use navigator_analysis::roh::{RohConfig, RohPattern, RohResult, RohSegment, RohSummary};
pub use navigator_analysis::sex::{Confidence as SexConfidence, InferredSex, SexInferenceResult};
pub use navigator_analysis::sv::types::{SvAnalysisResult, SvCall, SvType};
pub use navigator_analysis::unified::UnifiedMetricsResult;
pub use navigator_domain::ancestry::{
    side_label_default, AncestryResult, AncestrySegment, ConfidenceInterval, PaintingResult, PopulationComponent,
    SuperPopulationSummary,
};
// The format of the ancestry panel. This crate exports it again, so a panel tool and a test depend
// on navigator-app alone.
pub use navigator_analysis::ancestry::{AncestryPanel, PanelSite as AncestryPanelSite};

/// One haplogroup assignment. It holds the candidates in their order. For the terminal node that it
/// reports, it also holds each child branch with the evidence of each SNP.
///
/// That evidence shows the reason that the descent stopped. A split with no support shows an
/// ancestral SNP. A split with no answer shows a no-call.
#[derive(Debug, Clone)]
pub struct HaploAssignment {
    pub ranked: Vec<ScoredHaplogroup>,
    pub branches: Vec<BranchEvidence>,
    /// The evidence of each SNP along the lineage, from the root to the terminal node. The list
    /// holds each mutation that defines a node, and the state of the sample at that mutation. A
    /// state has one of three values: Derived, Ancestral, or NoCall.
    ///
    /// The variant **profile**, which pools many sources, reconciles this set.
    ///
    /// This list is not `branches`. That field holds the child branches that the descent did not
    /// take. It shows the reason that the descent stopped, so most of its states are Ancestral or
    /// NoCall.
    pub lineage: Vec<SnpEvidence>,
}

/// A descent report for one lineage, which is the Y lineage or the mtDNA lineage. The report has the
/// shape of a YFull YReport.
///
/// It holds the path from the root to the terminal node. Each node holds the SNPs that define it,
/// and the call state of the subject at each SNP.
///
/// The type takes a [`DnaType`] parameter, so the Y-DNA tab and the mtDNA tab share one model and
/// one renderer. [`App::descent_report`] builds it.
#[derive(Debug, Clone)]
pub struct DescentReport {
    pub dna: DnaType,
    /// The reported terminal haplogroup name (e.g. "R-FGC29071", "U5a1b1g").
    pub terminal: String,
    /// The nodes from the root to the terminal node. Each node holds the SNPs that define it and
    /// the state of the sample, in a `NodeEvidence` value.
    pub nodes: Vec<NodeEvidence>,
}

/// The **block tree** of the cohort of one project.
///
/// The tree is the part of the haplotree that covers the terminal haplogroups of the members. Each
/// node is a *block* of SNPs that the tree treats as equivalent. Each member appears below its own
/// terminal node.
///
/// This type is the group-project form of [`DescentReport`]. That report draws the path of one
/// subject from the root to its terminal node. This tree draws the place of each member of a cohort
/// against the other members.
///
/// [`App::project_block_tree`] builds it. See `documents/design/project-block-tree.md`.
///
/// This view **reads** a placement and never makes one. So it can add no placement error.
#[derive(Debug, Clone)]
pub struct ProjectBlockTree {
    pub dna: DnaType,
    /// Induced-subtree blocks in pre-order (a parent always precedes its children).
    pub blocks: Vec<Block>,
    /// The members with no placement, and the members whose terminal node this tree does not hold.
    ///
    /// The report names them and does not remove them. In a cohort from many laboratories, a
    /// difference between providers and builds is normal. Without those members, the reader can not
    /// see how much of the project the tree covers.
    pub unplaced: Vec<UnplacedMember>,
    /// The tree of this view, which is `"decodingus"` or `"ftdna"`. The `Block::loci` values belong
    /// to that tree.
    pub provider: String,
    /// The coordinate space of each position in `Block::loci`.
    ///
    /// The node names and the shape of the tree do not depend on the build. Only the positions
    /// depend on it. So the view carries the one build key that the code parsed it under, which is
    /// the most frequent build of the cohort.
    pub build_key: String,
    /// The count of groups of shared private variants that the code **removed** for a conflict.
    ///
    /// The member set of such a group shares some members with an accepted group, and the accepted
    /// group does not hold it. Two such sets can not both be a branch of one tree.
    ///
    /// The report gives this count and does not hide it. A count above zero shows recurrent calls,
    /// or a real conflict in the phylogeny of the cohort. The reader needs that fact.
    pub candidate_conflicts: usize,
    /// The count of positions that the code refused as **recurrent**. Each one would define a
    /// candidate branch below more than one parent block. So it occurred more than one time, and it
    /// can not mark a new branch.
    ///
    /// The report gives this count and does not hide it. A high value shows that the private calls
    /// of the cohort hold noise from the same cause.
    pub candidate_recurrent: usize,
}

/// One block of a [`ProjectBlockTree`]. It holds a branch and the run of SNPs that define that
/// branch. The tree treats those SNPs as equivalent. Each member below the branch carries each of
/// them, and no observation in this cohort separates them.
#[derive(Debug, Clone)]
pub struct Block {
    pub node_id: i64,
    pub name: String,
    /// Parent within the induced subtree (`None` at a root).
    pub parent: Option<i64>,
    /// The depth of this block in the subtree. The root has the value 0. The layout uses this value
    /// as its x coordinate.
    pub depth: usize,
    /// The equivalent SNPs of this block. After a collapse, the list holds the loci of each absorbed
    /// branch, with the loci nearest the root first. Inside *this cohort*, those branches are one
    /// block that nothing divides.
    pub loci: Vec<Locus>,
    /// Members whose terminal *is* this block.
    pub members: Vec<BlockMember>,
    /// The count of members at this block and below it. A collapsed branch shows that count.
    pub subtree_members: usize,
    /// Names of the member-less branches this block absorbed when collapsed (root-most first).
    /// Empty for an ordinary block. Kept so the UI can still name what it folded away.
    pub collapsed: Vec<String>,
    /// True when this block is a **candidate branch**. The published tree holds no such node. The
    /// code makes the group from the private variants that two members or more share, and those
    /// variants have no name.
    ///
    /// For a candidate, `node_id` is a value that the code made, and it is negative. The `name`
    /// field is empty, because the view supplies the label in the language of the user.
    ///
    /// A published tree can not give this answer, and the app can. The branch is real in the data,
    /// and nobody has named it yet.
    pub candidate: bool,
    /// For a candidate branch, the evidence of each carrier at each shared position. A reader can
    /// then judge the branch. The list is empty on a named block, because the tree states those SNPs
    /// and this app does not.
    pub evidence: Vec<CandidateEvidence>,
}

/// The evidence of one carrier at one shared position of a candidate branch.
///
/// A candidate is a deduction, and the statement "three men share 1 SNP" is not enough to judge it.
///
/// The read evidence behind the call of each carrier decides between a branch and a mapping
/// artefact. That evidence is the depth, and the share of the reads that hold the derived allele. A
/// cell holds one copy of that chromosome, so a true call has almost no other allele.
///
/// The aggregate carries this evidence, so a reader can judge the branch and does not trust it
/// without data.
#[derive(Debug, Clone)]
pub struct CandidateEvidence {
    pub guid: SampleGuid,
    /// Display name of the carrier.
    pub member: String,
    pub position: i64,
    pub reference: char,
    pub alternate: char,
    /// Read depth at the site; `0` when the source reported none.
    pub depth: u32,
    /// The count of reads that hold the derived allele.
    pub alt_depth: u32,
    /// The share of the reads that hold the derived allele. A cell holds one copy of chrY, so a
    /// true call gives a value near 1.0.
    pub allele_fraction: f64,
    /// Whether this call clears the federation publish gate.
    pub publishable: bool,
}

/// A project member placed at a [`Block`].
#[derive(Debug, Clone)]
pub struct BlockMember {
    pub guid: SampleGuid,
    /// The name for the display. The value is the donor identifier, and then the guid. The leaf of
    /// the tree shows it.
    pub name: String,
    /// The count of private variants below the terminal node of this member. Those variants have no
    /// name.
    ///
    /// The value is `None` until the app calculates the private-Y data of the subject. That state is
    /// not the same as `Some(0)`, which means that the app calculated the data and found no variant.
    ///
    /// Phase 3 writes this value. See `documents/design/project-block-tree.md` §9. Before that
    /// phase, the value is `None`.
    pub private_novel: Option<usize>,
    /// The part of the count above that the app can **publish**. Each such variant is new, is in
    /// unique sequence, has almost no other allele, and has enough reads. [`PublishGate`] holds
    /// those rules.
    ///
    /// The app makes a claim about a branch only from this count, and the block tree takes the mean
    /// of it. The `private_novel` count above uses weaker rules. It is a value for work in progress,
    /// and it is not a result.
    pub private_publishable: Option<usize>,
    pub private_total: Option<usize>,
}

/// A project member the block tree could not place.
#[derive(Debug, Clone)]
pub struct UnplacedMember {
    pub guid: SampleGuid,
    pub name: String,
    /// The terminal that failed to resolve against this tree; `None` = no placement at all.
    pub terminal: Option<String>,
}

/// A branch report with one row for each marker. It holds the genotype of the sample at each marker
/// that defines a node in the **subtree below** a tree node that the user chose. The lineage is the Y
/// lineage or the mtDNA lineage.
///
/// A researcher uses this report to check a placement, and to send observations to another
/// researcher.
///
/// [`App::branch_report`] builds it, and [`crate::export::branch_report_tsv`] writes it to a file.
///
/// This report is not a [`DescentReport`]. That report reads the ancestors of the placement from the
/// stored profile, from the root to the terminal node. This report genotypes the subtree again. So
/// it also holds each branch off the path where the sample is *ancestral*.
#[derive(Debug, Clone)]
pub struct BranchReport {
    pub dna: DnaType,
    /// The queried root node's haplogroup name (e.g. `R-FGC29071`).
    pub root: String,
    pub contig: String,
    /// True when the bases and the evidence came from the GVCF sidecar file of the sample. False
    /// when they came from the pileup.
    pub gvcf_backed: bool,
    pub rows: Vec<BranchRow>,
}

/// One marker that defines a branch in a [`BranchReport`]. It holds the call of the sample and the
/// evidence for that call.
#[derive(Debug, Clone)]
pub struct BranchRow {
    pub node: String,
    pub parent: String,
    pub marker: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub observed_base: Option<char>,
    pub state: CallState,
    /// The depth of the reference allele and the depth of the alternate allele, as a pair. The value
    /// is `None` on a reference block, and on the pileup path.
    pub ad: Option<(u32, u32)>,
    pub dp: Option<u32>,
    pub gq: Option<u32>,
    /// Evidence origin: `gvcf_variant` | `gvcf_refblock` | `gvcf` (uncovered) | `pileup`.
    pub source: &'static str,
    /// Human note (indel/MNV, hom-ref block, no call, …); empty for a clean call.
    pub note: String,
}

impl BranchReport {
    /// `(derived, ancestral, no-call)` marker tallies over the subtree.
    pub fn counts(&self) -> (usize, usize, usize) {
        let (mut d, mut a, mut n) = (0, 0, 0);
        for r in &self.rows {
            match r.state {
                CallState::Derived => d += 1,
                CallState::Ancestral => a += 1,
                CallState::NoCall => n += 1,
            }
        }
        (d, a, n)
    }
}

impl DescentReport {
    /// The count of SNPs on the path that define a node and that the sample carries. Each one has a
    /// derived state.
    pub fn derived(&self) -> usize {
        self.nodes
            .iter()
            .flat_map(|n| &n.snps)
            .filter(|s| s.state == CallState::Derived)
            .count()
    }

    /// The count of SNPs on the path that define a node, in each state.
    pub fn total(&self) -> usize {
        self.nodes.iter().map(|n| n.snps.len()).sum()
    }
}

/// How a private (off-backbone) variant relates to the tree.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum PrivateClass {
    /// A known SNP of the tree that is not on the path that the code assigned. It supports a finer
    /// branch, or a branch beside the assigned one.
    OffPathKnown(String),
    /// The tree does not hold this SNP. It is a candidate for a new branch.
    Novel,
}

/// A derived variant the sample carries that the haplogroup placement does not explain.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PrivateVariant {
    pub position: i64,
    pub reference: char,
    pub alternate: char,
    pub depth: u32,
    /// The count of reads that hold the derived allele, which is the alternate allele. The value
    /// `alt_depth / depth` is about equal to `allele_fraction`.
    ///
    /// The field is explicit, so the publish gate can set a minimum count of such reads. An older
    /// cached bucket holds no such field, and `serde(default)` reads it as 0. The app then
    /// calculates that bucket again.
    #[serde(default)]
    pub alt_depth: u32,
    pub allele_fraction: f64,
    pub class: PrivateClass,
    /// The structural class of this position on chrY in CHM13, from the curated list. The classes
    /// are a palindrome, an amplicon, and an AZF-DYZ region.
    ///
    /// Such a region holds paralogs, and a short read maps there without reliability. So the call is
    /// doubtful. This value is an annotation only, and the code removes no call for it.
    ///
    /// A value of `None` means unique sequence, or a build that is not CHM13.
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
    /// The calls in a curated structural region of chrY, which holds paralogs. Such a call is
    /// doubtful. A report must give it less weight, and it must not show it as a confident new
    /// variant.
    pub fn in_structural_region(&self) -> usize {
        self.variants.iter().filter(|v| v.region.is_some()).count()
    }
    /// The new calls in *unique* sequence, with no structural-region mark. These calls are the
    /// candidates for a new branch with high confidence, and they are separate from the noise of a
    /// paralog region.
    pub fn novel_in_unique_sequence(&self) -> usize {
        self.variants
            .iter()
            .filter(|v| v.class == PrivateClass::Novel && v.region.is_none())
            .count()
    }

    /// The part that the app can **publish**. Each such call is new, is in unique sequence, and also
    /// passes the strict `gate` for a new marker. That gate needs a high share of derived reads and
    /// a minimum count of reads.
    ///
    /// The app sends this set to the AppView as a set of single candidates that nobody verified.
    ///
    /// The rules are much stricter than the placement rules of the caller. So a call from a paralog,
    /// from contamination, or with little evidence never becomes a claim.
    pub fn publishable(&self, gate: PublishGate) -> Vec<&PrivateVariant> {
        self.variants.iter().filter(|v| gate.admits(v)).collect()
    }

    /// Count of [`publishable`](Self::publishable) variants under `gate`.
    pub fn publishable_count(&self, gate: PublishGate) -> usize {
        self.variants.iter().filter(|v| gate.admits(v)).count()
    }

    /// A QC banner when the novel-in-unique count is implausibly high for one sample (contamination /
    /// low coverage / reference mismatch), else `None`. See [`private_y_qc_banner`].
    pub fn qc_banner(&self) -> Option<String> {
        navigator_domain::results_context::private_y_qc_banner(self.novel_in_unique_sequence())
    }
}

/// The limits that decide which private variants the app can **publish** to the AppView as
/// candidates for a new branch.
///
/// A cell holds one copy of chrY, so a call there has almost no second allele. The code refuses an
/// allele fraction between 0.5 and 0.9, which marks a mixture or a paralog. The placement caller
/// still accepts such a fraction.
///
/// The code also refuses a call with too few reads to trust as a real single variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PublishGate {
    /// Minimum derived-allele fraction (haploid → expect ≈1.0).
    pub min_allele_fraction: f64,
    /// The minimum count of reads that hold the derived allele.
    pub min_alt_depth: u32,
}

impl Default for PublishGate {
    /// The default values for a short-read WGS sample. The call needs almost no second allele, and
    /// it needs 10 reads or more.
    fn default() -> Self {
        Self {
            min_allele_fraction: 0.9,
            min_alt_depth: 10,
        }
    }
}

impl PublishGate {
    /// The gate for the mean read length of the sample.
    ///
    /// A HiFi read, and each other long read, gives a confident haploid observation from many fewer
    /// reads. [`adaptive_min_depth`] uses the same reason. So the minimum count of reads becomes 3.
    /// The rule for the allele fraction does not change.
    pub fn for_read_len(read_len: f64) -> Self {
        let mut g = Self::default();
        if read_len > 1000.0 {
            g.min_alt_depth = 3;
        }
        g
    }

    /// Shows whether a variant passes the gate. Such a variant must be new and have no name. It
    /// must also be in unique sequence, have almost no second allele, and have enough reads.
    pub fn admits(&self, v: &PrivateVariant) -> bool {
        v.class == PrivateClass::Novel
            && v.region.is_none()
            && v.allele_fraction >= self.min_allele_fraction
            && v.alt_depth >= self.min_alt_depth
    }
}

#[cfg(test)]
mod publish_gate_tests {
    use super::*;
    use navigator_analysis::mask::YRegionClass;

    fn var(class: PrivateClass, region: Option<YRegionClass>, alt_depth: u32, af: f64) -> PrivateVariant {
        PrivateVariant {
            position: 1_000,
            reference: 'A',
            alternate: 'G',
            depth: alt_depth + 1,
            alt_depth,
            allele_fraction: af,
            class,
            region,
        }
    }

    #[test]
    fn publish_gate_admits_only_confident_unique_novels() {
        let g = PublishGate::default(); // af >= 0.9, alt_depth >= 10

        // The one that should publish: novel, unique, homozygous, deep.
        assert!(g.admits(&var(PrivateClass::Novel, None, 30, 1.0)));
        // Off-path-known is informational, never a novel-branch claim.
        assert!(!g.admits(&var(PrivateClass::OffPathKnown("M269".into()), None, 30, 1.0)));
        // Paralog-prone structural region → rejected even when deep/homozygous.
        assert!(!g.admits(&var(PrivateClass::Novel, Some(YRegionClass::Palindrome), 30, 1.0)));
        // This site has mixed alleles. The placement caller accepts a fraction of 0.5, and a
        // publish must not.
        assert!(!g.admits(&var(PrivateClass::Novel, None, 30, 0.6)));
        // A short-read sample needs more reads than this call holds.
        assert!(!g.admits(&var(PrivateClass::Novel, None, 4, 1.0)));
    }

    #[test]
    fn hifi_gate_relaxes_depth_not_fraction() {
        let g = PublishGate::for_read_len(15_000.0); // HiFi
        assert_eq!(g.min_alt_depth, 3);
        assert!(g.admits(&var(PrivateClass::Novel, None, 3, 1.0))); // 3 reads OK for HiFi
        assert!(!g.admits(&var(PrivateClass::Novel, None, 3, 0.7))); // fraction still enforced
    }

    #[test]
    fn publishable_subset_and_count_agree() {
        let bucket = PrivateBucket {
            terminal: "R-FGC29071".into(),
            variants: vec![
                var(PrivateClass::Novel, None, 30, 1.0),                         // publishable
                var(PrivateClass::Novel, Some(YRegionClass::Amplicon), 30, 1.0), // structural → no
                var(PrivateClass::OffPathKnown("Z".into()), None, 30, 1.0),      // off-path → no
                var(PrivateClass::Novel, None, 2, 1.0),                          // shallow → no
            ],
        };
        let g = PublishGate::default();
        assert_eq!(bucket.publishable_count(g), 1);
        assert_eq!(bucket.publishable(g).len(), 1);
        assert_eq!(bucket.novel_in_unique_sequence(), 2); // DISPLAY keeps the shallow one
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
pub use maintenance::{Chore, ChoreOutcome, ChoreSurvey, PrivateYRefresh, TreeReplace};
pub use navigator_domain::identity::{ExternalId, FtdnaMember, Lineage, Mdka};
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
/// [`use_os_keychain`] opts this process into the real OS keychain. The shipped binary calls it once
/// at the top of `main`; nothing else may. Without it every secret lives in an in-memory map, which
/// is what keeps the test suite off the developer's login keychain (see
/// `navigator_sync::secret_store`). [`os_keychain_enabled`] reports the current state, so tests can
/// assert they are *not* on the real keychain.
pub use navigator_sync::{os_keychain_enabled, use_os_keychain};
pub use navigator_sync::{
    AlignmentRecord, AncestralOriginRecord, BiosampleRecord, ContigMetrics, FeedPostRecord, OriginExternalId,
    PdsClient, PopulationBreakdownRecord, PrivateVariantsRecord, RecordRef, SequenceRunRecord, VariantCallEntry,
    ANCESTRAL_ORIGIN_COLLECTION, NS_ALIGNMENT, NS_BIOSAMPLE, NS_FEED_POST, NS_POPULATION_BREAKDOWN, NS_SEQUENCERUN,
    PRIVATE_VARIANTS_COLLECTION,
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
    /// The count of sites with a call in **both** samples. That count is the true size of the
    /// comparison.
    ///
    /// A small overlap makes a call on a short segment weak. Two chips give such an overlap, and a
    /// chip against a WGS sample also gives one, because the chip sites limit it. The report gives
    /// this count and does not hide it.
    pub overlapping_sites: usize,
}

/// One sample of an IBD comparison. It is a WGS **alignment** in a CRAM file, which the code
/// genotypes at the IBD-panel sites. It can also be a **chip** profile that a user imported, which
/// the code re-keys to the same CHM13 sites.
///
/// Both forms give a dosage at each site of the canonical IBD panel. So the comparison does not
/// depend on the kind of data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IbdSource {
    Alignment(i64),
    Chip(i64),
    /// A genome-wide imported variant set (a WGS VCF / CompleteGenomics masterVar), resolved to
    /// panel dosages with unlisted sites taken as homozygous reference.
    VariantSet(i64),
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

/// A federated-IBD candidate from the match engine of the AppView. The candidate has no name.
///
/// The `suggested_sample_guid` value is the opaque handle of the AppView for the other person. It is
/// not a DID, and it holds no personal data. The app uses it to ask for an introduction.
///
/// The `signals` field names each source behind the `score` value, such as `POPULATION_OVERLAP`,
/// `HAPLOGROUP`, and `SHARED_MATCH`.
///
/// The `target_sample_guid` value is the handle of the AppView for **our own** sample, and the
/// engine ranked the candidate against that sample. We already own it, so it gives away nothing.
///
/// But a client that publishes its own records has no other way to learn its handle on the server.
/// [`App::ibd_attest`] can not report a complete comparison without it. The value is `None` on an
/// AppView from before that field.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdSuggestion {
    pub target_sample_guid: Option<String>,
    pub suggested_sample_guid: String,
    pub suggestion_type: String,
    pub score: f64,
    pub signals: Vec<String>,
}

/// The level of trust in the score of a federated-IBD candidate. The "Genetic relatives" card of
/// Simple mode uses these words.
///
/// This value reads the evidence, and it does not only draw it. The line between "strong" and
/// "possible" decides which of three statements the app makes about the relationship between a
/// stranger and the user.
///
/// This code is beside [`IbdSuggestion`] and not in the card that draws it. So the rule has one
/// home. A later change to it changes the reading of the evidence, and not a widget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchStrength {
    /// The signals agree strongly; presented as a likely relative.
    Strong,
    /// The two samples agree enough for the user to act.
    Likely,
    /// Weak or single-signal evidence; presented as a possibility only.
    Possible,
}

impl IbdSuggestion {
    /// Change the `score` of this candidate into the level that the UI names.
    ///
    /// The AppView gives a score from 0 to 1, and that score joins each of the `signals` values. So
    /// the limits here are careful. The method names a candidate strong only when the evidence is
    /// far above the middle of that range.
    ///
    /// A statement that is too strong makes a user write to a stranger on weak evidence.
    pub fn strength(&self) -> MatchStrength {
        if self.score >= 0.8 {
            MatchStrength::Strong
        } else if self.score >= 0.5 {
            MatchStrength::Likely
        } else {
            MatchStrength::Possible
        }
    }
}

/// The result of a request for an introduction to a candidate. It holds the request URI of the
/// AppView and the status of that request. The first status is `PENDING`, and the two parties then
/// exchange their consent.
///
/// The server chooses the `purpose` value from the strongest signal of the suggestion. The values
/// are `IBD_AUTOSOMAL`, `IBD_Y`, and `IBD_MT`.
///
/// That value decides the genomic region of a later attestation. So the app records it at the
/// introduction, and it does not wait for the session to give it.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdIntroResult {
    pub request_uri: String,
    pub status: String,
    pub purpose: String,
}

/// An exchange request that arrived and that needs the consent of this account. The view is
/// **symmetric-blind**: the app does not see the sender until both parties agree. The value comes
/// from `GET /api/v1/exchange/incoming`.
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

/// The result of `POST /api/v1/exchange/consent`. The value is `CONSENTED`, with the `session_id`
/// of the new session. It can also be `DECLINED`. It can also be `PENDING`, which means that the
/// server recorded our answer and waits for the other party.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsentOutcome {
    pub status: String,
    pub session_id: Option<String>,
}

/// Who opened a matching conversation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchingDirection {
    /// This account asked for the introduction.
    Outbound,
    /// Another person asked for an introduction to this account.
    Inbound,
}

impl MatchingDirection {
    /// Stable ledger token (independent of display strings).
    pub fn as_str(self) -> &'static str {
        match self {
            MatchingDirection::Outbound => "OUTBOUND",
            MatchingDirection::Inbound => "INBOUND",
        }
    }
    fn parse(s: &str) -> Self {
        match s {
            "INBOUND" => MatchingDirection::Inbound,
            _ => MatchingDirection::Outbound,
        }
    }
}

/// The state of a matching conversation. The value records only what this device can *know*.
///
/// The broker is symmetric-blind. So a partner who declined looks the same as a partner who did not
/// answer, and both stay at [`MatchingStatus::Requested`]. The value
/// [`MatchingStatus::Declined`] means that **this account** declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchingStatus {
    /// We asked; the counterpart has not consented (or has not answered).
    Requested,
    /// The request arrived, and this account must decide.
    AwaitingConsent,
    /// We declined.
    Declined,
    /// Both parties agreed. A session is open, and the encrypted exchange can run.
    Ready,
    /// The exchange ran, and the store holds a result.
    Exchanged,
    /// The exchange ran and failed. The `last_error` field gives the reason.
    Failed,
}

impl MatchingStatus {
    /// Stable ledger token (independent of display strings).
    pub fn as_str(self) -> &'static str {
        match self {
            MatchingStatus::Requested => "REQUESTED",
            MatchingStatus::AwaitingConsent => "AWAITING_CONSENT",
            MatchingStatus::Declined => "DECLINED",
            MatchingStatus::Ready => "READY",
            MatchingStatus::Exchanged => "EXCHANGED",
            MatchingStatus::Failed => "FAILED",
        }
    }
    fn parse(s: &str) -> Self {
        match s {
            "AWAITING_CONSENT" => MatchingStatus::AwaitingConsent,
            "DECLINED" => MatchingStatus::Declined,
            "READY" => MatchingStatus::Ready,
            "EXCHANGED" => MatchingStatus::Exchanged,
            "FAILED" => MatchingStatus::Failed,
            _ => MatchingStatus::Requested,
        }
    }
    /// True when the conversation has no more work. It gave a result, or this account declined it.
    /// The UI keeps such a conversation out of the list of actions.
    pub fn is_terminal(self) -> bool {
        matches!(self, MatchingStatus::Exchanged | MatchingStatus::Declined)
    }
}

/// One matching conversation, as the UI reads it: the durable ledger row plus the exchange result
/// once there is one. Assembled by [`App::matching_entries`].
#[derive(Debug, Clone)]
pub struct MatchingEntry {
    pub request_uri: String,
    pub direction: MatchingDirection,
    pub purpose: String,
    pub status: MatchingStatus,
    /// The server gives this value only after both parties agree. It is `None` while the request
    /// stays blind.
    pub partner_did: Option<String>,
    pub session_id: Option<String>,
    /// The local subject whose dosages this conversation exchanges.
    pub biosample_guid: Option<SampleGuid>,
    /// The two sample handles of the AppView, ours and theirs. An attestation uses them as its
    /// key.
    pub my_sample_ref: Option<String>,
    pub partner_sample_ref: Option<String>,
    /// Our own consent decision; `None` until we make one.
    pub consent_given: Option<bool>,
    /// True once the AppView accepted our attestation for this comparison.
    pub attested: bool,
    pub last_error: Option<String>,
    pub created_at: String,
    pub updated_at: String,
    /// The stored exchange result, present once `status` is [`MatchingStatus::Exchanged`].
    pub result: Option<StoredIbdExchange>,
}

/// A pulled relay envelope: the opaque ciphertext `blob`, its route fields (`from_did`/`seq`), and
/// the broker `id` to ack. From `GET /api/v1/exchange/relay/pull`.
#[derive(Debug, Clone, PartialEq)]
pub struct RelayEnvelope {
    pub id: i64,
    pub from_did: String,
    pub seq: i32,
    pub blob: String,
}

/// A live exchange session with a derived shared key, ready to seal/open payloads. It holds key
/// material. So it is deliberately not `Debug`/`Serialize`, and you must keep it in memory only.
#[derive(Clone)]
pub struct EstablishedSession {
    pub session_id: String,
    pub partner_did: String,
    key: [u8; 32],
}

/// Jetstream-ingest retry budget for a device key that the app published a moment ago. A 403
/// immediately after the publish means the AppView has not ingested our `deviceKey` record yet.
/// Exponential backoff of 1+2+4+8 s, about 15 s in total, before the app stops.
const DEVICE_KEY_INGEST_RETRIES: u32 = 4;

/// Poll rounds (≈1s each) an IBD exchange waits for the partner's dosages / attestation.
const EXCHANGE_POLL_ROUNDS: u32 = 30;

/// The fed-record collections this client publishes. A PULL reconcile scans each one, and the set
/// mirrors the `publish_*` NSIDs. The app tracks the derived-summary collections, but it does not
/// overwrite them locally.
const PUBLISHED_COLLECTIONS: &[&str] = &[
    NS_BIOSAMPLE,
    NS_ALIGNMENT,
    NS_POPULATION_BREAKDOWN,
    NS_SEQUENCERUN,
    PRIVATE_VARIANTS_COLLECTION,
    HAPLOGROUP_RECONCILIATION_COLLECTION,
];

/// PDS collection NSID for a published IBD match attestation (the AppView indexes these through
/// Jetstream).
const IBD_ATTESTATION_COLLECTION: &str = "com.decodingus.atmosphere.ibdAttestation";

/// Above this many sites, the app decimates the exchanged dosage vector to fit the relay's 1 MiB
/// envelope.
const EXCHANGE_SITE_BUDGET: usize = 100_000;
/// Decimation stride when over budget: keep sites at `position % N == 0`. The rule is
/// **position-based** and not index-based, so both peers keep the *same physical sites*. This keeps
/// the IBD intersection even when their panels differ in size (WGS against chip). The result is
/// about 1/N of the canonical panel.
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

/// Parse the AppView's `/api/v1/ibd/suggestions` body into [`IbdSuggestion`]s. The parser accepts
/// both field casings (camel and snake) and both `signals` shapes (an object map or an array). A
/// small change to the contract then loses only some fields, and it does not drop every candidate.
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
            // Optional: older AppViews omit it, and the row is still usable for everything except
            // an attestation. So a missing value must not drop the candidate.
            let target_sample_guid = it
                .get("targetSampleGuid")
                .or_else(|| it.get("target_sample_guid"))
                .and_then(|v| v.as_str())
                .map(str::to_string);
            let score = it.get("score").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let signals = it
                .get("metadata")
                .and_then(|m| m.get("signals"))
                .or_else(|| it.get("signals"))
                .map(parse_ibd_signals)
                .unwrap_or_default();
            Some(IbdSuggestion {
                target_sample_guid,
                suggested_sample_guid,
                suggestion_type,
                score,
                signals,
            })
        })
        .collect()
}

/// Signal names. The AppView emits an array of plain strings, as in
/// `["POPULATION_OVERLAP", "HAPLOGROUP"]`. The parser also accepts an array of `{name|source}`
/// objects, or an object map whose keys are the names. A change to the contract then loses less.
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

pub use navigator_analysis::ibd_attest::{IbdAttestation, IbdExchangeMsg, IbdSite};
use navigator_domain::bisdna;
pub use navigator_domain::brief::{
    AncestryBrief, AncientComponent, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief,
};
use navigator_domain::chipprofile::{self, ChipProfile, NewChipProfile};
pub use navigator_domain::consensus::{DiploidSourceObs, DiploidVariant};
use navigator_domain::filetype;
pub use navigator_domain::filetype::DetectedData;
use navigator_domain::mtdna::{self, MtdnaSequence, NewMtdnaSequence};
use navigator_domain::reconciliation::{self, CallProvenance, RunHaplogroupCall};
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
pub use navigator_store::ibd_request::StoredIbdRequest;
pub use navigator_store::source_file::SourceFile;
use navigator_store::{
    alignment, ancestry_result, artifact, biosample, biosample_project, chip_profile, consensus_profile,
    haplogroup_call, mdka, mtdna as mtdna_store, project, reconciliation as recon_store, sequence_run, sig_cache,
    source_file, str_profile, sync_history, sync_outbox, sync_state, variant_set, variant_set_genotype,
    variant_set_private_y, Store, StoreError,
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
pub use update::UpdateInfo;

/// Artifact kind for de-novo calls, keyed by contig so different contigs do not
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
        .unwrap_or_else(|| navigator_domain::paths::decodingus_dir().join("trees"));
    dir.join(file)
}

/// Sidecar path that holds the HTTP `ETag` of a cached haplotree (`<cache>.etag`).
/// [`App::fetch_tree`] sends it back as `If-None-Match` on a refresh. An unchanged tree then
/// returns a small `304` instead of the full 60 to 127 MB body.
fn tree_etag_path(cache_path: &Path) -> PathBuf {
    let mut p = cache_path.as_os_str().to_owned();
    p.push(".etag");
    PathBuf::from(p)
}

/// How long the app trusts a cached haplotree before [`App::fetch_tree`] downloads it again. The
/// AppView's curated tree changes slowly (curator review, periodic builds). So a weekly refresh
/// keeps placements current, and it does not touch the network on every run. Override with
/// `NAVIGATOR_TREE_TTL_DAYS` (0 = always refetch).
const TREE_CACHE_TTL_DAYS_DEFAULT: u64 = 7;

/// Whole-request timeout for a haplotree download. reqwest's `.timeout()` bounds the *entire*
/// request, and it includes the body read. The trees are large: the DecodingUs Y tree is about
/// 60 MB, and the FTDNA Y tree about 127 MB. A short cap stops the body read part of the way
/// through. reqwest reports that as "error decoding response body". A refresh can then never
/// complete, and the cache stays stale for ever.
///
/// The value is large on purpose. A cache that is present makes any failure fall back immediately
/// (see [`App::fetch_tree`]). So only a first fetch with no cache can wait this long.
const TREE_DOWNLOAD_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(300);

/// Is the cached tree at `path` still within its TTL (default 7 days; `NAVIGATOR_TREE_TTL_DAYS`
/// overrides)? Unknown mtime or unreadable metadata → not fresh, which forces a refresh.
fn tree_cache_is_fresh(path: &Path) -> bool {
    let days = std::env::var("NAVIGATOR_TREE_TTL_DAYS")
        .ok()
        .and_then(|v| v.trim().parse::<u64>().ok())
        .or_else(|| AppSettings::load().tree_ttl_days)
        .unwrap_or(TREE_CACHE_TTL_DAYS_DEFAULT);
    brief::cache_is_fresh(path, days)
}

/// Score a tree against the sample calls and attach the terminal's child-branch evidence.
///
/// The Kulczynski `score` ranks the candidates by proportional similarity, and it supplies the
/// alternatives list. But two steps choose the *reported terminal*.
///
/// The first step takes the best-ranked candidate that the path-supported parsimony guard admits.
/// The guard rejects a candidate whose lineage tunnels through a branch that the sample
/// contradicts, which is the distal-Y paralog artifact.
///
/// The second step calls [`haplo::deepen_terminal`], which descends into any child that the sample
/// clearly entered. This corrects an under-call at an **unsplit tree node**, where a half-ancestral
/// SNP block scores below its parent.
///
/// The function moves the chosen node to the front, so every `ranked.first()` consumer gets it.
/// See `documents/design/PangenomeExpansion.md`.
/// Pool every source's vote into one consensus map by a `SourceType`-weighted majority.
///
/// The key `K` is the SNP **name** for Y, which is portable across builds. For mt it is the rCRS
/// **position**. The value `V` is a **state** for each Y SNP, which is independent of strand and
/// build. CHM13 and GRCh38 can flip a base, but neither changes whether the sample carries the
/// derived allele. For mt the value is a **base**, because mt has one coordinate system.
///
/// The weight matches the `SourceType` term of the variant reconcile, in
/// [`navigator_domain::consensus::obs_weight`]. The value with the highest weight wins for each
/// key. The function places the pooled set on the tree **once**, at the genome level. It does not
/// vote among the terminal labels of each run.
fn pool_votes<K, V>(sources: &[(SourceType, HashMap<K, V>)]) -> HashMap<K, V>
where
    K: std::hash::Hash + Eq + Clone,
    V: std::hash::Hash + Eq + Clone + Ord,
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
                // The highest weight wins. On a tie, break by the allele itself, so the pooled
                // call is deterministic. Without that tie-break the `HashMap` iteration order
                // picked the winner at random, which flipped the placed terminal between runs
                // over identical genotypes.
                .max_by(|a, b| {
                    a.1.partial_cmp(&b.1)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then(a.0.cmp(&b.0))
                })
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

/// Terminal selection for **named Y-SNP panel** data (BISDNA chip), against the alignment-tuned
/// [`assemble_assignment`].
///
/// Such panels give confident but sparse genotype calls. A few recurrent or mis-probed ancestral
/// calls on backbone nodes can make the strict `path_admissible` guard veto the genuine deep
/// lineage. The call then drops to a shallow node such as A1. That guard exists to kill distal
/// tunnel artifacts in alignment data with limited coverage.
///
/// With confident chip calls that failure mode dominates. So here the code trusts the top of the
/// proportional Kulczynski rank, which survives a few stray calls. It then calls
/// [`deepen_terminal`] on the children that the sample clearly entered.
///
/// Checked: this kit's chromo2 export gives R-S1121 on both the DecodingUs/hs1 tree and the
/// FTDNA/GRCh38 tree, on the lineage to its WGS-confirmed R-FGC29071.
fn assemble_assignment_robust(
    tree: &navigator_analysis::haplo::HaploTree,
    calls: &HashMap<i64, char>,
) -> HaploAssignment {
    use navigator_analysis::haplo;
    let mut ranked = haplo::score(tree, calls);
    if let Some(top_id) = ranked.first().map(|r| r.id) {
        let terminal_id = haplo::deepen_terminal(tree, calls, top_id);
        // Parsimony back-off: do not report a terminal deeper than the evidence supports. Trim any
        // net-contradicted tail of the lineage, which a sparse panel or damaged aDNA can make too
        // deep. A lone contradiction that deeper derived support outweighs still reaches the deep
        // terminal.
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

/// The root→`target` path of node ids (inclusive), or empty if `target` is not reachable.
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

/// Root→`name` lineage of haplogroup names from the tree (empty if the name is not found). Used to
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

/// Back off an over-deep terminal to the node with the maximum cumulative support along its
/// lineage.
///
/// The walk goes from the root to the terminal. Each node contributes
/// `(covered derived − covered ancestral)` over the SNPs that define it and that the sample has a
/// call for. The chosen terminal is the deepest node at which that cumulative balance peaks.
///
/// The function trims a net-contradicted tail, which has more ancestral than derived calls. A
/// sparse chip or a degraded aDNA sample makes such a tail when it tunnels into a wrong sub-clade.
/// But the function keeps a tail whose deeper derived calls outweigh a shallow contradiction. A
/// tie favours the deeper node, which keeps the "survive a lone backbone contradiction" behaviour.
///
/// Returns `terminal_id` unchanged when the code can not trace its lineage.
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
        // derived-supported. A contradiction that a deeper derived call recovers then still
        // reaches the deep terminal. A net-negative tail, or a flat run of nodes with no marker,
        // goes back to the last node with positive support. That flat run is the sparse-panel or
        // aDNA tunnel.
        if balance > best_balance || (balance == best_balance && node_derived) {
            best_balance = balance;
            best_id = id;
        }
    }
    best_id
}

/// Reconcile chip genotype calls to a haplotree's strand.
///
/// Consumer arrays report alleles on the reference plus strand. But some sites sit on the opposite
/// strand from the ancestral/derived convention of the tree.
///
/// For each call at a tree position, keep the observed base when it already equals the ancestral
/// or the derived allele. If not, use its complement when *that* base matches. If neither matches,
/// keep the observed base, and the score counts it against the branch. A position that the tree
/// does not have passes through unchanged, and it changes no score.
///
/// This does nothing for BISDNA calls that the dictionary reconciled, because their base is always
/// the derived allele. So it is safe on the shared chip-placement path.
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

/// Map GVCF-decoded bases at *lifted* positions back to tree positions. This is the GVCF form of
/// [`App::build_calls_from_lifted`].
///
/// A variant base wins. If there is none, a callable hom-ref lifted site takes the **reference
/// base** at that lifted position. A minus-strand lift reverse-complements both. Any other
/// position is a no-call.
///
/// `ref_base` uses the lifted position as its key, which is the GVCF or reference coordinate, and
/// not the tree position.
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

/// Minimum callable depth, adapted to the read technology. The default of 4 is a short-read
/// assumption: about 4 reads to call a base with confidence.
///
/// Long, accurate reads (HiFi, with a mean read length above 1 kb) give a confident haploid
/// observation from a *single* read. So a HiFi sample at about 4x is callable at 1x. A
/// floor clamped at 2 threw away half of its already shallow coverage for no gain.
///
/// ONT long reads are less accurate. Look at this again if the code ever adapts by platform and
/// not by read length.
fn adaptive_min_depth(base: u32, read_len: f64) -> u32 {
    if read_len > 1000.0 {
        1
    } else {
        base
    }
}

/// Haploid-caller params adapted to the sample's read technology (see [`adaptive_min_depth`]).
/// The function samples the head of the BAM, and falls back to the defaults on any error. It
/// blocks, because it reads the BAM, so call it inside `spawn_blocking`.
fn adaptive_haploid_params(bam_path: &Path, reference: Option<&Path>) -> HaploidCallerParams {
    let mut params = HaploidCallerParams::default();
    if let Ok((read_len, _)) = coverage::estimate_molecule_lengths(bam_path, reference) {
        params.min_depth = adaptive_min_depth(params.min_depth, read_len);
    }
    params
}

/// Minimum genotyped sites for a reliable AIMs ancestry estimate (Scala `minSnpsAims`).
/// `$NAVIGATOR_ANCESTRY_MIN_SNPS` overrides it (tests use a small panel).
fn ancestry_min_snps() -> usize {
    std::env::var("NAVIGATOR_ANCESTRY_MIN_SNPS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3000)
}

/// Resolve an ancestry/IBD asset path under `<refgenome base>/ancestry/`: an `$<env_var>` override
/// (when non-empty) wins, else `<base>/ancestry/<stem>_<build>.<ext>`. The wrapper for each asset
/// below delegates here, so the override, join, and format pattern lives in one place.
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
/// `<refgenome base>/ancestry/ancestry_pca_<build>.bin`. Optional. When it is absent, the
/// AF-likelihood estimate runs with no PCA coordinates.
fn ancestry_pca_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_PCA", "ancestry_pca", build, "bin")
}

/// The fine-population frequency asset path (`$NAVIGATOR_ANCESTRY_FREQ` override, else
/// `<base>/ancestry/ancestry_freq_global_<build>.bin`). Optional. Without it, the app skips fine
/// admixture.
fn ancestry_freq_global_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_FREQ", "ancestry_freq_global", build, "bin")
}

/// The phased-haplotype reference asset (`$NAVIGATOR_ANCESTRY_HAPS` override, else
/// `<base>/ancestry/ancestry_haps_<build>.bin`): the phased 1000G haplotypes the statistical phaser
/// copies from, for the parent-split chromosome painter. Optional. When it is absent, the painter
/// falls back to the unphased diploid path, which gives two sorted copies and not parental sides.
fn ancestry_haps_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_HAPS", "ancestry_haps", build, "bin")
}

/// The **ancient** deep-source frequency asset (`$NAVIGATOR_ANCESTRY_FREQ_ANCIENT` override, else
/// `<base>/ancestry/ancestry_freq_ancient_<build>.bin`): WHG/ANF/Steppe alt-allele frequencies for
/// each site, which `panelbuild ancient-panel` builds from the AADR. Optional. Without it, the app
/// skips deep ancestry. It covers the AIM panel's own sites, so one genotyping pass supplies it.
///
/// This replaces the old `ancestry_pca_ancient_<build>.bin`, which nothing reads now. Ancient
/// centroids that a PCA projects collapse onto the modern cloud and carry no ancient signal.
fn ancestry_freq_ancient_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_FREQ_ANCIENT", "ancestry_freq_ancient", build, "bin")
}

/// The qpAdm deep-ancestry panel asset path (`$NAVIGATOR_ANCESTRY_QPADM` override, else
/// `<base>/ancestry/ancestry_qpadm_<build>.bin`). The full-1240k Patterson-2022 config, with
/// WHG/EEF/Steppe sources and sister outgroups. See
/// documents/design/ancient-ancestry-rebuild.md §7.14.
fn ancestry_qpadm_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ANCESTRY_QPADM", "ancestry_qpadm", build, "bin")
}

/// The archaic (Neanderthal / Denisovan) marker panel asset path
/// (`$NAVIGATOR_ARCHAIC_MARKERS` override, else `<base>/ancestry/archaic_markers_<build>.bin`).
/// `panelbuild archaic-panel` builds it. See documents/design/ArchaicAncestry_Design.md §4.
fn archaic_markers_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ARCHAIC_MARKERS", "archaic_markers", build, "bin")
}

/// The archaic percentile reference asset path (`$NAVIGATOR_ARCHAIC_DIST` override, else
/// `<base>/ancestry/archaic_marker_dist_<build>.bin`).
fn archaic_marker_dist_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ARCHAIC_DIST", "archaic_marker_dist", build, "bin")
}

/// Tier B: positions variable in the African outgroup, which let the code remove shared variants.
/// Cache signature for a Tier B segment result: the alignment it came from, plus the genotype
/// version of the caller. So a newer caller invalidates the result.
pub(crate) fn archaic_segment_sig(alignment_id: i64, called_contigs: &[String]) -> String {
    // Three things make a result stale, and all three are in the key.
    //
    // The METHOD version, because Tier B changed from a private-variant density model to a match
    // against the archaic genomes. Without it, a workspace that holds output from the withdrawn
    // caller would continue to serve that output.
    //
    // The CONTIGS THE CODE CALLED, because the result covers only those. Take a subject called on
    // chr21 alone, and later called genome-wide. Without this term the subject would keep the
    // two-chromosome answer for ever. A genome-wide check showed exactly that: a report of 1.94 Mb
    // over 2 contigs, while 22 more sat in the cache and ready.
    let mut contigs: Vec<&str> = called_contigs.iter().map(String::as_str).collect();
    contigs.sort_unstable();
    format!(
        "aln{alignment_id}:gt{}:m{}:c[{}]",
        navigator_analysis::caller::GENOTYPE_VERSION,
        navigator_analysis::archaic_match::METHOD_VERSION,
        contigs.join(",")
    )
}

/// The autosomes that have cached de-novo diploid calls for `alignment_id`, as bare contig names.
pub(crate) async fn called_diploid_contigs(
    store: &navigator_store::Store,
    alignment_id: i64,
) -> Result<Vec<String>, AppError> {
    const PREFIX: &str = "diploid_denovo:";
    Ok(navigator_store::artifact::list_kinds(store.pool(), alignment_id)
        .await?
        .into_iter()
        .filter_map(|k| k.strip_prefix(PREFIX).map(str::to_string))
        .collect())
}

fn archaic_outgroup_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ARCHAIC_OUTGROUP", "archaic_outgroup_af", build, "bin")
}

/// Tier B: genome-wide archaic diagnostic sites, which give a label to each called segment.
fn archaic_classify_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ARCHAIC_CLASSIFY", "archaic_classify", build, "bin")
}

/// Tier B: the count of callable bases in each window. Without this the segment HMM calls mapping
/// artifacts.
fn archaic_callable_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_ARCHAIC_CALLABLE", "archaic_callable", build, "bin")
}

/// The chip-compatible IBD panel asset path (`$NAVIGATOR_IBD_PANEL` override, else
/// `<base>/ancestry/ibd_panel_<build>.bin`).
fn ibd_panel_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_IBD_PANEL", "ibd_panel", build, "bin")
}

/// The ancestry/IBD reference assets for the analysis build (CHM13). Each one carries its presence
/// and its manifest check, which is the "data sources" transparency line. This looks only at the
/// file system, and does no analysis.
pub fn ancestry_asset_status() -> Vec<AssetStatus> {
    let build = ReferenceBuild::Chm13v2;
    let manifest = load_asset_manifest(build);
    [
        ("super-pop panel", ancestry_panel_path(build)),
        ("PCA (modern)", ancestry_pca_path(build)),
        ("fine frequencies", ancestry_freq_global_path(build)),
        ("ancient frequencies", ancestry_freq_ancient_path(build)),
        ("genetic map", genetic_map_path(build)),
        ("IBD panel", ibd_panel_path(build)),
        ("archaic markers", archaic_markers_path(build)),
        ("archaic percentiles", archaic_marker_dist_path(build)),
        ("archaic outgroup", archaic_outgroup_path(build)),
        ("archaic classify", archaic_classify_path(build)),
        ("archaic callable", archaic_callable_path(build)),
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

/// Copy every regular file in `src_dir` into `dest_dir` that is not already present there. It never
/// overwrites a file that exists, because an asset refreshed from the CDN must win over the bundled
/// one. It creates `dest_dir`. A `src_dir` that is missing or unreadable does nothing, and returns
/// the empty summary. The function is pure over the two directories, with no globals, so a unit
/// test can drive it.
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
        // Skip hidden files, such as the `.staged` marker, and documents such as a bundled README.
        // The seed copies only real data assets.
        if name.to_string_lossy().starts_with('.') || path.extension().and_then(|e| e.to_str()) == Some("md") {
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

/// Locate the bundled ancestry-asset resource directory inside the installed image:
/// `$NAVIGATOR_BUNDLED_ASSETS` (override), else candidates relative to the current executable for
/// each packaged layout (macOS `.app` Resources, Linux `usr/lib|share/<app>`, Windows alongside).
/// `None` for a dev `target/` build with no bundle.
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

/// Locate the bundled chrY-mask resource directory (`masks/`), which holds the private-Y filter
/// assets: the callable mask and the cohort-shared exclude list. The resolution is the same as
/// [`bundled_assets_dir`], plus a compile-time repo path. So a `cargo run` dev build seeds
/// straight from the `assets/masks/` in the repo. These masks are small enough to live in git,
/// unlike the ancestry `.bin` files.
fn bundled_masks_dir() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_BUNDLED_MASKS") {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
    }
    // Dev build: the checked-in repo assets (exists on a developer's machine, not in a package).
    let repo_assets = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../assets/masks"));
    if repo_assets.is_dir() {
        return Some(repo_assets);
    }
    let exe = std::env::current_exe().ok()?;
    let dir = exe.parent()?;
    [
        dir.join("../Resources/masks"),         // macOS .app/Contents/MacOS → ../Resources
        dir.join("masks"),                      // Windows (alongside) / portable
        dir.join("../lib/DUNavigator/masks"),   // Linux .deb/AppImage usr/bin → usr/lib/<app>
        dir.join("../share/DUNavigator/masks"), // Linux usr/share/<app>
        dir.join("resources/masks"),            // generic
    ]
    .into_iter()
    .find(|c| c.is_dir())
}

/// Seed the bundled chrY private-Y masks into `<cache base>/masks/` on first run (gzipped BEDs;
/// [`RegionMask::from_bed`](navigator_analysis::mask::RegionMask::from_bed) reads them transparently).
/// Never overwrites, so a user's own uncompressed override wins. Best-effort ⇒ empty summary when no
/// source is present.
pub fn seed_bundled_masks() -> SeedSummary {
    let Some(src) = bundled_masks_dir() else {
        return SeedSummary::default();
    };
    let dest = refgenome_cache::base_dir().join("masks");
    seed_assets_from(&src, &dest).unwrap_or_default()
}

/// Locate the bundled STR-reference resource directory (`str/`), which holds the HipSTR reference
/// BEDs. The resolution is the same as [`bundled_masks_dir`], but with **no** compile-time repo
/// fallback. The GRCh38 HipSTR BED is about 20 MB, which is too large for git. So it ships in the
/// installer image, staged there from `~/.decodingus/str/` or from the asset release.
fn bundled_str_dir() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_BUNDLED_STR") {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
    }
    let exe = std::env::current_exe().ok()?;
    let dir = exe.parent()?;
    [
        dir.join("../Resources/str"),       // macOS .app/Contents/MacOS → ../Resources
        dir.join("str"),                    // Windows (alongside) / portable
        dir.join("../lib/DUNavigator/str"), // Linux .deb/AppImage usr/bin → usr/lib/<app>
        dir.join("../share/DUNavigator/str"),
        dir.join("resources/str"),
    ]
    .into_iter()
    .find(|c| c.is_dir())
}

/// Seed the bundled HipSTR reference BEDs into `<cache base>/str/` on the first run, so STR calling
/// works immediately for the shipped builds. It never overwrites, because a user's own override
/// wins. It is best-effort: it returns the empty summary when no bundle is present, as in a lean
/// dev build.
pub fn seed_bundled_str() -> SeedSummary {
    let Some(src) = bundled_str_dir() else {
        return SeedSummary::default();
    };
    let dest = refgenome_cache::base_dir().join("str");
    seed_assets_from(&src, &dest).unwrap_or_default()
}

/// Seed every bundled asset directory into the cache: the ancestry panels, the chrY masks, and the
/// STR references. The function is idempotent, because it skips a destination that already exists.
/// So only a first run copies.
pub fn seed_bundled_all() -> SeedSummary {
    let mut total = SeedSummary::default();
    for s in [seed_bundled_assets(), seed_bundled_masks(), seed_bundled_str()] {
        total.copied += s.copied;
        total.skipped += s.skipped;
    }
    total
}

/// The in-flight background seed started by [`spawn_bundled_seed`], if any.
static BUNDLED_SEED: std::sync::Mutex<Option<std::thread::JoinHandle<SeedSummary>>> = std::sync::Mutex::new(None);

/// Run [`seed_bundled_all`] on a background thread, so a first-run copy does not sit in front of
/// the GUI's first frame. The GRCh38 HipSTR BED alone is about 20 MB, and that delay lands exactly
/// when a new user looks for the window. Anything that reads the assets must call
/// [`await_bundled_assets`] first. [`App::open`] does, which covers every data path.
///
/// A headless run seeds synchronously instead, because its analysis starts immediately and it has
/// no window to show. So a headless run never calls this, and the await below does nothing.
pub fn spawn_bundled_seed() {
    let mut slot = BUNDLED_SEED.lock().unwrap_or_else(|e| e.into_inner());
    if slot.is_none() {
        *slot = Some(std::thread::spawn(seed_bundled_all));
    }
}

/// Block until a [`spawn_bundled_seed`] has finished. It returns immediately when no seed started,
/// or when a caller already awaited it. The function holds the lock across the join on purpose. A
/// second caller then waits for the first. It does not race past a half-copied asset directory.
pub fn await_bundled_assets() -> SeedSummary {
    let mut slot = BUNDLED_SEED.lock().unwrap_or_else(|e| e.into_inner());
    match slot.take() {
        Some(h) => h.join().unwrap_or_default(),
        None => SeedSummary::default(),
    }
}

/// The asset integrity manifest path for a build (`<base>/ancestry/ancestry_manifest_<build>.json`).
fn ancestry_manifest_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("", "ancestry_manifest", build, "json")
}

/// Load the build's asset manifest, if the build has one. `None` means the manifest is absent or
/// unparseable, and the app then skips the integrity checks. They are advisory.
fn load_asset_manifest(build: ReferenceBuild) -> Option<navigator_analysis::manifest::AssetManifest> {
    std::fs::read_to_string(ancestry_manifest_path(build))
        .ok()
        .and_then(|s| navigator_analysis::manifest::AssetManifest::from_json(&s).ok())
}

/// Read an asset file (`None` if absent), and check its SHA-256 against the build manifest when the
/// build has one. A **checksum mismatch is a hard error**: refuse a corrupt or truncated asset, and
/// do not analyze against it. A manifest that is missing, or a file the manifest does not list,
/// passes through with no check.
///
/// The code also quarantines a mismatched file in the **managed cache**, with a rename to
/// `.corrupt`, so the next [`App::ensure_ancestry_asset`] downloads a good copy. Without that, a
/// message that tells the user to "re-download it" asks for something the app gives them no way to
/// do. The error then repeats for ever. The code never touches a file at a user-specified
/// `$NAVIGATOR_*` override. That file is theirs, and a mismatch there most probably means their
/// manifest is stale, and not their asset bad.
fn read_verified_asset(build: ReferenceBuild, path: &Path) -> Result<Option<Vec<u8>>, AppError> {
    let Ok(bytes) = std::fs::read(path) else {
        return Ok(None);
    };
    if let Some(manifest) = load_asset_manifest(build) {
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if let Err((expected, got)) = manifest.verify(name, &bytes) {
                let managed = refgenome_cache::base_dir().join("ancestry").join(name);
                let quarantined = path == managed && std::fs::rename(path, path.with_extension("corrupt")).is_ok();
                return Err(AppError::Import(format!(
                    "asset {name} failed its integrity check (manifest sha256 {expected}, file {got}){}",
                    if quarantined {
                        " — moved aside; retry to download a fresh copy"
                    } else {
                        " — re-download it"
                    }
                )));
            }
        }
    }
    Ok(Some(bytes))
}

/// Base URL of the GitHub release that hosts the prebuilt ancestry/IBD assets for `build`. For
/// example:
/// `https://github.com/JamesKane/decodingus-navigator/releases/download/assets-chm13v2.0`.
/// `NAVIGATOR_ASSET_REPO` and `NAVIGATOR_ASSET_RELEASE` override it. The offline packager uses the
/// same two variables.
fn asset_release_base_url(build: ReferenceBuild) -> String {
    let repo = std::env::var("NAVIGATOR_ASSET_REPO").unwrap_or_else(|_| "JamesKane/decodingus-navigator".to_string());
    let tag = std::env::var("NAVIGATOR_ASSET_RELEASE").unwrap_or_else(|_| format!("assets-{}", build.as_str()));
    format!("https://github.com/{repo}/releases/download/{tag}")
}

/// The genetic-map asset path for a build (`$NAVIGATOR_GENETIC_MAP` override, else
/// `<base>/ancestry/genetic_map_<build>.bin`). Optional. Without it, IBD falls back to a uniform
/// map.
fn genetic_map_path(build: ReferenceBuild) -> PathBuf {
    ancestry_asset_path("NAVIGATOR_GENETIC_MAP", "genetic_map", build, "bin")
}

/// Load the real recombination map for IBD if the asset is present. If not, fall back to a
/// uniform 1 cM/Mb map over `lengths`, and log that. `lengths` is the observed `(chromosome, max_bp)` of
/// the compared samples. Only the uniform fallback uses it.
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

/// Map a computed [`AncestryResult`] onto the shared federated wire record. The record takes the
/// analysis method verbatim from the estimator that produced the result, and never infers it. So
/// the published `analysisMethod` always matches the composition on the screen.
/// How many outbox rows one [`App::drain_outbox`] pass tries.
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
    /// Whole-genome diploid variant calls (SNV and indel) for an alignment, as a VCF. It is heavy,
    /// because it walks the BAM again for each primary chromosome. The app caches the result.
    DiploidVcf(i64),
    /// The subject-level **consensus** diploid VCF: the joint genotype across the subject's
    /// alignments on the same build. It is heavy, with a call and a force-call for each alignment.
    ConsensusDiploidVcf(SampleGuid),
    /// The plain-language "DNA Story" brief for a subject, as a self-contained HTML document.
    SubjectBriefHtml(SampleGuid),
    /// The Y-DNA / mtDNA descent report as TSV: the root→terminal lineage with the call state of
    /// each SNP.
    DescentTsv(SampleGuid, DnaType),
}

impl ExportRequest {
    /// File extension (no dot) for the save dialog + filter.
    pub fn extension(&self) -> &'static str {
        match self {
            ExportRequest::CoverageHtml(_) | ExportRequest::AncestryHtml(_) | ExportRequest::SubjectBriefHtml(_) => {
                "html"
            }
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
            ExportRequest::SubjectBriefHtml(_) => "DNA story (HTML)",
            ExportRequest::DescentTsv(_, _) => "descent report (TSV)",
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
            ExportRequest::SubjectBriefHtml(_) => return format!("dna_story.{}", self.extension()),
            ExportRequest::DescentTsv(_, dna) => {
                let kind = match dna {
                    DnaType::Y => "y",
                    DnaType::Mt => "mt",
                };
                return format!("{kind}_descent.{}", self.extension());
            }
        };
        format!("{stem}_{id}.{}", self.extension())
    }
}

/// One row of the Y-STR concordance view. It holds a marker called from sequence, with its
/// FTDNA-convention value and its calibration status. Beside it are the subject's imported vendor
/// value and whether the two agree.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StrConcordanceRow {
    pub marker: String,
    /// Called FTDNA-convention value, or `None` if the marker was not called from sequence.
    pub called: Option<i32>,
    /// Calibration status: `Reliable` | `ConventionOffset` | `Excluded` | `Uncalibrated` | `NotCalled`.
    pub status: String,
    /// Whether the corpus calibrated the called value (Reliable/ConventionOffset), which makes it
    /// comparable.
    pub calibrated: bool,
    /// Imported vendor value (e.g. `"13"`, or a multi-copy `"11-15"`), or `None` if not in the profile.
    pub imported: Option<String>,
    pub depth: u32,
    /// Calibrated call whose value matches the imported single value.
    pub agree: bool,
}

/// Outcome of a batch [`App::add_data_batch`] run over more than one file or folder. It feeds the
/// import summary that the GUI shows after a multi-file Add Data or a drag-and-drop.
#[derive(Debug, Clone, Default)]
pub struct BatchImportSummary {
    /// `(filename, detected-type description)` for each successfully imported file.
    pub imported: Vec<(String, String)>,
    /// `(filename, reason)` for each file skipped or errored (unrecognized / import failure).
    pub skipped: Vec<(String, String)>,
}

/// Outcome of one staged **sample directory** that [`App::add_sample_dir`] ingested. This is the
/// CLI `ingest` fast path for the D2C bulk side-load. It records the alignment from the header,
/// then places Y/mt, sex, metrics, and coverage from the pipeline sidecars, with no CRAM decode.
#[derive(Debug, Clone, Default)]
pub struct SampleDirSummary {
    /// Alignment records created (BAM/CRAM newly recorded onto the subject's run).
    pub alignments_created: usize,
    /// Alignment files already recorded (idempotent re-ingest).
    pub alignments_skipped: usize,
    /// Variant files (VCFs) imported from the directory.
    pub variants_imported: usize,
    /// Whether the sidecar fast path ran (a haplogroup GVCF was present and attached).
    pub sidecars_ingested: bool,
    /// Y haplogroup placed from the chrY GVCF, if any.
    pub y_haplogroup: Option<String>,
    /// mt haplogroup placed from the chrM GVCF, if any.
    pub mt_haplogroup: Option<String>,
    /// Sex filled from the `.sex` sidecar (`M`/`F`/`U`), if any.
    pub sex: Option<String>,
    /// Read metrics filled from a `stats.txt`/flagstat/Picard sidecar.
    pub read_metrics: bool,
    /// Lite coverage filled from a `coverage.txt`/Picard WGS-metrics sidecar.
    pub lite_coverage: bool,
    /// Loose-bundle fallback (no alignment/variant/GVCF): one `(filename, description)` for each
    /// import.
    pub imported: Vec<(String, String)>,
    /// Loose-bundle fallback skips, plus any files that failed: `(filename, reason)`.
    pub skipped: Vec<(String, String)>,
    /// Non-fatal fast-path errors (a sidecar that failed while the rest proceeded).
    pub errors: Vec<String>,
}

/// A file's low-cost signature (`mtime_secs:size`) for analysis-cache staleness. It reads no
/// content. `None` if the file is missing, or if the code can not stat it.
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

/// Whether a cached artifact is still fresh for the current source signature.
///
/// A stale entry fails, which means the code knows both signatures and they are not equal. An unknown
/// stored sig passes, which is a legacy entry or a source that is not a file. An unknown current
/// sig also passes, because the file is gone and there is nothing to compute against.
fn artifact_is_fresh(stored: Option<&str>, current: Option<&str>) -> bool {
    match (stored, current) {
        (Some(s), Some(c)) => s == c,
        _ => true,
    }
}

/// A recognized data-file extension. This is the pre-filter for directory expansion, because the
/// code walks a dropped folder for these. `add_data` sniffs a text file (csv/tsv/txt) again, to
/// route chip, STR, or variant data.
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
        // CompleteGenomics masterVar dumps ship gzip/bzip2-compressed TSV.
        ".tsv.gz",
        ".tsv.bz2",
    ]
    .iter()
    .any(|e| n.ends_with(e))
}

/// Collect recognized data files from `path`: a file gives itself, if the code recognizes it, and
/// the code walks a directory recursively. The depth and the file count have bounds, so a large
/// tree can not run away.
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

/// The immediate child directories of `root` that gave at least one collected file. This is the
/// "this folder holds more than one sample" signal.
///
/// The folder of a single sample fans data into one subdirectory at most. FTDNA gives
/// `<sample>/<kit>/<uuid>.bam` plus a top-level results CSV, which is the kit directory alone. A
/// *parent* of many sample folders spreads the data across some of them.
///
/// A file directly in `root` does not count, because it belongs to the picked folder itself.
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

/// The result of one [`App::drain_outbox`] pass, which the UI reports and shows in its indicator.
#[derive(Debug, Clone, Default)]
pub struct DrainOutcome {
    /// `(kind, at-uri)` of each row published this pass.
    pub published: Vec<(String, String)>,
    /// Rows that hit a non-transient error, which the code then marked FAILED.
    pub failed: usize,
    /// Whether a transient failure rescheduled a row (i.e. we are likely offline).
    pub retry_scheduled: usize,
    /// Rows that still wait for a successful push after this pass.
    pub pending: i64,
}

/// The result of a [`App::pull_sync`] reconcile pass over the account's PDS records.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PullOutcome {
    /// Records unchanged since our last push.
    pub in_sync: usize,
    /// Records changed on the PDS and applied locally (where applicable).
    pub applied: usize,
    /// Remote records with no local mapping. These are PII-free summaries, and the app tracks them
    /// but does not rebuild them.
    pub adopted: usize,
    /// Locally-published records that the PDS does not have. The app flags them for re-publish.
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

/// Build a community feed-post record. This is the shared
/// `com.decodingus.atmosphere.feed.post` contract that the AppView mirrors into its feed, with a
/// top-level `createdAt` and optional topic or reply pointers. It holds no PII beyond the text the
/// user chose to publish. `reply` is the `(root_uri, parent_uri)` pair on a threaded reply, and
/// `None` for a top-level post.
fn feed_post_record(content: &str, topic: Option<&str>, reply: Option<(&str, &str)>) -> serde_json::Value {
    let mut rec = FeedPostRecord::new(content, Utc::now().to_rfc3339()).with_topic(topic.map(str::to_string));
    if let Some((root, parent)) = reply {
        rec = rec.with_reply(root, parent);
    }
    // A struct of plain strings always serializes. Show a build bug loudly, and do not drop the
    // post without a word.
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

/// Ordinal rank of an artifact's `completeness` for downgrade checks: `full` (2) beats `partial`
/// (1) beats anything else (0). Used by [`App::save_analysis_no_downgrade`] so a fast-path sidecar
/// result never overwrites an equal-or-fuller stored one.
fn completeness_rank(completeness: &str) -> u8 {
    match completeness {
        "full" => 2,
        "partial" => 1,
        _ => 0,
    }
}

/// Deterministic PDS record key for a subject's biosample record.
///
/// The key comes from the subject GUID, which is a random local UUID and holds no donor PII. So
/// the code knows the record's at:// URI *before* it publishes the record. Every child record
/// (coverage, ancestry, sequence-run) can then reference it, and a second publish overwrites in
/// place instead of a duplicate. The device-key record uses the same deterministic-rkey pattern.
fn biosample_rkey(guid: SampleGuid) -> String {
    format!("bio-{}", guid.0.simple())
}

/// Deterministic PDS record key for a sequence-run record (stable within the account's repo).
fn seqrun_rkey(run_id: i64) -> String {
    format!("run-{run_id}")
}

/// Deterministic PDS record key for an alignment (coverage-summary) record.
///
/// The key is stable for each alignment, so a second publish overwrites in place and makes no
/// duplicate. Without it the record took the `create` path, with a fresh PDS-assigned TID each
/// time. Two concurrent drains of the same alignment then raced into two records, which is the
/// federated "2 samples" duplicate. This mirrors the deterministic rkeys that the biosample and
/// sequence-run records already use.
fn alignment_rkey(alignment_id: i64) -> String {
    format!("aln-{alignment_id}")
}

/// The at:// URI that a published biosample record has in `did`'s repo. Child records link to this
/// anchor.
fn biosample_at_uri(did: &str, guid: SampleGuid) -> String {
    format!("at://{did}/{NS_BIOSAMPLE}/{}", biosample_rkey(guid))
}

/// What an ancestral-origin batch did. `considered` is every row that the consent predicate
/// allowed. `refused` is how many of those a field gate then rejected, which is a normal outcome.
/// It is the number to watch, because a jump in it means the MDKA data changed shape.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize)]
pub struct OriginPublishReport {
    pub considered: usize,
    pub publishable: usize,
    pub refused: usize,
    /// Of the publishable, how many carry a place (an ancestor with a birth year) …
    pub with_place: usize,
    /// … and how many the precision ladder reduced to their country.
    pub country_only: usize,
}

/// Deterministic PDS record key for one lineage's ancestral-origin record. It is stable for each
/// (subject, lineage) pair, so a corrected MDKA overwrites and makes no duplicate ancestor.
fn origin_rkey(guid: SampleGuid, lineage: &str) -> String {
    format!("origin-{}-{}", lineage.to_lowercase(), guid.0.simple())
}

/// The at:// URI a published sequence-run record has in `did`'s repo.
fn seqrun_at_uri(did: &str, run_id: i64) -> String {
    format!("at://{did}/{NS_SEQUENCERUN}/{}", seqrun_rkey(run_id))
}

/// Fallback reference build for a batch import when neither the header nor the filename
/// identifies one. This app's analysis reference is CHM13v2.0, so the code binds an unlabeled file
/// to it and does not leave it unresolved. The project folders on the NAS are CHM13v2.0.
const DEFAULT_IMPORT_BUILD: &str = "chm13v2.0";

/// Detect an alignment's reference build for batch import. The function reads the BAM/CRAM
/// **header** and no more, which is fast and needs no reference FASTA. It prefers the `@SQ`/`@PG`
/// signal, and falls back to the filename heuristic. It returns `"unknown"` only when both are
/// silent. The IO blocks, so call this inside `spawn_blocking`. It returns `(build, source)`, where
/// `source` names how the code decided the build, for the import diagnostics.
fn detect_build_for(path: &Path) -> (String, &'static str) {
    match navigator_analysis::probe::probe_alignment(path) {
        Ok(probe) => match probe.reference_build {
            Some(b) => (b, "header probe"),
            None => (reference_build_for(path), "filename"),
        },
        Err(e) => {
            eprintln!(
                "project import: header probe failed for {} ({e}); falling back to the filename",
                path.display()
            );
            (reference_build_for(path), "filename (probe failed)")
        }
    }
}

/// A fast look at a VCF header: the joined `##` meta block, plus the contig names from
/// `##contig=<ID=…>`. It reads only the header, and stops at the first data line. It decompresses
/// a gzip/BGZF (`.vcf.gz`) input, the same as the import parser, so it also classifies a bgzipped
/// vendor VCF.
fn peek_vcf_header(path: &Path) -> (String, Vec<String>) {
    use std::io::BufRead;
    let Ok(reader) = navigator_analysis::gzio::open_maybe_gz(path) else {
        return (String::new(), Vec::new());
    };
    let mut meta = String::new();
    let mut contigs = Vec::new();
    for line in reader.lines().map_while(Result::ok) {
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

/// Parse a VCF into the **subject's** SNP calls, and obey the genotype.
///
/// A vendor VCF (FTDNA Big Y, YSEQ, …) lists a site's REF/ALT even where the sample is
/// homozygous-reference (`GT 0/0`). An example is `chrY 2781955 C T … 0/0`, where the sample is C
/// and not T. A parser that takes `ALT[0]` blindly, as a sites-only parser does, records that T as
/// a derived call. A Big Y export carries thousands of such reference sites, and the placement then
/// goes deeper into branches that the sample does not carry.
///
/// So when the VCF has a genotyped sample column, the parser reads its `GT`. It keeps a
/// single-base ALT only when the genotype selects it, which is the first non-zero allele, and it
/// handles a multi-allelic row. It drops a `0/0` row and a `./.` row.
///
/// A VCF with no FORMAT or sample column is a sites-only list. It keeps its old sense: every
/// listed ALT is one of the subject's variants.
fn parse_vcf_subject_snps(path: &Path) -> Result<Vec<variants::VariantCall>, AppError> {
    use std::io::BufRead;
    let reader = navigator_analysis::gzio::open_maybe_gz(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
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
        let sample_field = |key: &str| -> Option<&str> {
            (f.len() >= 10)
                .then(|| {
                    f[8].split(':')
                        .position(|k| k == key)
                        .and_then(|i| f[9].split(':').nth(i))
                })
                .flatten()
                .filter(|v| !v.is_empty() && *v != ".")
        };
        let gt = sample_field("GT");
        // `alt_index` is the ALT that the genotype selected. The code needs it to pick the correct
        // AD entry on a multi-allelic row, where AD is [ref, alt1, alt2, …].
        let (alt, genotype, alt_index) = match gt {
            Some(gt) => {
                // The first non-zero allele index selects the carried ALT. All-zero (0/0) or
                // no-call (./.) means the subject is reference here, so skip it.
                match gt
                    .split(['/', '|'])
                    .filter_map(|a| a.parse::<usize>().ok())
                    .find(|&a| a > 0)
                {
                    Some(idx) => match alts.get(idx - 1) {
                        Some(&a) => (a, Some(gt.to_string()), idx),
                        None => continue,
                    },
                    None => continue,
                }
            }
            None => (alts[0], None, 1), // sites-only VCF: the listed variant is the subject's
        };

        // Evidence the source supplies. Every field stays `None` when it is absent. A missing DP
        // means "the vendor did not say", and a stored 0 would make a good call look unsupported.
        let ad: Option<Vec<u32>> = sample_field("AD").map(|v| v.split(',').map(|x| x.parse().unwrap_or(0)).collect());
        let evidence = variants::CallEvidence {
            qual: f.get(5).and_then(|q| q.parse::<f64>().ok()),
            // The code stores only failures. `.` and `PASS` are almost all of the rows, and they
            // mean nothing.
            filter: f
                .get(6)
                .filter(|v| !v.is_empty() && **v != "." && !v.eq_ignore_ascii_case("PASS"))
                .map(|v| v.to_string()),
            dp: sample_field("DP").and_then(|v| v.parse().ok()),
            gq: sample_field("GQ").and_then(|v| v.parse().ok()),
            ad_ref: ad.as_ref().and_then(|a| a.first().copied()),
            ad_alt: ad.as_ref().and_then(|a| a.get(alt_index).copied()),
        };
        if let Some(call) =
            variants::snp_call_with_evidence(chrom, pos, reference, alt, Some(id.to_string()), genotype, evidence)
        {
            out.push(call);
        }
    }
    Ok(out)
}

/// The reference build of an external autosomal gVCF, as a token that the IBD panel's `locus()`
/// accepts (`GRCh38` / `GRCh37` / `chm13`). `NAVIGATOR_CALLSET_BUILD` wins, else the VCF `##` meta
/// (`detect_vcf_build`), else the `chr1`/`1` contig length in the header, else GRCh38, which is the
/// GATK4 WGS default. The code normalizes the value, so a `chm13v2.0`, `hg38`, or `b37` form
/// still resolves.
fn callset_build_for(path: &Path) -> String {
    let raw = std::env::var("NAVIGATOR_CALLSET_BUILD")
        .ok()
        .filter(|s| !s.trim().is_empty())
        .or_else(|| {
            let (meta, _) = peek_vcf_header(path);
            detect_vcf_build(&meta).or_else(|| {
                // Fall back to the chr1/1 contig length.
                meta.lines().find_map(|line| {
                    let after = line.strip_prefix("##contig=<ID=")?;
                    let id: String = after.chars().take_while(|&c| c != ',' && c != '>').collect();
                    if id != "chr1" && id != "1" {
                        return None;
                    }
                    let lp = line.find("length=")?;
                    let len: String = line[lp + 7..].chars().take_while(|c| c.is_ascii_digit()).collect();
                    Some(match len.as_str() {
                        "249250621" => "GRCh37".to_string(),
                        "248387328" => "chm13".to_string(),
                        _ => "GRCh38".to_string(),
                    })
                })
            })
        })
        .unwrap_or_else(|| "GRCh38".to_string());
    let l = raw.to_ascii_lowercase();
    if l.contains("chm13") || l.contains("t2t") || l.contains("hs1") {
        "chm13".to_string()
    } else if l.contains("38") {
        "GRCh38".to_string()
    } else if l.contains("37") || l.contains("19") || l.contains("b37") {
        "GRCh37".to_string()
    } else {
        "GRCh38".to_string()
    }
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

/// Read a `readme.txt` in the same directory, if it is present. An FTDNA Big Y bundle puts one
/// beside `variants.vcf`.
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

/// A label context that makes a vendor VCF distinct. It is the parent directory when the file name
/// is the generic vendor name (`variants.vcf`), and the file name itself in all other cases.
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

/// SHA-256 of a file's content, in hex, computed off the async runtime. It streams the file, which
/// is safe for a large alignment.
async fn sha256_file_async(path: PathBuf) -> Result<String, AppError> {
    let hash = tokio::task::spawn_blocking(move || du_bio::hash::sha256_file(&path)).await??;
    Ok(hash)
}

/// SHA-256 of an in-memory string, in hex, for a tree JSON or other small content.
fn sha256_str(s: &str) -> String {
    du_bio::hash::sha256_hex(s.as_bytes())
}

/// Rebuild a small [`HaploAssignment`] from a recorded call: the terminal and the lineage, with no
/// full ranked list and no branch evidence. A cache hit on the score returns this. The recorded
/// call is the source of truth, and only a fresh score needs the detail.
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

// Watson–Crick complement of a base, for a reverse-strand lift. This is the shared helper in
// navigator-domain, which the chip dosage and BISDNA QC paths also use.
use navigator_domain::seq::complement_base;

/// The build that a haplotree's positions are in, by contig. The FTDNA Y tree is GRCh38. mtDNA
/// (`chrM`) is rCRS and stays a direct query with no chain, so it returns `None`.
/// Whether a stored reference-build string means GRCh38, which is the FTDNA Y tree's native
/// coordinate space. `None` means GRCh38, the vendor-Y-VCF import default. The FTDNA-provider Y
/// consensus uses this to admit only GRCh38 vendor sets, because others do not match the GRCh38
/// tree positions.
fn is_grch38_build(build: &Option<String>) -> bool {
    match build {
        None => true,
        Some(b) => {
            let b = b.to_ascii_lowercase();
            b.contains("grch38") || b.contains("hg38") || b == "38" || b == "b38"
        }
    }
}

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

/// Whether a trusted external caller wins over Navigator's own genotyping. Such a caller is a GATK4
/// GVCF or a 1240K call set that the sidecar fast path imported.
///
/// When this is on, an external haplogroup call wins the reconcile, and Navigator's internal caller
/// does not walk that alignment again. The env var `NAVIGATOR_PREFER_EXTERNAL_CALLS` wins, then the
/// settings, then the **default of on**. Never dilute a call that the user deliberately produced.
///
/// This is a pure resolver, split out so a test can drive it.
fn resolve_prefer_external(env: Option<&str>, settings: Option<bool>) -> bool {
    match env.map(str::trim) {
        Some(v) if v.eq_ignore_ascii_case("false") || v == "0" || v.eq_ignore_ascii_case("off") => false,
        Some(v) if v.eq_ignore_ascii_case("true") || v == "1" || v.eq_ignore_ascii_case("on") => true,
        _ => settings.unwrap_or(true),
    }
}

pub(crate) fn prefer_external_calls() -> bool {
    let env = std::env::var("NAVIGATOR_PREFER_EXTERNAL_CALLS").ok();
    resolve_prefer_external(env.as_deref(), AppSettings::load().prefer_external_calls)
}

/// The painter's built-in knob defaults, as plain scalars, for the Settings UI to show when a knob
/// has no value. The UI does not depend on `navigator-analysis`. A second copy of the literals
/// there would let the two sets become different. That is exactly how a stale calibration goes
/// back into the settings.
pub struct LaiKnobDefaults {
    pub recomb_per_cm: f64,
    pub max_ref_haps: u32,
    pub min_ancestry: f64,
    pub switch_per_cm: f64,
    pub min_segment_cm: f64,
    pub size_normalize: f64,
    pub mismatch: f64,
}

/// [`CopyingLaiParams::default`], flattened for the UI.
pub fn lai_knob_defaults() -> LaiKnobDefaults {
    let d = navigator_analysis::lai::CopyingLaiParams::default();
    LaiKnobDefaults {
        recomb_per_cm: d.recomb_per_cm,
        max_ref_haps: d.max_ref_haps as u32,
        min_ancestry: d.min_ancestry,
        switch_per_cm: d.switch_per_cm,
        min_segment_cm: d.min_segment_cm,
        size_normalize: d.size_normalize,
        mismatch: d.mismatch,
    }
}

/// The copying-LAI chromosome-painter parameters, resolved from [`AppSettings`] with the built-in
/// [`CopyingLaiParams::default`] for any knob with no value. The code reads them at paint time, so
/// an edit applies immediately.
pub(crate) fn copying_lai_params() -> navigator_analysis::lai::CopyingLaiParams {
    let s = AppSettings::load();
    let d = navigator_analysis::lai::CopyingLaiParams::default();
    navigator_analysis::lai::CopyingLaiParams {
        recomb_per_cm: s.lai_recomb_per_cm.unwrap_or(d.recomb_per_cm),
        max_ref_haps: s.lai_max_ref_haps.map(|v| v as usize).unwrap_or(d.max_ref_haps),
        min_ancestry: s.lai_min_ancestry.unwrap_or(d.min_ancestry),
        switch_per_cm: s.lai_switch_per_cm.unwrap_or(d.switch_per_cm),
        min_segment_cm: s.lai_min_segment_cm.unwrap_or(d.min_segment_cm),
        size_normalize: s.lai_size_normalize.unwrap_or(d.size_normalize),
        mismatch: s.lai_mismatch.unwrap_or(d.mismatch),
        min_ref_haps: d.min_ref_haps,
    }
}

/// `haplogroup_call.source_key` for an alignment's **external** (sidecar fast-path) Y call. An
/// external call and an internal (CRAM-walk) call use different keys, and the walk uses
/// `aln:{id}` or `aln:{id}:mt`. So neither upsert can ever overwrite the other.
pub(crate) fn external_y_source_key(alignment_id: i64) -> String {
    format!("aln:{alignment_id}:ext")
}

/// `haplogroup_call.source_key` for an alignment's **external** (sidecar fast-path) mtDNA call.
pub(crate) fn external_mt_source_key(alignment_id: i64) -> String {
    format!("aln:{alignment_id}:ext:mt")
}

/// One DNA type's caller comparison for an alignment: the trusted external terminal, from an
/// imported GVCF, against Navigator's own internal-walk terminal. [`App::compare_callers`] produces
/// it. That is the "Compare callers" diagnostic. It shows a difference between GATK and
/// Navigator, such as ancient-DNA damage, even when the external caller is the preferred one.
#[derive(Debug, Clone)]
pub struct CallerComparison {
    pub dna_type: DnaType,
    /// Terminal from the trusted external caller (sidecar GVCF), if the store has one.
    pub external: Option<String>,
    /// Terminal from Navigator's own genotyping. The code forces it, whatever the prefer-external
    /// policy says.
    pub navigator: Option<String>,
}

impl CallerComparison {
    /// True when both callers produced a terminal and the two match. `false` when one is missing,
    /// or when the two are different. The second case is the one to show.
    pub fn agree(&self) -> bool {
        matches!((&self.external, &self.navigator), (Some(a), Some(b)) if a == b)
    }
}

/// Application interface mode. `Simple` is a casual, single-person experience with plain-language
/// briefs. `Advanced` is the full power-user UI, with projects and analysis of each source. This is
/// app-level UI state, held in [`AppSettings`].
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
    // A recognized env value wins. The code ignores a value it does not recognize, and falls
    // through to the settings.
    env.and_then(UiMode::parse).or_else(|| settings.and_then(UiMode::parse))
}

/// The configured UI mode. `NAVIGATOR_UI_MODE` wins over the stored setting. `None` when neither
/// has a value, which is the first run, where the UI takes a default from a workspace heuristic.
pub fn configured_ui_mode() -> Option<UiMode> {
    let env = std::env::var("NAVIGATOR_UI_MODE").ok();
    let settings = AppSettings::load().ui_mode;
    resolve_ui_mode(env.as_deref(), settings.as_deref())
}

/// Store the chosen UI mode, and keep all other settings.
pub fn persist_ui_mode(mode: UiMode) -> std::io::Result<()> {
    let mut s = AppSettings::load();
    s.ui_mode = Some(mode.as_str().to_string());
    s.save()
}

/// Production DecodingUs AppView (the `/api/v1/*` backend: trees, lab lookup, IBD, exchange,
/// social). This is the default when neither `DECODINGUS_APPVIEW_URL` nor the `appview_url` setting
/// has a value. Override either one for a local dev backend.
const DEFAULT_APPVIEW_URL: &str = "https://decoding-us.org";

/// Resolve the AppView base URL. The function is pure. The env var wins, then the settings, then
/// the production default. It trims a slash at the end, and ignores a blank value.
fn resolve_appview_url(env: Option<String>, settings: Option<String>) -> String {
    env.or(settings)
        .map(|s| s.trim_end_matches('/').to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| DEFAULT_APPVIEW_URL.to_string())
}

fn decodingus_appview_url() -> String {
    resolve_appview_url(
        std::env::var("DECODINGUS_APPVIEW_URL").ok(),
        AppSettings::load().appview_url,
    )
}

/// Production OAuth client identity: the URL of Navigator's **native/public** client-metadata
/// document at decoding-us.org (loopback redirect, PKCE, `token_endpoint_auth_method: none`). It
/// is different from the confidential *web* client at `/oauth/client-metadata.json`, because a
/// desktop app can not hold a signing key or receive a server-side redirect.
const DEFAULT_OAUTH_CLIENT_ID: &str = "https://decoding-us.org/oauth/navigator-client-metadata.json";

/// OAuth scope Navigator requests: identity (`atproto`) **plus** transitional generic write access
/// (`transition:generic`). A publish of federated records to the user's PDS needs write scope, and
/// not identity alone. This must match the `scope` in the published client-metadata document.
const OAUTH_SCOPE: &str = "atproto transition:generic";

/// Resolve Navigator's OAuth client config. The function is pure.
///
/// `DECODINGUS_OAUTH_CLIENT_ID` overrides the hosted default. The literal `loopback` selects the
/// atproto dev loopback client, for a login against a local or test PDS that has not registered
/// the production document. The code reads any other non-blank value as a hosted client-metadata
/// URL.
fn resolve_oauth_config(env_client_id: Option<String>) -> OAuthConfig {
    match env_client_id.map(|s| s.trim().to_string()).filter(|s| !s.is_empty()) {
        Some(v) if v.eq_ignore_ascii_case("loopback") => OAuthConfig::loopback(OAUTH_SCOPE),
        Some(url) => OAuthConfig::hosted(url, OAUTH_SCOPE),
        None => OAuthConfig::hosted(DEFAULT_OAUTH_CLIENT_ID, OAUTH_SCOPE),
    }
}

fn decodingus_oauth_config() -> OAuthConfig {
    resolve_oauth_config(std::env::var("DECODINGUS_OAUTH_CLIENT_ID").ok())
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
    /// Provenance for each source: which tests contributed, and how many variants each one gave.
    #[serde(default)]
    pub sources: Vec<YSourceSummary>,
}

/// The Y-DNA view of a [`ConsensusProfile`]. The Y adapter is the first consumer of the generic
/// consensus aggregate. The name stays as it is for the Y-DNA tab and the worker contract.
pub type YProfile = ConsensusProfile;

/// One source that contributes to a [`ConsensusProfile`], for the provenance display.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct YSourceSummary {
    pub label: String,
    pub source_type: SourceType,
    pub variant_count: usize,
}

/// A subject's **autosomal** multi-source consensus profile. It is the diploid (0/1/2) equivalent
/// of [`ConsensusProfile`], over the canonical CHM13 IBD-panel sites. It goes in the same
/// `consensus_profile` table, under `dna_type='Auto'`. It has no lineage label, because an autosome
/// has no haplogroup.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DiploidProfile {
    pub variants: Vec<DiploidVariant>,
    pub summary: YProfileSummary,
    /// Provenance for each source: which tests contributed, and how many sites each one gave.
    #[serde(default)]
    pub sources: Vec<YSourceSummary>,
}

/// `ancestry_result.alignment_id` sentinel for a result that comes from the subject's autosomal
/// **consensus**, which pools all sources, and not from one sequencing alignment.
pub const CONSENSUS_SOURCE_ID: i64 = 0;

/// One row of the deep-ancestry stability diagnostic ([`App::ancient_ancestry_stability`]). It is
/// the ancient mixture fitted over one *view* of a subject: the pooled consensus, one source alone,
/// or a thinned site set. It is a diagnostic only. The app never stores it and never publishes it.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct AncientFitRow {
    /// Which view of the subject this fit came from.
    pub label: String,
    /// Panel sites with a genotype in this view.
    pub sites: usize,
    /// `(population_code, percentage)`, ranked.
    pub components: Vec<(String, f64)>,
    /// Model dispersion (≈1 at the noise floor; grows as the sample leaves the sources' span).
    pub dispersion: f64,
    /// European share (%) by the modern super-pop estimate, which is the deep model's scope check.
    pub european: f64,
    /// Whether the released estimator would report this fit, or suppress it as inapplicable.
    pub reported: bool,
}

/// Ancient (deep) ancestry. **Rebuilt, but still OFF: it fails the stability gate.**
///
/// The original PCA-centroid implementation made up numbers (§1–2), and the team disabled it. Some
/// rebuild tries then looked like failures: a frequency-mixture EM, ascertainment, and
/// pseudo-haploid. Those "walls" were in fact an outgroup mistake, 16k-SNP imprecision, and a
/// genotype-label bug in a diagnostic tool (docs §7.9–7.13).
///
/// The method that works, and that is **now enabled**, is the qpAdm f4 estimator
/// ([`navigator_analysis::ancestry::estimate_qpadm_ancestry`]) over the full-1240k
/// **Patterson-2022 sister-outgroup panel** ([`ancestry_qpadm_path`]). The sources are WHG, EEF,
/// and Steppe. Each one has a sister outgroup, which gives the fit the power to separate them
/// (§7.14). [`App::estimate_deep_ancestry`] computes it, and genotypes the subject's CHM13
/// alignments at about 1.15M sites.
///
/// It passes the **stability** gate of §5.4, the one most diagnostic test, which every earlier
/// try failed. Subject `huF98AFD` gives **WHG 14.6 / EEF 44.8 / Steppe 40.6 from his WGS**. His
/// 23andMe chip gives **WHG 14.3 / EEF 44.9 / Steppe 40.8**. That is the same person by two
/// independent means, which agree to 0.3% or better, and the fit accepts both models (±1–2% SE). A
/// cross-check against real `admixtools2` agrees to 0.2%, over a 99.84% same-person genotype
/// concordance. This is a literature-grade British breakdown.
///
/// The gate covers the computation, the display, **and** the publish. The estimator returns `None`,
/// and stores nothing, for any sample outside the model's scope or with a rejected fit. So an
/// inapplicable breakdown never reaches the UI or the PDS.
///
/// See `documents/design/ancient-ancestry-rebuild.md` (start at §7.14).
pub const ANCIENT_ANCESTRY_ENABLED: bool = true;

/// Whether the app computes, reads back, or shows Tier B **archaic segments**, which are the
/// introgressed-tract caller and its chromosome browser.
///
/// **ON**, for the rebuilt caller ([`navigator_analysis::archaic_match`]), and specifically as a
/// **within-population** measure. The history below is the first implementation, which the team
/// withdrew.
///
/// What changed: the caller no longer counts private-variant density, and it now matches the
/// archaic genomes directly. A hold-out of 30 Europeans that the fit never saw gives
/// **r = +0.710 (p < 0.0001)** for the extent in each individual. The old caller reached −0.018
/// (p = 0.94) on the same test. Every individual of 90 scores above their own random-placement
/// null. A concordance filter takes precision from 54 % to 90 %.
///
/// **You must not use it to compare people of different ancestries.**
///
/// East Asian tracts match our four sequenced archaic genomes less well than European tracts do.
/// So the extent is under-called for them. The reported figure then orders the two populations in the wrong direction against the
/// truth. That is a property of which archaic genomes have a sequence, and not a threshold. So the
/// UI states the limit, and does not suggest a universal percentage.
pub const ARCHAIC_SEGMENTS_ENABLED: bool = true;

// The rationale of the earlier gate, kept because it is the reason the current caller exists.
//
// **Off**: a check against an external truth set for each individual showed that the caller
// carries no signal for one person.
//
// Tier B shipped on the strength of one number. Its total extent landed at 1.01x the hmmix
// European mean. A proper check, against hmmix's own calls **for the same individuals** (n=20
// Europeans, chr21+22), showed that this number is all there is.
//
// **The locations disagree.** For HG00096, we also call 2.1 % of hmmix's archaic bases, against
// an expectation of 5.0 % (p95 9.4 %). That expectation is for segments of our own lengths at
// *random* positions in the same span. So the result is below chance.
//
// It is not a coordinate artefact. The curve of overlap against shift is flat across +/-2 Mb with
// no peak. And 70.7 % of the truth lies inside our callable territory, so the tracts were
// reachable.
//
// **The amounts do not track the individual.** Across the 20, Pearson r = -0.018 (p = 0.94) and
// Spearman rho = -0.020, against a true range of 1.19-2.97 Mb. Our own spread is 0.63x the spread
// of the truth. The two individuals with the least archaic ancestry drew our two highest calls.
//
// The mean ratio is about 0.92. The caller reproduces the cohort average and nothing else, which
// is exactly what its three fitted parameters make it do. An honest report needs a measurement of
// *this person*. So the app withholds the feature, and does not show it with a caveat.
//
// The machinery stays behind this flag, with its unit tests. This is the same discipline as
// `attribute_lineage` and the ancient-ancestry precedent above. To enable it again needs a change
// of method, not a change of threshold. The design records Skov-2020 haplotype matching as the
// path. It also needs a re-run of the check harness that produced these numbers.
//
// Tier A is the marker **count** and percentile. It is a different method on a different asset,
// and this flag does **not** gate it.
//
// See `documents/design/ArchaicAncestry_Design.md`, "Tier B validation (2026-07-30)" and
// "Why it failed" for the full record. Historical only: the live gate is the `true` above.

/// The stored method name of the deep-ancestry breakdown. It is re-exported, so the UI reads the
/// rebuilt method by name and can never fall back to a retired one.
pub use navigator_analysis::ancestry::ANCIENT_ADMIXTURE;

/// Bridge the autosomal consensus to the genotype carrier that the ancestry estimators and the IBD
/// detector consume. Each reconciled site becomes a [`SiteGenotype`] with the consensus dosage,
/// which is the count of the CHM13 ALT. That is the canonical orientation, and the AIM freq and PCA
/// assets use it as their key. A no-call has dosage -1, it passes through, and the code downstream
/// ignores it like any missing genotype.
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

/// Flatten a placement's branch SNP evidence into one observation for each SNP. The name dedupes
/// them: a SNP defines one branch, but guard against a duplicate. `in_tree` is true for a SNP that
/// defines a tree node.
///
/// Build one observation for each SNP of the multi-source consensus profile, from a placement's
/// **lineage**. That lineage holds the root→terminal mutations that the sample carries. Do not use
/// its child branches, which are the deeper splits that the descent did not take, and which are
/// ancestral or no-call by construction. With `branches`, a single-source profile read as
/// all-no-call even when the terminal placed cleanly.
fn snp_obs_from_assignment(assignment: &HaploAssignment, in_tree: bool) -> Vec<YObsInput> {
    let mut by_name: std::collections::HashMap<String, YObsInput> = std::collections::HashMap::new();
    for snp in &assignment.lineage {
        by_name.entry(snp.name.clone()).or_insert_with(|| {
            // Carry the observed base, and not the state alone. The code imputes the state from
            // the base. It keeps the base so that the profile can take a new imputation against a
            // corrected polarity, with no second genotype pass (see `consensus::reproject`).
            YObsInput::observed(
                snp.name.clone(),
                snp.position,
                snp.ancestral.clone(),
                snp.derived.clone(),
                snp.base,
                in_tree,
            )
        });
    }
    by_name.into_values().collect()
}

/// Reference-genome status and override for each build, for the Settings UI.
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
mod appview;
pub use analysis::AnalysisStep;
mod auth;
mod blocktree;
pub use blocktree::COLLAPSE_MIN_RUN;
mod brief;
mod commands;
mod dm;
mod fastpath;
mod ftdna_import;
mod haplogroup;
mod ibd_exchange;
mod import_profiles;
mod import_unified;
pub mod llm;
pub use llm::{ChatTurn, NarratedBrief};
pub use navigator_domain::results_context::SignalKind;
mod maintenance;
mod matching;
mod publish;
mod queries;
mod realign;
/// This is re-exported alone, and the module stays closed. The UI must name the build that it
/// offers to realign to, and nothing else in the module is its business.
pub use realign::{is_target_build, DEFAULT_TARGET_BUILD};
pub mod realign_job;
mod recruitment;
mod social;
mod sync;
pub mod update;

impl App {
    /// Reference-genome settings and cache status, with one row for each supported build.
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

    /// Set the local-FASTA override and the auto-download flag for a build, and store
    /// `reference_sources.json`. It applies on the next reference resolve. No restart is necessary,
    /// because the gateway reads the file again.
    pub fn set_reference_override(
        &self,
        build: &str,
        local_path: Option<String>,
        auto_download: bool,
    ) -> Result<(), AppError> {
        self.set_reference_overrides(&[ReferenceOverrideInput {
            build: build.to_string(),
            local_path,
            auto_download,
        }])
    }

    /// Store **all** reference-source overrides in ONE load-change-save.
    ///
    /// The Settings "References" table has one row for each build. The old code sent a separate
    /// `SetReferenceOverride` command for each row. Every worker command goes through
    /// `tokio::spawn`, so those N concurrent load-change-`save` cycles raced the shared
    /// `reference_sources.json`. That lost rows, *and* tore the file into corrupt JSON, with the
    /// head of one write and the tail of another (issue #26).
    ///
    /// This applies every row against a single load, and writes once, atomically, through
    /// [`UserConfig::save`]. That removes both failure modes.
    pub fn set_reference_overrides(&self, rows: &[ReferenceOverrideInput]) -> Result<(), AppError> {
        let path = self.gateway.config_path();
        let mut cfg = navigator_refgenome::UserConfig::load(&path);
        for row in rows {
            let key = canonical_build(&row.build)
                .map(|b| b.as_str().to_string())
                .unwrap_or_else(|| row.build.clone());
            let entry = cfg.references.entry(key).or_default();
            entry.local_path = row
                .local_path
                .as_ref()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty());
            entry.auto_download = row.auto_download;
        }
        cfg.save(&path)?;
        Ok(())
    }
}

/// One reference-source override to store. It is a row of the Settings "References" table, with a
/// build's local-FASTA path and its auto-download flag. [`App::set_reference_overrides`] batches
/// them.
#[derive(Debug, Clone)]
pub struct ReferenceOverrideInput {
    pub build: String,
    pub local_path: Option<String>,
    pub auto_download: bool,
}

/// One instrument→lab association from the AppView `sequencer` endpoints (D8). It mirrors the
/// `SequencerLabDto` shape, and it accepts extra fields.
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

/// The DecodingUs tree's own coordinate space, and the only one that carries almost all of it.
/// The hs1 build has coordinates for 99.8% of the tree's ~204k variants, against 86.5% for
/// GRCh38. The
/// reason is that the project called most DecodingUs-discovered (`DU`-named) SNPs in CHM13. It
/// mapped only a few hundred back to the older references.
///
/// Parse under this build whenever the code joins the tree to data by SNP **name**. There a locus
/// position is for display, and not for lookup. Any narrower build drops the variants it lacks, and
/// says nothing. Placement is the opposite case. It queries an alignment by position, so it must
/// use *that* alignment's build ([`decodingus_build_key`]).
pub(crate) const DECODINGUS_NATIVE_BUILD: &str = "hs1";

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

/// Whether an alignment's reference build matches a GVCF name's build token, such as `chm13`. The
/// comparison uses the canonical build, so `chm13`, `chm13v2`, and `hs1` all agree. A token that
/// does not resolve to a known build counts as a non-match, and the code falls back to the first
/// alignment.
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

/// Read the first 64 KiB of a file as lossy UTF-8. That is enough to find a text file's type,
/// and it does not read a multi-MB chip export into memory.
fn read_head(path: &Path) -> Result<String, AppError> {
    use std::io::Read;
    // Transparently decompress gzip/BGZF/bzip2 so the fingerprint sees text even for a compressed
    // dump (e.g. a CompleteGenomics `var-*.tsv.bz2`). Plain files read straight through.
    let mut reader = navigator_analysis::gzio::open_maybe_compressed(path)?;
    let mut buf = vec![0u8; 64 * 1024];
    let mut filled = 0;
    // Decoders satisfy a read in small chunks; loop until the head buffer is full or EOF.
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    buf.truncate(filled);
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Group the genotypes of each site into a dosage array for each chromosome, sorted by position,
/// for the IBD detector.
/// Artifact kind for an alignment's cached IBD-panel genotypes.
///
/// The `2` suffix retires every genotype cached before GRCh37/GRCh38 support landed (`3cf4956`,
/// 2026-07-13). Before that, the code genotyped those builds at the *CHM13* panel coordinates on a
/// BAM that is not CHM13. The positions were wrong, so the calls were near-random and never
/// heterozygous. The code fix did not change the cache key, which is the panel-manifest salt plus
/// `GENOTYPE_VERSION`. So the corrupt dosages continued to feed the autosomal consensus (IBD,
/// ancestry, identity). A new stem forces a one-time re-genotype with the build-aware path.
const IBD_PANEL_KIND: &str = "ibd_panel_genotypes2";

/// The IBD-panel genotype cache kind, salted with the CHM13 panel asset's manifest sha256 (the
/// first 16 hex characters). So a new build of the panel invalidates a stale genotype of an
/// alignment. It falls back to the bare kind when the build publishes no manifest, and
/// `GENOTYPE_VERSION` is then the only key.
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

/// Cache kind for the archaic-panel genotypes of one alignment.
const ARCHAIC_PANEL_KIND: &str = "archaic_panel_genotypes";

/// The archaic-panel genotype cache kind, salted with the panel asset's manifest sha256, exactly as
/// [`ibd_panel_cache_kind`] does it. The archaic panel's site list changes at every recalibration
/// of its thresholds. Genotypes taken over an older site set would corrupt the count, with no
/// warning.
fn archaic_panel_cache_kind() -> String {
    let build = ReferenceBuild::Chm13v2;
    let name = archaic_markers_path(build)
        .file_name()
        .and_then(|n| n.to_str())
        .map(String::from);
    match (load_asset_manifest(build), name) {
        (Some(m), Some(n)) => match m.assets.get(&n) {
            Some(e) => format!("{ARCHAIC_PANEL_KIND}:{}", &e.sha256[..16.min(e.sha256.len())]),
            None => ARCHAIC_PANEL_KIND.to_string(),
        },
        _ => ARCHAIC_PANEL_KIND.to_string(),
    }
}

/// Count of the sites called, with a dosage within the ploidy, in **both** samples. This is the
/// effective size of the IBD comparison. The app shows it, so that a sparse chip↔chip or chip↔WGS
/// overlap does not look like a confident result.
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
/// the count of the shared sites. The alignment-pair path and the chip-or-WGS compare path both
/// use this.
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

/// IBD detection over two [`IbdSite`] dosage vectors. This is the federated-exchange path, where
/// the partner's dosages arrive as `IbdSite` and not as [`SiteGenotype`]. It mirrors
/// [`detect_ibd`], but it groups directly from the compact wire type.
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
    // The called sites that both samples share (dosage 0..=2 in both).
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

/// Coarse analysis state for each subject, for the Subjects list. `Pending` means at least one of
/// the subject's alignments has no full `coverage` artifact, as with a new file that no analysis
/// has touched. `Complete` means the app analyzed every alignment. A subject with no alignments at
/// all is absent from the status map, and the list shows it with no status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubjectAnalysisStatus {
    Pending,
    Complete,
}

/// The count for each field from [`App::backfill_read_profiles`], for the CLI `backfill-profiles`
/// summary.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct ReadProfileBackfill {
    /// Runs looked at.
    pub runs_examined: usize,
    /// Runs that gained `total_bases` from a cached read-metrics artifact.
    pub total_bases_filled: usize,
    /// Runs that gained `read_type` with no file I/O, from a low-cost platform or test-type
    /// inference, or from the cached mean-read-length evidence.
    pub read_type_filled: usize,
    /// Runs that gained `read_type` from a bounded read-name rescan (`--rescan`).
    pub read_type_rescanned: usize,
    /// Runs still missing `read_type` (generic-WGS PacBio without `--rescan`, or no readable file).
    pub read_type_unresolved: usize,
}

/// Outcome of [`App::backfill_catalog_ids`], for the CLI `backfill-catalog-ids` summary.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct CatalogBackfill {
    /// Whether the code wrote the ids (`false` = dry run).
    pub applied: bool,
    /// Subjects scanned.
    pub subjects_examined: usize,
    /// Subjects with at least one derivable public-catalog id.
    pub subjects_matched: usize,
    /// Derivable ids not already present (the dry-run "would add" count).
    pub ids_to_add: usize,
    /// Ids added (0 on a dry run).
    pub ids_added: usize,
    /// Ids whose `(namespace, value)` already belongs to a different subject, from a duplicate
    /// import. The code skips them.
    pub conflicts: usize,
}

/// Outcome of [`App::backfill_accessions`], for the CLI `backfill-accessions` summary.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct AccessionBackfill {
    /// Whether the code wrote the accessions (`false` = dry run).
    pub applied: bool,
    /// Subjects queried against the samples API.
    pub examined: usize,
    /// Subjects the API returned a record for.
    pub resolved: usize,
    /// Subjects the API had no record for (404, because the alias is not yet in the catalog).
    pub not_found: usize,
    /// Network/parse failures (counted, not fatal).
    pub errors: usize,
    /// External ids that the code can attach, and that are not already present. These are the
    /// catalog **name** id (IGSR/HGDP) and the authoritative INSDC **accession**. This is the
    /// dry-run "would add" count.
    pub ids_to_add: usize,
    /// External ids attached (0 on a dry run).
    pub ids_added: usize,
    /// Local `sample_accession` placeholders corrected to the authoritative value.
    pub accession_updated: usize,
    /// Ids whose `(namespace, value)` already belongs to another subject. The code skips them.
    pub conflicts: usize,
    /// A few `alias → accession` examples for the report.
    pub examples: Vec<String>,
}

/// Outcome of [`App::prune_orphan_alignments`], for the CLI `prune-orphans` summary.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct PruneReport {
    /// Whether the code applied the deletes (`false` = dry run).
    pub applied: bool,
    /// Alignment records seen in the PDS repo.
    pub examined: usize,
    /// Orphan records deleted (0 on a dry run).
    pub deleted: usize,
    /// The orphan rkeys. The code deletes them when `applied` is true. If not, this is the dry-run
    /// list of what it would remove.
    pub orphans: Vec<String>,
}

/// One row of a project's report for each sample: the coverage roll-up and the haplogroup
/// consensus. A coverage field is `None` until the app computes the coverage. A haplogroup field is
/// `None` until the store holds the calls (this slice defers that).
#[derive(Debug, Clone)]
pub struct ProjectSampleReport {
    pub biosample: Biosample,
    /// An alignment to drive "recompute coverage" from: the one that carries the coverage if there
    /// is one, else the first. `None` if the sample has no alignments.
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
    /// The coverage on the screen is a `partial` result, from a lite sidecar, and a deep walk can
    /// improve it. `false` when the coverage is full, or when there is none yet.
    pub coverage_partial: bool,
    /// The last analysis try on the primary alignment failed, as with a corrupt CRAM that the
    /// code can not decode. This field holds the failure message. `Some` separates a sample that
    /// failed
    /// from one that no analysis has touched. The report then shows "Failed", and not an empty
    /// cell.
    pub decode_error: Option<String>,
}

/// A stored marker that a Navigator walk failed for an alignment, as with a corrupt CRAM that the
/// code can not decode. It goes in the `error`/`"1"` artifact (see
/// [`App::record_analysis_error`]). So the report shows a failure. It does not look the same as a
/// sample that no analysis has touched.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AnalysisError {
    /// Which step failed (`metrics`, `Y`, …).
    pub step: String,
    /// The (truncated) error message.
    pub message: String,
}

/// Artifact key for [`AnalysisError`] markers.
pub(crate) const ERROR_KIND: &str = "error";
pub(crate) const ERROR_VERSION: &str = "1";

/// Artifact key for the pipeline sidecar paths that an alignment came from. This record lets
/// [`App::replace_against_current_tree`] replay the fast path. Without it the code must scan again
/// for files that it found only at import time.
pub(crate) const SIDECARS_KIND: &str = "sidecars";
pub(crate) const SIDECARS_VERSION: &str = "1";

/// One member row of a project's FTDNA-style Y-DNA STR overview. It holds the identity columns, the
/// subject's consensus STR marker values (normalized marker name → value), and the terminal Y
/// haplogroup. The query omits a member with no STR profile.
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
    /// True when SNP evidence supports the haplogroup. The app places from SNP evidence alone, so
    /// it confirms every haplogroup that it places. `false` leaves room for a future STR-predicted label.
    pub y_confirmed: bool,
    /// The highest STR panel or tier reached, such as "Y-111" or "Alpha". This is the "Test" column.
    pub test: Option<String>,
    /// Consensus STR values keyed by uppercase marker name (DYS393 → "13", DYS385 → "11-15").
    pub markers: std::collections::HashMap<String, String>,
}

/// One cell of the project Y-STR chart, ready to draw. It holds the marker value text and its
/// deviation from the subgroup's modal value, which gives the colour. The code computes it in
/// advance, so the UI does no work in a frame.
#[derive(Debug, Clone)]
pub struct StrChartCell {
    pub text: String,
    pub dev: navigator_domain::strchart::Deviation,
}

/// What a [`StrChartRow`] represents. This decides how the UI styles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrRowKind {
    /// A subgroup banner (haplogroup heading); `cells` are empty.
    Group,
    /// The subgroup's MIN row, with one value for each marker.
    Min,
    /// The subgroup's MAX row, with one value for each marker.
    Max,
    /// The subgroup's MODE row, with one value for each marker.
    Mode,
    /// One member's marker values.
    Member,
}

/// One fully-prepared row of the project Y-STR overview, in display order. The code computes the
/// whole table once, off the UI thread, and the renderer only iterates over these.
#[derive(Debug, Clone)]
pub struct StrChartRow {
    pub kind: StrRowKind,
    /// Depth for the indentation. A subgroup sits under its tree ancestors.
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
    /// One cell for each marker, aligned to [`ProjectStrChart::markers`]. It is empty for a
    /// non-member row that does not fill every column.
    pub cells: Vec<StrChartCell>,
}

/// The FTDNA-style project Y-DNA STR overview, computed in advance. It holds the marker column
/// order, and a flat, ordered list of rows that are ready to draw: subgroup banners, MIN/MAX/MODE,
/// and members. The rows group by assigned Y haplogroup, and follow the tree topology, from basal
/// to derived, with children under their parents.
#[derive(Debug, Clone, Default)]
pub struct ProjectStrChart {
    pub markers: Vec<String>,
    pub rows: Vec<StrChartRow>,
    pub member_count: usize,
    pub group_count: usize,
}

/// A reference build that an import needs but does not have in the cache. The app shows it, so
/// the UI can prompt for it and download it before a second try.
#[derive(Debug, Clone)]
pub struct BuildNeed {
    pub build: String,
    pub url: String,
    pub est_bytes: u64,
}

/// Outcome of a project-wide analyze pass: coverage and Y haplogroup for each sample.
#[derive(Debug, Clone)]
pub struct AnalyzeSummary {
    pub project_id: i64,
    pub samples: usize,
    pub coverage_done: usize,
    pub y_done: usize,
    pub sex_done: usize,
    pub metrics_done: usize,
    /// Failures for one sample. The pass is best-effort: one error does not stop the rest.
    pub errors: Vec<String>,
}

/// What the deep analyze pass filled for one biosample, or skipped because it was already there.
/// `had_alignment` is false when the sample has no alignment with a BAM to walk, and the caller
/// then skips it and does not count it. Each `*_done` is true when that artifact is now present,
/// whether the pass computed it or the cache held it. A failure goes in `errors`.
#[derive(Debug, Clone, Default)]
pub struct SampleAnalyzeOutcome {
    pub had_alignment: bool,
    pub coverage_done: bool,
    pub y_done: bool,
    pub sex_done: bool,
    pub metrics_done: bool,
    pub errors: Vec<String>,
}

/// Outcome of a BISDNA chromo2 Y-SNP import: the variant set it created, plus a tally for each
/// category. The UI and the CLI show the coverage and any name the dictionary could not place.
#[derive(Debug, Clone)]
pub struct BisdnaImportSummary {
    pub variant_set: VariantSet,
    /// The reference build that the caller used for the calls, `"hs1"` for example.
    pub build: String,
    /// Total marker rows parsed from the file.
    pub total_markers: usize,
    /// Positive (derived) calls resolved to a locus and emitted as variant calls.
    pub derived_calls: usize,
    /// Negative (ancestral) markers. They are not variants, so the import counts them only.
    pub ancestral: usize,
    /// `no_call` markers (genotype `00`).
    pub no_call: usize,
    /// Back-mutated markers. The import flags them and keeps them out of the placement.
    pub back_mutated: usize,
    /// Markers whose name the dictionary does not have on this build. The code can not place them.
    pub unresolved: usize,
    /// A sample of unresolved names for diagnostics (capped).
    pub unresolved_names: Vec<String>,
    /// Positive calls whose genotype disagreed with the dictionary alleles on either strand. This
    /// is a QC signal. The import still emits the call, and trusts the file.
    pub strand_mismatches: usize,
}

/// Outcome of a batch project-directory import. It is idempotent, and counts only what is new.
#[derive(Debug, Clone)]
pub struct ProjectImportSummary {
    pub project: Project,
    pub samples_total: usize,
    pub samples_created: usize,
    pub alignments_created: usize,
    pub alignments_skipped: usize,
    /// Sample ids whose alignment had no index (.crai/.bai) beside it. Coverage needs one.
    pub missing_index: Vec<String>,
    /// Failures for one sample that the code skipped, so the rest of the batch could import
    /// (`"<sample>: <detail>"`). A non-empty list means the import completed *partially*.
    pub sample_errors: Vec<String>,
    /// Human-readable reference-resolution notes, with one line for each distinct build that the
    /// batch found, such as `"chm13v2.0: 27 alignment(s) → …/chm13v2.0.fa (header probe)"`.
    /// The app shows these, so the user can see which reference the code bound each file to.
    pub reference_notes: Vec<String>,
    /// Roll-up of the fast-path sidecar ingest across the imported samples.
    pub fast_path: FastPathSummary,
}

/// What the fast-path sidecar ingest filled across a project import, with one tally for each
/// result kind. The import then returns immediately, with the report already full.
#[derive(Debug, Clone, Default)]
pub struct FastPathSummary {
    /// Samples that had pipeline sidecars to ingest.
    pub samples_with_sidecars: usize,
    pub y_placed: usize,
    pub mt_placed: usize,
    pub sex_filled: usize,
    pub metrics_filled: usize,
    pub coverage_filled: usize,
    /// Ingest errors for one sample (`"<sample>: <detail>"`), non-fatal.
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
        // Reload the account from the last launch. A keychain error means "nobody".
        let active = tokens.active().ok().flatten();
        Auth {
            tokens,
            config: decodingus_oauth_config(),
            http: dev_http_client(),
            active: Arc::new(Mutex::new(active)),
            online: Arc::new(AtomicBool::new(true)),
        }
    }
}

/// The application. It is low-cost to clone, because the store wraps a connection pool.
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
        // The GUI seeds bundled assets on a background thread, so the window paints first. Wait
        // for it here, on the worker thread, so no analysis can read a half-seeded cache. This
        // does nothing for a headless run, which seeds synchronously first, or for a test, which
        // never spawns it.
        await_bundled_assets();
        Ok(App::new(Store::open(path).await?))
    }
}

/// Render a project report as CSV, with one header row and one row for each sample. A coverage or
/// haplogroup that the app has not computed gets an empty cell. The code formats the text itself,
/// to keep out a CSV dependency, and it quotes a value that holds a comma or a quote.
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

        // A file directly in the picked folder does not count as a subdir.
        let flat = [PathBuf::from("/data/FTDNA/a.csv"), PathBuf::from("/data/FTDNA/b.bam")];
        assert!(contributing_subdirs(root, &flat).is_empty());
    }

    // A six-node spine for the back-off tests: root(1) → A(2,@146) → B(3,@263) → C(4,@750)
    // → D(5,@1000) → F(6,@1100). One SNP defines each node.
    const SPINE6: &str = r#"{ "allNodes": {
      "1": {"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
      "2": {"haplogroupId":2,"name":"A","isRoot":false,"variants":[{"variant":"a","position":146,"ancestral":"A","derived":"G"}],"children":[3]},
      "3": {"haplogroupId":3,"name":"B","isRoot":false,"variants":[{"variant":"b","position":263,"ancestral":"A","derived":"G"}],"children":[4]},
      "4": {"haplogroupId":4,"name":"C","isRoot":false,"variants":[{"variant":"c","position":750,"ancestral":"C","derived":"T"}],"children":[5]},
      "5": {"haplogroupId":5,"name":"D","isRoot":false,"variants":[{"variant":"d","position":1000,"ancestral":"G","derived":"A"}],"children":[6]},
      "6": {"haplogroupId":6,"name":"F","isRoot":false,"variants":[{"variant":"f","position":1100,"ancestral":"C","derived":"T"}],"children":[]}
    }}"#;

    /// The consensus-profile observations come from the placed **lineage**, which holds the
    /// derived mutations from the root to the terminal. They do not come from the child branches
    /// that the descent did not take. This is the regression test for the all-no-call profile. A
    /// clean single-source placement at C must give Derived obs for a/b/c. It must NOT carry the
    /// SNP `d` of child branch D, which is a no-call because the descent stopped.
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

    /// The parsimony back-off trims a net-contradicted deep tail, which a sparse panel or aDNA
    /// makes too deep. It trims to the node where the cumulative (derived − ancestral) support
    /// peaks. But it keeps a clean deep path, and it accepts a lone contradiction that deeper
    /// support outweighs.
    #[test]
    fn support_backoff_trims_net_negative_tail_but_keeps_supported_depth() {
        let tree = parse_ftdna_json(SPINE6).unwrap();
        // Derived A+B, with the peak at B. Below B: ancestral@750, contradiction@1000 (G≠der A),
        // and a lone derived@1100, so the tail is net −1. It must back off F(6) → B(3).
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

    /// Genome-level pool. Take a sparse source that stops shallow on its own, and a dense source
    /// that confirms the deep branches. Together they place the *pooled* call set, by position on
    /// one tree, at the deep terminal. A sparse call from one run no longer pulls the consensus
    /// shallow.
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

    /// A higher-weight source wins the vote at a position: WGS (0.85) derived outvotes a Chip
    /// ancestral call at the same SNP.
    #[test]
    fn pool_vote_prefers_the_higher_weight_source() {
        let wgs: HashMap<i64, char> = [(750, 'T')].into_iter().collect(); // derived
        let chip: HashMap<i64, char> = [(750, 'C')].into_iter().collect(); // ancestral
        let pooled = pool_votes(&[(SourceType::WgsShortRead, wgs), (SourceType::Chip, chip)]);
        assert_eq!(pooled.get(&750), Some(&'T'), "WGS derived outweighs chip ancestral");
    }

    /// The code flips a chip allele on the tree's opposite strand to the ancestral or derived
    /// allele that matches. It does not touch an in-tree match or an out-of-tree position. A
    /// flipped derived call then places as deep as the forward one would.
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

    // root → A(146) → B(263) → C(750) → D(1000). One SNP defines each node.
    const TREE: &str = r#"{ "allNodes": {
      "1": {"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
      "2": {"haplogroupId":2,"name":"A","isRoot":false,"variants":[{"variant":"a","position":146,"ancestral":"A","derived":"G"}],"children":[3]},
      "3": {"haplogroupId":3,"name":"B","isRoot":false,"variants":[{"variant":"b","position":263,"ancestral":"A","derived":"G"}],"children":[4]},
      "4": {"haplogroupId":4,"name":"C","isRoot":false,"variants":[{"variant":"c","position":750,"ancestral":"C","derived":"T"}],"children":[5]},
      "5": {"haplogroupId":5,"name":"D","isRoot":false,"variants":[{"variant":"d","position":1000,"ancestral":"G","derived":"A"}],"children":[]}
    }}"#;

    /// A deep lineage with a single stray ancestral call on a backbone node (C). This is the
    /// sparse-chip failure mode. Strict selection vetoes the whole lineage and stops shallow, at
    /// B. Robust selection trusts the proportional top and reaches the deep terminal, at D.
    #[test]
    fn robust_selection_survives_a_backbone_contradiction() {
        let tree = parse_ftdna_json(TREE).unwrap();
        // Derived at 146, 263, and 1000, but ANCESTRAL (C) at 750: a lone contradiction on C.
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

    /// The GVCF fast path rebuilds exactly the `calls` that a pileup would give. A fully derived
    /// path, where every SNP of a node is a variant, places to the deep terminal D.
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

    /// A hom-ref tree SNP, which is callable with no variant, rebuilds as the **reference base**.
    /// On a real reference that base can be the *derived* allele (CHM13 Y = J1). Here the
    /// reference base at position 750 is the derived T. So the evidence supports node C, and the
    /// placement reaches D. The old "assume ancestral" logic got this case wrong, and stopped at B.
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

    /// Lifted assembly maps GVCF observations back to tree positions, and reverse-complements a
    /// minus-strand lift. A hom-ref lifted site takes the reference base at the lifted position.
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
mod external_precedence_tests {
    use super::*;
    use navigator_store::{haplogroup_call, Store};

    fn ycall(label: &str, hg: &str, score: f64) -> RunHaplogroupCall {
        RunHaplogroupCall {
            source_label: label.into(),
            haplogroup: hg.into(),
            lineage: vec!["R".into(), hg.into()],
            score,
            matched: 0,
            expected: 0,
        }
    }

    /// PRJEB37976 idempotence gate. An external (GATK4 GVCF) Y call and Navigator's internal walk
    /// on the *same* alignment stay as two distinct rows, and neither overwrites the other. With
    /// the default "prefer external caller" policy, the external terminal wins the consensus, even
    /// though the damaged ancient-DNA walk scored higher. Before the fix, the walk shared the
    /// `aln:1` key and overwrote the external call.
    #[tokio::test]
    async fn external_call_is_not_clobbered_and_wins_consensus() {
        // Deterministic policy regardless of the developer's settings.json (env wins in the resolver).
        std::env::set_var("NAVIGATOR_PREFER_EXTERNAL_CALLS", "true");
        let app = App::new(Store::open_in_memory().await.unwrap());
        let bio = app.add_biosample(None, "ANCIENT-1", None, None).await.unwrap();

        // External placement: distinct `:ext` key + External provenance, as the sidecar fast path writes it.
        haplogroup_call::upsert(
            app.store.pool(),
            bio.guid,
            DnaType::Y,
            &external_y_source_key(1),
            &ycall("gvcf", "R-EXT", 0.60),
            CallProvenance::External,
            Some("gv:abc"),
        )
        .await
        .unwrap();

        // Internal walk on the SAME alignment: a higher score, a DIFFERENT (wrong, aDNA-damaged) terminal.
        app.record_haplogroup_call(bio.guid, DnaType::Y, "aln:1", &ycall("walk", "R-WRONG", 0.95))
            .await
            .unwrap();

        // No clobber: both rows survive under their distinct keys.
        assert!(
            haplogroup_call::get_one(app.store.pool(), bio.guid, DnaType::Y, &external_y_source_key(1))
                .await
                .unwrap()
                .is_some()
        );
        assert!(
            haplogroup_call::get_one(app.store.pool(), bio.guid, DnaType::Y, "aln:1")
                .await
                .unwrap()
                .is_some()
        );

        // Default policy prefers external → external terminal wins despite the walk's higher score.
        let c = app.haplogroup_consensus(bio.guid, DnaType::Y).await.unwrap().unwrap();
        assert_eq!(c.haplogroup, "R-EXT");

        std::env::remove_var("NAVIGATOR_PREFER_EXTERNAL_CALLS");
    }
}

#[cfg(test)]
mod publish_tests {
    use super::*;
    use navigator_domain::workspace::NewSequenceRun;
    use navigator_store::Store;

    /// The published sequence-run record carries the inferred `instrumentId` (camelCase) so the
    /// AppView can crowd-source the instrument→lab map. This is a regression guard: the field
    /// held a hardcoded `None` during the work that restored the lab inference.
    #[tokio::test]
    async fn sequence_run_record_publishes_instrument_id() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                instrument_model: Some("NovaSeq".into()),
                ..NewSequenceRun::new(b.guid, "ILLUMINA", "WGS")
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
            Some("SHORT"),
        )
        .await
        .unwrap();
        // Exact yield → the standardized label's Gbases figure.
        sequence_run::set_read_stats(
            app.store.pool(),
            run.id,
            Some(300_000_000),
            Some(150.0),
            None,
            None,
            Some(45_000_000_000),
        )
        .await
        .unwrap();
        sequence_run::set_facility(app.store.pool(), run.id, "Dante Labs")
            .await
            .unwrap();
        let reloaded = sequence_run::get(app.store.pool(), run.id).await.unwrap().unwrap();

        let value = app.sequence_run_record("did:plc:test", &reloaded).await.unwrap();
        assert_eq!(value.get("instrumentId").and_then(|v| v.as_str()), Some("A00182"));
        // The record publishes the known sequencing lab, so the AppView can show it. Its
        // instrument→lab map does not cover every serial, PacBio for example.
        assert_eq!(
            value.get("sequencingFacility").and_then(|v| v.as_str()),
            Some("Dante Labs")
        );
        // The record publishes the read-profile fields behind the standardized label.
        assert_eq!(value.get("totalBases").and_then(|v| v.as_i64()), Some(45_000_000_000));
        assert_eq!(value.get("readType").and_then(|v| v.as_str()), Some("SHORT"));
        assert_eq!(value.get("$type").and_then(|v| v.as_str()), Some(NS_SEQUENCERUN));
        // Links back to its subject's biosample record through the deterministic at:// URI.
        assert_eq!(
            value.get("biosampleRef").and_then(|v| v.as_str()),
            Some(biosample_at_uri("did:plc:test", reloaded.biosample_guid).as_str())
        );
    }

    /// The published biosample record carries the subject's external identifiers as
    /// `externalIds[{namespace,value}]`. That is the deterministic dedup anchor that the AppView
    /// keys on. Both vendor kits and public catalog ids go in plaintext, and the AppView gates
    /// their visibility by namespace.
    #[tokio::test]
    async fn biosample_record_publishes_external_ids() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "huF98AFD", None, None).await.unwrap();
        app.add_external_id(b.guid, "YSEQ", "229").await.unwrap();
        app.add_external_id(b.guid, "PGP", "huF98AFD").await.unwrap();

        let value = app.biosample_record("did:plc:test", b.guid).await.unwrap();
        let ids = value
            .get("externalIds")
            .and_then(|v| v.as_array())
            .expect("externalIds present");
        let mut pairs: Vec<(String, String)> = ids
            .iter()
            .map(|e| {
                (
                    e.get("namespace")
                        .and_then(|v| v.as_str())
                        .unwrap_or_default()
                        .to_string(),
                    e.get("value").and_then(|v| v.as_str()).unwrap_or_default().to_string(),
                )
            })
            .collect();
        pairs.sort();
        // Pure field rename source→namespace, external_id→value; both a vendor kit and a public id.
        assert_eq!(
            pairs,
            vec![("PGP".into(), "huF98AFD".into()), ("YSEQ".into(), "229".into())]
        );
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
    /// subject's match profile from cached data, and rank by Y-STR genetic distance, lowest first.
    /// There is no Y tree and there are no SNP profiles, so this drives the offline path.
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

    /// The code drops a subject with no comparable Y data. An empty workspace gives no matches.
    #[tokio::test]
    async fn y_matches_drops_incomparable_and_self() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let q = app.add_biosample(None, "Query", None, None).await.unwrap();
        seed_str(&app, q.guid, &[("DYS393", "13")]).await;
        // A subject with no STR / Y data at all.
        let _empty = app.add_biosample(None, "Empty", None, None).await.unwrap();
        // A subject whose markers do not overlap the query's.
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

    /// A signed attestation checks out on the same code path that the AppView runs. A change to
    /// the bytes breaks it.
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

    /// Two peers that compute the same summary produce the same agreement hash. Two different
    /// summaries do not.
    #[test]
    fn summary_hash_drives_agreement() {
        use navigator_analysis::ibd_attest::summary_hash;
        assert_eq!(summary_hash(&summary(42.0)), summary_hash(&summary(42.0)));
        assert_ne!(summary_hash(&summary(42.0)), summary_hash(&summary(7.0)));
    }

    /// The store keeps an exchange result and reads it back for each subject (the UI's
    /// saved-results path).
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

        // The ledger adopts a conversation whose open it never saw. So the completed exchange
        // still reads as one entry with its result attached, and not as an orphan row.
        app.mark_matching_exchanged(b.guid, &session, "exchange:r")
            .await
            .unwrap();
        let entries = app.matching_entries().await.unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].status, MatchingStatus::Exchanged);
        assert_eq!(entries[0].partner_did.as_deref(), Some("did:key:zB"));
        assert_eq!(entries[0].biosample_guid, Some(b.guid));
        assert_eq!(entries[0].result.as_ref().map(|r| r.total_shared_cm), Some(75.0));
    }

    /// A failed exchange must read as failed, with its reason, and must not sit at READY for
    /// ever. To forget a conversation drops it locally.
    #[tokio::test]
    async fn matching_failure_and_forget() {
        use navigator_store::Store;
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, None).await.unwrap();
        let session = EstablishedSession {
            session_id: "sess-y".into(),
            partner_did: "did:key:zC".into(),
            key: [0u8; 32],
        };
        app.mark_matching_exchanged(b.guid, &session, "exchange:f")
            .await
            .unwrap();

        app.record_matching_failure("exchange:f", "relay timeout")
            .await
            .unwrap();
        let e = app.matching_entry("exchange:f").await.unwrap();
        assert_eq!(e.status, MatchingStatus::Failed);
        assert_eq!(e.last_error.as_deref(), Some("relay timeout"));
        assert!(!e.status.is_terminal(), "a failure stays retryable");

        // A record against an unknown request does nothing, and it is not an error.
        app.record_matching_failure("exchange:nope", "x").await.unwrap();
        assert!(app.matching_entry("exchange:nope").await.is_err());

        app.forget_matching_request("exchange:f").await.unwrap();
        assert!(app.matching_entries().await.unwrap().is_empty());
    }

    /// A gate protects the attestation. Without the AppView sample handles there is nothing to
    /// file, and the app must not report a disputed summary as a match. Neither case is an error.
    #[tokio::test]
    async fn attest_is_skipped_without_handles_or_agreement() {
        use navigator_store::Store;
        let app = App::new(Store::open_in_memory().await.unwrap());
        let b = app.add_biosample(None, "S1", None, None).await.unwrap();
        let session = EstablishedSession {
            session_id: "sess-z".into(),
            partner_did: "did:key:zD".into(),
            key: [0u8; 32],
        };
        let mut result = IbdExchangeResult {
            summary: summary(75.0),
            segments: vec![],
            overlapping_sites: 100,
            my_attestation: IbdAttestation::unsigned(
                "exchange:a",
                "sess-z",
                "did:key:zA",
                Some(b.guid.to_string()),
                None,
                &summary(75.0),
                "t",
            ),
            partner_attestation: IbdAttestation::unsigned(
                "exchange:a",
                "sess-z",
                "did:key:zD",
                None,
                Some(b.guid.to_string()),
                &summary(75.0),
                "t",
            ),
            agreed: true,
        };
        app.record_ibd_exchange(b.guid, &session, "exchange:a", &result)
            .await
            .unwrap();
        app.mark_matching_exchanged(b.guid, &session, "exchange:a")
            .await
            .unwrap();

        // No sample handles (a direct request never carries them) → nothing to attest, no network.
        assert!(!app.attest_exchange_if_possible("exchange:a").await.unwrap());

        // With handles but a disputed summary, we still file nothing.
        app.set_matching_sample_refs("exchange:a", Some("s-mine"), Some("s-theirs"))
            .await
            .unwrap();
        result.agreed = false;
        app.record_ibd_exchange(b.guid, &session, "exchange:a", &result)
            .await
            .unwrap();
        assert!(!app.attest_exchange_if_possible("exchange:a").await.unwrap());
        assert!(!app.matching_entry("exchange:a").await.unwrap().attested);
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

    /// The tier boundaries decide which of three claims the app makes about a stranger's
    /// relatedness. So this code holds them, and the widget that draws the card does not.
    #[test]
    fn match_strength_tiers_are_conservative() {
        let at = |score: f64| {
            IbdSuggestion {
                target_sample_guid: None,
                suggested_sample_guid: "h".into(),
                suggestion_type: "SHARED_MATCH".into(),
                score,
                signals: Vec::new(),
            }
            .strength()
        };
        // Boundaries are inclusive at the bottom of each tier.
        assert_eq!(at(1.0), MatchStrength::Strong);
        assert_eq!(at(0.8), MatchStrength::Strong);
        assert_eq!(at(0.79), MatchStrength::Likely);
        assert_eq!(at(0.5), MatchStrength::Likely);
        assert_eq!(at(0.49), MatchStrength::Possible);
        assert_eq!(at(0.0), MatchStrength::Possible);
        // A missing score parses as 0.0, which must read as the weakest claim, never the strongest.
        assert_eq!(
            at(f64::NAN),
            MatchStrength::Possible,
            "an unusable score must not overstate"
        );
    }
}

#[cfg(test)]
mod export_tests {
    use super::*;

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
        let app = App::new(Store::open_in_memory().await.unwrap());
        assert!(matches!(app.publish_coverage(1).await, Err(AppError::NotAuthenticated)));
        // No account → nothing enqueued, and the accessors degrade gracefully.
        assert_eq!(app.outbox_pending_count().await.unwrap(), 0);
        assert!(app.outbox_entries().await.unwrap().is_empty());
        assert!(app.sync_history(10).await.unwrap().is_empty());
        // A drain with no account does nothing, and it does no harm.
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
        assert_eq!(resolve_appview_url(None, None), DEFAULT_APPVIEW_URL);
        // the code ignores a blank value, and falls through to the default
        assert_eq!(resolve_appview_url(Some("".into()), None), DEFAULT_APPVIEW_URL);
    }

    #[test]
    fn deterministic_pds_refs_are_stable() {
        let g = SampleGuid(Uuid::nil());
        // Deterministic rkeys: the code knows the biosample's at:// URI before it publishes. So
        // child records can reference it, and a second publish overwrites and makes no duplicate.
        assert_eq!(biosample_rkey(g), "bio-00000000000000000000000000000000");
        assert_eq!(seqrun_rkey(7), "run-7");
        assert_eq!(
            biosample_at_uri("did:plc:abc", g),
            "at://did:plc:abc/com.decodingus.atmosphere.biosample/bio-00000000000000000000000000000000"
        );
        assert_eq!(
            seqrun_at_uri("did:plc:abc", 7),
            "at://did:plc:abc/com.decodingus.atmosphere.sequencerun/run-7"
        );
    }

    #[test]
    fn oauth_config_defaults_hosted_and_env_overrides() {
        let redirect = "http://127.0.0.1:5001/callback";

        // No override → the hosted production native client + write scope.
        let prod = resolve_oauth_config(None);
        assert_eq!(prod.client_id(redirect), DEFAULT_OAUTH_CLIENT_ID);
        assert_eq!(prod.scope, OAUTH_SCOPE);
        assert!(
            prod.scope.contains("transition:generic"),
            "publishing needs write scope"
        );

        // The code ignores a blank value, so this is still the hosted default.
        assert_eq!(
            resolve_oauth_config(Some("  ".into())).client_id(redirect),
            DEFAULT_OAUTH_CLIENT_ID
        );

        // `loopback` selects the dev client (client_id derived from the loopback redirect).
        let dev = resolve_oauth_config(Some("loopback".into()));
        assert!(dev.client_id(redirect).starts_with("http://localhost?redirect_uri="));

        // Any other value is a hosted client-metadata URL (e.g. a local dev document).
        assert_eq!(
            resolve_oauth_config(Some("https://dev.example/cm.json".into())).client_id(redirect),
            "https://dev.example/cm.json"
        );
    }

    #[test]
    fn app_settings_serde_round_trip_and_defaults() {
        let s = AppSettings {
            y_tree_provider: Some("ftdna".into()),
            prefer_external_calls: Some(false),
            appview_url: Some("https://av.example".into()),
            tree_ttl_days: Some(3),
            theme: Some("light".into()),
            prompt_before_download: Some(false),
            ui_scale: Some(1.5),
            ui_mode: Some("simple".into()),
            llm_enabled: Some(true),
            llm_base_url: Some("http://localhost:1234/v1".into()),
            llm_model: Some("llama-3.1-8b-instruct".into()),
            llm_max_tokens: Some(8192),
            check_for_updates: Some(false),
            skip_update_version: Some("0.2.0-alpha".into()),
            window_size: Some([1280.0, 800.0]),
            last_nav: Some("subjects".into()),
            last_subject: Some("11111111-1111-1111-1111-111111111111".into()),
            last_detail_tab: Some("ancestry".into()),
            lai_recomb_per_cm: Some(0.15),
            lai_max_ref_haps: Some(40),
            lai_min_ancestry: Some(0.03),
            lai_switch_per_cm: Some(0.04),
            lai_min_segment_cm: Some(3.0),
            lai_size_normalize: Some(0.4),
            lai_mismatch: Some(0.01),
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
        // the code ignores a token it does not recognize
        assert_eq!(resolve_ui_mode(Some("expert"), Some("simple")), Some(UiMode::Simple));
        assert_eq!(resolve_ui_mode(Some("expert"), None), None);
    }
}

#[cfg(test)]
mod vcf_genotype_tests {
    use super::vcf_genotypes_at;
    use std::collections::HashSet;

    fn genotypes(name: &str, body: &str, targets: &[i64]) -> std::collections::HashMap<i64, char> {
        let dir = std::env::temp_dir().join(format!("nav-vcf-gt-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{name}.vcf"));
        std::fs::write(&path, body).unwrap();
        let t: HashSet<i64> = targets.iter().copied().collect();
        let out = vcf_genotypes_at(&path, "chrY", &t).unwrap();
        let _ = std::fs::remove_file(&path);
        out
    }

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

    #[test]
    fn hom_ref_rows_become_ancestral_evidence() {
        // The reason this exists: a 0/0 row says the donor is *ancestral* here. An import of the
        // non-reference rows alone loses that, and placement can not tell ancestral from
        // uncovered.
        let g = genotypes(
            "homref",
            &format!(
                "{HEADER}chrY\t100\t.\tA\tG\t99\tPASS\t.\tGT\t1/1\n\
                 chrY\t200\t.\tC\tT\t99\tPASS\t.\tGT\t0/0\n"
            ),
            &[100, 200],
        );
        assert_eq!(g.get(&100), Some(&'G'), "derived → the ALT base");
        assert_eq!(g.get(&200), Some(&'C'), "hom-ref → the REF base");
    }

    #[test]
    fn a_no_call_is_not_ancestral() {
        // `./.` also yields no non-zero allele index, but it is not evidence of the reference.
        let g = genotypes(
            "nocall",
            &format!("{HEADER}chrY\t300\t.\tA\tG\t99\tPASS\t.\tGT\t./.\n"),
            &[300],
        );
        assert!(g.is_empty(), "a no-call must stay absent, not read as ancestral");
    }

    #[test]
    fn only_target_positions_and_the_right_contig_are_reported() {
        let g = genotypes(
            "targets",
            &format!(
                "{HEADER}chrY\t100\t.\tA\tG\t99\tPASS\t.\tGT\t1/1\n\
                 chrY\t999\t.\tA\tG\t99\tPASS\t.\tGT\t1/1\n\
                 chr1\t100\t.\tA\tG\t99\tPASS\t.\tGT\t1/1\n"
            ),
            &[100],
        );
        assert_eq!(g.len(), 1, "off-target and off-contig rows are ignored");
        assert_eq!(g.get(&100), Some(&'G'));
    }

    #[test]
    fn the_called_alt_is_used_on_a_multiallelic_row() {
        let g = genotypes(
            "multi",
            &format!("{HEADER}chrY\t400\t.\tA\tG,T\t99\tPASS\t.\tGT\t2/2\n"),
            &[400],
        );
        assert_eq!(g.get(&400), Some(&'T'), "GT 2 selects ALT[1]");
    }

    #[test]
    fn indels_and_sites_only_rows_are_skipped() {
        let g = genotypes(
            "skip",
            &format!(
                "{HEADER}chrY\t500\t.\tA\tAT\t99\tPASS\t.\tGT\t1/1\n\
                 chrY\t600\t.\tAT\tA\t99\tPASS\t.\tGT\t0/0\n"
            ),
            &[500, 600],
        );
        assert!(g.is_empty(), "no single observed base to report for an indel");

        // A sites-only row has no sample column, so it says nothing about *this* donor.
        let sites = genotypes(
            "sitesonly",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrY\t700\t.\tA\tG\t99\t.\t.\n",
            &[700],
        );
        assert!(sites.is_empty());
    }

    #[test]
    fn contig_naming_is_matched_leniently() {
        // GRCh37-style bare `Y` must match a `chrY` query (see contig::bare_upper).
        let g = genotypes(
            "bareY",
            &format!("{HEADER}Y\t800\t.\tA\tG\t99\tPASS\t.\tGT\t1/1\n"),
            &[800],
        );
        assert_eq!(g.get(&800), Some(&'G'));
    }
}

#[cfg(test)]
mod vcf_evidence_tests {
    use super::parse_vcf_subject_snps;

    /// `name` keeps each test on its own file. The tests run in parallel, and a shared path lets
    /// one test delete the file that another still reads.
    fn parse(name: &str, body: &str) -> Vec<navigator_domain::variants::VariantCall> {
        let dir = std::env::temp_dir().join(format!("nav-vcf-ev-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("{name}.vcf"));
        std::fs::write(&path, body).unwrap();
        let out = parse_vcf_subject_snps(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        out
    }

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

    #[test]
    fn captures_qual_depth_and_allele_depths() {
        // The sample column is ONE colon-separated field, keyed by FORMAT.
        let v = parse(
            "qual",
            &format!("{HEADER}chrY\t100\trs1\tA\tG\t512.7\tPASS\t.\tGT:AD:DP:GQ\t1/1:2,38:40:99\n"),
        );
        assert_eq!(v.len(), 1);
        let e = &v[0].evidence;
        assert_eq!(e.qual, Some(512.7));
        assert_eq!(e.dp, Some(40));
        assert_eq!(e.gq, Some(99));
        assert_eq!(e.ad_ref, Some(2));
        assert_eq!(e.ad_alt, Some(38));
        assert_eq!(e.allele_fraction(), Some(0.95));
        assert!(!e.is_filtered(), "PASS is not a failure");
    }

    #[test]
    fn ad_alt_follows_the_called_allele_on_a_multiallelic_row() {
        // GT 2 selects ALT[1] = T, so the code must read AD at index 2, and not at index 1.
        let v = parse(
            "multi",
            &format!("{HEADER}chrY\t200\t.\tA\tG,T\t99\t.\t.\tGT:AD\t2/2:1,3,30\n"),
        );
        assert_eq!(v[0].alternate, "T", "the genotype-selected ALT is kept");
        assert_eq!(v[0].evidence.ad_ref, Some(1));
        assert_eq!(v[0].evidence.ad_alt, Some(30), "AD index follows the ALT index");
    }

    #[test]
    fn a_failing_filter_is_recorded_but_pass_and_dot_are_not() {
        let v = parse(
            "filter",
            &format!(
                "{HEADER}chrY\t300\t.\tA\tG\t10\tLowQual\t.\tGT\t1/1\n\
                 chrY\t301\t.\tA\tG\t10\t.\t.\tGT\t1/1\n"
            ),
        );
        assert_eq!(v[0].evidence.filter.as_deref(), Some("LowQual"));
        assert!(v[0].evidence.is_filtered());
        assert_eq!(v[1].evidence.filter, None, "'.' is not a failure");
    }

    #[test]
    fn absent_evidence_stays_absent_rather_than_becoming_zero() {
        // A sites-only VCF: no FORMAT/sample columns, QUAL '.'.
        let v = parse(
            "sitesonly",
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrY\t400\t.\tA\tG\t.\t.\t.\n",
        );
        assert_eq!(v.len(), 1);
        let e = &v[0].evidence;
        assert!(e.is_empty(), "nothing captured → empty, so the set tags as BASIC");
        assert_eq!(e.dp, None, "missing DP must not read as 0 supporting reads");
        assert_eq!(e.allele_fraction(), None, "no AD → no fraction, not 0.0");
    }

    #[test]
    fn a_reference_or_nocall_genotype_is_still_skipped() {
        let v = parse(
            "refcall",
            &format!(
                "{HEADER}chrY\t500\t.\tA\tG\t99\tPASS\t.\tGT:DP\t0/0:40\n\
                 chrY\t501\t.\tA\tG\t99\tPASS\t.\tGT:DP\t./.:40\n"
            ),
        );
        assert!(v.is_empty(), "the subject carries no ALT at either site");
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
    use super::{seed_assets_from, tree_etag_path, SeedSummary};
    use std::path::Path;

    #[test]
    fn etag_sidecar_appends_not_replaces_extension() {
        // `<cache>.etag`: append, so the tree keeps its own `.json` extension. A `with_extension`
        // would drop it, and two trees with the same stem could then collide.
        assert_eq!(
            tree_etag_path(Path::new("/t/decodingus-ytree.json")),
            Path::new("/t/decodingus-ytree.json.etag")
        );
    }

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
        // `b.bin` already exists in the cache, a copy refreshed from the CDN for example. Keep it.
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

/// Genotype a VCF **at a fixed set of target positions**. This is the file-based counterpart of a
/// walk over a BAM/CRAM at the tree's sites ([`App::base_calls`]).
///
/// Returns `position → observed base` for every target that the VCF has something to say about:
///
/// - a non-reference genotype gives the **called ALT** base, and the donor carries the derived
///   allele here;
/// - an explicit hom-ref (`0/0` / `0|0`) gives the **REF** base, and the donor is *ancestral*;
/// - no record, or a no-call (`./.`), gives nothing, which is a genuine no-call.
///
/// That middle case is the whole point. An import of the non-reference rows alone leaves the
/// workspace unable to separate "ancestral" from "not covered". Every backbone node then scores as
/// a no-call, and the placement runs on a few dozen sites. A vendor Y export already carries the
/// hom-ref rows: an aengine Big Y is about 218k PASS records, most of them `0/0`. This reads them.
///
/// The code uses only single-base REF/ALT rows. An indel has no single observed base to report,
/// and the targets describe the tree's SNP loci. The `contig` match is lenient (`chrY` == `Y`).
fn vcf_genotypes_at(
    path: &Path,
    contig: &str,
    targets: &std::collections::HashSet<i64>,
) -> Result<HashMap<i64, char>, AppError> {
    use std::io::BufRead;
    let want = navigator_domain::contig::bare_upper(contig);
    let reader = navigator_analysis::gzio::open_maybe_gz(path)?;
    let mut out = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let Ok(pos) = f[1].parse::<i64>() else { continue };
        if !targets.contains(&pos) || navigator_domain::contig::bare_upper(f[0]) != want {
            continue;
        }
        let (reference, alt_field) = (f[3], f[4]);
        let ref_base = match single_base(reference) {
            Some(b) => b,
            None => continue,
        };
        // Genotype selects which allele the donor carries; without a sample column the row is a
        // sites-only listing and says nothing about *this* donor's state.
        let gt = (f.len() >= 10)
            .then(|| {
                f[8].split(':')
                    .position(|k| k == "GT")
                    .and_then(|i| f[9].split(':').nth(i))
            })
            .flatten();
        let Some(gt) = gt else { continue };
        let idx = gt
            .split(['/', '|'])
            .filter_map(|a| a.parse::<usize>().ok())
            .find(|&a| a > 0);
        match idx {
            // Derived: the ALT that the genotype selected, and not ALT[0].
            Some(i) => {
                if let Some(b) = alt_field.split(',').nth(i - 1).and_then(single_base) {
                    out.insert(pos, b);
                }
            }
            // Ancestral, but only for an explicit hom-ref. `./.` also parses to no indices, and a
            // no-call is not evidence of the reference allele.
            None if gt.split(['/', '|']).any(|a| a == "0") => {
                out.insert(pos, ref_base);
            }
            None => {}
        }
    }
    Ok(out)
}

/// The single upper-cased base of a one-character A/C/G/T allele, else `None` (indel/symbolic).
fn single_base(allele: &str) -> Option<char> {
    let mut cs = allele.chars();
    let c = cs.next()?.to_ascii_uppercase();
    (cs.next().is_none() && matches!(c, 'A' | 'C' | 'G' | 'T')).then_some(c)
}

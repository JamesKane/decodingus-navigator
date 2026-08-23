//! Multi-source variant **consensus engine**, with no dependence on the DNA type.
//!
//! The input is a set of sources: the placement of a WGS alignment, a chip or BISDNA panel, a
//! private bucket, and others. Each source gives calls at each variant, with the **name** as the
//! key. The name is independent of the build, because M269 is M269 whether the source aligned to
//! GRCh37 or to GRCh38.
//!
//! [`reconcile`] groups those calls and weight-votes the consensus state. It classifies each
//! variant as confirmed, novel, conflict, or single-source, and it calculates a quality-weighted
//! confidence. This mirrors the Scala `YVariantConcordance`.
//!
//! This engine is the shared foundation for the Y-DNA profile, through the [`crate::yprofile`]
//! adapter today. By design it is also the foundation for the future mtDNA consumer (variants vs
//! rCRS) and the future autosomal consumer. It holds nothing specific to a DNA type. A caller
//! collects the observations and gives the variant identity. The DNA type, and the consensus label
//! (a haplogroup, where that applies), live at the persistence layer and the app layer.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::variants::SourceType;

/// One source's call state at a variant position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusState {
    /// Carries the derived (mutant) allele, so it is positive for the branch of the variant. For
    /// mtDNA this is "differs from rCRS". For autosomes, a future adapter maps a diploid genotype
    /// onto this axis.
    Derived,
    /// Carries the ancestral (reference) allele.
    Ancestral,
    /// No confident call.
    NoCall,
}

/// Cross-source status of a variant after reconciliation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusStatus {
    /// ≥2 sources agree on the consensus state and the variant is a known tree/reference variant.
    Confirmed,
    /// Derived but not a known tree variant (private / off-path).
    Novel,
    /// Sources disagree (weighted minority > 30%).
    Conflict,
    /// Only one source reports the variant.
    SingleSource,
    /// It has data, but the weighted confidence is below the confirmation threshold, and it does
    /// not reach the conflict line. This is rare, and it stays for parity with the Scala
    /// `YVariantConcordance`.
    Pending,
    /// No source made a confident call (every observation was NoCall).
    NoCoverage,
}

/// How callable each position of an observation is. It scales the concordance weight of that
/// observation, because a base in a region with no coverage, or with poor mapping, carries little
/// confidence. Mirrors the Scala `YCallableState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CallableState {
    Callable,
    LowCoverage,
    ExcessiveCoverage,
    PoorMappingQuality,
    NoCoverage,
    RefN,
}

impl CallableState {
    /// Confidence multiplier (Scala weights): full for CALLABLE, none for NO_COVERAGE / REF_N.
    pub fn weight(self) -> f64 {
        match self {
            CallableState::Callable => 1.0,
            CallableState::LowCoverage => 0.5,
            CallableState::ExcessiveCoverage | CallableState::PoorMappingQuality => 0.3,
            CallableState::NoCoverage | CallableState::RefN => 0.0,
        }
    }
}

/// One source's observation of a variant (for provenance display).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceObs {
    pub label: String,
    pub source_type: SourceType,
    pub state: ConsensusState,
    /// The **observed base** (allele) this source called at the variant. `None` for a no-call, and
    /// for a source or legacy profile that carries only a state. The store keeps the base, and not
    /// only the derived or ancestral reading of it. [`reproject`] can then run [`impute_state`]
    /// again against a corrected tree polarity, or a different one, and it never reads the BAM or
    /// CRAM again.
    #[serde(default)]
    pub base: Option<String>,
}

/// A reconciled variant across the subject's sources.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConsensusVariant {
    /// Variant name (e.g. "M269"); for unnamed/novel calls this is a `@<position>` placeholder.
    pub name: String,
    /// A representative position (from the consensus-side sources; builds may differ).
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    /// The **consensus observed base**: the weighted-majority nucleotide over the sources, with the
    /// strand normalized to the alleles of this SNP. This is the primary observation.
    /// [`consensus`](Self::consensus) is its derived or ancestral reading against the tree. `None`
    /// means no source made a call. A base that matches neither allele is a genuine third allele,
    /// and it survives here as itself.
    #[serde(default)]
    pub consensus_base: Option<String>,
    pub consensus: ConsensusState,
    pub status: ConsensusStatus,
    /// Sources matching the consensus state.
    pub support: usize,
    /// Sources with any call (excludes NoCall).
    pub total: usize,
    /// Whether the variant is a known reference/haplotree variant (vs a private/novel call).
    pub in_tree: bool,
    /// Weighted confidence in the consensus = consensusWeight / totalWeight (0 when no call).
    #[serde(default)]
    pub confidence_score: f64,
    pub sources: Vec<SourceObs>,
}

/// Counts for each status, for the profile header.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct ConsensusSummary {
    pub total: usize,
    pub confirmed: usize,
    pub novel: usize,
    pub conflict: usize,
    pub single_source: usize,
    /// Overall profile confidence: `(confirmed + 0.7·novel − 0.5·conflict) / total`, clamped [0,1].
    #[serde(default)]
    pub overall_confidence: f64,
}

// ---------------------------------------------------------------------------------------------
// Observation-first storage. A persisted profile holds only OBSERVATIONS: for each SNP and each
// source, the observed base, the quality, and the identity. It never holds a fixed derived or
// ancestral reading. [`interpret`] calculates the state, the vote, the status, the support and the
// summary on demand, against the polarity of the CURRENT tree. So a fix to the tree polarity, or a
// switch of provider, corrects every view, and nothing genotypes again. This is the type the code
// writes to the `consensus_profile` payload.
// ---------------------------------------------------------------------------------------------

fn one() -> f64 {
    1.0
}
fn schema_v1() -> u8 {
    1
}

/// The raw observation of a variant by one source: the **observed base**, plus the quality inputs
/// to the concordance weight. It carries no derived or ancestral state, because [`impute_state`]
/// makes that at read time. The store keeps depth, MQ, callable and region. The in-memory tally of
/// `reconcile` used those and never stored them. Now [`interpret`] can weight exactly, which fixes
/// the loss of quality that the old `reproject` warned about.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ObservedSource {
    pub label: String,
    pub source_type: SourceType,
    /// Observed allele; `None` = no confident call at this position for this source.
    #[serde(default)]
    pub base: Option<String>,
    #[serde(default)]
    pub depth: Option<u32>,
    #[serde(default)]
    pub mapq: Option<f64>,
    #[serde(default)]
    pub callable: Option<CallableState>,
    #[serde(default = "one")]
    pub region_modifier: f64,
}

/// A variant that the sources of the subject observed: the identity, plus the observed base from
/// each source. The derived and ancestral polarity comes from the current tree at [`interpret`]
/// time, by name. The stored `ref_allele` and `alt_allele` are the polarity fallback. That fallback
/// serves a novel or private call off the tree, and an mtDNA mutation the tree map does not hold.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ObservedVariant {
    /// Variant name (e.g. "M269"); empty for a novel/unnamed call (then keyed by position).
    pub name: String,
    pub position: i64,
    pub in_tree: bool,
    #[serde(default)]
    pub ref_allele: Option<String>,
    #[serde(default)]
    pub alt_allele: Option<String>,
    pub sources: Vec<ObservedSource>,
}

/// The provenance of one source that contributes (label, type, count). It reads nothing into the
/// data, and it is here for display.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceSummary {
    pub label: String,
    pub source_type: SourceType,
    pub variant_count: usize,
}

/// The persisted, observation-only profile (the `consensus_profile` payload). Interpreted into the
/// display view (`ConsensusVariant` + `ConsensusSummary`) on demand by [`interpret`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ObservedProfile {
    /// Payload schema tag. If it is there, this is not a legacy fixed `ConsensusProfile` JSON.
    #[serde(default = "schema_v1")]
    pub schema_version: u8,
    pub variants: Vec<ObservedVariant>,
    #[serde(default)]
    pub sources: Vec<SourceSummary>,
    /// The terminal haplogroup label of the placement. It is an output of the placement, and not a
    /// reading of each SNP.
    #[serde(default)]
    pub terminal_hint: Option<String>,
}

/// The call of one source at a variant, for [`reconcile`]. The quality fields refine the
/// concordance weight (see [`obs_weight`]). A source that does not carry them (a chip, a tree
/// placement) leaves them `None` or `1.0`, and falls back to the plain source-type weight.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusObs {
    pub name: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub state: ConsensusState,
    /// The observed base (allele) this source called, when known. Carried through to
    /// [`SourceObs::base`] so the state can be re-imputed later without the BAM. `None` keeps the
    /// supplied `state` authoritative (sources that carry only a state, e.g. private calls).
    pub base: Option<String>,
    /// Whether this variant is a known tree/reference variant (true for placement SNPs, false for
    /// private calls).
    pub in_tree: bool,
    /// Read depth at the call (sequencing sources): a `√depth/10` bonus, with a cap of +1.0.
    pub depth: Option<u32>,
    /// Mean mapping quality: an `MQ/60` factor, with a cap of 1.0.
    pub mapq: Option<f64>,
    /// How callable the position is. It scales the weight (`NoCoverage`/`RefN` → 0).
    pub callable: Option<CallableState>,
    /// Region-confidence modifier (e.g. <1 in palindrome/amplicon zones), clamped [0.1, 1.0].
    pub region_modifier: f64,
}

impl ConsensusObs {
    /// A SNP or variant observation with no quality data for the call, so the weight is the
    /// source-type weight. A caller can set the quality fields after this, for a source that
    /// carries them (for example sequencing depth). `base` stays `None`. For an observation that
    /// carries its called allele, use [`ConsensusObs::observed`].
    pub fn snp(
        name: impl Into<String>,
        position: i64,
        ancestral: impl Into<String>,
        derived: impl Into<String>,
        state: ConsensusState,
        in_tree: bool,
    ) -> Self {
        ConsensusObs {
            name: name.into(),
            position,
            ancestral: ancestral.into(),
            derived: derived.into(),
            state,
            base: None,
            in_tree,
            depth: None,
            mapq: None,
            callable: None,
            region_modifier: 1.0,
        }
    }

    /// A SNP or variant observation with the **observed base**. [`impute_state`] makes the state
    /// from that base against the polarity of the variant, and the base stays for a later
    /// imputation ([`reproject`]). `base = None` means a no-call (`NoCall`).
    pub fn observed(
        name: impl Into<String>,
        position: i64,
        ancestral: impl Into<String>,
        derived: impl Into<String>,
        base: Option<char>,
        in_tree: bool,
    ) -> Self {
        let ancestral = ancestral.into();
        let derived = derived.into();
        let state = impute_state(base, &ancestral, &derived);
        ConsensusObs {
            name: name.into(),
            position,
            ancestral,
            derived,
            state,
            base: base.map(|b| b.to_string()),
            in_tree,
            depth: None,
            mapq: None,
            callable: None,
            region_modifier: 1.0,
        }
    }
}

/// Watson-Crick complement of a single base (non-ACGT passes through unchanged).
fn complement_base(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// True when the two alleles of a SNP are strand-ambiguous (`A↔T` / `C↔G`). The complement of one
/// allele is the other, so the observed base does not give the strand.
fn strand_ambiguous(a: char, d: char) -> bool {
    let mut pair = [a.to_ascii_uppercase(), d.to_ascii_uppercase()];
    pair.sort_unstable();
    pair == ['A', 'T'] || pair == ['C', 'G']
}

/// Sentinel observed "base" for an indel locus the sample **carries** (derived). The indel
/// genotyper (`navigator_analysis::caller::call_indels_at`) writes it. Mirrors
/// `navigator_analysis::haplo::INDEL_DERIVED`.
pub const INDEL_DERIVED: char = '+';
/// Sentinel for an indel locus the sample does not carry (ancestral). Mirrors `haplo::INDEL_ANCESTRAL`.
pub const INDEL_ANCESTRAL: char = '-';

/// Make a [`ConsensusState`] from an observed `base`, against the `ancestral` and `derived` alleles
/// of a variant. This is the canonical projection that turns a stored base back into derived or
/// ancestral. Genotyping applies it ([`ConsensusObs::observed`]), and [`reproject`] applies it again
/// against a corrected polarity. It accepts the strand-complement of the alleles, because some
/// trees record a SNP on the strand opposite to the reference. For a strand-ambiguous SNP it keeps
/// a literal match. A base that matches neither strand of either allele, and no base at all, is
/// `NoCall`.
///
/// Mirrors `navigator_analysis::haplo::locus_state`, which works on the analysis `CallState` and
/// `Locus` types. Keep the two in step.
pub fn impute_state(base: Option<char>, ancestral: &str, derived: &str) -> ConsensusState {
    // Indel or MNP (an allele of more than one character). One *base* can not evaluate it. But the
    // indel genotyper resolves it and passes its verdict as a sentinel, so obey that first.
    match base {
        Some(INDEL_DERIVED) => return ConsensusState::Derived,
        Some(INDEL_ANCESTRAL) => return ConsensusState::Ancestral,
        _ => {}
    }
    // If not, nothing can evaluate a multi-base allele from a raw base observation. An insertion
    // or a deletion shares its anchor base, so a compare of the first base would read every sample
    // as derived.
    if ancestral.chars().count() > 1 || derived.chars().count() > 1 {
        return ConsensusState::NoCall;
    }
    let Some(d) = derived.chars().next().map(|c| c.to_ascii_uppercase()) else {
        return ConsensusState::NoCall;
    };
    let a = ancestral.chars().next().map(|c| c.to_ascii_uppercase());
    let Some(b) = base.map(|c| c.to_ascii_uppercase()) else {
        return ConsensusState::NoCall;
    };
    if b == d {
        return ConsensusState::Derived;
    }
    if Some(b) == a {
        return ConsensusState::Ancestral;
    }
    let ambiguous = a.is_some_and(|a| strand_ambiguous(a, d));
    if !ambiguous {
        let bc = complement_base(b);
        if bc == d {
            return ConsensusState::Derived;
        }
        if Some(bc) == a {
            return ConsensusState::Ancestral;
        }
    }
    ConsensusState::NoCall
}

/// Concordance weight for one observation (Scala `YVariantConcordance.calculateWeight`):
/// `snp_weight · (1 + min(√depth/10, 1)) · min(MQ/60, 1) · callableWeight · clamp(region, 0.1, 1)`.
/// Missing depth → no bonus; missing MQ/callable → factor 1.0.
pub fn obs_weight(
    source_type: SourceType,
    depth: Option<u32>,
    mapq: Option<f64>,
    callable: Option<CallableState>,
    region_modifier: f64,
) -> f64 {
    let method = source_type.snp_weight();
    let depth_bonus = depth
        .filter(|&d| d > 0)
        .map(|d| ((d as f64).sqrt() / 10.0).min(1.0))
        .unwrap_or(0.0);
    let mapq_factor = mapq.filter(|&q| q > 0.0).map(|q| (q / 60.0).min(1.0)).unwrap_or(1.0);
    let callable_factor = callable.map(|c| c.weight()).unwrap_or(1.0);
    let region_factor = region_modifier.clamp(0.1, 1.0);
    method * (1.0 + depth_bonus) * mapq_factor * callable_factor * region_factor
}

/// The weighted share of support that does not agree, above which a variant is a conflict.
const CONFLICT_FRACTION: f64 = 0.30;
/// Consensus confidence at or above which a variant from more than one source, with no conflict,
/// counts as confirmed.
const CONFIRMATION_FRACTION: f64 = 0.70;

/// Key a variant, so that it groups across sources and across builds. Use the name when there is
/// one, which is independent of the build. If not, use the position: a novel or unnamed call only
/// ever matches the same position on the same build.
fn group_key(obs: &ConsensusObs) -> String {
    if obs.name.trim().is_empty() {
        format!("@{}", obs.position)
    } else {
        obs.name.trim().to_uppercase()
    }
}

/// Normalize the strand of an observed base into the allele space of this SNP. If the base matches
/// an allele, keep it. If the SNP is not strand-ambiguous, and the complement of the base matches
/// an allele, use the complement. That is a read on the opposite strand. If neither, keep the
/// base as it is: it is a genuine third allele, it survives the vote as itself, and nothing drops
/// it.
///
/// This compares the first base of each allele, because a SNP is one base. An indel allele falls
/// through to a literal compare.
fn canonicalize_base(base: char, ancestral: &str, derived: &str) -> String {
    let b = base.to_ascii_uppercase();
    let a = ancestral.chars().next().map(|c| c.to_ascii_uppercase());
    let d = derived.chars().next().map(|c| c.to_ascii_uppercase());
    if Some(b) == a || Some(b) == d {
        return b.to_string();
    }
    let ambiguous = matches!((a, d), (Some(a), Some(d)) if strand_ambiguous(a, d));
    if !ambiguous {
        let bc = complement_base(b);
        if Some(bc) == a || Some(bc) == d {
            return bc.to_string();
        }
    }
    b.to_string()
}

/// The voted outcome over the **observed bases** of one variant, from each source.
struct BaseTally {
    /// The weighted-majority base (already strand-normalized), or `None` when no source called.
    consensus_base: Option<String>,
    /// Sources whose base equals the consensus base.
    support: usize,
    /// Sources with a call (base present).
    total: usize,
    /// The weight of the base that won, divided by the total weight.
    confidence_score: f64,
}

/// Weight-vote the **canonical bases** of a variant, from each source, into a consensus base. This
/// is the observation-first core. The consensus is the real nucleotide the sources agree on, which
/// is any of A/C/G/T, and a third allele too. It is not a collapse into a binary derived or
/// ancestral value. [`impute_state`] makes the state afterward, from the consensus base against the
/// tree.
fn tally_bases(obs: &[(Option<String>, f64)]) -> BaseTally {
    let mut weights: BTreeMap<String, f64> = BTreeMap::new();
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut total = 0usize;
    for (base, weight) in obs {
        if let Some(base) = base {
            *weights.entry(base.clone()).or_default() += *weight;
            *counts.entry(base.clone()).or_default() += 1;
            total += 1;
        }
    }
    if total == 0 {
        return BaseTally {
            consensus_base: None,
            support: 0,
            total: 0,
            confidence_score: 0.0,
        };
    }
    // Argmax by weight. A tie goes to the base with more raw sources behind it, then to a stable
    // lexical order.
    let consensus_base = weights
        .iter()
        .max_by(|a, b| {
            a.1.partial_cmp(b.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| counts[a.0].cmp(&counts[b.0]))
                .then_with(|| b.0.cmp(a.0))
        })
        .map(|(k, _)| k.clone())
        .expect("total > 0 implies a winner");
    let total_weight: f64 = weights.values().sum();
    let confidence_score = if total_weight > 0.0 {
        weights[&consensus_base] / total_weight
    } else {
        0.0
    };
    let support = counts[&consensus_base];
    BaseTally {
        consensus_base: Some(consensus_base),
        support,
        total,
        confidence_score,
    }
}

/// The cross-source status of a variant given its consensus state, tree membership, coverage, and
/// agreement. Shared taxonomy with the Scala `YVariantConcordance`.
fn status_of(state: ConsensusState, in_tree: bool, total: usize, confidence_score: f64) -> ConsensusStatus {
    let minority_fraction = 1.0 - confidence_score;
    if total == 0 {
        ConsensusStatus::NoCoverage
    } else if minority_fraction > CONFLICT_FRACTION {
        ConsensusStatus::Conflict
    } else if state == ConsensusState::Derived && !in_tree {
        // A derived call that is off the tree is novel or private, even from one source, which is
        // the common case.
        ConsensusStatus::Novel
    } else if total == 1 {
        ConsensusStatus::SingleSource
    } else if confidence_score >= CONFIRMATION_FRACTION {
        ConsensusStatus::Confirmed
    } else {
        ConsensusStatus::Pending
    }
}

/// Group the [`ConsensusObs`] of each source into an [`ObservedProfile`], which is the persisted
/// form that holds observations only. It groups by name, which is independent of the build, and by
/// position when there is no name. It keeps the observed base and the quality of each source. It
/// does NOT store the state, because [`interpret`] makes that on read. The ancestral and derived
/// alleles of the representative become the `ref_allele` and `alt_allele` polarity fallback of the
/// variant. That fallback serves a novel call off the tree, and an mtDNA mutation the tree map does
/// not hold.
pub fn to_observed(sources: &[(String, SourceType, Vec<ConsensusObs>)]) -> ObservedProfile {
    struct Acc {
        repr: ConsensusObs,
        sources: Vec<ObservedSource>,
    }
    let mut groups: BTreeMap<String, Acc> = BTreeMap::new();
    for (label, source_type, observations) in sources {
        for o in observations {
            let key = group_key(o);
            let acc = groups.entry(key).or_insert_with(|| Acc {
                repr: o.clone(),
                sources: Vec::new(),
            });
            if acc.repr.name.trim().is_empty() && !o.name.trim().is_empty() {
                acc.repr = o.clone();
            }
            acc.sources.push(ObservedSource {
                label: label.clone(),
                source_type: *source_type,
                base: o.base.clone(),
                depth: o.depth,
                mapq: o.mapq,
                callable: o.callable,
                region_modifier: o.region_modifier,
            });
        }
    }
    let variants = groups
        .into_values()
        .map(|acc| ObservedVariant {
            name: acc.repr.name,
            position: acc.repr.position,
            in_tree: acc.repr.in_tree,
            ref_allele: Some(acc.repr.ancestral),
            alt_allele: Some(acc.repr.derived),
            sources: acc.sources,
        })
        .collect();
    let source_summaries = sources
        .iter()
        .map(|(label, source_type, obs)| SourceSummary {
            label: label.clone(),
            source_type: *source_type,
            variant_count: obs.len(),
        })
        .collect();
    ObservedProfile {
        schema_version: schema_v1(),
        variants,
        sources: source_summaries,
        terminal_hint: None,
    }
}

/// Interpret an [`ObservedProfile`] against a `polarity` map (`SNP name → (ancestral, derived)`,
/// for example from the current DecodingUs, FTDNA or rCRS tree). The output is the display view:
/// the reconciled [`ConsensusVariant`] list and the [`ConsensusSummary`]. This is the whole point of
/// observation-first storage. The state, the status, the support and the consensus come fresh from
/// the **observed base** of each source, against the **current** polarity. So a corrected tree
/// changes every view, and nothing genotypes again.
///
/// For each variant, resolve the polarity from the map by the upper-case name. If it is not there,
/// fall back to the stored `ref_allele` and `alt_allele`. That covers a novel or private call, and
/// an mtDNA mutation the map does not hold. [`impute_state`] makes the state of each source from
/// its base, and a source with no base gives `NoCall`. [`obs_weight`] weights it over the persisted
/// quality, and [`tally_states`] then counts it.
pub fn interpret(
    observed: &ObservedProfile,
    polarity: &BTreeMap<String, (String, String)>,
) -> (Vec<ConsensusVariant>, ConsensusSummary) {
    let upper: BTreeMap<String, &(String, String)> =
        polarity.iter().map(|(k, v)| (k.trim().to_uppercase(), v)).collect();

    let mut out: Vec<ConsensusVariant> = observed
        .variants
        .iter()
        .map(|v| {
            // Polarity: the current tree by name, else the stored ref/alt fallback.
            let (ancestral, derived) = upper
                .get(v.name.trim().to_uppercase().as_str())
                .map(|(a, d)| ((*a).clone(), (*d).clone()))
                .unwrap_or_else(|| {
                    (
                        v.ref_allele.clone().unwrap_or_default(),
                        v.alt_allele.clone().unwrap_or_default(),
                    )
                });
            let mut source_obs = Vec::with_capacity(v.sources.len());
            let mut bases = Vec::with_capacity(v.sources.len());
            for s in &v.sources {
                let base = s.base.as_deref().and_then(|b| b.chars().next());
                // Vote the real nucleotide, with the strand normalized to the alleles of this SNP.
                // Do not collapse it into a binary derived or ancestral value, so that a
                // multiallelic call and a third-allele call survive.
                let canonical = base.map(|b| canonicalize_base(b, &ancestral, &derived));
                let weight = obs_weight(s.source_type, s.depth, s.mapq, s.callable, s.region_modifier);
                bases.push((canonical, weight));
                source_obs.push(SourceObs {
                    label: s.label.clone(),
                    source_type: s.source_type,
                    state: impute_state(base, &ancestral, &derived),
                    base: s.base.clone(),
                });
            }
            let t = tally_bases(&bases);
            // The state is the interpretation of the consensus *base* against the tree polarity.
            let consensus_char = t.consensus_base.as_deref().and_then(|b| b.chars().next());
            let state = impute_state(consensus_char, &ancestral, &derived);
            let status = status_of(state, v.in_tree, t.total, t.confidence_score);
            ConsensusVariant {
                name: v.name.clone(),
                position: v.position,
                ancestral,
                derived,
                consensus_base: t.consensus_base,
                consensus: state,
                status,
                support: t.support,
                total: t.total,
                in_tree: v.in_tree,
                confidence_score: t.confidence_score,
                sources: source_obs,
            }
        })
        .collect();

    // Conflicts first (most actionable), then novel, then by name.
    out.sort_by(|a, b| {
        status_rank(a.status)
            .cmp(&status_rank(b.status))
            .then_with(|| a.name.cmp(&b.name))
    });
    let summary = summarize(&out);
    (out, summary)
}

/// Reconcile the variant observations of each source into the display view. This is a convenience
/// wrapper: it groups into an [`ObservedProfile`], then it runs [`interpret`] against the **own**
/// stored polarity of each variant. An empty map gives the `ancestral` and `derived` of the
/// observations. New code that persists must call [`to_observed`] and interpret against the current
/// tree, so that a polarity fix reaches everything. [`impute_state`] makes the state of each source
/// from its observed base, and an observation with no base is a `NoCall`.
pub fn reconcile(sources: &[(String, SourceType, Vec<ConsensusObs>)]) -> Vec<ConsensusVariant> {
    interpret(&to_observed(sources), &BTreeMap::new()).0
}

fn status_rank(s: ConsensusStatus) -> u8 {
    match s {
        ConsensusStatus::Conflict => 0,
        ConsensusStatus::Novel => 1,
        ConsensusStatus::Pending => 2,
        ConsensusStatus::SingleSource => 3,
        ConsensusStatus::Confirmed => 4,
        ConsensusStatus::NoCoverage => 5,
    }
}

/// Counts for each status, and the overall confidence, over a reconciled variant list.
pub fn summarize(variants: &[ConsensusVariant]) -> ConsensusSummary {
    let mut s = ConsensusSummary {
        total: variants.len(),
        ..Default::default()
    };
    for v in variants {
        match v.status {
            ConsensusStatus::Confirmed => s.confirmed += 1,
            ConsensusStatus::Novel => s.novel += 1,
            ConsensusStatus::Conflict => s.conflict += 1,
            ConsensusStatus::SingleSource => s.single_source += 1,
            // Pending / NoCoverage are not headline counts; they fold into `total` only.
            ConsensusStatus::Pending | ConsensusStatus::NoCoverage => {}
        }
    }
    // Scala profile confidence: (confirmed + 0.7·novel − 0.5·conflict) / total, clamped [0,1].
    s.overall_confidence = if s.total == 0 {
        0.0
    } else {
        ((s.confirmed as f64 + 0.7 * s.novel as f64 - 0.5 * s.conflict as f64) / s.total as f64).clamp(0.0, 1.0)
    };
    s
}

// ---------------------------------------------------------------------------------------------
// Diploid (autosomal) reconciler. It has the same quality weights, the same status taxonomy, and
// the same summary. But it votes a genotype of three classes (alt-allele dosage 0/1/2), and not a
// binary derived or ancestral state. The autosomal adapter genotypes each source over a fixed site
// panel, and reconciles here.
// ---------------------------------------------------------------------------------------------

/// The diploid call of one source at an autosomal site, for [`reconcile_diploid`]. `dosage` is the
/// alt-allele count 0/1/2, or -1 for a no-call. `depth` drives the weight bonus of the call.
#[derive(Debug, Clone, PartialEq)]
pub struct DiploidObs {
    pub name: String,
    pub contig: String,
    pub position: i64,
    pub reference: String,
    pub alternate: String,
    pub dosage: i8,
    pub depth: Option<u32>,
}

/// One source's diploid observation at a site (for provenance display).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiploidSourceObs {
    pub label: String,
    pub source_type: SourceType,
    pub dosage: i8,
}

/// A reconciled autosomal site over the sources of the subject: a voted diploid genotype.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiploidVariant {
    pub name: String,
    pub contig: String,
    pub position: i64,
    pub reference: String,
    pub alternate: String,
    /// Consensus alt-allele dosage 0/1/2, or -1 when no source made a call.
    pub consensus_dosage: i8,
    pub status: ConsensusStatus,
    /// Sources matching the consensus dosage.
    pub support: usize,
    /// Sources with any call (excludes no-calls).
    pub total: usize,
    #[serde(default)]
    pub confidence_score: f64,
    pub sources: Vec<DiploidSourceObs>,
}

/// Reconcile the diploid genotype calls of each source into one profile. The key is the site name,
/// an rsID, which is independent of the build. Mirrors [`reconcile`], but it votes a genotype
/// of three classes (dosage 0/1/2), and not a binary derived or ancestral state. `Novel` never
/// applies, because every site is a known panel site.
pub fn reconcile_diploid(sources: &[(String, SourceType, Vec<DiploidObs>)]) -> Vec<DiploidVariant> {
    struct ObsRec {
        label: String,
        source_type: SourceType,
        dosage: i8,
        weight: f64,
    }
    struct Acc {
        repr: DiploidObs,
        obs: Vec<ObsRec>,
    }
    let mut groups: BTreeMap<String, Acc> = BTreeMap::new();

    for (label, source_type, observations) in sources {
        for o in observations {
            let key = o.name.trim().to_uppercase(); // rsID — panel sites are always named
            let acc = groups.entry(key).or_insert_with(|| Acc {
                repr: o.clone(),
                obs: Vec::new(),
            });
            // Depth-bonus only (chips have depth 0 → bare method weight; deep WGS earns the bonus).
            let weight = obs_weight(*source_type, o.depth, None, None, 1.0);
            acc.obs.push(ObsRec {
                label: label.clone(),
                source_type: *source_type,
                dosage: o.dosage,
                weight,
            });
        }
    }

    let mut out: Vec<DiploidVariant> = groups
        .into_values()
        .map(|acc| {
            let repr = acc.repr;
            // Weighted vote over the three dosage classes {0,1,2}; no-calls (-1) excluded.
            let mut w = [0.0f64; 3];
            let mut counts = [0usize; 3];
            let mut total = 0usize;
            for o in &acc.obs {
                if (0..=2).contains(&o.dosage) {
                    w[o.dosage as usize] += o.weight;
                    counts[o.dosage as usize] += 1;
                    total += 1;
                }
            }
            // argmax weight. A tie goes to more raw sources behind it, then to the lower dosage.
            let mut best = 0usize;
            for d in 1..3 {
                if w[d] > w[best] || (w[d] == w[best] && counts[d] > counts[best]) {
                    best = d;
                }
            }
            let consensus_dosage: i8 = if total == 0 { -1 } else { best as i8 };

            let total_weight: f64 = w.iter().sum();
            let confidence_score = if total_weight > 0.0 {
                w[best] / total_weight
            } else {
                0.0
            };
            let minority_fraction = 1.0 - confidence_score;

            let status = if total == 0 {
                ConsensusStatus::NoCoverage
            } else if minority_fraction > CONFLICT_FRACTION {
                ConsensusStatus::Conflict
            } else if total == 1 {
                ConsensusStatus::SingleSource
            } else if confidence_score >= CONFIRMATION_FRACTION {
                ConsensusStatus::Confirmed
            } else {
                ConsensusStatus::Pending
            };

            let support = acc
                .obs
                .iter()
                .filter(|o| consensus_dosage >= 0 && o.dosage == consensus_dosage)
                .count();
            let sources = acc
                .obs
                .iter()
                .map(|o| DiploidSourceObs {
                    label: o.label.clone(),
                    source_type: o.source_type,
                    dosage: o.dosage,
                })
                .collect();

            DiploidVariant {
                name: repr.name,
                contig: repr.contig,
                position: repr.position,
                reference: repr.reference,
                alternate: repr.alternate,
                consensus_dosage,
                status,
                support,
                total,
                confidence_score,
                sources,
            }
        })
        .collect();

    out.sort_by(|a, b| {
        status_rank(a.status)
            .cmp(&status_rank(b.status))
            .then_with(|| a.name.cmp(&b.name))
    });
    out
}

/// Counts for each status, and the overall confidence, over a reconciled diploid variant list.
/// `Novel` does not apply to an autosomal site, so the confidence is
/// `(confirmed − 0.5·conflict) / total`.
pub fn summarize_diploid(variants: &[DiploidVariant]) -> ConsensusSummary {
    let mut s = ConsensusSummary {
        total: variants.len(),
        ..Default::default()
    };
    for v in variants {
        match v.status {
            ConsensusStatus::Confirmed => s.confirmed += 1,
            ConsensusStatus::Conflict => s.conflict += 1,
            ConsensusStatus::SingleSource => s.single_source += 1,
            ConsensusStatus::Novel | ConsensusStatus::Pending | ConsensusStatus::NoCoverage => {}
        }
    }
    s.overall_confidence = if s.total == 0 {
        0.0
    } else {
        ((s.confirmed as f64 - 0.5 * s.conflict as f64) / s.total as f64).clamp(0.0, 1.0)
    };
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build an observation with a base that agrees with the wanted state (anc=A, der=G). The
    // observation-first path (`to_observed` → `interpret`) then derives that state from the base.
    fn obs(name: &str, pos: i64, state: ConsensusState, in_tree: bool) -> ConsensusObs {
        let base = match state {
            ConsensusState::Derived => Some('G'),
            ConsensusState::Ancestral => Some('A'),
            ConsensusState::NoCall => None,
        };
        ConsensusObs::observed(name, pos, "A", "G", base, in_tree)
    }

    #[test]
    fn impute_state_literal_complement_and_palindrome() {
        // Literal matches.
        assert_eq!(impute_state(Some('G'), "A", "G"), ConsensusState::Derived);
        assert_eq!(impute_state(Some('A'), "A", "G"), ConsensusState::Ancestral);
        // A read on the opposite strand matches through the complement (non-ambiguous A>C: comp
        // T/G).
        assert_eq!(impute_state(Some('G'), "A", "C"), ConsensusState::Derived); // comp(G)=C=derived
        assert_eq!(impute_state(Some('T'), "A", "C"), ConsensusState::Ancestral); // comp(T)=A=ancestral

        // Strand-ambiguous C/G: complement of derived G is ancestral C → keep literal only.
        assert_eq!(impute_state(Some('C'), "C", "G"), ConsensusState::Ancestral);
        assert_eq!(impute_state(Some('A'), "C", "G"), ConsensusState::NoCall); // genuine third allele

        // No base → no call.
        assert_eq!(impute_state(None, "A", "G"), ConsensusState::NoCall);
    }

    #[test]
    fn impute_state_indel_is_nocall() {
        // An indel shares its anchor base between the alleles (G against GAGC), so one observed
        // base can not evaluate it. It must be a no-call, and not a false derived.
        assert_eq!(impute_state(Some('G'), "G", "GAGC"), ConsensusState::NoCall); // insertion
        assert_eq!(impute_state(Some('G'), "GAGC", "G"), ConsensusState::NoCall); // deletion
        assert_eq!(impute_state(Some('A'), "AT", "GC"), ConsensusState::NoCall);
        // MNP
    }

    #[test]
    fn observed_constructor_imputes_and_keeps_base() {
        let o = ConsensusObs::observed("PF1016", 100, "C", "T", Some('T'), true);
        assert_eq!(o.state, ConsensusState::Derived);
        assert_eq!(o.base.as_deref(), Some("T"));
    }

    #[test]
    fn interpret_flips_state_against_corrected_polarity_from_stored_base() {
        // One source observed base T. Against an FTDNA-style inverted polarity (anc=T, der=C) it
        // reads Ancestral. Against the true DecodingUs polarity (anc=C, der=T) the SAME stored base
        // reads Derived. The code calculates this live, and nothing genotypes again. The consensus
        // *base* is T in both, and only the reading of it changes.
        let observed = to_observed(&[(
            "aln #1".into(),
            SourceType::WgsShortRead,
            vec![ConsensusObs::observed("PF1016", 100, "T", "C", Some('T'), true)],
        )]);

        // Empty map → falls back to the stored ref/alt polarity (T>C) → Ancestral.
        let (v0, _) = interpret(&observed, &BTreeMap::new());
        assert_eq!(v0[0].consensus, ConsensusState::Ancestral);
        assert_eq!(v0[0].consensus_base.as_deref(), Some("T"));

        // Corrected polarity C>T → the same base is now Derived.
        let polarity: BTreeMap<String, (String, String)> = [("PF1016".to_string(), ("C".to_string(), "T".to_string()))]
            .into_iter()
            .collect();
        let (v1, _) = interpret(&observed, &polarity);
        assert_eq!(v1[0].consensus, ConsensusState::Derived);
        assert_eq!(v1[0].consensus_base.as_deref(), Some("T"));
        assert_eq!(v1[0].ancestral, "C");
        assert_eq!(v1[0].derived, "T");
        assert_eq!(v1[0].sources[0].state, ConsensusState::Derived);
    }

    #[test]
    fn multiallelic_third_allele_survives_the_vote() {
        // At an A>G SNP, two sources read a genuine third allele T on the *forward* strand. It is
        // not the A allele and not the G allele, and comp(T)=A is ancestral. So the code treats T
        // as an ancestral read on the opposite strand. A cleaner third-allele case is next: a
        // strand-ambiguous A/T with a C read stays C.
        let observed = to_observed(&[
            (
                "a".into(),
                SourceType::WgsShortRead,
                vec![ConsensusObs::observed("S1", 1, "A", "T", Some('C'), true)],
            ),
            (
                "b".into(),
                SourceType::WgsShortRead,
                vec![ConsensusObs::observed("S1", 1, "A", "T", Some('C'), true)],
            ),
        ]);
        let (v, _) = interpret(&observed, &BTreeMap::new());
        // A/T is strand-ambiguous, so a C read matches no allele, and it stays as itself. The
        // consensus base is the real third allele C, and nothing folds it away.
        assert_eq!(v[0].consensus_base.as_deref(), Some("C"));
        assert_eq!(v[0].consensus, ConsensusState::NoCall); // C is neither ancestral nor derived
        assert_eq!(v[0].total, 2);
    }

    #[test]
    fn to_observed_preserves_per_call_quality() {
        // Depth and region go into the stored observation, so that interpret can weight exactly.
        let mut o = ConsensusObs::observed("M269", 100, "A", "G", Some('G'), true);
        o.depth = Some(100);
        o.region_modifier = 0.4;
        let observed = to_observed(&[("aln".into(), SourceType::WgsShortRead, vec![o])]);
        let s = &observed.variants[0].sources[0];
        assert_eq!(s.depth, Some(100));
        assert!((s.region_modifier - 0.4).abs() < 1e-9);
    }

    #[test]
    fn obs_weight_applies_depth_mapq_callable() {
        // No quality data → bare source-type weight.
        assert!((obs_weight(SourceType::WgsShortRead, None, None, None, 1.0) - 0.85).abs() < 1e-9);
        // depth 100 → bonus min(√100/10,1)=1.0 → ×2; MQ 60 → ×1; callable → ×1.
        assert!(
            (obs_weight(
                SourceType::WgsShortRead,
                Some(100),
                Some(60.0),
                Some(CallableState::Callable),
                1.0
            ) - 1.7)
                .abs()
                < 1e-9
        );
        // Low coverage halves; a region modifier <1 scales down further.
        let w = obs_weight(
            SourceType::WgsShortRead,
            None,
            None,
            Some(CallableState::LowCoverage),
            0.5,
        );
        assert!((w - 0.85 * 0.5 * 0.5).abs() < 1e-9);
        // NoCoverage callability zeroes the weight.
        assert_eq!(
            obs_weight(
                SourceType::Sanger,
                Some(50),
                Some(60.0),
                Some(CallableState::NoCoverage),
                1.0
            ),
            0.0
        );
    }

    #[test]
    fn confidence_score_and_overall() {
        let v = reconcile(&[
            (
                "a".into(),
                SourceType::WgsShortRead,
                vec![obs("M269", 1, ConsensusState::Derived, true)],
            ),
            (
                "b".into(),
                SourceType::Chip,
                vec![obs("M269", 1, ConsensusState::Derived, true)],
            ),
        ]);
        assert!((v[0].confidence_score - 1.0).abs() < 1e-9); // unanimous → full confidence
        let s = summarize(&v);
        assert!((s.overall_confidence - 1.0).abs() < 1e-9); // 1 confirmed / 1 total
    }

    #[test]
    fn two_sources_agree_in_tree_is_confirmed() {
        let v = reconcile(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![obs("M269", 100, ConsensusState::Derived, true)],
            ),
            (
                "consumer".into(),
                SourceType::Chip,
                vec![obs("M269", 200, ConsensusState::Derived, true)],
            ),
        ]);
        assert_eq!(v.len(), 1); // grouped by name across differing positions/builds
        assert_eq!(v[0].name, "M269");
        assert_eq!(v[0].consensus, ConsensusState::Derived);
        assert_eq!(v[0].status, ConsensusStatus::Confirmed);
        assert_eq!(v[0].support, 2);
        assert_eq!(v[0].total, 2);
    }

    #[test]
    fn derived_not_in_tree_is_novel() {
        let v = reconcile(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![obs("FT1", 100, ConsensusState::Derived, false)],
            ),
            (
                "aln #2".into(),
                SourceType::WgsShortRead,
                vec![obs("FT1", 100, ConsensusState::Derived, false)],
            ),
        ]);
        assert_eq!(v[0].status, ConsensusStatus::Novel);
    }

    #[test]
    fn comparable_weight_disagreement_is_conflict() {
        // WGS (0.85) derived vs Chip (0.5) ancestral → minority 0.5/1.35 ≈ 0.37 > 0.30 → conflict.
        let v = reconcile(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![obs("M269", 100, ConsensusState::Derived, true)],
            ),
            (
                "consumer".into(),
                SourceType::Chip,
                vec![obs("M269", 100, ConsensusState::Ancestral, true)],
            ),
        ]);
        assert_eq!(v[0].status, ConsensusStatus::Conflict);
        assert_eq!(v[0].consensus, ConsensusState::Derived); // higher weight wins the consensus
    }

    #[test]
    fn dominant_weight_disagreement_is_not_conflict() {
        // Sanger (1.0) derived vs Manual (0.3) ancestral → minority 0.3/1.3 ≈ 0.23 ≤ 0.30 → confirmed.
        let v = reconcile(&[
            (
                "sanger".into(),
                SourceType::Sanger,
                vec![obs("M269", 100, ConsensusState::Derived, true)],
            ),
            (
                "manual".into(),
                SourceType::Manual,
                vec![obs("M269", 100, ConsensusState::Ancestral, true)],
            ),
        ]);
        assert_eq!(v[0].consensus, ConsensusState::Derived);
        assert_eq!(v[0].status, ConsensusStatus::Confirmed);
        assert_eq!(v[0].support, 1); // only the Sanger source matches the derived consensus
    }

    #[test]
    fn single_source_is_single_source() {
        let v = reconcile(&[(
            "aln #1".into(),
            SourceType::WgsShortRead,
            vec![obs("M269", 100, ConsensusState::Derived, true)],
        )]);
        assert_eq!(v[0].status, ConsensusStatus::SingleSource);
    }

    #[test]
    fn nocall_excluded_from_vote() {
        let v = reconcile(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![obs("M269", 100, ConsensusState::Derived, true)],
            ),
            (
                "aln #2".into(),
                SourceType::WgsShortRead,
                vec![obs("M269", 100, ConsensusState::NoCall, true)],
            ),
        ]);
        assert_eq!(v[0].total, 1); // NoCall not counted
        assert_eq!(v[0].status, ConsensusStatus::SingleSource);
        assert_eq!(v[0].sources.len(), 2); // but still shown for provenance
    }

    #[test]
    fn summary_counts_by_status() {
        let v = reconcile(&[
            (
                "a".into(),
                SourceType::WgsShortRead,
                vec![
                    obs("M269", 1, ConsensusState::Derived, true),
                    obs("FT1", 2, ConsensusState::Derived, false),
                ],
            ),
            (
                "b".into(),
                SourceType::Chip,
                vec![obs("M269", 1, ConsensusState::Derived, true)],
            ),
        ]);
        let s = summarize(&v);
        assert_eq!(s.total, 2);
        assert_eq!(s.confirmed, 1); // M269 (2 sources agree, in tree)
        assert_eq!(s.novel, 1); // FT1 (derived, not in tree → novel even single-source)
        assert_eq!(s.single_source, 0);
    }

    fn dobs(name: &str, dosage: i8, depth: Option<u32>) -> DiploidObs {
        DiploidObs {
            name: name.into(),
            contig: "chr1".into(),
            position: 100,
            reference: "A".into(),
            alternate: "G".into(),
            dosage,
            depth,
        }
    }

    #[test]
    fn diploid_two_sources_agree_is_confirmed() {
        let v = reconcile_diploid(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![dobs("rs1", 1, Some(30))],
            ),
            ("chip".into(), SourceType::Chip, vec![dobs("rs1", 1, None)]),
        ]);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].consensus_dosage, 1); // both het
        assert_eq!(v[0].status, ConsensusStatus::Confirmed);
        assert!((v[0].confidence_score - 1.0).abs() < 1e-9);
        assert_eq!((v[0].support, v[0].total), (2, 2));
    }

    #[test]
    fn diploid_comparable_weight_disagreement_is_conflict() {
        // Two equal-weight WGS sources, hom-ref vs hom-alt → minority 0.5 > 0.30 → conflict; the
        // tie resolves to the lower dosage.
        let v = reconcile_diploid(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![dobs("rs1", 0, None)]),
            ("aln #2".into(), SourceType::WgsShortRead, vec![dobs("rs1", 2, None)]),
        ]);
        assert_eq!(v[0].status, ConsensusStatus::Conflict);
        assert_eq!(v[0].consensus_dosage, 0);
    }

    #[test]
    fn diploid_single_source_and_nocall() {
        // One het call + one no-call → counted as a single source (no-call excluded from the vote),
        // but the no-call is still shown for provenance.
        let v = reconcile_diploid(&[
            (
                "aln #1".into(),
                SourceType::WgsShortRead,
                vec![dobs("rs1", 1, Some(30))],
            ),
            ("chip".into(), SourceType::Chip, vec![dobs("rs1", -1, None)]),
        ]);
        assert_eq!(v[0].total, 1);
        assert_eq!(v[0].status, ConsensusStatus::SingleSource);
        assert_eq!(v[0].consensus_dosage, 1);
        assert_eq!(v[0].sources.len(), 2);
    }

    #[test]
    fn diploid_summary_counts_and_confidence() {
        let v = reconcile_diploid(&[
            (
                "a".into(),
                SourceType::WgsShortRead,
                vec![dobs("rs1", 2, None), dobs("rs2", 0, None)],
            ),
            (
                "b".into(),
                SourceType::WgsShortRead,
                vec![dobs("rs1", 2, None), dobs("rs2", 2, None)],
            ),
        ]);
        let s = summarize_diploid(&v);
        assert_eq!(s.total, 2);
        assert_eq!(s.confirmed, 1); // rs1 (both hom-alt)
        assert_eq!(s.conflict, 1); // rs2 (0 vs 2)

        // (1 confirmed − 0.5·1 conflict) / 2 = 0.25
        assert!((s.overall_confidence - 0.25).abs() < 1e-9);
    }
}

//! Multi-source variant **consensus engine** — DNA-type-agnostic.
//!
//! Given a set of sources (a WGS alignment's placement, a chip/BISDNA panel, a private bucket, …),
//! each contributing per-variant calls keyed **by name** (build-independent — M269 is M269 whether
//! the source aligned to GRCh37 or GRCh38), [`reconcile`] groups them and weight-votes the consensus
//! state, classifying each variant as confirmed / novel / conflict / single-source and computing a
//! quality-weighted confidence. Mirrors the Scala `YVariantConcordance`.
//!
//! This engine is the shared foundation for the Y-DNA profile (the [`crate::yprofile`] adapter today)
//! and — by design — the future mtDNA (variants vs rCRS) and autosomal consumers. It carries no
//! DNA-type specifics: callers gather observations and supply the variant identity; the DNA type and
//! consensus label (haplogroup, where applicable) live at the persistence / app layer.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::variants::SourceType;

/// One source's call state at a variant position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusState {
    /// Carries the derived (mutant) allele — positive for the variant's branch. For mtDNA this is
    /// "differs from rCRS"; for autosomes a future adapter maps a diploid genotype onto this axis.
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
    /// Has data but the weighted confidence is below the confirmation threshold without crossing the
    /// conflict line (rare — kept for parity with the Scala `YVariantConcordance`).
    Pending,
    /// No source made a confident call (every observation was NoCall).
    NoCoverage,
}

/// Per-position callability of a source's observation — scales its concordance weight (a base in a
/// no-coverage / poor-mapping region carries little confidence). Mirrors the Scala `YCallableState`.
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
    /// The **observed base** (allele) this source called at the variant — `None` for a no-call, or
    /// for sources/legacy profiles that carry only a state. Persisting the base (not just the
    /// derived/ancestral interpretation) lets the state be re-[`impute_state`]d against a corrected
    /// or different tree polarity via [`reproject`] — without re-reading the BAM/CRAM.
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
    /// The **consensus observed base** — the weighted-majority nucleotide across sources (strand-
    /// normalized to this SNP's alleles). This is the primary observation; [`consensus`](Self::consensus)
    /// is its derived/ancestral interpretation against the tree. `None` = no source made a call. A base
    /// matching neither allele (a genuine third allele) survives here as itself.
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

/// Per-status counts for the profile header.
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
// Observation-first storage. A persisted profile holds only OBSERVATIONS — per-SNP, per-source
// observed bases + quality + identity — never a baked derived/ancestral interpretation. The state,
// vote, status, support, and summary are computed on demand by [`interpret`] against the CURRENT
// tree's polarity, so a tree-polarity fix (or provider switch) corrects every view with no
// re-genotyping. This is the type actually written to the `consensus_profile` payload.
// ---------------------------------------------------------------------------------------------

fn one() -> f64 {
    1.0
}
fn schema_v1() -> u8 {
    1
}

/// One source's raw observation of a variant — the **observed base** plus the quality inputs to the
/// concordance weight. Carries no derived/ancestral state: that is [`impute_state`]d at read time.
/// (Persisting depth/MQ/callable/region — which `reconcile`'s in-memory tally used but never stored —
/// lets [`interpret`] re-weight exactly, fixing the quality-loss the old `reproject` warned about.)
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

/// A variant observed across the subject's sources — identity + each source's observed base. The
/// derived/ancestral polarity comes from the current tree at [`interpret`] time (by name); for an
/// off-tree novel/private call — and for mtDNA mutations absent from the tree map — the stored
/// `ref_allele`/`alt_allele` are the polarity fallback.
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

/// One contributing source's provenance (label, type, count) — non-interpretive, carried for display.
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
    /// Payload schema tag — presence distinguishes this from a legacy baked `ConsensusProfile` JSON.
    #[serde(default = "schema_v1")]
    pub schema_version: u8,
    pub variants: Vec<ObservedVariant>,
    #[serde(default)]
    pub sources: Vec<SourceSummary>,
    /// The placement's terminal haplogroup label (a placement output, not a per-SNP interpretation).
    #[serde(default)]
    pub terminal_hint: Option<String>,
}

/// One source's call at a variant, fed into [`reconcile`]. Quality fields refine the concordance
/// weight (see [`obs_weight`]); sources that don't carry them (chip, tree placement) leave them
/// `None` / `1.0` and fall back to the plain source-type weight.
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
    /// Read depth at the call (sequencing sources) — a `√depth/10` bonus, capped at +1.0.
    pub depth: Option<u32>,
    /// Mean mapping quality — an `MQ/60` factor, capped at 1.0.
    pub mapq: Option<f64>,
    /// Callability of the position — scales the weight (`NoCoverage`/`RefN` → 0).
    pub callable: Option<CallableState>,
    /// Region-confidence modifier (e.g. <1 in palindrome/amplicon zones), clamped [0.1, 1.0].
    pub region_modifier: f64,
}

impl ConsensusObs {
    /// A SNP/variant observation with no per-call quality data (weight = the source-type weight).
    /// Quality fields can be set afterward for sources that carry them (e.g. sequencing depth).
    /// `base` is left `None`; for an observation that carries its called allele use
    /// [`ConsensusObs::observed`].
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

    /// A SNP/variant observation carrying the **observed base**; the state is imputed from the base
    /// against the variant's polarity ([`impute_state`]) and the base is retained for later
    /// re-imputation ([`reproject`]). `base = None` means a no-call (`NoCall`).
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

/// Whether a SNP's two alleles are strand-ambiguous (`A↔T` / `C↔G`): the complement of one allele
/// equals the other, so strand can't be inferred from the observed base.
fn strand_ambiguous(a: char, d: char) -> bool {
    let mut pair = [a.to_ascii_uppercase(), d.to_ascii_uppercase()];
    pair.sort_unstable();
    pair == ['A', 'T'] || pair == ['C', 'G']
}

/// Impute a [`ConsensusState`] from an observed `base` against a variant's `ancestral`/`derived`
/// alleles. The canonical projection that turns a stored base back into derived/ancestral — applied
/// at genotyping time ([`ConsensusObs::observed`]) and re-applied against corrected polarity by
/// [`reproject`]. Accepts the strand-complement of the alleles (some trees record a SNP on the
/// opposite strand from the reference) except for strand-ambiguous SNPs, where literal matching is
/// kept. A base matching neither strand of either allele, or no base, is `NoCall`.
///
/// Mirrors `navigator_analysis::haplo::locus_state` (which operates on the analysis `CallState` /
/// `Locus` types); keep the two in step.
pub fn impute_state(base: Option<char>, ancestral: &str, derived: &str) -> ConsensusState {
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

/// Fraction of disagreeing (weighted) support above which a variant is a conflict.
const CONFLICT_FRACTION: f64 = 0.30;
/// Consensus confidence at or above which a multi-source, non-conflicting variant is confirmed.
const CONFIRMATION_FRACTION: f64 = 0.70;

/// Key a variant for cross-source/cross-build grouping: by name when present (build-independent),
/// else by position (a novel/unnamed call only ever matches the same build's same position).
fn group_key(obs: &ConsensusObs) -> String {
    if obs.name.trim().is_empty() {
        format!("@{}", obs.position)
    } else {
        obs.name.trim().to_uppercase()
    }
}

/// Strand-normalize an observed base to this SNP's allele space: if it matches an allele keep it;
/// if (for a non-strand-ambiguous SNP) its complement matches an allele, use the complement (an
/// opposite-strand read); otherwise keep it as-is — a genuine third allele that survives the vote as
/// itself rather than being discarded. Compares the first base of each allele (SNPs are single-base;
/// indel alleles fall through to a literal compare).
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

/// The voted outcome over one variant's per-source **observed bases**.
struct BaseTally {
    /// The weighted-majority base (already strand-normalized), or `None` when no source called.
    consensus_base: Option<String>,
    /// Sources whose base equals the consensus base.
    support: usize,
    /// Sources with a call (base present).
    total: usize,
    /// Winning base weight / total weight.
    confidence_score: f64,
}

/// Weight-vote a variant's per-source **canonicalized bases** into a consensus base. This is the
/// observation-first core: the consensus is the actual nucleotide the sources agree on (any of
/// A/C/G/T, incl. a third allele), not a binary derived/ancestral collapse — the state is derived
/// afterward by [`impute_state`]ing the consensus base against the tree.
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
    // Argmax by weight; ties broken by more raw supporting sources, then a stable lexical order.
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
        // Derived off-tree call is novel/private — even from a single source (the common case).
        ConsensusStatus::Novel
    } else if total == 1 {
        ConsensusStatus::SingleSource
    } else if confidence_score >= CONFIRMATION_FRACTION {
        ConsensusStatus::Confirmed
    } else {
        ConsensusStatus::Pending
    }
}

/// Group per-source [`ConsensusObs`] into an [`ObservedProfile`] — the persisted, observation-only
/// form. Groups by name (build-independent) else position, keeping each source's observed base +
/// quality; the state is NOT stored (it is [`interpret`]ed on read). The representative's
/// ancestral/derived become the variant's `ref_allele`/`alt_allele` polarity fallback (used for
/// off-tree novel calls and mtDNA mutations absent from the tree map).
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

/// Interpret an [`ObservedProfile`] against a `polarity` map (`SNP name → (ancestral, derived)`, e.g.
/// from the current DecodingUs/FTDNA/rCRS tree) into the display view — the reconciled
/// [`ConsensusVariant`]s + [`ConsensusSummary`]. This is the whole point of observation-first
/// storage: state/status/support/consensus are derived here, fresh, from each source's **observed
/// base** against the **current** polarity — so a corrected tree flips every view with no
/// re-genotyping.
///
/// Per variant: resolve polarity from the map by upper-cased name, else fall back to the stored
/// `ref_allele`/`alt_allele` (novel/private, and mtDNA mutations not in the map). Each source's state
/// is [`impute_state`]d from its base (base-less sources → `NoCall`), weighted by [`obs_weight`] over
/// the persisted quality, then [`tally_states`]d.
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
                // Vote the actual nucleotide (strand-normalized to this SNP's alleles), not a binary
                // derived/ancestral collapse — so multiallelic / third-allele calls survive.
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

/// Reconcile per-source variant observations into the display view. Convenience wrapper: group into
/// an [`ObservedProfile`] then [`interpret`] against each variant's **own** stored polarity (empty
/// map → the observations' `ancestral`/`derived`). New code that persists should call [`to_observed`]
/// and interpret against the current tree, so a polarity fix propagates. Each source's state is
/// imputed from its observed base — a base-less observation is a `NoCall`.
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

/// Per-status counts + overall confidence over a reconciled variant list.
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
            // Pending / NoCoverage aren't headline counts; they fold into `total` only.
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
// Diploid (autosomal) reconciler — the same quality-weighting + status taxonomy + summary, but
// voting a three-class genotype (alt-allele dosage 0/1/2) instead of a binary derived/ancestral
// state. The autosomal adapter genotypes each source over a fixed site panel and reconciles here.
// ---------------------------------------------------------------------------------------------

/// One source's diploid call at an autosomal site, fed into [`reconcile_diploid`]. `dosage` is the
/// alt-allele count 0/1/2, or -1 for a no-call. `depth` drives the per-call weight bonus.
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

/// A reconciled autosomal site across the subject's sources — a voted diploid genotype.
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

/// Reconcile per-source diploid genotype calls into one profile, keyed by site name (rsID —
/// build-independent). Mirrors [`reconcile`] but votes a three-class genotype (dosage 0/1/2)
/// instead of a binary derived/ancestral state. `Novel` never applies — every site is a known
/// panel site.
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
            // argmax weight; tie → more raw supporting sources, then the lower dosage.
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

/// Per-status counts + overall confidence over a reconciled diploid variant list. `Novel` doesn't
/// apply to autosomal sites, so the confidence is `(confirmed − 0.5·conflict) / total`.
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

    // Build an observation carrying a base consistent with the desired state (anc=A, der=G), so the
    // observation-first path (`to_observed` → `interpret`) re-derives that state from the base.
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
        // Opposite-strand reads match via the complement (non-ambiguous A>C: comp T/G).
        assert_eq!(impute_state(Some('G'), "A", "C"), ConsensusState::Derived); // comp(G)=C=derived
        assert_eq!(impute_state(Some('T'), "A", "C"), ConsensusState::Ancestral); // comp(T)=A=ancestral
        // Strand-ambiguous C/G: complement of derived G is ancestral C → keep literal only.
        assert_eq!(impute_state(Some('C'), "C", "G"), ConsensusState::Ancestral);
        assert_eq!(impute_state(Some('A'), "C", "G"), ConsensusState::NoCall); // genuine third allele
        // No base → no call.
        assert_eq!(impute_state(None, "A", "G"), ConsensusState::NoCall);
    }

    #[test]
    fn observed_constructor_imputes_and_keeps_base() {
        let o = ConsensusObs::observed("PF1016", 100, "C", "T", Some('T'), true);
        assert_eq!(o.state, ConsensusState::Derived);
        assert_eq!(o.base.as_deref(), Some("T"));
    }

    #[test]
    fn interpret_flips_state_against_corrected_polarity_from_stored_base() {
        // One source observed base T. Interpreting against an FTDNA-style inverted polarity (anc=T,
        // der=C) reads it Ancestral; against the true DecodingUs polarity (anc=C, der=T) the SAME
        // stored base reads Derived — computed live, no re-genotyping. The consensus *base* is T in
        // both; only its interpretation flips.
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
        let polarity: BTreeMap<String, (String, String)> =
            [("PF1016".to_string(), ("C".to_string(), "T".to_string()))].into_iter().collect();
        let (v1, _) = interpret(&observed, &polarity);
        assert_eq!(v1[0].consensus, ConsensusState::Derived);
        assert_eq!(v1[0].consensus_base.as_deref(), Some("T"));
        assert_eq!(v1[0].ancestral, "C");
        assert_eq!(v1[0].derived, "T");
        assert_eq!(v1[0].sources[0].state, ConsensusState::Derived);
    }

    #[test]
    fn multiallelic_third_allele_survives_the_vote() {
        // At an A>G SNP, two sources read a genuine third allele T on the *forward* strand (not the
        // A/G alleles, and comp(T)=A is ancestral — so T is treated as an opposite-strand ancestral
        // read here). A cleaner third-allele case: strand-ambiguous A/T with a C read stays C.
        let observed = to_observed(&[
            ("a".into(), SourceType::WgsShortRead, vec![ConsensusObs::observed("S1", 1, "A", "T", Some('C'), true)]),
            ("b".into(), SourceType::WgsShortRead, vec![ConsensusObs::observed("S1", 1, "A", "T", Some('C'), true)]),
        ]);
        let (v, _) = interpret(&observed, &BTreeMap::new());
        // A/T is strand-ambiguous, so a C read matches no allele and is kept as itself — the
        // consensus base is the actual third allele C, not folded away.
        assert_eq!(v[0].consensus_base.as_deref(), Some("C"));
        assert_eq!(v[0].consensus, ConsensusState::NoCall); // C is neither ancestral nor derived
        assert_eq!(v[0].total, 2);
    }

    #[test]
    fn to_observed_preserves_per_call_quality() {
        // Depth/region are carried into the stored observation so interpret can weight exactly.
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

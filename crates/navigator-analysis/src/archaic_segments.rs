//! Tier B — archaic **segment** calling (design `documents/design/ArchaicAncestry_Design.md` §5).
//!
//! Where Tier A counts marker copies, this finds the actual introgressed **tracts**: an hmmix-style
//! (Skov et al. 2018) two-state HMM over the density of *private* derived variants — the subject's
//! variants that no African outgroup individual carries. Anything Africans also carry is not
//! evidence of introgression, so stripping them is what makes the remaining density informative.
//!
//! The HMM **can not tell Neanderthal from Denisovan** — the two lineages coalesce before either
//! meets modern humans (§3) — so it finds segments and a downstream pass labels them by counting
//! derived-allele matches against the archaic genomes (`ArchaicClassify`).
//!
//! Log-space Viterbi + forward/backward, cM-scaled transitions: the same idiom as
//! [`crate::roh`] and the chromosome painter. Only the emission differs — a Poisson point process
//! over private-variant counts per window, rather than het/hom.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::archaic::{ArchaicCallable, ArchaicClassify, ArchaicOutgroup, DiagnosticClass};
use crate::caller::SiteGenotype;
use crate::ibd::GeneticMap;

/// Which archaic lineage a called segment was attributed to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArchaicSource {
    Neanderthal,
    Denisovan,
    /// Archaic by density, but the diagnostic sites in it do not favour either lineage — the
    /// honest label for a segment we can not attribute, and a substantial share in real data
    /// (Skov 2020 reported ~12 % unknown on Icelanders).
    Unknown,
}

/// One called archaic tract.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSegment {
    pub contig: String,
    /// 1-based inclusive bounds.
    pub start: i64,
    pub end: i64,
    /// Mean posterior probability of the archaic state across the segment's windows.
    pub posterior: f64,
    /// Private derived variants inside the segment — the evidence the call rests on.
    pub n_private: usize,
    pub source: ArchaicSource,
    /// Diagnostic-site matches supporting each lineage (the basis for `source`).
    pub neanderthal_matches: usize,
    pub denisovan_matches: usize,
}

impl ArchaicSegment {
    pub fn length_mb(&self) -> f64 {
        (self.end - self.start).max(0) as f64 / 1_000_000.0
    }
}

/// Genome-level totals.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSummary {
    pub total_mb: f64,
    /// Archaic share of the **callable** span, as a percentage. Reported against callable rather
    /// than nominal genome length: a partially-covered genome would otherwise read as having less
    /// archaic ancestry purely because it was sequenced less.
    pub pct_callable: f64,
    pub callable_mb: f64,
    pub neanderthal_mb: f64,
    pub denisovan_mb: f64,
    pub unknown_mb: f64,
    pub n_segments: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSegmentResult {
    pub segments: Vec<ArchaicSegment>,
    pub summary: ArchaicSummary,
}

/// Tuning knobs, same shape as `RohConfig` / `PaintParams`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArchaicConfig {
    /// Window size in bp for the density process.
    pub window_bp: i64,
    /// Expected private variants per window outside archaic tracts. `None` = estimate from the
    /// sample's own genome-wide private rate, which adapts to depth and to how aggressively the
    /// outgroup stripped.
    pub rate_background: Option<f64>,
    /// Multiple of the background rate expected inside an archaic tract. Introgressed haplotypes
    /// carry an excess of derived variants absent from Africans; this is that excess.
    pub archaic_rate_multiple: f64,
    /// Prior probability of the archaic state.
    pub prior_archaic: f64,
    /// Expected state switches per centimorgan — the recombination-scaled transition.
    pub switches_per_cm: f64,
    /// Discard called tracts shorter than this.
    pub min_segment_bp: i64,
    /// Discard tracts whose mean posterior is below this.
    pub min_posterior: f64,
    /// A segment is attributed to a lineage only when its **enrichment** over the expected
    /// carrier rate favours that lineage by at least this ratio; otherwise it is `Unknown`.
    pub min_lineage_ratio: f64,
    /// Expected fraction of Neanderthal-diagnostic sites at which a non-archaic-specific genome
    /// carries the derived allele, and the same for Denisovan-diagnostic sites.
    ///
    /// These base rates are the reason raw match counts can not attribute a lineage. Measured on the
    /// ground-truth European: 4.3 % at Neanderthal-diagnostic sites versus 3.9 % at
    /// Denisovan-diagnostic ones — a ratio of 1.10, essentially no discrimination. Carrying a
    /// "Denisovan-diagnostic" allele mostly reflects ordinary shared ancestry, not Denisovan
    /// introgression, so attribution must compare observed against expected rather than
    /// Neanderthal-count against Denisovan-count.
    pub base_rate_neanderthal: f64,
    pub base_rate_denisovan: f64,
    /// Minimum diagnostic matches of the winning lineage before any attribution is made. A handful
    /// of matches on a sub-megabase segment is noise at these base rates.
    pub min_lineage_matches: usize,
    /// Whether to attempt per-segment lineage attribution at all. **Default `false`.**
    ///
    /// Attribution is built and unit-tested but is NOT validated, and on the ground-truth European
    /// it produces the opposite of the known pattern: 0.00 Mb Neanderthal against 0.48 Mb
    /// Denisovan, where §7 expects essentially all Neanderthal and no Denisovan. The cause is
    /// visible in the base rates — a European carries the derived allele at 4.3 % of
    /// Neanderthal-diagnostic sites and 3.9 % of Denisovan-diagnostic ones, a ratio of 1.10, so
    /// there is almost no discriminating signal to work with at segment scale.
    ///
    /// Shipping a lineage split on that basis would manufacture exactly the Denisovan-in-Europeans
    /// claim §7 forbids, so segments are reported as archaic-but-unattributed until the method is
    /// validated. This is the same discipline that gated ancient ancestry.
    pub attribute_lineage: bool,
    /// Minimum callable fraction for a window to be modelled at all. Windows below this are
    /// **excluded**, not down-weighted: a mostly-uncallable window's variant density is dominated by
    /// mapping error, and including it is what made the first real run call 3.62 % archaic out of
    /// repetitive sequence.
    pub min_callable_fraction: f64,
}

impl Default for ArchaicConfig {
    fn default() -> Self {
        ArchaicConfig {
            window_bp: 1_000,
            rate_background: None,
            // Calibrated against the hmmix 1000G callset (Zenodo, CC BY 4.0) on chr21+chr22 of the
            // ground-truth European: 45 segments / 2.01 Mb against their EUR target of 43 / 2.09 Mb,
            // both inside the p10-p90 spread (35-51 segments, 1.51-2.65 Mb).
            archaic_rate_multiple: 6.0,
            prior_archaic: 0.02,
            // 5, not 1. The transition rate IS the tract-length prior: at 1.0 with the 1 cM/Mb
            // fallback a 1 kb window switches with probability ~0.001, implying ~1 Mb tracts against
            // a real median of 31 kb. That single parameter was why the caller produced a third as
            // many segments, each several times too long.
            switches_per_cm: 5.0,
            // 5 kb, not 50 kb. hmmix's median European tract on these chromosomes is 31 kb and its
            // p10 is 7 kb, so a 50 kb floor discarded more than half of all real segments by
            // construction.
            min_segment_bp: 5_000,
            min_posterior: 0.70,
            min_lineage_ratio: 2.0,
            base_rate_neanderthal: 0.043,
            base_rate_denisovan: 0.039,
            min_lineage_matches: 10,
            attribute_lineage: false,
            min_callable_fraction: 0.5,
        }
    }
}

/// ln(k!) for the small counts a 1 kb window produces.
fn ln_factorial(k: u32) -> f64 {
    (2..=k).map(|i| (i as f64).ln()).sum()
}

/// Poisson log-pmf.
fn ln_poisson(k: u32, lambda: f64) -> f64 {
    let l = lambda.max(1e-9);
    k as f64 * l.ln() - l - ln_factorial(k)
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m.is_infinite() {
        return m;
    }
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// cM between two bp positions, falling back to 1 cM/Mb when the map lacks the contig.
fn span_cm(gmap: &GeneticMap, chr: &str, start_bp: i64, end_bp: i64) -> f64 {
    gmap.interval_cm(chr, start_bp as i32, end_bp as i32)
        .unwrap_or_else(|| (end_bp - start_bp).max(0) as f64 / 1_000_000.0)
}

/// Call archaic tracts from a subject's genome-wide diploid calls.
///
/// `calls` should be the de-novo diploid variant calls for one alignment (Tier B is gated to
/// WGS/VCF input — a chip can not supply the density this needs).
pub fn call_archaic_segments(
    calls: &[SiteGenotype],
    outgroup: &ArchaicOutgroup,
    classify: &ArchaicClassify,
    callable: &ArchaicCallable,
    gmap: &GeneticMap,
    cfg: &ArchaicConfig,
) -> ArchaicSegmentResult {
    // Per-contig: the subject's variant positions (carrying a non-reference allele), plus the
    // callable extent. A dosage of 0 is a called reference site — real information about coverage,
    // but not a variant.
    let mut variants: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    let mut extent: BTreeMap<String, (i64, i64)> = BTreeMap::new();
    let mut alleles: BTreeMap<(String, i64), (char, char, i32)> = BTreeMap::new();
    for g in calls {
        if !(0..=g.ploidy as i32).contains(&g.dosage) {
            continue;
        }
        // Key by the EXACT contig name, never normalized. The Tier B assets are built from the
        // CHM13 VCFs and carry their naming; the subject's calls are on the same build and carry the
        // same naming. Normalizing one side silently breaks the join, and the failure mode is
        // invisible — retain_private finds no contig, returns nothing, and the genome reads as 0%
        // archaic, which looks like a perfectly plausible answer.
        let contig = g.contig.clone();
        let e = extent.entry(contig.clone()).or_insert((g.position, g.position));
        e.0 = e.0.min(g.position);
        e.1 = e.1.max(g.position);
        if g.dosage >= 1 {
            variants.entry(contig.clone()).or_default().push(g.position);
            if let (Some(r), Some(a)) = (g.reference_allele.chars().next(), g.alternate_allele.chars().next()) {
                alleles.insert((contig, g.position), (r, a, g.dosage));
            }
        }
    }

    // Strip everything the African outgroup also carries — the step that turns raw variant density
    // into an introgression signal.
    let mut private: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    let mut total_private = 0usize;
    let mut total_windows = 0f64;
    for (contig, mut positions) in variants {
        positions.sort_unstable();
        positions.dedup();
        let kept = outgroup.retain_private(&contig, &positions);
        // Count only callable windows toward the background rate — otherwise the denominator is
        // inflated by regions the model will not look at anyway.
        if let Some((lo, hi)) = extent.get(&contig) {
            let mut w = *lo;
            while w <= *hi {
                if callable.callable_fraction(&contig, w) >= cfg.min_callable_fraction {
                    total_windows += 1.0;
                }
                w += cfg.window_bp;
            }
        }
        total_private += kept
            .iter()
            .filter(|&&p| callable.callable_fraction(&contig, p) >= cfg.min_callable_fraction)
            .count();
        private.insert(contig, kept);
    }

    // Background rate: the sample's own genome-wide private density unless pinned. Adapts to depth
    // and to how much the outgroup stripped, both of which vary per sample.
    let background = cfg
        .rate_background
        .unwrap_or_else(|| (total_private as f64 / total_windows.max(1.0)).clamp(0.001, 10.0));
    let archaic_rate = background * cfg.archaic_rate_multiple.max(1.1);

    let mut segments = Vec::new();
    for (contig, positions) in &private {
        let Some(&(lo, hi)) = extent.get(contig) else { continue };
        if positions.is_empty() || hi <= lo {
            continue;
        }
        segments.extend(call_contig(
            contig,
            positions,
            lo,
            hi,
            gmap,
            cfg,
            background,
            archaic_rate,
            classify,
            callable,
            &alleles,
        ));
    }

    // Both sides of the ratio must be measured in the SAME units. A segment's span includes windows
    // the mask excluded, so counting span against callable megabases mixes them and inflates the
    // percentage — the first run read 4.80% of "96.4 Mb callable" while the callable track was
    // 44.6 Mb. Archaic extent is therefore summed over callable bases only.
    let seg_callable_mb = |s: &ArchaicSegment| -> f64 {
        let mut bp = 0.0;
        let mut w = s.start;
        while w <= s.end {
            bp += callable.callable_fraction(&s.contig, w) * cfg.window_bp as f64;
            w += cfg.window_bp;
        }
        bp / 1_000_000.0
    };
    let (mut nea, mut den, mut unk) = (0.0, 0.0, 0.0);
    for s in &segments {
        let mb = seg_callable_mb(s);
        match s.source {
            ArchaicSource::Neanderthal => nea += mb,
            ArchaicSource::Denisovan => den += mb,
            ArchaicSource::Unknown => unk += mb,
        }
    }
    let total_mb = nea + den + unk;
    // Denominator: callable bases within the analysed contigs, from the mask itself.
    let callable_mb: f64 = callable
        .contigs
        .iter()
        .filter(|c| private.contains_key(&c.contig))
        .flat_map(|c| c.callable_bp.iter())
        .map(|&b| b as f64)
        .sum::<f64>()
        / 1_000_000.0;
    let summary = ArchaicSummary {
        total_mb,
        pct_callable: if callable_mb > 0.0 {
            total_mb * 100.0 / callable_mb
        } else {
            0.0
        },
        callable_mb,
        neanderthal_mb: nea,
        denisovan_mb: den,
        unknown_mb: unk,
        n_segments: segments.len(),
    };
    ArchaicSegmentResult { segments, summary }
}

/// The per-contig HMM: window the callable span, run Viterbi + forward/backward, stitch archaic
/// runs, then attribute each to a lineage.
#[allow(clippy::too_many_arguments)]
fn call_contig(
    contig: &str,
    private: &[i64],
    lo: i64,
    hi: i64,
    gmap: &GeneticMap,
    cfg: &ArchaicConfig,
    background: f64,
    archaic_rate: f64,
    classify: &ArchaicClassify,
    callable: &ArchaicCallable,
    alleles: &BTreeMap<(String, i64), (char, char, i32)>,
) -> Vec<ArchaicSegment> {
    let n_windows = (((hi - lo) / cfg.window_bp) + 1) as usize;
    if n_windows < 2 {
        return Vec::new();
    }
    let mut counts = vec![0u32; n_windows];
    for &p in private {
        let idx = ((p - lo) / cfg.window_bp) as usize;
        if idx < n_windows {
            counts[idx] += 1;
        }
    }

    // Callable fraction per window. A window below the floor is uninformative: it emits nothing in
    // either state, so it neither supports nor breaks a segment.
    let frac: Vec<f64> = (0..n_windows)
        .map(|i| callable.callable_fraction(contig, lo + i as i64 * cfg.window_bp))
        .collect();
    let usable: Vec<bool> = frac.iter().map(|f| *f >= cfg.min_callable_fraction).collect();
    if !usable.iter().any(|u| *u) {
        return Vec::new();
    }

    let ln = |x: f64| x.max(1e-300).ln();
    let ln_pi = [ln(1.0 - cfg.prior_archaic), ln(cfg.prior_archaic)];
    // Expected counts scale with how much of the window is actually callable, so a half-callable
    // window is not mistaken for a variant-poor one.
    let emit = |i: usize, k: u32| -> [f64; 2] {
        if !usable[i] {
            return [0.0, 0.0];
        }
        let f = frac[i].max(1e-3);
        [ln_poisson(k, background * f), ln_poisson(k, archaic_rate * f)]
    };

    // Transition log-probabilities between adjacent windows, scaled by genetic distance: tracts
    // break at recombination, so a wide gap in cM means a switch is likelier.
    let trans = |i: usize| -> [[f64; 2]; 2] {
        let a = lo + i as i64 * cfg.window_bp;
        let b = a + cfg.window_bp;
        let cm = span_cm(gmap, contig, a, b);
        let sw = (1.0 - (-cfg.switches_per_cm * cm.max(0.0)).exp()).clamp(1e-6, 0.5);
        [[ln(1.0 - sw), ln(sw)], [ln(sw), ln(1.0 - sw)]]
    };

    // Forward.
    let mut fwd = vec![[f64::NEG_INFINITY; 2]; n_windows];
    let e0 = emit(0, counts[0]);
    fwd[0] = [ln_pi[0] + e0[0], ln_pi[1] + e0[1]];
    for i in 1..n_windows {
        let t = trans(i - 1);
        let e = emit(i, counts[i]);
        for s in 0..2 {
            fwd[i][s] = ln_sum_exp(fwd[i - 1][0] + t[0][s], fwd[i - 1][1] + t[1][s]) + e[s];
        }
    }
    // Backward.
    let mut bwd = vec![[0.0f64; 2]; n_windows];
    for i in (0..n_windows - 1).rev() {
        let t = trans(i);
        let e = emit(i + 1, counts[i + 1]);
        for s in 0..2 {
            bwd[i][s] = ln_sum_exp(t[s][0] + e[0] + bwd[i + 1][0], t[s][1] + e[1] + bwd[i + 1][1]);
        }
    }
    // Posterior of the archaic state per window.
    let post: Vec<f64> = (0..n_windows)
        .map(|i| {
            let a = fwd[i][0] + bwd[i][0];
            let b = fwd[i][1] + bwd[i][1];
            let m = a.max(b);
            let (ea, eb) = ((a - m).exp(), (b - m).exp());
            eb / (ea + eb)
        })
        .collect();

    // Stitch runs of archaic-posterior windows, then apply the length and confidence floors.
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < n_windows {
        if post[i] < 0.5 {
            i += 1;
            continue;
        }
        let start_w = i;
        while i < n_windows && post[i] >= 0.5 {
            i += 1;
        }
        let end_w = i - 1;
        let start = lo + start_w as i64 * cfg.window_bp;
        let end = (lo + (end_w as i64 + 1) * cfg.window_bp - 1).min(hi);
        if end - start < cfg.min_segment_bp {
            continue;
        }
        let mean_post = post[start_w..=end_w].iter().sum::<f64>() / (end_w - start_w + 1) as f64;
        if mean_post < cfg.min_posterior {
            continue;
        }
        // At least half the run's windows must be callable, so a tract can not be carried by a
        // stretch of uninformative windows riding on the transition prior.
        let callable_windows = usable[start_w..=end_w].iter().filter(|u| **u).count();
        if callable_windows * 2 < end_w - start_w + 1 {
            continue;
        }
        let n_private = private.iter().filter(|&&p| p >= start && p <= end).count();
        let (source, nea, den) = attribute(contig, start, end, classify, alleles, cfg);
        out.push(ArchaicSegment {
            contig: contig.to_string(),
            start,
            end,
            posterior: mean_post,
            n_private,
            source,
            neanderthal_matches: nea,
            denisovan_matches: den,
        });
    }
    out
}

/// Attribute a segment by counting the subject's derived-allele matches at diagnostic sites.
///
/// Requires a clear margin (`min_lineage_ratio`); a segment whose evidence is balanced, or which
/// has no diagnostic sites at all, is `Unknown`. Guessing here would manufacture exactly the
/// Denisovan-in-Europeans claim §7 forbids.
fn attribute(
    contig: &str,
    start: i64,
    end: i64,
    classify: &ArchaicClassify,
    alleles: &BTreeMap<(String, i64), (char, char, i32)>,
    cfg: &ArchaicConfig,
) -> (ArchaicSource, usize, usize) {
    // Count both MATCHES and the diagnostic sites available, per lineage: the enrichment is
    // matches / (sites x base rate), and without the site counts the comparison silently favours
    // whichever lineage happens to have more sites in this segment.
    let (mut nea, mut den) = (0usize, 0usize);
    let (mut nea_sites, mut den_sites) = (0usize, 0usize);
    for (pos, derived, class) in classify.in_range(contig, start, end) {
        match class {
            DiagnosticClass::Neanderthal => nea_sites += 1,
            DiagnosticClass::Denisovan => den_sites += 1,
            DiagnosticClass::SharedArchaic => {}
        }
        let Some(&(r, a, dosage)) = alleles.get(&(contig.to_string(), pos)) else {
            continue;
        };
        let d = derived.to_ascii_uppercase();
        // Does the subject actually carry the archaic-derived base here?
        let carries = (a.to_ascii_uppercase() == d && dosage >= 1) || (r.to_ascii_uppercase() == d && dosage <= 1);
        if !carries {
            continue;
        }
        match class {
            DiagnosticClass::Neanderthal => nea += 1,
            DiagnosticClass::Denisovan => den += 1,
            DiagnosticClass::SharedArchaic => {}
        }
    }
    // Enrichment over the expected carrier rate, not raw counts.
    let exp_nea = nea_sites as f64 * cfg.base_rate_neanderthal;
    let exp_den = den_sites as f64 * cfg.base_rate_denisovan;
    let enr_nea = if exp_nea > 0.0 { nea as f64 / exp_nea } else { 0.0 };
    let enr_den = if exp_den > 0.0 { den as f64 / exp_den } else { 0.0 };

    if !cfg.attribute_lineage {
        // Machinery preserved and unit-tested; the label is withheld until it validates.
        return (ArchaicSource::Unknown, nea, den);
    }
    let source = if nea >= cfg.min_lineage_matches && enr_nea >= enr_den * cfg.min_lineage_ratio {
        ArchaicSource::Neanderthal
    } else if den >= cfg.min_lineage_matches && enr_den >= enr_nea * cfg.min_lineage_ratio {
        ArchaicSource::Denisovan
    } else {
        // Not enough evidence to separate the lineages. For a European this is the expected
        // outcome and the honest one — §7 forbids manufacturing a Denisovan finding.
        ArchaicSource::Unknown
    };
    (source, nea, den)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archaic::{ArchaicCallable, CallableContig, ClassifyContig, PositionStream};

    fn gt(contig: &str, position: i64, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: contig.into(),
            position,
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

    /// Background variants every 5 kb across 2 Mb, plus a dense block (every 200 bp) in
    /// [1.0 Mb, 1.3 Mb] — an introgressed tract. None are in the outgroup.
    /// Fully-callable track over the synthetic 2 Mb contig.
    fn callable_all() -> ArchaicCallable {
        ArchaicCallable {
            build: "chm13v2.0".into(),
            window_bp: 1_000,
            contigs: vec![CallableContig {
                contig: "chr21".into(),
                start: 0,
                callable_bp: vec![1_000u16; 2_100],
            }],
        }
    }

    fn sample() -> (Vec<SiteGenotype>, ArchaicOutgroup) {
        let mut calls = Vec::new();
        let mut p = 1i64;
        while p < 2_000_000 {
            calls.push(gt("chr21", p, 1));
            p += 5_000;
        }
        let mut q = 1_000_000i64;
        while q < 1_300_000 {
            calls.push(gt("chr21", q, 1));
            q += 200;
        }
        let og = ArchaicOutgroup {
            build: "chm13v2.0".into(),
            min_allele_count: 1,
            contigs: vec![PositionStream::encode("chr21", &[])],
        };
        (calls, og)
    }

    #[test]
    fn finds_a_dense_private_block_and_ignores_background() {
        let (calls, og) = sample();
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &callable_all(),
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert_eq!(r.segments.len(), 1, "exactly the dense block should call");
        let s = &r.segments[0];
        assert!(s.start >= 950_000 && s.start <= 1_050_000, "start {} off", s.start);
        assert!(s.end >= 1_250_000 && s.end <= 1_350_000, "end {} off", s.end);
        assert!(s.posterior > 0.9);
        // No diagnostic sites supplied -> honestly Unknown, not guessed.
        assert_eq!(s.source, ArchaicSource::Unknown);
        assert!(r.summary.total_mb > 0.2 && r.summary.total_mb < 0.4);
    }

    #[test]
    fn stripping_the_outgroup_is_what_makes_the_signal() {
        // Same data, but now the outgroup carries every variant in the dense block: the density
        // excess vanishes and nothing should be called. This is the step that separates
        // introgression from ordinary variation.
        let (calls, _) = sample();
        let mut shared: Vec<i64> = (1_000_000..1_300_000).step_by(200).collect();
        shared.sort_unstable();
        let og = ArchaicOutgroup {
            build: "chm13v2.0".into(),
            min_allele_count: 1,
            contigs: vec![PositionStream::encode("chr21", &shared)],
        };
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &callable_all(),
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert!(r.segments.is_empty(), "outgroup-shared density must not call archaic");
    }

    #[test]
    fn an_uncallable_dense_region_is_not_called() {
        // The same dense block, but the mask says that region is not callable — the exact situation
        // that made the first real run report 3.62% archaic out of repetitive sequence.
        let (calls, og) = sample();
        let mut track = callable_all();
        for w in 950..1_350 {
            track.contigs[0].callable_bp[w] = 0;
        }
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &track,
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert!(
            r.segments.is_empty(),
            "density in an uncallable region must not be called archaic, got {:?}",
            r.segments
        );
    }

    #[test]
    fn attribution_compares_enrichment_not_raw_counts() {
        // 40 Neanderthal-diagnostic sites and 40 Denisovan ones in the segment. The subject matches
        // 12 of each. Raw counts say "tie"; against base rates of 4.3% vs 3.9% the expected counts
        // are 1.7 and 1.6, so both are enriched almost identically -> genuinely Unknown, which is
        // the honest answer rather than a coin flip.
        let mut positions: Vec<i64> = Vec::new();
        let mut derived: Vec<u8> = Vec::new();
        let mut classes: Vec<u8> = Vec::new();
        for i in 0..80i64 {
            positions.push(1_000_000 + i * 10);
            derived.push(b'G');
            classes.push(if i < 40 { 0 } else { 1 });
        }
        let cls = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &positions),
                derived,
                classes,
            }],
        };
        let mut alleles = BTreeMap::new();
        for i in 0..12i64 {
            alleles.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        for i in 40..52i64 {
            alleles.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        // Attribution is off by default (see `attribute_lineage`); this test exercises the logic.
        let cfg = ArchaicConfig {
            attribute_lineage: true,
            ..ArchaicConfig::default()
        };
        let (src, nea, den) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles, &cfg);
        assert_eq!((nea, den), (12, 12));
        assert_eq!(src, ArchaicSource::Unknown, "equal enrichment must not pick a lineage");

        // Now give Neanderthal a genuine excess and strip the Denisovan matches: it should call.
        let mut alleles2 = BTreeMap::new();
        for i in 0..30i64 {
            alleles2.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        let (src, nea, den) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles2, &cfg);
        assert_eq!((nea, den), (30, 0));
        assert_eq!(src, ArchaicSource::Neanderthal);

        // Too few matches to mean anything at these base rates -> Unknown, not a guess.
        let mut alleles3 = BTreeMap::new();
        for i in 0..3i64 {
            alleles3.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        let (src, _, _) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles3, &cfg);
        assert_eq!(src, ArchaicSource::Unknown);
    }
}

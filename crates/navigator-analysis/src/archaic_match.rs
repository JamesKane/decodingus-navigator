//! Tier B, second attempt — archaic tracts by **matching the archaic genomes**, not by counting
//! mutations.
//!
//! # Why this replaces the density caller
//!
//! [`crate::archaic_segments`] follows Skov 2018 (hmmix): strip variants Africans also carry, then
//! look for regions dense in what remains. That method exists for people who do **not** have archaic
//! reference genomes and must infer them indirectly. We have all four, and already ship
//! [`ArchaicClassify`] — 2,031,406 sites where the archaics carry a derived allele.
//!
//! Measured on a real European against hmmix's own calls for the same person, the difference is not
//! subtle. Both observables carry the same ~3x contrast, but they differ 30-fold in how much
//! evidence one tract holds, and that is what decides whether a tract can be called at all:
//!
//! | observable | evidence per 36 kb tract | sensitivity at 5 % false positives |
//! |---|---|---|
//! | private-variant density | ~1 variant | 14.3 % |
//! | archaic-allele matching (this) | ~30 sites | 95.1 % |
//!
//! Density does not reach 80 % sensitivity at **500 kb**; matching reaches 95 % at the real median
//! tract of 36 kb. See `documents/design/ArchaicAncestry_Design.md` § *Why it failed*.
//!
//! # The model
//!
//! An introgressed tract is a haplotype inherited intact from an archaic ancestor, so it carries the
//! archaic allele at a large share of the diagnostic sites it spans; elsewhere the subject carries
//! them only at the background rate. That is a two-state HMM whose observation is one **bit per
//! diagnostic site** — carried or not — with Bernoulli emissions, indexed **by site rather than by
//! base pair**.
//!
//! Indexing by site is what makes this robust where the density model was not. Diagnostic sites
//! become the denominator, so their uneven density cancels out: the mutation-rate map the density
//! model needed (and which no available proxy supplied — the best explained 38 % of a 14.6x
//! overdispersion) is simply not required here.
//!
//! Transitions stay recombination-scaled between consecutive sites, as in [`crate::roh`] and the
//! chromosome painter.
//!
//! # Validation
//!
//! Scored against hmmix's own calls for the same individuals, 60 Europeans on chr21+22, split 30
//! **train** / 30 **test** on a fixed seed. Thresholds were fitted on train only; every figure below
//! is the held-out half. The split exists because the previous caller was tuned until a cohort
//! statistic matched and the statistic was then reported as evidence.
//!
//! | | density caller | this, uncalibrated | this, calibrated |
//! |---|---|---|---|
//! | base-level F1 | — | 27.9 % | **34.5 %** |
//! | precision | 1.5 % | 20.2 % | **34.9 %** |
//! | extent ratio ours/theirs | 1.45 | 2.23 | **0.98** |
//! | per-individual extent `r` | −0.018 (p = 0.94) | +0.520 | **+0.710 (p < 0.0001)** |
//!
//! The extent ratio of 0.98 is the one to notice: the caller is no longer systematically
//! over-calling, which the emission-ratio sweep is what fixed.
//!
//! On locations, all 20 individuals of an earlier cohort scored above their own random-placement
//! null (mean 45.3 % sensitivity against a 7.1 % null); the density caller scored 2.1 % against a
//! 5.0 % null, i.e. below chance.
//!
//! **Still not enough to re-enable**, and the limits are specific rather than general unease:
//! precision is 34.9 %, so two thirds of called sequence is not in the reference callset; the cohort
//! is **European only** and **chr21+22 only**; and the reference callset is itself weakly supported
//! (hmmix's own tracts are enriched only 1.84x for their own archaic SNPs), so agreement with it
//! caps out well below 100 % even for a correct caller. `ARCHAIC_SEGMENTS_ENABLED` stays `false`
//! until this reproduces outside Europe — East Asians are the sharp test, since the truth predicts
//! ~1.18x more archaic sequence there and a caller merely tracking European structure would miss it.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::archaic::{ArchaicCallable, ArchaicClassify, DiagnosticClass};
use crate::archaic_segments::{ArchaicSegment, ArchaicSegmentResult, ArchaicSource, ArchaicSummary};
use crate::caller::SiteGenotype;
use crate::ibd::GeneticMap;

/// One diagnostic site, reduced to what the HMM consumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteObs {
    pub position: i64,
    /// Whether the subject carries the archaic-derived allele here.
    pub carries: bool,
    pub class: DiagnosticClass,
}

/// Tuning knobs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MatchConfig {
    /// Rate at which a **non-introgressed** genome carries the archaic allele at a diagnostic site.
    /// `None` estimates it from the subject's own genome-wide rate, which absorbs coverage, call
    /// behaviour and ancestry.
    ///
    /// Estimated rather than fixed because it is the denominator of the whole inference — but
    /// estimated *directly*, not by EM. Unconstrained Baum-Welch on the previous caller diverged to
    /// a degenerate fit (a 22x emission ratio and 9 kb tracts, calling 7x the truth), so parameters
    /// here are measured, not fitted.
    pub p_background: Option<f64>,
    /// Rate inside an introgressed tract. `None` derives it as `p_background * archaic_ratio`.
    pub p_archaic: Option<f64>,
    /// Multiple of the background rate expected inside a tract when `p_archaic` is `None`.
    /// Measured at 3.04x (39.5 % inside real tracts against 13.0 % elsewhere).
    pub archaic_ratio: f64,
    /// Expected state switches per centimorgan.
    pub switches_per_cm: f64,
    /// Discard tracts whose mean posterior is below this.
    pub min_posterior: f64,
    /// Minimum diagnostic sites in a tract. A tract resting on one or two sites is exactly the
    /// failure mode of the density caller, restated in a new observable.
    pub min_sites: usize,
    /// Discard tracts shorter than this.
    pub min_segment_bp: i64,
    /// Minimum callable fraction for a site's window to be used at all.
    pub min_callable_fraction: f64,
    /// Whether to attempt per-segment Neanderthal/Denisovan attribution. Default `false`, unchanged
    /// from the density caller: the lineage signal has not been shown to work, and this module does
    /// not by itself change that.
    pub attribute_lineage: bool,
}

impl Default for MatchConfig {
    fn default() -> Self {
        MatchConfig {
            p_background: None,
            p_archaic: None,
            // FITTED (not measured): the observed enrichment inside real tracts is 3.04x, but the
            // model separates best at 4.5x. That is not a contradiction — 3.04x is the *average*
            // over an external tract set that is itself only weakly supported, while the emission
            // ratio is what makes the HMM selective enough to place boundaries. Fitted on 30
            // Europeans, reported on 30 held-out ones; it is the parameter that removed the
            // over-calling (extent ratio 2.23 -> 0.98).
            archaic_ratio: 4.5,
            switches_per_cm: 1.0,
            // All three CALIBRATED on train, reported on held-out test (see the module docs).
            // Objective was base-level F1: sensitivity alone is bought by calling more sequence, and
            // the uncalibrated caller over-called 2.2x while still scoring 45 %.
            min_posterior: 0.98,
            min_sites: 16,
            // 5 kb, though the grid's argmax preferred 10 kb. Within the plateau the two differ by
            // 0.1 F1 points, the 5 kb floor is slightly BETTER on per-individual extent correlation
            // (+0.710 vs +0.706), and it discards half as many real tracts (8 % of the truth under
            // 5 kb against 16 % under 10 kb). An earlier sweep wanted 40 kb, which would have
            // discarded 61 %; the design records the same trap once before at 50 kb. Structural
            // exclusion of real tracts is not worth a tenth of a point.
            min_segment_bp: 5_000,
            min_callable_fraction: 0.5,
            attribute_lineage: false,
        }
    }
}

/// Reduce one contig's diagnostic sites to observations.
///
/// `ref_base` supplies the reference base at a position; sites where the archaic-derived allele
/// **is** the reference base are dropped. At such a site every reference-matching genome trivially
/// "carries" the derived allele, so it separates nothing and would dilute the contrast — and,
/// because the caller emits only variant records, a no-call there means the subject *does* carry it,
/// the opposite of what a no-call means everywhere else.
///
/// A site with no variant record is hom-reference, hence **not** carrying. Restricting instead to
/// sites where the subject happens to have a call is the trap that made an early version of this
/// analysis report an 80 % carrying rate against a known 4.3 % background: it samples only sites
/// where a variant already exists.
pub fn observations_for_contig(
    contig: &str,
    classify: &ArchaicClassify,
    calls_by_pos: &BTreeMap<i64, &SiteGenotype>,
    ref_base: impl Fn(i64) -> Option<u8>,
    callable: &ArchaicCallable,
    min_callable_fraction: f64,
) -> Vec<SiteObs> {
    let Some(c) = classify.contigs.iter().find(|c| c.positions.contig == contig) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for (i, pos) in c.positions.iter().enumerate() {
        let Some(&derived) = c.derived.get(i) else { continue };
        // Uninformative: the reference already carries the archaic allele.
        if ref_base(pos) == Some(derived) {
            continue;
        }
        if callable.callable_fraction(contig, pos) < min_callable_fraction {
            continue;
        }
        let carries = calls_by_pos.get(&pos).is_some_and(|g| {
            g.dosage > 0 && g.alternate_allele.as_bytes().first() == Some(&derived)
        });
        let class = match c.classes.get(i).copied().unwrap_or(2) {
            0 => DiagnosticClass::Neanderthal,
            1 => DiagnosticClass::Denisovan,
            _ => DiagnosticClass::SharedArchaic,
        };
        out.push(SiteObs {
            position: pos,
            carries,
            class,
        });
    }
    out
}

fn ln(x: f64) -> f64 {
    x.max(1e-300).ln()
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Posterior probability of the archaic state at each observation.
///
/// Log-space forward/backward with recombination-scaled transitions, as in [`crate::roh`]. Exposed
/// so the decoding can be tested against hand-computed posteriors without constructing assets.
pub fn posteriors(obs: &[SiteObs], contig: &str, gmap: &GeneticMap, p_bg: f64, p_arch: f64, switches_per_cm: f64) -> Vec<f64> {
    let n = obs.len();
    if n == 0 {
        return Vec::new();
    }
    let emit = |i: usize| -> [f64; 2] {
        if obs[i].carries {
            [ln(p_bg), ln(p_arch)]
        } else {
            [ln(1.0 - p_bg), ln(1.0 - p_arch)]
        }
    };
    // Switch probability between consecutive sites, from the genetic distance between them.
    let sw = |i: usize| -> f64 {
        let cm = gmap
            .interval_cm(contig, obs[i].position as i32, obs[i + 1].position as i32)
            .unwrap_or_else(|| (obs[i + 1].position - obs[i].position).max(0) as f64 / 1_000_000.0);
        (1.0 - (-switches_per_cm * cm.max(0.0)).exp()).clamp(1e-9, 0.5)
    };

    // Prior: the stationary share of the archaic state, from the rates themselves rather than a
    // tuned constant — with p_arch > p_bg the algebra puts it at a few percent, matching reality.
    let prior_arch = ((p_bg - (1.0 - p_arch) * 0.0) / p_arch).clamp(0.001, 0.5) * 0.1;
    let mut fwd = vec![[f64::NEG_INFINITY; 2]; n];
    let e0 = emit(0);
    fwd[0] = [ln(1.0 - prior_arch) + e0[0], ln(prior_arch) + e0[1]];
    for i in 1..n {
        let s = sw(i - 1);
        let (stay, go) = (ln(1.0 - s), ln(s));
        let e = emit(i);
        for st in 0..2 {
            let from0 = fwd[i - 1][0] + if st == 0 { stay } else { go };
            let from1 = fwd[i - 1][1] + if st == 1 { stay } else { go };
            fwd[i][st] = ln_sum_exp(from0, from1) + e[st];
        }
    }
    let mut bwd = vec![[0.0f64; 2]; n];
    for i in (0..n - 1).rev() {
        let s = sw(i);
        let (stay, go) = (ln(1.0 - s), ln(s));
        let e = emit(i + 1);
        for st in 0..2 {
            let to0 = bwd[i + 1][0] + e[0] + if st == 0 { stay } else { go };
            let to1 = bwd[i + 1][1] + e[1] + if st == 1 { stay } else { go };
            bwd[i][st] = ln_sum_exp(to0, to1);
        }
    }
    let total = ln_sum_exp(fwd[n - 1][0], fwd[n - 1][1]);
    (0..n)
        .map(|i| (fwd[i][1] + bwd[i][1] - total).exp().clamp(0.0, 1.0))
        .collect()
}

/// Call archaic tracts for one subject by matching the archaic genomes.
///
/// `observations` is per contig, already reduced by [`observations_for_contig`], so this function
/// does no I/O and no asset decoding — it is the model, and is unit-testable as such.
pub fn call_from_observations(
    observations: &BTreeMap<String, Vec<SiteObs>>,
    gmap: &GeneticMap,
    callable: &ArchaicCallable,
    cfg: &MatchConfig,
) -> ArchaicSegmentResult {
    let (carried, total): (usize, usize) = observations
        .values()
        .flatten()
        .fold((0, 0), |(c, t), o| (c + usize::from(o.carries), t + 1));
    if total == 0 {
        return ArchaicSegmentResult {
            segments: Vec::new(),
            summary: ArchaicSummary {
                total_mb: 0.0,
                pct_callable: 0.0,
                callable_mb: 0.0,
                neanderthal_mb: 0.0,
                denisovan_mb: 0.0,
                unknown_mb: 0.0,
                n_segments: 0,
            },
        };
    }
    // The genome-wide rate is dominated by non-archaic sequence (archaic tracts are a few percent
    // of it), so it estimates the background directly.
    let p_bg = cfg
        .p_background
        .unwrap_or((carried as f64 / total as f64).clamp(0.001, 0.5));
    let p_arch = cfg
        .p_archaic
        .unwrap_or((p_bg * cfg.archaic_ratio).clamp(p_bg * 1.1, 0.95));

    let mut segments = Vec::new();
    for (contig, obs) in observations {
        if obs.len() < cfg.min_sites {
            continue;
        }
        let post = posteriors(obs, contig, gmap, p_bg, p_arch, cfg.switches_per_cm);
        let mut i = 0usize;
        while i < post.len() {
            if post[i] < cfg.min_posterior {
                i += 1;
                continue;
            }
            let start = i;
            while i < post.len() && post[i] >= cfg.min_posterior {
                i += 1;
            }
            let end = i - 1;
            let n_sites = end - start + 1;
            let span = obs[end].position - obs[start].position;
            if n_sites < cfg.min_sites || span < cfg.min_segment_bp {
                continue;
            }
            let mean_post = post[start..=end].iter().sum::<f64>() / n_sites as f64;
            let (mut nea, mut den) = (0usize, 0usize);
            for o in &obs[start..=end] {
                if !o.carries {
                    continue;
                }
                match o.class {
                    DiagnosticClass::Neanderthal => nea += 1,
                    DiagnosticClass::Denisovan => den += 1,
                    DiagnosticClass::SharedArchaic => {}
                }
            }
            segments.push(ArchaicSegment {
                contig: contig.clone(),
                start: obs[start].position,
                end: obs[end].position,
                posterior: mean_post,
                n_private: obs[start..=end].iter().filter(|o| o.carries).count(),
                // Attribution stays off by default; the lineage signal is a separate question this
                // module does not answer (see `MatchConfig::attribute_lineage`).
                source: ArchaicSource::Unknown,
                neanderthal_matches: nea,
                denisovan_matches: den,
            });
        }
    }

    let callable_mb: f64 = callable
        .contigs
        .iter()
        .filter(|c| observations.contains_key(&c.contig))
        .flat_map(|c| c.callable_bp.iter())
        .map(|&b| b as f64)
        .sum::<f64>()
        / 1_000_000.0;
    let total_mb: f64 = segments.iter().map(|s| s.length_mb()).sum();
    let summary = ArchaicSummary {
        total_mb,
        pct_callable: if callable_mb > 0.0 { total_mb * 100.0 / callable_mb } else { 0.0 },
        callable_mb,
        neanderthal_mb: 0.0,
        denisovan_mb: 0.0,
        unknown_mb: total_mb,
        n_segments: segments.len(),
    };
    ArchaicSegmentResult { segments, summary }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archaic::{CallableContig, ClassifyContig, PositionStream};

    fn gmap(contig: &str, len: i32) -> GeneticMap {
        GeneticMap::uniform(1.0, &[(contig, len)])
    }

    fn callable(contig: &str, windows: usize) -> ArchaicCallable {
        ArchaicCallable {
            build: "chm13v2.0".into(),
            window_bp: 1_000,
            contigs: vec![CallableContig {
                contig: contig.into(),
                start: 0,
                callable_bp: vec![1_000u16; windows],
            }],
        }
    }

    fn obs(positions: &[(i64, bool)]) -> Vec<SiteObs> {
        positions
            .iter()
            .map(|&(position, carries)| SiteObs {
                position,
                carries,
                class: DiagnosticClass::SharedArchaic,
            })
            .collect()
    }

    /// A run of carried sites against a background of non-carried ones is what a real tract looks
    /// like, and is the thing this model exists to find.
    #[test]
    fn finds_a_run_of_carried_sites() {
        let mut sites: Vec<(i64, bool)> = (0..60).map(|i| (10_000 + i * 500, false)).collect();
        for s in sites.iter_mut().skip(20).take(20) {
            s.1 = true; // a 10 kb tract, 20 diagnostic sites, all carried
        }
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        // Thresholds pinned rather than inherited: this test is about whether the model finds a
        // run at all, and should not move when the calibrated defaults do. (It broke once when
        // `min_posterior` rose to 0.98 and trimmed the run's edges — correct behaviour, wrong
        // thing for this test to be sensitive to.)
        let cfg = MatchConfig {
            p_background: Some(0.13),
            p_archaic: Some(0.40),
            min_posterior: 0.80,
            min_sites: 5,
            min_segment_bp: 1_000,
            ..Default::default()
        };
        let r = call_from_observations(&m, &gmap("chr21", 60_000), &callable("chr21", 60), &cfg);
        assert_eq!(r.segments.len(), 1, "one tract expected, got {:?}", r.segments);
        let seg = &r.segments[0];
        assert!(seg.start >= 19_000 && seg.start <= 21_000, "start {} off", seg.start);
        assert!(seg.end >= 29_000 && seg.end <= 31_000, "end {} off", seg.end);
    }

    /// The failure that gated the density caller was calling tracts out of background noise. With
    /// no carried sites at all there is nothing to call, and the model must say so.
    #[test]
    fn calls_nothing_on_a_background_only_contig() {
        let sites: Vec<(i64, bool)> = (0..200).map(|i| (10_000 + i * 500, false)).collect();
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        let r = call_from_observations(
            &m,
            &gmap("chr21", 200_000),
            &callable("chr21", 200),
            &MatchConfig {
                p_background: Some(0.13),
                p_archaic: Some(0.40),
                ..Default::default()
            },
        );
        assert!(r.segments.is_empty(), "background should call nothing, got {:?}", r.segments);
    }

    /// Scattered carried sites at the background rate must not accumulate into a tract — the
    /// density caller's defining failure, restated in this observable.
    #[test]
    fn scattered_background_carriers_do_not_form_a_tract() {
        // 13 % carried, evenly spread: exactly the background rate, no run.
        let sites: Vec<(i64, bool)> = (0..300).map(|i| (10_000 + i * 500, i % 8 == 0)).collect();
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        let r = call_from_observations(
            &m,
            &gmap("chr21", 300_000),
            &callable("chr21", 300),
            &MatchConfig {
                p_background: Some(0.13),
                p_archaic: Some(0.40),
                ..Default::default()
            },
        );
        assert!(r.segments.is_empty(), "background-rate carriers formed {:?}", r.segments);
    }

    /// A site whose derived allele IS the reference base separates nothing, and a no-call there
    /// means the opposite of what it means elsewhere. Such sites must be dropped, not counted.
    #[test]
    fn observations_drop_sites_where_reference_is_derived() {
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &[1_000, 2_000, 3_000]),
                derived: vec![b'A', b'C', b'G'],
                classes: vec![0, 1, 2],
            }],
        };
        let calls: BTreeMap<i64, &SiteGenotype> = BTreeMap::new();
        // The reference carries the derived allele at 2_000 only.
        let out = observations_for_contig(
            "chr21",
            &classify,
            &calls,
            |p| if p == 2_000 { Some(b'C') } else { Some(b'T') },
            &callable("chr21", 10),
            0.5,
        );
        assert_eq!(out.len(), 2, "the reference-derived site must be dropped");
        assert!(out.iter().all(|o| o.position != 2_000));
        assert!(out.iter().all(|o| !o.carries), "no calls means nothing carried");
    }

    /// A no-call is hom-reference, i.e. NOT carrying. Conditioning on "has a call" instead is what
    /// made an early version of this analysis report ~80 % carrying against a 4.3 % background.
    #[test]
    fn a_missing_call_is_not_a_carrier() {
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &[1_000, 2_000]),
                derived: vec![b'A', b'A'],
                classes: vec![0, 0],
            }],
        };
        let carried = SiteGenotype {
            name: String::new(),
            contig: "chr21".into(),
            position: 1_000,
            reference_allele: "T".into(),
            alternate_allele: "A".into(),
            ploidy: 2,
            dosage: 1,
            gq: 60,
            depth: 30,
            ref_depth: 15,
            alt_depth: 15,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        };
        let mut calls: BTreeMap<i64, &SiteGenotype> = BTreeMap::new();
        calls.insert(1_000, &carried);
        let out = observations_for_contig("chr21", &classify, &calls, |_| Some(b'T'), &callable("chr21", 10), 0.5);
        assert_eq!(out.len(), 2);
        assert!(out[0].carries, "a called derived allele carries");
        assert!(!out[1].carries, "an absent call is hom-reference, not a carrier");
    }
}

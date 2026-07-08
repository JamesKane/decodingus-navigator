//! Haploid **local reassembly** resolver (private-Y Option B, phase 1) — pure Rust, no external
//! tools, Windows/MSVC-clean.
//!
//! Design: `docs/design/haploid-reassembly-caller.md`. This module owns **Stages B–E** over a single
//! active window: read selection (mapping-quality gate + fragment dedup), candidate haplotypes
//! (per-SNV in v1: reference vs reference-with-one-substitution), read↔haplotype likelihood via a
//! **base-quality-aware PairHMM** (`bio::stats::pairhmm`), and haploid genotyping by the aggregate
//! log-odds. `caller.rs` owns Stage A (active-region detection — it already tallies the per-position
//! counts) and Stage F (turning [`ReassemblyCall`]s into `VariantCall`s); this module is deliberately
//! **I/O-free** so it is unit-testable on synthetic windows.
//!
//! Why it exists: the pileup caller (`caller.rs`) rejects a position whose pileup is ~50/50 as a
//! suspected paralog artifact (`is_paralogous`). At Y segmental-duplication / ampliconic loci that
//! throws away *true* derived SNVs, because reads from a paralogous region mismap and carry the
//! reference base onto the site. GATK resolves these by local reassembly + a base-quality PairHMM;
//! this is the haploid-only equivalent. Proven on WGS229 (POC `examples/reassembly_probe.rs`): the
//! base-quality PairHMM recovers the misaligned-ref sites the crude match/mismatch pileup ties.
//!
//! v1 is **per-candidate-SNV** (one alternate haplotype per candidate position). Linked variants and
//! short indels via POA multi-haplotype assembly are the v2 extension (see the design doc); POA still
//! serves here as an optional cross-check for the caller.

use std::collections::HashMap;

use bio::stats::pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission};
use bio::stats::{LogProb, Prob};

/// Natural-log → Phred scale factor (`10 / ln 10`); `LogProb` is base-*e*.
const PHRED_PER_NAT: f64 = 4.342_944_819_032_518;

/// Tuning for the reassembly resolver. Defaults are the POC-validated starting points; the design
/// doc's §Open-questions flags τ / window size for calibration on the full truth set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReassemblyParams {
    /// Reads below this mapping quality are excluded (GATK default) — this is what drops the
    /// ambiguously-placed paralog reads that masquerade as high-base-quality reference support.
    pub min_mapping_quality: u8,
    /// Minimum aggregate log-odds (nats) for a haploid DERIVED call; symmetric for ANCESTRAL.
    pub min_log_odds: f64,
    /// A DERIVED call needs at least this many alt-supporting fragments (post-dedup).
    pub min_alt_fragments: u32,
}

impl Default for ReassemblyParams {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            min_log_odds: 2.0,
            min_alt_fragments: 2,
        }
    }
}

/// What a read observes at one candidate site (the base it carries there and that base's quality).
/// `None` in [`WindowRead::site_obs`] means the read does not span that candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteObs {
    pub base: u8,
    pub qual: u8,
}

/// A read projected onto the active window's reference frame. Construction (the CIGAR walk that
/// yields the window-frame sequence, per-base qualities, and per-candidate [`SiteObs`]) is the
/// caller's job; this module consumes the projection so it stays I/O-free and testable.
#[derive(Debug, Clone)]
pub struct WindowRead {
    /// Fragment identity (query name) — same name for a read and its mate, used for dedup.
    pub name: Vec<u8>,
    /// Window-frame bases (uppercase), for the whole-read PairHMM realignment.
    pub seq: Vec<u8>,
    /// Per-base Phred qualities, parallel to `seq`.
    pub quals: Vec<u8>,
    /// Mapping quality of the source record.
    pub mapq: u8,
    /// Observation at each candidate site; parallel to the `candidates` slice given to
    /// [`genotype_window`]. `None` = read does not cover that site.
    pub site_obs: Vec<Option<SiteObs>>,
}

/// A candidate variant position within the window (1-based reference coordinate).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Candidate {
    pub position: i64,
    pub ref_base: u8,
    pub alt_base: u8,
}

/// The haploid genotype the resolver assigns to a candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    /// The alternate haplotype explains the reads (a private variant).
    Derived,
    /// The reference haplotype explains the reads (drop it).
    Ancestral,
    /// Neither wins by `min_log_odds` — genuinely undecided (do not call).
    Ambiguous,
}

/// A genotyped candidate. The caller keeps [`Zygosity::Derived`] calls and turns them into
/// `VariantCall`s (Stage F); the others are returned so tests and diagnostics can see the decision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReassemblyCall {
    pub position: i64,
    pub ref_base: u8,
    pub alt_base: u8,
    /// Fragments spanning the site after the MAPQ gate + mate dedup.
    pub depth: u32,
    /// Spanning fragments whose site base is the alternate allele.
    pub alt_depth: u32,
    pub allele_fraction: f64,
    /// Aggregate `Σ ln P(read|alt) − ln P(read|ref)` (nats); >0 favours the alt haplotype.
    pub log_odds: f64,
    /// Phred-scaled confidence of the winning genotype (GQ-like).
    pub quality: f64,
    pub genotype: Zygosity,
}

/// Genotype every candidate in one active window against the reads projected onto it.
///
/// `ref_window` is the uppercase reference over the window; `window_start` is the 1-based reference
/// coordinate of `ref_window[0]`. Each read's `site_obs` must be parallel to `candidates`.
pub fn genotype_window(
    ref_window: &[u8],
    window_start: i64,
    candidates: &[Candidate],
    reads: &[WindowRead],
    params: &ReassemblyParams,
) -> Vec<ReassemblyCall> {
    let mut hmm = PairHMM::new(&GapParams);
    candidates
        .iter()
        .enumerate()
        .map(|(ci, cand)| genotype_candidate(&mut hmm, ref_window, window_start, ci, cand, reads, params))
        .collect()
}

/// Resolve a single candidate: build the alt haplotype, select + dedup spanning reads, score each
/// against ref vs alt with the PairHMM, and genotype by the aggregate log-odds.
fn genotype_candidate(
    hmm: &mut PairHMM,
    ref_window: &[u8],
    window_start: i64,
    ci: usize,
    cand: &Candidate,
    reads: &[WindowRead],
    params: &ReassemblyParams,
) -> ReassemblyCall {
    // Alternate haplotype = reference window with the single substitution at the candidate.
    let off = (cand.position - window_start) as usize;
    let mut alt_hap = ref_window.to_vec();
    if off < alt_hap.len() {
        alt_hap[off] = cand.alt_base;
    }

    // Stage B — select reads that clear the MAPQ gate and span this candidate, then collapse
    // overlapping mate pairs to one fragment (keep the record whose site base has higher quality).
    let mut by_fragment: HashMap<&[u8], usize> = HashMap::new();
    for (ri, read) in reads.iter().enumerate() {
        if read.mapq < params.min_mapping_quality {
            continue;
        }
        let Some(Some(obs)) = read.site_obs.get(ci) else {
            continue; // does not span the site
        };
        by_fragment
            .entry(read.name.as_slice())
            .and_modify(|kept| {
                let kept_q = reads[*kept].site_obs[ci].map(|o| o.qual).unwrap_or(0);
                if obs.qual > kept_q {
                    *kept = ri;
                }
            })
            .or_insert(ri);
    }

    // Stages D/E — per-fragment likelihood ratio and site-base vote.
    let mut log_odds = 0.0f64;
    let mut depth = 0u32;
    let mut alt_depth = 0u32;
    for &ri in by_fragment.values() {
        let read = &reads[ri];
        let lp_ref = hap_likelihood(hmm, &read.seq, &read.quals, ref_window);
        let lp_alt = hap_likelihood(hmm, &read.seq, &read.quals, &alt_hap);
        log_odds += *lp_alt - *lp_ref;
        depth += 1;
        if read.site_obs[ci].map(|o| o.base) == Some(cand.alt_base) {
            alt_depth += 1;
        }
    }

    let allele_fraction = if depth > 0 { alt_depth as f64 / depth as f64 } else { 0.0 };
    let genotype = if log_odds > params.min_log_odds && alt_depth >= params.min_alt_fragments {
        Zygosity::Derived
    } else if log_odds < -params.min_log_odds {
        Zygosity::Ancestral
    } else {
        Zygosity::Ambiguous
    };

    ReassemblyCall {
        position: cand.position,
        ref_base: cand.ref_base,
        alt_base: cand.alt_base,
        depth,
        alt_depth,
        allele_fraction,
        log_odds,
        quality: (log_odds.abs() * PHRED_PER_NAT).min(99.0),
        genotype,
    }
}

/// Log-probability that `read` (with `quals`) was produced by `hap`, marginalised over alignments.
fn hap_likelihood(hmm: &mut PairHMM, read: &[u8], quals: &[u8], hap: &[u8]) -> LogProb {
    hmm.prob_related(&ReadHapEmission { read, quals, hap }, &Semiglobal, None)
}

// ---- base-quality-aware PairHMM emission model (POC-validated) --------------------------------

/// Phred score → error probability, clamped to Q2–Q60 (never a certain match/mismatch).
fn phred_err(q: u8) -> f64 {
    let q = q.clamp(2, 60) as f64;
    10f64.powf(-q / 10.0)
}

/// Emission: `x` = read (carries per-base quality), `y` = candidate haplotype.
struct ReadHapEmission<'a> {
    read: &'a [u8],
    quals: &'a [u8],
    hap: &'a [u8],
}

impl EmissionParameters for ReadHapEmission<'_> {
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        let err = phred_err(self.quals[i]);
        if self.read[i] == self.hap[j] {
            XYEmission::Match(LogProb::from(Prob(1.0 - err)))
        } else {
            XYEmission::Mismatch(LogProb::from(Prob(err / 3.0)))
        }
    }
    fn prob_emit_x(&self, _i: usize) -> LogProb {
        LogProb::ln_one() // insertion in read: base is real; cost is the gap-open prob
    }
    fn prob_emit_y(&self, _j: usize) -> LogProb {
        LogProb::ln_one() // deletion (gap in read): hap base emitted against a gap
    }
    fn len_x(&self) -> usize {
        self.read.len()
    }
    fn len_y(&self) -> usize {
        self.hap.len()
    }
}

/// GATK-ish affine gap model (indels rare relative to substitutions).
struct GapParams;
impl GapParameters for GapParams {
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(Prob(1e-4))
    }
    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(Prob(1e-4))
    }
    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::from(Prob(0.1))
    }
    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::from(Prob(0.1))
    }
}

/// Semiglobal in the read: free leading/trailing offset so window-edge trimming isn't penalised.
struct Semiglobal;
impl StartEndGapParameters for Semiglobal {
    fn free_start_gap_x(&self) -> bool {
        true
    }
    fn free_end_gap_x(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A 40 bp reference window; candidate sits at its centre.
    const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const WIN_START: i64 = 1000;
    const CAND_OFF: usize = 20;

    fn candidate() -> Candidate {
        Candidate {
            position: WIN_START + CAND_OFF as i64,
            ref_base: REF[CAND_OFF],
            alt_base: b'T', // REF[20] is 'A'; flip A>T
        }
    }

    /// A read whose window sequence carries `site_base` at the candidate offset; uniform quality.
    fn read(name: &str, site_base: u8, qual: u8, mapq: u8) -> WindowRead {
        let mut seq = REF.to_vec();
        seq[CAND_OFF] = site_base;
        WindowRead {
            name: name.as_bytes().to_vec(),
            seq,
            quals: vec![qual; REF.len()],
            mapq,
            site_obs: vec![Some(SiteObs { base: site_base, qual })],
        }
    }

    fn call(reads: &[WindowRead]) -> ReassemblyCall {
        genotype_window(REF, WIN_START, &[candidate()], reads, &ReassemblyParams::default())
            .pop()
            .unwrap()
    }

    #[test]
    fn clean_derived_site_is_called() {
        // Twelve fragments all carrying the alt allele → strongly DERIVED.
        let reads: Vec<_> = (0..12).map(|i| read(&format!("r{i}"), b'T', 35, 60)).collect();
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Derived);
        assert_eq!(c.depth, 12);
        assert_eq!(c.alt_depth, 12);
        assert!(c.log_odds > 2.0, "log_odds {}", c.log_odds);
        assert_eq!(c.alt_base, b'T');
    }

    #[test]
    fn low_mapq_paralog_reference_reads_are_dropped_recovering_the_site() {
        // The misaligned-ref case: 8 clean alt fragments (MAPQ 60) + 6 paralog reference fragments
        // that carry the ref base but are ambiguously placed (MAPQ 5). The MAPQ gate excludes the
        // paralogs, so the site is recovered as DERIVED instead of rejected as ~50/50.
        let mut reads: Vec<_> = (0..8).map(|i| read(&format!("alt{i}"), b'T', 35, 60)).collect();
        reads.extend((0..6).map(|i| read(&format!("par{i}"), b'A', 35, 5)));
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Derived);
        assert_eq!(c.depth, 8, "paralog low-MAPQ reads must be excluded");
        assert_eq!(c.alt_depth, 8);
    }

    #[test]
    fn genuinely_balanced_high_quality_site_is_not_called() {
        // Specificity: an even split of high-quality, well-placed ref and alt fragments is truly
        // undecided — reassembly must NOT invent a call.
        let mut reads: Vec<_> = (0..6).map(|i| read(&format!("alt{i}"), b'T', 35, 60)).collect();
        reads.extend((0..6).map(|i| read(&format!("ref{i}"), b'A', 35, 60)));
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Ambiguous);
    }

    #[test]
    fn overlapping_mates_are_counted_once() {
        // Four distinct alt fragments plus a read and its mate (same name) both covering the site.
        // Fragment dedup must collapse the mate pair so depth is 5, not 6.
        let mut reads: Vec<_> = (0..4).map(|i| read(&format!("f{i}"), b'T', 35, 60)).collect();
        reads.push(read("pair", b'T', 20, 60)); // read
        reads.push(read("pair", b'T', 35, 60)); // its mate (higher qual → the kept one)
        let c = call(&reads);
        assert_eq!(c.depth, 5, "overlapping mates double-counted");
        assert_eq!(c.alt_depth, 5);
        assert_eq!(c.genotype, Zygosity::Derived);
    }

    #[test]
    fn all_reference_reads_are_ancestral() {
        let reads: Vec<_> = (0..10).map(|i| read(&format!("r{i}"), b'A', 35, 60)).collect();
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Ancestral);
        assert_eq!(c.alt_depth, 0);
        assert!(c.log_odds < -2.0, "log_odds {}", c.log_odds);
    }
}

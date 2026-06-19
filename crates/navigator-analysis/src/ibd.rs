//! IBD (identity-by-descent) detection + relationship estimation — port of the Scala
//! `ibd.engine` (PairwiseIbdDetector / GeneticMap / RelationshipEstimator). Pure math:
//! the network matching layer (crypto/protocol/relay) is out of scope here.
//!
//! Input is per-chromosome diploid dosage genotypes (0/1/2, -1 no-call) — exactly what
//! [`crate::caller::genotype_sites`] produces. The detector classifies IBS at shared
//! sites, finds high-IBS runs with a sliding window + error tolerance, converts spans to
//! centiMorgans via a [`GeneticMap`], filters by length, and merges nearby segments.
//! Total shared cM maps to a relationship category (Shared cM Project / ISOGG values).

use std::collections::{BTreeSet, HashMap};

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// Relationship category from total shared cM (mirrors the Scala enum).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RelationshipEstimate {
    ParentChild,
    FullSibling,
    HalfSibling,
    Grandparent,
    AuntUncle,
    FirstCousin,
    FirstCousinOnceRemoved,
    SecondCousin,
    SecondCousinOnceRemoved,
    ThirdCousin,
    FourthCousin,
    FifthCousin,
    Distant,
    Unknown,
}

/// Estimate relationship from total shared cM (thresholds match the Scala estimator).
pub fn estimate_relationship(total_shared_cm: f64) -> RelationshipEstimate {
    use RelationshipEstimate::*;
    if total_shared_cm >= 3400.0 {
        ParentChild
    } else if total_shared_cm >= 2550.0 {
        FullSibling
    } else if total_shared_cm >= 1700.0 {
        Grandparent
    } else if total_shared_cm >= 1200.0 {
        AuntUncle
    } else if total_shared_cm >= 680.0 {
        FirstCousin
    } else if total_shared_cm >= 400.0 {
        FirstCousinOnceRemoved
    } else if total_shared_cm >= 200.0 {
        SecondCousin
    } else if total_shared_cm >= 90.0 {
        SecondCousinOnceRemoved
    } else if total_shared_cm >= 50.0 {
        ThirdCousin
    } else if total_shared_cm >= 25.0 {
        FourthCousin
    } else if total_shared_cm >= 10.0 {
        FifthCousin
    } else if total_shared_cm >= 7.0 {
        RelationshipEstimate::Distant
    } else {
        Unknown
    }
}

/// A detected IBD segment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IbdSegment {
    pub chromosome: String,
    pub start_position: i64,
    pub end_position: i64,
    pub length_cm: f64,
    pub snp_count: Option<u32>,
    pub is_half_identical: Option<bool>,
}

/// Aggregate of shared segments between two individuals.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MatchSummary {
    pub total_shared_cm: f64,
    pub segment_count: usize,
    pub longest_segment_cm: f64,
    pub relationship: RelationshipEstimate,
}

fn round_cm(cm: f64) -> f64 {
    (cm * 100.0).round() / 100.0
}

impl MatchSummary {
    pub fn from_segments(segments: &[IbdSegment]) -> Self {
        if segments.is_empty() {
            return MatchSummary {
                total_shared_cm: 0.0,
                segment_count: 0,
                longest_segment_cm: 0.0,
                relationship: RelationshipEstimate::Unknown,
            };
        }
        let total = round_cm(segments.iter().map(|s| s.length_cm).sum());
        let longest = round_cm(segments.iter().map(|s| s.length_cm).fold(f64::MIN, f64::max));
        MatchSummary {
            total_shared_cm: total,
            segment_count: segments.len(),
            longest_segment_cm: longest,
            relationship: estimate_relationship(total),
        }
    }
}

/// Diploid dosage genotypes for one chromosome. `dosages`: 0 hom-ref, 1 het, 2 hom-alt,
/// -1 no-call. `positions` must be sorted ascending and the same length as `dosages`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChromosomeGenotypes {
    pub chromosome: String,
    pub positions: Vec<i32>,
    pub dosages: Vec<i8>,
}

impl ChromosomeGenotypes {
    pub fn size(&self) -> usize {
        self.positions.len()
    }
}

/// Detector configuration (defaults match the Scala `IbdDetectorConfig`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IbdDetectorConfig {
    pub min_segment_cm: f64,
    pub min_snp_count: usize,
    pub window_size: usize,
    pub ibs_threshold: f64,
    pub error_tolerance: f64,
    pub min_gap_bp: i64,
}

impl Default for IbdDetectorConfig {
    fn default() -> Self {
        IbdDetectorConfig {
            min_segment_cm: 7.0,
            min_snp_count: 100,
            window_size: 100,
            ibs_threshold: 0.70,
            error_tolerance: 0.01,
            min_gap_bp: 1_000_000,
        }
    }
}

/// Normalize a chromosome name to its bare form ("chr1" -> "1", "chrX" -> "X").
pub fn normalize_chromosome(chr: &str) -> String {
    let s = chr.to_lowercase();
    let s = s.strip_prefix("chr").unwrap_or(&s);
    match s.parse::<i64>() {
        Ok(n) => n.to_string(),
        Err(_) => s.to_uppercase(),
    }
}

/// Sorted position/cM arrays for one chromosome; linear interpolation with end extrapolation.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
struct ChromosomeMap {
    positions: Vec<i32>,
    cm: Vec<f64>,
}

impl ChromosomeMap {
    fn interpolate(&self, position: i32) -> f64 {
        match self.positions.binary_search(&position) {
            Ok(idx) => self.cm[idx],
            Err(ins) => {
                let n = self.positions.len();
                if ins == 0 {
                    if n >= 2 {
                        let rate = (self.cm[1] - self.cm[0]) / (self.positions[1] - self.positions[0]) as f64;
                        self.cm[0] + rate * (position - self.positions[0]) as f64
                    } else {
                        self.cm[0]
                    }
                } else if ins >= n {
                    if n >= 2 {
                        let rate =
                            (self.cm[n - 1] - self.cm[n - 2]) / (self.positions[n - 1] - self.positions[n - 2]) as f64;
                        self.cm[n - 1] + rate * (position - self.positions[n - 1]) as f64
                    } else {
                        self.cm[n - 1]
                    }
                } else {
                    let (lo, hi) = (ins - 1, ins);
                    let frac =
                        (position - self.positions[lo]) as f64 / (self.positions[hi] - self.positions[lo]) as f64;
                    self.cm[lo] + frac * (self.cm[hi] - self.cm[lo])
                }
            }
        }
    }
}

/// bp -> cM genetic map, keyed by normalized chromosome.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GeneticMap {
    maps: HashMap<String, ChromosomeMap>,
}

impl GeneticMap {
    /// Deserialize a built genetic-map asset (bincode), as written by `panelbuild genetic-map`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("genetic map decode: {e}")))
    }

    /// Serialize to the binary asset form (bincode).
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("genetic map encode: {e}")))
    }

    /// Build from per-chromosome `(name, positions, cm)` marker arrays.
    pub fn from_markers(markers: impl IntoIterator<Item = (String, Vec<i32>, Vec<f64>)>) -> Self {
        let maps = markers
            .into_iter()
            .map(|(chr, positions, cm)| (normalize_chromosome(&chr), ChromosomeMap { positions, cm }))
            .collect();
        GeneticMap { maps }
    }

    /// A uniform `cm_per_mb` map over the given `(chromosome, length_bp)` pairs.
    pub fn uniform(cm_per_mb: f64, lengths: &[(&str, i32)]) -> Self {
        let markers = lengths.iter().map(|&(chr, len)| {
            (
                chr.to_string(),
                vec![1, len],
                vec![0.0, len as f64 / 1_000_000.0 * cm_per_mb],
            )
        });
        GeneticMap::from_markers(markers)
    }

    pub fn position_to_cm(&self, chromosome: &str, position: i32) -> Option<f64> {
        self.maps
            .get(&normalize_chromosome(chromosome))
            .map(|m| m.interpolate(position))
    }

    pub fn interval_cm(&self, chromosome: &str, start_bp: i32, end_bp: i32) -> Option<f64> {
        let s = self.position_to_cm(chromosome, start_bp)?;
        let e = self.position_to_cm(chromosome, end_bp)?;
        Some((e - s).abs())
    }

    pub fn has_chromosome(&self, chromosome: &str) -> bool {
        self.maps.contains_key(&normalize_chromosome(chromosome))
    }
}

/// IBS state between two diploid dosages: 2 = both alleles shared (same genotype),
/// 1 = one shared, 0 = opposite homozygotes.
fn ibs_state(g1: i8, g2: i8) -> i8 {
    match (g1 - g2).abs() {
        0 => 2,
        1 => 1,
        _ => 0,
    }
}

/// Pairwise IBD detector.
pub struct PairwiseIbdDetector {
    pub config: IbdDetectorConfig,
}

impl PairwiseIbdDetector {
    pub fn new(config: IbdDetectorConfig) -> Self {
        PairwiseIbdDetector { config }
    }

    /// Detect IBD segments across all chromosomes shared by both samples.
    pub fn detect_segments(
        &self,
        sample1: &HashMap<String, ChromosomeGenotypes>,
        sample2: &HashMap<String, ChromosomeGenotypes>,
        genetic_map: &GeneticMap,
    ) -> Vec<IbdSegment> {
        let shared: BTreeSet<&String> = sample1.keys().filter(|k| sample2.contains_key(*k)).collect();
        let mut out = Vec::new();
        for chr in shared {
            out.extend(self.detect_chromosome_segments(&sample1[chr], &sample2[chr], genetic_map));
        }
        out
    }

    /// Detect IBD segments on a single chromosome.
    pub fn detect_chromosome_segments(
        &self,
        geno1: &ChromosomeGenotypes,
        geno2: &ChromosomeGenotypes,
        genetic_map: &GeneticMap,
    ) -> Vec<IbdSegment> {
        let (aligned_pos, g1, g2) = intersect_positions(geno1, geno2);
        if aligned_pos.len() < self.config.min_snp_count {
            return Vec::new();
        }
        let ibs: Vec<i8> = g1.iter().zip(&g2).map(|(&a, &b)| ibs_state(a, b)).collect();
        let candidates = self.find_candidate_segments(&aligned_pos, &ibs);
        let merged = self.merge_segments(candidates);

        let mut out = Vec::new();
        for (start, end, snp_count, _ibs2) in merged {
            let start_bp = aligned_pos[start];
            let end_bp = aligned_pos[end];
            if let Some(cm) = genetic_map.interval_cm(&geno1.chromosome, start_bp, end_bp) {
                if cm >= self.config.min_segment_cm && snp_count >= self.config.min_snp_count {
                    out.push(IbdSegment {
                        chromosome: normalize_chromosome(&geno1.chromosome),
                        start_position: start_bp as i64,
                        end_position: end_bp as i64,
                        length_cm: round_cm(cm),
                        snp_count: Some(snp_count as u32),
                        is_half_identical: Some(true), // IBS-based detection finds IBD1
                    });
                }
            }
        }
        out
    }

    /// Sliding-window candidate segments: `(start_idx, end_idx, snp_count, ibs2_count)`.
    fn find_candidate_segments(&self, positions: &[i32], ibs: &[i8]) -> Vec<(usize, usize, usize, usize)> {
        let n = positions.len();
        if n < self.config.window_size {
            return Vec::new();
        }
        let mut candidates = Vec::new();
        let (mut in_segment, mut seg_start) = (false, 0usize);
        let (mut ibs0, mut ibs2, mut seg_snp) = (0usize, 0usize, 0usize);
        let half = self.config.window_size / 2;

        for i in 0..n {
            let look_back = i.saturating_sub(half);
            let look_forward = (i + half).min(n - 1);
            let (mut local_ibs2, mut local_ibs0, mut local_total) = (0usize, 0usize, 0usize);
            for &s in &ibs[look_back..=look_forward] {
                match s {
                    2 => local_ibs2 += 1,
                    0 => local_ibs0 += 1,
                    _ => {}
                }
                local_total += 1;
            }
            let ibs_fraction = if local_total > 0 {
                (local_ibs2 as f64 + 0.5 * (local_total - local_ibs2 - local_ibs0) as f64) / local_total as f64
            } else {
                0.0
            };

            if ibs_fraction >= self.config.ibs_threshold {
                if !in_segment {
                    in_segment = true;
                    seg_start = i;
                    ibs0 = 0;
                    ibs2 = 0;
                    seg_snp = 0;
                }
                seg_snp += 1;
                if ibs[i] == 0 {
                    ibs0 += 1;
                }
                if ibs[i] == 2 {
                    ibs2 += 1;
                }
            } else if in_segment {
                let error_rate = if seg_snp > 0 { ibs0 as f64 / seg_snp as f64 } else { 1.0 };
                if error_rate <= self.config.error_tolerance || seg_snp < self.config.window_size {
                    seg_snp += 1;
                    if ibs[i] == 0 {
                        ibs0 += 1;
                    }
                    if ibs[i] == 2 {
                        ibs2 += 1;
                    }
                } else {
                    if seg_snp >= self.config.min_snp_count {
                        candidates.push((seg_start, i - 1, seg_snp, ibs2));
                    }
                    in_segment = false;
                }
            }
        }
        if in_segment && seg_snp >= self.config.min_snp_count {
            candidates.push((seg_start, n - 1, seg_snp, ibs2));
        }
        candidates
    }

    /// Merge candidate segments separated by gaps `<= window_size` (in SNP index).
    fn merge_segments(&self, candidates: Vec<(usize, usize, usize, usize)>) -> Vec<(usize, usize, usize, usize)> {
        if candidates.len() <= 1 {
            return candidates;
        }
        let mut merged: Vec<(usize, usize, usize, usize)> = vec![candidates[0]];
        for &(cur_start, cur_end, cur_snps, cur_ibs2) in &candidates[1..] {
            let last = *merged.last().unwrap();
            let (_, prev_end, prev_snps, prev_ibs2) = last;
            if cur_start.saturating_sub(prev_end) <= self.config.window_size {
                let merged_snps = prev_snps + cur_snps + (cur_start - prev_end);
                let i = merged.len() - 1;
                merged[i] = (last.0, cur_end, merged_snps, prev_ibs2 + cur_ibs2);
            } else {
                merged.push((cur_start, cur_end, cur_snps, cur_ibs2));
            }
        }
        merged
    }
}

/// Intersect two sorted genotype arrays to shared positions where both have valid calls.
fn intersect_positions(geno1: &ChromosomeGenotypes, geno2: &ChromosomeGenotypes) -> (Vec<i32>, Vec<i8>, Vec<i8>) {
    let (mut pos, mut g1, mut g2) = (Vec::new(), Vec::new(), Vec::new());
    let (mut i, mut j) = (0usize, 0usize);
    while i < geno1.size() && j < geno2.size() {
        let (p1, p2) = (geno1.positions[i], geno2.positions[j]);
        if p1 == p2 {
            if geno1.dosages[i] >= 0 && geno2.dosages[j] >= 0 {
                pos.push(p1);
                g1.push(geno1.dosages[i]);
                g2.push(geno2.dosages[j]);
            }
            i += 1;
            j += 1;
        } else if p1 < p2 {
            i += 1;
        } else {
            j += 1;
        }
    }
    (pos, g1, g2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relationship_thresholds() {
        assert_eq!(estimate_relationship(3500.0), RelationshipEstimate::ParentChild);
        assert_eq!(estimate_relationship(2600.0), RelationshipEstimate::FullSibling);
        assert_eq!(estimate_relationship(700.0), RelationshipEstimate::FirstCousin);
        assert_eq!(estimate_relationship(8.0), RelationshipEstimate::Distant);
        assert_eq!(estimate_relationship(3.0), RelationshipEstimate::Unknown);
    }

    #[test]
    fn ibs_state_classification() {
        assert_eq!(ibs_state(0, 0), 2); // AA vs AA
        assert_eq!(ibs_state(1, 1), 2); // AB vs AB
        assert_eq!(ibs_state(0, 1), 1); // AA vs AB
        assert_eq!(ibs_state(0, 2), 0); // AA vs BB
        assert_eq!(ibs_state(2, 0), 0);
    }

    #[test]
    fn normalize_chromosome_forms() {
        assert_eq!(normalize_chromosome("chr1"), "1");
        assert_eq!(normalize_chromosome("1"), "1");
        assert_eq!(normalize_chromosome("chrX"), "X");
        assert_eq!(normalize_chromosome("01"), "1");
    }

    #[test]
    fn genetic_map_interpolates_and_extrapolates() {
        let m = GeneticMap::from_markers([("1".to_string(), vec![1_000, 2_000, 3_000], vec![0.0, 1.0, 2.0])]);
        assert_eq!(m.position_to_cm("chr1", 2_000), Some(1.0)); // exact
        assert_eq!(m.position_to_cm("1", 1_500), Some(0.5)); // midpoint
        assert_eq!(m.interval_cm("1", 1_000, 3_000), Some(2.0));
        // extrapolation past the last marker uses the last rate (1 cM / 1000 bp).
        assert_eq!(m.position_to_cm("1", 4_000), Some(3.0));
        assert!(m.position_to_cm("9", 1).is_none());
    }

    #[test]
    fn genetic_map_round_trips_through_bincode() {
        let m = GeneticMap::from_markers([
            ("1".to_string(), vec![1_000, 2_000, 3_000], vec![0.0, 1.0, 2.5]),
            ("X".to_string(), vec![5_000, 9_000], vec![0.0, 4.0]),
        ]);
        let back = GeneticMap::from_bytes(&m.to_bytes().unwrap()).unwrap();
        assert_eq!(back, m);
        // The reloaded map interpolates identically.
        assert_eq!(back.position_to_cm("chr1", 1_500), Some(0.5));
        assert_eq!(back.position_to_cm("X", 7_000), Some(2.0));
    }

    /// Build a chromosome of `n` SNPs spaced `step` bp apart, dosages from `f(i)`.
    fn chrom(n: usize, step: i32, f: impl Fn(usize) -> i8) -> ChromosomeGenotypes {
        ChromosomeGenotypes {
            chromosome: "1".into(),
            positions: (0..n).map(|i| 1 + i as i32 * step).collect(),
            dosages: (0..n).map(f).collect(),
        }
    }

    #[test]
    fn detects_a_long_identical_segment_but_not_a_discordant_region() {
        // 300 SNPs at 50 kb spacing (15 Mb -> 15 cM at 1 cM/Mb). First 250 identical
        // (IBS-2), last 50 opposite homozygotes (IBS-0).
        let n = 300;
        let s1 = chrom(n, 50_000, |i| (i % 3) as i8); // 0,1,2 cycling
        let s2 = chrom(n, 50_000, |i| if i < 250 { (i % 3) as i8 } else { 2 - (i % 3) as i8 });
        let map = GeneticMap::uniform(1.0, &[("1", 16_000_000)]);
        let det = PairwiseIbdDetector::new(IbdDetectorConfig::default());

        let segs = det.detect_chromosome_segments(&s1, &s2, &map);
        assert_eq!(segs.len(), 1, "expected one segment, got {segs:?}");
        let seg = &segs[0];
        assert_eq!(seg.chromosome, "1");
        assert!(seg.length_cm >= 7.0, "length {} below threshold", seg.length_cm);
        assert!(seg.snp_count.unwrap() >= 100);
        assert_eq!(seg.is_half_identical, Some(true));

        // Two unrelated samples (always opposite homozygotes) -> no segments.
        let u1 = chrom(n, 50_000, |_| 0);
        let u2 = chrom(n, 50_000, |_| 2);
        assert!(det.detect_chromosome_segments(&u1, &u2, &map).is_empty());
    }

    #[test]
    fn match_summary_aggregates_and_estimates() {
        let seg = |cm| IbdSegment {
            chromosome: "1".into(),
            start_position: 1,
            end_position: 2,
            length_cm: cm,
            snp_count: Some(200),
            is_half_identical: Some(true),
        };
        let s = MatchSummary::from_segments(&[seg(120.0), seg(90.0)]);
        assert_eq!(s.segment_count, 2);
        assert_eq!(s.total_shared_cm, 210.0);
        assert_eq!(s.longest_segment_cm, 120.0);
        assert_eq!(s.relationship, RelationshipEstimate::SecondCousin); // 210 >= 200
        assert_eq!(
            MatchSummary::from_segments(&[]).relationship,
            RelationshipEstimate::Unknown
        );
    }
}

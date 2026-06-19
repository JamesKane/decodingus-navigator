//! Depth segmenter — port of the Scala `DepthSegmenter`. Z-score CNV detection: expected
//! reads/bin from coverage, Poisson-ish variance, greedy extension of aberrant runs
//! (tolerating short dips), size filter, then merge + conversion to SV calls.

use super::evidence::DepthSegment;
use super::types::{SvCall, SvCallerConfig, SvType};
use std::collections::BTreeMap;

/// Segment per-contig depth bins into CNV `DepthSegment`s.
pub fn segment(
    depth_bins: &BTreeMap<String, Vec<u32>>,
    contig_lengths: &BTreeMap<String, i64>,
    mean_coverage: f64,
    read_length: f64,
    config: &SvCallerConfig,
) -> Vec<DepthSegment> {
    let bin_size = config.bin_size;
    let min_z = config.min_depth_z_score;
    let expected_per_bin = mean_coverage * bin_size as f64 / read_length;

    let mut segments: Vec<DepthSegment> = Vec::new();

    for (contig, bins) in depth_bins {
        let contig_length = contig_lengths
            .get(contig)
            .copied()
            .unwrap_or(bins.len() as i64 * bin_size);

        let z_scores: Vec<f64> = bins
            .iter()
            .map(|&rc| {
                if expected_per_bin > 0.0 {
                    let variance = expected_per_bin.max(1.0);
                    (rc as f64 - expected_per_bin) / variance.sqrt()
                } else {
                    0.0
                }
            })
            .collect();

        let mut i = 0usize;
        while i < z_scores.len() {
            let z = z_scores[i];
            if z.abs() < min_z {
                i += 1;
                continue;
            }
            // Start of an aberrant run.
            let is_deletion = z < 0.0;
            let start_bin = i;
            let mut count = 1usize;
            let mut sum_z = z;
            let mut sum_depth = bins[i] as f64;

            // Extend while the next bin is on the same side, or short dips are bridged.
            loop {
                let end_bin = start_bin + count - 1;
                if end_bin + 1 >= z_scores.len() {
                    break;
                }
                let next_z = z_scores[end_bin + 1];
                let same_side = (is_deletion && next_z < -min_z * 0.5) || (!is_deletion && next_z > min_z * 0.5);
                if same_side {
                    count += 1;
                    sum_z += next_z;
                    sum_depth += bins[end_bin + 1] as f64;
                    continue;
                }
                // Look ahead up to 3 bins; bridge if >= 2 are aberrant.
                let look_ahead = (end_bin + 3).min(z_scores.len() - 1);
                let future_aberrant = ((end_bin + 1)..=look_ahead)
                    .filter(|&j| {
                        let fz = z_scores[j];
                        (is_deletion && fz < -min_z * 0.5) || (!is_deletion && fz > min_z * 0.5)
                    })
                    .count();
                if future_aberrant >= 2 {
                    count += 1;
                    sum_z += next_z;
                    sum_depth += bins[end_bin + 1] as f64;
                } else {
                    break;
                }
            }

            let end_bin = start_bin + count - 1;
            let segment_start = start_bin as i64 * bin_size;
            let segment_end = ((end_bin + 1) as i64 * bin_size).min(contig_length);
            if segment_end - segment_start >= config.min_cnv_size {
                let mean_depth = sum_depth / count as f64;
                let mean_z = sum_z / count as f64;
                let log2_ratio = if expected_per_bin > 0.0 {
                    (mean_depth / expected_per_bin).log2()
                } else {
                    0.0
                };
                segments.push(DepthSegment {
                    chrom: contig.clone(),
                    start: segment_start,
                    end: segment_end,
                    mean_depth,
                    log2_ratio,
                    z_score: mean_z,
                    num_bins: count as u32,
                    sv_type: if mean_z < 0.0 { SvType::Del } else { SvType::Dup },
                });
            }
            i = start_bin + count;
        }
    }

    segments.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
    segments
}

/// Merge nearby same-type segments (default gap 50 kb), weighting by bin count.
pub fn merge_nearby_segments(segments: &[DepthSegment], max_gap: i64) -> Vec<DepthSegment> {
    if segments.is_empty() {
        return Vec::new();
    }
    let mut sorted = segments.to_vec();
    sorted.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));

    let mut merged: Vec<DepthSegment> = Vec::new();
    let mut current = sorted[0].clone();
    for next in sorted.into_iter().skip(1) {
        if current.chrom == next.chrom && current.sv_type == next.sv_type && next.start - current.end <= max_gap {
            let total = current.num_bins + next.num_bins;
            let w = |a: f64, b: f64| (a * current.num_bins as f64 + b * next.num_bins as f64) / total as f64;
            current = DepthSegment {
                chrom: current.chrom.clone(),
                start: current.start,
                end: next.end,
                mean_depth: w(current.mean_depth, next.mean_depth),
                log2_ratio: w(current.log2_ratio, next.log2_ratio),
                z_score: w(current.z_score, next.z_score),
                num_bins: total,
                sv_type: current.sv_type,
            };
        } else {
            merged.push(current);
            current = next;
        }
    }
    merged.push(current);
    merged
}

/// Convert depth segments to SV calls (depth-only: no PE/SR support).
pub fn to_sv_calls(segments: &[DepthSegment], config: &SvCallerConfig) -> Vec<SvCall> {
    segments
        .iter()
        .enumerate()
        .map(|(idx, seg)| {
            let quality = (seg.z_score.abs() * 10.0).min(99.0);
            let ci_size = (config.bin_size / 2).max(100) as i32;
            let genotype = match seg.sv_type {
                SvType::Del if seg.log2_ratio < -0.9 => "1/1",
                SvType::Dup if seg.log2_ratio > 0.7 => "1/1",
                _ => "0/1",
            };
            let sv_len = if seg.sv_type == SvType::Del {
                -(seg.end - seg.start)
            } else {
                seg.end - seg.start
            };
            SvCall {
                id: format!("CNV_{}_{}_{}", seg.chrom, seg.start, idx),
                chrom: seg.chrom.clone(),
                start: seg.start,
                end: seg.end,
                sv_type: seg.sv_type,
                sv_len,
                ci_pos: (-ci_size, ci_size),
                ci_end: (-ci_size, ci_size),
                quality,
                paired_end_support: 0,
                split_read_support: 0,
                relative_depth: Some(2f64.powf(seg.log2_ratio)),
                mate_chrom: None,
                mate_pos: None,
                filter: if quality >= config.min_quality {
                    "PASS".into()
                } else {
                    "LowQual".into()
                },
                genotype: genotype.into(),
            }
        })
        .collect()
}

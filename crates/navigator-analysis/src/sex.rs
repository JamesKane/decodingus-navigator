//! Sex inference — Rust port of the Scala `SexInference`. Infers biological sex from
//! the chrX:autosome coverage ratio: males (XY) sit near 0.5×, females (XX) near 1.0×.
//! Drives per-contig ploidy for variant calling.
//!
//! Uses the **BAI index metadata** (per-reference aligned-record counts) — the Scala
//! fast path — so it is O(contigs), not a read scan. An unindexed BAM is an error,
//! matching the Scala behaviour.

use std::path::Path;

use noodles::bam;
use noodles::csi::binning_index::ReferenceSequence as _;

use crate::contig;
use crate::error::AnalysisError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InferredSex {
    Male,
    Female,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Confidence {
    High,
    Medium,
    Low,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SexInferenceResult {
    pub inferred_sex: InferredSex,
    pub x_autosome_ratio: f64,
    /// Autosome reads per 100 bp (the Scala "coverage" proxy).
    pub autosome_mean_coverage: f64,
    pub x_coverage: f64,
    pub confidence: Confidence,
}

// Thresholds (identical to the Scala constants).
const MALE_RATIO_THRESHOLD: f64 = 0.65;
const FEMALE_RATIO_THRESHOLD: f64 = 0.85;
const MIN_AUTOSOME_COVERAGE: f64 = 5.0;

/// Infer sex from an indexed BAM by comparing chrX to autosome read density.
pub fn infer_from_bam(bam_path: &Path) -> Result<SexInferenceResult, AnalysisError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .map_err(|e| AnalysisError::io(bam_path, e))?;
    let header = reader
        .read_header()
        .map_err(|e| AnalysisError::io(bam_path, e))?;

    // BAI reference sequences are in header order; zip names/lengths with counts.
    let bai_path = bam_path.with_extension("bam.bai");
    let index = bam::bai::read(&bai_path)
        .map_err(|e| AnalysisError::io(&bai_path, e))?;
    let counts: Vec<u64> = index
        .reference_sequences()
        .iter()
        .map(|rs| rs.metadata().map_or(0, |m| m.mapped_record_count()))
        .collect();

    let mut autosome_reads: u64 = 0;
    let mut autosome_length: u64 = 0;
    let mut x_reads: u64 = 0;
    let mut x_length: Option<u64> = None;

    for (i, (name_bytes, map)) in header.reference_sequences().iter().enumerate() {
        let name = String::from_utf8_lossy(name_bytes.as_ref());
        let length = map.length().get() as u64;
        let count = counts.get(i).copied().unwrap_or(0);
        if contig::is_autosome(&name) {
            autosome_reads += count;
            autosome_length += length;
        } else if contig::is_chr_x(&name) {
            x_reads += count;
            x_length = Some(length);
        }
    }

    if autosome_length == 0 {
        return Err(AnalysisError::Message(
            "no autosomal chromosomes found in BAM header".into(),
        ));
    }
    let Some(x_length) = x_length.filter(|&l| l > 0) else {
        return Err(AnalysisError::Message("chrX not found in BAM header".into()));
    };
    if autosome_reads == 0 {
        return Err(AnalysisError::Message(
            "no autosomal reads found - cannot infer sex".into(),
        ));
    }

    // Reads per 100 bp, then the chrX:autosome ratio.
    let autosome_coverage = autosome_reads as f64 / autosome_length as f64 * 100.0;
    let x_coverage = x_reads as f64 / x_length as f64 * 100.0;
    let ratio = if autosome_coverage > 0.0 {
        x_coverage / autosome_coverage
    } else {
        0.0
    };

    let (inferred_sex, confidence) = determine_sex(ratio, autosome_coverage);
    Ok(SexInferenceResult {
        inferred_sex,
        x_autosome_ratio: ratio,
        autosome_mean_coverage: autosome_coverage,
        x_coverage,
        confidence,
    })
}

/// Classify the ratio into sex + confidence (pure; mirrors the Scala `determineSex`).
pub fn determine_sex(ratio: f64, autosome_coverage: f64) -> (InferredSex, Confidence) {
    if autosome_coverage < MIN_AUTOSOME_COVERAGE {
        if ratio < MALE_RATIO_THRESHOLD {
            (InferredSex::Male, Confidence::Low)
        } else if ratio > FEMALE_RATIO_THRESHOLD {
            (InferredSex::Female, Confidence::Low)
        } else {
            (InferredSex::Unknown, Confidence::Low)
        }
    } else if ratio < MALE_RATIO_THRESHOLD {
        let conf = if ratio < 0.55 { Confidence::High } else { Confidence::Medium };
        (InferredSex::Male, conf)
    } else if ratio > FEMALE_RATIO_THRESHOLD {
        let conf = if ratio > 0.95 { Confidence::High } else { Confidence::Medium };
        (InferredSex::Female, conf)
    } else {
        (InferredSex::Unknown, Confidence::Low)
    }
}

/// Ploidy for a contig given inferred sex; `None` means skip the contig (chrY in
/// females). Mirrors the Scala `ploidyForContig`.
pub fn ploidy_for_contig(contig_name: &str, sex: InferredSex) -> Option<u32> {
    if contig::is_chr_x(contig_name) {
        match sex {
            InferredSex::Female => Some(2),
            InferredSex::Male => Some(1),
            InferredSex::Unknown => Some(2),
        }
    } else if contig::is_chr_y(contig_name) {
        match sex {
            InferredSex::Female => None,
            InferredSex::Male => Some(1),
            InferredSex::Unknown => Some(1),
        }
    } else if contig::is_chr_m(contig_name) {
        Some(1)
    } else {
        Some(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn determine_sex_covers_all_branches() {
        // High coverage, clear signals.
        assert_eq!(determine_sex(0.50, 30.0), (InferredSex::Male, Confidence::High));
        assert_eq!(determine_sex(0.60, 30.0), (InferredSex::Male, Confidence::Medium));
        assert_eq!(determine_sex(1.00, 30.0), (InferredSex::Female, Confidence::High));
        assert_eq!(determine_sex(0.90, 30.0), (InferredSex::Female, Confidence::Medium));
        assert_eq!(determine_sex(0.75, 30.0), (InferredSex::Unknown, Confidence::Low));
        // Low coverage -> always low confidence.
        assert_eq!(determine_sex(0.50, 2.0), (InferredSex::Male, Confidence::Low));
        assert_eq!(determine_sex(1.00, 2.0), (InferredSex::Female, Confidence::Low));
        assert_eq!(determine_sex(0.75, 2.0), (InferredSex::Unknown, Confidence::Low));
    }

    #[test]
    fn ploidy_follows_sex() {
        assert_eq!(ploidy_for_contig("chrX", InferredSex::Male), Some(1));
        assert_eq!(ploidy_for_contig("chrX", InferredSex::Female), Some(2));
        assert_eq!(ploidy_for_contig("chrY", InferredSex::Female), None);
        assert_eq!(ploidy_for_contig("chrY", InferredSex::Male), Some(1));
        assert_eq!(ploidy_for_contig("chrM", InferredSex::Female), Some(1));
        assert_eq!(ploidy_for_contig("chr7", InferredSex::Male), Some(2));
    }
}

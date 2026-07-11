//! Sex inference — Rust port of the Scala `SexInference`. Infers biological sex from
//! the chrX:autosome coverage ratio: males (XY) sit near 0.5×, females (XX) near 1.0×.
//! Drives per-contig ploidy for variant calling.
//!
//! For BAM, uses the **BAI index metadata** (per-reference aligned-record counts) — the
//! Scala fast path — so it is O(contigs), not a read scan; an unindexed BAM is an error.
//! CRAM indexes (`.crai`) carry no per-reference counts, so CRAM falls back to a single
//! record scan tallying mapped reads per chromosome (O(reads), `reference` required).

use std::path::Path;

use noodles::bam;
use noodles::csi::binning_index::ReferenceSequence as _;
use noodles::sam;

use serde::{Deserialize, Serialize};

use crate::contig;
use crate::error::AnalysisError;
use crate::reader::{self, Format};
use crate::readview::AlnRead;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InferredSex {
    Male,
    Female,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Confidence {
    High,
    Medium,
    Low,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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

/// Per-chromosome-class accumulators: (autosome reads, autosome length, chrX reads, chrX length).
type Tally = (u64, u64, u64, Option<u64>);

/// Infer sex from an indexed BAM or CRAM by comparing chrX to autosome read density.
/// `reference` is required for CRAM (the record-scan fallback decodes it).
pub fn infer_from_bam(bam_path: &Path, reference: Option<&Path>) -> Result<SexInferenceResult, AnalysisError> {
    let tally = match reader::detect_format(bam_path) {
        Format::Bam => tally_via_bai(bam_path)?,
        Format::Cram => tally_via_scan(bam_path, reference)?,
    };
    result_from_tally(tally)
}

/// Turn a per-class read/length tally `(autosome_reads, autosome_length, x_reads, x_length)`
/// into the inferred-sex result (reads-per-100bp → chrX:autosome ratio → classification).
/// Shared by every tally source (BAI fast path, CRAM scan, and the parallel walker's summed
/// per-contig counts).
pub(crate) fn result_from_tally(tally: Tally) -> Result<SexInferenceResult, AnalysisError> {
    let (autosome_reads, autosome_length, x_reads, x_length) = tally;

    if autosome_length == 0 {
        return Err(AnalysisError::Message(
            "no autosomal chromosomes found in alignment header".into(),
        ));
    }
    let Some(x_length) = x_length.filter(|&l| l > 0) else {
        return Err(AnalysisError::Message("chrX not found in alignment header".into()));
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

/// Per-chromosome-class read tally over a record stream, shared by the standalone CRAM
/// scan and the fused [`crate::unified`] walker (which already touches every record, so it
/// tallies sex directly rather than depending on BAI). Build with [`SexState::new`] from
/// the header, feed every record via [`SexState::accept`], then [`SexState::finish`].
pub(crate) struct SexState {
    /// ref_id -> class: 0 = other, 1 = autosome, 2 = chrX.
    class: Vec<u8>,
    autosome_length: u64,
    x_length: Option<u64>,
    autosome_reads: u64,
    x_reads: u64,
}

impl SexState {
    pub(crate) fn new(header: &sam::Header) -> Self {
        let mut class = Vec::with_capacity(header.reference_sequences().len());
        let (mut autosome_length, mut x_length) = (0u64, None);
        for (name_bytes, map) in header.reference_sequences() {
            let name = String::from_utf8_lossy(name_bytes.as_ref());
            let length = map.length().get() as u64;
            if contig::is_autosome(&name) {
                autosome_length += length;
                class.push(1u8);
            } else if contig::is_chr_x(&name) {
                x_length = Some(length);
                class.push(2u8);
            } else {
                class.push(0u8);
            }
        }
        SexState {
            class,
            autosome_length,
            x_length,
            autosome_reads: 0,
            x_reads: 0,
        }
    }

    pub(crate) fn accept(&mut self, record: &impl AlnRead) {
        if record.flags().is_unmapped() {
            return;
        }
        if let Some(id) = record.reference_sequence_id() {
            match self.class.get(id).copied().unwrap_or(0) {
                1 => self.autosome_reads += 1,
                2 => self.x_reads += 1,
                _ => {}
            }
        }
    }

    fn tally(&self) -> Tally {
        (self.autosome_reads, self.autosome_length, self.x_reads, self.x_length)
    }

    pub(crate) fn finish(&self) -> Result<SexInferenceResult, AnalysisError> {
        result_from_tally(self.tally())
    }
}

/// BAM fast path: per-reference mapped-record counts from the BAI metadata (O(contigs)).
fn tally_via_bai(bam_path: &Path) -> Result<Tally, AnalysisError> {
    let header = reader::read_header(bam_path, None)?;
    let bai_path = bam_path.with_extension("bam.bai");
    let index = bam::bai::fs::read(&bai_path).map_err(|e| AnalysisError::io(&bai_path, e))?;
    let counts: Vec<u64> = index
        .reference_sequences()
        .iter()
        .map(|rs| rs.metadata().map_or(0, |m| m.mapped_record_count()))
        .collect();

    let (mut autosome_reads, mut autosome_length, mut x_reads, mut x_length) = (0u64, 0u64, 0u64, None);
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
    Ok((autosome_reads, autosome_length, x_reads, x_length))
}

/// CRAM fallback: a single record scan tallying mapped reads per chromosome class (CRAI
/// has no per-reference counts). Lengths come from the header; reads from `reference_sequence_id`.
fn tally_via_scan(bam_path: &Path, reference: Option<&Path>) -> Result<Tally, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, reference)?;
    let mut state = SexState::new(&header);
    for result in reader.records_lazy(&header) {
        state.accept(&result?);
    }
    Ok(state.tally())
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
        let conf = if ratio < 0.55 {
            Confidence::High
        } else {
            Confidence::Medium
        };
        (InferredSex::Male, conf)
    } else if ratio > FEMALE_RATIO_THRESHOLD {
        let conf = if ratio > 0.95 {
            Confidence::High
        } else {
            Confidence::Medium
        };
        (InferredSex::Female, conf)
    } else {
        (InferredSex::Unknown, Confidence::Low)
    }
}

/// Minimum chrY reads before an alignment can be judged Y-scoped — guards against calling a
/// near-empty file "male" off a handful of stray reads.
const Y_SCOPED_MIN_Y_READS: u64 = 1_000;
/// chrY reads must exceed this multiple of the autosome + chrX read count for an alignment to
/// count as Y-scoped. A whole-genome male carries ~100× MORE autosome than chrY reads (the
/// autosomes are ~40× the sequence and diploid), so a genuine Y-only extract — chrY in the
/// millions, autosomes a few dozen mismapped reads — clears this by orders of magnitude while
/// WGS and genuine females never approach it.
const Y_SCOPED_DOMINANCE: u64 = 8;

/// Does this alignment's per-contig read distribution look **Y-scoped** — chrY carrying
/// essentially all the reads while the autosomes and chrX hold only a trace of mismapped ones?
/// That is the shape of a Y-only extract (e.g. GRCh38 chrY reads realigned to hs1) or a Y-Elite /
/// Big Y capture. The chrX:autosome ratio [`determine_sex`] relies on is meaningless for such data
/// and can read as **female** — which would silently disable the entire Y pipeline (the Y
/// haplogroup step skips females *before* it ever downloads the tree). A `true` here means the
/// donor sequenced his Y, so callers should treat him as male regardless of the ratio.
///
/// The discriminator is read *counts*, not depth: pass `(contig_name, mapped_reads)` per contig.
pub fn is_y_scoped<'a>(per_contig_reads: impl IntoIterator<Item = (&'a str, u64)>) -> bool {
    let (mut y_reads, mut other_reads) = (0u64, 0u64);
    for (name, reads) in per_contig_reads {
        if contig::is_chr_y(name) {
            y_reads += reads;
        } else if contig::is_autosome(name) || contig::is_chr_x(name) {
            other_reads += reads;
        }
    }
    y_reads >= Y_SCOPED_MIN_Y_READS && y_reads > other_reads.saturating_mul(Y_SCOPED_DOMINANCE)
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
    fn y_scoped_detects_y_only_extracts() {
        // chrY in the millions, autosomes only a few dozen mismapped reads, no chrX → Y-scoped.
        assert!(is_y_scoped([("chrY", 3_000_000), ("chr1", 30), ("chr2", 24), ("chr7", 12)]));
        // A pure chrY-only alignment (nothing elsewhere) → Y-scoped.
        assert!(is_y_scoped([("chrY", 2_000_000)]));
        // chrY + chrM only (the chrYM.cram shape) → Y-scoped (chrM is neither autosome nor chrX).
        assert!(is_y_scoped([("chrY", 2_000_000), ("chrM", 50_000)]));
    }

    #[test]
    fn y_scoped_rejects_wgs_and_females() {
        // Male WGS: autosomes dwarf chrY → not Y-scoped (the ratio walk handles these).
        assert!(!is_y_scoped([("chr1", 200_000_000), ("chrX", 5_000_000), ("chrY", 3_000_000)]));
        // Female WGS: chrY only a trace of mismapping → not Y-scoped.
        assert!(!is_y_scoped([("chr1", 200_000_000), ("chrX", 10_000_000), ("chrY", 300)]));
        // Near-empty alignment: a handful of chrY reads is not enough to judge.
        assert!(!is_y_scoped([("chrY", 50)]));
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

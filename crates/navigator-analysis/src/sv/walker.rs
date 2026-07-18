//! SV evidence walker — port of the Scala `SvEvidenceWalker`. Single pass over the alignment
//! collecting per-bin read depth (CNV), discordant read pairs (BreakDancer-style), and
//! split reads from the SA tag (Pindel-style).
//!
//! Reads through [`crate::reader::open_seq`], so it handles BAM and CRAM alike. It walks records as
//! [`AlnRead`] views rather than a concrete record type: the BAM path stays on the lazy, zero-copy
//! `bam::Record` (this is a whole-genome pass, so a per-read owned copy would be costly), while the
//! CRAM path gets the decoded record it has no cheaper form of.

use std::collections::BTreeMap;
use std::path::Path;

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;

use super::evidence::{DiscordantPair, DiscordantReason, SplitRead, SvEvidenceCollection};
use super::types::SvCallerConfig;
use crate::error::AnalysisError;
use crate::readview::AlnRead;

const SA_TAG: Tag = Tag::new(b'S', b'A');

/// Collect SV evidence in a single pass. `contig_lengths` selects which contigs get
/// depth bins (and their sizes); `expected_insert_size`/`insert_size_sd` come from
/// read-metrics. `reference` is required for CRAM (ignored for BAM) — SV evidence never consults
/// reference *bases*, but decoding a CRAM record at all does. Mirrors the Scala walker's
/// thresholds and filters.
pub fn collect_evidence(
    bam_path: &Path,
    reference: Option<&Path>,
    contig_lengths: &BTreeMap<String, i64>,
    expected_insert_size: f64,
    insert_size_sd: f64,
    config: &SvCallerConfig,
) -> Result<SvEvidenceCollection, AnalysisError> {
    let (header, mut reader) = crate::reader::open_seq(bam_path, reference)?;

    // Reference id -> name (header order).
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .collect();

    // Depth bins per requested contig.
    let mut depth_bins: BTreeMap<String, Vec<u32>> = contig_lengths
        .iter()
        .map(|(c, &len)| {
            let num_bins = ((len + config.bin_size - 1) / config.bin_size).max(0) as usize;
            (c.clone(), vec![0u32; num_bins])
        })
        .collect();

    let insert_max = expected_insert_size + config.insert_size_z_threshold * insert_size_sd;
    let insert_min = (expected_insert_size - config.insert_size_z_threshold * insert_size_sd).max(0.0);

    let mut discordant_pairs = Vec::new();
    let mut split_reads = Vec::new();

    for result in reader.records_lazy(&header) {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() {
            continue;
        }
        let Some(ref_id) = record.reference_sequence_id() else {
            continue;
        };
        let Some(contig) = names.get(ref_id) else { continue };
        let Some(start) = record.alignment_start().map(|p| p as i64) else {
            continue;
        };
        let mapq = record.mapping_quality().unwrap_or(255);
        let secondary_or_supp = flags.is_secondary() || flags.is_supplementary();

        // 1. Depth (primary, non-supplementary only).
        if !secondary_or_supp {
            if let Some(bins) = depth_bins.get_mut(contig) {
                let bin = (start / config.bin_size) as usize;
                if bin < bins.len() {
                    bins[bin] += 1;
                }
            }
        }

        // 2. Discordant pairs (primary, paired only).
        if !secondary_or_supp && flags.is_segmented() {
            if let Some(dp) = detect_discordant_pair(&record, contig, &names, mapq, insert_min, insert_max, config) {
                discordant_pairs.push(dp);
            }
        }

        // 3. Split reads (SA tag).
        if mapq >= config.min_mapq {
            if let Some(sr) = extract_split_read(&record, contig, mapq, config) {
                split_reads.push(sr);
            }
        }
    }

    Ok(SvEvidenceCollection {
        discordant_pairs,
        split_reads,
        depth_bins,
        sample_name: "unknown".to_string(), // RG SM only feeds deferred VCF/summary output
        expected_insert_size,
        insert_size_sd,
    })
}

/// The record's name as a `String`, empty when unset or not UTF-8.
fn read_name(record: &impl AlnRead) -> String {
    record
        .name()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .unwrap_or_default()
}

#[allow(clippy::too_many_arguments)]
fn detect_discordant_pair(
    record: &impl AlnRead,
    contig: &str,
    names: &[String],
    mapq: u8,
    insert_min: f64,
    insert_max: f64,
    config: &SvCallerConfig,
) -> Option<DiscordantPair> {
    let flags = record.flags();
    if flags.is_mate_unmapped() || mapq < config.min_mapq {
        return None;
    }
    let read_strand = if flags.is_reverse_complemented() { '-' } else { '+' };
    let mate_strand = if flags.is_mate_reverse_complemented() { '-' } else { '+' };

    let ref_id = record.reference_sequence_id();
    let mate_ref_id = record.mate_reference_sequence_id();
    let mate_chrom = mate_ref_id.and_then(|i| names.get(i).cloned()).unwrap_or_default();
    let mate_pos = record.mate_alignment_start().map_or(0, |p| p as i64);
    let pos1 = record.alignment_start().map_or(0, |p| p as i64);
    let read_name = read_name(record);

    let mk = |insert_size: i32, reason: DiscordantReason| DiscordantPair {
        read_name: read_name.clone(),
        chrom1: contig.to_string(),
        pos1,
        strand1: read_strand,
        chrom2: mate_chrom.clone(),
        pos2: mate_pos,
        strand2: mate_strand,
        insert_size,
        mapq,
        reason,
    };

    // Inter-chromosomal takes precedence.
    if ref_id != mate_ref_id {
        return Some(mk(0, DiscordantReason::InterChromosomal));
    }
    let insert_size = record.template_length().abs();
    if insert_size as f64 > insert_max || (insert_size > 0 && (insert_size as f64) < insert_min) {
        return Some(mk(insert_size, DiscordantReason::InsertSizeOutlier));
    }
    if !is_expected_orientation(record, pos1, mate_pos) {
        return Some(mk(insert_size, DiscordantReason::WrongOrientation));
    }
    None
}

/// Standard Illumina FR orientation check (mirrors the Scala logic).
fn is_expected_orientation(record: &impl AlnRead, pos1: i64, mate_pos: i64) -> bool {
    let flags = record.flags();
    let read_neg = flags.is_reverse_complemented();
    let mate_neg = flags.is_mate_reverse_complemented();
    if read_neg == mate_neg {
        false // tandem
    } else if pos1 < mate_pos {
        !read_neg && mate_neg // read upstream -> read on +
    } else {
        read_neg && !mate_neg // read downstream -> read on -
    }
}

/// Parse the first SA-tag alignment into a [`SplitRead`]; clip length is the read's own
/// soft/hard-clip total.
fn extract_split_read(
    record: &impl AlnRead,
    contig: &str,
    mapq: u8,
    config: &SvCallerConfig,
) -> Option<SplitRead> {
    let sa = record.string_tag(SA_TAG)?;
    if sa.is_empty() {
        return None;
    }
    let first = sa.split(';').next().unwrap_or("");
    let parts: Vec<&str> = first.split(',').collect();
    if parts.len() < 5 {
        return None;
    }
    let (Ok(supp_pos), Some(supp_strand), Ok(supp_mapq)) =
        (parts[1].parse::<i64>(), parts[2].chars().next(), parts[4].parse::<u8>())
    else {
        return None;
    };

    // Clip length from this read's CIGAR (sum of S/H ops).
    let clip_length: i32 = record.cigar_with(|ops| {
        ops.filter(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
            .map(|(_, len)| len as i32)
            .sum()
    });

    if supp_mapq >= config.min_mapq && clip_length >= 10 {
        Some(SplitRead {
            read_name: read_name(record),
            primary_chrom: contig.to_string(),
            primary_pos: record.alignment_start().map_or(0, |p| p as i64),
            primary_strand: if record.flags().is_reverse_complemented() {
                '-'
            } else {
                '+'
            },
            supp_chrom: parts[0].to_string(),
            supp_pos,
            supp_strand,
            clip_length,
            mapq: mapq.min(supp_mapq),
        })
    } else {
        None
    }
}

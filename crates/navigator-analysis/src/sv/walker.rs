//! SV evidence walker — port of the Scala `SvEvidenceWalker`. Single pass over the BAM
//! collecting per-bin read depth (CNV), discordant read pairs (BreakDancer-style), and
//! split reads from the SA tag (Pindel-style).

use std::collections::BTreeMap;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::value::Value;
use noodles::sam::alignment::record::data::field::Tag;

use super::evidence::{DiscordantPair, DiscordantReason, SplitRead, SvEvidenceCollection};
use super::types::SvCallerConfig;
use crate::error::AnalysisError;

const SA_TAG: Tag = Tag::new(b'S', b'A');

/// Collect SV evidence in a single pass. `contig_lengths` selects which contigs get
/// depth bins (and their sizes); `expected_insert_size`/`insert_size_sd` come from
/// read-metrics. Mirrors the Scala walker's thresholds and filters.
pub fn collect_evidence(
    bam_path: &Path,
    contig_lengths: &BTreeMap<String, i64>,
    expected_insert_size: f64,
    insert_size_sd: f64,
    config: &SvCallerConfig,
) -> Result<SvEvidenceCollection, AnalysisError> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .map_err(|e| AnalysisError::io(bam_path, e))?;
    let header = reader
        .read_header()
        .map_err(|e| AnalysisError::io(bam_path, e))?;

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

    for result in reader.records() {
        let record = result.map_err(|e| AnalysisError::io(bam_path, e))?;
        let flags = record.flags();
        if flags.is_unmapped() {
            continue;
        }
        let Some(ref_id) = opt_usize(record.reference_sequence_id(), bam_path)? else {
            continue;
        };
        let Some(contig) = names.get(ref_id) else { continue };
        let start = match record.alignment_start() {
            Some(p) => p.map_err(|e| AnalysisError::io(bam_path, e))?.get() as i64,
            None => continue,
        };
        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
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
            if let Some(dp) = detect_discordant_pair(
                &record, contig, &names, mapq, insert_min, insert_max, config, bam_path,
            )? {
                discordant_pairs.push(dp);
            }
        }

        // 3. Split reads (SA tag).
        if mapq >= config.min_mapq {
            if let Some(sr) = extract_split_read(&record, contig, mapq, config, bam_path)? {
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

#[allow(clippy::too_many_arguments)]
fn detect_discordant_pair(
    record: &bam::Record,
    contig: &str,
    names: &[String],
    mapq: u8,
    insert_min: f64,
    insert_max: f64,
    config: &SvCallerConfig,
    path: &Path,
) -> Result<Option<DiscordantPair>, AnalysisError> {
    let flags = record.flags();
    if flags.is_mate_unmapped() || mapq < config.min_mapq {
        return Ok(None);
    }
    let read_strand = if flags.is_reverse_complemented() { '-' } else { '+' };
    let mate_strand = if flags.is_mate_reverse_complemented() { '-' } else { '+' };

    let ref_id = opt_usize(record.reference_sequence_id(), path)?;
    let mate_ref_id = opt_usize(record.mate_reference_sequence_id(), path)?;
    let mate_chrom = mate_ref_id
        .and_then(|i| names.get(i).cloned())
        .unwrap_or_default();
    let mate_pos = record
        .mate_alignment_start()
        .transpose()
        .map_err(|e| AnalysisError::io(path, e))?
        .map_or(0, |p| p.get() as i64);
    let pos1 = record
        .alignment_start()
        .transpose()
        .map_err(|e| AnalysisError::io(path, e))?
        .map_or(0, |p| p.get() as i64);
    let read_name = record.name().map(|n| n.to_string()).unwrap_or_default();

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
        return Ok(Some(mk(0, DiscordantReason::InterChromosomal)));
    }
    let insert_size = record.template_length().abs();
    if insert_size as f64 > insert_max || (insert_size > 0 && (insert_size as f64) < insert_min) {
        return Ok(Some(mk(insert_size, DiscordantReason::InsertSizeOutlier)));
    }
    if !is_expected_orientation(record, pos1, mate_pos) {
        return Ok(Some(mk(insert_size, DiscordantReason::WrongOrientation)));
    }
    Ok(None)
}

/// Standard Illumina FR orientation check (mirrors the Scala logic).
fn is_expected_orientation(record: &bam::Record, pos1: i64, mate_pos: i64) -> bool {
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
    record: &bam::Record,
    contig: &str,
    mapq: u8,
    config: &SvCallerConfig,
    path: &Path,
) -> Result<Option<SplitRead>, AnalysisError> {
    let data = record.data();
    let sa = match data.get(&SA_TAG) {
        Some(r) => match r.map_err(|e| AnalysisError::io(path, e))? {
            Value::String(s) => s.to_string(),
            _ => return Ok(None),
        },
        None => return Ok(None),
    };
    if sa.is_empty() {
        return Ok(None);
    }
    let first = sa.split(';').next().unwrap_or("");
    let parts: Vec<&str> = first.split(',').collect();
    if parts.len() < 5 {
        return Ok(None);
    }
    let (Ok(supp_pos), Some(supp_strand), Ok(supp_mapq)) = (
        parts[1].parse::<i64>(),
        parts[2].chars().next(),
        parts[4].parse::<u8>(),
    ) else {
        return Ok(None);
    };

    // Clip length from this read's CIGAR (sum of S/H ops).
    let mut clip_length: i32 = 0;
    for op in record.cigar().iter() {
        let op = op.map_err(|e| AnalysisError::io(path, e))?;
        if matches!(op.kind(), Kind::SoftClip | Kind::HardClip) {
            clip_length += op.len() as i32;
        }
    }

    if supp_mapq >= config.min_mapq && clip_length >= 10 {
        let pos1 = record
            .alignment_start()
            .transpose()
            .map_err(|e| AnalysisError::io(path, e))?
            .map_or(0, |p| p.get() as i64);
        Ok(Some(SplitRead {
            read_name: record.name().map(|n| n.to_string()).unwrap_or_default(),
            primary_chrom: contig.to_string(),
            primary_pos: pos1,
            primary_strand: if record.flags().is_reverse_complemented() { '-' } else { '+' },
            supp_chrom: parts[0].to_string(),
            supp_pos,
            supp_strand,
            clip_length,
            mapq: mapq.min(supp_mapq),
        }))
    } else {
        Ok(None)
    }
}

fn opt_usize(
    v: Option<std::io::Result<usize>>,
    path: &Path,
) -> Result<Option<usize>, AnalysisError> {
    match v {
        Some(r) => Ok(Some(r.map_err(|e| AnalysisError::io(path, e))?)),
        None => Ok(None),
    }
}

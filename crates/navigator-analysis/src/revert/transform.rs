//! One alignment record → one original read, or a reason it was dropped.
//!
//! Everything an aligner did to a read has to be undone here, because the mapper downstream is
//! going to redo it against a different reference. The subtle one is orientation: SAM stores SEQ
//! and QUAL in *reference* orientation, so a read that mapped to the reverse strand is held
//! reverse-complemented relative to what the sequencer produced. Emitting that to FASTQ without
//! flipping it back would hand the mapper a read that never existed.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::RecordBuf;

use super::{HardClipPolicy, RevertParams, RevertStats};

/// `OQ` — the original base qualities, before any recalibration overwrote `QUAL`.
const OQ_TAG: Tag = Tag::new(b'O', b'Q');

/// Phred used when a record carries no qualities at all. See [`RevertStats::qualities_synthesized`]
/// — the reads are worth keeping, but the qualities are invented and must be counted as such.
const SYNTHETIC_PHRED: u8 = 40;

/// Which end of a template a read is, from the segment flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Mate {
    /// Not part of a pair (`0x1` clear), or paired but with neither/both segment bits set — a
    /// record whose flags contradict themselves can not be placed in an R1/R2 file.
    Unpaired,
    /// First segment (`0x40`).
    One,
    /// Last segment (`0x80`).
    Two,
}

/// Why a record produced no read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Skipped {
    /// Secondary alignment (`0x100`) — the read's full sequence lives on its primary record.
    Secondary,
    /// Supplementary alignment (`0x800`) — likewise, and usually hard-clipped besides.
    Supplementary,
    /// `SEQ` was `*`; there is no read here to recover.
    NoSequence,
    /// A hard-clipped primary under [`HardClipPolicy::Skip`].
    HardClipped,
    /// No read name, so the record can never be paired with anything.
    NoName,
}

/// An original, unaligned read: exactly what the sequencer emitted, as far as the record allows.
///
/// Ordered by `(name, mate)` so a sort brings mates together with R1 ahead of R2 — which is the
/// entire mechanism [`super::collate`] relies on.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct RevertedRead {
    pub name: Vec<u8>,
    pub mate: Mate,
    /// ASCII bases in sequencer orientation.
    pub sequence: Vec<u8>,
    /// Raw phred scores (no `+33` offset), same orientation and length as `sequence`.
    pub qualities: Vec<u8>,
}

impl RevertedRead {
    /// Rough heap footprint, for the sort's memory budget. The two `Vec`s dominate; the constant
    /// covers the struct itself and allocator overhead closely enough to size a buffer by.
    pub fn heap_bytes(&self) -> usize {
        // The read sits inline in the collator's `Vec`, so its own size counts; the three vectors
        // each cost an allocation on top of their contents. The flat 64 this replaced was smaller
        // than the struct alone, which was harmless against a fixed 256 MB budget and is not
        // against a budget sized from the machine (see `navigator_resource::spill_budget`).
        std::mem::size_of::<Self>()
            + self.name.len()
            + self.sequence.len()
            + self.qualities.len()
            + 3 * navigator_resource::ALLOCATION_OVERHEAD
    }
}

/// Record the reason a record was dropped.
pub fn count_skip(skipped: Skipped, stats: &mut RevertStats) {
    match skipped {
        Skipped::Secondary => stats.secondary_dropped += 1,
        Skipped::Supplementary => stats.supplementary_dropped += 1,
        Skipped::NoSequence | Skipped::NoName => stats.no_sequence_dropped += 1,
        // Already counted by `revert_record` before it consulted the policy — the stat means
        // "records affected by hard clipping" under both policies, so counting it again here
        // would double it on the skip path only.
        Skipped::HardClipped => {}
    }
}

/// Strip a record back to the read it was made from.
///
/// Alignment state (position, CIGAR, MAPQ, mate fields, aligner tags like `NM`/`MD`/`AS`) needs no
/// explicit clearing: none of it is carried into [`RevertedRead`], so it is dropped by
/// construction rather than by a list of fields someone has to remember to extend.
pub fn revert_record(
    record: &RecordBuf,
    params: &RevertParams,
    stats: &mut RevertStats,
) -> Result<RevertedRead, Skipped> {
    let flags = record.flags();

    // Primaries only. Both of these carry a partial view of a read whose whole sequence is on
    // another record, so keeping them would duplicate reads and truncate some of them.
    if flags.is_secondary() {
        return Err(Skipped::Secondary);
    }
    if flags.is_supplementary() {
        return Err(Skipped::Supplementary);
    }

    let name = record.name().ok_or(Skipped::NoName)?;
    let sequence = record.sequence().as_ref();
    if sequence.is_empty() {
        return Err(Skipped::NoSequence);
    }

    // A hard clip on a *primary* means sequence was thrown away before we ever saw the record —
    // rare, but real, and undetectable downstream once it reaches a FASTQ.
    let hard_clipped = record.cigar().as_ref().iter().any(|op| op.kind() == Kind::HardClip);
    if hard_clipped {
        stats.hard_clipped += 1;
        if params.hard_clipped == HardClipPolicy::Skip {
            // Counted here rather than through `count_skip`, so the stat means "records affected
            // by hard clipping" under either policy instead of changing definition with the flag.
            return Err(Skipped::HardClipped);
        }
    }

    let mut qualities =
        original_qualities(record, params, stats).unwrap_or_else(|| record.quality_scores().as_ref().to_vec());

    // QUAL may legitimately be `*`. The read still maps, so keep it and flag the invention.
    if qualities.is_empty() {
        qualities = vec![SYNTHETIC_PHRED; sequence.len()];
        stats.qualities_synthesized += 1;
    }

    // A record whose QUAL disagrees with SEQ is malformed; trust SEQ (the mapper needs one
    // quality per base) and pad or truncate rather than emitting an unwritable FASTQ record.
    if qualities.len() != sequence.len() {
        qualities.resize(sequence.len(), SYNTHETIC_PHRED);
    }

    let mut sequence = sequence.to_vec();
    if flags.is_reverse_complemented() {
        reverse_complement(&mut sequence);
        qualities.reverse();
    }

    if flags.is_unmapped() {
        stats.unmapped_reads += 1;
    }

    Ok(RevertedRead {
        name: name.to_vec(),
        mate: mate_of(flags),
        sequence,
        qualities,
    })
}

/// The `OQ` tag's qualities, if present and wanted. `OQ` is stored as an ASCII phred+33 string —
/// the same encoding FASTQ uses — so it is decoded back to raw phred here to match `QUAL`.
fn original_qualities(record: &RecordBuf, params: &RevertParams, stats: &mut RevertStats) -> Option<Vec<u8>> {
    use noodles::sam::alignment::record_buf::data::field::Value;

    if !params.prefer_original_qualities {
        return None;
    }
    let Value::String(oq) = record.data().get(&OQ_TAG)? else {
        return None;
    };
    if oq.is_empty() {
        return None;
    }
    stats.original_qualities_used += 1;
    Some(oq.iter().map(|q| q.saturating_sub(33)).collect())
}

fn mate_of(flags: noodles::sam::alignment::record::Flags) -> Mate {
    if !flags.is_segmented() {
        return Mate::Unpaired;
    }
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, false) => Mate::One,
        (false, true) => Mate::Two,
        // Both or neither: the record claims to be paired but will not say which end. It can not go in
        // a synchronized R1/R2 file, so it becomes a singleton rather than corrupting the pairing.
        _ => Mate::Unpaired,
    }
}

/// Reverse-complement in place. `N` and any other non-ACGT byte is passed through unchanged rather
/// than normalized — the read is meant to come out of here exactly as the sequencer produced it.
fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        };
    }
}

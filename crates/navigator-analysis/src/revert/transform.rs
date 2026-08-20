//! One alignment record becomes one original read, or a reason that the code dropped it.
//!
//! This code must undo everything that an aligner did to a read, because the mapper after it does
//! that work again, against a different reference.
//!
//! The orientation is the part to watch. SAM stores SEQ and QUAL in the orientation of the
//! *reference*. A read that mapped to the reverse strand then sits there as the reverse complement
//! of what the sequencer gave. To put that into a FASTQ, and not turn it back, would give the
//! mapper a read that never existed.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::RecordBuf;

use super::{HardClipPolicy, RevertParams, RevertStats};

/// The `OQ` tag. It holds the original base qualities, from before a recalibration wrote over
/// `QUAL`.
const OQ_TAG: Tag = Tag::new(b'O', b'Q');

/// The phred value that the code uses when a record carries no qualities at all. See
/// [`RevertStats::qualities_synthesized`]. The reads are worth a place in the output, but the code
/// invents those qualities, and the count must say so.
const SYNTHETIC_PHRED: u8 = 40;

/// Which end of a template a read is, from the segment flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Mate {
    /// The record is not part of a pair, where `0x1` is clear. Or it is part of a pair, and it has
    /// neither segment bit set, or both. The flags of such a record contradict themselves, so the
    /// code can not put it into an R1 or R2 file.
    Unpaired,
    /// First segment (`0x40`).
    One,
    /// Last segment (`0x80`).
    Two,
}

/// Why a record produced no read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Skipped {
    /// A secondary alignment, at flag `0x100`. The full sequence of the read lives on its primary
    /// record.
    Secondary,
    /// A supplementary alignment, at flag `0x800`. The same holds, and it usually carries a hard
    /// clip as well.
    Supplementary,
    /// `SEQ` was `*`; there is no read here to recover.
    NoSequence,
    /// A hard-clipped primary under [`HardClipPolicy::Skip`].
    HardClipped,
    /// The record has no read name, so nothing can ever pair with it.
    NoName,
}

/// An original, unaligned read. It is exactly what the sequencer gave, as far as the record lets
/// the code recover.
///
/// The order is by `(name, mate)`. A sort then brings two mates together, with the R1 before the
/// R2. That is the whole mechanism that [`super::collate`] stands on.
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
        // The read sits inside the `Vec` of the collator, so its own size counts. The three
        // vectors each cost an allocation on top of what they hold.
        //
        // The flat 64 before this was smaller than the struct alone. Against a fixed budget of
        // 256 MB that did no harm. Against a budget that comes from the machine it does. See
        // `navigator_resource::spill_budget`.
        std::mem::size_of::<Self>()
            + self.name.len()
            + self.sequence.len()
            + self.qualities.len()
            + 3 * navigator_resource::ALLOCATION_OVERHEAD
    }
}

/// Note down the reason that the code dropped a record.
pub fn count_skip(skipped: Skipped, stats: &mut RevertStats) {
    match skipped {
        Skipped::Secondary => stats.secondary_dropped += 1,
        Skipped::Supplementary => stats.supplementary_dropped += 1,
        Skipped::NoSequence | Skipped::NoName => stats.no_sequence_dropped += 1,
        // `revert_record` already counted this, before it read the policy. The statistic means
        // "the records that a hard clip affected", under either policy. To count it again here
        // would make it twice too large, and only on the skip path.
        Skipped::HardClipped => {}
    }
}

/// Strip a record back to the read that made it.
///
/// The alignment state needs no explicit removal. That state is the position, the CIGAR, the MAPQ,
/// the mate fields, and an aligner tag such as `NM`, `MD` or `AS`. None of it goes into a
/// [`RevertedRead`], so the construction itself drops all of it. There is no list of fields that
/// somebody must remember to extend.
pub fn revert_record(
    record: &RecordBuf,
    params: &RevertParams,
    stats: &mut RevertStats,
) -> Result<RevertedRead, Skipped> {
    let flags = record.flags();

    // Primary records alone. The other two kinds each hold part of a read whose whole sequence
    // sits on another record. To keep them would put some reads in twice, and cut others short.
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

    // A hard clip on a *primary* record means that somebody threw sequence away before this code
    // saw the record. That is rare, and it is real. Once such a read reaches a FASTQ, no later
    // step can find it.
    let hard_clipped = record.cigar().as_ref().iter().any(|op| op.kind() == Kind::HardClip);
    if hard_clipped {
        stats.hard_clipped += 1;
        if params.hard_clipped == HardClipPolicy::Skip {
            // The count goes here, and not through `count_skip`. The statistic then says "the
            // records that a hard clip affected" under either policy, and the flag does not move
            // what it says.
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

    // A record whose QUAL does not agree with SEQ is not correct. Trust SEQ, because the mapper
    // needs one quality at each base. Then add to the qualities, or cut them, and do not emit a
    // FASTQ record that no writer can write.
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

/// The qualities of the `OQ` tag, when that tag is present and the caller wants it. `OQ` holds an
/// ASCII phred+33 string, which is the same encoding that FASTQ uses. So the code decodes it back
/// to a raw phred here, to match `QUAL`.
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
        // Both bits, or neither. The record says that it is part of a pair, and it does not say
        // which end. It can not go into an R1 or R2 file that must stay in step. So it becomes a
        // singleton, and it does not put the wrong reads together.
        _ => Mate::Unpaired,
    }
}

/// Take the reverse complement in place. An `N`, and any other byte that is not ACGT, goes through
/// unchanged. The code does not normalize it, because the read must come out of here exactly as
/// the sequencer gave it.
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

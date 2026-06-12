//! `AlnRead` — the minimal record view the quality-metrics walkers need, abstracted over the
//! **lazy** `bam::Record` (zero-copy, the hot path) and the owned `RecordBuf` (CRAM). The walkers
//! used to consume `RecordBuf`, which forced a per-read owned copy of the sequence, qualities,
//! CIGAR, name, and *every* optional tag — measured at ~half the per-read CPU on a WGS BAM. The
//! walkers touch only a handful of scalar fields plus the CIGAR ops and per-base qualities, all of
//! which `bam::Record` exposes as borrowed views, so this trait lets the same accumulator code run
//! over either record type with no allocation on the BAM path.
//!
//! Implementations map noodles' lazy `io::Result` accessors to plain `Option`/values, treating a
//! decode error as "absent" (skips that field) rather than aborting the walk — strictly more robust
//! than the old `RecordBuf` conversion, which would have errored the whole pass on a bad record.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::RecordBuf;

/// The fields the coverage / read-metrics / sex walkers read from an alignment record.
pub trait AlnRead {
    fn flags(&self) -> Flags;
    /// Reference-sequence id (`@SQ` index), or `None` if unmapped/unset/undecodable.
    fn reference_sequence_id(&self) -> Option<usize>;
    fn mate_reference_sequence_id(&self) -> Option<usize>;
    /// 1-based alignment start, or `None`.
    fn alignment_start(&self) -> Option<usize>;
    fn mate_alignment_start(&self) -> Option<usize>;
    /// Mapping quality (`None` == 255/unavailable).
    fn mapping_quality(&self) -> Option<u8>;
    fn template_length(&self) -> i32;
    fn sequence_len(&self) -> usize;
    /// Run `f` with the per-base phred qualities (raw, no +33; indexable by query offset) and an
    /// iterator of CIGAR `(kind, len)` ops. Via a callback so the lazy `bam::Record` views (which
    /// borrow the record's buffer through a temporary wrapper) stay alive for the duration — no
    /// per-read allocation. Undecodable CIGAR ops are skipped.
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T;
}

impl AlnRead for RecordBuf {
    fn flags(&self) -> Flags {
        RecordBuf::flags(self)
    }
    fn reference_sequence_id(&self) -> Option<usize> {
        RecordBuf::reference_sequence_id(self)
    }
    fn mate_reference_sequence_id(&self) -> Option<usize> {
        RecordBuf::mate_reference_sequence_id(self)
    }
    fn alignment_start(&self) -> Option<usize> {
        RecordBuf::alignment_start(self).map(|p| p.get())
    }
    fn mate_alignment_start(&self) -> Option<usize> {
        RecordBuf::mate_alignment_start(self).map(|p| p.get())
    }
    fn mapping_quality(&self) -> Option<u8> {
        RecordBuf::mapping_quality(self).map(|m| m.get())
    }
    fn template_length(&self) -> i32 {
        RecordBuf::template_length(self)
    }
    fn sequence_len(&self) -> usize {
        self.sequence().len()
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        use noodles::sam::alignment::record::Cigar as _; // RecordBuf's Cigar iterates via the trait
        let quals = self.quality_scores();
        let cigar = self.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(quals.as_ref(), &mut ops)
    }
}

impl AlnRead for noodles::bam::Record {
    fn flags(&self) -> Flags {
        noodles::bam::Record::flags(self)
    }
    fn reference_sequence_id(&self) -> Option<usize> {
        self.reference_sequence_id().and_then(|r| r.ok())
    }
    fn mate_reference_sequence_id(&self) -> Option<usize> {
        self.mate_reference_sequence_id().and_then(|r| r.ok())
    }
    fn alignment_start(&self) -> Option<usize> {
        self.alignment_start().and_then(|r| r.ok()).map(|p| p.get())
    }
    fn mate_alignment_start(&self) -> Option<usize> {
        self.mate_alignment_start().and_then(|r| r.ok()).map(|p| p.get())
    }
    fn mapping_quality(&self) -> Option<u8> {
        noodles::bam::Record::mapping_quality(self).map(|m| m.get())
    }
    fn template_length(&self) -> i32 {
        noodles::bam::Record::template_length(self)
    }
    fn sequence_len(&self) -> usize {
        self.sequence().len()
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        let quals = self.quality_scores();
        let cigar = self.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(quals.as_ref(), &mut ops)
    }
}

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
use noodles::sam::alignment::record::{Cigar as _, Flags, QualityScores as _, Sequence as _};
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

/// A borrowed view over a decoded **CRAM** record (`noodles::cram::Record`) paired with the header,
/// implementing [`AlnRead`] by delegating to the `sam::alignment::Record` trait. CRAM stores the
/// sequence as deltas against the reference, so a `cram::Record` already holds the per-read data in
/// borrowed/lightweight form — driving the walkers off it directly skips the per-read
/// `RecordBuf::try_from_alignment_record` copy (sequence + quals + cigar + name + *every* tag into
/// owned form), measured at ~1.74× the per-read decode cost on a 30× short-read WGS CRAM. The
/// header is only needed by the trait's `reference_sequence_id` accessors (which, for CRAM, ignore
/// it and return the record's stored id — but the signature requires one).
pub struct CramRead<'a, 'c> {
    pub rec: &'a noodles::cram::Record<'c>,
    pub header: &'a noodles::sam::Header,
}

impl AlnRead for CramRead<'_, '_> {
    fn flags(&self) -> Flags {
        use noodles::sam::alignment::Record as _;
        self.rec.flags().unwrap_or(Flags::UNMAPPED)
    }
    fn reference_sequence_id(&self) -> Option<usize> {
        use noodles::sam::alignment::Record as _;
        self.rec.reference_sequence_id(self.header).and_then(|r| r.ok())
    }
    fn mate_reference_sequence_id(&self) -> Option<usize> {
        use noodles::sam::alignment::Record as _;
        self.rec.mate_reference_sequence_id(self.header).and_then(|r| r.ok())
    }
    fn alignment_start(&self) -> Option<usize> {
        use noodles::sam::alignment::Record as _;
        self.rec.alignment_start().and_then(|r| r.ok()).map(|p| p.get())
    }
    fn mate_alignment_start(&self) -> Option<usize> {
        use noodles::sam::alignment::Record as _;
        self.rec.mate_alignment_start().and_then(|r| r.ok()).map(|p| p.get())
    }
    fn mapping_quality(&self) -> Option<u8> {
        use noodles::sam::alignment::Record as _;
        self.rec.mapping_quality().and_then(|r| r.ok()).map(|m| m.get())
    }
    fn template_length(&self) -> i32 {
        use noodles::sam::alignment::Record as _;
        self.rec.template_length().unwrap_or(0)
    }
    fn sequence_len(&self) -> usize {
        use noodles::sam::alignment::Record as _;
        self.rec.sequence().len()
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        use noodles::sam::alignment::Record as _;
        // CRAM exposes qualities only via an iterator (not a contiguous slice), so collect them
        // once per read — a small (~read-length) allocation, still far cheaper than the full
        // `RecordBuf` materialization the high-level reader would do. The cigar is a lazy view over
        // the record's features, iterated directly with no allocation.
        let quals: Vec<u8> = self.rec.quality_scores().iter().map(|r| r.unwrap_or(0)).collect();
        let cigar = self.rec.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(&quals, &mut ops)
    }
}

/// A record yielded by a **sequential** (whole-file, no index) walk over either format: the
/// **lazy, zero-copy** `bam::Record` on the BAM path (no owned `RecordBuf` decode/tag-parse — the
/// hot-path win) and the decoded `RecordBuf` on the CRAM path (CRAM has no cheaper lazy form). It
/// implements [`AlnRead`] by delegating to the per-type impls above, so the same accumulator code
/// (`CoverageState`/`ReadMetricsState`/`SexState`) drives both with no allocation on the BAM path —
/// the sequential counterpart to the indexed [`crate::reader::RecordSink`] fan-out.
pub enum SeqRecord {
    Bam(noodles::bam::Record),
    Cram(RecordBuf),
}

impl AlnRead for SeqRecord {
    fn flags(&self) -> Flags {
        match self {
            SeqRecord::Bam(r) => AlnRead::flags(r),
            SeqRecord::Cram(r) => AlnRead::flags(r),
        }
    }
    fn reference_sequence_id(&self) -> Option<usize> {
        match self {
            SeqRecord::Bam(r) => AlnRead::reference_sequence_id(r),
            SeqRecord::Cram(r) => AlnRead::reference_sequence_id(r),
        }
    }
    fn mate_reference_sequence_id(&self) -> Option<usize> {
        match self {
            SeqRecord::Bam(r) => AlnRead::mate_reference_sequence_id(r),
            SeqRecord::Cram(r) => AlnRead::mate_reference_sequence_id(r),
        }
    }
    fn alignment_start(&self) -> Option<usize> {
        match self {
            SeqRecord::Bam(r) => AlnRead::alignment_start(r),
            SeqRecord::Cram(r) => AlnRead::alignment_start(r),
        }
    }
    fn mate_alignment_start(&self) -> Option<usize> {
        match self {
            SeqRecord::Bam(r) => AlnRead::mate_alignment_start(r),
            SeqRecord::Cram(r) => AlnRead::mate_alignment_start(r),
        }
    }
    fn mapping_quality(&self) -> Option<u8> {
        match self {
            SeqRecord::Bam(r) => AlnRead::mapping_quality(r),
            SeqRecord::Cram(r) => AlnRead::mapping_quality(r),
        }
    }
    fn template_length(&self) -> i32 {
        match self {
            SeqRecord::Bam(r) => AlnRead::template_length(r),
            SeqRecord::Cram(r) => AlnRead::template_length(r),
        }
    }
    fn sequence_len(&self) -> usize {
        match self {
            SeqRecord::Bam(r) => AlnRead::sequence_len(r),
            SeqRecord::Cram(r) => AlnRead::sequence_len(r),
        }
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        match self {
            SeqRecord::Bam(r) => AlnRead::pileup_with(r, f),
            SeqRecord::Cram(r) => AlnRead::pileup_with(r, f),
        }
    }
}

//! `AlnRead` is the smallest view of a record that the quality-metrics walkers need. It covers two
//! types: the **lazy** `bam::Record`, which is zero-copy and the hot path, and the owned
//! `RecordBuf`, which a CRAM gives.
//!
//! The walkers once took a `RecordBuf`. That forced an owned copy, at each read, of the sequence,
//! the qualities, the CIGAR, the name, and *every* optional tag. A measurement on a WGS BAM put
//! that at about half of the CPU at each read.
//!
//! The walkers touch a few scalar fields, plus the CIGAR operations and the quality of each base.
//! `bam::Record` gives all of those as borrowed views. So this trait lets the same accumulator
//! code run over either record type, and it allocates nothing on the BAM path.
//!
//! An implementation maps the lazy `io::Result` accessors of noodles to a plain `Option` or a
//! value. It reads a decode error as "absent", and it skips that field. It does not stop the walk.
//! That is more robust than the old conversion to a `RecordBuf`, which would have failed the whole
//! pass on one bad record.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Cigar as _, Data as _, Flags, QualityScores as _, Sequence as _};
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
    /// The name of the read, which is its template name, as raw bytes. It is `None` when the
    /// record has none. This borrows the bytes, so the caller decides whether to pay for a
    /// `String`.
    fn name(&self) -> Option<&[u8]>;
    /// An auxiliary tag whose value is a string, such as `SA`. It is `None` when the tag is
    /// absent, when the code can not decode it, or when it holds another type.
    ///
    /// This owns the value, because the three record types each give their `Data` view a
    /// different shape. The one tag that this serves is `SA`, and it occurs on a small minority of
    /// the reads. So the allocation does not sit on the hot path.
    fn string_tag(&self, tag: Tag) -> Option<String>;
    /// Run `f` with an iterator over the CIGAR operations, as `(kind, len)`. The callback form
    /// keeps the borrowed view of the lazy `bam::Record` alive while `f` runs. It skips an
    /// operation that the code can not decode.
    ///
    /// Use this, and not [`AlnRead::pileup_with`], when you need the CIGAR alone. The CRAM version
    /// of `pileup_with` builds the quality scores, and this function does not.
    fn cigar_with<T>(&self, f: impl FnOnce(&mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T;
    /// Run `f` with two things: the phred quality of each base, and an iterator over the CIGAR
    /// operations as `(kind, len)`. The qualities are raw, with no +33, and you index them by the
    /// offset into the query.
    ///
    /// It uses a callback so that the views of the lazy `bam::Record` stay alive while `f` runs.
    /// Those views borrow the buffer of the record through a temporary wrapper. There is then no
    /// allocation at each read. It skips a CIGAR operation that the code can not decode.
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
    fn name(&self) -> Option<&[u8]> {
        RecordBuf::name(self).map(|n| &**n)
    }
    fn string_tag(&self, tag: Tag) -> Option<String> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        match RecordBuf::data(self).get(&tag)? {
            Value::String(s) => Some(s.to_string()),
            _ => None,
        }
    }
    fn cigar_with<T>(&self, f: impl FnOnce(&mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        use noodles::sam::alignment::record::Cigar as _; // RecordBuf's Cigar iterates via the trait
        let cigar = self.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(&mut ops)
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
    fn name(&self) -> Option<&[u8]> {
        noodles::bam::Record::name(self).map(|n| &**n)
    }
    fn string_tag(&self, tag: Tag) -> Option<String> {
        use noodles::sam::alignment::record::data::field::Value;
        match self.data().get(&tag)?.ok()? {
            Value::String(s) => Some(s.to_string()),
            _ => None,
        }
    }
    fn cigar_with<T>(&self, f: impl FnOnce(&mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        let cigar = self.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(&mut ops)
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        let quals = self.quality_scores();
        let cigar = self.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(quals.as_ref(), &mut ops)
    }
}

/// A borrowed view over a decoded **CRAM** record, which is a `noodles::cram::Record`, together
/// with the header. It implements [`AlnRead`], and it hands the work to the
/// `sam::alignment::Record` trait.
///
/// A CRAM stores the sequence as deltas against the reference. So a `cram::Record` already holds
/// the data of each read in a borrowed, light form. To drive the walkers from it directly
/// leaves out the `RecordBuf::try_from_alignment_record` copy at each read.
///
/// That copy takes the sequence, the qualities, the CIGAR, the name and *every* tag into owned
/// form. A measurement on a 30x short-read WGS CRAM put it at about 1.74 times the cost of the
/// decode of one read.
///
/// The header serves the `reference_sequence_id` accessors of the trait alone. For a CRAM those
/// ignore it, and they return the id that the record stores. But the signature needs one.
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
    fn name(&self) -> Option<&[u8]> {
        use noodles::sam::alignment::Record as _;
        self.rec.name().map(|n| &**n)
    }
    fn string_tag(&self, tag: Tag) -> Option<String> {
        use noodles::sam::alignment::record::data::field::Value;
        use noodles::sam::alignment::Record as _;
        match self.rec.data().get(&tag)?.ok()? {
            Value::String(s) => Some(s.to_string()),
            _ => None,
        }
    }
    fn cigar_with<T>(&self, f: impl FnOnce(&mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        use noodles::sam::alignment::Record as _;
        let cigar = self.rec.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(&mut ops)
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        use noodles::sam::alignment::Record as _;
        // A CRAM gives the qualities through an iterator alone, and not as one slice. So the code
        // collects them once at each read. That allocation is small, at about one read length, and
        // it still costs far less than the full `RecordBuf` that the high-level reader would
        // build. The cigar is a lazy view over the features of the record, and the code walks it
        // directly, with no allocation.
        let quals: Vec<u8> = self.rec.quality_scores().iter().map(|r| r.unwrap_or(0)).collect();
        let cigar = self.rec.cigar();
        let mut ops = cigar.iter().filter_map(|op| op.ok().map(|o| (o.kind(), o.len())));
        f(&quals, &mut ops)
    }
}

/// A record that a **sequential** walk gives, over either format. Such a walk covers the whole
/// file and uses no index.
///
/// On the BAM path it is the **lazy, zero-copy** `bam::Record`. There is then no decode into an
/// owned `RecordBuf`, and no parse of a tag, and that is the gain on the hot path. On the CRAM
/// path it is the decoded `RecordBuf`, because a CRAM has no cheaper lazy form.
///
/// It implements [`AlnRead`] by a call into the implementation of each type above. The same
/// accumulator code, which is `CoverageState`, `ReadMetricsState` and `SexState`, then drives
/// both, and it allocates nothing on the BAM path. This is the sequential counterpart of the
/// indexed [`crate::reader::RecordSink`] fan-out.
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
    fn name(&self) -> Option<&[u8]> {
        match self {
            SeqRecord::Bam(r) => AlnRead::name(r),
            SeqRecord::Cram(r) => AlnRead::name(r),
        }
    }
    fn string_tag(&self, tag: Tag) -> Option<String> {
        match self {
            SeqRecord::Bam(r) => AlnRead::string_tag(r, tag),
            SeqRecord::Cram(r) => AlnRead::string_tag(r, tag),
        }
    }
    fn cigar_with<T>(&self, f: impl FnOnce(&mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        match self {
            SeqRecord::Bam(r) => AlnRead::cigar_with(r, f),
            SeqRecord::Cram(r) => AlnRead::cigar_with(r, f),
        }
    }
    fn pileup_with<T>(&self, f: impl FnOnce(&[u8], &mut dyn Iterator<Item = (Kind, usize)>) -> T) -> T {
        match self {
            SeqRecord::Bam(r) => AlnRead::pileup_with(r, f),
            SeqRecord::Cram(r) => AlnRead::pileup_with(r, f),
        }
    }
}

//! Tests for the coordinate sort.
//!
//! The property that matters most is that the sort loses nothing. A sort that drops records, where
//! nobody sees it, would show later as a coverage that reads a little low. Almost nobody would
//! find that.
//!
//! The second property is the test that the spill path and the no-spill path agree. The result
//! must not depend on how much of the input fit in memory.

use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::{map, Map};
use noodles::{bam, sam};

use super::sort::*;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-sort-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn header(contigs: &[(&str, usize)]) -> sam::Header {
    let mut builder = sam::Header::builder();
    for (name, len) in contigs {
        builder = builder.add_reference_sequence(
            name.as_bytes(),
            Map::<map::ReferenceSequence>::new(std::num::NonZeroUsize::new(*len).unwrap()),
        );
    }
    builder.build()
}

/// A record at `(ref_id, pos)`, or unplaced when `ref_id` is `None`.
fn record(name: &str, ref_id: Option<usize>, pos: usize) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_sequence(Sequence::from(b"ACGTACGTAC".to_vec()))
        .set_quality_scores(QualityScores::from(vec![30; 10]));
    match ref_id {
        Some(id) => {
            builder = builder
                .set_flags(Flags::from(0u16))
                .set_reference_sequence_id(id)
                .set_alignment_start(noodles::core::Position::new(pos).unwrap());
        }
        None => builder = builder.set_flags(Flags::UNMAPPED),
    }
    builder.build()
}

fn write_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(header).unwrap();
    for record in records {
        writer.write_alignment_record(header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

/// `(name, ref_id, pos)` for every record, in file order.
fn read_back(path: &Path) -> (sam::Header, Vec<(String, Option<usize>, usize)>) {
    let file = std::fs::File::open(path).unwrap();
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.unwrap();
            (
                String::from_utf8_lossy(r.name().unwrap_or_default()).to_string(),
                r.reference_sequence_id(),
                r.alignment_start().map(|p| p.get()).unwrap_or(0),
            )
        })
        .collect();
    (header, records)
}

fn run(input: &Path, dir: &Path, params: SortParams) -> (PathBuf, SortStats) {
    let output = dir.join("sorted.bam");
    let stats = sort_alignment(
        input,
        &output,
        &dir.join("scratch"),
        &params,
        &CancelToken::none(),
        &mut |_| {},
    )
    .expect("sort should succeed");
    (output, stats)
}

/// Shuffled input across two contigs, plus unplaced reads.
fn unsorted_fixture(dir: &Path) -> (PathBuf, sam::Header, usize) {
    let hdr = header(&[("chr1", 100_000), ("chr2", 100_000)]);
    let mut records = Vec::new();
    for i in 0..60 {
        // Deliberately out of order, and interleaved between the two contigs.
        let pos = ((i * 4703) % 90_000) + 1;
        let ref_id = i % 2;
        records.push(record(&format!("r{i:03}"), Some(ref_id), pos));
    }
    for i in 0..5 {
        records.push(record(&format!("u{i}"), None, 0));
    }
    let input = dir.join("unsorted.bam");
    let total = records.len();
    write_bam(&input, &hdr, &records);
    (input, hdr, total)
}

// ---- the buffer estimate ---------------------------------------------------

/// The size of the buffer is now a fraction of the machine, and not a constant that somebody
/// chose. So the cost that the code charges a record against that buffer must be about correct. It
/// was not correct before: the tag dictionary cost nothing, and a mapped record carries a dozen
/// tags.
#[test]
fn the_buffer_estimate_counts_the_tag_dictionary() {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    let bare = record("r0", Some(0), 1);
    let mut tagged = bare.clone();
    for (tag, value) in [
        (Tag::ALIGNMENT_HIT_COUNT, Value::from(1i32)),
        (Tag::MISMATCHED_POSITIONS, Value::from("10")),
        (Tag::ALIGNMENT_SCORE, Value::from(60i32)),
    ] {
        tagged.data_mut().insert(tag, value);
    }

    let charged = heap_bytes(&tagged) - heap_bytes(&bare);
    assert_eq!(charged, 3 * TAG_ENTRY_BYTES, "three tags should cost three entries");
}

/// An upgrade of noodles that makes `Value` larger must fail here. Without this test, every buffer
/// would hold more than its budget says, and nobody would see it.
#[test]
fn a_tag_entry_is_not_larger_than_the_estimate_assumes() {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    assert!(
        std::mem::size_of::<(Tag, Value)>() <= TAG_ENTRY_BYTES,
        "a tag entry is {} bytes, which the estimate does not cover",
        std::mem::size_of::<(Tag, Value)>()
    );
}

/// The own size of the record is part of what it costs, because it lives inside the `Vec` of the
/// buffer.
#[test]
fn the_buffer_estimate_covers_the_record_itself() {
    let empty = RecordBuf::default();
    assert!(heap_bytes(&empty) >= std::mem::size_of::<RecordBuf>());
}

// ---- the properties that matter -------------------------------------------

/// Coordinate order, with unplaced reads at the end where SAM puts them.
#[test]
fn records_come_out_in_coordinate_order_with_unplaced_reads_last() {
    let dir = scratch("order");
    let (input, _, total) = unsorted_fixture(&dir);
    let (output, stats) = run(&input, &dir, SortParams::default());

    assert_eq!(stats.records as usize, total);
    assert_eq!(stats.unplaced, 5);

    let (_, records) = read_back(&output);
    assert_eq!(records.len(), total);

    let placed: Vec<_> = records.iter().filter(|r| r.1.is_some()).collect();
    let unplaced: Vec<_> = records.iter().filter(|r| r.1.is_none()).collect();
    assert_eq!(unplaced.len(), 5);
    assert_eq!(
        records[records.len() - 5..].iter().filter(|r| r.1.is_none()).count(),
        5,
        "the unplaced reads must be the final records, not scattered"
    );

    for pair in placed.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        assert!((a.1, a.2) <= (b.1, b.2), "out of order: {:?} before {:?}", a, b);
    }
}

/// A sort that loses records would show later as a coverage that is a little low, and almost
/// nobody would find that. The sort must drop nothing, and that includes the reads with no place,
/// which realignment exists to recover.
#[test]
fn the_sort_is_lossless() {
    let dir = scratch("lossless");
    let (input, _, total) = unsorted_fixture(&dir);
    let (output, stats) = run(&input, &dir, SortParams { buffer_bytes: 1 });

    let (_, records) = read_back(&output);
    assert_eq!(records.len(), total);
    assert_eq!(stats.records as usize, total);

    let mut names: Vec<&str> = records.iter().map(|r| r.0.as_str()).collect();
    names.sort_unstable();
    let mut expected: Vec<String> = (0..60).map(|i| format!("r{i:03}")).collect();
    expected.extend((0..5).map(|i| format!("u{i}")));
    expected.sort_unstable();
    assert_eq!(names, expected, "every input record must appear exactly once");
}

/// The result must not depend on how much of the input fit in memory. The external sort of the
/// revert stage makes the same claim. That is the reason you can rely on the path that uses the
/// disk.
#[test]
fn spilling_produces_the_same_order_as_sorting_in_memory() {
    let dir = scratch("spill");
    let (input, _, _) = unsorted_fixture(&dir);

    let in_memory_dir = dir.join("mem");
    std::fs::create_dir_all(&in_memory_dir).unwrap();
    let (a, a_stats) = run(
        &input,
        &in_memory_dir,
        SortParams {
            buffer_bytes: 64 * 1024 * 1024,
        },
    );

    let spilled_dir = dir.join("spill");
    std::fs::create_dir_all(&spilled_dir).unwrap();
    let (b, b_stats) = run(&input, &spilled_dir, SortParams { buffer_bytes: 1 });

    assert_eq!(a_stats.runs, 1, "everything fit");
    assert!(b_stats.runs > 1, "the budget forced multiple runs");
    assert_eq!(read_back(&a).1, read_back(&b).1);
}

/// An index is correct only for a file in coordinate order, and a reader reads this to decide
/// whether it may query a region. A sort that works, and that does not say so, makes every reader
/// scan the file again.
#[test]
fn the_output_header_declares_coordinate_order() {
    let dir = scratch("header");
    let (input, _, _) = unsorted_fixture(&dir);
    let (output, _) = run(&input, &dir, SortParams::default());

    let (header, _) = read_back(&output);
    let hd = header.header().expect("@HD must be present");
    use noodles::sam::header::record::value::map::header::tag;
    let so = hd.other_fields().get(&tag::SORT_ORDER).expect("@HD must carry SO");
    assert_eq!(so.as_slice(), b"coordinate");
}

/// Runs are the size of the alignment itself, so they must not outlive the sort.
#[test]
fn spilled_runs_are_cleaned_up() {
    let dir = scratch("cleanup");
    let (input, _, _) = unsorted_fixture(&dir);
    let _ = run(&input, &dir, SortParams { buffer_bytes: 1 });

    let leftovers: Vec<String> = std::fs::read_dir(dir.join("scratch"))
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("sort-run-"))
        .collect();
    assert!(leftovers.is_empty(), "run files left behind: {leftovers:?}");
}

/// Sorting a WGS is long enough to need a cancel that works, reported as itself.
#[test]
fn cancellation_stops_the_sort() {
    let dir = scratch("cancel");
    let (input, _, _) = unsorted_fixture(&dir);
    let cancel = CancelToken::new();
    cancel.cancel();

    let err = sort_alignment(
        &input,
        &dir.join("out.bam"),
        &dir.join("scratch"),
        &SortParams::default(),
        &cancel,
        &mut |_| {},
    )
    .unwrap_err();
    assert!(matches!(err, AnalysisError::Cancelled));
}

/// A file that is already sorted is a usual input. It comes from a second run, or from a mapper
/// that gave its output in order. It must come out unchanged, and the sort must not move a record
/// where nobody looks.
#[test]
fn an_already_sorted_file_is_unchanged() {
    let dir = scratch("idempotent");
    let (input, _, _) = unsorted_fixture(&dir);
    let (once, _) = run(&input, &dir, SortParams::default());

    let twice_dir = dir.join("again");
    std::fs::create_dir_all(&twice_dir).unwrap();
    let (twice, _) = run(&once, &twice_dir, SortParams::default());

    assert_eq!(read_back(&once).1, read_back(&twice).1);
}

// ---- the mark on a duplicate ----------------------------------------------

use super::markdup::{mark_duplicates, MarkDupParams};

/// A paired record. `pos` is the alignment start, `mate_pos` the mate's.
#[allow(clippy::too_many_arguments)]
fn pair_record(
    name: &str,
    ref_id: usize,
    pos: usize,
    reverse: bool,
    first: bool,
    mate_ref: usize,
    mate_pos: usize,
    mate_reverse: bool,
    cigar: &str,
) -> RecordBuf {
    let mut flags = 0x1u16 | if first { 0x40 } else { 0x80 };
    if reverse {
        flags |= 0x10;
    }
    if mate_reverse {
        flags |= 0x20;
    }
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::from(flags))
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(noodles::core::Position::new(pos).unwrap())
        .set_cigar(parse_cigar(cigar))
        .set_mate_reference_sequence_id(mate_ref)
        .set_mate_alignment_start(noodles::core::Position::new(mate_pos).unwrap())
        // SEQ must have exactly the length that the CIGAR takes from the query. A soft clip
        // counts, and a hard clip does not. noodles refuses the record if the two differ, and that
        // is how somebody found this.
        .set_sequence(Sequence::from(vec![b'A'; query_len(cigar)]))
        .set_quality_scores(QualityScores::from(vec![30; query_len(cigar)]))
        .build()
}

/// Bases of the read a CIGAR accounts for: `M`/`I`/`S`/`=`/`X` consume query, `D`/`N`/`H` do not.
fn query_len(spec: &str) -> usize {
    let mut total = 0;
    let mut digits = String::new();
    for c in spec.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            let n: usize = digits.parse().unwrap();
            digits.clear();
            if matches!(c, 'M' | 'I' | 'S' | '=' | 'X') {
                total += n;
            }
        }
    }
    total
}

/// A small CIGAR parser for the fixtures, for example "3S10M".
fn parse_cigar(spec: &str) -> noodles::sam::alignment::record_buf::Cigar {
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};
    let mut ops = Vec::new();
    let mut digits = String::new();
    for c in spec.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else {
            let n: usize = digits.parse().unwrap();
            digits.clear();
            let kind = match c {
                'M' => Kind::Match,
                'S' => Kind::SoftClip,
                'H' => Kind::HardClip,
                'D' => Kind::Deletion,
                'I' => Kind::Insertion,
                other => panic!("unhandled cigar op {other}"),
            };
            ops.push(Op::new(kind, n));
        }
    }
    noodles::sam::alignment::record_buf::Cigar::from(ops)
}

/// `(name, is_duplicate)` for every record, in file order.
fn duplicate_flags(path: &Path) -> Vec<(String, bool)> {
    let file = std::fs::File::open(path).unwrap();
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.unwrap();
            (
                String::from_utf8_lossy(r.name().unwrap_or_default()).to_string(),
                r.flags().is_duplicate(),
            )
        })
        .collect()
}

fn mark(input: &Path, dir: &Path, params: MarkDupParams) -> (PathBuf, super::MarkDupStats) {
    let output = dir.join("marked.bam");
    let stats =
        mark_duplicates(input, &output, &params, &CancelToken::none(), &mut |_| {}).expect("marking should succeed");
    (output, stats)
}

/// Two templates from the same molecule, with the same endpoints and the same strands. One of them
/// stays without a mark, and the code flags the other. A separate template at another position does
/// not change.
#[test]
fn identical_fragments_are_marked_and_one_representative_is_kept() {
    let dir = scratch("dupes");
    let hdr = header(&[("chr1", 100_000)]);
    let records = vec![
        pair_record("a", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("b", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("other", 0, 900, false, true, 0, 1500, true, "10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (output, stats) = mark(&input, &dir, MarkDupParams::default());
    let flags = duplicate_flags(&output);

    assert_eq!(stats.records, 3);
    assert_eq!(stats.duplicates, 1, "exactly one of the two copies is marked");
    assert_eq!(flags[0], ("a".into(), false), "the first seen is the representative");
    assert_eq!(flags[1], ("b".into(), true));
    assert_eq!(
        flags[2],
        ("other".into(), false),
        "an independent fragment is untouched"
    );
}

/// **The property the module docs argue for.** A template with one end marked and the other not
/// would show consumers half a pair. Both ends of a duplicate template must agree.
#[test]
fn both_ends_of_a_duplicate_template_are_marked_alike() {
    let dir = scratch("bothends");
    let hdr = header(&[("chr1", 100_000)]);
    // Two duplicate templates, each with both ends present, in coordinate order.
    let records = vec![
        pair_record("a", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("b", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("a", 0, 500, true, false, 0, 100, false, "10M"),
        pair_record("b", 0, 500, true, false, 0, 100, false, "10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (output, _) = mark(&input, &dir, MarkDupParams::default());
    let flags = duplicate_flags(&output);

    let verdicts = |name: &str| -> Vec<bool> { flags.iter().filter(|(n, _)| n == name).map(|(_, d)| *d).collect() };
    assert_eq!(verdicts("a"), vec![false, false], "both ends of 'a' agree");
    assert_eq!(verdicts("b"), vec![true, true], "both ends of 'b' agree");
}

/// Two copies of one molecule can carry different clips, and one mismatch near an end is enough to
/// cause that. A group on the alignment start would miss them. A group on the 5' position before
/// the clip finds them.
#[test]
fn differently_clipped_copies_of_one_fragment_are_still_duplicates() {
    let dir = scratch("clipping");
    let hdr = header(&[("chr1", 100_000)]);
    // Both molecules begin at 100. One aligns from 100 with no clip. The other carries a clip of
    // 3, so it *starts* at 103, and it covers the same fragment.
    let records = vec![
        pair_record("plain", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("clipped", 0, 103, false, true, 0, 500, true, "3S10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (output, stats) = mark(&input, &dir, MarkDupParams::default());
    assert_eq!(stats.duplicates, 1, "clipping must not hide a duplicate");
    assert!(duplicate_flags(&output)[1].1, "the clipped copy is the duplicate");
}

/// The same start, and a different mate. These are two separate molecules that share one endpoint
/// by chance. A signature that left out the mate would put them together, and that would delete
/// real coverage.
#[test]
fn fragments_sharing_one_end_but_not_the_other_are_not_duplicates() {
    let dir = scratch("mate");
    let hdr = header(&[("chr1", 100_000)]);
    let records = vec![
        pair_record("a", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("b", 0, 100, false, true, 0, 900, true, "10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (_, stats) = mark(&input, &dir, MarkDupParams::default());
    assert_eq!(stats.duplicates, 0, "different mates means different molecules");
}

/// Same endpoints on opposite strands are different molecules, not copies.
#[test]
fn opposite_strands_are_not_duplicates() {
    let dir = scratch("strand");
    let hdr = header(&[("chr1", 100_000)]);
    let records = vec![
        pair_record("fwd", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("rev", 0, 100, true, true, 0, 500, true, "10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (_, stats) = mark(&input, &dir, MarkDupParams::default());
    assert_eq!(stats.duplicates, 0);
}

/// A long-read library usually has no PCR step, and two long reads rarely share an endpoint by
/// chance. A mark on them would then throw away real coverage. The option that turns the mark off
/// must turn it off.
#[test]
fn marking_can_be_turned_off_for_long_reads() {
    let dir = scratch("disabled");
    let hdr = header(&[("chr1", 100_000)]);
    let records = vec![
        pair_record("a", 0, 100, false, true, 0, 500, true, "10M"),
        pair_record("b", 0, 100, false, true, 0, 500, true, "10M"),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (output, stats) = mark(
        &input,
        &dir,
        MarkDupParams {
            enabled: false,
            ..Default::default()
        },
    );
    assert!(stats.skipped);
    assert_eq!(stats.duplicates, 0);
    assert!(duplicate_flags(&output).iter().all(|(_, d)| !d));
}

/// An unmapped record, a secondary record and a supplementary record represent no molecule of
/// their own, because the primary record already does. They go through with no mark, and the count
/// puts them outside the set that the pass can mark.
#[test]
fn ineligible_records_pass_through_unmarked() {
    let dir = scratch("ineligible");
    let hdr = header(&[("chr1", 100_000)]);
    let mut secondary = pair_record("sec", 0, 100, false, true, 0, 500, true, "10M");
    *secondary.flags_mut() = Flags::from(0x1u16 | 0x40 | 0x100);
    let mut supplementary = pair_record("sup", 0, 100, false, true, 0, 500, true, "10M");
    *supplementary.flags_mut() = Flags::from(0x1u16 | 0x40 | 0x800);

    let records = vec![
        pair_record("a", 0, 100, false, true, 0, 500, true, "10M"),
        secondary,
        supplementary,
        record("unmapped", None, 0),
    ];
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &records);

    let (output, stats) = mark(&input, &dir, MarkDupParams::default());
    assert_eq!(stats.ineligible, 3);
    assert_eq!(stats.duplicates, 0, "none of them may be marked");
    assert!(duplicate_flags(&output).iter().all(|(_, d)| !d));
}

/// A second run, and an input that a vendor already marked, must both take the answer of this
/// pass. Neither may keep an answer that this pass did not reach.
#[test]
fn pre_existing_duplicate_flags_are_recomputed() {
    let dir = scratch("recompute");
    let hdr = header(&[("chr1", 100_000)]);
    // Flagged on input, but not a duplicate of anything.
    let mut stale = pair_record("stale", 0, 100, false, true, 0, 500, true, "10M");
    *stale.flags_mut() = Flags::from(0x1u16 | 0x40 | 0x400);
    let input = dir.join("in.bam");
    write_bam(&input, &hdr, &[stale]);

    let (output, stats) = mark(&input, &dir, MarkDupParams::default());
    assert_eq!(stats.duplicates, 0);
    assert!(
        !duplicate_flags(&output)[0].1,
        "a stale flag must be cleared, not carried"
    );
}

/// The mark never drops a record. It changes a flag and nothing else.
#[test]
fn marking_is_lossless() {
    let dir = scratch("mdlossless");
    let (input, _, total) = unsorted_fixture(&dir);
    let sorted_dir = dir.join("sorted");
    std::fs::create_dir_all(&sorted_dir).unwrap();
    let (sorted, _) = run(&input, &sorted_dir, SortParams::default());

    let (output, stats) = mark(&sorted, &dir, MarkDupParams::default());
    assert_eq!(stats.records as usize, total);
    assert_eq!(duplicate_flags(&output).len(), total);
}

// ---- CRAM emit and index --------------------------------------------------

use super::cram::{crai_path, index_cram, write_cram};

/// A CRAM stores a read as the difference from the reference, so a test needs a real reference.
/// That means bases that the records match, and the `.fai` that the repository reads. A stub
/// reference would either fail to build, or encode a mismatch at every base where nobody sees
/// it.
fn write_reference_fasta(dir: &Path, name: &str, len: usize) -> PathBuf {
    let path = dir.join("ref.fa");
    let bases: Vec<u8> = (0..len).map(|i| b"ACGT"[(i * 7 + 3) % 4]).collect();

    let mut text = Vec::new();
    text.extend_from_slice(format!(">{name}\n").as_bytes());
    let line_width = 60;
    for chunk in bases.chunks(line_width) {
        text.extend_from_slice(chunk);
        text.push(b'\n');
    }
    std::fs::write(&path, &text).unwrap();

    // A small `.fai`. Its five fields are the name, the length, the offset of the first base, the
    // bases in a line, and the bytes in a line.
    let offset = name.len() + 2; // ">name\n"
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!("{name}\t{len}\t{offset}\t{line_width}\t{}\n", line_width + 1),
    )
    .unwrap();
    path
}

fn reference_bases_at(len: usize) -> Vec<u8> {
    (0..len).map(|i| b"ACGT"[(i * 7 + 3) % 4]).collect()
}

/// A record whose sequence matches the reference at `pos`. A CRAM then has nothing to store but
/// the position. That is the case that this test needs, because a fixture with many mismatches
/// would not cover the compression against a reference at all.
fn matching_record(name: &str, pos: usize, len: usize, reference: &[u8]) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::from(0u16))
        .set_reference_sequence_id(0)
        .set_alignment_start(noodles::core::Position::new(pos).unwrap())
        .set_cigar(parse_cigar(&format!("{len}M")))
        .set_sequence(Sequence::from(reference[pos - 1..pos - 1 + len].to_vec()))
        .set_quality_scores(QualityScores::from(vec![30; len]))
        .build()
}

/// A sorted BAM, and the reference that the mapper aligned it to.
fn cram_fixture(dir: &Path, count: usize) -> (PathBuf, PathBuf, usize) {
    let contig_len = 10_000;
    let reference = write_reference_fasta(dir, "chr1", contig_len);
    let bases = reference_bases_at(contig_len);

    let hdr = header(&[("chr1", contig_len)]);
    let mut records: Vec<RecordBuf> = (0..count)
        .map(|i| matching_record(&format!("r{i:03}"), 1 + i * 50, 50, &bases))
        .collect();
    // Coordinate order is the condition. The fixture is already in that order. But write the
    // header the way that the sort would, so that the check under test sees what it expects.
    records.sort_by_key(|r| r.alignment_start().map(|p| p.get()).unwrap_or(0));

    let bam = dir.join("sorted.bam");
    write_bam_sorted(&bam, &hdr, &records);
    (bam, reference, count)
}

/// The same as [`write_bam`], and it also writes `@HD SO:coordinate`, as the sort does.
fn write_bam_sorted(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    use noodles::sam::header::record::value::map::header::tag;
    use noodles::sam::header::record::value::{map, Map};

    let mut header = header.clone();
    let mut hd = Map::<map::Header>::default();
    hd.other_fields_mut()
        .insert(tag::SORT_ORDER, b"coordinate".as_slice().into());
    *header.header_mut() = Some(hd);
    write_bam(path, &header, records);
}

/// The main property. A sorted BAM becomes a CRAM that reads back with the same records, and an
/// index sits beside it.
#[test]
fn cram_round_trips_every_record_and_writes_an_index() {
    let dir = scratch("cram");
    let (bam, reference, count) = cram_fixture(&dir, 40);
    let out = dir.join("out.cram");

    let result =
        write_cram(&bam, &out, &reference, &CancelToken::none(), &mut |_| {}).expect("CRAM emit should succeed");

    assert_eq!(result.records as usize, count);
    assert!(out.is_file());
    assert_eq!(result.index, crai_path(&out));
    assert!(result.index.is_file(), "the .crai must be written beside the CRAM");
    assert!(
        std::fs::metadata(&result.index).unwrap().len() > 0,
        "an empty index would leave every query scanning the whole file"
    );

    // Read it back through the same path that Navigator uses for a vendor CRAM. If no reader can
    // take the realigned output that way, then no analysis in the app can use it.
    let (header, mut reader) = crate::reader::open_seq(&out, Some(&reference)).unwrap();
    let names: Vec<String> = reader
        .records(&header)
        .map(|r| String::from_utf8_lossy(r.unwrap().name().unwrap_or_default()).to_string())
        .collect();
    assert_eq!(names.len(), count, "no record lost in compression");
    assert_eq!(names[0], "r000");
}

/// A CRAM builds the bases again from the reference. If that round trip were wrong, the sequences
/// would come back changed, and the read itself would not fail. So this test compares the bases
/// one by one.
#[test]
fn sequences_survive_reference_based_compression() {
    let dir = scratch("crambases");
    let (bam, reference, _) = cram_fixture(&dir, 10);
    let out = dir.join("out.cram");
    write_cram(&bam, &out, &reference, &CancelToken::none(), &mut |_| {}).unwrap();

    let expected = reference_bases_at(10_000);
    let (header, mut reader) = crate::reader::open_seq(&out, Some(&reference)).unwrap();
    for result in reader.records(&header) {
        let record = result.unwrap();
        let start = record.alignment_start().unwrap().get();
        let seq = record.sequence().as_ref().to_vec();
        assert_eq!(
            seq,
            expected[start - 1..start - 1 + seq.len()].to_vec(),
            "CRAM must give back the bases that went in"
        );
    }
}

/// Input in read order makes a file that is slow to write, larger than the BAM, and of no use for
/// a region query. To refuse it at the start is better than to find that out after hours.
#[test]
fn unsorted_input_is_refused_before_compressing() {
    let dir = scratch("cramunsorted");
    let contig_len = 10_000;
    let reference = write_reference_fasta(&dir, "chr1", contig_len);
    let bases = reference_bases_at(contig_len);

    // Written without the sort's @HD SO stamp.
    let hdr = header(&[("chr1", contig_len)]);
    let records = vec![matching_record("a", 1, 50, &bases)];
    let bam = dir.join("unsorted.bam");
    write_bam(&bam, &hdr, &records);

    let out = dir.join("out.cram");
    let err = write_cram(&bam, &out, &reference, &CancelToken::none(), &mut |_| {});
    assert!(err.is_err(), "unsorted input must be refused");
    match err.unwrap_err() {
        AnalysisError::Message(m) => assert!(m.contains("coordinate"), "unhelpful message: {m}"),
        other => panic!("expected a clear message, got {other:?}"),
    }
}

/// Indexing reads the finished file back, so it is also the repair path for a CRAM whose index was
/// lost or truncated. It must work standalone.
#[test]
fn an_index_can_be_rebuilt_from_a_finished_cram() {
    let dir = scratch("reindex");
    let (bam, reference, _) = cram_fixture(&dir, 20);
    let out = dir.join("out.cram");
    let result = write_cram(&bam, &out, &reference, &CancelToken::none(), &mut |_| {}).unwrap();

    let original = std::fs::read(&result.index).unwrap();
    std::fs::remove_file(&result.index).unwrap();

    let rebuilt = index_cram(&out).unwrap();
    assert_eq!(rebuilt, result.index);
    assert_eq!(
        std::fs::read(&rebuilt).unwrap(),
        original,
        "a rebuilt index must match the one written alongside"
    );
}

/// The name every reader looks for.
#[test]
fn the_index_sits_beside_the_cram() {
    assert_eq!(
        crai_path(Path::new("/data/sample.cram")),
        Path::new("/data/sample.cram.crai")
    );
}

/// Compression is long enough to need a cancel that works.
#[test]
fn cancellation_stops_cram_emission() {
    let dir = scratch("cramcancel");
    let (bam, reference, _) = cram_fixture(&dir, 10);
    let cancel = CancelToken::new();
    cancel.cancel();

    let err = write_cram(&bam, &dir.join("out.cram"), &reference, &cancel, &mut |_| {}).unwrap_err();
    assert!(matches!(err, AnalysisError::Cancelled));
}

/// A secondary alignment with `SEQ: *`. That is legal SAM, and it is what minimap2 gives, because
/// the primary alignment alone carries the bases. The code can not encode it as a difference from
/// the reference. So it drops that record and counts it. Without that, the writer would panic from
/// inside noodles.
#[test]
fn cram_drops_a_secondary_record_that_carries_no_sequence() {
    let dir = scratch("cram-secondary");
    let contig_len = 10_000;
    let reference = write_reference_fasta(&dir, "chr1", contig_len);
    let bases = reference_bases_at(contig_len);
    let hdr = header(&[("chr1", contig_len)]);

    let primary = matching_record("r000", 1, 65, &bases);
    let secondary = RecordBuf::builder()
        .set_name("r000")
        .set_flags(Flags::SECONDARY)
        .set_reference_sequence_id(0)
        .set_alignment_start(noodles::core::Position::new(500).unwrap())
        .set_cigar(parse_cigar("65M"))
        .build();

    let bam = dir.join("sorted.bam");
    write_bam_sorted(&bam, &hdr, &[primary, secondary]);

    let out = dir.join("out.cram");
    let result = write_cram(&bam, &out, &reference, &CancelToken::none(), &mut |_| {}).expect("CRAM emit");

    assert_eq!(result.records, 1, "the primary is written");
    assert_eq!(result.sequenceless_dropped, 1, "the secondary is dropped, and counted");
}

/// The same shape on a *primary* alignment means that a read goes missing. It is not a record that
/// says nothing new. So the code fails loudly there, and it does not put that record into a count
/// of drops.
#[test]
fn cram_refuses_a_primary_record_that_carries_no_sequence() {
    let dir = scratch("cram-primary-noseq");
    let contig_len = 10_000;
    let reference = write_reference_fasta(&dir, "chr1", contig_len);
    let hdr = header(&[("chr1", contig_len)]);

    let orphan = RecordBuf::builder()
        .set_name("r000")
        .set_flags(Flags::from(0u16))
        .set_reference_sequence_id(0)
        .set_alignment_start(noodles::core::Position::new(1).unwrap())
        .set_cigar(parse_cigar("65M"))
        .build();

    let bam = dir.join("sorted.bam");
    write_bam_sorted(&bam, &hdr, &[orphan]);

    let err = write_cram(
        &bam,
        &dir.join("out.cram"),
        &reference,
        &CancelToken::none(),
        &mut |_| {},
    )
    .expect_err("a primary with no SEQ must not be dropped silently");
    assert!(
        format!("{err}").contains("no SEQ"),
        "the error should name the cause: {err}"
    );
}

/// The last step moves the marked BAM into place and makes its index. This test covers *more than
/// one* contig, which is the shape that the code can not index as a CRAM. Every slice that crosses
/// a contig boundary holds more than one reference, and `cram::fs::index` decodes those against an
/// empty reference repository.
#[test]
fn finalizing_moves_the_bam_into_place_and_indexes_it() {
    use super::finalize::{bai_path, finalize_bam};

    let dir = scratch("finalize");
    let contig_len = 10_000;
    let bases = reference_bases_at(contig_len);
    let hdr = header(&[("chr1", contig_len), ("chr2", contig_len), ("chrM", contig_len)]);

    let mut records: Vec<RecordBuf> = Vec::new();
    for contig in 0..3usize {
        for i in 0..20usize {
            let mut record = matching_record(&format!("r{contig}{i:03}"), 1 + i * 50, 50, &bases);
            *record.reference_sequence_id_mut() = Some(contig);
            records.push(record);
        }
    }

    let marked = dir.join("marked.bam");
    write_bam_sorted(&marked, &hdr, &records);

    let out = dir.join("alignment-1.chm13v2.0.bam");
    let finalized = finalize_bam(&marked, &out).expect("finalize");

    assert!(!marked.exists(), "the intermediate is moved, not copied");
    assert_eq!(finalized.bam, out);
    assert_eq!(finalized.index, bai_path(&out));
    assert!(finalized.index.exists(), "a .bai sits beside the alignment");

    let (_, read_back) = {
        let file = std::fs::File::open(&out).unwrap();
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header().unwrap();
        let records: Vec<_> = reader.record_bufs(&header).map(|r| r.unwrap()).collect();
        (header, records)
    };
    assert_eq!(read_back.len(), 60, "every record survives finalising");
}

// ---- is the file complete -------------------------------------------------
//
// The realignment resume in `navigator-app` trusts `is_complete_bam`. That answer decides one
// thing. Can the app take up the 60 GB intermediate of a run that somebody killed, or must it
// derive that file again over some hours? A wrong answer in either direction costs much, so this
// section holds both directions.

#[test]
fn a_finished_bam_is_complete() {
    let dir = scratch("complete-finished");
    let hdr = header(&[("chr1", 1000)]);
    let path = dir.join("finished.bam");
    write_bam(&path, &hdr, &[record("a", Some(0), 10)]);

    assert!(crate::postprocess::is_complete_bam(&path));
}

/// This is the case that the check exists for: somebody killed a writer in the middle of its
/// stream. The file is large, it looks correct, and a reader can read it up to the cut. Only the
/// end-of-file marker separates it from a good file.
#[test]
fn a_truncated_bam_is_not_complete() {
    let dir = scratch("complete-truncated");
    let hdr = header(&[("chr1", 1000)]);
    let path = dir.join("whole.bam");
    write_bam(&path, &hdr, &[record("a", Some(0), 10), record("b", Some(0), 20)]);

    let cut = dir.join("cut.bam");
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&cut, &bytes[..bytes.len() - 8]).unwrap();

    assert!(!crate::postprocess::is_complete_bam(&cut));
}

#[test]
fn a_missing_or_empty_file_is_not_complete() {
    let dir = scratch("complete-missing");
    assert!(!crate::postprocess::is_complete_bam(&dir.join("nothing.bam")));

    let empty = dir.join("empty.bam");
    std::fs::write(&empty, b"").unwrap();
    assert!(!crate::postprocess::is_complete_bam(&empty));
}

/// The writers that the pipeline uses go through `bamio`, which now paces its flushes and its
/// syncs on the way out. Whatever that path does, the file that it leaves must read as
/// complete.
#[test]
fn what_the_pipeline_writes_reads_as_complete() {
    let dir = scratch("complete-pipeline");
    let (input, _, _) = unsorted_fixture(&dir);
    let (sorted, _) = run(&input, &dir.join("runs"), SortParams { buffer_bytes: 1 });

    assert!(
        crate::postprocess::is_complete_bam(&sorted),
        "the sort's own output must satisfy the predicate resume checks"
    );
}

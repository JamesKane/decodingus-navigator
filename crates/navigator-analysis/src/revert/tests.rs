//! Tests for the revert stage.
//!
//! Records are built by hand rather than read from BAM fixtures: the container layer is
//! [`crate::reader`]'s responsibility and is covered there, and everything interesting here — flag
//! handling, orientation, pairing, the external sort's spill/merge boundary — is independent of
//! how the bytes arrived. Building records directly also lets a test construct the malformed
//! inputs that matter (contradictory flags, absent qualities) which a real aligner rarely emits.

use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;

use super::*;

const OQ: Tag = Tag::new(b'O', b'Q');

/// Unique scratch dir per test, under the system temp dir (matching the convention elsewhere in
/// the crate — no `tempfile` dependency).
fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-revert-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A minimal record: name, flags, sequence, and qualities matching the sequence length.
fn record(name: &str, flags: u16, seq: &str, quals: &[u8]) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::from(flags))
        .set_sequence(Sequence::from(seq.as_bytes().to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build()
}

fn paired(name: &str, first: bool, reverse: bool, seq: &str, quals: &[u8]) -> RecordBuf {
    let mut flags = 0x1u16 | if first { 0x40 } else { 0x80 };
    if reverse {
        flags |= 0x10;
    }
    record(name, flags, seq, quals)
}

fn run(records: Vec<RecordBuf>, dir: &Path, params: &RevertParams) -> RevertOutput {
    revert_records(records.into_iter().map(Ok), dir, params, &CancelToken::none()).unwrap()
}

fn read_lines(path: &Path) -> Vec<String> {
    use std::io::BufRead;
    // The reverted FASTQ is gzipped (see writer.rs on why), so tests read it the same way the
    // mapper does rather than assuming plain text.
    crate::gzio::open_maybe_gz(path)
        .unwrap()
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

// ---- the transform --------------------------------------------------------

/// The core correctness property: a reverse-strand alignment stores the read flipped, and the
/// FASTQ has to carry what the sequencer produced, not what the aligner stored.
#[test]
fn a_reverse_strand_read_is_restored_to_sequencer_orientation() {
    let dir = scratch("revcomp");
    let out = run(
        vec![record("r", 0x10, "AACCG", &[10, 20, 30, 40, 50])],
        &dir,
        &RevertParams::default(),
    );

    let lines = read_lines(&out.singletons);
    assert_eq!(lines[1], "CGGTT", "reverse-complemented back");
    // Qualities travel with the bases they describe, so they reverse too (but are not complemented).
    assert_eq!(
        lines[3].bytes().map(|b| b - 33).collect::<Vec<_>>(),
        vec![50, 40, 30, 20, 10]
    );
}

/// A forward-strand read must be passed through untouched — the mirror of the test above, so a
/// bug that reverse-complements unconditionally cannot pass both.
#[test]
fn a_forward_strand_read_is_left_alone() {
    let dir = scratch("forward");
    let out = run(
        vec![record("r", 0x0, "AACCG", &[10, 20, 30, 40, 50])],
        &dir,
        &RevertParams::default(),
    );
    let lines = read_lines(&out.singletons);
    assert_eq!(lines[1], "AACCG");
    assert_eq!(
        lines[3].bytes().map(|b| b - 33).collect::<Vec<_>>(),
        vec![10, 20, 30, 40, 50]
    );
}

/// Secondary and supplementary records duplicate a read whose full sequence lives on the primary;
/// keeping them would emit the same read more than once and, for supplementaries, truncated.
#[test]
fn secondary_and_supplementary_records_are_dropped() {
    let dir = scratch("nonprimary");
    let out = run(
        vec![
            record("r", 0x0, "ACGT", &[30; 4]),
            record("r", 0x100, "ACG", &[30; 3]),
            record("r", 0x800, "CGT", &[30; 3]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.secondary_dropped, 1);
    assert_eq!(out.stats.supplementary_dropped, 1);
    assert_eq!(out.stats.reads_emitted, 1, "only the primary survives");
}

/// Unmapped reads are the realignment payoff — they must flow through, and be counted so the
/// payoff is measurable.
#[test]
fn unmapped_reads_are_kept_and_counted() {
    let dir = scratch("unmapped");
    let out = run(vec![record("u", 0x4, "ACGT", &[30; 4])], &dir, &RevertParams::default());
    assert_eq!(out.stats.unmapped_reads, 1);
    assert_eq!(out.stats.reads_emitted, 1);
}

/// `OQ` holds the qualities from before recalibration overwrote `QUAL`, ASCII-encoded. Preferring
/// it is the difference between reverting to the original read and reverting to a processed one.
#[test]
fn original_qualities_are_preferred_over_recalibrated_ones() {
    let dir = scratch("oq");
    let rec = RecordBuf::builder()
        .set_name("r")
        .set_flags(Flags::from(0x0u16))
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![2, 2, 2, 2]))
        // Phred 40 as FASTQ ASCII ('I' == 40 + 33).
        .set_data([(OQ, Value::String("IIII".into()))].into_iter().collect::<Data>())
        .build();

    let out = run(vec![rec], &dir, &RevertParams::default());
    assert_eq!(out.stats.original_qualities_used, 1);
    assert_eq!(read_lines(&out.singletons)[3], "IIII", "OQ won, decoded and re-encoded");
}

/// Opting out has to actually opt out, or the flag is decoration.
#[test]
fn original_qualities_can_be_declined() {
    let dir = scratch("oq-off");
    let rec = RecordBuf::builder()
        .set_name("r")
        .set_flags(Flags::from(0x0u16))
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![2, 2, 2, 2]))
        .set_data([(OQ, Value::String("IIII".into()))].into_iter().collect::<Data>())
        .build();

    let params = RevertParams {
        prefer_original_qualities: false,
        ..Default::default()
    };
    let out = run(vec![rec], &dir, &params);
    assert_eq!(out.stats.original_qualities_used, 0);
    assert_eq!(read_lines(&out.singletons)[3], "####", "QUAL 2 == '#', one per base");
}

/// A hard-clipped primary has already lost sequence. Skipping is the default because a dropped
/// read shows up in the stats and a truncated one does not.
#[test]
fn hard_clipped_primaries_are_skipped_by_default_and_emittable_on_request() {
    let rec = || {
        RecordBuf::builder()
            .set_name("h")
            .set_flags(Flags::from(0x0u16))
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; 4]))
            .set_cigar(Cigar::from(vec![Op::new(Kind::HardClip, 10), Op::new(Kind::Match, 4)]))
            .build()
    };

    let skipped = run(vec![rec()], &scratch("hardclip-skip"), &RevertParams::default());
    assert_eq!(skipped.stats.hard_clipped, 1);
    assert_eq!(skipped.stats.reads_emitted, 0);

    let params = RevertParams {
        hard_clipped: HardClipPolicy::Emit,
        ..Default::default()
    };
    let emitted = run(vec![rec()], &scratch("hardclip-emit"), &params);
    assert_eq!(emitted.stats.hard_clipped, 1, "counted under either policy");
    assert_eq!(emitted.stats.reads_emitted, 1);
}

/// `QUAL` of `*` is legal. The read is still mappable, so it is kept — but the qualities that come
/// out are invented, and the stat is the only thing that says so.
#[test]
fn missing_qualities_are_synthesized_and_counted() {
    let dir = scratch("noqual");
    let out = run(vec![record("r", 0x0, "ACGT", &[])], &dir, &RevertParams::default());
    assert_eq!(out.stats.qualities_synthesized, 1);
    let lines = read_lines(&out.singletons);
    assert_eq!(lines[3].len(), 4, "one quality per base, as FASTQ requires");
}

// ---- pairing --------------------------------------------------------------

/// The headline case: mates arrive far apart in coordinate order and must come back together.
#[test]
fn mates_separated_in_the_input_are_paired_in_the_output() {
    let dir = scratch("pairing");
    // Interleave three templates so no pair is adjacent on input.
    let out = run(
        vec![
            paired("a", true, false, "AAAA", &[30; 4]),
            paired("b", true, false, "BBBB", &[30; 4]),
            paired("c", true, false, "CCCC", &[30; 4]),
            paired("b", false, false, "TTTT", &[30; 4]),
            paired("c", false, false, "GGGG", &[30; 4]),
            paired("a", false, false, "CCCC", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.pairs, 3);
    assert_eq!(out.stats.singletons, 0);

    let r1 = read_lines(&out.read1);
    let r2 = read_lines(&out.read2);
    // Name order, and — the invariant the mapper depends on — R1 and R2 in lockstep.
    assert_eq!([r1[0].as_str(), r1[4].as_str(), r1[8].as_str()], ["@a", "@b", "@c"]);
    assert_eq!(r1.len(), r2.len(), "files stay the same length");
    for i in (0..r1.len()).step_by(4) {
        assert_eq!(r1[i], r2[i], "record {} pairs the same template", i / 4);
    }
}

/// A mate whose partner was dropped must not be written into `_1` — doing so would shift every
/// later pair by one and mis-pair the rest of the file.
#[test]
fn a_read_whose_mate_was_dropped_becomes_a_singleton() {
    let dir = scratch("orphan");
    let out = run(
        vec![
            paired("a", true, false, "AAAA", &[30; 4]),
            paired("a", false, false, "TTTT", &[30; 4]),
            // 'b' has only its first segment; its mate never appears.
            paired("b", true, false, "GGGG", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.pairs, 1);
    assert_eq!(out.stats.singletons, 1);
    assert_eq!(read_lines(&out.read1).len(), 4, "only the complete pair");
    assert_eq!(read_lines(&out.singletons)[0], "@b");
}

/// Flags that claim "paired" but not which end cannot be placed in a synchronized file.
#[test]
fn a_paired_record_with_contradictory_segment_flags_is_a_singleton() {
    let dir = scratch("contradictory");
    // 0x1 segmented, with both 0x40 and 0x80 set.
    let out = run(
        vec![record("x", 0x1 | 0x40 | 0x80, "ACGT", &[30; 4])],
        &dir,
        &RevertParams::default(),
    );
    assert_eq!(out.stats.pairs, 0);
    assert_eq!(out.stats.singletons, 1);
}

/// Two records claiming the same name *and* the same end is unresolvable — picking one would
/// silently drop a read and could mis-pair the template.
#[test]
fn duplicate_segment_bits_under_one_name_do_not_pair() {
    let dir = scratch("dupe-segment");
    let out = run(
        vec![
            paired("d", true, false, "AAAA", &[30; 4]),
            paired("d", true, false, "CCCC", &[30; 4]),
            paired("d", false, false, "GGGG", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );
    assert_eq!(out.stats.pairs, 0);
    assert_eq!(out.stats.singletons, 3, "all three go to singletons");
}

// ---- the external sort ----------------------------------------------------

/// The property the whole design rests on: the result must not depend on whether the input fit in
/// memory. A budget of 1 byte forces a spill per read and exercises the k-way merge; the output
/// has to be identical to the single-run case.
#[test]
fn spilling_to_disk_produces_the_same_output_as_sorting_in_memory() {
    let records = || {
        let mut v = Vec::new();
        for i in 0..50 {
            // Names deliberately out of lexical order relative to emission order.
            let name = format!("read{:03}", (i * 37) % 50);
            v.push(paired(&name, true, false, "ACGT", &[30; 4]));
            v.push(paired(&name, false, true, "ACGT", &[30; 4]));
        }
        v
    };

    let in_memory = run(
        records(),
        &scratch("sort-memory"),
        &RevertParams {
            sort_buffer_bytes: 64 * 1024 * 1024,
            ..Default::default()
        },
    );
    let spilled = run(
        records(),
        &scratch("sort-spill"),
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    assert_eq!(in_memory.stats.runs_spilled, 1, "everything fit");
    assert!(spilled.stats.runs_spilled > 1, "actually spilled multiple runs");
    assert_eq!(in_memory.stats.pairs, 50);
    assert_eq!(spilled.stats.pairs, in_memory.stats.pairs);
    assert_eq!(
        read_lines(&spilled.read1),
        read_lines(&in_memory.read1),
        "identical output regardless of how much fit in memory"
    );
    assert_eq!(read_lines(&spilled.read2), read_lines(&in_memory.read2));
}

/// Names must come out in sorted order after a merge across many runs — the grouping logic reads
/// runs of equal names, so an unsorted merge would silently split templates.
#[test]
fn the_merge_emits_names_in_sorted_order() {
    let dir = scratch("sorted");
    let mut records = Vec::new();
    for i in (0..30).rev() {
        records.push(record(&format!("r{i:03}"), 0x0, "ACGT", &[30; 4]));
    }
    let out = run(
        records,
        &dir,
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    let names: Vec<String> = read_lines(&out.singletons).into_iter().step_by(4).collect();
    let mut sorted = names.clone();
    sorted.sort();
    assert_eq!(names, sorted);
    assert_eq!(names.len(), 30);
}

/// Scratch files are large — tens of GB for a WGS — so they must not outlive the run.
#[test]
fn run_files_are_cleaned_up() {
    let dir = scratch("cleanup");
    let mut records = Vec::new();
    for i in 0..20 {
        records.push(record(&format!("r{i:03}"), 0x0, "ACGT", &[30; 4]));
    }
    let _ = run(
        records,
        &dir,
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("revert-run-"))
        .collect();
    assert!(leftovers.is_empty(), "spill files left behind: {leftovers:?}");
}

// ---- cancellation ---------------------------------------------------------

/// A revert is an hours-long job; an already-cancelled token must stop it rather than run to
/// completion, and must report itself as cancelled rather than as a failure.
#[test]
fn an_already_cancelled_token_stops_the_revert() {
    let dir = scratch("cancel");
    let cancel = CancelToken::new();
    cancel.cancel();

    let records = (0..10).map(|i| Ok(record(&format!("r{i}"), 0x0, "ACGT", &[30; 4])));
    let err = revert_records(records, &dir, &RevertParams::default(), &cancel).unwrap_err();
    assert!(matches!(err, AnalysisError::Cancelled));
}

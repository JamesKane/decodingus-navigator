//! End-to-end revert over real containers, against `paired.bam` / `paired.cram`.
//!
//! The unit tests in `src/revert/tests.rs` give the records directly. This test closes the gap
//! that they leave: the same pipeline must work when the records come out of a real BAM or CRAM
//! decode.
//!
//! `paired.bam` suits this by accident of how somebody built it. It holds two FR pairs, at chrM:1
//! and 31, and at chrM:5 and 25, in coordinate order. The file order is then pairA, pairB, pairB,
//! pairA, and the two mates of a template never sit beside each other. That is the exact condition
//! that the collation exists for.
//!
//! The `/2` records carry flag 147, and that flag holds `0x10`. So the restore of the reverse
//! complement also runs on real decoded records.

use std::path::{Path, PathBuf};

use navigator_analysis::cancel::CancelToken;
use navigator_analysis::revert::{revert_alignment, RevertParams};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-revert-it-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn lines(path: &Path) -> Vec<String> {
    use std::io::BufRead;
    // Gzipped output; read it the way the mapper does.
    navigator_analysis::gzio::open_maybe_gz(path)
        .unwrap()
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

#[test]
fn reverts_a_coordinate_sorted_bam_into_synchronized_pairs() {
    let dir = scratch("bam");
    let out = revert_alignment(
        &fixtures().join("paired.bam"),
        None,
        &dir,
        &RevertParams::default(),
        &CancelToken::none(),
    )
    .expect("revert should succeed");

    assert_eq!(out.stats.records_read, 4);
    assert_eq!(out.stats.pairs, 2, "both templates re-paired across the file");
    assert_eq!(out.stats.singletons, 0);
    assert_eq!(out.stats.reads_emitted, 4, "no read lost in the round trip");

    let r1 = lines(&out.read1);
    let r2 = lines(&out.read2);
    assert_eq!(r1.len(), 8, "two FASTQ records");
    assert_eq!(r1.len(), r2.len(), "R1 and R2 stay in lockstep");
    assert_eq!([r1[0].as_str(), r1[4].as_str()], ["@pairA", "@pairB"]);
    assert_eq!([r2[0].as_str(), r2[4].as_str()], ["@pairA", "@pairB"]);

    // The reads of the fixture are poly-A, and the /1 records store them forward. The /2 records
    // carry the reverse flag. What the file holds as poly-A comes back as poly-T, once the code
    // restores the orientation of the sequencer.
    assert_eq!(r1[1], "AAAAAAAAAA", "/1 is forward, passed through");
    assert_eq!(r2[1], "TTTTTTTTTT", "/2 is reverse-flagged, so it is complemented back");
}

#[test]
fn cram_reverts_identically_to_bam() {
    // The same reads, in a different container. The revert must not be able to see a difference. A
    // CRAM needs the reference to decode, and that is the one difference in the path that this
    // test must cover.
    let dir = fixtures();
    let from_bam = revert_alignment(
        &dir.join("paired.bam"),
        None,
        &scratch("parity-bam"),
        &RevertParams::default(),
        &CancelToken::none(),
    )
    .unwrap();
    let from_cram = revert_alignment(
        &dir.join("paired.cram"),
        Some(&dir.join("ref.fa")),
        &scratch("parity-cram"),
        &RevertParams::default(),
        &CancelToken::none(),
    )
    .unwrap();

    assert_eq!(from_cram.stats, from_bam.stats, "CRAM revert stats must equal BAM");
    assert_eq!(lines(&from_cram.read1), lines(&from_bam.read1));
    assert_eq!(lines(&from_cram.read2), lines(&from_bam.read2));
}

/// The same input, through a sort budget of one byte, must give the same FASTQ. This holds on real
/// decoded records, and not on synthetic ones alone.
#[test]
fn spilling_does_not_change_the_result_on_a_real_bam() {
    let bam = fixtures().join("paired.bam");
    let in_memory = revert_alignment(
        &bam,
        None,
        &scratch("spill-none"),
        &RevertParams::default(),
        &CancelToken::none(),
    )
    .unwrap();
    let spilled = revert_alignment(
        &bam,
        None,
        &scratch("spill-yes"),
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
        &CancelToken::none(),
    )
    .unwrap();

    assert!(spilled.stats.runs_spilled > 1, "the budget forced multiple runs");
    assert_eq!(lines(&spilled.read1), lines(&in_memory.read1));
    assert_eq!(lines(&spilled.read2), lines(&in_memory.read2));
}

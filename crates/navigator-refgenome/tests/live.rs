//! Live network tests. Each one carries `#[ignore]`, as the live-PDS tests do, and CI never runs
//! them. They go to the real public reference hosts. Run one of them explicitly, for example:
//!
//!   `cargo test -p navigator-refgenome --test live -- --ignored resolve_chm13 --nocapture`
//!
//! Note: the reference FASTA download is about 1 GB.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use navigator_refgenome::ReferenceGateway;

fn scratch() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-refgenome-live-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

#[tokio::test]
#[ignore = "downloads ~1 GB from the human-pangenomics bucket"]
async fn resolve_chm13() {
    let base = scratch();
    let g = ReferenceGateway::new(base.clone(), reqwest::Client::new());
    let last = AtomicU64::new(0);
    let mut progress = |received: u64, total: Option<u64>| {
        // Log roughly each 50 MB so a manual run shows movement.
        if received - last.load(Ordering::Relaxed) > 50_000_000 {
            last.store(received, Ordering::Relaxed);
            eprintln!("  {} MB / {:?}", received / 1_000_000, total.map(|t| t / 1_000_000));
        }
    };
    let fa = g.resolve_reference("chm13v2.0", &mut progress).await.expect("resolve");
    assert!(fa.exists());
    // The .fai must exist and hold chr1. CHM13 puts a `chr` prefix on a contig name.
    let fai = std::fs::read_to_string(fa.with_extension("fa.fai")).expect("fai");
    assert!(fai.lines().any(|l| l.starts_with("chr1\t")), "expected chr1 in .fai");
    let _ = std::fs::remove_dir_all(&base);
}

#[tokio::test]
#[ignore = "downloads a real liftover chain"]
async fn resolve_grch38_to_chm13_chain() {
    let base = scratch();
    let g = ReferenceGateway::new(base.clone(), reqwest::Client::new());
    let chain = g
        .resolve_chain("GRCh38", "chm13v2.0", &mut |_, _| {})
        .await
        .expect("resolve chain");
    assert!(chain.exists());
    // Parses as a UCSC chain and lifts at least one coordinate.
    let lo = g.load_liftover("GRCh38", "chm13v2.0").expect("parse");
    assert!(lo.lift("chr1", 1_000_000).is_some() || lo.lift("1", 1_000_000).is_some());
    let _ = std::fs::remove_dir_all(&base);
}

/// The CHM13 chrY structural masks help on GRCh38 and GRCh37 only if they lift. They are *safe*
/// only if the code drops a bad lift, and does not spread it across the chromosome. The real chain
/// decides both, so this is a live test.
#[tokio::test]
#[ignore = "downloads a real liftover chain"]
async fn lift_chry_intervals_chm13_to_grch38() {
    let base = scratch();
    let g = ReferenceGateway::new(base.clone(), reqwest::Client::new());
    g.resolve_chain("chm13v2.0", "GRCh38", &mut |_, _| {})
        .await
        .expect("resolve chain");

    // Two real CHM13 chrY palindrome spans plus one deliberate nonsense interval far past the end
    // of the chromosome, which must not survive.
    let src = [
        (6_000_000, 6_100_000),
        (18_000_000, 18_200_000),
        (200_000_000, 200_100_000),
    ];
    let (lifted, dropped) = g.lift_intervals("chm13v2.0", "GRCh38", "chrY", &src).expect("lift");

    assert!(dropped >= 1, "the off-chromosome interval must be dropped");
    for &(s, e) in &lifted {
        assert!(e > s, "a lifted interval keeps its orientation");
        assert!(e < 60_000_000, "GRCh38 chrY is ~57.2 Mb — a lift past that is a smear");
    }
    // Whatever survived must be no more than 2x its source span, the guard against a bad lift.
    let src_total: i64 = src.iter().map(|(s, e)| e - s).sum();
    let got_total: i64 = lifted.iter().map(|(s, e)| e - s).sum();
    assert!(got_total <= src_total * 2, "lifted span must stay bounded");
    let _ = std::fs::remove_dir_all(&base);
}

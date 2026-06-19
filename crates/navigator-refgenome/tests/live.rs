//! Live network tests — `#[ignore]` (like the live-PDS tests), never run in CI. They hit
//! the real public reference hosts. Run a single one explicitly, e.g.:
//!   cargo test -p navigator-refgenome --test live -- --ignored resolve_chm13 --nocapture
//! Note: the reference FASTA download is ~1 GB.

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
    // The .fai must have been built and carry chr1 (CHM13 uses chr-prefixed names).
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

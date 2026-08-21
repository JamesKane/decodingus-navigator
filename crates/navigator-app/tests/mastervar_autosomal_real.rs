//! A separate end-to-end check. It reads a real CompleteGenomics masterVar file, makes the
//! autosomal consensus, and then estimates the ancestry.
//!
//! The test uses an **in-memory** workspace. It does not touch the `~/.decodingus/navigator-rs.db`
//! file of the user.
//!
//! The test still reads the ancestry assets and the IBD-panel assets from `~/.decodingus/ancestry`,
//! or from the path in a `NAVIGATOR_*` variable. It reads and does not write them. So the test runs
//! only on a machine that holds those assets.
//!
//! The test has the `ignore` mark, because it needs the local dump and the assets. To run it:
//!
//! ```bash
//! MASTERVAR_TSV=/path/to/var-GS00253-DNA_A01_200_37-ASM.tsv.bz2 \
//!   cargo test -p navigator-app --test mastervar_autosomal_real -- --ignored --nocapture
//! ```

use std::path::Path;
use std::time::Instant;

use navigator_app::App;
use navigator_store::Store;

#[tokio::test]
#[ignore]
async fn mastervar_feeds_autosomal_and_ancestry() {
    let Ok(path) = std::env::var("MASTERVAR_TSV") else {
        eprintln!("set MASTERVAR_TSV to run");
        return;
    };

    // In-memory workspace: the real DB is never opened, created, or written.
    let app = App::new(Store::open_in_memory().await.unwrap());
    let subject = app.add_biosample(None, "CG-MASTERVAR-TEST", None, None).await.unwrap();

    let t = Instant::now();
    let detected = app.add_data(subject.guid, Path::new(&path)).await.unwrap();
    let sets = app.list_variant_sets(subject.guid).await.unwrap();
    let total_calls: usize = sets.iter().map(|s| s.calls.len()).sum();
    println!(
        "[{:>7.1?}] import: {:?} — {} variant set(s), {} calls, build {:?}",
        t.elapsed(),
        detected,
        sets.len(),
        total_calls,
        sets.first().and_then(|s| s.reference_build.clone())
    );
    assert_eq!(detected, navigator_app::DetectedData::CompleteGenomicsVar);
    assert!(total_calls > 1_000_000, "expected a genome-wide set");

    let t = Instant::now();
    let profile = app
        .build_autosomal_profile(subject.guid)
        .await
        .expect("autosomal consensus (needs the IBD panel asset)");
    println!(
        "[{:>7.1?}] autosomal consensus: {} source(s), {} reconciled sites",
        t.elapsed(),
        profile.sources.len(),
        profile.variants.len()
    );
    for s in &profile.sources {
        println!(
            "    source: {} ({:?}) — {} sites",
            s.label, s.source_type, s.variant_count
        );
    }
    assert!(
        !profile.sources.is_empty(),
        "the masterVar should be an autosomal source"
    );
    assert!(
        profile.variants.len() > 10_000,
        "a genome-wide source should densify to a large panel overlap, got {}",
        profile.variants.len()
    );

    let t = Instant::now();
    match app.estimate_ancestry_from_consensus(subject.guid).await {
        Ok(res) => {
            println!(
                "[{:>7.1?}] ancestry ({}): {} SNPs with genotype",
                t.elapsed(),
                res.method,
                res.snps_with_genotype
            );
            println!("  super-population composition:");
            for sp in &res.super_population_summary {
                if sp.percentage >= 0.5 {
                    println!("    {:<28} {:>5.1}%", sp.super_population, sp.percentage);
                }
            }
        }
        // An absent ancestry asset must not fail the full check. The autosomal consensus is the
        // part that this change adds, and it already shows that the masterVar file reaches the
        // pipeline.
        Err(e) => println!("[{:>7.1?}] ancestry skipped: {e}", t.elapsed()),
    }
}

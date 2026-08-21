//! A harness for development. It runs the Tier B segment caller on cached diploid calls.
//!
//! This tool once also swept the thresholds, and that sweep is gone. Its output sits in the design
//! document, under the M3 calibration, and the values that won are the defaults of
//! `ArchaicConfig`.
//!
//! ```text
//! archaic_segments_probe <calls.json> <outgroup.bin> <classify.bin> <callable.bin>
//! ```
use navigator_analysis::archaic::{ArchaicCallable, ArchaicClassify, ArchaicOutgroup};
use navigator_analysis::archaic_segments::{call_archaic_segments, ArchaicConfig};
use navigator_analysis::caller::SiteGenotype;
use navigator_analysis::ibd::GeneticMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let calls: Vec<SiteGenotype> = serde_json::from_str(&std::fs::read_to_string(a.next().unwrap())?)?;
    let og = ArchaicOutgroup::from_bytes(&std::fs::read(a.next().unwrap())?).map_err(|e| e.to_string())?;
    let cls = ArchaicClassify::from_bytes(&std::fs::read(a.next().unwrap())?).map_err(|e| e.to_string())?;
    let cal = ArchaicCallable::from_bytes(&std::fs::read(a.next().unwrap())?).map_err(|e| e.to_string())?;
    println!("calls {}   callable track {:.1} Mb", calls.len(), cal.callable_mb());
    let r = call_archaic_segments(
        &calls,
        &og,
        &cls,
        &cal,
        &GeneticMap::from_markers(Vec::new()),
        &ArchaicConfig::default(),
    );
    let s = &r.summary;
    println!(
        "segments {}  total {:.2} Mb  = {:.2}% of {:.1} Mb callable",
        s.n_segments, s.total_mb, s.pct_callable, s.callable_mb
    );
    println!(
        "  Neanderthal {:.2} Mb   Denisovan {:.2} Mb   Unknown {:.2} Mb",
        s.neanderthal_mb, s.denisovan_mb, s.unknown_mb
    );
    for seg in r.segments.iter().take(6) {
        println!(
            "  {} {}-{} ({:.2} Mb) post {:.2} private {} ({:.0}/Mb) {:?} nea{} den{}",
            seg.contig,
            seg.start,
            seg.end,
            seg.length_mb(),
            seg.posterior,
            seg.n_private,
            seg.n_private as f64 / seg.length_mb().max(1e-9),
            seg.source,
            seg.neanderthal_matches,
            seg.denisovan_matches
        );
    }
    Ok(())
}

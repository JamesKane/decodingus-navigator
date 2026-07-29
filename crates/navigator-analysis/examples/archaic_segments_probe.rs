//! Development harness: run the Tier B segment caller on cached diploid calls.
//!   archaic_segments_probe <calls.json> <outgroup.bin> <classify.bin> <callable.bin>
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
    // Optional sweep: prior,multiple,min_post,min_seg_kb (repeatable via env SWEEP=1)
    if std::env::var("SWEEP").is_ok() {
        println!("{:>8} {:>6} {:>6} {:>8} {:>9} {:>7} {:>8}", "sw/cM", "mult", "post", "minkb", "segments", "Mb", "medMb");
        println!("   TARGET (hmmix EUR chr21+22):  43 segments, 2.09 Mb, median 0.031 Mb");
        for sw in [5.0f64, 20.0, 50.0] {
            for mult in [4.0f64, 5.0, 6.0] {
                for post in [0.70f64, 0.80, 0.90] {
                    for minkb in [5i64] {
                        let cfg = ArchaicConfig {
                            switches_per_cm: sw, archaic_rate_multiple: mult,
                            min_posterior: post, min_segment_bp: minkb * 1000,
                            ..ArchaicConfig::default()
                        };
                        let r = call_archaic_segments(&calls, &og, &cls, &cal,
                            &GeneticMap::from_markers(Vec::new()), &cfg);
                        let medlen = { let mut v: Vec<f64> = r.segments.iter().map(|s| s.length_mb()).collect();
                                       v.sort_by(|a,b| a.total_cmp(b));
                                       if v.is_empty() {0.0} else {v[v.len()/2]} };
                        println!("{sw:>8.1} {mult:>6.1} {post:>6.2} {minkb:>8} {:>9} {:>7.2} {medlen:>8.3}",
                                 r.summary.n_segments, r.summary.total_mb);
                    }
                }
            }
        }
        return Ok(());
    }
    let r = call_archaic_segments(&calls, &og, &cls, &cal, &GeneticMap::from_markers(Vec::new()), &ArchaicConfig::default());
    let s = &r.summary;
    println!("segments {}  total {:.2} Mb  = {:.2}% of {:.1} Mb callable", s.n_segments, s.total_mb, s.pct_callable, s.callable_mb);
    println!("  Neanderthal {:.2} Mb   Denisovan {:.2} Mb   Unknown {:.2} Mb", s.neanderthal_mb, s.denisovan_mb, s.unknown_mb);
    for seg in r.segments.iter().take(6) {
        println!("  {} {}-{} ({:.2} Mb) post {:.2} private {} ({:.0}/Mb) {:?} nea{} den{}",
                 seg.contig, seg.start, seg.end, seg.length_mb(), seg.posterior, seg.n_private,
                 seg.n_private as f64 / seg.length_mb().max(1e-9), seg.source,
                 seg.neanderthal_matches, seg.denisovan_matches);
    }
    Ok(())
}

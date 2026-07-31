//! Run the reference-based archaic tract caller ([`archaic_match`]) on real cached calls.
//!
//! The unit tests prove the model behaves on synthetic runs; this is what shows whether it finds
//! REAL tracts. Emits segments as JSON for `scripts/archaic-validation/compare_locations.py`, which
//! scores them against an external callset and — critically — against the random-placement null the
//! density caller failed.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_match_probe -- \
//!   classify.bin callable.bin chm13v2.0.fa genetic_map.bin calls.chr21.json calls.chr22.json \
//!   > segments.json
//! ```

use std::collections::BTreeMap;

use navigator_analysis::archaic::{ArchaicCallable, ArchaicClassify};
use navigator_analysis::archaic_match::{call_from_observations, observations_for_contig, MatchConfig, SiteObs};
use navigator_analysis::caller::SiteGenotype;
use navigator_analysis::ibd::GeneticMap;
use navigator_analysis::reader::read_contig_sequence;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 5 {
        eprintln!(
            "usage: archaic_match_probe <classify.bin> <callable.bin> <reference.fa> \
             <genetic_map.bin|-> <calls.json> [calls.json ...]"
        );
        std::process::exit(2);
    }
    let classify = ArchaicClassify::from_bytes(&std::fs::read(&a[0])?).map_err(|e| e.to_string())?;
    let callable = ArchaicCallable::from_bytes(&std::fs::read(&a[1])?).map_err(|e| e.to_string())?;
    let reference = std::path::PathBuf::from(&a[2]);

    let mut calls: Vec<SiteGenotype> = Vec::new();
    for p in &a[4..] {
        let mut v: Vec<SiteGenotype> = serde_json::from_str(&std::fs::read_to_string(p)?)?;
        calls.append(&mut v);
    }
    let mut by_contig: BTreeMap<String, BTreeMap<i64, &SiteGenotype>> = BTreeMap::new();
    for c in &calls {
        by_contig.entry(c.contig.clone()).or_default().insert(c.position, c);
    }
    eprintln!("{} calls over {} contig(s)", calls.len(), by_contig.len());

    let mut observations: BTreeMap<String, Vec<SiteObs>> = BTreeMap::new();
    let mut lengths: Vec<(String, i32)> = Vec::new();
    for (contig, pos_map) in &by_contig {
        // The reference base decides which diagnostic sites are informative at all, so it is read
        // rather than assumed (see `observations_for_contig`).
        let seq = read_contig_sequence(&reference, contig)?;
        let obs = observations_for_contig(
            contig,
            &classify,
            pos_map,
            |p| seq.get((p - 1).max(0) as usize).copied().map(|b| b.to_ascii_uppercase()),
            &callable,
            0.5,
        );
        let carried = obs.iter().filter(|o| o.carries).count();
        eprintln!(
            "{contig}: {} informative diagnostic sites, {carried} carried ({:.1}%)",
            obs.len(),
            if obs.is_empty() { 0.0 } else { carried as f64 * 100.0 / obs.len() as f64 }
        );
        lengths.push((contig.clone(), seq.len() as i32));
        observations.insert(contig.clone(), obs);
    }

    let gmap = if a[3] == "-" {
        let pairs: Vec<(&str, i32)> = lengths.iter().map(|(c, l)| (c.as_str(), *l)).collect();
        eprintln!("genetic map: uniform 1 cM/Mb");
        GeneticMap::uniform(1.0, &pairs)
    } else {
        GeneticMap::from_bytes(&std::fs::read(&a[3])?).map_err(|e| e.to_string())?
    };

    // `ARCHAIC_RATIOS=2.0,2.5,3.04` sweeps the emission ratio in one process. It cannot be swept
    // post-hoc like the three thresholds — it changes the emissions, so the HMM must be re-decoded —
    // but the expensive part (reading the reference, walking the diagnostic sites) is per sample,
    // not per ratio, so doing it here costs one pass instead of one per value.
    if let Ok(spec) = std::env::var("ARCHAIC_RATIOS") {
        let mut out = serde_json::Map::new();
        for tok in spec.split(',').filter(|t| !t.trim().is_empty()) {
            let ratio: f64 = tok.trim().parse()?;
            let cfg = MatchConfig {
                archaic_ratio: ratio,
                ..Default::default()
            };
            let r = call_from_observations(&observations, &gmap, &callable, &cfg);
            eprintln!(
                "  ratio {ratio:>5.2} -> {} segments, {:.3} Mb",
                r.summary.n_segments, r.summary.total_mb
            );
            out.insert(tok.trim().to_string(), serde_json::to_value(&r)?);
        }
        println!("{}", serde_json::to_string(&out)?);
        return Ok(());
    }

    let result = call_from_observations(&observations, &gmap, &callable, &MatchConfig::default());
    let s = &result.summary;
    eprintln!(
        "SEGMENTS {}  total {:.3} Mb  = {:.2}% of {:.1} Mb callable",
        s.n_segments, s.total_mb, s.pct_callable, s.callable_mb
    );
    println!("{}", serde_json::to_string(&result)?);
    Ok(())
}

//! Per-window counts of African-outgroup segregating sites — a candidate local mutation-rate proxy.
//!
//! The Tier B emission model assumes one background rate genome-wide. Measured, the background
//! private-variant density varies 5.3x between its 10th and 90th percentile and is 14.6x
//! overdispersed relative to the Poisson it is modelled with, which is larger than the 2.89x
//! enrichment inside real archaic tracts — so the model calls its own upper tail archaic. hmmix
//! avoids this with a mutation-rate map; we have no such asset.
//!
//! The density of sites segregating in Africans is already in `archaic_outgroup_af_<build>.bin` and
//! is a direct measure of how variable a region is, for reasons that have nothing to do with
//! archaic introgression (mutation rate, reference quality, mappability). This dumps it so that
//! proxy can be tested as a normalizer before an asset is built for the purpose.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_outgroup_density -- \
//!   ~/.decodingus/ancestry/archaic_outgroup_af_chm13v2.0.bin 1000 chr21 chr22 > og_density.tsv
//! ```

use navigator_analysis::archaic::ArchaicOutgroup;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a.next().expect("usage: archaic_outgroup_density <outgroup.bin> [window_bp] [contig ...]");
    let window: i64 = a.next().and_then(|s| s.parse().ok()).unwrap_or(1000);
    let want: Vec<String> = a.collect();

    let og = ArchaicOutgroup::from_bytes(&std::fs::read(&path)?).map_err(|e| e.to_string())?;
    println!("contig\twindow_start\tn_outgroup_sites");
    for c in &og.contigs {
        if !want.is_empty() && !want.contains(&c.contig) {
            continue;
        }
        let mut counts: std::collections::BTreeMap<i64, u32> = Default::default();
        let mut n = 0u64;
        for p in c.iter() {
            *counts.entry(p / window * window).or_insert(0) += 1;
            n += 1;
        }
        eprintln!("{}: {n} outgroup sites in {} non-empty windows", c.contig, counts.len());
        for (w, k) in counts {
            println!("{}\t{}\t{}", c.contig, w, k);
        }
    }
    Ok(())
}

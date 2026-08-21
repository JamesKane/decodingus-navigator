//! The count of sites that vary in the African outgroup, in each window. It is a candidate proxy
//! for the local mutation rate.
//!
//! The Tier B emission model takes one background rate over the whole genome. A measurement shows
//! otherwise. The background density of private variants changes by 5.3x between its 10th and 90th
//! percentile. It is also 14.6x more spread out than the Poisson distribution that the model gives
//! it.
//!
//! That spread is larger than the 2.89x enrichment inside a real archaic tract. So the model calls
//! its own upper tail archaic. The hmmix tool avoids this with a map of the mutation rate, and
//! this project has no such asset.
//!
//! `archaic_outgroup_af_<build>.bin` already holds the density of the sites that vary in Africans.
//! That density measures directly how much a region varies. Its reasons have nothing to do with
//! archaic introgression. They are the mutation rate, the quality of the reference, and how well
//! reads map there. This tool writes it out, so that somebody can test it as a normalizer before anybody
//! builds an asset for the purpose.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_outgroup_density -- \
//!   ~/.decodingus/ancestry/archaic_outgroup_af_chm13v2.0.bin 1000 chr21 chr22 > og_density.tsv
//! ```

use navigator_analysis::archaic::ArchaicOutgroup;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a
        .next()
        .expect("usage: archaic_outgroup_density <outgroup.bin> [window_bp] [contig ...]");
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

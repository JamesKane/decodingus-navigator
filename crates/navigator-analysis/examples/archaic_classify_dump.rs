//! Dump the archaic **diagnostic** sites (position, derived base, lineage class) as TSV.
//!
//! Written to test a different observable for the Tier B HMM. The current model counts *all*
//! private variants per window, and that signal is weak: measured on a real European, archaic
//! tracts carry only 2.89x the background density while the background itself varies 5.3x between
//! its 10th and 90th percentile. Restricting the observable to sites where the derived allele is
//! actually known to be archaic should be far more specific.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_classify_dump -- \
//!   ~/.decodingus/ancestry/archaic_classify_chm13v2.0.bin chr21 chr22 > classify.tsv
//! ```

use navigator_analysis::archaic::ArchaicClassify;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a.next().expect("usage: archaic_classify_dump <classify.bin> [contig ...]");
    let want: Vec<String> = a.collect();

    let cls = ArchaicClassify::from_bytes(&std::fs::read(&path)?).map_err(|e| e.to_string())?;
    println!("contig\tposition\tderived\tclass");
    for c in &cls.contigs {
        let name = &c.positions.contig;
        if !want.is_empty() && !want.contains(name) {
            continue;
        }
        let mut n = 0usize;
        for (i, p) in c.positions.iter().enumerate() {
            let d = c.derived.get(i).copied().unwrap_or(b'N') as char;
            let k = c.classes.get(i).copied().unwrap_or(2);
            println!("{name}\t{p}\t{d}\t{k}");
            n += 1;
        }
        eprintln!("{name}: {n} diagnostic sites");
    }
    Ok(())
}

//! Write out the Tier A marker panel as a TSV, with the calls of **each archaic genome**.
//!
//! This is the independent evidence that decides a Tier B call. The segment caller
//! ([`navigator_analysis::archaic_match`]) reads `ArchaicClassify` alone, which holds a derived
//! base and a lineage class at each site. It never sees which archaic genome carries what. The
//! pattern over the genomes is information that nobody could have fitted the caller to, and that is
//! what makes it a referee.
//!
//! Here is why a referee is necessary. A measurement gave the precision against the callset of
//! hmmix. But a call that hmmix does not hold is not wrong by that fact alone. The own tracts of
//! hmmix show an enrichment of only 1.84x for their own archaic SNPs. So that callset is
//! incomplete by an amount that nobody knows.
//!
//! A score of a segment against the archaic genomes asks a direct question. Does this segment look
//! like an archaic haplotype that came down from an ancestor? It asks no other caller for an
//! opinion.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_panel_dump -- \
//!   ~/.decodingus/ancestry/archaic_markers_chm13v2.0.bin chr21 chr22 > panel.tsv
//! ```

use navigator_analysis::archaic::{ArchaicMarkerPanel, ARCHAIC_GENOMES};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a
        .next()
        .expect("usage: archaic_panel_dump <archaic_markers.bin> [contig ...]");
    let want: Vec<String> = a.collect();

    let panel = ArchaicMarkerPanel::from_bytes(&std::fs::read(&path)?).map_err(|e| e.to_string())?;
    eprintln!("panel: {} sites, build {}", panel.sites.len(), panel.build);

    // One column for each archaic genome. D means that the genome carries the derived allele. A
    // means that the caller positively called it homozygous-ancestral. A `.` means no call.
    //
    // The difference between A and `.` carries weight. To read a no-call as ancestral is the error
    // that gave about 19% Denisovan for a European, in an earlier pass.
    println!("contig\tposition\tderived\tclass\t{}", ARCHAIC_GENOMES.join("\t"));
    let mut n = 0usize;
    for s in &panel.sites {
        if !want.is_empty() && !want.contains(&s.contig) {
            continue;
        }
        let calls: Vec<&str> = s
            .calls
            .iter()
            .map(|c| {
                if c.carries_derived() {
                    "D"
                } else if matches!(c, navigator_analysis::archaic::ArchaicCall::HomAncestral) {
                    "A"
                } else {
                    "."
                }
            })
            .collect();
        println!(
            "{}\t{}\t{}\t{:?}\t{}",
            s.contig,
            s.position,
            s.archaic_derived_allele,
            s.diagnostic_class,
            calls.join("\t")
        );
        n += 1;
    }
    eprintln!("wrote {n} sites");
    Ok(())
}

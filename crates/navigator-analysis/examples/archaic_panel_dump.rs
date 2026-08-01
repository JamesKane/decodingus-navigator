//! Dump the Tier A marker panel with its **per-archaic-genome** calls, as TSV.
//!
//! This is the independent evidence for arbitrating Tier B calls. The segment caller
//! ([`navigator_analysis::archaic_match`]) reads only `ArchaicClassify` — a derived base and a
//! lineage class per site — and never sees which archaic genome carries what. So the per-genome
//! pattern is information the caller cannot have fitted to, which is what makes it usable as a
//! referee.
//!
//! Why a referee is needed: precision has been measured against hmmix's callset, but a call absent
//! from hmmix is not necessarily wrong — hmmix's own tracts are enriched only 1.84x for their own
//! archaic SNPs, so that callset is incomplete by an unknown amount. Scoring a segment against the
//! archaic genomes directly asks whether it looks like an inherited archaic haplotype, without
//! asking another caller's opinion.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_panel_dump -- \
//!   ~/.decodingus/ancestry/archaic_markers_chm13v2.0.bin chr21 chr22 > panel.tsv
//! ```

use navigator_analysis::archaic::{ArchaicMarkerPanel, ARCHAIC_GENOMES};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a.next().expect("usage: archaic_panel_dump <archaic_markers.bin> [contig ...]");
    let want: Vec<String> = a.collect();

    let panel = ArchaicMarkerPanel::from_bytes(&std::fs::read(&path)?).map_err(|e| e.to_string())?;
    eprintln!("panel: {} sites, build {}", panel.sites.len(), panel.build);

    // One column per archaic genome: D = carries the derived allele, A = positively called
    // homozygous-ancestral, . = no call. The A/. distinction is load-bearing — treating a no-call as
    // ancestral is the error that produced ~19 % Denisovan for a European in an earlier pass.
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

//! Genotype an alignment at every site in a 1240k-style BED (col4 = `rsid|ref|alt`) and dump
//! `rsid<TAB>contig<TAB>pos<TAB>dosage` (dosage 0/1/2, -1 = no-call). Used to genotype a target at
//! the full 1240k for the qpAdm rebuild (docs/design/ancient-ancestry-rebuild.md §7.11).
//!   genotype_bed <sites.bed> <bam_or_cram> <out.tsv> [reference.fa]
use navigator_analysis::caller::{genotype_sites_all_contigs, HaploidCallerParams, Site};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let bed = std::env::args().nth(1).expect("usage: genotype_bed <sites.bed> <bam> <out.tsv> [ref.fa]");
    let bam = PathBuf::from(std::env::args().nth(2).expect("bam"));
    let out = std::env::args().nth(3).expect("out.tsv");
    let reference = std::env::args().nth(4).map(PathBuf::from);

    let mut sites = Vec::new();
    let mut rsids = Vec::new();
    for line in std::fs::read_to_string(&bed)?.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let contig = f[0].to_string();
        let pos: i64 = f[2].parse()?; // BED end = 1-based position
        let name: Vec<&str> = f[3].split('|').collect();
        let (rsid, r, a) = (name[0], name.get(1).copied().unwrap_or("N"), name.get(2).copied().unwrap_or("N"));
        rsids.push(rsid.to_string());
        sites.push(Site {
            name: rsid.to_string(),
            contig,
            position: pos,
            reference_allele: r.to_string(),
            alternate_allele: a.to_string(),
        });
    }
    eprintln!("genotyping {} sites from {} ...", sites.len(), bam.display());
    let params = HaploidCallerParams::default();
    let gts = genotype_sites_all_contigs(&bam, &sites, 2, &params, reference.as_deref()).map_err(|e| anyhow::anyhow!("{e}"))?;

    // genotype_sites_all_contigs returns genotypes REORDERED (per-contig), so we must key each
    // returned genotype to its rsID by (contig,position) — NOT by input order. Zipping with `rsids`
    // would mislabel every genotype (only ~14% would land on the right rsID).
    let rsid_at: std::collections::HashMap<(&str, i64), &str> = sites
        .iter()
        .zip(&rsids)
        .map(|(s, r)| ((s.contig.as_str(), s.position), r.as_str()))
        .collect();
    let mut w = BufWriter::new(std::fs::File::create(&out)?);
    let mut called = 0usize;
    for g in &gts {
        if g.dosage >= 0 {
            called += 1;
        }
        let rsid = rsid_at.get(&(g.contig.as_str(), g.position)).copied().unwrap_or(".");
        writeln!(w, "{}\t{}\t{}\t{}", rsid, g.contig, g.position, g.dosage)?;
    }
    w.flush()?;
    eprintln!("wrote {out}: {} of {} sites called", called, gts.len());
    Ok(())
}

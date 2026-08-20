//! A throwaway tool that tests the Option A′ concept. It writes a copy of the ancient AF panel
//! that holds the sites a consumer array reaches, and no other.
//!
//! It reads the rsIDs that one consumer chip file, or more, assays. It maps each one to a CHM13
//! (contig,pos) pair through the IBD panel, which carries an rsid with its CHM13 locus. It then
//! keeps the ancient-panel sites at those positions alone. Usage:
//!   ascertain_chip <ancient.bin> <ibd_panel.bin> <out.bin> <chip1.txt> [chip2.txt ...]
use navigator_analysis::ancestry::AncestryPanel;
use navigator_analysis::ibd_panel::IbdPanel;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let ancient_path = args
        .next()
        .expect("usage: ascertain_chip <ancient.bin> <ibd.bin> <out.bin> <chip...>");
    let ibd_path = args.next().expect("ibd.bin");
    let out_path = args.next().expect("out.bin");
    let chip_files: Vec<String> = args.collect();
    anyhow::ensure!(!chip_files.is_empty(), "need at least one chip file");

    // rsIDs the consumer arrays assay (col 1 of each chip file; keep only rs-prefixed).
    let mut chip_rsids: HashSet<String> = HashSet::new();
    for f in &chip_files {
        let r = BufReader::new(std::fs::File::open(f)?);
        for line in r.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            if let Some(tok) = line.split(['\t', ',', ' ']).next() {
                if tok.starts_with("rs") {
                    chip_rsids.insert(tok.to_string());
                }
            }
        }
    }

    // Map each assayed rsID to a CHM13 (contig,pos) pair, through the IBD panel.
    let ibd = IbdPanel::from_bytes(&std::fs::read(&ibd_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let chip_loci: HashSet<(String, i64)> = ibd
        .sites
        .iter()
        .filter(|s| chip_rsids.contains(&s.rsid))
        .map(|s| (s.chm13.contig.clone(), s.chm13.position))
        .collect();

    let mut ancient = AncestryPanel::from_bytes(&std::fs::read(&ancient_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let before = ancient.sites.len();
    ancient
        .sites
        .retain(|s| chip_loci.contains(&(s.contig.clone(), s.position)));
    let after = ancient.sites.len();
    std::fs::write(&out_path, ancient.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;

    println!(
        "chip rsIDs: {} | IBD sites on those arrays: {} | ancient sites kept: {}/{} ({:.1}%) -> {}",
        chip_rsids.len(),
        chip_loci.len(),
        after,
        before,
        100.0 * after as f64 / before as f64,
        out_path
    );
    Ok(())
}

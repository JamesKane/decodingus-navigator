//! Diagnostic: overlap + allele orientation between two AncestryPanels (by contig,pos).
//!   panel_overlap <a.bin> <b.bin>
use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let a = AncestryPanel::from_bytes(&std::fs::read(std::env::args().nth(1).unwrap())?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let b = AncestryPanel::from_bytes(&std::fs::read(std::env::args().nth(2).unwrap())?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let bm: HashMap<(String, i64), (char, char)> = b
        .sites
        .iter()
        .map(|s| ((s.contig.clone(), s.position), (s.reference_allele, s.alternate_allele)))
        .collect();
    let (mut overlap, mut same, mut swapped, mut other) = (0, 0, 0, 0);
    for s in &a.sites {
        if let Some(&(r, al)) = bm.get(&(s.contig.clone(), s.position)) {
            overlap += 1;
            if (s.reference_allele, s.alternate_allele) == (r, al) {
                same += 1;
            } else if (s.reference_allele, s.alternate_allele) == (al, r) {
                swapped += 1;
            } else {
                other += 1;
            }
        }
    }
    println!("A {} sites, B {} sites", a.sites.len(), b.sites.len());
    println!("overlap {overlap}: same-orientation {same}, swapped {swapped}, other {other}");
    Ok(())
}

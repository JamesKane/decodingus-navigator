//! Throwaway: emit a copy of the ancient AF panel keeping only TRANSVERSION sites (ref/alt not a
//! transition A<->G or C<->T). aDNA post-mortem damage (cytosine deamination) corrupts transitions;
//! if restricting to transversions collapses the WGS/chip disagreement, damage is the cause.
use navigator_analysis::ancestry::AncestryPanel;

fn is_transition(r: char, a: char) -> bool {
    matches!((r, a), ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C'))
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: filter_tv <ancient.bin> <out.bin>");
    let out = std::env::args().nth(2).expect("out.bin");
    let mut panel = AncestryPanel::from_bytes(&std::fs::read(&path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let before = panel.sites.len();
    panel.sites.retain(|s| !is_transition(s.reference_allele, s.alternate_allele));
    let after = panel.sites.len();
    std::fs::write(&out, panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    println!("transversions only: kept {after}/{before} sites ({:.1}%) -> {out}", 100.0 * after as f64 / before as f64);
    Ok(())
}

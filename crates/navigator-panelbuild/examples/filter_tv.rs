//! A throwaway tool. It writes a copy of the ancient AF panel that holds the TRANSVERSION sites
//! alone. At such a site the ref and alt are not a transition, which is A to G, or C to T.
//!
//! Post-mortem damage in aDNA, which is cytosine deamination, corrupts a transition. So if the
//! transversions alone close the disagreement between the WGS answer and the chip answer, that
//! damage is the cause.
use navigator_analysis::ancestry::AncestryPanel;

fn is_transition(r: char, a: char) -> bool {
    matches!((r, a), ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C'))
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: filter_tv <ancient.bin> <out.bin>");
    let out = std::env::args().nth(2).expect("out.bin");
    let mut panel = AncestryPanel::from_bytes(&std::fs::read(&path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let before = panel.sites.len();
    panel
        .sites
        .retain(|s| !is_transition(s.reference_allele, s.alternate_allele));
    let after = panel.sites.len();
    std::fs::write(&out, panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    println!(
        "transversions only: kept {after}/{before} sites ({:.1}%) -> {out}",
        100.0 * after as f64 / before as f64
    );
    Ok(())
}

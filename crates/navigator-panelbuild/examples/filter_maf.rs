//! A throwaway tool. It writes a copy of the ancient AF panel that holds only the sites whose
//! modern EUR MAF reaches a threshold. That MAF comes from the super-pop panel, joined by contig
//! and pos.
//!
//! It checks the ascertainment-floor fix, and nobody has to build the panel again from the AADR.
//! Point NAVIGATOR_ANCESTRY_FREQ_ANCIENT at the output, and run the stability gate again.
use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let a = std::env::args()
        .nth(1)
        .expect("usage: filter_maf <ancient.bin> <super.bin> <min_maf> <out.bin>");
    let s = std::env::args().nth(2).expect("super.bin");
    let t: f64 = std::env::args().nth(3).expect("min_maf").parse()?;
    let out = std::env::args().nth(4).expect("out.bin");

    let mut ancient = AncestryPanel::from_bytes(&std::fs::read(&a)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let sup = AncestryPanel::from_bytes(&std::fs::read(&s)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let eur = sup.populations.iter().position(|p| p == "EUR").expect("no EUR");
    let maf: HashMap<(String, i64), f64> = sup
        .sites
        .iter()
        .filter_map(|x| {
            x.freqs
                .get(eur)
                .map(|&f| ((x.contig.clone(), x.position), (f as f64).min(1.0 - f as f64)))
        })
        .collect();

    let before = ancient.sites.len();
    ancient
        .sites
        .retain(|x| maf.get(&(x.contig.clone(), x.position)).is_some_and(|&m| m >= t));
    let after = ancient.sites.len();
    std::fs::write(&out, ancient.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    println!("min_maf={t}: kept {after}/{before} sites -> {out}");
    Ok(())
}

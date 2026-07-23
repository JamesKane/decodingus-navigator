//! Print the SUPER-POP admixture (the prior the chromosome painter uses) for a dosage TSV against
//! the super-pop panel — to check the painting's global-composition gate threshold.
//!   score_superpop_from_tsv <ancestry_panel.bin> <dosage.tsv>
use navigator_analysis::ancestry::{estimate_admixture, AncestryPanel};
use navigator_analysis::caller::SiteGenotype;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let panel_path = std::env::args().nth(1).expect("usage: score_superpop_from_tsv <panel.bin> <dosage.tsv>");
    let tsv = std::env::args().nth(2).expect("dosage.tsv");
    let panel = AncestryPanel::from_bytes(&std::fs::read(&panel_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    let mut dosage: HashMap<(String, i64), i32> = HashMap::new();
    for line in std::fs::read_to_string(&tsv)?.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (c, p, d) = match f.len() {
            n if n >= 4 => (f[1], f[2], f[3]),
            3 => (f[0], f[1], f[2]),
            _ => continue,
        };
        let (Ok(p), Ok(d)) = (p.trim().parse::<i64>(), d.trim().parse::<i32>()) else { continue };
        dosage.insert((c.to_string(), p), d);
    }
    let gts: Vec<SiteGenotype> = panel
        .sites
        .iter()
        .filter_map(|s| {
            dosage.get(&(s.contig.clone(), s.position)).map(|&d| SiteGenotype {
                name: String::new(),
                contig: s.contig.clone(),
                position: s.position,
                reference_allele: s.reference_allele.to_string(),
                alternate_allele: s.alternate_allele.to_string(),
                ploidy: 2,
                dosage: d,
                gq: 40,
                depth: 20,
                ref_depth: 10,
                alt_depth: 10,
                pls: vec![],
                gt: None,
                allele_depths: None,
            })
        })
        .collect();
    let r = estimate_admixture(&gts, &panel, "chm13v2.0");
    println!("super-pop panel: {} pops, {} sites used", panel.populations.len(), r.snps_with_genotype);
    let mut comps = r.components.clone();
    comps.sort_by(|a, b| b.percentage.total_cmp(&a.percentage));
    for c in &comps {
        let gated = if c.percentage < 2.0 { "  <- dropped by 2% gate" } else { "" };
        println!("  {:<5} {:>6.2} %{}", c.population_code, c.percentage, gated);
    }
    Ok(())
}

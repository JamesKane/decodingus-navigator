//! qpAdm from an external caller's genotypes (docs/design/ancient-ancestry-rebuild.md §7): is the
//! WGS batch effect our native caller, or the method? Reads a `contig<TAB>pos<TAB>dosage` TSV
//! (dosage 0/1/2, -1 = no-call) — e.g. from a GATK4 VCF via
//! `bcftools query -f '%CHROM\t%POS\t[%GT]\n'` mapped to dosages — matches it to the panel by
//! (contig,pos), and runs qpadm_fit. Compare the weights to our caller's qpadm_check result on the
//! same alignment.
//!   qpadm_from_tsv <qpadm_panel.bin> <dosage.tsv>
use navigator_analysis::ancestry::{qpadm_fit, AncestryPanel, F4_BLOCK_BP};
use navigator_analysis::caller::SiteGenotype;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let panel_path = std::env::args().nth(1).expect("usage: qpadm_from_tsv <panel.bin> <dosage.tsv>");
    let tsv = std::env::args().nth(2).expect("dosage.tsv");
    let panel = AncestryPanel::from_bytes(&std::fs::read(&panel_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    // (contig,pos) -> dosage.
    let mut dosage: HashMap<(String, i64), i32> = HashMap::new();
    for line in std::fs::read_to_string(&tsv)?.lines() {
        let mut it = line.split('\t');
        let (Some(c), Some(p), Some(d)) = (it.next(), it.next(), it.next()) else { continue };
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
    let called = gts.iter().filter(|g| g.dosage >= 0).count();
    eprintln!("{} panel sites, {} matched in TSV, {called} called", panel.sites.len(), gts.len());

    let src_codes = ["WHG", "ANF", "Steppe"];
    let sources: Vec<usize> = src_codes.iter().map(|c| panel.populations.iter().position(|p| p == c).unwrap()).collect();
    let outgroups: Vec<usize> = (0..panel.populations.len()).filter(|i| !sources.contains(i)).collect();

    let fit = qpadm_fit(&gts, &panel, &sources, &outgroups, F4_BLOCK_BP)
        .ok_or_else(|| anyhow::anyhow!("qpadm_fit returned None"))?;
    println!("\nsites {}  blocks {}  dof {}  chi2 {:.2}  p {:.4}", fit.n_sites, fit.n_blocks, fit.dof, fit.chi2, fit.p_value);
    for (code, i) in src_codes.iter().zip(0..) {
        println!("  {code:<8} {:>6.1} %   (SE {:.1})", fit.weights[i] * 100.0, fit.std_errors[i] * 100.0);
    }
    println!(
        "\nmodel {} at p=0.05; weights {}",
        if fit.p_value >= 0.05 { "ACCEPTED" } else { "REJECTED" },
        if fit.weights_feasible(0.02) { "feasible" } else { "INFEASIBLE" }
    );
    Ok(())
}

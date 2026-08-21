//! Score a candidate MODERN ancestry asset, which is a fine-admixture panel and an optional PCA,
//! against a subject's CHM13-oriented dosages. This is the offline check harness for the sweep over
//! panel depth (documents/design/ancient-ancestry-rebuild.md). It is the companion of
//! `qpadm_from_tsv`.
//!
//! It reads a dosage TSV in the form `[rsid\t]contig\tpos\tdosage`, where a dosage is 0, 1, or 2,
//! and -1 is a no-call. The `genotype_bed` example, or a chip resolver, can produce one. It matches
//! that file to each candidate asset by (contig,pos), runs `estimate_fine_admixture`, and prints
//! the breakdown with the count of sites that it used. That count is the quantity that the sweep
//! must grow.
//!
//! Point it at a 20k, 100k, or 200k candidate. That shows whether the extra depth changes the
//! estimate, and where it stops. Run it on the same subject's WGS dosages and chip dosages, to check
//! that the two agree.
//!   score_modern_from_tsv <fine_panel.bin> <dosage.tsv> [pca.bin]
use navigator_analysis::ancestry::{estimate_fine_admixture, project_pca, AncestryPanel, PcaLoadings};
use navigator_analysis::caller::SiteGenotype;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let fine_path = std::env::args()
        .nth(1)
        .expect("usage: score_modern_from_tsv <fine_panel.bin> <dosage.tsv> [pca.bin]");
    let tsv = std::env::args().nth(2).expect("dosage.tsv");
    let pca_path = std::env::args().nth(3).filter(|s| !s.is_empty());

    let fine = AncestryPanel::from_bytes(&std::fs::read(&fine_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    // (contig,pos) -> dosage. Accepts genotype_bed's 4-col `rsid contig pos dosage` and a plain
    // 3-col `contig pos dosage`.
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
        let (Ok(p), Ok(d)) = (p.trim().parse::<i64>(), d.trim().parse::<i32>()) else {
            continue;
        };
        dosage.insert((c.to_string(), p), d);
    }

    // Build the genotypes over the fine panel's sites. `estimate_admixture` keys by (contig,pos),
    // and reads the alt frequency of each population from the panel. So the caller must orient a
    // dosage to CHM13 and to the panel first. The `genotype_bed` and `resolve_chip` producers
    // promise that.
    let gts: Vec<SiteGenotype> = fine
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

    let result = estimate_fine_admixture(&gts, &fine, "chm13v2.0");
    println!(
        "fine panel {} sites | {} matched in TSV | {} sites used (called & informative)",
        fine.sites.len(),
        gts.len(),
        result.snps_with_genotype,
    );

    println!("\nsuper-population roll-up:");
    let mut sup = result.super_population_summary.clone();
    sup.sort_by(|a, b| b.percentage.total_cmp(&a.percentage));
    for s in sup.iter().filter(|s| s.percentage >= 0.05) {
        println!("  {:<6} {:>6.1} %", s.super_population, s.percentage);
    }

    println!("\nfine populations (>= 0.5%):");
    let mut comps = result.components.clone();
    comps.sort_by(|a, b| b.percentage.total_cmp(&a.percentage));
    for c in comps.iter().filter(|c| c.percentage >= 0.5) {
        println!(
            "  {:<5} {:<22} {:>6.1} %",
            c.population_code, c.population_name, c.percentage
        );
    }

    if let Some(pca_path) = pca_path {
        let pca = PcaLoadings::from_bytes(&std::fs::read(&pca_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
        let coords = project_pca(&gts, &pca);
        let shown: Vec<String> = coords.iter().take(4).map(|c| format!("{c:.2}")).collect();
        println!("\nPCA ({} sites): PC[1..4] = [{}]", pca.sites.len(), shown.join(", "));
    }
    Ok(())
}

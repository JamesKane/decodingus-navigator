//! Lever-2 viability probe (docs/design/ancient-ancestry-rebuild.md §7): does the qpAdm f4 estimator
//! pull a WGS sample out of the fabricated ~80% Steppe band into the sane NW-European range
//! (Steppe 40–55 / ANF 25–40 / WHG 10–25), and does the model-fit p-value behave?
//!
//! Genotypes one alignment at the qpAdm panel sites (sources WHG/ANF/Steppe + outgroups), then runs
//! qpadm_fit. Sources = the panel populations named WHG/ANF/Steppe; every other population is an
//! outgroup. Compare across sources (WGS vs chip) for the stability gate.
//!   qpadm_check <qpadm_panel.bin> <bam_or_cram> [reference.fa]
use navigator_analysis::ancestry::{qpadm_fit, AncestryPanel, F4_BLOCK_BP};
use navigator_analysis::caller::{genotype_sites_all_contigs, HaploidCallerParams, Site};
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let panel_path = std::env::args().nth(1).expect("usage: qpadm_check <qpadm_panel.bin> <bam> [ref.fa]");
    let bam = PathBuf::from(std::env::args().nth(2).expect("bam"));
    let reference = std::env::args().nth(3).map(PathBuf::from);
    let panel = AncestryPanel::from_bytes(&std::fs::read(&panel_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    // Sources = WHG/ANF/Steppe (in that order); outgroups = every other population.
    let src_codes = ["WHG", "EEF", "Steppe"];
    let sources: Vec<usize> = src_codes
        .iter()
        .map(|c| panel.populations.iter().position(|p| p == c).unwrap_or_else(|| panic!("panel missing source {c}")))
        .collect();
    let outgroups: Vec<usize> = (0..panel.populations.len()).filter(|i| !sources.contains(i)).collect();
    eprintln!(
        "panel: {} sites, sources {:?}, outgroups {:?}",
        panel.sites.len(),
        src_codes,
        outgroups.iter().map(|&i| panel.populations[i].as_str()).collect::<Vec<_>>()
    );

    let sites: Vec<Site> = panel
        .sites
        .iter()
        .map(|s| Site {
            name: format!("{}:{}", s.contig, s.position),
            contig: s.contig.clone(),
            position: s.position,
            reference_allele: s.reference_allele.to_string(),
            alternate_allele: s.alternate_allele.to_string(),
        })
        .collect();

    eprintln!("genotyping {} sites from {} ...", sites.len(), bam.display());
    let params = HaploidCallerParams::default();
    let gts = genotype_sites_all_contigs(&bam, &sites, 2, &params, reference.as_deref()).map_err(|e| anyhow::anyhow!("{e}"))?;
    let called = gts.iter().filter(|g| g.dosage >= 0).count();
    eprintln!("genotyped: {called} of {} sites called", gts.len());

    let fit = qpadm_fit(&gts, &panel, &sources, &outgroups, F4_BLOCK_BP)
        .ok_or_else(|| anyhow::anyhow!("qpadm_fit returned None (too few sites/blocks or singular system)"))?;

    println!("\n(reference: chip ~58% Steppe · old frequency-EM on WGS ~80% Steppe · NW-Eur band 40–55)\n");
    println!("sites {}  blocks {}  dof {}  chi2 {:.2}  p {:.4}", fit.n_sites, fit.n_blocks, fit.dof, fit.chi2, fit.p_value);
    for (code, i) in src_codes.iter().zip(0..) {
        println!("  {code:<8} {:>6.1} %   (SE {:.1})", fit.weights[i] * 100.0, fit.std_errors[i] * 100.0);
    }
    println!(
        "\nmodel {} at p=0.05; weights {}",
        if fit.p_value >= 0.05 { "ACCEPTED" } else { "REJECTED" },
        if fit.weights_feasible(0.02) { "feasible" } else { "INFEASIBLE (outside [0,1])" }
    );
    Ok(())
}

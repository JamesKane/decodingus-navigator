//! Lever-1 prototype (docs/design/ancient-ancestry-rebuild.md §5.1): does harmonizing the WGS target
//! to the pseudo-haploid ancient references collapse the WGS-vs-chip stability split?
//!
//! Genotypes one CHM13-native alignment at the ancient-panel sites (with allele depths), then fits
//! the ancient admixture four ways: diploid vs read-level pseudo-haploid (sample one allele per site
//! with P(alt)=alt_depth/depth), each over all sites vs transversions-only. Compare to the diploid
//! consensus (~75% Steppe) and the chip (~58%).
//!   phaploid_fit <ancient.bin> <bam_or_cram> [reference.fa]
use navigator_analysis::ancestry::{ancient_admixture_fit, AncestryPanel};
use navigator_analysis::caller::{genotype_sites_all_contigs, HaploidCallerParams, Site, SiteGenotype};
use std::path::PathBuf;

fn is_transition(r: &str, a: &str) -> bool {
    matches!((r, a), ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C"))
}

// Deterministic per-site draw of one allele, P(alt) = alt_depth/depth.
fn draw_alt(contig: &str, pos: i64, ref_d: u32, alt_d: u32) -> Option<bool> {
    let total = ref_d + alt_d;
    if total == 0 {
        return None;
    }
    let mut h = (pos as u64).wrapping_mul(0x9E3779B97F4A7C15);
    for b in contig.bytes() {
        h = (h ^ b as u64).wrapping_mul(0x100000001B3);
    }
    Some((h % total as u64) < alt_d as u64)
}

fn fit(label: &str, gts: &[SiteGenotype], panel: &AncestryPanel) {
    match ancient_admixture_fit(gts, panel, "chm13v2.0") {
        Some(r) => {
            let get = |c: &str| r.components.iter().find(|x| x.population_code == c).map_or(0.0, |x| x.percentage);
            println!(
                "{:<26} {:>6}   WHG {:>5.1}  ANF {:>5.1}  Steppe {:>5.1}   disp {:>5.2}",
                label,
                r.snps_with_genotype,
                get("WHG"),
                get("ANF"),
                get("Steppe"),
                r.fit_distance.unwrap_or(f64::NAN),
            );
        }
        None => println!("{label:<26} (fit returned None — too few sites)"),
    }
}

fn main() -> anyhow::Result<()> {
    let panel_path = std::env::args().nth(1).expect("usage: phaploid_fit <ancient.bin> <bam> [ref.fa]");
    let bam = PathBuf::from(std::env::args().nth(2).expect("bam"));
    let reference = std::env::args().nth(3).map(PathBuf::from);
    let panel = AncestryPanel::from_bytes(&std::fs::read(&panel_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

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

    eprintln!("genotyping {} ancient sites from {} ...", sites.len(), bam.display());
    let params = HaploidCallerParams::default();
    let gts = genotype_sites_all_contigs(&bam, &sites, 2, &params, reference.as_deref(), &navigator_analysis::CancelToken::none()).map_err(|e| anyhow::anyhow!("{e}"))?;
    let called = gts.iter().filter(|g| g.dosage >= 0).count();
    eprintln!("genotyped: {} sites called (of {})", called, gts.len());

    let tv = |g: &SiteGenotype| !is_transition(&g.reference_allele, &g.alternate_allele);

    // Diploid, as-genotyped.
    let dip: Vec<SiteGenotype> = gts.iter().filter(|g| g.dosage >= 0).cloned().collect();
    let dip_tv: Vec<SiteGenotype> = dip.iter().filter(|g| tv(g)).cloned().collect();

    // Read-level pseudo-haploid: one allele per site → homozygous representation (0 or 2).
    let phap: Vec<SiteGenotype> = gts
        .iter()
        .filter(|g| g.dosage >= 0)
        .filter_map(|g| {
            draw_alt(&g.contig, g.position, g.ref_depth, g.alt_depth).map(|alt| SiteGenotype {
                dosage: if alt { 2 } else { 0 },
                ..g.clone()
            })
        })
        .collect();
    let phap_tv: Vec<SiteGenotype> = phap.iter().filter(|g| tv(g)).cloned().collect();

    println!("\n(reference: diploid consensus ~75% Steppe · chip ~58% Steppe)\n");
    fit("diploid (all)", &dip, &panel);
    fit("diploid (transversions)", &dip_tv, &panel);
    fit("pseudo-haploid (all)", &phap, &panel);
    fit("pseudo-haploid (transv.)", &phap_tv, &panel);
    Ok(())
}

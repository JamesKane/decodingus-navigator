//! qpAdm control (docs/design/ancient-ancestry-rebuild.md §7): can the *real* panel + outgroups
//! recover a KNOWN European mixture drawn from the panel's own source frequencies? This isolates the
//! estimator/outgroup adequacy from any target-side data batch effect — no BAM is genotyped.
//!
//! For each site, target frequency = Σ wᵢ·sourceᵢ(site) with a known w, drawn as a diploid genome.
//! qpadm_fit must return ~w, feasible, not rejected. If it can't recover a self-consistent mixture,
//! the outgroups don't resolve WHG/ANF/Steppe (a construction problem, not a target problem).
//!   qpadm_selftest <qpadm_panel.bin>
use navigator_analysis::ancestry::{qpadm_fit, AncestryPanel, F4_BLOCK_BP};
use navigator_analysis::caller::SiteGenotype;

struct Lcg(u64);
impl Lcg {
    fn f(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn dosage(&mut self, p: f64) -> i32 {
        (self.f() < p) as i32 + (self.f() < p) as i32
    }
}

fn sg(contig: &str, pos: i64, dosage: i32) -> SiteGenotype {
    SiteGenotype {
        name: String::new(),
        contig: contig.to_string(),
        position: pos,
        reference_allele: "A".into(),
        alternate_allele: "G".into(),
        ploidy: 2,
        dosage,
        gq: 40,
        depth: 20,
        ref_depth: 10,
        alt_depth: 10,
        pls: vec![],
        gt: None,
        allele_depths: None,
    }
}

fn run(panel: &AncestryPanel, label: &str, truth: [f64; 3], seed: u64) {
    let idx = |c: &str| panel.populations.iter().position(|p| p == c).unwrap();
    let (wi, ai, si) = (idx("WHG"), idx("ANF"), idx("Steppe"));
    let sources = vec![wi, ai, si];
    let outgroups: Vec<usize> = (0..panel.populations.len()).filter(|i| !sources.contains(i)).collect();
    let mut rng = Lcg(seed);
    let gts: Vec<SiteGenotype> = panel
        .sites
        .iter()
        .map(|s| {
            let f = truth[0] * s.freqs[wi] as f64 + truth[1] * s.freqs[ai] as f64 + truth[2] * s.freqs[si] as f64;
            sg(&s.contig, s.position, rng.dosage(f.clamp(0.0, 1.0)))
        })
        .collect();
    match qpadm_fit(&gts, panel, &sources, &outgroups, F4_BLOCK_BP) {
        Some(f) => {
            println!(
                "{label:<28} WHG {:>6.1} ANF {:>6.1} Steppe {:>6.1}  | SE {:.0}/{:.0}/{:.0}  p {:.3}  {}",
                f.weights[0] * 100.0,
                f.weights[1] * 100.0,
                f.weights[2] * 100.0,
                f.std_errors[0] * 100.0,
                f.std_errors[1] * 100.0,
                f.std_errors[2] * 100.0,
                f.p_value,
                if f.weights_feasible(0.02) { "feasible" } else { "INFEASIBLE" },
            );
        }
        None => println!("{label:<28} qpadm_fit -> None"),
    }
}

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: qpadm_selftest <qpadm_panel.bin>");
    let panel = AncestryPanel::from_bytes(&std::fs::read(&path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    println!(
        "panel {} sites, outgroups {:?}\ntruth order [WHG, ANF, Steppe]:\n",
        panel.sites.len(),
        panel.populations[3..].iter().collect::<Vec<_>>()
    );
    run(&panel, "NW-Euro-ish 20/30/50", [0.20, 0.30, 0.50], 1);
    run(&panel, "pure Steppe 0/0/100", [0.0, 0.0, 1.0], 2);
    run(&panel, "pure WHG 100/0/0", [1.0, 0.0, 0.0], 3);
    run(&panel, "pure ANF 0/100/0", [0.0, 1.0, 0.0], 4);
    run(&panel, "Sardinian-ish 10/70/20", [0.10, 0.70, 0.20], 5);
    Ok(())
}

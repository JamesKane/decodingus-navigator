//! The §3.4 validation gates for the ancient (deep-ancestry) frequency panel.
//!
//! The previous ancient implementation shipped fabricated numbers, and it shipped them because
//! nobody ran the reference populations back through it. This does exactly that: it **simulates
//! individuals** from a known reference population's allele frequencies (drawing each genotype
//! `g ~ Binomial(2, f)` at the ancient panel's own sites), runs them through the shipping estimator
//! [`navigator_analysis::ancestry::estimate_ancient_admixture`], and prints what came back.
//!
//! Simulating rather than using real subjects is deliberate: it gives a *ground truth we control*
//! for every population — including the ones no workspace has a sample of (Sardinian, Yoruba) — and
//! it tests the panel and the model without the confound of the read pipeline. The gates:
//!
//! 1. **Sanity band.** A NW-European must land near Steppe 40–55 / ANF 25–40 / WHG 10–25. A
//!    Sardinian must be ANF-dominant with little Steppe. A Yoruba or Han must be rejected outright
//!    (dispersion above threshold → the estimator returns `None`), not handed a confident number.
//! 2. **Density.** `--downsample` re-runs on a random half of the sites; the answer must barely
//!    move, or the fit is ill-conditioned.
//! 3. **Fit residual.** The dispersion is reported for every population, so the separation between
//!    "fits" and "doesn't fit" is visible rather than asserted.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ancestry::{
    ancient_admixture_fit, estimate_admixture, estimate_ancient_admixture, west_eurasian_share, AncestryPanel,
};
use navigator_analysis::caller::SiteGenotype;

#[derive(Parser)]
pub struct ValidateAncientArgs {
    /// The ancient source panel to validate (`ancestry_freq_ancient_<build>.bin`).
    #[arg(long)]
    ancient: PathBuf,
    /// Reference frequency panel to simulate individuals FROM (`ancestry_freq_global_<build>.bin`).
    #[arg(long)]
    reference: PathBuf,
    /// The modern super-population panel (`ancestry_panel_<build>.bin`). Deep ancestry is scoped by
    /// the modern estimate, so the validator must score both models or it would be validating a
    /// policy the app doesn't actually run.
    #[arg(long)]
    panel: PathBuf,
    /// Reference populations to simulate, comma-separated.
    ///
    /// **Use 1000G populations only.** The `freq_global` asset's SGDP-derived columns (Sardinian,
    /// Basque, French, Orcadian, Han, …) record `0.0` where a population had no called sample —
    /// indistinguishable from a true zero — so 60%+ of their sites are fake zeros and an individual
    /// simulated from them is not that population, it's noise. The 1000G columns are called
    /// essentially everywhere and are the only trustworthy simulation source in that asset.
    /// (This is the same no-data-as-zero defect the ancient panel exists to avoid; see
    /// `pca::build_ancient_panel`.)
    #[arg(long, default_value = "GBR,CEU,FIN,TSI,IBS,PJL,CHB,JPT,YRI,LWK")]
    pops: String,
    /// Simulated individuals per population.
    #[arg(long, default_value_t = 20)]
    replicates: usize,
    /// Keep only this fraction of panel sites (the density gate). 1.0 = all sites.
    #[arg(long, default_value_t = 1.0)]
    downsample: f64,
    /// RNG seed (the simulation is fully deterministic given the seed).
    #[arg(long, default_value_t = 42)]
    seed: u64,
}

/// A deterministic, dependency-free PRNG (xorshift64*). Reproducibility matters more than quality
/// here: a validation run must give the same numbers twice.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Uniform in [0, 1).
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// A diploid alt-allele dosage drawn from HWE at frequency `f`.
    fn draw_dosage(&mut self, f: f64) -> i32 {
        (self.next_f64() < f) as i32 + (self.next_f64() < f) as i32
    }
}

/// Known mixtures to round-trip through the estimator, as `(label, weights over the panel axis)`.
/// The pure sources test conditioning (three nearly-collinear sources would smear a pure WHG across
/// all three); the mixture tests that the EM recovers proportions it was never told.
const RECOVERY_CASES: [(&str, [f64; 3]); 5] = [
    ("pure WHG", [1.0, 0.0, 0.0]),
    ("pure ANF", [0.0, 1.0, 0.0]),
    ("pure Steppe", [0.0, 0.0, 1.0]),
    ("20/30/50", [0.20, 0.30, 0.50]),
    ("50/25/25", [0.50, 0.25, 0.25]),
];

/// Gate 0 — **recovery**. Simulate individuals whose true ancestry we set ourselves, by drawing
/// genotypes from the panel's own source frequencies, and check the estimator gets them back.
///
/// This is the gate that would have caught the old implementation on its own: the sources here are
/// *exactly* the model's, so a correct estimator must return the input mixture, and the dispersion
/// must sit at ≈1.0 — that number is the model's noise floor, and it is what the thresholds for
/// real data are measured against.
fn recovery_gate(ancient: &AncestryPanel, super_panel: &AncestryPanel, args: &ValidateAncientArgs) -> Result<()> {
    anyhow::ensure!(
        ancient.populations.len() == 3,
        "recovery cases assume a three-source panel, got {:?}",
        ancient.populations
    );
    println!("── gate 0: recovery (simulated FROM the ancient panel — the true answer is known)\n");
    print!("{:<14}{:>10}", "true mixture", "");
    for p in &ancient.populations {
        print!("{:>14}", p);
    }
    println!("{:>18}{:>10}", "dispersion (worst)", "reported");

    for (label, truth) in RECOVERY_CASES {
        let mut sums = [0.0f64; 3];
        let mut disp = 0.0;
        let mut reported = 0usize;
        for rep in 0..args.replicates {
            let mut rng = Rng::new(args.seed ^ ((rep as u64) << 8) ^ (label.len() as u64));
            let genotypes: Vec<SiteGenotype> = ancient
                .sites
                .iter()
                .map(|site| {
                    // An admixed individual draws each allele from a source chosen by `truth`, so
                    // the alt-allele probability is exactly the mixture frequency Σ q_k·p_k.
                    let f: f64 = (0..3).map(|k| truth[k] * site.freqs[k] as f64).sum();
                    site_genotype(site, rng.draw_dosage(f))
                })
                .collect();
            let modern = estimate_admixture(&genotypes, super_panel, "validate");
            if let Some(r) = estimate_ancient_admixture(&genotypes, ancient, &modern, "validate") {
                reported += 1;
                disp += r.fit_distance.unwrap_or(f64::NAN);
                for (k, code) in ancient.populations.iter().enumerate() {
                    sums[k] += r
                        .components
                        .iter()
                        .find(|c| &c.population_code == code)
                        .map_or(0.0, |c| c.percentage);
                }
            }
        }
        let truth_s: Vec<String> = truth.iter().map(|t| format!("{:.0}", t * 100.0)).collect();
        print!("{:<14}{:>10}", label, format!("[{}]", truth_s.join("/")));
        for sum in sums {
            if reported > 0 {
                print!("{:>14.1}", sum / reported as f64);
            } else {
                print!("{:>14}", "—");
            }
        }
        println!(
            "{:>12}{:>10}",
            if reported > 0 {
                format!("{:.2}", disp / reported as f64)
            } else {
                "—".into()
            },
            format!("{reported}/{}", args.replicates)
        );
    }
    println!();
    Ok(())
}

pub fn validate_ancient(args: ValidateAncientArgs) -> Result<()> {
    let ancient = read_panel(&args.ancient)?;
    let reference = read_panel(&args.reference)?;
    let super_panel = read_panel(&args.panel)?;
    recovery_gate(&ancient, &super_panel, &args)?;
    println!("── gates 1–3: modern reference populations (simulated from the 1000G frequencies)");

    // The ancient panel's sites are the coordinate space we simulate in; look each one up in the
    // reference panel to get the source population's frequency there.
    let ref_idx: HashMap<(&str, i64), usize> = reference
        .sites
        .iter()
        .enumerate()
        .map(|(i, s)| ((s.contig.as_str(), s.position), i))
        .collect();

    let pops: Vec<&str> = args.pops.split(',').map(str::trim).filter(|p| !p.is_empty()).collect();
    println!(
        "ancient panel: {} sites × {:?}",
        ancient.sites.len(),
        ancient.populations
    );
    println!(
        "simulating {} individuals per population from {}{}\n",
        args.replicates,
        args.reference.display(),
        if args.downsample < 1.0 {
            format!(" (density gate: keeping {:.0}% of sites)", args.downsample * 100.0)
        } else {
            String::new()
        }
    );

    let width = ancient.populations.iter().map(|p| p.len().max(6)).collect::<Vec<_>>();
    print!("{:<12}{:>7}", "population", "sites");
    for (p, w) in ancient.populations.iter().zip(&width) {
        print!("{:>w$}", p, w = w + 8);
    }
    println!("{:>12}{:>10}", "dispersion", "reported");

    for pop in &pops {
        let Some(&pop_col) = reference.populations.iter().position(|p| p == pop).as_ref() else {
            println!("{pop:<12}  — not in the reference panel, skipped");
            continue;
        };

        let mut sums = vec![0.0f64; ancient.populations.len()];
        let mut sumsq = vec![0.0f64; ancient.populations.len()];
        let mut disp_sum = 0.0f64;
        let mut disp_sq = 0.0f64;
        let mut disp_max = f64::MIN;
        let mut eur_sum = 0.0f64;
        let mut reported = 0usize;
        let mut sites_used = 0usize;

        for rep in 0..args.replicates {
            // Seed per (population, replicate) so a run is reproducible and populations are
            // independent of each other's draw order.
            let mut rng = Rng::new(
                args.seed
                    .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                    .wrapping_add((pop.len() as u64) << 32)
                    .wrapping_add(pop.bytes().map(u64::from).sum::<u64>())
                    .wrapping_add(rep as u64),
            );
            let mut genotypes = Vec::with_capacity(ancient.sites.len());
            for site in &ancient.sites {
                if args.downsample < 1.0 && rng.next_f64() >= args.downsample {
                    continue;
                }
                let Some(&ri) = ref_idx.get(&(site.contig.as_str(), site.position)) else {
                    continue;
                };
                let f = reference.sites[ri].freqs[pop_col] as f64;
                genotypes.push(site_genotype(site, rng.draw_dosage(f)));
            }
            sites_used = genotypes.len();

            // Always fit, so the dispersion is reported whether or not the gate accepts it — the
            // threshold is only defensible if the separation it rests on is visible.
            let Some(fit) = ancient_admixture_fit(&genotypes, &ancient, "validate") else {
                continue;
            };
            let modern = estimate_admixture(&genotypes, &super_panel, "validate");
            eur_sum += west_eurasian_share(&modern);
            let d = fit.fit_distance.unwrap_or(f64::NAN);
            disp_sum += d;
            disp_sq += d * d;
            disp_max = disp_max.max(d);
            for (i, code) in ancient.populations.iter().enumerate() {
                let pct = fit
                    .components
                    .iter()
                    .find(|c| &c.population_code == code)
                    .map_or(0.0, |c| c.percentage);
                sums[i] += pct;
                sumsq[i] += pct * pct;
            }
            // …and separately record what the *shipping* estimator decided to report.
            if estimate_ancient_admixture(&genotypes, &ancient, &modern, "validate").is_some() {
                reported += 1;
            }
        }

        let n = args.replicates as f64;
        print!("{pop:<12}{sites_used:>7}");
        for (i, w) in width.iter().enumerate() {
            let mean = sums[i] / n;
            let var = (sumsq[i] / n - mean * mean).max(0.0);
            print!("{:>w$}", format!("{mean:.1}±{:.1}", var.sqrt()), w = w + 8);
        }
        // Report the worst individual as well as the mean: an applicability threshold set on means
        // would flicker — accepting most of a population and silently dropping its tail.
        let mean = disp_sum / n;
        let sd = (disp_sq / n - mean * mean).max(0.0).sqrt();
        println!(
            "{:>8.0}{:>18}{:>10}",
            eur_sum / n,
            format!("{mean:.2}±{sd:.2} (≤{disp_max:.2})"),
            if reported == 0 {
                "REJECTED".to_string()
            } else {
                format!("{reported}/{}", args.replicates)
            }
        );
    }

    println!(
        "\nPercentages are the RAW fit — for a REJECTED row they are exactly the numbers the\n\
         applicability gate exists to suppress. Expected: NW-European ≈ Steppe 40–55 / ANF 25–40 /\n\
         WHG 10–25 · non-European REJECTED (by EUR% scope or by dispersion) · dispersion ≈1 at the\n\
         model's noise floor."
    );
    Ok(())
}

fn site_genotype(site: &navigator_analysis::ancestry::PanelSite, dosage: i32) -> SiteGenotype {
    SiteGenotype {
        name: String::new(),
        contig: site.contig.clone(),
        position: site.position,
        reference_allele: site.reference_allele.to_string(),
        alternate_allele: site.alternate_allele.to_string(),
        ploidy: 2,
        dosage,
        gq: 99,
        depth: 30,
        ref_depth: 15,
        alt_depth: 15,
        pls: Vec::new(),
        gt: None,
        allele_depths: None,
    }
}

fn read_panel(path: &PathBuf) -> Result<AncestryPanel> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    AncestryPanel::from_bytes(&bytes).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
}

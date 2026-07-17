//! Throwaway: does the ancient panel's per-source zero-frequency rate depend on modern MAF?
//! Joins the ancient AF panel (WHG/ANF/Steppe) to the super-pop panel's EUR frequency by
//! (contig, pos), bins ancient sites by EUR MAF, and reports the per-source "freq==0" rate and
//! mean freq per bin. Hypothesis: at low modern MAF, the sparse sources (esp. WHG) collapse to 0
//! while Steppe stays nonzero, so the EM over-weights Steppe there.
use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let ancient_path = std::env::args().nth(1).expect("usage: ascertainment <ancient.bin> <super.bin>");
    let super_path = std::env::args().nth(2).expect("usage: ascertainment <ancient.bin> <super.bin>");
    let ancient = AncestryPanel::from_bytes(&std::fs::read(&ancient_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let sup = AncestryPanel::from_bytes(&std::fs::read(&super_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    let eur = sup.populations.iter().position(|p| p == "EUR").expect("super panel has no EUR");
    let eur_maf: HashMap<(String, i64), f64> = sup
        .sites
        .iter()
        .filter_map(|s| s.freqs.get(eur).map(|&f| ((s.contig.clone(), s.position), (f as f64).min(1.0 - f as f64))))
        .collect();

    // Bin edges on EUR MAF. Last bin captures the "common" (chip-like) sites.
    let edges = [0.0, 0.01, 0.05, 0.10, 0.20, 0.50001];
    let nb = edges.len() - 1;
    let bin_of = |m: f64| -> usize { edges.windows(2).position(|w| m >= w[0] && m < w[1]).unwrap_or(nb - 1) };

    let mut n = vec![0usize; nb];
    let mut zero = [vec![0usize; nb], vec![0usize; nb], vec![0usize; nb]]; // WHG, ANF, Steppe
    let mut sum = [vec![0.0f64; nb], vec![0.0f64; nb], vec![0.0f64; nb]];
    let mut joined = 0usize;

    for s in &ancient.sites {
        let Some(&m) = eur_maf.get(&(s.contig.clone(), s.position)) else { continue };
        joined += 1;
        let b = bin_of(m);
        n[b] += 1;
        for c in 0..3 {
            let f = s.freqs.get(c).copied().unwrap_or(0.0) as f64;
            if f == 0.0 {
                zero[c][b] += 1;
            }
            sum[c][b] += f;
        }
    }

    println!(
        "ancient sites={} joined to EUR MAF={} ({:.1}%)\n",
        ancient.sites.len(),
        joined,
        100.0 * joined as f64 / ancient.sites.len() as f64
    );
    println!(
        "{:<14}{:>7}   {:>6} {:>6} {:>6}   {:>7} {:>7} {:>7}",
        "EUR MAF bin", "sites", "WHG=0", "ANF=0", "Stp=0", "meanWHG", "meanANF", "meanStp"
    );
    for b in 0..nb {
        let nn = n[b].max(1) as f64;
        println!(
            "[{:.2},{:.2}){:>6}   {:>5.1}% {:>5.1}% {:>5.1}%   {:>7.3} {:>7.3} {:>7.3}",
            edges[b],
            edges[b + 1].min(0.5),
            n[b],
            100.0 * zero[0][b] as f64 / nn,
            100.0 * zero[1][b] as f64 / nn,
            100.0 * zero[2][b] as f64 / nn,
            sum[0][b] / nn,
            sum[1][b] / nn,
            sum[2][b] / nn,
        );
    }
    Ok(())
}

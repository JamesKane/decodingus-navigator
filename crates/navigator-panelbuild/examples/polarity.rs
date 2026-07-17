//! Throwaway: is the ancient AF panel's ref/alt orientation consistent with the (correctly
//! oriented) super-pop panel? estimate_admixture joins sample dosage to panel freqs by (contig,pos)
//! ONLY — no allele check — so a site where the ancient panel's alt is the super panel's ref feeds
//! an inverted dosage into the mixture. Tabulate aligned / swapped / mismatched, split by EUR MAF.
use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashMap;

fn main() -> anyhow::Result<()> {
    let a = std::env::args().nth(1).expect("usage: polarity <ancient.bin> <super.bin>");
    let s = std::env::args().nth(2).expect("super.bin");
    let ancient = AncestryPanel::from_bytes(&std::fs::read(&a)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let sup = AncestryPanel::from_bytes(&std::fs::read(&s)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let eur = sup.populations.iter().position(|p| p == "EUR").expect("no EUR");

    // super site index by (contig,pos) -> (ref, alt, eur_maf)
    let sup_idx: HashMap<(String, i64), (char, char, f64)> = sup
        .sites
        .iter()
        .filter_map(|x| {
            x.freqs
                .get(eur)
                .map(|&f| ((x.contig.clone(), x.position), (x.reference_allele, x.alternate_allele, (f as f64).min(1.0 - f as f64))))
        })
        .collect();

    let edges = [0.0, 0.05, 0.10, 0.20, 0.50001];
    let nb = edges.len() - 1;
    let bin = |m: f64| edges.windows(2).position(|w| m >= w[0] && m < w[1]).unwrap_or(nb - 1);
    let (mut aligned, mut swapped, mut other) = (vec![0usize; nb], vec![0usize; nb], vec![0usize; nb]);

    for x in &ancient.sites {
        let Some(&(sr, sa, m)) = sup_idx.get(&(x.contig.clone(), x.position)) else { continue };
        let b = bin(m);
        let (ar, aa) = (x.reference_allele, x.alternate_allele);
        if ar == sr && aa == sa {
            aligned[b] += 1;
        } else if ar == sa && aa == sr {
            swapped[b] += 1;
        } else {
            other[b] += 1;
        }
    }

    println!("{:<14}{:>8}{:>9}{:>9}{:>9}", "EUR MAF bin", "aligned", "SWAPPED", "other", "%swap");
    for b in 0..nb {
        let tot = (aligned[b] + swapped[b] + other[b]).max(1);
        println!(
            "[{:.2},{:.2}){:>8}{:>9}{:>9}{:>8.1}%",
            edges[b],
            edges[b + 1].min(0.5),
            aligned[b],
            swapped[b],
            other[b],
            100.0 * swapped[b] as f64 / tot as f64
        );
    }
    let ta: usize = aligned.iter().sum();
    let ts: usize = swapped.iter().sum();
    let to: usize = other.iter().sum();
    println!("\ntotal: aligned={ta} swapped={ts} other={to}  ({:.1}% swapped)", 100.0 * ts as f64 / (ta + ts + to).max(1) as f64);
    Ok(())
}

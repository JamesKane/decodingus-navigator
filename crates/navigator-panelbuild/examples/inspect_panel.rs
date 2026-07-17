//! Throwaway: dump an AncestryPanel's populations + per-population AF summary and pairwise Fst.
use navigator_analysis::ancestry::AncestryPanel;

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: inspect_panel <panel.bin>");
    let bytes = std::fs::read(&path)?;
    let panel = AncestryPanel::from_bytes(&bytes).map_err(|e| anyhow::anyhow!("{e}"))?;
    println!("build={} sites={} pops={}", panel.build, panel.sites.len(), panel.populations.len());
    println!("populations: {:?}", panel.populations);

    let k = panel.populations.len();
    // Per-pop: mean AF, and how many sites are exactly 0.0 (the builder's "no data" sentinel).
    let mut sum = vec![0.0f64; k];
    let mut zeros = vec![0usize; k];
    let mut fixed = vec![0usize; k]; // 0.0 or 1.0
    for s in &panel.sites {
        for i in 0..k.min(s.freqs.len()) {
            let f = s.freqs[i] as f64;
            sum[i] += f;
            if f == 0.0 {
                zeros[i] += 1;
            }
            if f == 0.0 || f == 1.0 {
                fixed[i] += 1;
            }
        }
    }
    let n = panel.sites.len() as f64;
    println!("\npop            meanAF   %AF==0   %fixed");
    for i in 0..k {
        println!(
            "{:<14} {:>6.3}  {:>6.1}%  {:>6.1}%",
            panel.populations[i],
            sum[i] / n,
            100.0 * zeros[i] as f64 / n,
            100.0 * fixed[i] as f64 / n
        );
    }

    // Pairwise Nei Fst between a few populations of interest.
    let want = ["WHG", "ANF", "Steppe", "EHG", "CHG", "Iran_N", "GBR", "CEU", "TSI", "YRI", "Han"];
    let idx: Vec<(usize, &str)> = want
        .iter()
        .filter_map(|w| panel.populations.iter().position(|p| p == w).map(|i| (i, *w)))
        .collect();
    println!("\npairwise Nei Fst (sites where BOTH pops have data, i.e. not both-fixed-at-0):");
    print!("{:<10}", "");
    for (_, b) in &idx {
        print!("{:>8}", b);
    }
    println!();
    for (ia, a) in &idx {
        print!("{:<10}", a);
        for (ib, _) in &idx {
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            let mut used = 0usize;
            for s in &panel.sites {
                if s.freqs.len() != k {
                    continue;
                }
                let (p1, p2) = (s.freqs[*ia] as f64, s.freqs[*ib] as f64);
                let hs = (p1 * (1.0 - p1) + p2 * (1.0 - p2)) / 2.0;
                let pbar = (p1 + p2) / 2.0;
                let ht = pbar * (1.0 - pbar);
                if ht > 0.0 {
                    num += ht - hs;
                    den += ht;
                    used += 1;
                }
            }
            let _ = used;
            print!("{:>8.3}", if den > 0.0 { num / den } else { f64::NAN });
        }
        println!();
    }
    Ok(())
}

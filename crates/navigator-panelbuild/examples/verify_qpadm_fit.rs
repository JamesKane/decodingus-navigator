//! Verify our own `qpadm_fit` (the app-path pooled-frequency estimator) reproduces the admixtools2
//! Patterson-2022 result on James's real WGS. Builds an in-memory CHM13 frequency panel from the AADR
//! `.traw` (per-individual, hg19) — reconciling each population's frequency to the CHM13 bed_alt
//! orientation — plus James's genotypes, and runs `qpadm_fit`. Expect ~WHG 15 / EEF 45 / Steppe 41.
//!   verify_qpadm_fit <pat.traw> <pat.ind> <rs_chm13.tsv> <rs_bedalt.tsv> <james.tsv>
//! rs_chm13.tsv: `rsid<TAB>contig<TAB>pos` (CHM13); james.tsv: `rsid<TAB>contig<TAB>pos<TAB>dosage`.
use navigator_analysis::ancestry::{qpadm_fit, AncestryPanel, PanelSite, Pop, Quartet, F4_BLOCK_BP};
use navigator_analysis::caller::SiteGenotype;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

fn comp(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        x => x,
    }
}

fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (traw, ind, rschm, rsba, james) = (&a[1], &a[2], &a[3], &a[4], &a[5]);

    // Population order: sources first, then outgroups.
    let pops = ["WHG", "EEF", "Steppe", "AnatoliaOG", "Afanasievo", "IronGates", "African"];
    let pop_idx: HashMap<&str, usize> = pops.iter().enumerate().map(|(i, &p)| (p, i)).collect();

    // .ind → label per traw sample column (skip the appended Target row).
    let labels: Vec<Option<usize>> = std::fs::read_to_string(ind)?
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').nth(2).and_then(|lab| pop_idx.get(lab).copied()))
        .collect();

    // rsID → CHM13 (contig,pos); rsID → bed_alt; rsID → James dosage(count bed_alt).
    let mut chm: HashMap<String, (String, i64)> = HashMap::new();
    for l in std::fs::read_to_string(rschm)?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() >= 3 {
            chm.insert(f[0].into(), (f[1].into(), f[2].parse().unwrap_or(0)));
        }
    }
    let mut bedalt: HashMap<String, u8> = HashMap::new();
    for l in std::fs::read_to_string(rsba)?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() >= 2 && !f[1].is_empty() {
            bedalt.insert(f[0].into(), f[1].as_bytes()[0]);
        }
    }
    let mut jd: HashMap<String, i32> = HashMap::new();
    for l in std::fs::read_to_string(james)?.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() >= 4 {
            jd.insert(f[0].into(), f[3].parse().unwrap_or(-1));
        }
    }

    let k = pops.len();
    let mut sites = Vec::new();
    let mut genos = Vec::new();
    let file = std::fs::File::open(traw)?;
    for (li, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if li == 0 {
            continue; // header
        }
        let f: Vec<&str> = line.split('\t').collect();
        let rsid = f[1];
        let Some(&(ref contig, pos)) = chm.get(rsid) else { continue };
        let Some(&ba) = bedalt.get(rsid) else { continue };
        let counted = f[4].as_bytes()[0];
        let alt = f[5].as_bytes()[0];
        // Per-population sum of COUNTED-allele copies and called count.
        let mut sum = vec![0.0f64; k];
        let mut n = vec![0u32; k];
        for (c, cell) in f[6..].iter().enumerate() {
            if let Some(p) = labels.get(c).copied().flatten() {
                if *cell != "NA" {
                    sum[p] += cell.parse::<f64>().unwrap_or(0.0);
                    n[p] += 1;
                }
            }
        }
        if n.contains(&0) {
            continue; // require all pops present at the site
        }
        // Orient frequency to count bed_alt (matches James's raw dosage), reconciling COUNTED/ALT.
        let flip = if ba == counted || comp(ba) == counted {
            false
        } else if ba == alt || comp(ba) == alt {
            true
        } else {
            continue;
        };
        let freqs: Vec<f32> = (0..k)
            .map(|p| {
                let f_counted = sum[p] / (2.0 * n[p] as f64);
                (if flip { 1.0 - f_counted } else { f_counted }) as f32
            })
            .collect();
        let Some(&d) = jd.get(rsid) else { continue };
        if d < 0 {
            continue;
        }
        sites.push(PanelSite {
            contig: contig.clone(),
            position: pos,
            reference_allele: 'A',
            alternate_allele: 'G',
            freqs,
        });
        genos.push(SiteGenotype {
            name: String::new(),
            contig: contig.clone(),
            position: pos,
            reference_allele: "A".into(),
            alternate_allele: "G".into(),
            ploidy: 2,
            dosage: d,
            gq: 40,
            depth: 20,
            ref_depth: 10,
            alt_depth: 10,
            pls: vec![],
            gt: None,
            allele_depths: None,
        });
    }
    eprintln!("panel {} sites, {} James genotypes", sites.len(), genos.len());
    let panel = AncestryPanel {
        build: "chm13v2.0".into(),
        populations: pops.iter().map(|s| s.to_string()).collect(),
        sites,
    };
    // Cross-check f4 wiring compiles for this panel (unused sanity ref).
    let _ = Quartet::new(Pop::Target, Pop::Ref(0), Pop::Ref(3), Pop::Ref(4));

    let sources = [0usize, 1, 2];
    let outgroups = [3usize, 4, 5, 6];
    let fit = qpadm_fit(&genos, &panel, &sources, &outgroups, F4_BLOCK_BP)
        .ok_or_else(|| anyhow::anyhow!("qpadm_fit returned None"))?;
    println!("\n== our qpadm_fit — James (Patterson config) ==");
    println!("sites {}  blocks {}  dof {}  chi2 {:.2}  p {:.4}", fit.n_sites, fit.n_blocks, fit.dof, fit.chi2, fit.p_value);
    for (c, i) in ["WHG", "EEF", "Steppe"].iter().zip(0..) {
        println!("  {c:<8} {:>6.1} %  (SE {:.1})", fit.weights[i] * 100.0, fit.std_errors[i] * 100.0);
    }
    println!(
        "model {} at p=0.05; weights {}",
        if fit.p_value >= 0.05 { "ACCEPTED" } else { "REJECTED" },
        if fit.weights_feasible(0.02) { "feasible" } else { "INFEASIBLE" }
    );
    Ok(())
}

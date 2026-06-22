//! One-off placement debug for a vendor Y VCF (FTDNA Big Y / YSEQ / Full Genomes): parse the VCF's
//! chrY calls, place them on the FTDNA tree, and print the terminal + lineage. Mirrors the app's
//! `place_chip_panel` (score → path_admissible → deepen_terminal) without the store/network.
//!
//!   cargo run --release --example place_y_vcf -p navigator-analysis -- <variants.vcf> <tree.json>

use std::collections::HashMap;
use std::path::Path;

use navigator_analysis::haplo;

/// Genotype-aware chrY SNP calls from a VCF: honor the sample GT (drop 0/0 / ./.), single-base only.
fn chr_y_calls(path: &Path) -> HashMap<i64, char> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = HashMap::new();
    let mut rows = 0usize;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 || !(f[0].eq_ignore_ascii_case("chrY") || f[0].eq_ignore_ascii_case("y")) {
            continue;
        }
        rows += 1;
        let Ok(pos) = f[1].parse::<i64>() else { continue };
        if f[4] == "." {
            continue;
        }
        let alts: Vec<&str> = f[4].split(',').collect();
        let gt = (f.len() >= 10)
            .then(|| {
                f[8].split(':')
                    .position(|k| k == "GT")
                    .and_then(|i| f[9].split(':').nth(i))
            })
            .flatten();
        let allele = match gt {
            Some(gt) => match gt
                .split(['/', '|'])
                .filter_map(|a| a.parse::<usize>().ok())
                .find(|&a| a > 0)
            {
                Some(idx) => match alts.get(idx - 1) {
                    Some(&a) => a,
                    None => continue,
                },
                None => continue,
            },
            None => alts[0],
        };
        if allele.len() == 1 {
            if let Some(b) = allele.chars().next() {
                out.insert(pos, b.to_ascii_uppercase());
            }
        }
    }
    eprintln!("chrY rows={rows} derived calls={}", out.len());
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 3 {
        eprintln!("usage: place_y_vcf <variants.vcf> <tree.json>");
        std::process::exit(2);
    }
    let (vcf, tree_path) = (Path::new(&a[1]), &a[2]);

    let tree = haplo::parse_ftdna_json(&std::fs::read_to_string(tree_path).unwrap()).unwrap();
    let calls = chr_y_calls(vcf);

    let ranked = haplo::score(&tree, &calls);
    println!("\n=== top 10 by Kulczynski ===");
    for r in ranked.iter().take(10) {
        println!(
            "  {:<16} score={:.4} depth={:>2} matched={}/{} admissible={}",
            r.name,
            r.score,
            r.depth,
            r.matched,
            r.expected,
            haplo::path_admissible(&tree, &calls, r.id),
        );
    }

    let start = ranked
        .iter()
        .find(|r| haplo::path_admissible(&tree, &calls, r.id))
        .map(|r| (r.id, r.name.clone()));
    if let Some((sid, sname)) = start {
        let term = haplo::deepen_terminal(&tree, &calls, sid);
        let tname = tree.nodes.get(&term).map(|n| n.name.as_str()).unwrap_or("?");
        println!("\nfirst admissible start = {sname}; deepen_terminal -> {tname}");
    } else {
        println!("\nno admissible placement");
    }
}

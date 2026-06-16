//! STR convention calibration: run the caller on a corpus of Big Y kits (BAM + FTDNA DYS CSV each)
//! and tabulate, per marker, the offset (ftdna − caller) distribution across kits — to classify each
//! marker reliable / convention-offset / variable-exclude and measure callability.
//! Usage: cargo run --release -p navigator-analysis --example str_calibrate -- <kits_dir> <hipstr.bed.gz>
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use navigator_analysis::strcaller::{genotype_str_loci, StrCallerParams, StrConfidence};
use navigator_analysis::strref::{load_hipstr_contig, StrLocus};

/// Find (csv, bam) under a kit folder: a `*_DYS_Results*.csv` (any level) + a `*.bam` (any level).
fn kit_pair(dir: &Path) -> Option<(PathBuf, PathBuf)> {
    let mut csv = None;
    let mut bam = None;
    for e in walkdir(dir) {
        let n = e.file_name().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        if n.ends_with(".csv") && n.contains("dys") {
            csv = Some(e.clone());
        } else if n.ends_with(".bam") {
            bam = Some(e.clone());
        }
    }
    Some((csv?, bam?))
}

fn walkdir(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(rd) = std::fs::read_dir(dir) {
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                out.extend(walkdir(&p));
            } else {
                out.push(p);
            }
        }
    }
    out
}

/// Parse the wide FTDNA DYS CSV (row1 names, row2 values; quoted, leading spaces, "-" = no call).
fn parse_ftdna(csv: &Path) -> HashMap<String, String> {
    let text = std::fs::read_to_string(csv).unwrap_or_default();
    let mut lines = text.lines();
    let names: Vec<String> = lines.next().unwrap_or("").split(',').map(|s| s.trim().trim_matches('"').trim().to_string()).collect();
    let vals: Vec<String> = lines.next().unwrap_or("").split(',').map(|s| s.trim().trim_matches('"').trim().to_string()).collect();
    names.into_iter().zip(vals).collect()
}

fn base_marker(n: &str) -> String {
    let n = n.split('/').next().unwrap_or(n);
    let n = n.split('_').next().unwrap_or(n);
    n.split('.').next().unwrap_or(n).to_string()
}

fn main() {
    let mut a = std::env::args().skip(1);
    let kits_dir = PathBuf::from(a.next().expect("kits dir"));
    let bed = PathBuf::from(a.next().expect("hipstr bed"));

    let loci: Vec<StrLocus> = load_hipstr_contig(&bed, "chrY", 2).expect("load bed");
    eprintln!("{} chrY loci", loci.len());

    // marker -> Vec<(offset, kit)> for single-copy comparisons; and callability counts.
    let mut offsets: HashMap<String, Vec<i32>> = HashMap::new();
    let mut callable: HashMap<String, usize> = HashMap::new();
    let mut n_kits = 0;

    let mut folders: Vec<PathBuf> = std::fs::read_dir(&kits_dir).unwrap().flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
    folders.sort();
    for folder in folders {
        let Some((csv, bam)) = kit_pair(&folder) else { continue };
        let ftdna = parse_ftdna(&csv);
        if ftdna.is_empty() {
            continue;
        }
        let genos = match genotype_str_loci(&bam, "chrY", &loci, 1, &StrCallerParams::default(), None) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("  {}: caller error {e}", folder.display());
                continue;
            }
        };
        n_kits += 1;
        eprintln!("  [{n_kits}] {} -> {} genotypes", folder.file_name().unwrap().to_string_lossy(), genos.len());
        // One called value per base marker (high-confidence, single-allele = single-copy).
        for g in genos.iter().filter(|g| g.confidence != StrConfidence::Low && g.alleles.len() == 1) {
            let m = base_marker(&g.name);
            *callable.entry(m.clone()).or_default() += 1;
            if let Some(fv) = ftdna.get(&m) {
                if let Ok(f) = fv.parse::<i32>() {
                    offsets.entry(m).or_default().push(f - g.alleles[0]);
                }
            }
        }
    }

    // Per-marker classification.
    println!("\n# marker  n  callable_kits  modal_offset  agreement%  class");
    let mut markers: Vec<&String> = offsets.keys().collect();
    markers.sort();
    let (mut reliable, mut convention, mut variable) = (0, 0, 0);
    for m in markers {
        let o = &offsets[m];
        let mut counts: BTreeMap<i32, usize> = BTreeMap::new();
        for &d in o {
            *counts.entry(d).or_default() += 1;
        }
        let (modal, mc) = counts.iter().max_by_key(|(_, c)| **c).map(|(d, c)| (*d, *c)).unwrap();
        let agree = mc as f64 / o.len() as f64;
        let class = if agree < 0.7 {
            variable += 1;
            "VARIABLE/exclude"
        } else if modal == 0 {
            reliable += 1;
            "reliable"
        } else if modal.abs() <= 3 {
            convention += 1;
            "convention-offset"
        } else {
            variable += 1;
            "large/exclude"
        };
        println!("{m:<14} {:>2} {:>13} {:>+12} {:>9.0}%  {class}", o.len(), callable.get(m).copied().unwrap_or(0), modal, agree * 100.0);
    }
    println!("\n# {n_kits} kits · {reliable} reliable(offset 0) · {convention} convention-offset(±1-3) · {variable} variable/large(exclude)");
}

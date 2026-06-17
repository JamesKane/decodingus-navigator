//! STR convention calibration: run the caller on a corpus of Big Y kits (BAM + FTDNA DYS CSV each)
//! and tabulate, per marker, the offset (ftdna − caller) distribution across kits — to classify each
//! marker reliable / convention-offset / variable-exclude and measure callability.
//! Usage: cargo run --release -p navigator-analysis --example str_calibrate -- <kits_dir> <hipstr.bed.gz>
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use navigator_analysis::strcaller::{genotype_str_loci, StrCallerParams, StrConfidence};
use navigator_analysis::strmarker::{to_ftdna, MarkerStatus};
use navigator_analysis::strref::{load_hipstr_contig, StrLocus};

/// Find (csv, alignment) under a kit folder: a `*_DYS_Results*.csv` + an alignment, preferring the
/// CHM13-realigned `.cram` (chrYM) over the original `.bam`.
fn kit_pair(dir: &Path) -> Option<(PathBuf, PathBuf)> {
    let mut csv = None;
    let mut cram = None;
    let mut bam = None;
    for e in walkdir(dir) {
        let n = e.file_name().and_then(|s| s.to_str()).unwrap_or("").to_ascii_lowercase();
        if n.ends_with(".csv") && n.contains("dys") {
            csv = Some(e.clone());
        } else if n.ends_with(".cram") {
            cram = Some(e.clone());
        } else if n.ends_with(".bam") {
            bam = Some(e.clone());
        }
    }
    Some((csv?, cram.or(bam)?))
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
    let reference = a.next().map(PathBuf::from); // required for CRAM decode (CHM13 FASTA)

    let loci: Vec<StrLocus> = load_hipstr_contig(&bed, "chrY", 2).expect("load bed");
    eprintln!("{} chrY loci · reference {:?}", loci.len(), reference);

    // marker -> Vec<offset> for single-copy comparisons; and callability counts.
    let mut offsets: HashMap<String, Vec<i32>> = HashMap::new();
    let mut callable: HashMap<String, usize> = HashMap::new();
    let (mut n_kits, mut n_skipped_swap, mut n_panic) = (0, 0, 0);
    // One bad CRAM (an unsupported noodles compression codec, etc.) must not abort the whole corpus.
    std::panic::set_hook(Box::new(|_| {}));

    let mut folders: Vec<PathBuf> = std::fs::read_dir(&kits_dir).unwrap().flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
    folders.sort();
    for folder in folders {
        let Some((csv, aln)) = kit_pair(&folder) else { continue };
        let ftdna = parse_ftdna(&csv);
        if ftdna.is_empty() {
            continue;
        }
        let kit = folder.file_name().unwrap().to_string_lossy().into_owned();
        let called = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            genotype_str_loci(&aln, "chrY", &loci, 1, &StrCallerParams::default(), reference.as_deref())
        }));
        let genos = match called {
            Ok(Ok(g)) => g,
            Ok(Err(e)) => {
                eprintln!("  {kit}: caller error {e}");
                continue;
            }
            Err(_) => {
                n_panic += 1;
                eprintln!("  {kit}: PANIC (unsupported CRAM codec?) — skipped");
                continue;
            }
        };
        let single: Vec<_> = genos.iter().filter(|g| g.confidence != StrConfidence::Low && g.alleles.len() == 1).collect();

        // QC: end-user-provided data may have a swapped BAM↔CSV. Fingerprint the kit on the markers
        // strmarker already calls Reliable (offset 0): if ≥10 are comparable but <70% match this
        // kit's own CSV, the alignment and the CSV are likely different people — exclude it.
        let (mut rel_ok, mut rel_tot) = (0, 0);
        for g in &single {
            let cm = to_ftdna(&g.name, g.alleles[0]);
            if cm.status == MarkerStatus::Reliable {
                if let Some(f) = ftdna.get(&cm.marker).and_then(|v| v.parse::<i32>().ok()) {
                    rel_tot += 1;
                    if f == cm.value {
                        rel_ok += 1;
                    }
                }
            }
        }
        if rel_tot >= 10 && (rel_ok as f64 / rel_tot as f64) < 0.7 {
            n_skipped_swap += 1;
            eprintln!("  SKIP {kit} — likely BAM/CSV swap ({rel_ok}/{rel_tot} reliable agree)");
            continue;
        }

        n_kits += 1;
        eprintln!("  [{n_kits}] {kit} -> {} genotypes ({rel_ok}/{rel_tot} reliable ok)", genos.len());
        // Re-derive offsets from RAW caller value vs CSV (independent of strmarker's current table).
        for g in &single {
            let m = base_marker(&g.name);
            *callable.entry(m.clone()).or_default() += 1;
            if let Some(f) = ftdna.get(&m).and_then(|v| v.parse::<i32>().ok()) {
                offsets.entry(m).or_default().push(f - g.alleles[0]);
            }
        }
    }
    eprintln!("calibrated on {n_kits} kits ({n_skipped_swap} skipped as likely swaps, {n_panic} CRAM-decode panics)");

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

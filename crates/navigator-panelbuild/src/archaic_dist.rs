//! Build the archaic percentile reference (`archaic_marker_dist_<build>.bin`), which is Asset 4 of
//! `documents/design/ArchaicAncestry_Design.md` §4.
//!
//! Tier A reports a **count**, and a count on its own means nothing without a cohort to set it
//! against. Design §1 gives the reason: 23andMe shows a percentile beside its count for exactly
//! this purpose.
//!
//! So this builder scores every reference sample through the same `count_archaic_markers`
//! arithmetic that the app runs on a subject. The percentile then compares like with like.
//!
//! **The asset holds a count for each fine population, and does not reduce them to
//! super-populations first** (design §9 Q3).
//!
//! Version 1 draws the percentile against a super-population. To key it to the user's inferred
//! fine ancestry would let an error in that ancestry move the archaic figure, with no warning. But
//! the fine-grained counts stay, so a fine-pop percentile later needs only a new key, and not a new
//! asset.
//!
//! The input is a `bcftools query` matrix at the panel sites, plus the 1kGP `sample/pop/super_pop`
//! panel. That is how every other stage feeds this crate.

use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::archaic::{ArchaicCountDistribution, ArchaicMarkerPanel, CohortCounts};

use crate::pca::{open_maybe_gz, write_bin};

#[derive(Parser)]
pub struct ArchaicDistArgs {
    /// The marker panel that a build produced (`archaic_markers_<build>.bin`). The distribution
    /// means something only against the panel that scored it, so the output records that panel's
    /// site count.
    #[arg(long)]
    panel: PathBuf,
    /// Genotype matrix at the panel sites: `CHROM POS REF ALT [GT...]` from `bcftools query`.
    #[arg(long)]
    matrix: PathBuf,
    /// The sample ids, one on each line, in the column order of the matrix (`bcftools query -l`).
    #[arg(long)]
    samples: PathBuf,
    /// A panel in the 1kGP style: `sample<TAB>pop<TAB>super_pop[...]`, with a header line. The
    /// builder drops a sample that this panel does not hold. That is also how a related sample
    /// stays out, because the unrelated-2504 panel is the one that carries the labels.
    #[arg(long)]
    pops: PathBuf,
    /// Extra site sets to measure variance inflation on, as `contig<TAB>pos` files (repeatable).
    ///
    /// This matters for a real user. A random subset of the panel, and a real capture or chip
    /// subset of the SAME size, behave very differently. The content of an array sits in clusters
    /// along the genome, so it holds far more linkage.
    ///
    /// Measured here: a random 3 % subset gives EUR an inflation of 2.0x, and the real 1240k
    /// intersection, which is 2.6 % of the panel, gives 5.3x.
    ///
    /// So give this option the site sets that users have, such as 1240k, or the manifest of a
    /// consumer array. The interpolation then stays honest. Without them it inflates too little,
    /// and every percentile becomes over-confident.
    #[arg(long = "subset-sites")]
    subset_sites: Vec<PathBuf>,
    /// Output (bincode `ArchaicCountDistribution`).
    #[arg(long)]
    out: PathBuf,
}

/// `sample -> (pop, super_pop)`.
fn load_pops(path: &Path) -> Result<HashMap<String, (String, String)>> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = HashMap::new();
    for (i, line) in rdr.lines().enumerate() {
        let line = line?;
        let f: Vec<&str> = line.split(['\t', ' ']).filter(|s| !s.is_empty()).collect();
        if f.len() < 3 {
            continue;
        }
        // Skip a header row, and do not depend on its exact text.
        if i == 0 && f[0].eq_ignore_ascii_case("sample") {
            continue;
        }
        out.insert(f[0].to_string(), (f[1].to_string(), f[2].to_string()));
    }
    Ok(out)
}

fn load_samples(path: &Path) -> Result<Vec<String>> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(rdr
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect())
}

pub fn build_archaic_dist(args: ArchaicDistArgs) -> Result<()> {
    let bytes = std::fs::read(&args.panel).with_context(|| format!("reading {}", args.panel.display()))?;
    let panel = ArchaicMarkerPanel::from_bytes(&bytes).map_err(|e| anyhow::anyhow!("{e}"))?;
    let samples = load_samples(&args.samples)?;
    let pops = load_pops(&args.pops)?;
    eprintln!(
        "panel {} sites, {} matrix samples, {} labelled samples",
        panel.len(),
        samples.len(),
        pops.len()
    );

    // The derived allele at each site, keyed by position. The matrix and the panel agree on the
    // CHM13 coordinates, but the matrix carries its own REF and ALT. So the match uses the derived
    // BASE. A dosage would invert at every site where the panel's orientation puts the derived
    // allele on REF.
    let derived: HashMap<(&str, i64), char> = panel
        .sites
        .iter()
        .map(|s| {
            (
                (s.contig.as_str(), s.position),
                s.archaic_derived_allele.to_ascii_uppercase(),
            )
        })
        .collect();

    // The panel-order index of each site, so a frequency lands in the order that the runtime uses
    // as its index.
    let site_index: HashMap<(&str, i64), usize> = panel
        .sites
        .iter()
        .enumerate()
        .map(|(i, s)| ((s.contig.as_str(), s.position), i))
        .collect();

    let mut totals = vec![0u32; samples.len()];
    let mut scored = vec![0u32; samples.len()];
    // The derived-allele count of each (site, super-population) pair, for the frequency table.
    let mut super_of: Vec<Option<String>> = Vec::with_capacity(samples.len());
    for smp in &samples {
        super_of.push(pops.get(smp).map(|(_, sup)| sup.clone()));
    }
    let mut supers: Vec<String> = super_of.iter().flatten().cloned().collect();
    supers.sort();
    supers.dedup();
    let super_idx: HashMap<&str, usize> = supers.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut ac = vec![vec![0u32; panel.len()]; supers.len()];
    let mut an = vec![vec![0u32; panel.len()]; supers.len()];

    // The code measures the variance inflation at a LADDER of densities, and not once.
    //
    // That inflation depends strongly on how many sites a subset takes from each linked haplotype
    // block. It is 52x on the full panel, and 5.3x on a subset that is 2.6% of the same panel.
    //
    // So one factor, applied to a chip, would widen the spread about 3 times too much, and push
    // every percentile toward 50. The rungs are deterministic subsets by site index, one inside the
    // next, and this pass accumulates them with everything else.
    const DENSITY_LADDER: [f32; 5] = [1.0, 0.3, 0.1, 0.03, 0.01];
    let in_rung = |site_idx: usize, rung: usize| -> bool {
        // The rungs are deterministic, and one sits inside the next. Rung r keeps every site whose
        // index falls in the first `density` part of a fixed stride. So a sparser rung is a subset
        // of a denser one.
        let d = DENSITY_LADDER[rung];
        ((site_idx as u64).wrapping_mul(2654435761) % 1_000_000) < (d as f64 * 1_000_000.0) as u64
    };
    let mut rung_tot = vec![vec![0u32; samples.len()]; DENSITY_LADDER.len()];

    // Real site sets (1240k, chip manifests) measured alongside the synthetic rungs.
    let mut subset_masks: Vec<(String, Vec<bool>)> = Vec::new();
    for path in &args.subset_sites {
        let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
        let mut keys: std::collections::HashSet<(String, i64)> = std::collections::HashSet::new();
        for line in rdr.lines() {
            let line = line?;
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() >= 2 {
                if let Ok(pos) = f[1].parse::<i64>() {
                    keys.insert((f[0].to_string(), pos));
                }
            }
        }
        let mask: Vec<bool> = panel
            .sites
            .iter()
            .map(|s| keys.contains(&(s.contig.clone(), s.position)))
            .collect();
        let n = mask.iter().filter(|b| **b).count();
        let label = path
            .file_stem()
            .and_then(|x| x.to_str())
            .unwrap_or("subset")
            .to_string();
        eprintln!("subset {label}: {n} of {} panel sites", panel.len());
        if n > 0 {
            subset_masks.push((label, mask));
        }
    }
    let mut subset_tot = vec![vec![0u32; samples.len()]; subset_masks.len()];
    let (mut rows, mut matched, mut skipped) = (0usize, 0usize, 0usize);

    let rdr = open_maybe_gz(&args.matrix).with_context(|| format!("opening {}", args.matrix.display()))?;
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        rows += 1;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let Ok(pos) = f[1].parse::<i64>() else { continue };
        let Some(&d) = derived.get(&(f[0], pos)) else {
            skipped += 1;
            continue;
        };
        let (Some(r), Some(a)) = (f[2].chars().next(), f[3].chars().next()) else {
            continue;
        };
        // Which allele index carries the derived base? A site whose alleles match neither one is a
        // disagreement between the matrix and the panel. The code skips it, and makes no guess.
        let (r, a) = (r.to_ascii_uppercase(), a.to_ascii_uppercase());
        let derived_idx = if a == d {
            1u8
        } else if r == d {
            0u8
        } else {
            skipped += 1;
            continue;
        };
        matched += 1;

        for (col, gt) in f[4..].iter().enumerate() {
            if col >= samples.len() {
                break;
            }
            let mut called = false;
            let mut copies = 0u32;
            for allele in gt.split(['/', '|']) {
                match allele.trim() {
                    "0" => {
                        called = true;
                        copies += u32::from(derived_idx == 0);
                    }
                    "1" => {
                        called = true;
                        copies += u32::from(derived_idx == 1);
                    }
                    _ => {}
                }
            }
            if called {
                scored[col] += 1;
                totals[col] += copies;
                if let (Some(sup), Some(&si)) = (super_of[col].as_deref(), site_index.get(&(f[0], pos))) {
                    if let Some(pi) = super_idx.get(sup) {
                        ac[*pi][si] += copies;
                        an[*pi][si] += 2;
                    }
                    for (r, tot) in rung_tot.iter_mut().enumerate() {
                        if in_rung(si, r) {
                            tot[col] += copies;
                        }
                    }
                    for (m, tot) in subset_tot.iter_mut().enumerate() {
                        if subset_masks[m].1[si] {
                            tot[col] += copies;
                        }
                    }
                }
            }
        }
    }
    eprintln!("matrix rows {rows}, matched to panel {matched}, skipped {skipped}");
    anyhow::ensure!(
        matched > 0,
        "no matrix row matched the panel — wrong build or wrong matrix?"
    );

    // Group by fine population, and keep the counts sorted, so a percentile is one quick scan.
    let mut by_pop: BTreeMap<(String, String), Vec<u32>> = BTreeMap::new();
    let mut unlabelled = 0usize;
    for (i, s) in samples.iter().enumerate() {
        let Some((pop, sup)) = pops.get(s) else {
            unlabelled += 1;
            continue;
        };
        // A sample with no called sites would drag the distribution toward zero for a reason that
        // has nothing to do with its archaic ancestry.
        if scored[i] == 0 {
            continue;
        }
        by_pop.entry((pop.clone(), sup.clone())).or_default().push(totals[i]);
    }
    eprintln!("{unlabelled} matrix samples had no population label (dropped — related/unlabelled set)");

    let cohorts: Vec<CohortCounts> = by_pop
        .into_iter()
        .map(|((population, super_population), mut counts)| {
            counts.sort_unstable();
            CohortCounts {
                population,
                super_population,
                counts,
            }
        })
        .collect();
    anyhow::ensure!(!cohorts.is_empty(), "no labelled sample survived — check --pops");

    let mut by_super: BTreeMap<&str, (usize, u64)> = BTreeMap::new();
    for c in &cohorts {
        let e = by_super.entry(c.super_population.as_str()).or_insert((0, 0));
        e.0 += c.counts.len();
        e.1 += c.counts.iter().map(|&v| v as u64).sum::<u64>();
    }
    eprintln!("cohorts: {} populations", cohorts.len());
    for (sup, (n, sum)) in &by_super {
        eprintln!("  {sup:<5} n={n:<5} mean archaic copies {:.0}", *sum as f64 / *n as f64);
    }

    // The derived-allele frequency at each site, for each super-population. A site with no called
    // sample in a population gets 0.0, and adds nothing to the mean or to the variance.
    let site_freqs: Vec<Vec<f32>> = ac
        .iter()
        .zip(&an)
        .map(|(a, n)| {
            a.iter()
                .zip(n)
                .map(|(&c, &t)| if t == 0 { 0.0 } else { c as f32 / t as f32 })
                .collect()
        })
        .collect();

    // The LD variance inflation: how much wider the cohort's real spread is than a model of
    // independent sites predicts.
    //
    // An archaic allele travels in a linked block, so the binomial sum states too small a variance.
    // Without this correction, every percentile moves toward an extreme, and an ordinary person
    // looks remarkable. The code measures it on the full panel, and applies it to any subset.
    let mut variance_inflation: Vec<Vec<(f32, f32)>> = Vec::with_capacity(supers.len());
    for (pi, sup) in supers.iter().enumerate() {
        let mut ladder: Vec<(f32, f32)> = Vec::new();
        for (r, rung) in rung_tot.iter().enumerate() {
            // Predicted (independent-sites) variance over just this rung's sites.
            let predicted: f64 = site_freqs[pi]
                .iter()
                .enumerate()
                .filter(|(si, _)| in_rung(*si, r))
                .map(|(_, &f)| 2.0 * f as f64 * (1.0 - f as f64))
                .sum();
            let vals: Vec<f64> = samples
                .iter()
                .enumerate()
                .filter(|(i, _)| super_of[*i].as_deref() == Some(sup.as_str()) && scored[*i] > 0)
                .map(|(i, _)| rung[i] as f64)
                .collect();
            if vals.len() < 2 || predicted <= 0.0 {
                continue;
            }
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let observed = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
            // Actual density realised by the hash, not the nominal target.
            let realised = site_freqs[pi]
                .iter()
                .enumerate()
                .filter(|(si, _)| in_rung(*si, r))
                .count() as f32
                / panel.len().max(1) as f32;
            ladder.push((realised, (observed / predicted).max(1.0) as f32));
        }
        // A real site set: the same computation, but over the real mask, and not a hash rung.
        for (m, (label, mask)) in subset_masks.iter().enumerate() {
            let predicted: f64 = site_freqs[pi]
                .iter()
                .enumerate()
                .filter(|(si, _)| mask[*si])
                .map(|(_, &f)| 2.0 * f as f64 * (1.0 - f as f64))
                .sum();
            let vals: Vec<f64> = samples
                .iter()
                .enumerate()
                .filter(|(i, _)| super_of[*i].as_deref() == Some(sup.as_str()) && scored[*i] > 0)
                .map(|(i, _)| subset_tot[m][i] as f64)
                .collect();
            if vals.len() < 2 || predicted <= 0.0 {
                continue;
            }
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let observed = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
            let realised = mask.iter().filter(|b| **b).count() as f32 / panel.len().max(1) as f32;
            let inflation = (observed / predicted).max(1.0) as f32;
            eprintln!(
                "  {sup:<5} real subset {label}: density {:.2}% inflation {inflation:.1}x",
                realised * 100.0
            );
            ladder.push((realised, inflation));
        }
        // A real subset TAKES THE PLACE OF the synthetic rungs, and does not join them. With both
        // in one ladder, the ladder no longer climbs in order: EUR measures 5.3x on the real 2.6%
        // subset, but only 2.0x on a random 3% one.
        //
        // An interpolation across both would then *lower* the inflation as the density rises, and
        // give a chip user an over-confident percentile. So keep the full-panel rung, which is real
        // by construction, plus every subset that the code measured.
        if !subset_masks.is_empty() {
            let full = ladder.first().copied();
            let mut kept: Vec<(f32, f32)> = ladder.drain(DENSITY_LADDER.len().min(ladder.len())..).collect();
            if let Some(f) = full {
                kept.push(f);
            }
            ladder = kept;
        }
        ladder.sort_by(|a, b| a.0.total_cmp(&b.0));
        let shown: Vec<String> = ladder
            .iter()
            .map(|(d, i)| format!("{:.1}%:{:.1}x", d * 100.0, i))
            .collect();
        eprintln!("  {sup:<5} variance inflation by density  {}", shown.join("  "));
        variance_inflation.push(ladder);
    }

    // Use the shared hash helper, and not a local sha2 implementation. The codebase brought every
    // sha256 onto du_bio::hash on purpose.
    let panel_fingerprint = navigator_analysis::manifest::sha256_hex(&bytes);

    let dist = ArchaicCountDistribution {
        build: panel.build.clone(),
        panel_sites: panel.len(),
        cohorts,
        populations: supers,
        site_freqs,
        variance_inflation,
        panel_fingerprint,
    };
    let bytes = dist.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    write_bin(&args.out, &bytes)?;
    eprintln!("wrote {}", args.out.display());
    Ok(())
}

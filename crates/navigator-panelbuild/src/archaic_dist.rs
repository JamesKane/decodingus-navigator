//! Build the archaic percentile reference (`archaic_marker_dist_<build>.bin`) — Asset 4 of
//! `documents/design/ArchaicAncestry_Design.md` §4.
//!
//! Tier A reports a **count**, and a bare count means nothing without a cohort to place it against
//! (design §1: 23andMe pairs its count with a percentile for exactly this reason). This builder
//! scores every reference sample through the same `count_archaic_markers` arithmetic the app will
//! run on a subject, so the percentile compares like with like.
//!
//! **Counts are stored per fine population, not pre-reduced to super-populations** (design §9 Q3).
//! v1 renders the percentile against a super-population — keying it to the user's inferred fine
//! ancestry would let an ancestry error silently move the archaic headline — but keeping the
//! fine-grained counts means a fine-pop percentile later is a re-keying, not an asset rebuild.
//!
//! Input is a `bcftools query` matrix at the panel sites plus the 1kGP `sample/pop/super_pop`
//! panel, mirroring how every other stage feeds this crate.

use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::archaic::{ArchaicCountDistribution, ArchaicMarkerPanel, CohortCounts};

use crate::pca::{open_maybe_gz, write_bin};

#[derive(Parser)]
pub struct ArchaicDistArgs {
    /// The built marker panel (`archaic_markers_<build>.bin`). The distribution is only meaningful
    /// against the panel it was scored with, so its site count is recorded in the output.
    #[arg(long)]
    panel: PathBuf,
    /// Genotype matrix at the panel sites: `CHROM POS REF ALT [GT...]` from `bcftools query`.
    #[arg(long)]
    matrix: PathBuf,
    /// Sample ids, one per line, in the matrix's column order (`bcftools query -l`).
    #[arg(long)]
    samples: PathBuf,
    /// 1kGP-style panel: `sample<TAB>pop<TAB>super_pop[...]`, with a header line. Samples absent
    /// here are dropped — which is also how related samples are excluded, since the unrelated-2504
    /// panel is the one that carries labels.
    #[arg(long)]
    pops: PathBuf,
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
        // Skip a header row without needing to know its exact spelling.
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

    // Derived allele per site, keyed by position — the matrix and the panel agree on CHM13
    // coordinates, but the matrix states its own REF/ALT, so the derived BASE is what we match on
    // (a dosage would invert wherever the panel's orientation put the derived allele on REF).
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

    let mut totals = vec![0u32; samples.len()];
    let mut scored = vec![0u32; samples.len()];
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
        // Which allele index carries the derived base? A site whose alleles match neither is a
        // matrix/panel disagreement and is skipped rather than guessed.
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
            }
        }
    }
    eprintln!("matrix rows {rows}, matched to panel {matched}, skipped {skipped}");
    anyhow::ensure!(
        matched > 0,
        "no matrix row matched the panel — wrong build or wrong matrix?"
    );

    // Group per fine population, keeping counts sorted so a percentile is a cheap scan.
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

    let dist = ArchaicCountDistribution {
        build: panel.build.clone(),
        panel_sites: panel.len(),
        cohorts,
    };
    let bytes = dist.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    write_bin(&args.out, &bytes)?;
    eprintln!("wrote {}", args.out.display());
    Ok(())
}

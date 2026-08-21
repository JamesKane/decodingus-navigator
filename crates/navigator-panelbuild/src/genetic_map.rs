//! Build the IBD genetic-map asset (`genetic_map_<build>.bin`) from a recombination map.
//!
//! The length of an IBD segment in cM, and so the relationship band that the app reports, is only
//! as good as the map. The app has a flat 1 cM/Mb stand-in. This replaces it with a real
//! sex-averaged map, from deCODE 2019 or HapMap II, that somebody has **already lifted to CHM13**.
//!
//! That lift moves coordinates alone, and no alleles, so a stage-2 CrossMap BED lift of the map's
//! positions is enough. This step parses the lifted text and serializes it to the bincode
//! [`navigator_analysis::ibd::GeneticMap`] that the app loads.
//!
//! Whitespace or a tab separates the input columns, which are
//! `chromosome  position(bp)  …  cumulative_cM`. The **last** column holds the cumulative genetic
//! position. That matches HapMap's `Chromosome Position(bp) Rate(cM/Mb) Map(cM)`, and the simple
//! `chrom pos cM` form.
//!
//! A first data row that is not numeric counts as a header, and the parser skips it. Every position
//! must be a CHM13 coordinate.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ibd::GeneticMap;

#[derive(Parser)]
pub struct GeneticMapArgs {
    /// Recombination-map text (CHM13 coordinates): `chromosome  position(bp)  …  cumulative_cM`.
    #[arg(long)]
    pub input: PathBuf,
    /// Output asset (bincode), e.g. `~/.decodingus/ancestry/genetic_map_chm13v2.0.bin`.
    #[arg(long)]
    pub out: PathBuf,
}

pub fn build_genetic_map(args: GeneticMapArgs) -> Result<()> {
    let file = File::open(&args.input).with_context(|| format!("open {}", args.input.display()))?;
    let mut by_chrom: BTreeMap<String, Vec<(i32, f64)>> = BTreeMap::new();
    let mut parsed = 0usize;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        // Column 1 holds the position, and the last column holds the cumulative cM. The parser
        // skips a header row, which is not numeric.
        let (Ok(pos), Ok(cm)) = (f[1].parse::<i32>(), f[f.len() - 1].parse::<f64>()) else {
            continue;
        };
        by_chrom.entry(f[0].to_string()).or_default().push((pos, cm));
        parsed += 1;
    }
    if parsed == 0 {
        anyhow::bail!("no (chrom, pos, cM) rows parsed from {}", args.input.display());
    }
    let n_chrom = by_chrom.len();

    let markers = by_chrom.into_iter().map(|(chrom, mut rows)| {
        rows.sort_by_key(|(p, _)| *p);
        rows.dedup_by_key(|(p, _)| *p); // collapse duplicate positions (keep first after sort)
        let positions = rows.iter().map(|(p, _)| *p).collect::<Vec<_>>();
        let cm = rows.iter().map(|(_, c)| *c).collect::<Vec<_>>();
        (chrom, positions, cm)
    });
    let map = GeneticMap::from_markers(markers);
    let bytes = map.to_bytes().map_err(|e| anyhow::anyhow!(e.to_string()))?;
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).ok();
    }
    File::create(&args.out)?.write_all(&bytes)?;
    eprintln!(
        "wrote {} ({parsed} markers across {n_chrom} chromosomes)",
        args.out.display()
    );
    Ok(())
}

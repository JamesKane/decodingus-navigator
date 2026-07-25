//! Numeric validation gates for the copying (chromosome-painter) local-ancestry model —
//! [`navigator_analysis::lai::paint_copying_lai`].
//!
//! The painter's calibration knobs ([`CopyingLaiParams`]) were tuned by *looking* at one kit's
//! painted chromosomes, which cannot distinguish "the smear is gone" from "the smear moved". This
//! runs the shipping painter against ground truth we control and prints the numbers:
//!
//! 1. **Leave-one-out gate** — take a real reference individual (both of its haplotypes), remove
//!    them from the reference so it cannot copy itself, paint, and score every site against the
//!    individual's known population. A NW-European reference individual painted as Finnish is the
//!    exact defect the recent recalibration commits were chasing, and here it is a number.
//! 2. **Simulated-admixture gate** — splice held-out donor haplotypes from two (or more)
//!    populations into a mosaic with exponential tract lengths for `g` generations, so the truth is
//!    known *per site* including the breakpoints. Scores accuracy and composition error.
//!
//! Scoring is **phase-insensitive**: at each site the two called labels are matched to the two true
//! labels as unordered multisets, so a phase switch between the sides is not counted as an ancestry
//! error (phasing accuracy is a separate concern, exercised with `--phase`). Sites are weighted by
//! the base-pair span they represent, so the reported composition matches what the painter draws.
//!
//! With `--sweep` the whole case set is re-run over a grid of knob values, one row per combination —
//! the tuning loop the visual check was standing in for. Example:
//!
//! ```text
//! navigator-panelbuild validate-lai --replicates 3 \
//!   --sweep "recomb_per_cm=0.05,0.1,0.3;max_ref_haps=30,50,100" --tsv sweep.tsv
//! ```
//!
//! Caveat worth keeping in mind when reading the output: the reference panel unions 1000G and HGDP,
//! which carry **duplicate labels** for the same population (`SRD`/`Sardinian`, `FRN`/`French`, …).
//! A call landing on the sibling label is a labelling artefact, not an ancestry error, so the report
//! carries a "group" accuracy column that merges the known duplicates alongside the raw fine one.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ancestry::HaplotypeReference;
use navigator_analysis::caller::SiteGenotype;
use navigator_analysis::ibd::GeneticMap;
use navigator_analysis::lai::{paint_copying_lai, CopyingLaiParams};
use navigator_analysis::phasing::{PhaseParams, PhasedGenotypes, PhasedSite, Phaser, ReferencePhaser};
use navigator_domain::ancestry::{population_super, AncestrySegment};
use rayon::prelude::*;

#[derive(Parser)]
pub struct LaiValidateArgs {
    /// Phased-haplotype reference to validate (`ancestry_haps_<build>.bin`).
    #[arg(long, default_value = "~/.decodingus/ancestry/ancestry_haps_chm13v2.0.bin")]
    haps: PathBuf,
    /// Genetic map asset (`genetic_map_<build>.bin`). Falls back to a uniform 1 cM/Mb map, which
    /// changes the answers — the copying model's switch costs are in cM.
    #[arg(long, default_value = "~/.decodingus/ancestry/genetic_map_chm13v2.0.bin")]
    map: PathBuf,
    /// Populations to hold out and paint, comma-separated; `all` uses every population with at
    /// least `--min-individuals` individuals in the reference.
    #[arg(long, default_value = "GBR,CEU,FIN,TSI,IBS,YRI,CHB,PJL")]
    pops: String,
    /// Individuals held out and painted per population.
    #[arg(long, default_value_t = 3)]
    replicates: usize,
    /// A population needs this many individuals to be testable (one is always held out).
    #[arg(long, default_value_t = 10)]
    min_individuals: usize,
    /// Restrict to these contigs (comma-separated, e.g. `chr1,chr2`). Default: the whole panel.
    #[arg(long)]
    contigs: Option<String>,
    /// Keep only every N-th reference site before painting — the marker-density knob. Sweep it
    /// (`--sweep "thin=1,2,4,8"` is not a painter knob, so pass repeated runs) to measure how much
    /// of the painter's accuracy is density-limited rather than model-limited.
    #[arg(long, default_value_t = 1)]
    thin: usize,
    /// Simulated admixture cases, comma-separated `POP+POP[+POP]` (e.g. `GBR+YRI,CEU+CHB`).
    /// Sources are mixed in equal proportions; empty disables the gate.
    #[arg(long, default_value = "")]
    admixed: String,
    /// Generations since admixture for the simulated mosaics (higher → shorter true tracts).
    #[arg(long, default_value_t = 8)]
    generations: usize,
    /// Donor individuals per source population in an admixture case (all of them are held out).
    #[arg(long, default_value_t = 4)]
    donors: usize,
    /// Run the production phaser first (statistical phasing of the collapsed genotypes) instead of
    /// feeding the true haplotype sides. Measures the end-to-end path; much slower.
    #[arg(long)]
    phase: bool,
    /// Genome-wide composition prior handed to the painter — the app passes the `estimate_admixture`
    /// result, which for a single-population sample is its own continent. `truth` mimics that;
    /// `flat` spreads the prior evenly over all super-populations and `none` passes an empty prior,
    /// both of which disable the gate's help and so measure whether the *copying model itself* gets
    /// the continent right (with `truth`, super% is ≈100% by construction).
    #[arg(long, default_value = "truth")]
    prior: PriorMode,
    /// Populations counted as drifted isolates for the over-call metric (the failure mode the
    /// recalibration commits were chasing).
    #[arg(
        long,
        default_value = "FIN,RUS,SRD,BSQ,ORC,Russian,Sardinian,Basque,Orcadian,Adygei,Kalash"
    )]
    isolates: String,
    /// Grid of knob values to sweep, `name=v1,v2;name2=v3,v4`. Names: mismatch, recomb_per_cm,
    /// switch_per_cm, min_ref_haps, max_ref_haps, min_segment_cm, min_ancestry.
    #[arg(long)]
    sweep: Option<String>,
    /// Write the per-case (or per-sweep-row) numbers here as TSV.
    #[arg(long)]
    tsv: Option<PathBuf>,
    /// RNG seed — a run is fully reproducible given the seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    // ── Painter knobs (default = the shipping `CopyingLaiParams::default()`) ──
    #[arg(long)]
    mismatch: Option<f64>,
    #[arg(long)]
    recomb_per_cm: Option<f64>,
    #[arg(long)]
    switch_per_cm: Option<f64>,
    #[arg(long)]
    min_ref_haps: Option<usize>,
    #[arg(long)]
    max_ref_haps: Option<usize>,
    #[arg(long)]
    min_segment_cm: Option<f64>,
    #[arg(long)]
    min_ancestry: Option<f64>,
    #[arg(long)]
    size_normalize: Option<f64>,
}

/// Where the painter's genome-wide composition prior comes from in a validation run.
#[derive(Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum PriorMode {
    /// The case's true super-population mix (what a working `estimate_admixture` would return).
    Truth,
    /// Even weight on every super-population the reference carries.
    Flat,
    /// No prior at all — the global-composition gate is disabled.
    None,
}

impl LaiValidateArgs {
    /// The painter parameters this run starts from: the shipping defaults with any CLI override
    /// applied (so a validation run always describes a configuration the app could actually have).
    fn base_params(&self) -> CopyingLaiParams {
        let d = CopyingLaiParams::default();
        CopyingLaiParams {
            mismatch: self.mismatch.unwrap_or(d.mismatch),
            recomb_per_cm: self.recomb_per_cm.unwrap_or(d.recomb_per_cm),
            switch_per_cm: self.switch_per_cm.unwrap_or(d.switch_per_cm),
            min_ref_haps: self.min_ref_haps.unwrap_or(d.min_ref_haps),
            max_ref_haps: self.max_ref_haps.unwrap_or(d.max_ref_haps),
            min_segment_cm: self.min_segment_cm.unwrap_or(d.min_segment_cm),
            min_ancestry: self.min_ancestry.unwrap_or(d.min_ancestry),
            size_normalize: self.size_normalize.unwrap_or(d.size_normalize),
        }
    }
}

/// Deterministic, dependency-free PRNG (xorshift64*) — the same one the ancient-panel validator
/// uses. Reproducibility matters more than statistical quality: a validation run must give the same
/// numbers twice.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }
}

/// 1000G and HGDP label the same population differently; a call landing on the sibling label is a
/// labelling artefact of the unioned panel, not an ancestry error. Reported as a separate "group"
/// accuracy so both readings are visible. Maps a code to its canonical group.
const GROUP_ALIASES: [(&str, &str); 12] = [
    ("Sardinian", "SRD"),
    ("Basque", "BSQ"),
    ("French", "FRN"),
    ("Orcadian", "ORC"),
    ("Russian", "RUS"),
    ("Tuscan", "TSI"),
    ("Bergamo", "ITN"),
    ("Italian", "ITN"),
    ("Han", "CHB"),
    ("Japanese", "JPT"),
    ("Yoruba", "YRI"),
    ("Mandenka", "GWD"),
];

/// Regional clusters — the granularity between a fine population and a whole continent. Fine calls
/// scatter across the populations *within* a cluster (a real British genome comes back CEU/GBR/FRN
/// in varying mixtures), so the cluster is the level at which a call may actually be reportable;
/// `region%` measures whether that is true. Codes with no defensible cluster fall back to their
/// super-population.
const REGIONS: [(&str, &str); 20] = [
    // NW Europe
    ("GBR", "NWE"),
    ("CEU", "NWE"),
    ("FRN", "NWE"),
    ("French", "NWE"),
    ("ORC", "NWE"),
    ("Orcadian", "NWE"),
    // South / SW Europe
    ("IBS", "SEU"),
    ("TSI", "SEU"),
    ("Tuscan", "SEU"),
    ("SRD", "SEU"),
    ("Sardinian", "SEU"),
    ("ITN", "SEU"),
    ("Bergamo", "SEU"),
    ("Italian", "SEU"),
    ("BSQ", "SEU"),
    ("Basque", "SEU"),
    // NE Europe / Baltic-Finnic
    ("FIN", "NEE"),
    ("RUS", "NEE"),
    ("Russian", "NEE"),
    ("Adygei", "CAU"),
];

fn region_of(code: &str) -> &str {
    REGIONS
        .iter()
        .find(|(from, _)| *from == code)
        .map(|(_, to)| *to)
        .unwrap_or_else(|| super_of(code))
}

fn group_of(code: &str) -> &str {
    GROUP_ALIASES
        .iter()
        .find(|(from, _)| *from == code)
        .map(|(_, to)| *to)
        .unwrap_or(code)
}

fn super_of(code: &str) -> &str {
    population_super(code).unwrap_or(code)
}

/// One reference individual: its two haplotype rows and population.
#[derive(Clone, Copy)]
struct Individual {
    hap0: usize,
    hap1: usize,
    pop: usize,
}

/// A validation case: the haplotypes to withhold from the reference, the two sides to paint, and
/// the true population label at every selected site on each side.
struct Case {
    label: String,
    kind: &'static str,
    hold_out: Vec<usize>,
    sides: [Vec<u8>; 2],
    truth: [Vec<String>; 2],
}

/// What one painted case scored. All shares are base-pair weighted over the selected sites.
struct CaseMetrics {
    label: String,
    kind: &'static str,
    truth_composition: BTreeMap<String, f64>,
    fine_acc: f64,
    group_acc: f64,
    region_acc: f64,
    super_acc: f64,
    /// Called fine-label composition (over called sites), for the mis-call breakdown.
    called: BTreeMap<String, f64>,
    /// Accuracy a painter that picked uniformly among the labels it actually used would score —
    /// the line `fine%` has to clear to mean anything. Fine accuracy near it says the copying model
    /// is not resolving sub-populations at all, however clean the segments look.
    chance: f64,
    /// Weight fraction of sides covered by any segment (should be ≈1).
    covered: f64,
    n_segments: usize,
    mean_segment_mb: f64,
    /// Mean absolute error of the called super-population composition vs the truth, in points.
    composition_err: f64,
    /// Is the true population the single largest called label? This — not per-site accuracy — is
    /// what the UI's per-side population list claims, so it is scored separately: a configuration
    /// can raise per-site accuracy while pushing the true population out of the top of the list.
    top1: f64,
    /// …and is it anywhere in the top three called labels?
    top3: f64,
}

pub fn validate_lai(args: LaiValidateArgs) -> Result<()> {
    let reference = load_reference(&expand(&args.haps))?.thin_sites(args.thin);
    if args.thin > 1 {
        println!("thinned to every {}th site", args.thin);
    }
    println!(
        "reference: {} haplotypes × {} sites, {} populations, build {}",
        reference.n_haplotypes,
        reference.n_sites,
        reference.populations.len(),
        reference.build
    );

    let sel = select_sites(&reference, args.contigs.as_deref());
    anyhow::ensure!(!sel.is_empty(), "no reference sites match --contigs");
    let weights = site_weights(&reference, &sel);
    let map = load_map(&expand(&args.map), &reference);
    println!(
        "painting {} sites over {} contigs{}\n",
        sel.len(),
        contig_count(&reference, &sel),
        if args.phase {
            " (through the production phaser)"
        } else {
            ""
        }
    );

    let individuals = individuals(&reference);
    let cases = build_cases(&args, &reference, &individuals, &sel, &map)?;
    anyhow::ensure!(!cases.is_empty(), "no cases to run — check --pops / --admixed");
    println!("{} cases: {}\n", cases.len(), case_summary(&cases));

    let isolates: HashSet<&str> = args
        .isolates
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();

    match &args.sweep {
        None => {
            let params = args.base_params();
            let metrics = run_all(
                &cases, &reference, &map, &sel, &weights, &params, args.phase, args.prior,
            );
            report_detail(&metrics, &isolates);
            if let Some(path) = &args.tsv {
                write_case_tsv(path, &metrics)?;
            }
        }
        Some(spec) => {
            let grid = expand_grid(&args.base_params(), spec)?;
            println!(
                "sweeping {} parameter combinations × {} cases\n",
                grid.len(),
                cases.len()
            );
            let mut rows = Vec::new();
            for params in &grid {
                let metrics = run_all(&cases, &reference, &map, &sel, &weights, params, args.phase, args.prior);
                rows.push((params.clone(), summarize(&metrics, &isolates)));
            }
            report_sweep(&rows);
            if let Some(path) = &args.tsv {
                write_sweep_tsv(path, &rows)?;
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_all(
    cases: &[Case],
    reference: &HaplotypeReference,
    map: &GeneticMap,
    sel: &[usize],
    weights: &[f64],
    params: &CopyingLaiParams,
    phase: bool,
    prior: PriorMode,
) -> Vec<CaseMetrics> {
    cases
        .par_iter()
        .map(|case| run_case(case, reference, map, sel, weights, params, phase, prior))
        .collect()
}

/// The genome-wide composition prior the app would hand the painter, per [`PriorMode`].
fn build_prior(case: &Case, reference: &HaplotypeReference, mode: PriorMode) -> Vec<(String, f64)> {
    let mut prior: BTreeMap<String, f64> = BTreeMap::new();
    match mode {
        PriorMode::None => return Vec::new(),
        PriorMode::Truth => {
            for side in &case.truth {
                for code in side {
                    *prior.entry(super_of(code).to_string()).or_default() += 1.0;
                }
            }
        }
        PriorMode::Flat => {
            for pop in &reference.populations {
                prior.entry(super_of(pop).to_string()).or_insert(1.0);
            }
        }
    }
    let total: f64 = prior.values().sum::<f64>().max(1.0);
    prior.into_iter().map(|(k, v)| (k, v / total)).collect()
}

/// Paint one case against the reference with the case's own haplotypes withheld, and score it.
#[allow(clippy::too_many_arguments)]
fn run_case(
    case: &Case,
    reference: &HaplotypeReference,
    map: &GeneticMap,
    sel: &[usize],
    weights: &[f64],
    params: &CopyingLaiParams,
    phase: bool,
    prior_mode: PriorMode,
) -> CaseMetrics {
    let reduced = reference.without_haplotypes(&case.hold_out);
    let prior = build_prior(case, reference, prior_mode);

    let phased = if phase {
        let genotypes: Vec<SiteGenotype> = sel
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                let site = &reference.sites[c];
                SiteGenotype {
                    name: String::new(),
                    contig: site.contig.clone(),
                    position: site.position,
                    reference_allele: site.reference_allele.to_string(),
                    alternate_allele: site.alternate_allele.to_string(),
                    ploidy: 2,
                    dosage: (case.sides[0][i] + case.sides[1][i]) as i32,
                    gq: 99,
                    depth: 30,
                    ref_depth: 15,
                    alt_depth: 15,
                    pls: Vec::new(),
                    gt: None,
                    allele_depths: None,
                }
            })
            .collect();
        ReferencePhaser::new(&reduced, map, PhaseParams::default()).phase(&genotypes)
    } else {
        PhasedGenotypes {
            sites: sel
                .iter()
                .enumerate()
                .map(|(i, &c)| PhasedSite {
                    contig: reference.sites[c].contig.clone(),
                    position: reference.sites[c].position,
                    side0: case.sides[0][i],
                    side1: case.sides[1][i],
                    confidence: 1.0,
                })
                .collect(),
        }
    };

    let segments = paint_copying_lai(&phased, &reduced, map, &prior, params);
    score(case, &segments, reference, sel, weights)
}

/// Score painted segments against the case's per-site truth. Phase-insensitive: the two called
/// labels at a site are matched to the two true labels as unordered multisets.
fn score(
    case: &Case,
    segments: &[AncestrySegment],
    reference: &HaplotypeReference,
    sel: &[usize],
    weights: &[f64],
) -> CaseMetrics {
    // Positions of the selected sites per contig, with their index into `sel` — for mapping a
    // segment's [start, end] back onto sites.
    let mut by_contig: HashMap<&str, Vec<(i64, usize)>> = HashMap::new();
    for (i, &c) in sel.iter().enumerate() {
        by_contig
            .entry(reference.sites[c].contig.as_str())
            .or_default()
            .push((reference.sites[c].position, i));
    }
    for v in by_contig.values_mut() {
        v.sort_unstable();
    }

    let mut called: [Vec<Option<&str>>; 2] = [vec![None; sel.len()], vec![None; sel.len()]];
    let mut seg_span_mb = 0.0f64;
    for seg in segments {
        let Some(sites) = by_contig.get(seg.contig.as_str()) else {
            continue;
        };
        let label = seg
            .fine_population_code
            .as_deref()
            .unwrap_or(seg.population_code.as_str());
        let lo = sites.partition_point(|(p, _)| *p < seg.start);
        let hi = sites.partition_point(|(p, _)| *p <= seg.end);
        for &(_, i) in &sites[lo..hi] {
            called[seg.copy.min(1) as usize][i] = Some(label);
        }
        seg_span_mb += (seg.end - seg.start).max(0) as f64 / 1e6;
    }

    let (mut fine_hit, mut group_hit, mut region_hit, mut super_hit, mut denom, mut covered) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut called_comp: BTreeMap<String, f64> = BTreeMap::new();
    let mut truth_comp: BTreeMap<String, f64> = BTreeMap::new();
    let mut called_super: BTreeMap<String, f64> = BTreeMap::new();
    let mut truth_super: BTreeMap<String, f64> = BTreeMap::new();
    for i in 0..sel.len() {
        let w = weights[i];
        denom += 2.0 * w;
        let c = [called[0][i], called[1][i]];
        let t = [case.truth[0][i].as_str(), case.truth[1][i].as_str()];
        fine_hit += w * unordered_hits(c, t, |x| x);
        group_hit += w * unordered_hits(c, t, group_of);
        region_hit += w * unordered_hits(c, t, region_of);
        super_hit += w * unordered_hits(c, t, super_of);
        for code in t {
            *truth_comp.entry(code.to_string()).or_default() += w;
            *truth_super.entry(super_of(code).to_string()).or_default() += w;
        }
        for code in c.into_iter().flatten() {
            covered += w;
            *called_comp.entry(code.to_string()).or_default() += w;
            *called_super.entry(super_of(code).to_string()).or_default() += w;
        }
    }

    let norm = |m: BTreeMap<String, f64>, total: f64| -> BTreeMap<String, f64> {
        m.into_iter()
            .map(|(k, v)| (k, if total > 0.0 { v / total } else { 0.0 }))
            .collect()
    };
    let called_total: f64 = called_comp.values().sum();
    let called_super = norm(called_super, called_total);
    let truth_super = norm(truth_super, denom);
    let mut keys: HashSet<&String> = called_super.keys().collect();
    keys.extend(truth_super.keys());
    let composition_err = keys
        .iter()
        .map(|k| (called_super.get(*k).copied().unwrap_or(0.0) - truth_super.get(*k).copied().unwrap_or(0.0)).abs())
        .sum::<f64>()
        / 2.0
        * 100.0;

    // Top-1 / top-3 containment, weighted by each truth label's share of the genome.
    let mut ranked: Vec<(&String, f64)> = called_comp.iter().map(|(k, v)| (k, *v)).collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    let rank_of = |code: &str| ranked.iter().position(|(k, _)| k.as_str() == code);
    let truth_total: f64 = truth_comp.values().sum();
    let (mut top1, mut top3) = (0.0, 0.0);
    for (code, w) in &truth_comp {
        let share = if truth_total > 0.0 { w / truth_total } else { 0.0 };
        match rank_of(code) {
            Some(0) => {
                top1 += share;
                top3 += share;
            }
            Some(r) if r < 3 => top3 += share,
            _ => {}
        }
    }
    let chance = if called_comp.is_empty() {
        0.0
    } else {
        1.0 / called_comp.len() as f64
    };
    CaseMetrics {
        label: case.label.clone(),
        kind: case.kind,
        chance,
        truth_composition: norm(truth_comp, denom),
        fine_acc: fine_hit / denom,
        group_acc: group_hit / denom,
        region_acc: region_hit / denom,
        super_acc: super_hit / denom,
        called: norm(called_comp, called_total),
        covered: covered / denom,
        top1,
        top3,
        n_segments: segments.len(),
        mean_segment_mb: if segments.is_empty() {
            0.0
        } else {
            seg_span_mb / segments.len() as f64
        },
        composition_err,
    }
}

/// How many of the two called labels match the two true labels as unordered multisets (0, 1 or 2),
/// after mapping both through `canon` (identity for fine, group/super roll-up otherwise).
fn unordered_hits(called: [Option<&str>; 2], truth: [&str; 2], canon: fn(&str) -> &str) -> f64 {
    let mut remaining = vec![canon(truth[0]), canon(truth[1])];
    let mut hits = 0.0;
    for maybe in called.into_iter().flatten() {
        let c = canon(maybe);
        if let Some(p) = remaining.iter().position(|t| *t == c) {
            remaining.remove(p);
            hits += 1.0;
        }
    }
    hits
}

// ── Case construction ────────────────────────────────────────────────────────────────────────

/// Reference haplotypes come in per-sample consecutive pairs (the builder pushes both sides of each
/// sample in order); pair them back into individuals, keeping only pairs that agree on population.
fn individuals(reference: &HaplotypeReference) -> Vec<Individual> {
    (0..reference.n_haplotypes / 2)
        .filter_map(|i| {
            let (h0, h1) = (2 * i, 2 * i + 1);
            (reference.hap_pop[h0] == reference.hap_pop[h1]).then_some(Individual {
                hap0: h0,
                hap1: h1,
                pop: reference.hap_pop[h0] as usize,
            })
        })
        .collect()
}

fn build_cases(
    args: &LaiValidateArgs,
    reference: &HaplotypeReference,
    individuals: &[Individual],
    sel: &[usize],
    map: &GeneticMap,
) -> Result<Vec<Case>> {
    // Individuals per population code.
    let mut by_pop: BTreeMap<&str, Vec<Individual>> = BTreeMap::new();
    for ind in individuals {
        by_pop
            .entry(reference.populations[ind.pop].as_str())
            .or_default()
            .push(*ind);
    }

    let wanted: Vec<String> = if args.pops.trim() == "all" {
        by_pop
            .iter()
            .filter(|(_, v)| v.len() >= args.min_individuals)
            .map(|(k, _)| k.to_string())
            .collect()
    } else {
        args.pops
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    };

    let mut cases = Vec::new();
    for pop in &wanted {
        let Some(pool) = by_pop.get(pop.as_str()) else {
            println!("  {pop}: not in the reference — skipped");
            continue;
        };
        if pool.len() < args.min_individuals {
            println!(
                "  {pop}: only {} individuals (< --min-individuals) — skipped",
                pool.len()
            );
            continue;
        }
        // Evenly spaced picks so replicates aren't all neighbours in the panel's sample order.
        let step = (pool.len() / args.replicates.max(1)).max(1);
        for r in 0..args.replicates.min(pool.len()) {
            let ind = pool[(r * step) % pool.len()];
            cases.push(Case {
                label: format!("{pop}#{r}"),
                kind: "held-out",
                hold_out: vec![ind.hap0, ind.hap1],
                sides: [
                    sel.iter().map(|&c| reference.allele(ind.hap0, c)).collect(),
                    sel.iter().map(|&c| reference.allele(ind.hap1, c)).collect(),
                ],
                truth: [vec![pop.clone(); sel.len()], vec![pop.clone(); sel.len()]],
            });
        }
    }

    for spec in args.admixed.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let sources: Vec<&str> = spec.split('+').map(str::trim).collect();
        anyhow::ensure!(
            sources.len() >= 2,
            "--admixed case `{spec}` needs at least two populations"
        );
        let mut pools: Vec<(&str, Vec<Individual>)> = Vec::new();
        for src in &sources {
            let Some(pool) = by_pop.get(src) else {
                println!("  {spec}: {src} not in the reference — skipped");
                pools.clear();
                break;
            };
            // Donors are held out wholesale, so a source needs donors + spare individuals.
            let donors: Vec<Individual> = pool.iter().take(args.donors).copied().collect();
            anyhow::ensure!(!donors.is_empty(), "{src} has no individuals to donate");
            pools.push((src, donors));
        }
        if pools.is_empty() {
            continue;
        }
        cases.push(simulate_admixed(spec, &pools, reference, sel, map, args));
    }
    Ok(cases)
}

/// Splice donor haplotypes into a two-sided mosaic with exponential tract lengths for
/// `generations` generations of recombination (crossover rate `g` per Morgan = `g/100` per cM), so
/// the true ancestry — and its breakpoints — are known at every site.
fn simulate_admixed(
    spec: &str,
    pools: &[(&str, Vec<Individual>)],
    reference: &HaplotypeReference,
    sel: &[usize],
    map: &GeneticMap,
    args: &LaiValidateArgs,
) -> Case {
    let hold_out: Vec<usize> = pools
        .iter()
        .flat_map(|(_, ds)| ds.iter().flat_map(|d| [d.hap0, d.hap1]))
        .collect();
    let mut sides = [vec![0u8; sel.len()], vec![0u8; sel.len()]];
    let mut truth = [vec![String::new(); sel.len()], vec![String::new(); sel.len()]];
    let rate = args.generations as f64 / 100.0; // switches per cM

    for (side, (alleles, labels)) in sides.iter_mut().zip(truth.iter_mut()).enumerate() {
        let mut rng = Rng::new(args.seed ^ ((side as u64 + 1) << 40) ^ spec.bytes().map(u64::from).sum::<u64>());
        let mut source = rng.below(pools.len());
        let mut donor = pick_donor(&mut rng, &pools[source].1);
        let (mut prev_contig, mut prev_pos) = ("", 0i64);
        for (i, &col) in sel.iter().enumerate() {
            let site = &reference.sites[col];
            // New contig → independent assortment: re-draw the source and the donor haplotype.
            let switch = if site.contig != prev_contig {
                true
            } else {
                let d = map
                    .interval_cm(&site.contig, prev_pos as i32, site.position as i32)
                    .unwrap_or(0.0)
                    .max(0.0);
                rng.next_f64() < 1.0 - (-d * rate).exp()
            };
            if switch {
                source = rng.below(pools.len());
                donor = pick_donor(&mut rng, &pools[source].1);
            }
            alleles[i] = reference.allele(donor, col);
            labels[i] = pools[source].0.to_string();
            prev_contig = site.contig.as_str();
            prev_pos = site.position;
        }
    }
    Case {
        label: format!("{spec} g{}", args.generations),
        kind: "admixed",
        hold_out,
        sides,
        truth,
    }
}

fn pick_donor(rng: &mut Rng, donors: &[Individual]) -> usize {
    let ind = donors[rng.below(donors.len())];
    if rng.next_f64() < 0.5 {
        ind.hap0
    } else {
        ind.hap1
    }
}

// ── Reporting ────────────────────────────────────────────────────────────────────────────────

/// Aggregate numbers across all cases — the headline row a sweep compares configurations on.
struct Summary {
    fine: f64,
    top1: f64,
    top3: f64,
    chance: f64,
    group: f64,
    region: f64,
    super_: f64,
    composition_err: f64,
    /// Weighted share of calls that landed on a drifted-isolate population when the truth was not
    /// one — the "Finnish/Sardinian over-call" failure mode, as a number.
    isolate_overcall: f64,
    segments: f64,
    covered: f64,
}

fn summarize(metrics: &[CaseMetrics], isolates: &HashSet<&str>) -> Summary {
    let n = metrics.len().max(1) as f64;
    let mean = |f: fn(&CaseMetrics) -> f64| metrics.iter().map(f).sum::<f64>() / n;
    let mut over = 0.0;
    let mut over_cases = 0.0;
    for m in metrics {
        if m.truth_composition.keys().any(|t| isolates.contains(t.as_str())) {
            continue;
        }
        over += m
            .called
            .iter()
            .filter(|(k, _)| isolates.contains(k.as_str()))
            .map(|(_, v)| v)
            .sum::<f64>();
        over_cases += 1.0;
    }
    Summary {
        fine: mean(|m| m.fine_acc),
        top1: mean(|m| m.top1),
        top3: mean(|m| m.top3),
        chance: mean(|m| m.chance),
        group: mean(|m| m.group_acc),
        region: mean(|m| m.region_acc),
        super_: mean(|m| m.super_acc),
        composition_err: mean(|m| m.composition_err),
        isolate_overcall: if over_cases > 0.0 { over / over_cases } else { 0.0 },
        segments: mean(|m| m.n_segments as f64),
        covered: mean(|m| m.covered),
    }
}

fn report_detail(metrics: &[CaseMetrics], isolates: &HashSet<&str>) {
    println!(
        "{:<16}{:>10}{:>8}{:>9}{:>8}{:>9}{:>8}{:>9}{:>8}{:>9}  top calls",
        "case", "kind", "fine%", "chance%", "group%", "region%", "super%", "compΔ", "segs", "meanMb"
    );
    for m in metrics {
        let mut top: Vec<(&String, &f64)> = m.called.iter().collect();
        top.sort_by(|a, b| b.1.total_cmp(a.1));
        let top: Vec<String> = top
            .iter()
            .take(4)
            .map(|(k, v)| format!("{k} {:.0}%", *v * 100.0))
            .collect();
        println!(
            "{:<16}{:>10}{:>8.1}{:>9.1}{:>8.1}{:>9.1}{:>8.1}{:>9.1}{:>8}{:>9.1}  {}",
            m.label,
            m.kind,
            m.fine_acc * 100.0,
            m.chance * 100.0,
            m.group_acc * 100.0,
            m.region_acc * 100.0,
            m.super_acc * 100.0,
            m.composition_err,
            m.n_segments,
            m.mean_segment_mb,
            top.join(" · ")
        );
    }
    let s = summarize(metrics, isolates);
    println!(
        "\nmean: fine {:.1}% (chance {:.1}%) · TOP-1 {:.1}% · TOP-3 {:.1}% · group {:.1}% · \
         region {:.1}% · super {:.1}% · composition error {:.1} pts · isolate over-call {:.1}% · \
         {:.0} segments · {:.1}% of sites covered",
        s.fine * 100.0,
        s.chance * 100.0,
        s.top1 * 100.0,
        s.top3 * 100.0,
        s.group * 100.0,
        s.region * 100.0,
        s.super_ * 100.0,
        s.composition_err,
        s.isolate_overcall * 100.0,
        s.segments,
        s.covered * 100.0
    );
    println!(
        "\nfine% = called sub-population correct · chance% = what picking uniformly among the labels\n\
         the painter used would score (fine% must clear it to mean anything) · group% merges the\n\
         panel's duplicate 1000G/HGDP labels (SRD/Sardinian, …) · region% rolls Europe up to NW/S/NE\n\
         clusters (the granularity a report could honestly use) · super% = continent correct · compΔ\n\
         = mean absolute error of the super-population composition, in points · isolate over-call =\n\
         share of calls landing on a drifted isolate ({}) for cases that are not one.",
        {
            let mut v: Vec<&str> = isolates.iter().copied().collect();
            v.sort_unstable();
            v.join(",")
        }
    );
}

fn report_sweep(rows: &[(CopyingLaiParams, Summary)]) {
    println!(
        "{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>8}{:>7}{:>8}{:>7}{:>7}{:>9}{:>7}{:>8}{:>9}{:>8}{:>8}",
        "mismatch",
        "recomb",
        "switch",
        "minHaps",
        "maxHaps",
        "minSegCM",
        "minAnc",
        "sizeN",
        "fine%",
        "top1%",
        "top3%",
        "chance%",
        "lift",
        "group%",
        "region%",
        "super%",
        "isol%"
    );
    for (p, s) in rows {
        println!(
            "{:>10.3}{:>10.3}{:>10.3}{:>10}{:>10}{:>10.1}{:>8.2}{:>7.2}{:>8.1}{:>7.1}{:>7.1}{:>9.1}{:>7.1}{:>8.1}{:>9.1}{:>8.1}{:>8.1}",
            p.mismatch,
            p.recomb_per_cm,
            p.switch_per_cm,
            p.min_ref_haps,
            p.max_ref_haps,
            p.min_segment_cm,
            p.min_ancestry,
            p.size_normalize,
            s.fine * 100.0,
            s.top1 * 100.0,
            s.top3 * 100.0,
            s.chance * 100.0,
            (s.fine - s.chance) * 100.0,
            s.group * 100.0,
            s.region * 100.0,
            s.super_ * 100.0,
            s.isolate_overcall * 100.0
        );
    }
    // Ranked on lift over chance, not raw accuracy: a knob that makes the painter use fewer labels
    // raises fine% by shrinking the space it can be wrong in, which is not an improvement.
    if let Some((best, s)) = rows
        .iter()
        .max_by(|a, b| (a.1.fine - a.1.chance).total_cmp(&(b.1.fine - b.1.chance)))
    {
        println!(
            "\nbest lift over chance ({:+.1} pts, isolate over-call {:.1}%) at: recomb_per_cm={} \
             switch_per_cm={} max_ref_haps={} min_ref_haps={} mismatch={} min_segment_cm={} \
             min_ancestry={}",
            (s.fine - s.chance) * 100.0,
            s.isolate_overcall * 100.0,
            best.recomb_per_cm,
            best.switch_per_cm,
            best.max_ref_haps,
            best.min_ref_haps,
            best.mismatch,
            best.min_segment_cm,
            best.min_ancestry
        );
    }
}

fn write_case_tsv(path: &Path, metrics: &[CaseMetrics]) -> Result<()> {
    let mut out = String::from(
        "case\tkind\tfine\tchance\tgroup\tregion\tsuper\tcomposition_err\tsegments\tmean_segment_mb\tcovered\ttop_calls\n",
    );
    for m in metrics {
        let mut top: Vec<(&String, &f64)> = m.called.iter().collect();
        top.sort_by(|a, b| b.1.total_cmp(a.1));
        let top: Vec<String> = top.iter().take(6).map(|(k, v)| format!("{k}={:.4}", *v)).collect();
        out.push_str(&format!(
            "{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}\t{:.2}\t{:.4}\t{}\n",
            m.label,
            m.kind,
            m.fine_acc,
            m.chance,
            m.group_acc,
            m.region_acc,
            m.super_acc,
            m.composition_err,
            m.n_segments,
            m.mean_segment_mb,
            m.covered,
            top.join(";")
        ));
    }
    std::fs::write(path, out).with_context(|| format!("writing {}", path.display()))?;
    println!("\nwrote {}", path.display());
    Ok(())
}

fn write_sweep_tsv(path: &Path, rows: &[(CopyingLaiParams, Summary)]) -> Result<()> {
    let mut out = String::from(
        "mismatch\trecomb_per_cm\tswitch_per_cm\tmin_ref_haps\tmax_ref_haps\tmin_segment_cm\tmin_ancestry\t\
         fine\tchance\tgroup\tregion\tsuper\tcomposition_err\tisolate_overcall\tsegments\n",
    );
    for (p, s) in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{:.4}\t{:.1}\n",
            p.mismatch,
            p.recomb_per_cm,
            p.switch_per_cm,
            p.min_ref_haps,
            p.max_ref_haps,
            p.min_segment_cm,
            p.min_ancestry,
            s.fine,
            s.chance,
            s.group,
            s.region,
            s.super_,
            s.composition_err,
            s.isolate_overcall,
            s.segments
        ));
    }
    std::fs::write(path, out).with_context(|| format!("writing {}", path.display()))?;
    println!("\nwrote {}", path.display());
    Ok(())
}

// ── Inputs ───────────────────────────────────────────────────────────────────────────────────

fn expand(path: &Path) -> PathBuf {
    match path.to_str().and_then(|s| s.strip_prefix("~/")) {
        Some(rest) => match std::env::var("HOME") {
            Ok(home) => PathBuf::from(home).join(rest),
            Err(_) => path.to_path_buf(),
        },
        None => path.to_path_buf(),
    }
}

fn load_reference(path: &Path) -> Result<HaplotypeReference> {
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let reference = HaplotypeReference::from_bytes(&bytes).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))?;
    anyhow::ensure!(!reference.is_empty(), "{} is empty", path.display());
    Ok(reference)
}

/// The real genetic map if present — the copying model's switch costs are in cM, so a uniform
/// fallback quietly changes every number here; say so loudly rather than reporting it as truth.
fn load_map(path: &Path, reference: &HaplotypeReference) -> GeneticMap {
    let mut lengths: BTreeMap<&str, i32> = BTreeMap::new();
    for s in &reference.sites {
        let e = lengths.entry(s.contig.as_str()).or_insert(0);
        *e = (*e).max(s.position as i32);
    }
    let pairs: Vec<(&str, i32)> = lengths.into_iter().collect();
    match std::fs::read(path).ok().and_then(|b| GeneticMap::from_bytes(&b).ok()) {
        Some(m) => m,
        None => {
            println!(
                "!! genetic map {} not readable — falling back to uniform 1 cM/Mb",
                path.display()
            );
            GeneticMap::uniform(1.0, &pairs)
        }
    }
}

fn select_sites(reference: &HaplotypeReference, contigs: Option<&str>) -> Vec<usize> {
    let keep: Option<HashSet<&str>> = contigs.map(|c| c.split(',').map(str::trim).filter(|s| !s.is_empty()).collect());
    (0..reference.n_sites)
        .filter(|&i| {
            keep.as_ref()
                .map_or(true, |k| k.contains(reference.sites[i].contig.as_str()))
        })
        .collect()
}

fn contig_count(reference: &HaplotypeReference, sel: &[usize]) -> usize {
    sel.iter()
        .map(|&i| reference.sites[i].contig.as_str())
        .collect::<HashSet<_>>()
        .len()
}

/// Base-pair weight of each selected site: half the distance to each neighbour within its contig,
/// so a sparse region doesn't count the same as a dense one and the reported composition matches
/// the painted (bp-proportioned) chromosomes.
fn site_weights(reference: &HaplotypeReference, sel: &[usize]) -> Vec<f64> {
    let pos = |i: usize| reference.sites[sel[i]].position as f64;
    let contig = |i: usize| reference.sites[sel[i]].contig.as_str();
    (0..sel.len())
        .map(|i| {
            let left = (i > 0 && contig(i - 1) == contig(i)).then(|| (pos(i) - pos(i - 1)) / 2.0);
            let right = (i + 1 < sel.len() && contig(i + 1) == contig(i)).then(|| (pos(i + 1) - pos(i)) / 2.0);
            match (left, right) {
                (Some(l), Some(r)) => l + r,
                (Some(l), None) => l * 2.0,
                (None, Some(r)) => r * 2.0,
                (None, None) => 1.0, // a lone site on its contig
            }
            .max(1.0)
        })
        .collect()
}

fn case_summary(cases: &[Case]) -> String {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for c in cases {
        *counts.entry(c.kind).or_default() += 1;
    }
    counts
        .iter()
        .map(|(k, v)| format!("{v} {k}"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Expand a `name=v1,v2;name2=v3,v4` sweep spec into the cartesian product of parameter sets,
/// starting from `base` (so unswept knobs keep their configured value).
fn expand_grid(base: &CopyingLaiParams, spec: &str) -> Result<Vec<CopyingLaiParams>> {
    let mut grid = vec![base.clone()];
    for axis in spec.split(';').map(str::trim).filter(|s| !s.is_empty()) {
        let (name, values) = axis
            .split_once('=')
            .ok_or_else(|| anyhow::anyhow!("sweep axis `{axis}` is not `name=v1,v2`"))?;
        let values: Vec<f64> = values
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|v| v.parse::<f64>().with_context(|| format!("sweep value `{v}`")))
            .collect::<Result<_>>()?;
        anyhow::ensure!(!values.is_empty(), "sweep axis `{axis}` has no values");
        let mut next = Vec::with_capacity(grid.len() * values.len());
        for p in &grid {
            for &v in &values {
                let mut p = p.clone();
                apply_knob(&mut p, name.trim(), v)?;
                next.push(p);
            }
        }
        grid = next;
    }
    Ok(grid)
}

fn apply_knob(p: &mut CopyingLaiParams, name: &str, v: f64) -> Result<()> {
    match name {
        "mismatch" => p.mismatch = v,
        "recomb_per_cm" => p.recomb_per_cm = v,
        "switch_per_cm" => p.switch_per_cm = v,
        "min_ref_haps" => p.min_ref_haps = v as usize,
        "max_ref_haps" => p.max_ref_haps = v as usize,
        "min_segment_cm" => p.min_segment_cm = v,
        "min_ancestry" => p.min_ancestry = v,
        "size_normalize" => p.size_normalize = v,
        other => anyhow::bail!("unknown painter knob `{other}`"),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unordered_hits_is_phase_insensitive() {
        // Both sides right, in the other order → still 2.
        assert_eq!(unordered_hits([Some("GBR"), Some("YRI")], ["YRI", "GBR"], |x| x), 2.0);
        // One right, one wrong → 1. An uncalled side scores 0 for that side.
        assert_eq!(unordered_hits([Some("GBR"), Some("CHB")], ["GBR", "YRI"], |x| x), 1.0);
        assert_eq!(unordered_hits([Some("GBR"), None], ["GBR", "GBR"], |x| x), 1.0);
        // The same label twice must not double-count a single true copy.
        assert_eq!(unordered_hits([Some("GBR"), Some("GBR")], ["GBR", "YRI"], |x| x), 1.0);
        // Duplicate panel labels agree at group level, not at fine level.
        assert_eq!(
            unordered_hits([Some("Sardinian"), Some("GBR")], ["SRD", "GBR"], |x| x),
            1.0
        );
        assert_eq!(
            unordered_hits([Some("Sardinian"), Some("GBR")], ["SRD", "GBR"], group_of),
            2.0
        );
        // Wrong sub-population inside the right continent: super level still scores it.
        assert_eq!(
            unordered_hits([Some("FIN"), Some("FIN")], ["GBR", "GBR"], super_of),
            2.0
        );
    }

    #[test]
    fn expand_grid_takes_the_cartesian_product_over_the_base() {
        let base = CopyingLaiParams::default();
        let grid = expand_grid(&base, "recomb_per_cm=0.1,0.3;max_ref_haps=25,50,100").unwrap();
        assert_eq!(grid.len(), 6);
        // Unswept knobs keep the base value.
        assert!(grid.iter().all(|p| p.mismatch == base.mismatch));
        assert!(grid.iter().any(|p| p.recomb_per_cm == 0.3 && p.max_ref_haps == 100));
        assert!(expand_grid(&base, "bogus=1").is_err());
    }
}

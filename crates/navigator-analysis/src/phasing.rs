//! Statistical haplotype phasing — the "which allele came from which parent" step the chromosome
//! painter needs to split ancestry into two internally-consistent parental sides.
//!
//! Navigator analyses **one subject at a time**, so there is no cohort to phase against; we phase
//! against a bundled panel of **phased reference haplotypes** ([`HaplotypeReference`]) using the
//! Li & Stephens copying model — the same reference-based mode EAGLE2/Beagle use. Each of the
//! sample's two haplotypes is modelled as a mosaic of reference haplotypes; the ordered pair of
//! copied reference haplotypes at each site implies the phase.
//!
//! The exact diploid HMM has `K²` states (K = number of reference haplotypes, ~5000), which is
//! infeasible. [`ReferencePhaser`] uses a **beam search** over ordered pair-states, with switch
//! targets restricted to a per-site candidate set of the reference haplotypes sharing the longest
//! IBS run with the sample (a cheap, PBWT-like heuristic). That makes it `O(N · B · M)` in the
//! number of sites `N`, beam width `B`, and candidate count `M`.
//!
//! The [`Phaser`] trait keeps the seam stable so a Mendelian [`TrioPhaser`] (used when a parent
//! sample is in the workspace) or a full PBWT phaser can drop in without touching callers.

use std::collections::HashMap;

use crate::ancestry::HaplotypeReference;
use crate::caller::SiteGenotype;
use crate::ibd::GeneticMap;

/// One phased site: the coordinate plus the allele placed on each of the two sides (`0` = ref
/// allele, `1` = alt). `side0`/`side1` are consistent across the whole chromosome (a genuine
/// parental split), not sorted per site. `confidence` is the phase confidence (1.0 at homozygous
/// sites, which are unambiguous; lower at heterozygous sites the model is unsure about).
#[derive(Debug, Clone, PartialEq)]
pub struct PhasedSite {
    pub contig: String,
    pub position: i64,
    pub side0: u8,
    pub side1: u8,
    pub confidence: f32,
}

/// A sample's phased genotypes: the sites that were both genotyped and present in the reference,
/// each with an allele on side 0 and side 1.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PhasedGenotypes {
    pub sites: Vec<PhasedSite>,
}

/// The phasing strategy. Reference-based statistical phasing by default; a Mendelian trio phaser
/// when a parent is available (see [`TrioPhaser`]).
pub trait Phaser {
    /// Phase the sample's genotypes into two consistent parental sides.
    fn phase(&self, genotypes: &[SiteGenotype]) -> PhasedGenotypes;
}

/// Tuning knobs for [`ReferencePhaser`].
#[derive(Debug, Clone)]
pub struct PhaseParams {
    /// Beam width `B`: ordered pair-states kept per site. Larger → more accurate, slower.
    pub beam: usize,
    /// Per-site switch-target candidates `M`: the reference haplotypes with the longest current
    /// IBS run with the sample that a side may recombine onto.
    pub candidates: usize,
    /// Copying mutation/mismatch rate μ: probability a copied reference allele is observed flipped.
    pub mutation: f64,
    /// Recombination intensity (expected copy switches per centiMorgan). Sets the distance-scaled
    /// switch probability `1 - exp(-d_cM · rate)`.
    pub recomb_per_cm: f64,
}

impl Default for PhaseParams {
    fn default() -> Self {
        Self {
            beam: 128,
            candidates: 64,
            mutation: 1e-3,
            recomb_per_cm: 0.5,
        }
    }
}

/// Reference-based diploid Li & Stephens phaser (beam search over ordered reference-haplotype
/// pairs). Borrows the phased reference and a genetic map (for distance-scaled switch costs).
pub struct ReferencePhaser<'a> {
    reference: &'a HaplotypeReference,
    map: &'a GeneticMap,
    params: PhaseParams,
}

/// A site usable for phasing: reference-hap column index, position, and observed dosage 0/1/2.
struct UsableSite {
    ref_col: usize,
    position: i64,
    dosage: u8,
}

impl<'a> ReferencePhaser<'a> {
    pub fn new(reference: &'a HaplotypeReference, map: &'a GeneticMap, params: PhaseParams) -> Self {
        Self { reference, map, params }
    }

    /// ln P(observe genotype `g` | the two copied reference alleles `c0`, `c1`) under independent
    /// per-copy mutation at rate μ. `c*`/`o*` are 0/1 alleles; `g` is the unordered dosage.
    fn emit_ln(&self, g: u8, c0: u8, c1: u8) -> f64 {
        let mu = self.params.mutation;
        // P(observe o | copied c): faithful with prob 1-μ, flipped with μ.
        let p = |c: u8, o: u8| if c == o { 1.0 - mu } else { mu };
        let prob = match g {
            0 => p(c0, 0) * p(c1, 0),
            2 => p(c0, 1) * p(c1, 1),
            // Heterozygous: one copy emits 0, the other 1 (either assignment).
            _ => p(c0, 0) * p(c1, 1) + p(c0, 1) * p(c1, 0),
        };
        prob.max(1e-300).ln()
    }

    /// Given the MAP copied alleles `(c0, c1)` at a heterozygous site, the ordered `(side0, side1)`
    /// alleles and a confidence. When the copies disagree (one 0, one 1) the phase is determined;
    /// when they agree the het is explained by a mutation and phase is ambiguous (low confidence).
    fn resolve_het(c0: u8, c1: u8) -> (u8, u8, f32) {
        match (c0, c1) {
            (0, 1) => (0, 1, 1.0),
            (1, 0) => (1, 0, 1.0),
            // Both copies same at a het site → a mutation; default ordering, low confidence.
            _ => (0, 1, 0.4),
        }
    }

    /// Phase one contig's usable sites. Returns `(side0, side1, confidence)` per site.
    fn phase_contig(&self, contig: &str, sites: &[UsableSite]) -> Vec<(u8, u8, f32)> {
        let n = sites.len();
        let k = self.reference.n_haplotypes;
        if n == 0 || k == 0 {
            return Vec::new();
        }
        let allele = |col: usize, hap: usize| self.reference.allele(hap, col);

        // Running IBS match length per reference haplotype: consecutive recent sites where the
        // haplotype's allele is consistent with the observed genotype (homozygous sites only
        // discriminate; heterozygous sites are consistent with every haplotype). Used to pick the
        // per-site switch candidates — the reference haplotypes sharing the longest tract.
        let mut match_len = vec![0u32; k];

        // Beam state: (side0 copied hap, side1 copied hap, ln prob, backpointer into prev beam).
        #[derive(Clone)]
        struct Bs {
            x: u32,
            y: u32,
            lp: f64,
            bp: u32,
        }

        let candidates_at = |match_len: &[u32], m: usize| -> Vec<usize> {
            // Top-m haplotypes by current match length (ties broken by index for determinism).
            let mut idx: Vec<usize> = (0..match_len.len()).collect();
            if idx.len() > m {
                idx.select_nth_unstable_by(m, |&a, &b| match_len[b].cmp(&match_len[a]).then(a.cmp(&b)));
                idx.truncate(m);
            }
            idx.sort_unstable();
            idx
        };

        let update_match = |match_len: &mut [u32], col: usize, g: u8| {
            for (h, ml) in match_len.iter_mut().enumerate() {
                let a = allele(col, h);
                let consistent = match g {
                    0 => a == 0,
                    2 => a == 1,
                    _ => true, // het: any haplotype is consistent
                };
                *ml = if consistent { *ml + 1 } else { 0 };
            }
        };

        // Prune a beam to its top-`width` states by ln prob (unstable partial sort).
        let prune_beam = |beam: &mut Vec<Bs>, width: usize| {
            if beam.len() > width {
                beam.select_nth_unstable_by(width, |a, b| b.lp.total_cmp(&a.lp));
                beam.truncate(width);
            }
        };

        // Initialise the beam at site 0 from the leading candidate set (uniform prior over pairs).
        update_match(&mut match_len, sites[0].ref_col, sites[0].dosage);
        let cand0 = candidates_at(&match_len, self.params.candidates);
        let mut beam: Vec<Bs> = Vec::new();
        for &x in &cand0 {
            let cx = allele(sites[0].ref_col, x);
            for &y in &cand0 {
                let cy = allele(sites[0].ref_col, y);
                beam.push(Bs {
                    x: x as u32,
                    y: y as u32,
                    lp: self.emit_ln(sites[0].dosage, cx, cy),
                    bp: 0,
                });
            }
        }
        prune_beam(&mut beam, self.params.beam);
        let mut trellis: Vec<Vec<Bs>> = Vec::with_capacity(n);
        trellis.push(beam);

        for i in 1..n {
            let d_cm = self
                .map
                .interval_cm(contig, sites[i - 1].position as i32, sites[i].position as i32)
                .unwrap_or(0.0)
                .max(0.0);
            let sw = 1.0 - (-d_cm * self.params.recomb_per_cm).exp();
            let sw = sw.clamp(1e-6, 0.999);
            let stay_ln = (1.0 - sw).ln();
            // A recombination lands on a specific candidate with prob sw/K; we only enumerate the
            // strong candidates but charge the per-target sw/K mass.
            let jump_ln = (sw / k as f64).max(1e-300).ln();

            update_match(&mut match_len, sites[i].ref_col, sites[i].dosage);
            let cand = candidates_at(&match_len, self.params.candidates);

            let prev = &trellis[i - 1];
            let col = sites[i].ref_col;
            let g = sites[i].dosage;

            // Collect candidate successor states, keyed by (x,y), keeping the best incoming lp.
            let mut next: HashMap<(u32, u32), (f64, u32)> = HashMap::new();
            let consider =
                |x: u32, y: u32, base_lp: f64, trans_ln: f64, bp: u32, next: &mut HashMap<(u32, u32), (f64, u32)>| {
                    let lp = base_lp + trans_ln + self.emit_ln(g, allele(col, x as usize), allele(col, y as usize));
                    let e = next.entry((x, y)).or_insert((f64::NEG_INFINITY, 0));
                    if lp > e.0 {
                        *e = (lp, bp);
                    }
                };

            for (bi, s) in prev.iter().enumerate() {
                let bp = bi as u32;
                // Stay on both.
                consider(s.x, s.y, s.lp, stay_ln + stay_ln, bp, &mut next);
                // Side 0 recombines to a candidate; side 1 stays.
                for &x in &cand {
                    let x = x as u32;
                    if x == s.x {
                        continue;
                    }
                    consider(x, s.y, s.lp, jump_ln + stay_ln, bp, &mut next);
                }
                // Side 1 recombines; side 0 stays.
                for &y in &cand {
                    let y = y as u32;
                    if y == s.y {
                        continue;
                    }
                    consider(s.x, y, s.lp, stay_ln + jump_ln, bp, &mut next);
                }
            }

            let mut beam: Vec<Bs> = next.into_iter().map(|((x, y), (lp, bp))| Bs { x, y, lp, bp }).collect();
            prune_beam(&mut beam, self.params.beam);
            trellis.push(beam);
        }

        // Backtrace from the best final state.
        let mut path = vec![(0u32, 0u32); n];
        let last = trellis[n - 1]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.lp.total_cmp(&b.1.lp))
            .map(|(i, _)| i)
            .unwrap_or(0);
        let mut cur = last;
        for i in (0..n).rev() {
            let s = &trellis[i][cur];
            path[i] = (s.x, s.y);
            cur = s.bp as usize;
        }

        path.iter()
            .enumerate()
            .map(|(i, &(x, y))| {
                let col = sites[i].ref_col;
                let c0 = allele(col, x as usize);
                let c1 = allele(col, y as usize);
                match sites[i].dosage {
                    0 => (0, 0, 1.0),
                    2 => (1, 1, 1.0),
                    _ => Self::resolve_het(c0, c1),
                }
            })
            .collect()
    }
}

impl<'a> Phaser for ReferencePhaser<'a> {
    fn phase(&self, genotypes: &[SiteGenotype]) -> PhasedGenotypes {
        // Reference-site lookup: (contig, position) → column index.
        let mut ref_col: HashMap<(&str, i64), usize> = HashMap::with_capacity(self.reference.sites.len());
        for (i, s) in self.reference.sites.iter().enumerate() {
            ref_col.insert((s.contig.as_str(), s.position), i);
        }

        // Group usable sites by contig, position-sorted.
        let mut by_contig: HashMap<String, Vec<UsableSite>> = HashMap::new();
        for g in genotypes {
            if g.dosage < 0 || g.dosage > 2 {
                continue;
            }
            let Some(&col) = ref_col.get(&(g.contig.as_str(), g.position)) else {
                continue;
            };
            by_contig.entry(g.contig.clone()).or_default().push(UsableSite {
                ref_col: col,
                position: g.position,
                dosage: g.dosage as u8,
            });
        }

        let mut out = PhasedGenotypes::default();
        let mut contigs: Vec<String> = by_contig.keys().cloned().collect();
        contigs.sort();
        for contig in contigs {
            let sites = by_contig.get_mut(&contig).unwrap();
            sites.sort_by_key(|s| s.position);
            let phased = self.phase_contig(&contig, sites);
            for (site, (side0, side1, confidence)) in sites.iter().zip(phased) {
                out.sites.push(PhasedSite {
                    contig: contig.clone(),
                    position: site.position,
                    side0,
                    side1,
                    confidence,
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ancestry::HapSite;

    /// Build a two-contig-free single-contig reference from allele rows (each `Vec<u8>` a haplotype).
    fn make_reference(n_sites: usize, rows: &[Vec<u8>]) -> HaplotypeReference {
        let sites: Vec<HapSite> = (0..n_sites)
            .map(|i| HapSite {
                contig: "chr1".to_string(),
                position: (i as i64 + 1) * 10_000,
                reference_allele: 'A',
                alternate_allele: 'G',
            })
            .collect();
        let hap_pop = vec![0u16; rows.len()];
        HaplotypeReference::from_rows("t2t".to_string(), sites, vec!["REF".to_string()], hap_pop, rows)
    }

    fn genotypes_from(contig: &str, hap_a: &[u8], hap_b: &[u8]) -> Vec<SiteGenotype> {
        hap_a
            .iter()
            .zip(hap_b)
            .enumerate()
            .map(|(i, (&a, &b))| SiteGenotype {
                name: String::new(),
                contig: contig.to_string(),
                position: (i as i64 + 1) * 10_000,
                reference_allele: "A".to_string(),
                alternate_allele: "G".to_string(),
                ploidy: 2,
                dosage: (a + b) as i32,
                gq: 60,
                depth: 30,
                ref_depth: 15,
                alt_depth: 15,
                pls: vec![],
                gt: None,
                allele_depths: None,
            })
            .collect()
    }

    fn uniform_map() -> GeneticMap {
        GeneticMap::uniform(1.0, &[("chr1", 250_000_000)])
    }

    #[test]
    fn packing_roundtrips() {
        let rows = vec![vec![0u8, 1, 1, 0], vec![1, 1, 0, 0], vec![0, 0, 0, 1]];
        let r = make_reference(4, &rows);
        for (h, row) in rows.iter().enumerate() {
            for (s, &a) in row.iter().enumerate() {
                assert_eq!(r.allele(h, s), a, "hap {h} site {s}");
            }
        }
        // bincode round-trip.
        let bytes = r.to_bytes().unwrap();
        let back = HaplotypeReference::from_bytes(&bytes).unwrap();
        assert_eq!(back, r);
    }

    /// The sample is a perfect mosaic of two distinct reference haplotypes: side 0 == ref-hap 0,
    /// side 1 == ref-hap 1. The phaser should recover both sides with near-zero switch error.
    #[test]
    fn recovers_known_phase_from_two_donors() {
        // Two clearly distinct donor haplotypes over 60 sites, plus a few decoys.
        let n = 60;
        let hap0: Vec<u8> = (0..n).map(|i| (i % 2) as u8).collect(); // 0,1,0,1,...
        let hap1: Vec<u8> = (0..n).map(|i| ((i / 3) % 2) as u8).collect(); // blocks of 3
        let decoy1: Vec<u8> = (0..n).map(|i| ((i + 1) % 2) as u8).collect();
        let decoy2: Vec<u8> = (0..n).map(|i| ((i / 5) % 2) as u8).collect();
        let rows = vec![hap0.clone(), hap1.clone(), decoy1, decoy2];
        let reference = make_reference(n, &rows);
        let map = uniform_map();
        let phaser = ReferencePhaser::new(&reference, &map, PhaseParams::default());

        let genos = genotypes_from("chr1", &hap0, &hap1);
        let phased = phaser.phase(&genos);
        assert_eq!(phased.sites.len(), n);

        // Count het sites resolved to the correct ordering (allowing a global side-swap, since
        // side labels are arbitrary until anchored).
        let mut agree_direct = 0;
        let mut agree_swapped = 0;
        let mut hets = 0;
        for (i, s) in phased.sites.iter().enumerate() {
            if hap0[i] == hap1[i] {
                // homozygous: both sides must equal the allele
                assert_eq!(s.side0, hap0[i]);
                assert_eq!(s.side1, hap1[i]);
                continue;
            }
            hets += 1;
            if s.side0 == hap0[i] && s.side1 == hap1[i] {
                agree_direct += 1;
            }
            if s.side0 == hap1[i] && s.side1 == hap0[i] {
                agree_swapped += 1;
            }
        }
        assert!(hets > 0);
        let best = agree_direct.max(agree_swapped);
        // Allow a couple of switch errors but demand the phase is overwhelmingly recovered.
        assert!(
            best as f64 >= 0.9 * hets as f64,
            "recovered {best}/{hets} het sites (direct {agree_direct}, swapped {agree_swapped})"
        );
    }

    #[test]
    fn homozygous_sites_are_certain() {
        let n = 20;
        let hap0 = vec![1u8; n];
        let hap1 = vec![1u8; n]; // all-alt homozygous everywhere
        let rows = vec![hap0.clone(), vec![0u8; n]];
        let reference = make_reference(n, &rows);
        let map = uniform_map();
        let phaser = ReferencePhaser::new(&reference, &map, PhaseParams::default());
        let genos = genotypes_from("chr1", &hap0, &hap1);
        let phased = phaser.phase(&genos);
        for s in &phased.sites {
            assert_eq!((s.side0, s.side1), (1, 1));
            assert_eq!(s.confidence, 1.0);
        }
    }
}

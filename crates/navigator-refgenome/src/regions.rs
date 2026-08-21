//! Genome-region metadata (gap §7). For each chromosome it gives the centromere, the telomere
//! caps, the cytoband ideogram, and the chrY pseudoautosomal regions. It parses them from the
//! authoritative UCSC `cytoBand` table, and serves them through a cache of two layers (see
//! [`crate::gateway::ReferenceGateway::genome_regions`]).
//!
//! The data is coordinate context, for QC and for display. It does not feed variant placement.
//!
//! Every interval is **0-based half-open**, as `[start, end)`, which is the UCSC and BED
//! convention. A query position is 1-based, and the code converts it.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::registry::Build;

/// Schema/source version stamped into the cached JSON; bump to invalidate stale caches when the
/// parser or overlay changes.
pub const REGIONS_VERSION: &str = "cytoband-v1";

/// One cytogenetic band (ideogram unit) from the UCSC `cytoBand` table.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cytoband {
    /// Band name within the chromosome, e.g. `p36.33`, `q11.1`.
    pub name: String,
    pub start: i64,
    pub end: i64,
    /// Giemsa stain class: `gneg`, `gpos25/50/75/100`, `acen` (centromere), `gvar`, `stalk`.
    pub stain: String,
}

/// Region metadata for one chromosome.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChromosomeRegions {
    /// Chromosome length (max cytoband end).
    pub length: i64,
    /// Centromere span (merged `acen` bands), if present.
    pub centromere: Option<(i64, i64)>,
    /// Nominal p-arm telomere cap (`[0, 10kb)`).
    pub telomere_p: Option<(i64, i64)>,
    /// Nominal q-arm telomere cap (`[length-10kb, length)`).
    pub telomere_q: Option<(i64, i64)>,
    /// Pseudoautosomal regions (chrY only; from build constants, not cytoBand).
    #[serde(default)]
    pub par: Vec<(i64, i64)>,
    /// The full cytoband ideogram for this chromosome.
    #[serde(default)]
    pub cytobands: Vec<Cytoband>,
}

impl ChromosomeRegions {
    /// The cytoband that holds the 1-based `position`, if there is one.
    pub fn cytoband_at(&self, position: i64) -> Option<&Cytoband> {
        let b = position - 1;
        self.cytobands.iter().find(|c| c.start <= b && b < c.end)
    }
}

/// All region metadata for a build, keyed by the chromosome name as it appears in the source
/// (`chr1`, `chrX`, …).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenomeRegions {
    pub build: String,
    pub version: String,
    pub chromosomes: BTreeMap<String, ChromosomeRegions>,
}

/// What [`GenomeRegions::annotate`] reports for a position.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RegionAnnotation {
    pub in_centromere: bool,
    pub in_telomere: bool,
    pub in_par: bool,
    /// The cytoband name at the position, if known.
    pub cytoband: Option<String>,
}

impl GenomeRegions {
    /// Parse the UCSC `cytoBand` table into the regions of each chromosome. The table columns are
    /// `chrom  start  end  name  gieStain`, and it is 0-based half-open.
    ///
    /// The regions are the length, which is the maximum end, and the centromere, which merges the
    /// `acen` bands. They also hold nominal 10kb telomere caps, and the full band list. PAR comes
    /// from somewhere else, because cytoBand does not carry it.
    pub fn from_cytoband(build: &str, text: &str) -> Self {
        let mut bands: BTreeMap<String, Vec<Cytoband>> = BTreeMap::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut f = line.split('\t');
            let (Some(chrom), Some(s), Some(e), name, stain) = (
                f.next(),
                f.next(),
                f.next(),
                f.next().unwrap_or(""),
                f.next().unwrap_or(""),
            ) else {
                continue;
            };
            if let (Ok(start), Ok(end)) = (s.parse::<i64>(), e.parse::<i64>()) {
                if end > start {
                    bands.entry(chrom.to_string()).or_default().push(Cytoband {
                        name: name.to_string(),
                        start,
                        end,
                        stain: stain.to_string(),
                    });
                }
            }
        }

        let mut chromosomes = BTreeMap::new();
        for (chrom, mut cbs) in bands {
            cbs.sort_by_key(|c| c.start);
            let length = cbs.iter().map(|c| c.end).max().unwrap_or(0);
            let acen: Vec<&Cytoband> = cbs.iter().filter(|c| c.stain == "acen").collect();
            let centromere = if acen.is_empty() {
                None
            } else {
                Some((
                    acen.iter().map(|c| c.start).min().unwrap(),
                    acen.iter().map(|c| c.end).max().unwrap(),
                ))
            };
            let telomere_p = (length > 0).then_some((0, 10_000.min(length)));
            let telomere_q = (length > 10_000).then_some((length - 10_000, length));
            chromosomes.insert(
                chrom,
                ChromosomeRegions {
                    length,
                    centromere,
                    telomere_p,
                    telomere_q,
                    par: Vec::new(),
                    cytobands: cbs,
                },
            );
        }

        let mut regions = GenomeRegions {
            build: build.to_string(),
            version: REGIONS_VERSION.to_string(),
            chromosomes,
        };
        regions.overlay_par(build);
        regions
    }

    /// Overlay the chrY pseudoautosomal regions for the build (PAR is not in cytoBand). Best-known
    /// constants for the builds we resolve; other builds get none.
    fn overlay_par(&mut self, build: &str) {
        let par = crate::registry::canonical_build(build)
            .map(par_regions)
            .unwrap_or_default();
        if par.is_empty() {
            return;
        }
        for key in ["chrY", "Y"] {
            if let Some(c) = self.chromosomes.get_mut(key) {
                c.par = par.clone();
                break;
            }
        }
    }

    /// The regions for `contig`. It accepts a `chr` prefix on one side only (`chr1` against `1`).
    pub fn chromosome(&self, contig: &str) -> Option<&ChromosomeRegions> {
        if let Some(c) = self.chromosomes.get(contig) {
            return Some(c);
        }
        let bare = navigator_domain::contig::bare(contig);
        let alt = if bare.len() < contig.len() {
            bare.to_string()
        } else {
            format!("chr{bare}")
        };
        self.chromosomes.get(&alt)
    }

    /// Annotate a 1-based `position` on `contig` with the region context that covers it.
    pub fn annotate(&self, contig: &str, position: i64) -> RegionAnnotation {
        let Some(c) = self.chromosome(contig) else {
            return RegionAnnotation::default();
        };
        let b = position - 1;
        let within = |iv: &Option<(i64, i64)>| iv.is_some_and(|(s, e)| s <= b && b < e);
        RegionAnnotation {
            in_centromere: within(&c.centromere),
            in_telomere: within(&c.telomere_p) || within(&c.telomere_q),
            in_par: c.par.iter().any(|&(s, e)| s <= b && b < e),
            cytoband: c.cytoband_at(position).map(|cb| cb.name.clone()),
        }
    }
}

/// The chrY pseudoautosomal regions for a build, 0-based half-open. These are well-documented
/// constants for the builds that we resolve. For any other build the list is empty, because a
/// guess would be worse.
fn par_regions(build: Build) -> Vec<(i64, i64)> {
    match build.nuclear() {
        // CHM13v2.0 chrY PAR1 chrY:1–2,458,320 / PAR2 chrY:62,122,809–62,460,029.
        Build::Chm13v2 => vec![(0, 2_458_320), (62_122_808, 62_460_029)],
        // GRCh38 chrY PAR1 chrY:10,001–2,781,479 / PAR2 chrY:56,887,903–57,217,415.
        Build::Grch38 => vec![(10_000, 2_781_479), (56_887_902, 57_217_415)],
        // GRCh37 chrY PAR1 chrY:10,001–2,649,520 / PAR2 chrY:59,034,050–59,373,566.
        Build::Grch37 => vec![(10_000, 2_649_520), (59_034_049, 59_373_566)],
        Build::Chm13v2MaskedRcrs => unreachable!("nuclear() collapses the masked variant"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CYTOBAND: &str = "\
chr1\t0\t2300000\tp36.33\tgneg
chr1\t2300000\t5300000\tp36.32\tgpos25
chr1\t121500000\t125000000\tp11.1\tacen
chr1\t125000000\t128900000\tq11\tacen
chr1\t128900000\t248956422\tq12\tgvar
chrY\t0\t300000\tp11.32\tgneg
chrY\t10300000\t10400000\tp11.1\tacen
chrY\t10400000\t57227415\tq11\tgpos50
";

    #[test]
    fn parses_length_centromere_and_telomeres() {
        let g = GenomeRegions::from_cytoband("GRCh38", CYTOBAND);
        let c1 = g.chromosome("chr1").unwrap();
        assert_eq!(c1.length, 248_956_422);
        assert_eq!(c1.centromere, Some((121_500_000, 128_900_000))); // merged acen p+q
        assert_eq!(c1.telomere_p, Some((0, 10_000)));
        assert_eq!(c1.telomere_q, Some((248_946_422, 248_956_422)));
        assert_eq!(c1.cytobands.len(), 5);
        // chr/no-chr tolerance.
        assert!(g.chromosome("1").is_some());
    }

    #[test]
    fn overlays_chr_y_par() {
        let g = GenomeRegions::from_cytoband("GRCh38", CYTOBAND);
        let y = g.chromosome("chrY").unwrap();
        assert_eq!(y.par, vec![(10_000, 2_781_479), (56_887_902, 57_217_415)]);
    }

    #[test]
    fn annotates_position_context() {
        let g = GenomeRegions::from_cytoband("GRCh38", CYTOBAND);
        // chr1 centromere span (1-based 121,500,001 .. 128,900,000).
        let a = g.annotate("chr1", 124_000_000);
        assert!(a.in_centromere && !a.in_telomere);
        assert_eq!(a.cytoband.as_deref(), Some("p11.1"));
        // p-terminal telomere cap.
        assert!(g.annotate("chr1", 1).in_telomere);
        // chrY PAR1.
        assert!(g.annotate("chrY", 1_000_000).in_par);
        // Unique mid-arm position.
        let u = g.annotate("chr1", 3_000_000);
        assert!(!u.in_centromere && !u.in_telomere && !u.in_par);
        assert_eq!(u.cytoband.as_deref(), Some("p36.32"));
    }

    #[test]
    fn round_trips_through_json() {
        let g = GenomeRegions::from_cytoband("GRCh38", CYTOBAND);
        let json = serde_json::to_string(&g).unwrap();
        let back: GenomeRegions = serde_json::from_str(&json).unwrap();
        assert_eq!(g, back);
        assert_eq!(back.version, REGIONS_VERSION);
    }
}

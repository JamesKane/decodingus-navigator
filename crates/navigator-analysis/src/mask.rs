//! The mask of the callable regions, from a BED file. It holds the variant calls to the regions
//! that the code can trust, such as the callable-Y mask of Poznik and 1KG, in `b38_sites.bed`.
//! Without it, artifacts from a palindrome, from heterochromatin and from a repeat control a
//! de-novo sweep over the whole chrY.
//!
//! BED is 0-based and half-open, as `[start, end)`. The positions in this project are 1-based.
//! The code loads the intervals of the contig that the caller asked for, sorts them, and joins the
//! ones that touch. [`RegionMask::contains`] is then a binary search.

use std::io::BufRead;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// Merged, sorted callable intervals (0-based half-open) for one contig.
#[derive(Debug, Clone)]
pub struct RegionMask {
    intervals: Vec<(i64, i64)>,
}

impl RegionMask {
    /// Load the intervals of `contig` from a BED file. It ignores every other contig. It also
    /// decompresses a BED that came through gzip or BGZF, and the caller sees no difference. It
    /// finds the compression from the content of the file. A large mask that ships with the app
    /// can then stay compressed, and that includes a bgzipped file of more than one block.
    pub fn from_bed(path: &Path, contig: &str) -> Result<Self, AnalysisError> {
        let reader = crate::gzio::open_maybe_gz(path).map_err(|e| AnalysisError::io(path, e))?;
        let mut intervals = Vec::new();
        for line in reader.lines() {
            let line = line.map_err(|e| AnalysisError::io(path, e))?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let mut f = line.split('\t');
            let (Some(c), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else {
                continue;
            };
            if c != contig {
                continue;
            }
            if let (Ok(s), Ok(e)) = (s.parse::<i64>(), e.parse::<i64>()) {
                if e > s {
                    intervals.push((s, e));
                }
            }
        }
        Ok(Self::from_intervals(intervals))
    }

    /// Build from raw `[start, end)` intervals (sorts + coalesces overlaps/adjacencies).
    pub fn from_intervals(mut intervals: Vec<(i64, i64)>) -> Self {
        intervals.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(intervals.len());
        for (s, e) in intervals {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        RegionMask { intervals: merged }
    }

    /// The `[start, end)` intervals, after the code joined the ones that touch. Use this when a
    /// caller must change them, for example to lift a mask to another build. Membership is not the
    /// only question you can ask.
    pub fn intervals(&self) -> &[(i64, i64)] {
        &self.intervals
    }

    /// Total callable bases.
    pub fn covered(&self) -> i64 {
        self.intervals.iter().map(|(s, e)| e - s).sum()
    }

    /// Whether this mask has no intervals.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// Return a new mask with `extra` `[start, end)` intervals added (re-sorted + coalesced).
    pub fn union(&self, extra: &[(i64, i64)]) -> Self {
        let mut all = self.intervals.clone();
        all.extend_from_slice(extra);
        Self::from_intervals(all)
    }

    /// Is the 1-based `position` inside a callable interval?
    pub fn contains(&self, position: i64) -> bool {
        let b = position - 1; // 0-based
        let idx = self.intervals.partition_point(|iv| iv.0 <= b);
        idx > 0 && {
            let (s, e) = self.intervals[idx - 1];
            s <= b && b < e
        }
    }
}

/// The structural class of a chrY region. It gives a Y call *less weight*, and it does not drop
/// that call. The weight follows how well a short read maps there.
///
/// Each class carries a **quality modifier** in `(0, 1]`. That ladder is the port of the Scala
/// `YRegionAnnotator`. Unique sequence, and X-degenerate sequence, carry the full weight: they have
/// no class, and their modifier is 1.0. A zone that is prone to a paralog carries less weight, and
/// so does a repeat zone. The more difficult the zone, the less it carries. A position in no class
/// counts as unique, at a modifier of 1.0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YRegionClass {
    /// A pseudoautosomal region, which recombines with X. The modifier is 0.5.
    Par,
    /// A palindrome or an inverted repeat. Gene conversion and a wrong mapping are both a risk.
    /// The modifier is 0.4.
    Palindrome,
    /// An X-transposed region. It is about 99% the same as X, so contamination is a risk. The
    /// modifier is 0.3.
    Xtr,
    /// An ampliconic block. Its repeat copies are almost the same, so a paralog is a large risk.
    /// The modifier is 0.3.
    Amplicon,
    /// A short-tandem-repeat region. recLOH and stutter are both a risk. The modifier is 0.25.
    Str,
    /// A centromeric region. Almost nothing maps there. The modifier is 0.1.
    Centromere,
    /// Yq12 heterochromatin, and the AZF-DYZ satellite. Nothing maps there. The modifier is 0.1.
    /// The old name of this was `AzfDyz`.
    #[serde(alias = "AzfDyz")]
    Heterochromatin,
}

impl YRegionClass {
    pub fn as_str(self) -> &'static str {
        match self {
            YRegionClass::Par => "par",
            YRegionClass::Palindrome => "palindrome",
            YRegionClass::Xtr => "xtr",
            YRegionClass::Amplicon => "amplicon",
            YRegionClass::Str => "str",
            YRegionClass::Centromere => "centromere",
            YRegionClass::Heterochromatin => "heterochromatin",
        }
    }

    /// The quality modifier, in `(0, 1]`. It says how much a call here counts when the code scores
    /// the concordance of a haplogroup. A lower value means more doubt about a paralog or a
    /// mapping. Unique sequence, which has no class, is 1.0.
    pub fn modifier(self) -> f64 {
        match self {
            YRegionClass::Par => 0.5,
            YRegionClass::Palindrome => 0.4,
            YRegionClass::Xtr | YRegionClass::Amplicon => 0.3,
            YRegionClass::Str => 0.25,
            YRegionClass::Centromere | YRegionClass::Heterochromatin => 0.1,
        }
    }
}

/// PAR1 of chrY in CHM13v2.0, at `chrY:1–2,458,320`, as a 0-based half-open interval. A
/// pseudoautosomal region recombines with X, so a Y-SNP placement never sits here. The code flags
/// such a position for QC, and it does not drop it.
const CHM13_PAR1: (i64, i64) = (0, 2_458_320);
/// CHM13v2.0 chrY PAR2 (`chrY:62,122,809–62,460,029`), 0-based half-open.
const CHM13_PAR2: (i64, i64) = (62_122_808, 62_460_029);
/// The bound of the Yq12 heterochromatin on chrY in CHM13v2.0, at `chrY:26,637,971–62,122,809`,
/// 0-based and half-open. This constant came over from the Scala port, and a check confirmed it.
/// The region is mostly satellite, and nothing maps there.
const CHM13_YQ12_HET: (i64, i64) = (26_637_970, 62_122_809);

/// The chrY PAR and heterochromatin bounds of each build, 0-based and half-open.
///
/// The code does **not lift** these. It takes them from each build natively, for two reasons. Each
/// assembly documents them exactly. And a chain is least reliable in exactly these places: chrX
/// shares a pseudoautosomal region, and Yq12 is satellite. A lift through either one is as likely
/// to be wrong as it is to be absent.
///
/// The palindromes and the amplicons *do* go through a lift. They sit in the euchromatin that only
/// males have, and the chain is reliable there.
///
/// GRCh38 (GCA_000001405.15) has PAR1 at `chrY:10,001–2,781,479`, and PAR2 at
/// `chrY:56,887,903–57,217,415`. The euchromatin that only males have ends at about 26.6 Mb. The
/// heterochromatic arm begins there, and it goes to the end of the assembly.
const GRCH38_PAR1: (i64, i64) = (10_000, 2_781_479);
const GRCH38_PAR2: (i64, i64) = (56_887_902, 57_217_415);
const GRCH38_YQ12_HET: (i64, i64) = (26_600_000, 57_227_415);

/// GRCh37/hg19: PAR1 `chrY:10,001–2,649,520`, PAR2 `chrY:59,034,050–59,363,566`; the
/// heterochromatic arm runs from ~28.8 Mb to the start of PAR2.
const GRCH37_PAR1: (i64, i64) = (10_000, 2_649_520);
const GRCH37_PAR2: (i64, i64) = (59_034_049, 59_363_566);
const GRCH37_YQ12_HET: (i64, i64) = (28_800_000, 59_034_050);

/// One build's chrY geometry: the two pseudoautosomal regions and the heterochromatic arm, all
/// 0-based half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct YLandmarks {
    pub par: [(i64, i64); 2],
    pub heterochromatin: (i64, i64),
}

/// The PAR and heterochromatin constants for a build key (`hs1`/`chm13`, `GRCh38`, `GRCh37`).
/// `None` for a build with no registered chrY geometry.
pub fn y_landmarks(build: &str) -> Option<YLandmarks> {
    let b = build.to_ascii_lowercase();
    if b.contains("chm13") || b.contains("hs1") || b.contains("t2t") {
        Some(YLandmarks {
            par: [CHM13_PAR1, CHM13_PAR2],
            heterochromatin: CHM13_YQ12_HET,
        })
    } else if b.contains("38") {
        Some(YLandmarks {
            par: [GRCH38_PAR1, GRCH38_PAR2],
            heterochromatin: GRCH38_YQ12_HET,
        })
    } else if b.contains("37") || b.contains("hg19") {
        Some(YLandmarks {
            par: [GRCH37_PAR1, GRCH37_PAR2],
            heterochromatin: GRCH37_YQ12_HET,
        })
    } else {
        None
    }
}

/// The curated structural regions of chrY on CHM13, with a quality modifier at each one. They give
/// less weight to a call in a zone that is prone to a paralog, or where nothing maps.
///
/// [`classify`](Self::classify) returns the class that matters most at a position, which is the
/// class with the lowest modifier. [`quality_modifier`](Self::quality_modifier) returns the
/// modifier of that class, and 1.0 for unique sequence.
#[derive(Debug, Clone)]
pub struct YStructuralRegions {
    par: RegionMask,
    palindrome: RegionMask,
    amplicon: RegionMask,
    /// Yq12 / AZF-DYZ satellite + the hardcoded heterochromatin bound.
    heterochromatin: RegionMask,
}

impl YStructuralRegions {
    /// Load from the three chrY BEDs of CHM13: the amplicons, the inverted repeats and
    /// palindromes, and AZF/DYZ. It then adds the CHM13 PAR1, PAR2 and Yq12-heterochromatin
    /// constants that this file holds. The AZF/DYZ BED covers the satellite arrays, and the
    /// constant fills the wider heterochromatic q-arm.
    pub fn from_beds(amplicon: &Path, palindrome: &Path, azf_dyz: &Path) -> Result<Self, AnalysisError> {
        Ok(Self::from_masks(
            RegionMask::from_intervals(vec![CHM13_PAR1, CHM13_PAR2]),
            RegionMask::from_bed(palindrome, "chrY")?,
            RegionMask::from_bed(amplicon, "chrY")?,
            RegionMask::from_bed(azf_dyz, "chrY")?.union(&[CHM13_YQ12_HET]),
        ))
    }

    /// The palindrome mask and the amplicon mask, for a lift to another build. This does not give
    /// out the PAR mask or the heterochromatin mask for that purpose, and that is deliberate. See
    /// [`y_landmarks`].
    pub fn structural_masks(&self) -> (&RegionMask, &RegionMask) {
        (&self.palindrome, &self.amplicon)
    }

    /// The AZF/DYZ and heterochromatin mask, for a lift of the satellite-array intervals.
    pub fn heterochromatin_mask(&self) -> &RegionMask {
        &self.heterochromatin
    }

    /// Build from masks that the caller gives. The BED loader and the unit tests share this seam.
    /// Nobody has found a source for the XTR, STR and centromere masks yet. Those tiers exist in
    /// [`YRegionClass`], ready for the day that their data arrives.
    pub fn from_masks(
        par: RegionMask,
        palindrome: RegionMask,
        amplicon: RegionMask,
        heterochromatin: RegionMask,
    ) -> Self {
        YStructuralRegions {
            par,
            palindrome,
            amplicon,
            heterochromatin,
        }
    }

    /// The structural class that matters most at the 1-based `position`, which is the class with
    /// the lowest modifier. It is `None` when that position sits in unique sequence, which is
    /// sequence that maps reliably, or that is X-degenerate.
    pub fn classify(&self, position: i64) -> Option<YRegionClass> {
        // The checks run from the lowest modifier up, so the class that matters most comes first.
        // Where two classes overlap, the result is then the one that takes the most weight
        // away.
        // An amplicon inside the heterochromatic arm gives Heterochromatin.
        if self.heterochromatin.contains(position) {
            Some(YRegionClass::Heterochromatin)
        } else if self.amplicon.contains(position) {
            Some(YRegionClass::Amplicon)
        } else if self.palindrome.contains(position) {
            Some(YRegionClass::Palindrome)
        } else if self.par.contains(position) {
            Some(YRegionClass::Par)
        } else {
            None
        }
    }

    /// The quality modifier for the 1-based `position` (the most-impactful class's, or 1.0 if the
    /// position is in unique sequence).
    pub fn quality_modifier(&self, position: i64) -> f64 {
        self.classify(position).map_or(1.0, |c| c.modifier())
    }
}

#[cfg(test)]
mod y_landmark_tests {
    use super::*;

    /// The PAR and heterochromatin bounds are a constant of each build, and no lift makes them. A
    /// wrong one then masks a whole assembly wrongly, and nobody sees it. This test holds the
    /// documented coordinates.
    #[test]
    fn each_build_carries_its_own_chry_geometry() {
        let l = y_landmarks("hs1").expect("CHM13");
        assert_eq!(l.par[0], CHM13_PAR1);
        assert_eq!(l.heterochromatin, CHM13_YQ12_HET);

        // GRCh38 PAR1 chrY:10,001-2,781,479 → 0-based half-open.
        let l = y_landmarks("GRCh38").expect("GRCh38");
        assert_eq!(l.par[0], (10_000, 2_781_479));
        assert_eq!(l.par[1], (56_887_902, 57_217_415));
        assert!(
            l.heterochromatin.0 < l.par[1].0,
            "the heterochromatic arm ends before PAR2 begins"
        );

        // GRCh37 PAR1 chrY:10,001-2,649,520.
        let l = y_landmarks("GRCh37").expect("GRCh37");
        assert_eq!(l.par[0], (10_000, 2_649_520));
        assert_eq!(l.par[1], (59_034_049, 59_363_566));
        assert!(l.heterochromatin.0 < l.par[1].0);

        assert!(y_landmarks("hg19").is_some(), "the UCSC name resolves too");
        assert!(
            y_landmarks("GRCm39").is_none(),
            "no geometry invented for an unknown build"
        );
    }

    /// The three builds must not collide: a GRCh37 position inside GRCh38's PAR2 is not in PAR.
    #[test]
    fn the_builds_do_not_share_par_bounds() {
        let g38 = y_landmarks("GRCh38").unwrap().par;
        let g37 = y_landmarks("GRCh37").unwrap().par;
        assert_ne!(g38[1], g37[1], "PAR2 moved between assemblies");
        assert!(g37[1].0 > g38[1].1, "GRCh37 PAR2 starts past the end of GRCh38's");
    }

    #[test]
    fn intervals_round_trip_through_the_accessor() {
        let m = RegionMask::from_intervals(vec![(30, 40), (10, 20), (18, 25)]);
        assert_eq!(m.intervals(), &[(10, 25), (30, 40)], "coalesced and sorted");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_and_tests_membership() {
        // [10,20) and [15,25) coalesce to [10,25); [40,50) separate.
        let m = RegionMask::from_intervals(vec![(40, 50), (10, 20), (15, 25)]);
        assert_eq!(m.covered(), 15 + 10); // [10,25)=15, [40,50)=10

        // 1-based positions: base0 = pos-1.
        assert!(!m.contains(10)); // base0 9 < 10
        assert!(m.contains(11)); // base0 10 in [10,25)
        assert!(m.contains(25)); // base0 24 in [10,25)
        assert!(!m.contains(26)); // base0 25 == end, excluded
        assert!(!m.contains(40)); // base0 39 < 40
        assert!(m.contains(41)); // base0 40 in [40,50)
        assert!(!m.contains(60));
    }

    #[test]
    fn reads_gzipped_bed_transparently() {
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("dun-maskgz-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.bed.gz");
        let mut enc =
            flate2::write::GzEncoder::new(std::fs::File::create(&path).unwrap(), flate2::Compression::default());
        // The loader ignores chrX. There are two chrY intervals, and the code joins one of them
        // to its neighbour.
        enc.write_all(b"chrY\t100\t200\nchrX\t0\t50\nchrY\t150\t260\n").unwrap();
        enc.finish().unwrap();
        let m = RegionMask::from_bed(&path, "chrY").unwrap();
        assert_eq!(m.covered(), 160); // [100,260) after coalescing
        assert!(m.contains(101) && m.contains(260) && !m.contains(261));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn reads_multiblock_bgzf_bed_without_truncation() {
        // A bgzipped BED is a chain of independent gzip members. The GzDecoder path before this
        // one stopped after the first member, and it dropped every later interval where nobody saw
        // it. This test holds that the code reads every member, through MultiGzDecoder in
        // gzio::open_maybe_gz.
        use std::io::Write;
        let dir = std::env::temp_dir().join(format!("dun-maskbgzf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.bed.gz");
        let member = |bytes: &[u8]| {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            enc.write_all(bytes).unwrap();
            enc.finish().unwrap()
        };
        let mut blob = member(b"chrY\t100\t200\n"); // first block
        blob.extend(member(b"chrY\t300\t400\n")); // second block — dropped by a single-member decoder
        std::fs::write(&path, &blob).unwrap();
        let m = RegionMask::from_bed(&path, "chrY").unwrap();
        assert_eq!(m.covered(), 200); // [100,200) + [300,400); would be 100 if truncated
        assert!(m.contains(350));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn region_modifier_ladder() {
        // Unique sequence is full weight; the ladder descends to the most-suspect zones.
        assert_eq!(YRegionClass::Par.modifier(), 0.5);
        assert_eq!(YRegionClass::Palindrome.modifier(), 0.4);
        assert_eq!(YRegionClass::Amplicon.modifier(), 0.3);
        assert_eq!(YRegionClass::Str.modifier(), 0.25);
        assert_eq!(YRegionClass::Heterochromatin.modifier(), 0.1);
        assert_eq!(YRegionClass::Centromere.modifier(), 0.1);
    }

    #[test]
    fn azf_dyz_alias_still_deserializes() {
        // Cached private-Y blobs stored the old "AzfDyz" name → must load as Heterochromatin.
        let c: YRegionClass = serde_json::from_str("\"AzfDyz\"").unwrap();
        assert_eq!(c, YRegionClass::Heterochromatin);
    }

    #[test]
    fn y_structural_classifies_most_impactful_first() {
        // Synthetic regions that do not overlap. The BED loader keys on chrY, and it ignores a
        // chrX line. The constructor that production uses also holds the CHM13 PAR1 and PAR2, and
        // the Yq12 heterochromatin bound.
        let dir = std::env::temp_dir().join(format!("dun-ymask-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let write = |name: &str, body: &str| {
            let p = dir.join(name);
            std::fs::write(&p, body).unwrap();
            p
        };
        let amp = write("amp.bed", "chrY\t3000000\t3000200\tA1\nchrX\t0\t999\n");
        let pal = write("pal.bed", "chrY\t3000300\t3000400\tIR1\n");
        let azf = write("azf.bed", "chrY\t3000500\t3000600\tAZFa\n");
        let r = YStructuralRegions::from_beds(&amp, &pal, &azf).unwrap();

        // Positions are >PAR1 (2,458,320) so they isolate the BED regions.
        assert_eq!(r.classify(3_000_150), Some(YRegionClass::Amplicon));
        assert_eq!(r.classify(3_000_350), Some(YRegionClass::Palindrome));
        assert_eq!(r.classify(3_000_550), Some(YRegionClass::Heterochromatin)); // AZF/DYZ tier
        assert_eq!(r.classify(3_000_700), None); // unique sequence → full weight
        assert_eq!(r.quality_modifier(3_000_700), 1.0);

        // Hardcoded CHM13 constants: PAR1 → Par; the Yq12 arm → Heterochromatin.
        assert_eq!(r.classify(1_000_000), Some(YRegionClass::Par));
        assert_eq!(r.quality_modifier(1_000_000), 0.5);
        assert_eq!(r.classify(30_000_000), Some(YRegionClass::Heterochromatin));
        let _ = std::fs::remove_dir_all(&dir);
    }
}

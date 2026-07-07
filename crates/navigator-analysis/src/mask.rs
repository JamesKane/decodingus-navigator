//! Callable-region mask from a BED file — restricts variant calls to reliable regions
//! (e.g. the Poznik/1KG callable-Y mask, `b38_sites.bed`). Without it, a whole-chrY
//! de-novo sweep is dominated by palindrome/heterochromatin/repeat artifacts.
//!
//! BED is 0-based, half-open `[start, end)`; our positions are 1-based. Intervals for the
//! requested contig are loaded, sorted, and coalesced so [`RegionMask::contains`] is a
//! binary search.

use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// Merged, sorted callable intervals (0-based half-open) for one contig.
#[derive(Debug, Clone)]
pub struct RegionMask {
    intervals: Vec<(i64, i64)>,
}

impl RegionMask {
    /// Load the intervals for `contig` from a BED file (other contigs ignored). A `.gz` path is
    /// transparently gunzipped, so large bundled masks can ship compressed.
    pub fn from_bed(path: &Path, contig: &str) -> Result<Self, AnalysisError> {
        let file = std::fs::File::open(path).map_err(|e| AnalysisError::io(path, e))?;
        let reader: Box<dyn BufRead> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
            Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
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

/// The structural class of a chrY region — for *down-weighting* (not dropping) Y calls by how
/// reliably short reads map there. Each class carries a **quality modifier** in `(0, 1]` (a port of
/// the Scala `YRegionAnnotator` ladder): unique / X-degenerate sequence is full weight (no class,
/// modifier 1.0); paralog-prone and repeat zones get progressively lower weight. A position not in
/// any class is treated as unique (modifier 1.0).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YRegionClass {
    /// Pseudoautosomal region (recombines with X) — modifier 0.5.
    Par,
    /// Palindrome / inverted repeat (gene-conversion + mapping risk) — modifier 0.4.
    Palindrome,
    /// X-transposed region (~99% X-identical, contamination risk) — modifier 0.3.
    Xtr,
    /// Ampliconic block — near-identical repeat copies, high paralog risk — modifier 0.3.
    Amplicon,
    /// Short-tandem-repeat region (recLOH / stutter risk) — modifier 0.25.
    Str,
    /// Centromeric region (nearly unmappable) — modifier 0.1.
    Centromere,
    /// Yq12 heterochromatin / AZF-DYZ satellite (unmappable) — modifier 0.1. (Was `AzfDyz`.)
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

    /// Quality modifier in `(0, 1]`: how much a call here counts in haplogroup concordance scoring
    /// (lower = more paralog-/mapping-suspect). Unique sequence (no class) is 1.0.
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

/// CHM13v2.0 chrY PAR1 (`chrY:1–2,458,320`), as a 0-based half-open interval. The pseudoautosomal
/// regions recombine with X, so Y-SNP placement never lives here — flagged for QC, not dropped.
const CHM13_PAR1: (i64, i64) = (0, 2_458_320);
/// CHM13v2.0 chrY PAR2 (`chrY:62,122,809–62,460,029`), 0-based half-open.
const CHM13_PAR2: (i64, i64) = (62_122_808, 62_460_029);
/// CHM13v2.0 chrY Yq12 heterochromatin bound (`chrY:26,637,971–62,122,809`), 0-based half-open —
/// the validated constant carried over from the Scala port (mostly satellite, unmappable).
const CHM13_YQ12_HET: (i64, i64) = (26_637_970, 62_122_809);

/// Curated CHM13 chrY structural regions with quality modifiers, for down-weighting calls in
/// paralog-prone / unmappable zones. [`classify`](Self::classify) returns the **most-impactful**
/// (lowest-modifier) class containing a position; [`quality_modifier`](Self::quality_modifier)
/// returns its modifier (1.0 for unique sequence).
#[derive(Debug, Clone)]
pub struct YStructuralRegions {
    par: RegionMask,
    palindrome: RegionMask,
    amplicon: RegionMask,
    /// Yq12 / AZF-DYZ satellite + the hardcoded heterochromatin bound.
    heterochromatin: RegionMask,
}

impl YStructuralRegions {
    /// Load from the three CHM13 chrY BEDs (amplicons, inverted-repeats/palindromes, AZF/DYZ),
    /// adding the hardcoded CHM13 PAR1/PAR2 and Yq12-heterochromatin constants (the AZF/DYZ BED
    /// covers the satellite arrays; the constant fills the broader heterochromatic q-arm).
    pub fn from_beds(amplicon: &Path, palindrome: &Path, azf_dyz: &Path) -> Result<Self, AnalysisError> {
        Ok(Self::from_masks(
            RegionMask::from_intervals(vec![CHM13_PAR1, CHM13_PAR2]),
            RegionMask::from_bed(palindrome, "chrY")?,
            RegionMask::from_bed(amplicon, "chrY")?,
            RegionMask::from_bed(azf_dyz, "chrY")?.union(&[CHM13_YQ12_HET]),
        ))
    }

    /// Build from explicit masks (the seam the BED loader + unit tests share). XTR/STR/centromere
    /// masks aren't sourced yet — those tiers exist in [`YRegionClass`] for when their data lands.
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

    /// The most-impactful (lowest-modifier) structural class containing the 1-based `position`, or
    /// `None` if it is in unique (reliably-mappable / X-degenerate) sequence.
    pub fn classify(&self, position: i64) -> Option<YRegionClass> {
        // Checked in ascending-modifier (most-impactful-first) order so overlaps resolve to the
        // strongest down-weight (e.g. an amplicon inside the heterochromatic arm → Heterochromatin).
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
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(&path).unwrap(),
            flate2::Compression::default(),
        );
        // chrX ignored; two chrY intervals, one of them coalescing.
        enc.write_all(b"chrY\t100\t200\nchrX\t0\t50\nchrY\t150\t260\n").unwrap();
        enc.finish().unwrap();
        let m = RegionMask::from_bed(&path, "chrY").unwrap();
        assert_eq!(m.covered(), 160); // [100,260) after coalescing
        assert!(m.contains(101) && m.contains(260) && !m.contains(261));
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
        // Disjoint synthetic regions (the BED loader keys on chrY; a chrX line is ignored). The
        // production constructor also bakes in CHM13 PAR1/PAR2 + the Yq12 heterochromatin bound.
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

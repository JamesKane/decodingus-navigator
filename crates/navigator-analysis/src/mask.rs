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
    /// Load the intervals for `contig` from a BED file (other contigs ignored).
    pub fn from_bed(path: &Path, contig: &str) -> Result<Self, AnalysisError> {
        let file = std::fs::File::open(path).map_err(|e| AnalysisError::io(path, e))?;
        let mut intervals = Vec::new();
        for line in BufReader::new(file).lines() {
            let line = line.map_err(|e| AnalysisError::io(path, e))?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let mut f = line.split('\t');
            let (Some(c), Some(s), Some(e)) = (f.next(), f.next(), f.next()) else { continue };
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

/// The structural class of a chrY region — for *annotating* (not dropping) Y calls that fall
/// in paralog-prone zones, where short-read mapping is unreliable. Ordered most-specific first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YRegionClass {
    /// Ampliconic block — near-identical repeat copies, the highest paralog risk.
    Amplicon,
    /// Palindrome / inverted repeat (the ampliconic mirror structure).
    Palindrome,
    /// AZFa/b/c or DYZ heterochromatic / satellite region.
    AzfDyz,
}

impl YRegionClass {
    pub fn as_str(self) -> &'static str {
        match self {
            YRegionClass::Amplicon => "amplicon",
            YRegionClass::Palindrome => "palindrome",
            YRegionClass::AzfDyz => "azf_dyz",
        }
    }
}

/// Curated CHM13 chrY structural regions (the T2T palindrome / amplicon / AZF-DYZ BEDs),
/// for flagging calls in paralog-prone zones. [`classify`](Self::classify) returns the
/// most-specific class containing a position (amplicon ⊂ palindrome ⊂ the larger AZF/DYZ).
#[derive(Debug, Clone)]
pub struct YStructuralRegions {
    amplicon: RegionMask,
    palindrome: RegionMask,
    azf_dyz: RegionMask,
}

impl YStructuralRegions {
    /// Load from the three CHM13 chrY BEDs (amplicons, inverted-repeats/palindromes, AZF/DYZ).
    pub fn from_beds(amplicon: &Path, palindrome: &Path, azf_dyz: &Path) -> Result<Self, AnalysisError> {
        Ok(YStructuralRegions {
            amplicon: RegionMask::from_bed(amplicon, "chrY")?,
            palindrome: RegionMask::from_bed(palindrome, "chrY")?,
            azf_dyz: RegionMask::from_bed(azf_dyz, "chrY")?,
        })
    }

    /// The most-specific structural class containing the 1-based `position`, or `None` if it is
    /// in unique (reliably-mappable) sequence.
    pub fn classify(&self, position: i64) -> Option<YRegionClass> {
        if self.amplicon.contains(position) {
            Some(YRegionClass::Amplicon)
        } else if self.palindrome.contains(position) {
            Some(YRegionClass::Palindrome)
        } else if self.azf_dyz.contains(position) {
            Some(YRegionClass::AzfDyz)
        } else {
            None
        }
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
    fn y_structural_classifies_most_specific_first() {
        // Nested regions: amplicon ⊂ palindrome ⊂ azf_dyz, all on chrY (a chrX line is ignored).
        let dir = std::env::temp_dir().join(format!("dun-ymask-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let write = |name: &str, body: &str| {
            let p = dir.join(name);
            std::fs::write(&p, body).unwrap();
            p
        };
        let amp = write("amp.bed", "chrY\t100\t200\tA1\nchrX\t0\t999\n");
        let pal = write("pal.bed", "chrY\t100\t300\tIR1\n");
        let azf = write("azf.bed", "chrY\t100\t500\tAZFa\n");
        let r = YStructuralRegions::from_beds(&amp, &pal, &azf).unwrap();

        assert_eq!(r.classify(150), Some(YRegionClass::Amplicon)); // in all → most specific
        assert_eq!(r.classify(250), Some(YRegionClass::Palindrome)); // palindrome+azf only
        assert_eq!(r.classify(400), Some(YRegionClass::AzfDyz)); // azf only
        assert_eq!(r.classify(600), None); // unique sequence
        let _ = std::fs::remove_dir_all(&dir);
    }
}

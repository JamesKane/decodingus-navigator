//! Callable-region mask from a BED file — restricts variant calls to reliable regions
//! (e.g. the Poznik/1KG callable-Y mask, `b38_sites.bed`). Without it, a whole-chrY
//! de-novo sweep is dominated by palindrome/heterochromatin/repeat artifacts.
//!
//! BED is 0-based, half-open `[start, end)`; our positions are 1-based. Intervals for the
//! requested contig are loaded, sorted, and coalesced so [`RegionMask::contains`] is a
//! binary search.

use std::io::{BufRead, BufReader};
use std::path::Path;

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
}

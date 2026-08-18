//! Runs-of-homozygosity domain types.
//!
//! The pattern read is a *classification*, not a rendering: it is computed once by
//! `navigator_analysis::roh` (which re-exports this enum) and then consumed by both the Advanced ROH
//! chart and the Simple-mode brief. It lives here, below the analysis engine, so the brief builder
//! in [`crate::brief`] can switch on the canonical verdict instead of re-deriving its own.

/// Coarse pattern read from the ROH length distribution. Heuristic — for narration, not diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum RohPattern {
    /// Little total ROH — outbred.
    Outbred,
    /// ROH mass dominated by short segments — background relatedness / endogamous population.
    Endogamy,
    /// ROH mass dominated by long segments — recent consanguinity in the pedigree.
    RecentConsanguinity,
    /// Substantial ROH across all classes.
    Mixed,
}

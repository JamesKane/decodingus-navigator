//! Runs-of-homozygosity domain types.
//!
//! The pattern read is a *classification*, and not a way to draw it. `navigator_analysis::roh`
//! computes it one time, and re-exports this enum. Both the Advanced ROH chart and the Simple-mode
//! brief then read it. It lives here, below the analysis engine. The brief builder in
//! [`crate::brief`] can then switch on the canonical verdict, and does not derive its own.

/// Coarse pattern read from the ROH length distribution. It is a heuristic, for narration and not
/// for diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum RohPattern {
    /// Little total ROH: outbred.
    Outbred,
    /// Short segments hold most of the ROH mass: background relatedness, or an endogamous
    /// population.
    Endogamy,
    /// Long segments hold most of the ROH mass: recent consanguinity in the pedigree.
    RecentConsanguinity,
    /// Large ROH across all classes.
    Mixed,
}

//! Ancestry estimation result types — a sample's population-proportion estimate plus the
//! reference population catalog. Pure types; the estimator lives in `navigator-analysis`
//! (which builds these), persistence in `navigator-store`/the app.
//!
//! Phase 1 works at **super-population** granularity (AFR/AMR/EAS/EUR/SAS), the resolution
//! the 1000G-on-CHM13 INFO allele counts give us directly. The fine-grained 26/33-population
//! catalog (and PCA coordinates) is deferred to phase 2 — the `pca_coordinates` field is
//! already carried so the result shape doesn't change when PCA lands.

use serde::{Deserialize, Serialize};

/// A reference (super-)population: its code, display name, and a hex display color.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Population {
    pub code: String,
    pub name: String,
    pub color: String,
}

/// The five 1000 Genomes super-populations, the phase-1 reference set. Order matches the
/// `populations` axis written by the panel builder so callers can zip by index.
pub fn super_populations() -> Vec<Population> {
    [
        ("AFR", "African", "#FF6600"),
        ("AMR", "Admixed American", "#CC0066"),
        ("EAS", "East Asian", "#00CC00"),
        ("EUR", "European", "#0066CC"),
        ("SAS", "South Asian", "#9900CC"),
    ]
    .into_iter()
    .map(|(code, name, color)| Population {
        code: code.to_string(),
        name: name.to_string(),
        color: color.to_string(),
    })
    .collect()
}

/// Display name for a (super-)population code, falling back to the code itself.
pub fn population_name(code: &str) -> String {
    super_populations()
        .into_iter()
        .find(|p| p.code == code)
        .map(|p| p.name)
        .unwrap_or_else(|| code.to_string())
}

/// Hex display color for a (super-)population code, falling back to a neutral grey.
pub fn population_color(code: &str) -> String {
    super_populations()
        .into_iter()
        .find(|p| p.code == code)
        .map(|p| p.color)
        .unwrap_or_else(|| "#888888".to_string())
}

/// Confidence-interval bounds (percent) on a component estimate.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ConfidenceInterval {
    pub lower: f64,
    pub upper: f64,
}

/// One population's share of the estimate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PopulationComponent {
    pub population_code: String,
    pub population_name: String,
    /// 0.0–100.0
    pub percentage: f64,
    pub confidence_interval: ConfidenceInterval,
    /// 1 = highest share.
    pub rank: usize,
}

/// A super-population (continental) summary. With the phase-1 super-population panel this is
/// 1:1 with the components; it stays distinct so the fine-grained phase-2 panel can roll its
/// constituent populations up here without changing the result shape.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SuperPopulationSummary {
    pub super_population: String,
    pub percentage: f64,
    pub populations: Vec<String>,
}

/// A sample's ancestry estimate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AncestryResult {
    /// "aims" | "genome-wide".
    pub panel_type: String,
    /// Total SNPs in the panel.
    pub snps_analyzed: usize,
    /// SNPs with a non-missing genotype.
    pub snps_with_genotype: usize,
    /// SNPs with a no-call.
    pub snps_missing: usize,
    pub components: Vec<PopulationComponent>,
    pub super_population_summary: Vec<SuperPopulationSummary>,
    /// Overall confidence (0–1) from data completeness.
    pub confidence_level: f64,
    pub pipeline_version: String,
    pub reference_version: String,
    /// First N PCA coordinates for visualization (phase 2; `None` for the AF-likelihood path).
    pub pca_coordinates: Option<Vec<f64>>,
}

impl AncestryResult {
    /// The top-ranked (super-)population, if any.
    pub fn primary(&self) -> Option<&PopulationComponent> {
        self.components.first()
    }

    /// Whether more than one super-population exceeds `threshold` percent (admixed signal).
    pub fn is_admixed(&self, threshold: f64) -> bool {
        self.super_population_summary
            .iter()
            .filter(|s| s.percentage >= threshold)
            .count()
            > 1
    }
}

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
        // Added with the SGDP diversity panel (continents 1000G doesn't cover).
        ("MEA", "Middle Eastern", "#996633"),
        ("CAS", "Central Asian & Siberian", "#66CCCC"),
        ("OCE", "Oceanian", "#009999"),
    ]
    .into_iter()
    .map(|(code, name, color)| Population {
        code: code.to_string(),
        name: name.to_string(),
        color: color.to_string(),
    })
    .collect()
}

/// The 26 fine-grained 1000 Genomes populations: `(code, name, super-population)`. Used for
/// fine-grained ancestry; each rolls up to one of [`super_populations`] for the summary.
const FINE_POPULATIONS: [(&str, &str, &str); 29] = [
    // African
    ("YRI", "Yoruba (Nigeria)", "AFR"),
    ("LWK", "Luhya (Kenya)", "AFR"),
    ("GWD", "Gambian", "AFR"),
    ("MSL", "Mende (Sierra Leone)", "AFR"),
    ("ESN", "Esan (Nigeria)", "AFR"),
    ("ASW", "African-American (SW US)", "AFR"),
    ("ACB", "African-Caribbean (Barbados)", "AFR"),
    // Admixed American
    ("MXL", "Mexican (LA)", "AMR"),
    ("PUR", "Puerto Rican", "AMR"),
    ("CLM", "Colombian", "AMR"),
    ("PEL", "Peruvian", "AMR"),
    // East Asian
    ("CHB", "Han Chinese (Beijing)", "EAS"),
    ("JPT", "Japanese (Tokyo)", "EAS"),
    ("CHS", "Southern Han Chinese", "EAS"),
    ("CDX", "Dai Chinese", "EAS"),
    ("KHV", "Kinh (Vietnam)", "EAS"),
    // European
    ("CEU", "NW European (Utah)", "EUR"),
    ("TSI", "Tuscan (Italy)", "EUR"),
    ("FIN", "Finnish", "EUR"),
    ("GBR", "British", "EUR"),
    ("IBS", "Iberian (Spain)", "EUR"),
    // South Asian
    ("GIH", "Gujarati", "SAS"),
    ("PJL", "Punjabi", "SAS"),
    ("BEB", "Bengali", "SAS"),
    ("STU", "Sri Lankan Tamil", "SAS"),
    ("ITU", "Indian Telugu", "SAS"),
    // SGDP-backed continents (each is a single reference group = its own super-population).
    ("MEA", "Middle Eastern", "MEA"),
    ("CAS", "Central Asian & Siberian", "CAS"),
    ("OCE", "Oceanian", "OCE"),
];

/// The super-population a (fine or super) population code belongs to, or `None` if unknown.
pub fn population_super(code: &str) -> Option<&'static str> {
    if super_populations().iter().any(|p| p.code == code) {
        // A super-population code maps to itself.
        return SUPER_CODES.iter().copied().find(|&c| c == code);
    }
    FINE_POPULATIONS.iter().find(|(c, _, _)| *c == code).map(|(_, _, sp)| *sp)
}

const SUPER_CODES: [&str; 8] = ["AFR", "AMR", "EAS", "SAS", "EUR", "MEA", "CAS", "OCE"];

/// Display name for a fine or super population code, falling back to the code itself.
pub fn population_name(code: &str) -> String {
    if let Some((_, name, _)) = FINE_POPULATIONS.iter().find(|(c, _, _)| *c == code) {
        return name.to_string();
    }
    super_populations()
        .into_iter()
        .find(|p| p.code == code)
        .map(|p| p.name)
        .unwrap_or_else(|| code.to_string())
}

/// Hex display color for a population code. Fine populations inherit their super-population's
/// color (so PCA clusters read by continent); falls back to neutral grey.
pub fn population_color(code: &str) -> String {
    let key = population_super(code).unwrap_or(code);
    super_populations()
        .into_iter()
        .find(|p| p.code == key)
        .map(|p| p.color)
        .unwrap_or_else(|| "#888888".to_string())
}

/// Approximate `(longitude, latitude)` of a population's homeland, for the geographic map.
/// Representative points (degrees); fine populations are placed in-country, super-only groups
/// (MEA/CAS/OCE) at a regional centroid.
pub fn population_lonlat(code: &str) -> Option<(f32, f32)> {
    let p = match code {
        // African
        "YRI" => (8.0, 8.0),
        "ESN" => (6.0, 6.5),
        "GWD" => (-15.5, 13.5),
        "MSL" => (-12.0, 8.5),
        "LWK" => (37.0, 0.3),
        "ASW" => (-90.0, 33.0),
        "ACB" => (-59.5, 13.2),
        "AFR" => (20.0, 2.0),
        // Admixed American
        "MXL" => (-102.0, 23.0),
        "PUR" => (-66.5, 18.2),
        "CLM" => (-74.3, 4.6),
        "PEL" => (-77.0, -12.0),
        "AMR" => (-85.0, 5.0),
        // East Asian
        "CHB" => (116.4, 39.9),
        "JPT" => (139.7, 35.7),
        "CHS" => (113.3, 28.2),
        "CDX" => (101.0, 22.0),
        "KHV" => (106.7, 16.5),
        "EAS" => (112.0, 35.0),
        // European
        "CEU" => (5.0, 52.0),
        "GBR" => (-2.0, 54.0),
        "FIN" => (26.0, 64.0),
        "IBS" => (-4.0, 40.0),
        "TSI" => (11.0, 43.5),
        "EUR" => (10.0, 50.0),
        // South Asian
        "GIH" => (72.0, 23.0),
        "PJL" => (74.3, 31.5),
        "BEB" => (90.4, 23.8),
        "STU" => (81.0, 7.0),
        "ITU" => (79.0, 16.5),
        "SAS" => (78.0, 22.0),
        // SGDP-added regions
        "MEA" => (45.0, 31.0),
        "CAS" => (95.0, 62.0),
        "OCE" => (143.0, -7.0),
        _ => return None,
    };
    Some(p)
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

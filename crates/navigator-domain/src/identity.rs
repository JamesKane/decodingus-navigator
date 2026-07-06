//! Vendor-neutral Subject identity + FTDNA-specific member/MDKA types (FTDNA project-import
//! design §4). Pure types, no IO.
//!
//! **Privacy:** [`ExternalId`], [`FtdnaMember`], and [`Mdka`] are **PII / never-federated** — they
//! must not be derived into a public PDS `fed` record nor put in an AppView-bound payload. They may
//! only ever enter the encrypted Edge-to-Edge tier. Keep distinct from our own computed haplogroup
//! calls (those live in `RunHaplogroupCall`).

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// A Subject's membership in a project (the M:N join, design §4.1).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProjectMembership {
    pub biosample_guid: SampleGuid,
    pub project_id: i64,
    /// Optional subgroup/branch label within the project.
    pub role: Option<String>,
    /// ISO-8601.
    pub added_at: String,
}

/// A vendor identifier for a Subject. `(source, external_id)` is the global cross-project dedup
/// anchor the matching engine keys on (design §4.2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExternalId {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// `FTDNA` | `YSEQ` | `NEBULA` | `WGS` | `MANUAL` | … — see [`IdSource`] for the well-known set.
    pub source: String,
    /// Kit number / vendor id.
    pub external_id: String,
}

/// Well-known [`ExternalId::source`] values. Stored as plain strings (open set — new vendors are
/// just a new value), but the common ones get constants to avoid typos at call sites.
pub struct IdSource;
impl IdSource {
    // ── vendor kits (background-only on the AppView — never surfaced publicly) ──
    pub const FTDNA: &'static str = "FTDNA";
    pub const YSEQ: &'static str = "YSEQ";
    pub const NEBULA: &'static str = "NEBULA";
    pub const DANTE: &'static str = "DANTE";
    pub const FGC: &'static str = "FGC";
    pub const WGS: &'static str = "WGS";
    pub const MANUAL: &'static str = "MANUAL";
    /// The Big Y variant/BAM package's internal sample UUID (links BAM ↔ variants; design §5).
    pub const FTDNA_BIGY_UUID: &'static str = "FTDNA_BIGY_UUID";

    // ── public / open-consent catalog ids (the AppView surfaces these) ──
    // These namespace tokens MUST match the AppView's `is_public` set exactly — it derives
    // displayability from the namespace, so a typo silently demotes a public id to background-only.
    pub const PGP: &'static str = "PGP";
    pub const IGSR: &'static str = "IGSR";
    pub const THOUSAND_GENOMES: &'static str = "1000G";
    pub const ENA: &'static str = "ENA";
    pub const SRA: &'static str = "SRA";
    pub const BIOSAMPLE: &'static str = "BIOSAMPLE";
    pub const HGDP: &'static str = "HGDP";
    pub const SGDP: &'static str = "SGDP";

    /// Whether a namespace is a public/open-consent catalog id (surfaced by the AppView) rather than
    /// a vendor kit (kept off every public surface). Mirrors the AppView's `is_public` policy so the
    /// two ends agree; an unrecognized namespace is treated as private (the safe default).
    pub fn is_public(source: &str) -> bool {
        matches!(
            source,
            Self::PGP | Self::IGSR | Self::THOUSAND_GENOMES | Self::ENA | Self::SRA | Self::BIOSAMPLE | Self::HGDP | Self::SGDP
        )
    }
}

/// Public/open-consent catalog identifiers derivable **purely from a sample's local provenance** —
/// used to seed the AppView-visible `external_ids` for bulk-imported public datasets so they match
/// their existing catalog rows. Deterministic pattern match only (no network/manifest lookup):
///
/// - a 1000 Genomes / IGSR sample name (`HG#####` / `NA#####`) → `(IGSR, name)`;
/// - an HGDP catalog id (`HGDP#####`) → `(HGDP, name)`;
/// - a genuine INSDC **sample** accession in `sample_accession` (`SAM*` → BIOSAMPLE, `ERS…` → ENA,
///   `SRS…` → SRA).
///
/// A dataset-specific friendly name (the common case in ancient-DNA / population sets, where the
/// accession is just a copy of the label) yields nothing — we never guess a namespace, because a
/// wrong token silently fails the AppView's `(namespace, value)` dedup. GIAB `HG00x` (< 5 digits)
/// is intentionally excluded to avoid colliding with build names.
pub fn catalog_ids_from_provenance(donor_identifier: &str, sample_accession: Option<&str>) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let d = donor_identifier.trim();
    if is_igsr_name(d) {
        out.push((IdSource::IGSR.to_string(), d.to_string()));
    } else if is_hgdp_name(d) {
        out.push((IdSource::HGDP.to_string(), d.to_string()));
    }
    if let Some(acc) = sample_accession.map(str::trim).filter(|s| !s.is_empty()) {
        if let Some(ns) = insdc_sample_namespace(acc) {
            out.push((ns.to_string(), acc.to_string()));
        }
    }
    out
}

/// `HG#####` / `NA#####` — a 1000 Genomes / IGSR sample name (≥ 5 digits after the prefix).
fn is_igsr_name(s: &str) -> bool {
    let rest = s.strip_prefix("HG").or_else(|| s.strip_prefix("NA"));
    matches!(rest, Some(r) if r.len() >= 5 && r.bytes().all(|b| b.is_ascii_digit()))
}

/// `HGDP#####` (optionally `HGDP_#####`) — an HGDP catalog id.
fn is_hgdp_name(s: &str) -> bool {
    let rest = s.strip_prefix("HGDP").map(|r| r.strip_prefix('_').unwrap_or(r));
    matches!(rest, Some(r) if !r.is_empty() && r.bytes().all(|b| b.is_ascii_digit()))
}

/// The INSDC **sample**-accession namespace for `acc`, if it's a real one (not a friendly name):
/// `SAM*` → BIOSAMPLE, `ERS…` → ENA, `SRS…` → SRA. `None` for anything else (a plain friendly name).
/// Used both by [`catalog_ids_from_provenance`] and by the API-driven accession backfill.
pub fn insdc_sample_namespace(acc: &str) -> Option<&'static str> {
    let u = acc.to_ascii_uppercase();
    let digits_after = |p: &str| u.strip_prefix(p).is_some_and(|r| !r.is_empty() && r.bytes().all(|b| b.is_ascii_digit()));
    if u.starts_with("SAMN") || u.starts_with("SAMEA") || u.starts_with("SAMD") {
        Some(IdSource::BIOSAMPLE)
    } else if digits_after("ERS") {
        Some(IdSource::ENA)
    } else if digits_after("SRS") {
        Some(IdSource::SRA)
    } else {
        None
    }
}

/// FTDNA-reported member labels only (the batch-file metadata we don't otherwise model). Computed
/// haplogroups stay in the haplogroup-call store — different provenance (design §4.2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FtdnaMember {
    pub biosample_guid: SampleGuid,
    pub member_name: Option<String>,
    pub y_haplogroup_ftdna: Option<String>,
    pub mt_haplogroup_ftdna: Option<String>,
    /// `predicted` | `confirmed`.
    pub haplo_status: Option<String>,
    /// `Advanced` | `Limited` | `None` — the pose-as gate, which also determines the reachable Big Y
    /// data tier (design §3.5).
    pub access_granted: Option<String>,
    /// `Publicly Share DNA Results` consent flag — gates whether this Subject may federate.
    pub publicly_shares: Option<bool>,
}

/// Lineage a [`Mdka`] (or haplogroup) pertains to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Lineage {
    /// Paternal (Y).
    Y,
    /// Maternal (mtDNA).
    Mt,
    /// Autosomal / general.
    Auto,
}

impl Lineage {
    pub fn as_str(self) -> &'static str {
        match self {
            Lineage::Y => "Y",
            Lineage::Mt => "Mt",
            Lineage::Auto => "Auto",
        }
    }

    pub fn parse(s: &str) -> Option<Lineage> {
        match s {
            "Y" => Some(Lineage::Y),
            "Mt" => Some(Lineage::Mt),
            "Auto" => Some(Lineage::Auto),
            _ => None,
        }
    }
}

/// Most Distant Known Ancestor on a lineage (design §4.3). One per Subject per lineage.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mdka {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// `Y` | `Mt` | `Auto` (see [`Lineage`]).
    pub lineage: String,
    pub ancestor_name: Option<String>,
    pub birth_year: Option<i32>,
    pub death_year: Option<i32>,
    pub origin_place: Option<String>,
    pub origin_country: Option<String>,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
    pub source: Option<String>,
    pub notes: Option<String>,
    /// ISO-8601.
    pub updated_at: String,
}

/// Insert/update payload for an MDKA row (the store stamps `id`/`updated_at`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NewMdka {
    pub lineage: String,
    pub ancestor_name: Option<String>,
    pub birth_year: Option<i32>,
    pub death_year: Option<i32>,
    pub origin_place: Option<String>,
    pub origin_country: Option<String>,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
    pub source: Option<String>,
    pub notes: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_ids_from_clean_names_only() {
        // 1000G / IGSR sample names → IGSR.
        assert_eq!(
            catalog_ids_from_provenance("HG00096", None),
            vec![("IGSR".to_string(), "HG00096".to_string())]
        );
        assert_eq!(
            catalog_ids_from_provenance("NA12878", None),
            vec![("IGSR".to_string(), "NA12878".to_string())]
        );
        // HGDP catalog id (with or without underscore) → HGDP, no collision with the HG prefix.
        assert_eq!(
            catalog_ids_from_provenance("HGDP00521", None),
            vec![("HGDP".to_string(), "HGDP00521".to_string())]
        );
        assert_eq!(catalog_ids_from_provenance("HGDP_00521", None)[0].0, "HGDP");
        // Dataset friendly names (the bulk-set common case) → nothing; we never guess.
        assert!(catalog_ids_from_provenance("Ale22", Some("Ale22")).is_empty());
        assert!(catalog_ids_from_provenance("BulgarianB4", Some("BulgarianB4")).is_empty());
        // GIAB HG002 (< 5 digits) is deliberately excluded to avoid build-name collisions.
        assert!(catalog_ids_from_provenance("HG002", None).is_empty());
    }

    #[test]
    fn catalog_ids_from_real_insdc_accessions() {
        assert_eq!(
            catalog_ids_from_provenance("Ale22", Some("SAMEA3302884")),
            vec![("BIOSAMPLE".to_string(), "SAMEA3302884".to_string())]
        );
        assert_eq!(catalog_ids_from_provenance("x", Some("ERS1234567"))[0].0, "ENA");
        assert_eq!(catalog_ids_from_provenance("x", Some("SRS999999"))[0].0, "SRA");
        // A clean name + a real accession yields both.
        let both = catalog_ids_from_provenance("HG00096", Some("SAMN12345678"));
        assert_eq!(both.len(), 2);
        assert!(both.contains(&("IGSR".to_string(), "HG00096".to_string())));
        assert!(both.contains(&("BIOSAMPLE".to_string(), "SAMN12345678".to_string())));
    }

    #[test]
    fn is_public_namespaces() {
        assert!(IdSource::is_public("IGSR"));
        assert!(IdSource::is_public("PGP"));
        assert!(IdSource::is_public("BIOSAMPLE"));
        assert!(!IdSource::is_public("FTDNA"));
        assert!(!IdSource::is_public("YSEQ"));
        assert!(!IdSource::is_public("WGS")); // vendor/background — the WGS229-style reconciliation case
        assert!(!IdSource::is_public("SOMETHING_NEW")); // unknown → private (safe default)
    }
}

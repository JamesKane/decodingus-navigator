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
    pub const FTDNA: &'static str = "FTDNA";
    pub const YSEQ: &'static str = "YSEQ";
    pub const NEBULA: &'static str = "NEBULA";
    pub const WGS: &'static str = "WGS";
    pub const MANUAL: &'static str = "MANUAL";
    /// The Big Y variant/BAM package's internal sample UUID (links BAM ↔ variants; design §5).
    pub const FTDNA_BIGY_UUID: &'static str = "FTDNA_BIGY_UUID";
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

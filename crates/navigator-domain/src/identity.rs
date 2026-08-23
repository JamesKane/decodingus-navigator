//! Vendor-neutral Subject identity, and the FTDNA-specific member and MDKA types (FTDNA
//! project-import design §4). Pure types, no IO.
//!
//! **Privacy:** [`ExternalId`], [`FtdnaMember`] and [`Mdka`] are **PII, and never federated**. No
//! code may derive them into a public PDS `fed` record, or put them in a payload bound for the
//! AppView. They may enter the encrypted Edge-to-Edge tier, and nothing else. Keep them separate
//! from the haplogroup calls we compute, which live in `RunHaplogroupCall`.

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
    /// `FTDNA` | `YSEQ` | `NEBULA` | `WGS` | `MANUAL` | … See [`IdSource`] for the well-known set.
    pub source: String,
    /// Kit number / vendor id.
    pub external_id: String,
}

/// Well-known [`ExternalId::source`] values. The store holds them as plain strings, because the
/// set is open and a new vendor is only a new value. The common ones have constants, to stop a typo
/// at a call site.
pub struct IdSource;
impl IdSource {
    // ── vendor kits (background-only on the AppView, and never public) ──
    pub const FTDNA: &'static str = "FTDNA";
    pub const YSEQ: &'static str = "YSEQ";
    pub const NEBULA: &'static str = "NEBULA";
    pub const DANTE: &'static str = "DANTE";
    pub const FGC: &'static str = "FGC";
    pub const WGS: &'static str = "WGS";
    pub const MANUAL: &'static str = "MANUAL";
    /// The Big Y variant/BAM package's internal sample UUID (links BAM ↔ variants; design §5).
    pub const FTDNA_BIGY_UUID: &'static str = "FTDNA_BIGY_UUID";

    // ── public / open-consent catalog ids (the AppView shows these) ──
    // These namespace tokens MUST match the `is_public` set of the AppView exactly. The AppView
    // decides what to display from the namespace, so a typo demotes a public id to
    // background-only, and gives no message.
    pub const PGP: &'static str = "PGP";
    pub const IGSR: &'static str = "IGSR";
    pub const THOUSAND_GENOMES: &'static str = "1000G";
    pub const ENA: &'static str = "ENA";
    pub const SRA: &'static str = "SRA";
    pub const BIOSAMPLE: &'static str = "BIOSAMPLE";
    pub const HGDP: &'static str = "HGDP";
    pub const SGDP: &'static str = "SGDP";

    /// True when a namespace is a public open-consent catalog id, which the AppView shows, and
    /// not a vendor kit, which stays off every public surface. This mirrors the `is_public` policy
    /// of the AppView, so that the two ends agree. A namespace it does not recognize is private,
    /// which is the safe default.
    pub fn is_public(source: &str) -> bool {
        matches!(
            source,
            Self::PGP
                | Self::IGSR
                | Self::THOUSAND_GENOMES
                | Self::ENA
                | Self::SRA
                | Self::BIOSAMPLE
                | Self::HGDP
                | Self::SGDP
        )
    }
}

/// Public open-consent catalog identifiers that come **only from the local provenance of a
/// sample**. They seed the `external_ids` the AppView can see, for a public dataset that a bulk
/// import brought in. The dataset then matches the catalog rows that already exist. This is a
/// deterministic pattern match only, with no network or manifest lookup:
///
/// - a 1000 Genomes or IGSR sample name (`HG#####` / `NA#####`) → `(IGSR, name)`;
/// - an HGDP catalog id (`HGDP#####`) → `(HGDP, name)`;
/// - a genuine INSDC **sample** accession in `sample_accession` (`SAM*` → BIOSAMPLE, `ERS…` → ENA,
///   `SRS…` → SRA).
///
/// A friendly name that belongs to one dataset gives nothing. That is the common case in
/// ancient-DNA and population sets, where the accession is only a copy of the label. We never guess
/// a namespace, because a wrong token fails the `(namespace, value)` dedup of the AppView, and
/// gives no message. GIAB `HG00x` (< 5 digits) stays out on purpose, so that it does not collide
/// with a build name.
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

/// `HG#####` / `NA#####`: a 1000 Genomes or IGSR sample name (≥ 5 digits after the prefix).
fn is_igsr_name(s: &str) -> bool {
    let rest = s.strip_prefix("HG").or_else(|| s.strip_prefix("NA"));
    matches!(rest, Some(r) if r.len() >= 5 && r.bytes().all(|b| b.is_ascii_digit()))
}

/// `HGDP#####` (also `HGDP_#####`): an HGDP catalog id.
fn is_hgdp_name(s: &str) -> bool {
    let rest = s.strip_prefix("HGDP").map(|r| r.strip_prefix('_').unwrap_or(r));
    matches!(rest, Some(r) if !r.is_empty() && r.bytes().all(|b| b.is_ascii_digit()))
}

/// The INSDC **sample**-accession namespace for `acc`, if it is a real one (not a friendly name):
/// `SAM*` → BIOSAMPLE, `ERS…` → ENA, `SRS…` → SRA. `None` for anything else (a plain friendly name).
/// Used both by [`catalog_ids_from_provenance`] and by the API-driven accession backfill.
pub fn insdc_sample_namespace(acc: &str) -> Option<&'static str> {
    let u = acc.to_ascii_uppercase();
    let digits_after = |p: &str| {
        u.strip_prefix(p)
            .is_some_and(|r| !r.is_empty() && r.bytes().all(|b| b.is_ascii_digit()))
    };
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

/// FTDNA-reported member labels only: the batch-file metadata we do not model in another place. A
/// computed haplogroup stays in the haplogroup-call store, because its provenance is different
/// (design §4.2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FtdnaMember {
    pub biosample_guid: SampleGuid,
    pub member_name: Option<String>,
    pub y_haplogroup_ftdna: Option<String>,
    pub mt_haplogroup_ftdna: Option<String>,
    /// `predicted` | `confirmed`.
    pub haplo_status: Option<String>,
    /// `Advanced` | `Limited` | `None`: the pose-as gate. It also sets which Big Y data tier the
    /// code can reach (design §3.5).
    pub access_granted: Option<String>,
    /// `Publicly Share DNA Results` consent flag. It gates whether this Subject may federate.
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

/// Most Distant Known Ancestor on a lineage (design §4.3). One for each Subject and each lineage.
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

/// Name particles that belong to the surname, and not to a given name. These are the tokens a
/// surname may correctly start with. Without them, [`surname_of`] would refuse `van der Berg` and
/// `de la Cruz`, which are surnames. It still refuses `Thomas Michael Kane`, which is not one.
const NAME_PARTICLES: &[&str] = &[
    "van", "von", "der", "den", "de", "del", "della", "di", "da", "dos", "du", "la", "le", "les", "mac", "mc", "st",
    "st.", "saint", "ter", "ten", "af", "av", "al", "bin", "ibn", "ap", "ó", "ni", "nic", "mag", "fitz", "o", "o'",
];

/// A token that starts the biographical tail, and does not continue the name: a year, a date, or
/// the word before one. Everything from here on is annotation.
fn is_annotation(token: &str) -> bool {
    let t = token.trim_matches(|c: char| !c.is_alphanumeric()).to_lowercase();
    if matches!(
        t.as_str(),
        "b" | "d" | "born" | "died" | "abt" | "about" | "circa" | "c" | "of" | "in"
    ) {
        return true;
    }
    // Any run of four digits that reads as a year (1000-2099). This covers `1770`, `~1770`,
    // `1919-1996` and `11/25/1843`.
    token.as_bytes().windows(4).any(|w| {
        w.iter().all(u8::is_ascii_digit) && {
            let y: i32 = std::str::from_utf8(w).unwrap_or("0").parse().unwrap_or(0);
            (1000..=2099).contains(&y)
        }
    })
}

/// Decode the HTML entities the FTDNA CSV importer leaves in place (`L&#225;ire`, `Died&#160;26`).
///
/// The root cause is the importer, and not this function. But this is the last point before a name
/// goes out, and to send `mac L&#225;ire` as a surname is worse than to decode it here. In the
/// reference corpus, 61 of the 6,218 names carry one.
fn decode_entities(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        let tail = &rest[i..];
        let Some(end) = tail.find(';').filter(|&e| e <= 8) else {
            out.push('&');
            rest = &tail[1..];
            continue;
        };
        let body = &tail[1..end];
        let ch = match body {
            "amp" => Some('&'),
            "lt" => Some('<'),
            "gt" => Some('>'),
            "quot" => Some('"'),
            "apos" => Some('\''),
            "nbsp" => Some(' '),
            _ => body
                .strip_prefix('#')
                .and_then(|n| match n.strip_prefix(['x', 'X']) {
                    Some(hex) => u32::from_str_radix(hex, 16).ok(),
                    None => n.parse::<u32>().ok(),
                })
                .and_then(char::from_u32),
        };
        match ch {
            Some(c) => {
                out.push(c);
                rest = &tail[end + 1..];
            }
            None => {
                out.push('&');
                rest = &tail[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// Generational suffixes. This code drops them before it takes the surname.
const NAME_SUFFIXES: &[&str] = &["jr", "jr.", "sr", "sr.", "i", "ii", "iii", "iv", "v", "esq", "esq."];

/// Reduce a most-distant-known-ancestor's full name to a **surname**.
///
/// This is a privacy gate, and not a matter of style. The surname, origin and dates of an MDKA are
/// genealogical context that may go out (`proposals/ancestral-origin-icicle.md` §2 in the AppView
/// repo). A given name may not, and a given name is what turns a published record into a named
/// individual. So the split happens here, at the edge, before anything becomes bytes. The full name
/// never leaves the workspace.
///
/// The rule is conservative. It takes the **last** token, plus any particle directly before it, and
/// it returns `None` when there is nothing usable. A wrong split leaks a forename, so the AppView
/// checks again, on its own, what arrives. This is the first of two gates, and not the only one.
///
/// ```text
/// "Thomas Michael Kane"      → "Kane"
/// "Kane"                     → "Kane"
/// "Pieter van der Berg"      → "van der Berg"
/// "Juan de la Cruz"          → "de la Cruz"
/// "Patrick O'Brien Jr."      → "O'Brien"
/// "Kane, Thomas"             → "Kane"          (comma-first form)
/// ```
pub fn surname_of(full_name: &str) -> Option<String> {
    let decoded = decode_entities(full_name);
    // `Surname, Given`. Genealogy files are full of it, and the plain last-token rule would take
    // the given name.
    let head = match decoded.split_once(',') {
        Some((last, _)) if !last.trim().is_empty() => last.to_string(),
        _ => decoded,
    };
    let mut tokens: Vec<&str> = head
        .split_whitespace()
        .filter(|t| !t.trim_matches(|c: char| !c.is_alphanumeric()).is_empty())
        .collect();
    // Cut at the first biographical annotation. A third of the reference corpus adds dates or a
    // birthplace to the name: `William Macaulay ~1770 of Balnicol`, `Michael OConnell b1854 d1928
    // St Louis`. A plain last-token rule takes `Balnicol` and `Louis` as surnames.
    if let Some(cut) = tokens.iter().position(|t| is_annotation(t)) {
        tokens.truncate(cut);
    }
    // Drop a generational suffix at the end (`Jr.`, `III`).
    while tokens.last().is_some_and(|t| {
        NAME_SUFFIXES.contains(&t.to_lowercase().trim_end_matches('.').to_string().as_str())
            || NAME_SUFFIXES.contains(&t.to_lowercase().as_str())
    }) {
        tokens.pop();
    }
    let last = tokens.pop()?;
    // Walk back over particles so multi-token surnames survive intact.
    let mut parts = vec![last];
    while let Some(prev) = tokens.last() {
        if NAME_PARTICLES.contains(&prev.to_lowercase().as_str()) {
            parts.push(tokens.pop().unwrap());
        } else {
            break;
        }
    }
    parts.reverse();
    let surname = parts.join(" ");
    (!surname.is_empty()).then_some(surname)
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
        // GIAB HG002 (< 5 digits) stays out on purpose, so that it does not collide with a build
        // name.
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

    #[test]
    fn surname_drops_the_given_name() {
        assert_eq!(surname_of("Thomas Michael Kane").as_deref(), Some("Kane"));
        assert_eq!(surname_of("Thomas Kane").as_deref(), Some("Kane"));
        assert_eq!(surname_of("Kane").as_deref(), Some("Kane"));
    }

    /// A surname that truly holds a space must survive whole. To refuse those would damage Dutch,
    /// Spanish and Gaelic lines, with no message, while the English ones passed.
    #[test]
    fn surname_keeps_its_particles() {
        assert_eq!(surname_of("Pieter van der Berg").as_deref(), Some("van der Berg"));
        assert_eq!(surname_of("Juan de la Cruz").as_deref(), Some("de la Cruz"));
        assert_eq!(surname_of("Angus Mac Donald").as_deref(), Some("Mac Donald"));
        assert_eq!(surname_of("Seán Ó Súilleabháin").as_deref(), Some("Ó Súilleabháin"));
    }

    /// `Surname, Given` is everywhere in genealogy exports, and the last-token rule would take
    /// exactly the wrong half of it.
    #[test]
    fn surname_handles_the_comma_first_form() {
        assert_eq!(surname_of("Kane, Thomas Michael").as_deref(), Some("Kane"));
        assert_eq!(surname_of("van der Berg, Pieter").as_deref(), Some("van der Berg"));
    }

    #[test]
    fn surname_drops_generational_suffixes() {
        assert_eq!(surname_of("Patrick O'Brien Jr.").as_deref(), Some("O'Brien"));
        assert_eq!(surname_of("John Smith III").as_deref(), Some("Smith"));
        assert_eq!(surname_of("John Smith Sr").as_deref(), Some("Smith"));
    }

    /// A third of the reference corpus appends dates or a birthplace to the name. Without the cut
    /// the last-token rule published `Balnicol` and `Louis` as surnames.
    #[test]
    fn surname_ignores_the_biographical_tail() {
        assert_eq!(
            surname_of("William Macaulay ~1770 of Balnicol, Uig, Lewis").as_deref(),
            Some("Macaulay")
        );
        assert_eq!(
            surname_of("Michael OConnell b1854 d1928 St Louis").as_deref(),
            Some("OConnell")
        );
        assert_eq!(
            surname_of("John A Moynihan 11/25/1843- 2/22/1920").as_deref(),
            Some("Moynihan")
        );
        assert_eq!(surname_of("Elizabeth Van Brunt b 1687").as_deref(), Some("Van Brunt"));
        assert_eq!(surname_of("Rosa Lee Aiken, 1919-1996").as_deref(), Some("Aiken"));
        // A particle surname must survive the cut intact.
        assert_eq!(
            surname_of("John De Lacy 1648/1722 - Ballingarry").as_deref(),
            Some("De Lacy")
        );
    }

    /// The FTDNA importer leaves HTML entities in the value. To send `mac L&#225;ire` out as a
    /// surname is worse than to decode it at the last point before it goes.
    #[test]
    fn surname_decodes_the_importers_html_entities() {
        assert_eq!(surname_of("Conall Corc mac L&#225;ire").as_deref(), Some("mac Láire"));
        assert_eq!(surname_of("Diarmaid &#211; Drisceoil").as_deref(), Some("Ó Drisceoil"));
        assert_eq!(surname_of("Jos&#233; de Mello").as_deref(), Some("de Mello"));
        // A bare ampersand stays as it is, and does not consume the rest of the string.
        assert_eq!(surname_of("Smith & Sons").as_deref(), Some("Sons"));
    }

    /// Nothing usable gives nothing. It never gives a stray fragment that would go out as a name.
    #[test]
    fn surname_of_nothing_is_none() {
        assert_eq!(surname_of(""), None);
        assert_eq!(surname_of("   "), None);
        assert_eq!(surname_of("Jr."), None, "a suffix alone is not a surname");
        assert_eq!(surname_of("?"), None);
        // Junk that the corpus holds. Nothing is better than a fragment that goes out as a name.
        assert_eq!(surname_of("trees.ancestry.com/tree/49418381/family"), None);
        assert_eq!(surname_of("1846 Duplin County, NC"), None);
        assert_eq!(surname_of("ABT. 1769 • Kilmalkedar, Co Kerry, Ireland"), None);
    }
}

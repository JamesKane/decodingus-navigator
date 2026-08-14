//! Which mapper preset a sample's reads need.
//!
//! minimap2's presets are not interchangeable tunings of one algorithm — `sr` and `map-ont` differ
//! in k-mer size, chaining, and gap costs by more than an order of magnitude in effect. Mapping
//! long reads under `sr` does not fail loudly; it produces plausible-looking, wrong alignments. So
//! the inference here refuses rather than guesses, which is the behaviour the design asks for
//! ("refuse (or warn loudly) on mixed/unknown technology rather than guessing").

use crate::error::AlignError;

/// A minimap2 preset, restricted to the ones this module maps reads under.
///
/// The assembly and splice presets exist upstream but have no meaning for realigning a
/// resequenced human sample, so they are deliberately not modelled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Preset {
    /// Illumina / short-read WGS and WES.
    ShortRead,
    /// PacBio HiFi.
    MapHifi,
    /// Oxford Nanopore.
    MapOnt,
}

impl Preset {
    /// The preset's minimap2 name, as `set_opt` and the CLI spell it.
    pub fn as_str(self) -> &'static str {
        match self {
            Preset::ShortRead => "sr",
            Preset::MapHifi => "map-hifi",
            Preset::MapOnt => "map-ont",
        }
    }

    /// Whether reads under this preset come in pairs. Only short-read data is mapped as pairs;
    /// long-read presets are single-end, and this is also what decides duplicate marking later
    /// (stage C marks duplicates for short reads only).
    pub fn is_paired(self) -> bool {
        matches!(self, Preset::ShortRead)
    }

    /// Parse an explicit user override.
    pub fn parse(name: &str) -> Result<Self, AlignError> {
        match name.trim().to_ascii_lowercase().as_str() {
            "sr" | "short" | "shortread" => Ok(Preset::ShortRead),
            "map-hifi" | "hifi" => Ok(Preset::MapHifi),
            "map-ont" | "ont" | "nanopore" => Ok(Preset::MapOnt),
            other => Err(AlignError::UnknownTechnology {
                what: format!("preset {other:?}"),
            }),
        }
    }

    /// Choose a preset from what the workspace already inferred about the run.
    ///
    /// `test_type` is `SequenceRun.test_type` (`WGS`, `WGS_HIFI`, `WGS_NANOPORE`, `WES`,
    /// `BIG_Y_700`, …) and `platform` is `SequenceRun.platform_name` (from `@RG PL`). The test type
    /// is consulted first because `testtype.rs` has already combined platform and read-length
    /// evidence to produce it; the platform is only a fallback for runs that never got one.
    ///
    /// Targeted panels (Big Y and friends) map under the chemistry that produced them, which for
    /// every such product Navigator ingests is short-read.
    pub fn infer(test_type: Option<&str>, platform: Option<&str>) -> Result<Self, AlignError> {
        if let Some(t) = test_type {
            let t = t.trim().to_ascii_uppercase();
            if !t.is_empty() {
                return match t.as_str() {
                    "WGS_HIFI" => Ok(Preset::MapHifi),
                    "WGS_NANOPORE" => Ok(Preset::MapOnt),
                    // WGS / WES / BIG_Y_* / any other targeted NGS product: short read.
                    "WGS" | "WES" => Ok(Preset::ShortRead),
                    other if other.starts_with("BIG_Y") => Ok(Preset::ShortRead),
                    other => Err(AlignError::UnknownTechnology {
                        what: format!("test type {other:?}"),
                    }),
                };
            }
        }

        // No test type recorded — fall back to the raw platform string.
        let p = platform.unwrap_or_default().to_ascii_uppercase();
        if p.contains("PACBIO") {
            // PacBio without a test type is ambiguous between CLR and HiFi. HiFi is what every
            // consumer PacBio product Navigator sees actually is, but the guess is worth flagging
            // rather than burying, so callers can surface it.
            return Ok(Preset::MapHifi);
        }
        if p.contains("NANOPORE") || p == "ONT" || p.contains("OXFORD") {
            return Ok(Preset::MapOnt);
        }
        if p.contains("ILLUMINA") || p.contains("BGI") || p.contains("MGI") || p.contains("DNBSEQ") {
            return Ok(Preset::ShortRead);
        }

        Err(AlignError::UnknownTechnology {
            what: if p.is_empty() {
                "a run with no test type or platform".to_string()
            } else {
                format!("platform {p:?}")
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_read_test_types_pick_their_own_presets() {
        assert_eq!(Preset::infer(Some("WGS_HIFI"), None).unwrap(), Preset::MapHifi);
        assert_eq!(Preset::infer(Some("WGS_NANOPORE"), None).unwrap(), Preset::MapOnt);
    }

    #[test]
    fn short_read_products_all_map_under_sr() {
        for t in ["WGS", "WES", "BIG_Y_700", "BIG_Y_500"] {
            assert_eq!(Preset::infer(Some(t), None).unwrap(), Preset::ShortRead, "{t}");
        }
    }

    /// The test type is the workspace's own considered inference, so it beats the raw platform
    /// string — a HiFi run still reports `PACBIO` as its platform.
    #[test]
    fn test_type_wins_over_platform() {
        assert_eq!(
            Preset::infer(Some("WGS_HIFI"), Some("ILLUMINA")).unwrap(),
            Preset::MapHifi
        );
    }

    #[test]
    fn platform_is_the_fallback_when_no_test_type_was_recorded() {
        assert_eq!(Preset::infer(None, Some("ILLUMINA")).unwrap(), Preset::ShortRead);
        assert_eq!(Preset::infer(Some(""), Some("DNBSEQ")).unwrap(), Preset::ShortRead);
        assert_eq!(Preset::infer(None, Some("OXFORD_NANOPORE")).unwrap(), Preset::MapOnt);
    }

    /// The property that matters: an unrecognized technology must stop the job. Mapping under the
    /// wrong preset yields wrong alignments quietly, which is worse than not running.
    #[test]
    fn an_unknown_technology_is_an_error_not_a_guess() {
        assert!(matches!(
            Preset::infer(None, None),
            Err(AlignError::UnknownTechnology { .. })
        ));
        assert!(matches!(
            Preset::infer(Some("WGS_SOMETHING_NEW"), None),
            Err(AlignError::UnknownTechnology { .. })
        ));
        assert!(matches!(
            Preset::infer(None, Some("SOLID")),
            Err(AlignError::UnknownTechnology { .. })
        ));
    }

    #[test]
    fn only_short_reads_are_mapped_as_pairs() {
        assert!(Preset::ShortRead.is_paired());
        assert!(!Preset::MapHifi.is_paired());
        assert!(!Preset::MapOnt.is_paired());
    }

    #[test]
    fn explicit_overrides_parse_by_name_and_alias() {
        assert_eq!(Preset::parse("sr").unwrap(), Preset::ShortRead);
        assert_eq!(Preset::parse("  MAP-ONT ").unwrap(), Preset::MapOnt);
        assert_eq!(Preset::parse("hifi").unwrap(), Preset::MapHifi);
        assert!(Preset::parse("asm5").is_err(), "not a resequencing preset");
    }
}

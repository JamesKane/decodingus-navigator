//! The DNA-test catalog — the kinds of test a subject can have (a `SequenceRun.test_type`
//! holds one of these codes). Ported from the Scala `test_types.conf` defaults: code,
//! display name, and the genomic region the test targets (which downstream gates Y/mt/
//! autosomal analysis). Static here; can move to a config file later as the Scala app does.

use serde::{Deserialize, Serialize};

/// What part of the genome a test covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TargetType {
    WholeGenome,
    YChromosome,
    MtDna,
    Autosomal,
    XChromosome,
    Mixed,
}

impl TargetType {
    pub fn label(self) -> &'static str {
        match self {
            TargetType::WholeGenome => "whole genome",
            TargetType::YChromosome => "Y chromosome",
            TargetType::MtDna => "mtDNA",
            TargetType::Autosomal => "autosomal",
            TargetType::XChromosome => "X chromosome",
            TargetType::Mixed => "mixed",
        }
    }
}

/// A known test type: a stable `code` (stored on the run), a human `display_name`, and the
/// `target` region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TestType {
    pub code: &'static str,
    pub display_name: &'static str,
    pub target: TargetType,
}

use TargetType::*;

/// The catalog (ported from `test_types.conf`). Order is presentation order.
pub const CATALOG: &[TestType] = &[
    TestType {
        code: "WGS",
        display_name: "Whole Genome Sequencing",
        target: WholeGenome,
    },
    TestType {
        code: "WGS_LOW_PASS",
        display_name: "Low-Pass WGS",
        target: WholeGenome,
    },
    TestType {
        code: "WGS_HIFI",
        display_name: "PacBio HiFi WGS",
        target: WholeGenome,
    },
    TestType {
        code: "WGS_NANOPORE",
        display_name: "Nanopore WGS",
        target: WholeGenome,
    },
    TestType {
        code: "WGS_CLR",
        display_name: "PacBio CLR WGS",
        target: WholeGenome,
    },
    TestType {
        code: "WES",
        display_name: "Whole Exome Sequencing",
        target: Autosomal,
    },
    TestType {
        code: "BIG_Y_500",
        display_name: "FTDNA Big Y-500",
        target: YChromosome,
    },
    TestType {
        code: "BIG_Y_700",
        display_name: "FTDNA Big Y-700",
        target: YChromosome,
    },
    TestType {
        code: "Y_ELITE",
        display_name: "Full Genomes Y Elite",
        target: YChromosome,
    },
    TestType {
        code: "Y_PRIME",
        display_name: "YSEQ Y Prime",
        target: YChromosome,
    },
    // Targeted tests recognized by coverage shape when the vendor can't be pinned down (see
    // `navigator-analysis::testtype::infer_test_type`) — honest generics, not a guessed product.
    TestType {
        code: "TARGETED_Y",
        display_name: "Targeted Y (vendor unknown)",
        target: YChromosome,
    },
    TestType {
        code: "TARGETED_MT",
        display_name: "Targeted mtDNA (vendor unknown)",
        target: MtDna,
    },
    TestType {
        code: "MT_FULL_SEQUENCE",
        display_name: "mtDNA Full Sequence",
        target: MtDna,
    },
    TestType {
        code: "MT_PLUS",
        display_name: "FTDNA mtDNA Plus",
        target: MtDna,
    },
    TestType {
        code: "MT_CR_ONLY",
        display_name: "mtDNA Control Region (HVR1/HVR2)",
        target: MtDna,
    },
    TestType {
        code: "YDNA_SNP_PACK_FTDNA",
        display_name: "FTDNA SNP Pack",
        target: YChromosome,
    },
    TestType {
        code: "YDNA_PANEL_YSEQ",
        display_name: "YSEQ Panel",
        target: YChromosome,
    },
    TestType {
        code: "YDNA_SNP_PANEL",
        display_name: "Y-DNA SNP Panel",
        target: YChromosome,
    },
    TestType {
        code: "ARRAY_BISDNA",
        display_name: "BISDNA Array",
        target: YChromosome,
    },
    TestType {
        code: "ARRAY_23ANDME_V5",
        display_name: "23andMe v5 Chip",
        target: Mixed,
    },
    TestType {
        code: "ARRAY_23ANDME_V4",
        display_name: "23andMe v4 Chip",
        target: Mixed,
    },
    TestType {
        code: "ARRAY_ANCESTRY_V2",
        display_name: "AncestryDNA v2",
        target: Mixed,
    },
    TestType {
        code: "ARRAY_FTDNA_FF",
        display_name: "FTDNA Family Finder",
        target: Autosomal,
    },
    TestType {
        code: "ARRAY_MYHERITAGE",
        display_name: "MyHeritage DNA",
        target: Mixed,
    },
    TestType {
        code: "ARRAY_LIVINGDNA",
        display_name: "LivingDNA",
        target: Mixed,
    },
];

/// Look a test type up by its code.
pub fn by_code(code: &str) -> Option<&'static TestType> {
    CATALOG.iter().find(|t| t.code == code)
}

/// Classify a stored `test_type` into its [`TargetType`] — tolerant of values that aren't a
/// canonical [`by_code`] code. A bulk import or a `--test-type` override may store a human label
/// like `"Big Y"` rather than `BIG_Y_500`/`BIG_Y_700`; without recognizing it the targeted-Y
/// scoping is lost and coverage walks the whole genome (slow on a targeted multi-reference CRAM).
/// Matches, in order: exact code, exact display name, then a small set of well-known vendor labels.
/// Returns `None` when nothing matches (caller treats that as whole-genome/unknown).
pub fn target_of(test_type: &str) -> Option<TargetType> {
    if let Some(t) = by_code(test_type) {
        return Some(t.target);
    }
    if let Some(t) = CATALOG.iter().find(|t| t.display_name.eq_ignore_ascii_case(test_type)) {
        return Some(t.target);
    }
    let s = test_type.trim().to_ascii_lowercase();
    const Y: &[&str] = &["big y", "big-y", "bigy", "y elite", "y-elite", "y prime", "y-prime", "targeted y"];
    const MT: &[&str] = &["mt full", "mtfull", "mt-full", "full mtdna", "full mitochondrial", "mtdna", "targeted mt"];
    if Y.iter().any(|p| s.contains(p)) {
        Some(TargetType::YChromosome)
    } else if MT.iter().any(|p| s.contains(p)) {
        Some(TargetType::MtDna)
    } else {
        None
    }
}

/// The display name for a code, falling back to the code itself if unknown.
pub fn display_name(code: &str) -> &str {
    by_code(code).map(|t| t.display_name).unwrap_or(code)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_by_code() {
        assert_eq!(by_code("BIG_Y_700").unwrap().target, YChromosome);
        assert_eq!(display_name("WGS"), "Whole Genome Sequencing");
        assert_eq!(display_name("NOPE"), "NOPE"); // unknown falls back to the code
    }

    #[test]
    fn target_of_tolerates_labels() {
        // Canonical code + display name.
        assert_eq!(target_of("BIG_Y_700"), Some(YChromosome));
        assert_eq!(target_of("FTDNA Big Y-500"), Some(YChromosome));
        assert_eq!(target_of("WGS"), Some(WholeGenome));
        // Human label a bulk import / --test-type override may have stored.
        assert_eq!(target_of("Big Y"), Some(YChromosome));
        assert_eq!(target_of("mtFull"), Some(MtDna));
        // Unknown → None (treated as whole-genome/unknown by callers).
        assert_eq!(target_of("Sanger panel"), None);
    }

    #[test]
    fn codes_are_unique() {
        let mut codes: Vec<_> = CATALOG.iter().map(|t| t.code).collect();
        let n = codes.len();
        codes.sort_unstable();
        codes.dedup();
        assert_eq!(codes.len(), n, "duplicate test-type code in CATALOG");
    }
}

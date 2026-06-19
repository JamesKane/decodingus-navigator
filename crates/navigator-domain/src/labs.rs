//! Catalog of known labs, sequencing centers, and genotyping vendors — a Rust port of the Scala
//! `LabsConfig`/`labs.conf`. Provides display names, ≤6-char abbreviations, categories, and
//! capabilities for the Data Sources UI (lab chips + the sequence-run lab dropdown), and
//! case-insensitive lookup by id / display name / alias for matching an inferred facility.
//!
//! Static data (no config file): the set is small and stable; adding a lab is a code change.

/// A lab / sequencing center / vendor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Lab {
    /// Stable internal id (the old HOCON key).
    pub id: &'static str,
    pub display_name: &'static str,
    /// Short code (≤6 chars) for UI chips.
    pub abbreviation: &'static str,
    /// `commercial-lab` | `consumer-vendor` | `sequencing-platform` | `academic`.
    pub category: &'static str,
    pub capabilities: &'static [&'static str],
    pub website: Option<&'static str>,
    pub aliases: &'static [&'static str],
}

/// Category constants.
pub mod category {
    pub const COMMERCIAL_LAB: &str = "commercial-lab";
    pub const CONSUMER_VENDOR: &str = "consumer-vendor";
    pub const SEQUENCING_PLATFORM: &str = "sequencing-platform";
    pub const ACADEMIC: &str = "academic";
}

/// The full catalog (25 labs), ported from `labs.conf`.
pub const CATALOG: &[Lab] = &[
    // Commercial DNA testing labs
    Lab {
        id: "familytreedna",
        display_name: "FamilyTreeDNA",
        abbreviation: "FTDNA",
        category: category::COMMERCIAL_LAB,
        capabilities: &["y-dna", "mt-dna", "str", "chip", "vcf-download", "bam-download"],
        website: Some("https://www.familytreedna.com"),
        aliases: &["FTDNA", "Family Tree DNA"],
    },
    Lab {
        id: "yseq",
        display_name: "YSEQ",
        abbreviation: "YSEQ",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs", "y-dna", "mt-dna", "str", "vcf-download", "bam-download"],
        website: Some("https://www.yseq.net"),
        aliases: &[],
    },
    Lab {
        id: "full-genomes",
        display_name: "Full Genomes Corporation",
        abbreviation: "FGC",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs", "y-dna", "mt-dna", "vcf-download"],
        website: Some("https://www.fullgenomes.com"),
        aliases: &["Full Genomes", "FGC"],
    },
    Lab {
        id: "nebula",
        display_name: "Nebula Genomics",
        abbreviation: "NEBULA",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs", "vcf-download"],
        website: Some("https://nebula.org"),
        aliases: &["Nebula"],
    },
    Lab {
        id: "dante",
        display_name: "Dante Labs",
        abbreviation: "DANTE",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs", "vcf-download"],
        website: Some("https://www.dantelabs.com"),
        aliases: &["Dante"],
    },
    Lab {
        id: "sequencing-com",
        display_name: "Sequencing.com",
        abbreviation: "SEQ",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs"],
        website: Some("https://sequencing.com"),
        aliases: &["Sequencing"],
    },
    Lab {
        id: "bisdna",
        display_name: "BISDNA",
        abbreviation: "BISDNA",
        category: category::COMMERCIAL_LAB,
        capabilities: &["y-dna", "chip"],
        website: None,
        aliases: &["BIS-DNA", "BIS DNA"],
    },
    Lab {
        id: "yoogene",
        display_name: "YooGene",
        abbreviation: "YOOGN",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs"],
        website: Some("https://www.yoogene.com"),
        aliases: &[],
    },
    Lab {
        id: "invitae",
        display_name: "Invitae",
        abbreviation: "INVIT",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs"],
        website: Some("https://www.invitae.com"),
        aliases: &[],
    },
    Lab {
        id: "genedx",
        display_name: "GeneDx",
        abbreviation: "GENEDX",
        category: category::COMMERCIAL_LAB,
        capabilities: &["wgs"],
        website: Some("https://www.genedx.com"),
        aliases: &[],
    },
    Lab {
        id: "yfull",
        display_name: "YFull",
        abbreviation: "YFULL",
        category: category::COMMERCIAL_LAB,
        capabilities: &["y-dna", "mt-dna"],
        website: Some("https://www.yfull.com"),
        aliases: &[],
    },
    // Consumer genotyping vendors
    Lab {
        id: "23andme",
        display_name: "23andMe",
        abbreviation: "23&ME",
        category: category::CONSUMER_VENDOR,
        capabilities: &["chip", "y-dna", "mt-dna"],
        website: Some("https://www.23andme.com"),
        aliases: &["23 and Me", "23andMe"],
    },
    Lab {
        id: "ancestrydna",
        display_name: "AncestryDNA",
        abbreviation: "ANCST",
        category: category::CONSUMER_VENDOR,
        capabilities: &["chip", "y-dna", "mt-dna"],
        website: Some("https://www.ancestry.com/dna"),
        aliases: &["Ancestry DNA", "AncestryDNA", "Ancestry"],
    },
    Lab {
        id: "myheritage",
        display_name: "MyHeritage",
        abbreviation: "MYHRT",
        category: category::CONSUMER_VENDOR,
        capabilities: &["chip", "y-dna", "mt-dna"],
        website: Some("https://www.myheritage.com"),
        aliases: &["MyHeritage DNA", "My Heritage"],
    },
    Lab {
        id: "livingdna",
        display_name: "LivingDNA",
        abbreviation: "LVDNA",
        category: category::CONSUMER_VENDOR,
        capabilities: &["chip", "y-dna", "mt-dna"],
        website: Some("https://livingdna.com"),
        aliases: &["Living DNA"],
    },
    // Sequencing platform vendors
    Lab {
        id: "illumina",
        display_name: "Illumina",
        abbreviation: "ILLUM",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://www.illumina.com"),
        aliases: &[],
    },
    Lab {
        id: "pacbio",
        display_name: "PacBio",
        abbreviation: "PACB",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://www.pacb.com"),
        aliases: &["Pacific Biosciences"],
    },
    Lab {
        id: "oxford-nanopore",
        display_name: "Oxford Nanopore",
        abbreviation: "ONT",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://nanoporetech.com"),
        aliases: &["Nanopore", "ONT"],
    },
    Lab {
        id: "mgi",
        display_name: "MGI Tech",
        abbreviation: "MGI",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://www.mgi-tech.com"),
        aliases: &["MGI", "BGI Genomics"],
    },
    Lab {
        id: "element",
        display_name: "Element Biosciences",
        abbreviation: "ELEM",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://www.elementbiosciences.com"),
        aliases: &["Element"],
    },
    Lab {
        id: "ultima",
        display_name: "Ultima Genomics",
        abbreviation: "ULTIMA",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://www.ultimagenomics.com"),
        aliases: &["Ultima"],
    },
    Lab {
        id: "singular",
        display_name: "Singular Genomics",
        abbreviation: "SINGLR",
        category: category::SEQUENCING_PLATFORM,
        capabilities: &["wgs"],
        website: Some("https://singulargenomics.com"),
        aliases: &[],
    },
    // Academic / research institutions
    Lab {
        id: "broad-institute",
        display_name: "Broad Institute",
        abbreviation: "BROAD",
        category: category::ACADEMIC,
        capabilities: &["wgs"],
        website: Some("https://www.broadinstitute.org"),
        aliases: &["Broad"],
    },
    Lab {
        id: "sanger",
        display_name: "Wellcome Sanger Institute",
        abbreviation: "SANGR",
        category: category::ACADEMIC,
        capabilities: &["wgs"],
        website: Some("https://www.sanger.ac.uk"),
        aliases: &["Sanger", "Wellcome Trust Sanger"],
    },
    Lab {
        id: "nhgri",
        display_name: "NHGRI",
        abbreviation: "NHGRI",
        category: category::ACADEMIC,
        capabilities: &["wgs"],
        website: Some("https://www.genome.gov"),
        aliases: &["National Human Genome Research Institute"],
    },
];

/// Find a lab by id, display name, or alias (case-insensitive), then by a fuzzy substring match
/// over display names and aliases. `None` if nothing matches.
pub fn find(identifier: &str) -> Option<&'static Lab> {
    let q = identifier.trim().to_ascii_lowercase();
    if q.is_empty() {
        return None;
    }
    CATALOG
        .iter()
        .find(|l| {
            l.id == q
                || l.display_name.to_ascii_lowercase() == q
                || l.aliases.iter().any(|a| a.to_ascii_lowercase() == q)
        })
        .or_else(|| {
            CATALOG.iter().find(|l| {
                l.display_name.to_ascii_lowercase().contains(&q)
                    || l.aliases.iter().any(|a| a.to_ascii_lowercase().contains(&q))
            })
        })
}

/// The ≤`max_chars`-char abbreviation for a lab name (configured, else the name's first chars
/// upper-cased).
pub fn abbreviation(name: &str, max_chars: usize) -> String {
    match find(name) {
        Some(lab) => lab.abbreviation.chars().take(max_chars).collect(),
        None => name.chars().take(max_chars).collect::<String>().to_uppercase(),
    }
}

/// The canonical display name for an identifier (returns the input unchanged if unknown).
pub fn display_name(identifier: &str) -> String {
    find(identifier)
        .map(|l| l.display_name.to_string())
        .unwrap_or_else(|| identifier.to_string())
}

/// Display names offered in the sequence-run lab dropdown: testing labs + sequencing platforms +
/// academic institutions (not consumer-array vendors), sorted.
pub fn sequence_run_lab_names() -> Vec<&'static str> {
    let mut names: Vec<&'static str> = CATALOG
        .iter()
        .filter(|l| {
            matches!(
                l.category,
                category::COMMERCIAL_LAB | category::SEQUENCING_PLATFORM | category::ACADEMIC
            )
        })
        .map(|l| l.display_name)
        .collect();
    names.sort_unstable();
    names
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_by_id_display_alias_and_fuzzy() {
        assert_eq!(find("yseq").unwrap().display_name, "YSEQ");
        assert_eq!(find("Full Genomes Corporation").unwrap().id, "full-genomes");
        assert_eq!(find("FGC").unwrap().id, "full-genomes"); // alias
        assert_eq!(find("ftdna").unwrap().display_name, "FamilyTreeDNA"); // alias, case-insensitive
        assert_eq!(find("dante labs").unwrap().id, "dante"); // exact display (lower)
        assert_eq!(find("nebula genomics").unwrap().id, "nebula");
        assert!(find("totally-unknown-lab").is_none());
    }

    #[test]
    fn abbreviation_falls_back_to_uppercased_prefix() {
        assert_eq!(abbreviation("YSEQ", 6), "YSEQ");
        assert_eq!(abbreviation("Full Genomes Corporation", 6), "FGC");
        assert_eq!(abbreviation("Mystery Lab", 4), "MYST");
    }

    #[test]
    fn run_dropdown_excludes_consumer_vendors_and_is_sorted() {
        let names = sequence_run_lab_names();
        assert!(names.contains(&"YSEQ") && names.contains(&"Illumina") && names.contains(&"Broad Institute"));
        assert!(
            !names.contains(&"23andMe"),
            "consumer array vendors aren't sequence-run labs"
        );
        assert!(names.windows(2).all(|w| w[0] <= w[1]), "sorted");
    }
}

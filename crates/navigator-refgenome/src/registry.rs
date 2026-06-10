//! The known-build registry: canonical build names + their source URLs (reference FASTAs
//! and liftover chains), with an optional user JSON override. Defaults come from
//! `docs/chm13-reference-resources.md` (CHM13 assets on the public human-pangenomics bucket)
//! plus the Broad public GRCh38/GRCh37 references.

use std::collections::HashMap;
use std::path::Path;

use serde::Deserialize;

/// A reference assembly Navigator can resolve.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Build {
    Grch38,
    Grch37,
    Chm13v2,
    /// The CHM13v2.0 analysis set with the Y PAR hard-masked and the mitochondrion replaced
    /// by rCRS — the recommended short-read calling reference (PAR-masking removes X/Y
    /// multi-mapping artifacts; rCRS chrM matches the haplotree coordinates). Its **nuclear**
    /// coordinates are identical to [`Build::Chm13v2`] (see [`Build::nuclear`]); only chrY is
    /// N-masked in the PAR and chrM is swapped — so it reuses CHM13's liftover chains.
    Chm13v2MaskedRcrs,
}

impl Build {
    /// Canonical label, also the cache filename stem and the value `reference_build_for`
    /// stamps on alignments.
    pub fn as_str(self) -> &'static str {
        match self {
            Build::Grch38 => "GRCh38",
            Build::Grch37 => "GRCh37",
            Build::Chm13v2 => "chm13v2.0",
            Build::Chm13v2MaskedRcrs => "chm13v2.0_maskedY_rCRS",
        }
    }

    /// The build whose **nuclear coordinate system** this one shares — itself for the plain
    /// assemblies, and [`Build::Chm13v2`] for the masked+rCRS variant (PAR-masking is N's, not
    /// a coordinate change). Liftover chains key off this, so the masked variant reuses
    /// CHM13's chains instead of duplicating them. (chrM differs — masked chrM is rCRS — but
    /// mtDNA is never lifted via a chain; it is a direct rCRS query.)
    pub fn nuclear(self) -> Build {
        match self {
            Build::Chm13v2MaskedRcrs => Build::Chm13v2,
            other => other,
        }
    }

    /// Provenance of this reference's haploid sequences — a standing reminder that the
    /// reference allele is a *coordinate system*, never a source of ancestral/derived
    /// polarity. See [`ReferencePolarity`].
    pub fn reference_polarity(self) -> ReferencePolarity {
        const RCRS_M: &str =
            "rCRS (NC_012920.1, haplogroup H2a2a1) — itself derived from the RSRS root, not ancestral";
        match self {
            Build::Chm13v2 => ReferencePolarity {
                chr_y: "HG002 Y, haplogroup J — the reference base is the DERIVED allele at many Y-SNP sites",
                chr_m: "CHM13's own mitochondrion (NOT rCRS) — handled via the rotation-aware rCRS↔chrM map",
            },
            Build::Chm13v2MaskedRcrs => ReferencePolarity {
                chr_y: "HG002 Y, haplogroup J (PAR hard-masked) — the reference base is the DERIVED allele at many Y-SNP sites",
                chr_m: RCRS_M,
            },
            Build::Grch38 => ReferencePolarity {
                chr_y: "GRCh38 chrY — a specific donor's Y, not the ancestral root",
                chr_m: RCRS_M,
            },
            Build::Grch37 => ReferencePolarity {
                chr_y: "GRCh37 chrY — a specific donor's Y, not the ancestral root",
                chr_m: RCRS_M,
            },
        }
    }
}

/// The provenance of a reference's haploid (chrY / chrM) sequences. It exists to make one
/// invariant explicit and discoverable: **ancestral/derived polarity must always come from
/// the haplotree, compared against the sample's own called base — never from "is the sample's
/// base equal to the reference (REF) or not (ALT)."** The canonical trap is CHM13v2.0, whose
/// chrY is HG002 (a haplogroup-J Y): the reference base is the *derived* allele at many Y-SNP
/// sites, so a REF-as-ancestral assumption would invert those calls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferencePolarity {
    /// What the reference's chrY is, and why its allele is not a polarity source.
    pub chr_y: &'static str,
    /// What the reference's mitochondrion is, and its polarity caveat.
    pub chr_m: &'static str,
}

/// Map any common spelling of a build to its canonical [`Build`] (case-insensitive).
pub fn canonical_build(name: &str) -> Option<Build> {
    let n = name.to_ascii_lowercase();
    match n.as_str() {
        "grch38" | "hg38" | "b38" | "grch38.p14" => Some(Build::Grch38),
        "grch37" | "hg19" | "b37" | "grch37.p13" => Some(Build::Grch37),
        "chm13v2.0_maskedy_rcrs" | "chm13v2_maskedy_rcrs" | "chm13_maskedy_rcrs"
        | "chm13v2.0-maskedy-rcrs" | "chm13v2.0_masked_rcrs" => Some(Build::Chm13v2MaskedRcrs),
        "chm13" | "chm13v2" | "chm13v2.0" | "t2t" | "hs1" | "t2t-chm13v2.0" => Some(Build::Chm13v2),
        _ => None,
    }
}

/// Where a reference FASTA is fetched from, with a rough size for the download prompt.
#[derive(Debug, Clone)]
pub struct ReferenceSource {
    pub build: Build,
    pub url: String,
    pub est_bytes: u64,
}

/// Where a liftover chain (UCSC `.chain`, 1:1) is fetched from.
#[derive(Debug, Clone)]
pub struct ChainSource {
    pub from: Build,
    pub to: Build,
    pub url: String,
}

const GB: u64 = 1_000_000_000;
const CHM13_FA: &str =
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz";
const CHM13_MASKED_RCRS_FA: &str =
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz";
const GRCH38_FA: &str =
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta";
const GRCH37_FA: &str =
    "https://storage.googleapis.com/genomics-public-data/references/hg19/v0/Homo_sapiens_assembly19.fasta.gz";
const CHAIN_BASE: &str = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo";

/// Per-build user override loaded from `reference_sources.json`.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct BuildOverride {
    /// Use this local FASTA as-is (already decompressed + indexed); never download.
    pub local_path: Option<String>,
    /// Override the download URL.
    pub url: Option<String>,
}

/// The optional user config at `~/.decodingus/config/reference_sources.json`.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct UserConfig {
    #[serde(default)]
    pub references: HashMap<String, BuildOverride>,
}

impl UserConfig {
    /// Load the config if present; a missing or unreadable file yields the empty default
    /// (overrides are advisory, never fatal).
    pub fn load(path: &Path) -> Self {
        std::fs::read_to_string(path)
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or_default()
    }

    fn for_build(&self, build: Build) -> Option<&BuildOverride> {
        self.references.get(build.as_str())
    }
}

/// The resolved set of sources (defaults merged with any user overrides).
#[derive(Debug, Clone)]
pub struct Registry {
    config: UserConfig,
}

impl Registry {
    pub fn new(config: UserConfig) -> Self {
        Registry { config }
    }

    /// A user-pinned local FASTA for this build, if configured.
    pub fn local_override(&self, build: Build) -> Option<&str> {
        self.config.for_build(build).and_then(|o| o.local_path.as_deref())
    }

    /// The download source for a build (user URL override else the built-in default).
    pub fn reference_source(&self, build: Build) -> ReferenceSource {
        let (default_url, est_bytes) = match build {
            Build::Grch38 => (GRCH38_FA, 3 * GB + GB / 10),
            Build::Grch37 => (GRCH37_FA, 9 * GB / 10),
            Build::Chm13v2 => (CHM13_FA, GB),
            Build::Chm13v2MaskedRcrs => (CHM13_MASKED_RCRS_FA, GB),
        };
        let url = self
            .config
            .for_build(build)
            .and_then(|o| o.url.clone())
            .unwrap_or_else(|| default_url.to_string());
        ReferenceSource { build, url, est_bytes }
    }

    /// The liftover chain source for a build pair, if one is registered. Builds are
    /// normalized to their nuclear coordinate system first, so the masked+rCRS variant reuses
    /// CHM13's chains (its nuclear coordinates are identical).
    pub fn chain_source(&self, from: Build, to: Build) -> Option<ChainSource> {
        let (from, to) = (from.nuclear(), to.nuclear());
        let file = match (from, to) {
            (Build::Grch38, Build::Chm13v2) => "grch38-chm13v2.chain",
            (Build::Chm13v2, Build::Grch38) => "chm13v2-grch38.chain",
            (Build::Grch37, Build::Chm13v2) => "hg19-chm13v2.chain",
            (Build::Chm13v2, Build::Grch37) => "chm13v2-hg19.chain",
            _ => return None,
        };
        Some(ChainSource { from, to, url: format!("{CHAIN_BASE}/{file}") })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_build_accepts_common_aliases() {
        assert_eq!(canonical_build("chm13v2.0"), Some(Build::Chm13v2));
        assert_eq!(canonical_build("CHM13"), Some(Build::Chm13v2));
        assert_eq!(canonical_build("hs1"), Some(Build::Chm13v2));
        assert_eq!(canonical_build("hg38"), Some(Build::Grch38));
        assert_eq!(canonical_build("GRCh37"), Some(Build::Grch37));
        assert_eq!(canonical_build("b37"), Some(Build::Grch37));
        assert_eq!(canonical_build("unknown"), None);
    }

    #[test]
    fn default_sources_resolve() {
        let reg = Registry::new(UserConfig::default());
        assert!(reg.reference_source(Build::Chm13v2).url.ends_with("chm13v2.0.fa.gz"));
        assert!(reg.local_override(Build::Chm13v2).is_none());
        let chain = reg.chain_source(Build::Grch38, Build::Chm13v2).unwrap();
        assert!(chain.url.ends_with("grch38-chm13v2.chain"));
        assert!(reg.chain_source(Build::Grch38, Build::Grch37).is_none());
    }

    #[test]
    fn masked_rcrs_is_a_resolvable_cacheable_build() {
        // Aliases canonicalize to the masked variant, and its as_str round-trips (so the cache
        // file is `chm13v2.0_maskedY_rCRS.fa` — distinct from plain chm13).
        for alias in ["chm13v2.0_maskedY_rCRS", "chm13_maskedY_rcrs", "CHM13V2.0_MASKEDY_RCRS"] {
            assert_eq!(canonical_build(alias), Some(Build::Chm13v2MaskedRcrs), "alias {alias}");
        }
        assert_eq!(Build::Chm13v2MaskedRcrs.as_str(), "chm13v2.0_maskedY_rCRS");
        assert_eq!(canonical_build(Build::Chm13v2MaskedRcrs.as_str()), Some(Build::Chm13v2MaskedRcrs));
        // Plain chm13 spellings still map to plain chm13.
        assert_eq!(canonical_build("chm13v2.0"), Some(Build::Chm13v2));

        let reg = Registry::new(UserConfig::default());
        assert!(reg.reference_source(Build::Chm13v2MaskedRcrs).url.ends_with("chm13v2.0_maskedY_rCRS.fa.gz"));
    }

    #[test]
    fn reference_polarity_records_the_chm13_y_is_j_trap() {
        // CHM13 (and the masked variant) carry HG002's haplogroup-J Y: the reference base is
        // derived, so the metadata must flag it as not-ancestral.
        for b in [Build::Chm13v2, Build::Chm13v2MaskedRcrs] {
            let p = b.reference_polarity();
            assert!(p.chr_y.contains("HG002") && p.chr_y.contains('J'), "{}: {}", b.as_str(), p.chr_y);
            assert!(p.chr_y.contains("DERIVED"));
        }
        // The masked variant's mito is rCRS; plain CHM13's is its own (not rCRS).
        assert!(Build::Chm13v2MaskedRcrs.reference_polarity().chr_m.contains("rCRS"));
        assert!(Build::Chm13v2.reference_polarity().chr_m.contains("NOT rCRS"));
        // Every reference documents an mt polarity caveat (rCRS is itself derived).
        for b in [Build::Grch38, Build::Grch37, Build::Chm13v2MaskedRcrs] {
            assert!(b.reference_polarity().chr_m.contains("rCRS"));
        }
    }

    #[test]
    fn masked_rcrs_shares_chm13_nuclear_coords_and_chains() {
        // Nuclear coordinate system is CHM13's; chrM (rCRS) is never chain-lifted.
        assert_eq!(Build::Chm13v2MaskedRcrs.nuclear(), Build::Chm13v2);
        assert_eq!(Build::Chm13v2.nuclear(), Build::Chm13v2);

        // So the masked variant reuses CHM13's chains — same file, normalized endpoints (no
        // duplicate keyed by the masked name).
        let reg = Registry::new(UserConfig::default());
        let direct = reg.chain_source(Build::Grch38, Build::Chm13v2).unwrap();
        let masked = reg.chain_source(Build::Grch38, Build::Chm13v2MaskedRcrs).unwrap();
        assert_eq!(masked.url, direct.url);
        assert_eq!(masked.to, Build::Chm13v2); // normalized for cache-key reuse
        assert!(reg.chain_source(Build::Chm13v2MaskedRcrs, Build::Grch38).unwrap().url.ends_with("chm13v2-grch38.chain"));
    }

    #[test]
    fn user_override_wins() {
        let mut references = HashMap::new();
        references.insert(
            "chm13v2.0".to_string(),
            BuildOverride { local_path: Some("/data/chm13.fa".into()), url: None },
        );
        let reg = Registry::new(UserConfig { references });
        assert_eq!(reg.local_override(Build::Chm13v2), Some("/data/chm13.fa"));
    }
}

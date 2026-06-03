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
}

impl Build {
    /// Canonical label, also the cache filename stem and the value `reference_build_for`
    /// stamps on alignments.
    pub fn as_str(self) -> &'static str {
        match self {
            Build::Grch38 => "GRCh38",
            Build::Grch37 => "GRCh37",
            Build::Chm13v2 => "chm13v2.0",
        }
    }
}

/// Map any common spelling of a build to its canonical [`Build`] (case-insensitive).
pub fn canonical_build(name: &str) -> Option<Build> {
    let n = name.to_ascii_lowercase();
    match n.as_str() {
        "grch38" | "hg38" | "b38" | "grch38.p14" => Some(Build::Grch38),
        "grch37" | "hg19" | "b37" | "grch37.p13" => Some(Build::Grch37),
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
        };
        let url = self
            .config
            .for_build(build)
            .and_then(|o| o.url.clone())
            .unwrap_or_else(|| default_url.to_string());
        ReferenceSource { build, url, est_bytes }
    }

    /// The liftover chain source for a build pair, if one is registered.
    pub fn chain_source(&self, from: Build, to: Build) -> Option<ChainSource> {
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

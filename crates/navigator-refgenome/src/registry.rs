//! The known-build registry: canonical build names + their source URLs (reference FASTAs
//! and liftover chains), with an optional user JSON override. Defaults come from
//! `docs/chm13-reference-resources.md` (CHM13 assets on the public human-pangenomics bucket)
//! plus the Broad public GRCh38/GRCh37 references.

use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

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
    /// Every supported build, in display order.
    pub fn all() -> &'static [Build] {
        &[Build::Grch38, Build::Grch37, Build::Chm13v2, Build::Chm13v2MaskedRcrs]
    }

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
        const RCRS_M: &str = "rCRS (NC_012920.1, haplogroup H2a2a1) — itself derived from the RSRS root, not ancestral";
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
        "chm13v2.0_maskedy_rcrs"
        | "chm13v2_maskedy_rcrs"
        | "chm13_maskedy_rcrs"
        | "chm13v2.0-maskedy-rcrs"
        | "chm13v2.0_masked_rcrs" => Some(Build::Chm13v2MaskedRcrs),
        "chm13" | "chm13v2" | "chm13v2.0" | "t2t" | "hs1" | "t2t-chm13v2.0" => Some(Build::Chm13v2),
        _ => None,
    }
}

/// Where a reference FASTA is fetched from, with a rough size for the download prompt and an
/// optional pinned SHA-256 of the downloaded artifact (publisher's hash, when known) used to
/// verify the download before it's accepted. `None` = no authoritative hash to pin against yet.
#[derive(Debug, Clone)]
pub struct ReferenceSource {
    pub build: Build,
    pub url: String,
    pub est_bytes: u64,
    pub sha256: Option<String>,
}

/// Where a liftover chain (UCSC `.chain`, 1:1) is fetched from, with an optional pinned SHA-256.
#[derive(Debug, Clone)]
pub struct ChainSource {
    pub from: Build,
    pub to: Build,
    pub url: String,
    pub sha256: Option<String>,
}

/// Where a named annotation-mask BED is fetched from (e.g. the curated CHM13 Y structural
/// regions). `name` is the cache key / filename stem. Optional pinned SHA-256.
#[derive(Debug, Clone)]
pub struct MaskSource {
    pub name: String,
    pub url: String,
    pub sha256: Option<String>,
}

/// The curated CHM13v2.0 chrY structural-region BEDs (marbl/CHM13, Rhie et al. 2023) — the
/// paralog-prone zones used to flag unreliable Y calls. Keyed by cache-stable name. All are
/// on the human-pangenomics bucket alongside the references and chains.
pub const Y_STRUCTURAL_MASKS: &[(&str, &str)] = &[
    ("chm13v2.0Y_inverted_repeats_v1", "chm13v2.0Y_inverted_repeats_v1.bed"),
    ("chm13v2.0Y_amplicons_v1", "chm13v2.0Y_amplicons_v1.bed"),
    ("chm13v2.0Y_AZF_DYZ_v1", "chm13v2.0Y_AZF_DYZ_v1.bed"),
];

const GB: u64 = 1_000_000_000;
const CHM13_FA: &str =
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz";
const CHM13_MASKED_RCRS_FA: &str =
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz";
// GRCh38: the community-standard "no-ALT analysis set" (Heng Li,
// https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use), served bgzipped by NCBI's
// HTTPS mirror. ~873 MB compressed vs. ~3.25 GB for the Broad plain FASTA we used before — a ~3.7x
// win on a slow/cellular connection (the `.gz` is decompressed + indexed locally by
// `decompress_and_index`). `chr`-prefixed contig names, rCRS chrM — matches the app's GRCh38
// convention. No ALT/decoy/HLA (Heng Li's recommended analysis set); a CRAM aligned to a
// Broad-specific decoy/ALT contig would need a `reference_sources.json` URL override.
const GRCH38_FA: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz";
// GRCh37: hs37-1kg (1000 Genomes phase-1 reference, Heng Li's recommendation), gzipped on the EBI
// 1000genomes HTTPS mirror. ~892 MB compressed vs. ~3.14 GB for the Broad plain FASTA. Bare contig
// names (1/X/MT), rCRS MT — matches the app's GRCh37 convention.
const GRCH37_FA: &str = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz";
const CHAIN_BASE: &str = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo";
const ANNOTATION_BASE: &str = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation";

/// Built-in authoritative SHA-256 (lowercase hex) of a reference FASTA's **downloaded artifact**
/// (the `.fa.gz` / `.fasta` exactly as served), when a publisher checksum has been confirmed.
/// `None` until verified — the integrity machinery ships ready and pins fill in over time. Add
/// values here (or via the per-build `sha256` override in `reference_sources.json`).
fn default_reference_sha(build: Build) -> Option<&'static str> {
    match build {
        // Awaiting confirmed publisher checksums (T2T human-pangenomics bucket / Broad references).
        Build::Grch38 | Build::Grch37 | Build::Chm13v2 | Build::Chm13v2MaskedRcrs => None,
    }
}

/// Per-build user override loaded from `reference_sources.json`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, PartialEq, Eq)]
pub struct BuildOverride {
    /// Use this local FASTA as-is (already decompressed + indexed); never download.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub local_path: Option<String>,
    /// Override the download URL.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub url: Option<String>,
    /// Pin an authoritative SHA-256 (lowercase hex) of the downloaded artifact; the download is
    /// rejected if it doesn't match. Lets a user supply a publisher checksum we don't ship.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub sha256: Option<String>,
    /// Whether a missing reference may be auto-downloaded for this build (default `true`).
    #[serde(default = "default_true")]
    pub auto_download: bool,
}

fn default_true() -> bool {
    true
}

/// The optional user config at `~/.decodingus/config/reference_sources.json`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, PartialEq, Eq)]
pub struct UserConfig {
    #[serde(default)]
    pub references: HashMap<String, BuildOverride>,
}

impl UserConfig {
    /// Load the config if present; a missing or unreadable file yields the empty default (overrides
    /// are advisory, never fatal — a novice with no config just gets the self-managed auto-download).
    ///
    /// A file that **exists but doesn't parse** also falls back to defaults, but is **warned about**:
    /// silently dropping it is how a power user's `local_path` override vanishes and the app surprises
    /// them with a full reference download (issue #26 — the config had been corrupted by a racing
    /// non-atomic write; see [`crate::cache::atomic_write`]). Say so instead of reverting in silence.
    pub fn load(path: &Path) -> Self {
        let Ok(text) = std::fs::read_to_string(path) else {
            return Self::default(); // absent / unreadable → empty (the normal no-config case)
        };
        match serde_json::from_str(&text) {
            Ok(cfg) => cfg,
            Err(e) => {
                eprintln!(
                    "reference_sources.json at {} is invalid ({e}) — ignoring it and using the default \
                     (auto-download) sources. Your reference overrides are NOT being applied; fix or \
                     delete the file to restore them.",
                    path.display()
                );
                Self::default()
            }
        }
    }

    /// Persist to `path` (creating the parent `config/` dir), pretty-printed. Written **atomically**
    /// (temp + rename, see [`crate::cache::atomic_write`]) — this file is rewritten from spawned worker
    /// tasks that can race, and a plain non-atomic write corrupts it into head-of-new + tail-of-old
    /// garbage. Callers should still avoid concurrent read-modify-write (prefer one bulk save) so an
    /// update isn't lost; atomicity only guarantees the file is never *torn*.
    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        let json = serde_json::to_string_pretty(self).map_err(std::io::Error::other)?;
        crate::cache::atomic_write(path, json.as_bytes())
    }

    /// The override for `build`, if any.
    pub fn for_build(&self, build: Build) -> Option<&BuildOverride> {
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
            Build::Grch38 => (GRCH38_FA, 873 * GB / 1000), // ~873 MB bgzipped no-ALT analysis set
            Build::Grch37 => (GRCH37_FA, 892 * GB / 1000), // ~892 MB gzipped hs37-1kg
            Build::Chm13v2 => (CHM13_FA, GB),
            Build::Chm13v2MaskedRcrs => (CHM13_MASKED_RCRS_FA, GB),
        };
        let ov = self.config.for_build(build);
        let url = ov
            .and_then(|o| o.url.clone())
            .unwrap_or_else(|| default_url.to_string());
        // User-pinned hash wins over the built-in (which is None until a publisher hash is confirmed).
        let sha256 = ov
            .and_then(|o| o.sha256.clone())
            .or_else(|| default_reference_sha(build).map(str::to_string));
        ReferenceSource {
            build,
            url,
            est_bytes,
            sha256,
        }
    }

    /// The liftover chain source for a build pair, if one is registered. Builds are
    /// normalized to their nuclear coordinate system first, so the masked+rCRS variant reuses
    /// CHM13's chains (its nuclear coordinates are identical).
    pub fn chain_source(&self, from: Build, to: Build) -> Option<ChainSource> {
        let (from, to) = (from.nuclear(), to.nuclear());
        // GRCh38↔GRCh37 use UCSC's gzipped over.chain (decompressed on download); the CHM13 pairs
        // are the curated uncompressed chains in the T2T bucket.
        let url = match (from, to) {
            (Build::Grch38, Build::Chm13v2) => format!("{CHAIN_BASE}/grch38-chm13v2.chain"),
            (Build::Chm13v2, Build::Grch38) => format!("{CHAIN_BASE}/chm13v2-grch38.chain"),
            (Build::Grch37, Build::Chm13v2) => format!("{CHAIN_BASE}/hg19-chm13v2.chain"),
            (Build::Chm13v2, Build::Grch37) => format!("{CHAIN_BASE}/chm13v2-hg19.chain"),
            (Build::Grch38, Build::Grch37) => {
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz".to_string()
            }
            (Build::Grch37, Build::Grch38) => {
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz".to_string()
            }
            _ => return None,
        };
        Some(ChainSource { from, to, url, sha256: None })
    }

    /// The UCSC `cytoBand` table URL for a build (gzipped) — the source for genome-region
    /// metadata (centromere/telomere/cytoband). A user URL override under `references["<build>:cytoband"]`
    /// is honored. `None` for builds without a known table.
    pub fn cytoband_source(&self, build: Build) -> Option<String> {
        let default = match build.nuclear() {
            Build::Grch38 => Some("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"),
            Build::Grch37 => Some("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"),
            Build::Chm13v2 => Some("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/database/cytoBandMapped.txt.gz"),
            Build::Chm13v2MaskedRcrs => unreachable!("nuclear() collapses the masked variant"),
        };
        let key = format!("{}:cytoband", build.as_str());
        self.config
            .references
            .get(&key)
            .and_then(|o| o.url.clone())
            .or_else(|| default.map(str::to_string))
    }

    /// The annotation-mask source for a registered name (see [`Y_STRUCTURAL_MASKS`]), or `None`
    /// if unknown. A user URL override under `references[name]` is honored.
    pub fn mask_source(&self, name: &str) -> Option<MaskSource> {
        let file = Y_STRUCTURAL_MASKS.iter().find(|(n, _)| *n == name).map(|(_, f)| *f)?;
        let ov = self.config.references.get(name);
        let url = ov
            .and_then(|o| o.url.clone())
            .unwrap_or_else(|| format!("{ANNOTATION_BASE}/{file}"));
        let sha256 = ov.and_then(|o| o.sha256.clone());
        Some(MaskSource {
            name: name.to_string(),
            url,
            sha256,
        })
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
        // GRCh38↔GRCh37 lift via UCSC's gzipped over.chain (decompressed on download).
        let g38_g37 = reg.chain_source(Build::Grch38, Build::Grch37).unwrap();
        assert!(g38_g37.url.ends_with("hg38ToHg19.over.chain.gz"));
        assert!(reg
            .chain_source(Build::Grch37, Build::Grch38)
            .unwrap()
            .url
            .ends_with("hg19ToHg38.over.chain.gz"));
    }

    #[test]
    fn masked_rcrs_is_a_resolvable_cacheable_build() {
        // Aliases canonicalize to the masked variant, and its as_str round-trips (so the cache
        // file is `chm13v2.0_maskedY_rCRS.fa` — distinct from plain chm13).
        for alias in ["chm13v2.0_maskedY_rCRS", "chm13_maskedY_rcrs", "CHM13V2.0_MASKEDY_RCRS"] {
            assert_eq!(canonical_build(alias), Some(Build::Chm13v2MaskedRcrs), "alias {alias}");
        }
        assert_eq!(Build::Chm13v2MaskedRcrs.as_str(), "chm13v2.0_maskedY_rCRS");
        assert_eq!(
            canonical_build(Build::Chm13v2MaskedRcrs.as_str()),
            Some(Build::Chm13v2MaskedRcrs)
        );
        // Plain chm13 spellings still map to plain chm13.
        assert_eq!(canonical_build("chm13v2.0"), Some(Build::Chm13v2));

        let reg = Registry::new(UserConfig::default());
        assert!(reg
            .reference_source(Build::Chm13v2MaskedRcrs)
            .url
            .ends_with("chm13v2.0_maskedY_rCRS.fa.gz"));
    }

    #[test]
    fn reference_polarity_records_the_chm13_y_is_j_trap() {
        // CHM13 (and the masked variant) carry HG002's haplogroup-J Y: the reference base is
        // derived, so the metadata must flag it as not-ancestral.
        for b in [Build::Chm13v2, Build::Chm13v2MaskedRcrs] {
            let p = b.reference_polarity();
            assert!(
                p.chr_y.contains("HG002") && p.chr_y.contains('J'),
                "{}: {}",
                b.as_str(),
                p.chr_y
            );
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
        assert!(reg
            .chain_source(Build::Chm13v2MaskedRcrs, Build::Grch38)
            .unwrap()
            .url
            .ends_with("chm13v2-grch38.chain"));
    }

    #[test]
    fn y_structural_mask_sources_resolve() {
        let reg = Registry::new(UserConfig::default());
        let m = reg.mask_source("chm13v2.0Y_amplicons_v1").unwrap();
        assert_eq!(m.name, "chm13v2.0Y_amplicons_v1");
        assert!(m.url.ends_with("/annotation/chm13v2.0Y_amplicons_v1.bed"), "{}", m.url);
        assert!(reg.mask_source("chm13v2.0Y_inverted_repeats_v1").is_some());
        assert!(reg.mask_source("chm13v2.0Y_AZF_DYZ_v1").is_some());
        assert!(reg.mask_source("not_a_mask").is_none());
        // All registered masks resolve.
        assert!(Y_STRUCTURAL_MASKS.iter().all(|(n, _)| reg.mask_source(n).is_some()));
    }

    #[test]
    fn user_override_wins() {
        let mut references = HashMap::new();
        references.insert(
            "chm13v2.0".to_string(),
            BuildOverride {
                local_path: Some("/data/chm13.fa".into()),
                url: None,
                sha256: None,
                auto_download: true,
            },
        );
        let reg = Registry::new(UserConfig { references });
        assert_eq!(reg.local_override(Build::Chm13v2), Some("/data/chm13.fa"));
    }

    #[test]
    fn user_config_round_trips_with_auto_download() {
        let dir = std::env::temp_dir().join(format!("dun-refcfg-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("config").join("reference_sources.json");

        let mut references = HashMap::new();
        references.insert(
            "GRCh38".to_string(),
            BuildOverride {
                local_path: Some("/refs/grch38.fa".into()),
                url: None,
                sha256: None,
                auto_download: false,
            },
        );
        let cfg = UserConfig { references };
        cfg.save(&path).unwrap();

        let loaded = UserConfig::load(&path);
        assert_eq!(loaded, cfg);
        let ov = loaded.for_build(Build::Grch38).unwrap();
        assert_eq!(ov.local_path.as_deref(), Some("/refs/grch38.fa"));
        assert!(!ov.auto_download);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn invalid_config_falls_back_to_default_not_a_panic() {
        // A corrupt/half-written config must not crash the app or apply garbage — it warns (stderr)
        // and yields the empty default, so resolution reverts to auto-download. (The #26 failure mode,
        // now impossible to *produce* thanks to atomic writes, but load must still tolerate one.)
        let dir = std::env::temp_dir().join(format!("dun-refcfg-bad-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("reference_sources.json");
        std::fs::create_dir_all(&dir).unwrap();
        // The exact torn shape from the bug report: short head + stale tail.
        std::fs::write(&path, "{\"references\":{\"chm13v2.0_maskedY_rCRS\":{\"auto_download\":true}}}eference/hs37d5.fa\",\"auto_download\":false}}}").unwrap();
        assert_eq!(UserConfig::load(&path), UserConfig::default());
        // An absent file is the same empty default (the normal no-config case).
        assert_eq!(UserConfig::load(&dir.join("nope.json")), UserConfig::default());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

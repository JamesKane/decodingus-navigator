//! The known-build registry: canonical build names + their source URLs (reference FASTAs
//! and liftover chains), with an optional user JSON override. Defaults come from
//! `documents/chm13-reference-resources.md` (CHM13 assets on the public human-pangenomics bucket)
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
    /// The CHM13v2.0 analysis set, with the Y PAR hard-masked and rCRS in place of the
    /// mitochondrion. This is the reference we recommend for a short-read call. The mask on the
    /// PAR removes the artifacts of a read that maps to both X and Y. An rCRS chrM matches the
    /// haplotree coordinates.
    ///
    /// Its **nuclear** coordinates are the same as [`Build::Chm13v2`], see [`Build::nuclear`]. Only
    /// two things change: chrY carries N in the PAR, and chrM is a different sequence. So it reuses
    /// CHM13's liftover chains.
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

    /// The build whose **nuclear coordinate system** this one shares. A plain assembly gives
    /// itself. The masked and rCRS variant gives [`Build::Chm13v2`], because a mask writes N and
    /// does not move a coordinate.
    ///
    /// The liftover chains key on this. So the masked variant reuses CHM13's chains, and the code
    /// holds no second copy of them. Its chrM is different, because the masked chrM is rCRS. But no
    /// chain ever lifts mtDNA: that is a direct rCRS query.
    pub fn nuclear(self) -> Build {
        match self {
            Build::Chm13v2MaskedRcrs => Build::Chm13v2,
            other => other,
        }
    }

    /// The provenance of this reference's haploid sequences. It is a permanent reminder that the
    /// reference allele is a *coordinate system*, and never a source of ancestral or derived
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

/// The provenance of a reference's haploid sequences, which are chrY and chrM.
///
/// It exists to make one rule plain, and easy to find. **The ancestral and derived polarity must
/// always come from the haplotree, against the sample's own called base.** It must never come from
/// the question "is the sample's base equal to the reference (REF), or not (ALT)?"
///
/// CHM13v2.0 is the trap. Its chrY is HG002, which is a haplogroup-J Y. At many Y-SNP sites the
/// reference base is the *derived* allele. So an assumption that REF is ancestral would invert
/// those calls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferencePolarity {
    /// What the reference's chrY is, and why its allele is not a polarity source.
    pub chr_y: &'static str,
    /// What the reference's mitochondrion is, and its polarity caveat.
    pub chr_m: &'static str,
}

/// Map any common form of a build name to its canonical [`Build`]. The match ignores case.
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

/// Where the app fetches a reference FASTA from. It carries a rough size, for the download prompt,
/// and an optional pinned SHA-256 of the downloaded artifact. That hash comes from the publisher,
/// when we know it, and the code checks the download against it before it accepts the file. `None`
/// means we have no authoritative hash to pin against yet.
#[derive(Debug, Clone)]
pub struct ReferenceSource {
    pub build: Build,
    pub url: String,
    pub est_bytes: u64,
    pub sha256: Option<String>,
}

/// Where the app fetches a liftover chain from. That is a UCSC `.chain`, 1:1. It can carry a
/// pinned SHA-256.
#[derive(Debug, Clone)]
pub struct ChainSource {
    pub from: Build,
    pub to: Build,
    pub url: String,
    pub sha256: Option<String>,
}

/// Where the app fetches a named annotation-mask BED from, such as the curated CHM13 Y structural
/// regions. `name` is the cache key, and the stem of the file name. It can carry a pinned
/// SHA-256.
#[derive(Debug, Clone)]
pub struct MaskSource {
    pub name: String,
    pub url: String,
    pub sha256: Option<String>,
}

/// The curated CHM13v2.0 chrY structural-region BEDs (marbl/CHM13, Rhie et al. 2023). These are
/// the zones where a paralog is likely, and the app uses them to flag a Y call that it can not
/// trust. The key is a name that stays stable in the cache. They all sit on the human-pangenomics
/// bucket, beside the references and the chains.
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
// GRCh38: the "no-ALT analysis set" that the community uses (Heng Li,
// https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use). NCBI's HTTPS mirror
// serves it bgzipped, at about 873 MB. The Broad plain FASTA that we used before is about 3.25 GB.
// That is about 3.7 times less to download on a slow or cellular connection, and
// `decompress_and_index` decompresses and indexes the `.gz` here.
//
// It puts a `chr` prefix on a contig name, and its chrM is rCRS, which matches the app's GRCh38
// convention. It carries no ALT, no decoy, and no HLA, which is what Heng Li recommends. A CRAM
// that aligns to a Broad-specific decoy or ALT contig needs a URL override in
// `reference_sources.json`.
const GRCH38_FA: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz";
// GRCh37: hs37-1kg, which is the 1000 Genomes phase-1 reference and Heng Li's recommendation. The
// EBI 1000genomes HTTPS mirror serves it gzipped, at about 892 MB. The Broad plain FASTA is about
// 3.14 GB. Its contig names are bare, as `1`, `X`, and `MT`, and its MT is rCRS. That matches the
// app's GRCh37 convention.
const GRCH37_FA: &str = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz";
const CHAIN_BASE: &str = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo";
const ANNOTATION_BASE: &str = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation";

/// The authoritative SHA-256 of a reference FASTA's **downloaded artifact**, in lowercase hex.
/// That artifact is the `.fa.gz` or `.fasta` exactly as the server sends it. A value appears here
/// once we confirm a publisher checksum.
///
/// It is `None` until then. The integrity machinery is ready, and the pins arrive over time. Add a
/// value here, or give one for a single build in the `sha256` override of
/// `reference_sources.json`.
fn default_reference_sha(build: Build) -> Option<&'static str> {
    match build {
        // We wait for confirmed publisher checksums, for the T2T human-pangenomics bucket and the
        // Broad references.
        Build::Grch38 | Build::Grch37 | Build::Chm13v2 | Build::Chm13v2MaskedRcrs => None,
    }
}

/// The user override for one build, from `reference_sources.json`.
#[derive(Debug, Clone, Default, Deserialize, Serialize, PartialEq, Eq)]
pub struct BuildOverride {
    /// Use this local FASTA as-is (already decompressed + indexed); never download.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub local_path: Option<String>,
    /// Override the download URL.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub url: Option<String>,
    /// Pin an authoritative SHA-256 of the downloaded artifact, in lowercase hex. The code refuses
    /// a download that does not match. This lets a user give a publisher checksum that we do not
    /// ship.
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
    /// Load the config when it is there. A file that is absent or unreadable gives the empty
    /// default. An override is advisory, and never fatal: a new user with no config gets the
    /// auto-download that the app manages itself.
    ///
    /// A file that **exists but does not parse** also falls back to the defaults, but the code
    /// **warns** about it. To drop such a file in silence is how a power user's `local_path`
    /// override goes away. The app then surprises them with a full reference download.
    ///
    /// That is issue #26, where a write that was not atomic, and that ran beside another write,
    /// corrupted the config. See [`crate::cache::atomic_write`]. Say so, and do not go back to the
    /// defaults without a word.
    pub fn load(path: &Path) -> Self {
        // Use `read_atomic`, and not `fs::read_to_string`. On Windows a save that runs beside this
        // read leaves the path delete-pending for a moment. An unreadable config here means the
        // user's overrides go away with no word. That is the symptom of issue #26, by another
        // route.
        let Ok(text) = crate::cache::read_atomic(path)
            .and_then(|b| String::from_utf8(b).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e)))
        else {
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

    /// Store to `path`, and make the parent `config/` directory. The output is pretty-printed.
    ///
    /// The write is **atomic**: a temp file and a rename, see [`crate::cache::atomic_write`].
    /// Spawned worker tasks rewrite this file, and two of them can run together. A write that is
    /// not atomic then corrupts it into a new head and an old tail.
    ///
    /// A caller must still keep two read-change-write cycles apart, and should prefer one bulk
    /// save. If not, an update can go missing. Atomicity promises only that the file is never
    /// *torn*.
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
        // A hash that the user pinned wins over the built-in one. That built-in is None until we
        // confirm a publisher hash.
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

    /// The liftover chain source for a pair of builds, if the registry holds one. The code first
    /// normalizes each build to its nuclear coordinate system. So the masked and rCRS variant
    /// reuses CHM13's chains, because its nuclear coordinates are the same.
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
        Some(ChainSource {
            from,
            to,
            url,
            sha256: None,
        })
    }

    /// The UCSC `cytoBand` table URL for a build, gzipped. It is the source for the genome-region
    /// metadata: the centromere, the telomere, and the cytoband. The code obeys a user URL override
    /// under `references["<build>:cytoband"]`. It gives `None` for a build with no known table.
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

    /// The annotation-mask source for a registered name, see [`Y_STRUCTURAL_MASKS`]. It gives
    /// `None` for a name that it does not know. The code obeys a user URL override under
    /// `references[name]`.
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
        // The GRCh38↔GRCh37 lift, through UCSC's gzipped over.chain. The download decompresses
        // it.
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
        // An alias canonicalizes to the masked variant, and its as_str round-trips. So the cache
        // file is `chm13v2.0_maskedY_rCRS.fa`, which is not the name of the plain chm13 file.
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
        // CHM13, and the masked variant, carry HG002's haplogroup-J Y. There the reference base is
        // the derived allele, so the metadata must mark it as not ancestral.
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

        // So the masked variant reuses CHM13's chains. It is the same file, with normalized
        // endpoints, and there is no second entry under the masked name.
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
        // A config that is corrupt, or half written, must not stop the app, and must not apply bad
        // values. The code warns on stderr and gives the empty default, so the resolve goes back to
        // the auto-download.
        //
        // That is the failure mode of issue #26. The atomic writes make it impossible to *produce*
        // now, but a load must still accept one.
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

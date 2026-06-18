//! The reference gateway: resolves a build name → a cached, decompressed, indexed FASTA
//! (fetching on a miss), and caches liftover chains for `du-bio` to parse. Cheap to clone
//! (the app holds one). Per-key locks prevent concurrent double-downloads.

use std::collections::HashMap;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex as StdMutex};

use tokio::sync::Mutex;

use crate::error::RefgenomeError;
use crate::registry::{canonical_build, Build, Registry, UserConfig};
use crate::{cache, download, index};

/// What [`ReferenceGateway::reference_status`] reports for a build (no download performed).
#[derive(Debug, Clone)]
pub enum RefStatus {
    /// Present in the cache (path is the indexed `.fa`).
    Cached(PathBuf),
    /// A user-pinned local FASTA (config `local_path`).
    LocalOverride(PathBuf),
    /// Not cached; would fetch `url` (~`est_bytes`).
    NeedsDownload { url: String, est_bytes: u64 },
    /// Unrecognized build name.
    Unknown,
}

/// Result of [`ReferenceGateway::verify_reference`] — re-hashing a cached reference against its
/// integrity sidecar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerifyOutcome {
    /// The cached file's SHA-256 matches its recorded sidecar.
    Verified,
    /// The cached file's hash differs from the sidecar — likely on-disk corruption.
    Mismatch { expected: String, got: String },
    /// Cached, but no sidecar to check against (e.g. a user-pinned local FASTA, or pre-dates this).
    NoSidecar,
    /// Nothing cached for this build.
    NotCached,
}

/// A tree position lifted to another build: the original `tree_pos` plus its `(contig, pos)`
/// in the target build (all 1-based). `reverse` is true when the target chain is on the minus
/// strand — the caller must reverse-complement the base it reads there (large tracts of the
/// CHM13 Y are inverted relative to GRCh38).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiftedPos {
    pub tree_pos: i64,
    pub contig: String,
    pub pos: i64,
    pub reverse: bool,
}

/// Lift/drop counts from [`ReferenceGateway::lift_hipstr_bed`].
#[derive(Debug, Default, Clone, Copy)]
pub struct LiftStats {
    pub total: usize,
    pub lifted: usize,
    /// An endpoint fell in a chain gap / non-syntenic region.
    pub dropped_unmapped: usize,
    /// The two endpoints lifted to different target contigs.
    pub dropped_split: usize,
    /// The lifted span was implausible vs the source (likely a bad lift through the repeat).
    pub dropped_span: usize,
}

#[derive(Clone)]
pub struct ReferenceGateway {
    base: PathBuf,
    http: reqwest::Client,
    locks: Arc<StdMutex<HashMap<String, Arc<Mutex<()>>>>>,
    /// In-memory (layer-1) cache of parsed genome regions, shared across clones.
    regions_cache: Arc<StdMutex<HashMap<Build, Arc<crate::regions::GenomeRegions>>>>,
}

impl ReferenceGateway {
    /// Build a gateway rooted at `base` (the cache dir).
    pub fn new(base: PathBuf, http: reqwest::Client) -> Self {
        ReferenceGateway {
            base,
            http,
            locks: Arc::new(StdMutex::new(HashMap::new())),
            regions_cache: Arc::new(StdMutex::new(HashMap::new())),
        }
    }

    /// Path to the user reference-source overrides file.
    pub fn config_path(&self) -> PathBuf {
        self.base.join("config").join("reference_sources.json")
    }

    /// A registry over the **current** on-disk overrides — re-read each call so edits made via the
    /// Settings UI apply without rebuilding the gateway. Resolution is per-analysis, so the small
    /// JSON read is negligible.
    fn registry(&self) -> Registry {
        Registry::new(UserConfig::load(&self.config_path()))
    }

    fn lock_for(&self, key: &str) -> Arc<Mutex<()>> {
        let mut m = self.locks.lock().unwrap();
        m.entry(key.to_string()).or_insert_with(|| Arc::new(Mutex::new(()))).clone()
    }

    /// Cache/override status of a build — no I/O beyond a stat, never downloads.
    pub fn reference_status(&self, build_name: &str) -> RefStatus {
        let Some(build) = canonical_build(build_name) else {
            return RefStatus::Unknown;
        };
        let registry = self.registry();
        if let Some(local) = registry.local_override(build) {
            let p = PathBuf::from(local);
            if cache::is_present(&p) {
                return RefStatus::LocalOverride(p);
            }
        }
        let fa = cache::reference_path(&self.base, build);
        if cache::is_present(&fa) && cache::is_present(&cache::reference_fai(&self.base, build)) {
            return RefStatus::Cached(fa);
        }
        let src = registry.reference_source(build);
        RefStatus::NeedsDownload { url: src.url, est_bytes: src.est_bytes }
    }

    /// The cached/overridden reference path, or `None` if a download would be required.
    pub fn cached_reference(&self, build_name: &str) -> Option<PathBuf> {
        match self.reference_status(build_name) {
            RefStatus::Cached(p) | RefStatus::LocalOverride(p) => Some(p),
            _ => None,
        }
    }

    /// Resolve a build to a usable indexed `.fa`, downloading + decompressing + indexing on a
    /// miss. `progress(received, total)` is called as bytes arrive during any download.
    pub async fn resolve_reference(
        &self,
        build_name: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, RefgenomeError> {
        let build = canonical_build(build_name).ok_or_else(|| RefgenomeError::UnknownBuild(build_name.to_string()))?;
        if let Some(p) = self.cached_reference(build_name) {
            return Ok(p);
        }

        let lock = self.lock_for(build.as_str());
        let _guard = lock.lock().await;
        if let Some(p) = self.cached_reference(build_name) {
            return Ok(p); // another caller finished while we waited
        }

        let src = self.registry().reference_source(build);
        let fa = cache::reference_path(&self.base, build);
        let dl = download_target(&fa, &src.url);
        let artifact_sha = download::download(&self.http, &src.url, &dl, progress).await?;
        // Pinned (publisher) verification, on the downloaded artifact exactly as served.
        verify_pinned(&dl, src.sha256.as_deref(), &artifact_sha)?;

        let fa_out = fa.clone();
        let fa_sha = tokio::task::spawn_blocking(move || index::decompress_and_index(&dl, &fa_out))
            .await
            .map_err(|e| RefgenomeError::Message(format!("indexing join error: {e}")))??;
        // TOFU sidecar of the decompressed reference (for later offline re-verification).
        write_sidecar(&fa, &fa_sha);
        Ok(fa)
    }

    /// Resolve a liftover chain to a cached `.chain` file, downloading on a miss.
    pub async fn resolve_chain(
        &self,
        from_name: &str,
        to_name: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, RefgenomeError> {
        let (from, to) = self.chain_builds(from_name, to_name)?;
        let path = cache::chain_path(&self.base, from, to);
        if cache::is_present(&path) {
            return Ok(path);
        }
        let lock = self.lock_for(&format!("chain:{}-{}", from.as_str(), to.as_str()));
        let _guard = lock.lock().await;
        if cache::is_present(&path) {
            return Ok(path);
        }
        let src = self
            .registry()
            .chain_source(from, to)
            .ok_or_else(|| RefgenomeError::NoChain { from: from.as_str().into(), to: to.as_str().into() })?;
        let sha = download::download(&self.http, &src.url, &path, progress).await?;
        verify_pinned(&path, src.sha256.as_deref(), &sha)?; // a chain is stored as-downloaded
        write_sidecar(&path, &sha);
        Ok(path)
    }

    /// Resolve a named annotation mask (see [`registry::Y_STRUCTURAL_MASKS`]) to a cached BED,
    /// downloading on a miss. Cached under `<base>/masks/<name>.bed` — pull once, use many.
    pub async fn resolve_mask(
        &self,
        name: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, RefgenomeError> {
        let path = cache::mask_path(&self.base, name);
        if cache::is_present(&path) {
            return Ok(path);
        }
        let lock = self.lock_for(&format!("mask:{name}"));
        let _guard = lock.lock().await;
        if cache::is_present(&path) {
            return Ok(path);
        }
        let src = self
            .registry()
            .mask_source(name)
            .ok_or_else(|| RefgenomeError::Message(format!("unknown mask {name}")))?;
        let sha = download::download(&self.http, &src.url, &path, progress).await?;
        verify_pinned(&path, src.sha256.as_deref(), &sha)?; // a mask BED is stored as-downloaded
        write_sidecar(&path, &sha);
        Ok(path)
    }

    /// Whether a named annotation mask is already cached (no I/O beyond a stat).
    pub fn cached_mask(&self, name: &str) -> Option<PathBuf> {
        let path = cache::mask_path(&self.base, name);
        cache::is_present(&path).then_some(path)
    }

    /// Resolve a build's genome-region metadata (centromere/telomere/cytoband/PAR) through a
    /// 2-layer cache: in-memory (parsed) over a disk JSON (`<base>/regions/<build>.json`), refreshed
    /// from the UCSC `cytoBand` table on a miss/expiry. If the refresh fails but a (possibly stale)
    /// disk copy exists, that copy is used — region data is stable, so stale beats nothing.
    pub async fn genome_regions(
        &self,
        build_name: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<Arc<crate::regions::GenomeRegions>, RefgenomeError> {
        let build = canonical_build(build_name)
            .ok_or_else(|| RefgenomeError::UnknownBuild(build_name.to_string()))?
            .nuclear();

        // Layer 1: in-memory.
        if let Some(r) = self.regions_cache.lock().unwrap().get(&build).cloned() {
            return Ok(r);
        }

        let lock = self.lock_for(&format!("regions:{}", build.as_str()));
        let _guard = lock.lock().await;
        if let Some(r) = self.regions_cache.lock().unwrap().get(&build).cloned() {
            return Ok(r); // another caller finished while we waited
        }

        let json_path = cache::regions_path(&self.base, build);
        // Layer 2: a fresh, version-matching disk copy.
        if let Some(age) = cache::age_days(&json_path) {
            if age < REGIONS_TTL_DAYS {
                if let Some(r) = load_regions_json(&json_path) {
                    return Ok(self.memo_regions(build, r));
                }
            }
        }

        // Refresh from UCSC cytoBand; fall back to a stale disk copy if the fetch fails.
        match self.fetch_regions(build, progress).await {
            Ok(regions) => {
                let json = serde_json::to_string(&regions).map_err(|e| RefgenomeError::Message(e.to_string()))?;
                if let Some(parent) = json_path.parent() {
                    let _ = std::fs::create_dir_all(parent);
                }
                std::fs::write(&json_path, json).map_err(|e| RefgenomeError::io(&json_path, e))?;
                Ok(self.memo_regions(build, regions))
            }
            Err(e) => match load_regions_json(&json_path) {
                Some(r) => Ok(self.memo_regions(build, r)), // stale, but usable offline
                None => Err(e),
            },
        }
    }

    /// Cached genome regions without any network — in-memory, else a disk copy (any age). `None`
    /// when neither is present.
    pub fn cached_genome_regions(&self, build_name: &str) -> Option<Arc<crate::regions::GenomeRegions>> {
        let build = canonical_build(build_name)?.nuclear();
        if let Some(r) = self.regions_cache.lock().unwrap().get(&build).cloned() {
            return Some(r);
        }
        let r = load_regions_json(&cache::regions_path(&self.base, build))?;
        Some(self.memo_regions(build, r))
    }

    fn memo_regions(&self, build: Build, regions: crate::regions::GenomeRegions) -> Arc<crate::regions::GenomeRegions> {
        let arc = Arc::new(regions);
        self.regions_cache.lock().unwrap().insert(build, arc.clone());
        arc
    }

    /// Download + gunzip + parse the UCSC cytoBand table for a build.
    async fn fetch_regions(
        &self,
        build: Build,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<crate::regions::GenomeRegions, RefgenomeError> {
        let url = self
            .registry()
            .cytoband_source(build)
            .ok_or_else(|| RefgenomeError::Message(format!("no cytoBand source for {}", build.as_str())))?;
        let gz = self.base.join("regions").join(format!("{}.cytoband.txt.gz", build.as_str()));
        download::download(&self.http, &url, &gz, progress).await?;
        let text = read_gz_to_string(&gz)?;
        Ok(crate::regions::GenomeRegions::from_cytoband(build.as_str(), &text))
    }

    /// Parse the cached chain for a build pair into a `du-bio` `Liftover` (call
    /// [`resolve_chain`](Self::resolve_chain) first to ensure it's present).
    pub fn load_liftover(&self, from_name: &str, to_name: &str) -> Result<du_bio::liftover::Liftover, RefgenomeError> {
        let (from, to) = self.chain_builds(from_name, to_name)?;
        let path = cache::chain_path(&self.base, from, to);
        if !cache::is_present(&path) {
            return Err(RefgenomeError::Message(format!(
                "liftover chain {}->{} not cached; resolve it first",
                from.as_str(),
                to.as_str()
            )));
        }
        let text = std::fs::read_to_string(&path).map_err(|e| RefgenomeError::io(&path, e))?;
        du_bio::liftover::Liftover::parse(&text).map_err(|e| RefgenomeError::Message(e.to_string()))
    }

    /// Whether a liftover chain is registered for this build pair (both names canonicalize
    /// and a chain source exists). No I/O.
    pub fn chain_available(&self, from: &str, to: &str) -> bool {
        match (canonical_build(from), canonical_build(to)) {
            (Some(f), Some(t)) => self.registry().chain_source(f, t).is_some(),
            _ => false,
        }
    }

    /// Lift 1-based `positions` on `contig` from build `from` to build `to`, using the cached
    /// chain (call [`resolve_chain`](Self::resolve_chain) first). Positions in gaps /
    /// non-syntenic regions are dropped. UCSC chains are 0-based half-open while genomic
    /// positions are 1-based, so we lift `p - 1` and return `q + 1`.
    pub fn lift_positions(
        &self,
        from: &str,
        to: &str,
        contig: &str,
        positions: &[i64],
    ) -> Result<Vec<LiftedPos>, RefgenomeError> {
        let lo = self.load_liftover(from, to)?;
        // Walk chains directly (rather than Liftover::lift) so we can capture the target
        // strand, which the base-reader needs to reverse-complement inverted lifts.
        Ok(positions
            .iter()
            .filter_map(|&p| {
                lo.chains.iter().filter(|c| c.t_name == contig).find_map(|c| {
                    c.lift(p - 1).map(|q| LiftedPos {
                        tree_pos: p,
                        contig: c.q_name.clone(),
                        pos: q + 1,
                        reverse: c.q_strand == '-',
                    })
                })
            })
            .collect())
    }

    /// Lift a HipSTR-format reference BED from `from` build to `to` build via the cached chain
    /// (call [`resolve_chain`](Self::resolve_chain) first), writing a new gzipped BED in the target
    /// coordinates. Each tract's endpoints (`[start+1, end+1]`, 1-based, end-inclusive — the HipSTR
    /// convention) are lifted; a locus is kept only when both endpoints map to the **same** target
    /// contig with a plausible span (0.5×–2× the source span — guards against bad lifts through the
    /// repeat). `ref_copies` is recomputed from the lifted span (the target assembly's own repeat
    /// count); period / name / motif carry over. `only_contig` (matched `chr`-prefix-insensitively)
    /// restricts the lift, e.g. `Some("chrY")` for a Y-only reference. Returns lift/drop counts.
    pub fn lift_hipstr_bed(
        &self,
        from: &str,
        to: &str,
        in_bed_gz: &Path,
        out_bed_gz: &Path,
        only_contig: Option<&str>,
    ) -> Result<LiftStats, RefgenomeError> {
        use std::io::{BufRead, BufReader, Write};

        let lo = self.load_liftover(from, to)?;
        let strip = |s: &str| s.strip_prefix("chr").unwrap_or(s).to_ascii_uppercase();
        let want = only_contig.map(strip);

        // Lift one 1-based position on a chr-prefixed source contig → (target contig, 1-based pos).
        let lift1 = |tname: &str, p: i64| -> Option<(String, i64)> {
            lo.chains
                .iter()
                .filter(|c| c.t_name == tname)
                .find_map(|c| c.lift(p - 1).map(|q| (c.q_name.clone(), q + 1)))
        };

        let file = std::fs::File::open(in_bed_gz).map_err(|e| RefgenomeError::io(in_bed_gz, e))?;
        let rd = BufReader::new(flate2::read::MultiGzDecoder::new(file));
        let out = std::fs::File::create(out_bed_gz).map_err(|e| RefgenomeError::io(out_bed_gz, e))?;
        let mut enc = flate2::write::GzEncoder::new(out, flate2::Compression::default());
        let mut stats = LiftStats::default();

        for line in rd.lines() {
            let line = line.map_err(|e| RefgenomeError::io(in_bed_gz, e))?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let (Ok(start), Ok(end)) = (f[1].parse::<i64>(), f[2].parse::<i64>()) else { continue };
            let contig = f[0];
            if let Some(w) = &want {
                if strip(contig) != *w {
                    continue;
                }
            }
            stats.total += 1;
            let period = f[3];
            let period_n: f64 = period.parse().unwrap_or(0.0);
            let name = f.get(5).copied().unwrap_or("");
            let motif = f.get(6).copied().unwrap_or("");

            let tname = format!("chr{}", strip(contig));
            let (Some(a), Some(b)) = (lift1(&tname, start + 1), lift1(&tname, end + 1)) else {
                stats.dropped_unmapped += 1;
                continue;
            };
            if a.0 != b.0 {
                stats.dropped_split += 1; // endpoints lifted to different contigs
                continue;
            }
            let (lo_pos, hi_pos) = (a.1.min(b.1), a.1.max(b.1));
            let (src_span, dst_span) = (end - start + 1, hi_pos - lo_pos + 1);
            if dst_span <= 0 || (dst_span as f64) < 0.5 * src_span as f64 || (dst_span as f64) > 2.0 * src_span as f64 {
                stats.dropped_span += 1; // implausible span — likely a bad lift through the repeat
                continue;
            }
            let ref_copies = if period_n > 0.0 { dst_span as f64 / period_n } else { 0.0 };
            // Back to BED (0-based-inclusive [lo-1, hi-1]); bare contig, matching the HipSTR format.
            writeln!(enc, "{}\t{}\t{}\t{period}\t{ref_copies}\t{name}\t{motif}", strip(&a.0), lo_pos - 1, hi_pos - 1)
                .map_err(|e| RefgenomeError::io(out_bed_gz, e))?;
            stats.lifted += 1;
        }
        enc.finish().map_err(|e| RefgenomeError::io(out_bed_gz, e))?;
        Ok(stats)
    }

    /// Re-hash a cached reference and compare to its integrity sidecar (TOFU, written at download
    /// time). Detects on-disk corruption of the cached `.fa`. Re-reads the whole FASTA, so call it
    /// from a blocking context (it's an explicit, user-triggered check, not the hot path). A
    /// user-pinned local FASTA has no sidecar → [`VerifyOutcome::NoSidecar`].
    pub fn verify_reference(&self, build_name: &str) -> Result<VerifyOutcome, RefgenomeError> {
        let fa = match self.reference_status(build_name) {
            RefStatus::Cached(p) | RefStatus::LocalOverride(p) => p,
            _ => return Ok(VerifyOutcome::NotCached),
        };
        let Some(expected) = read_sidecar(&fa) else { return Ok(VerifyOutcome::NoSidecar) };
        let got = index::hash_file(&fa)?;
        Ok(if expected.eq_ignore_ascii_case(&got) {
            VerifyOutcome::Verified
        } else {
            VerifyOutcome::Mismatch { expected, got }
        })
    }

    /// Resolve a `(from, to)` build-name pair for **chain** purposes, normalized to nuclear
    /// coordinates — so the masked+rCRS variant resolves to (and reuses the cache of) CHM13's
    /// chains rather than a duplicate keyed by its own name.
    fn chain_builds(&self, from: &str, to: &str) -> Result<(Build, Build), RefgenomeError> {
        let f = canonical_build(from).ok_or_else(|| RefgenomeError::UnknownBuild(from.to_string()))?;
        let t = canonical_build(to).ok_or_else(|| RefgenomeError::UnknownBuild(to.to_string()))?;
        Ok((f.nuclear(), t.nuclear()))
    }
}

/// Genome-region cache freshness window. Region metadata (cytoband/centromere) is effectively
/// static per assembly, so a long TTL avoids needless refetches; a stale copy is still used if a
/// refetch fails.
const REGIONS_TTL_DAYS: f64 = 90.0;

/// Load + deserialize a cached genome-regions JSON, dropping it if it predates the current schema
/// version (so a parser/overlay change invalidates stale caches).
fn load_regions_json(path: &Path) -> Option<crate::regions::GenomeRegions> {
    let text = std::fs::read_to_string(path).ok()?;
    let regions: crate::regions::GenomeRegions = serde_json::from_str(&text).ok()?;
    (regions.version == crate::regions::REGIONS_VERSION).then_some(regions)
}

/// gunzip a downloaded `.gz` into a UTF-8 string (cytoBand tables are small).
fn read_gz_to_string(path: &Path) -> Result<String, RefgenomeError> {
    use std::io::Read;
    let file = std::fs::File::open(path).map_err(|e| RefgenomeError::io(path, e))?;
    let mut dec = flate2::read::MultiGzDecoder::new(std::io::BufReader::new(file));
    let mut s = String::new();
    dec.read_to_string(&mut s).map_err(|e| RefgenomeError::io(path, e))?;
    Ok(s)
}

/// The `<file>.sha256` integrity-sidecar path for a cached artifact.
fn sidecar_path(file: &Path) -> PathBuf {
    let mut s: OsString = file.as_os_str().to_os_string();
    s.push(".sha256");
    PathBuf::from(s)
}

/// Write the TOFU integrity sidecar (best-effort: a missing sidecar just means "unverifiable",
/// never fatal).
fn write_sidecar(file: &Path, sha_hex: &str) {
    let _ = std::fs::write(sidecar_path(file), sha_hex);
}

/// Read the recorded sidecar digest for a cached file, if present (first whitespace-delimited token).
fn read_sidecar(file: &Path) -> Option<String> {
    let s = std::fs::read_to_string(sidecar_path(file)).ok()?;
    s.split_whitespace().next().map(str::to_string)
}

/// Compare a freshly-downloaded artifact's digest to a pinned (publisher) hash, if one is set.
/// On mismatch the partial download at `path` is removed and an [`RefgenomeError::Integrity`] is
/// returned; `None` pinned hash = nothing to check (TOFU only).
fn verify_pinned(path: &Path, pinned: Option<&str>, got: &str) -> Result<(), RefgenomeError> {
    if let Some(expected) = pinned {
        if !expected.eq_ignore_ascii_case(got) {
            let _ = std::fs::remove_file(path);
            return Err(RefgenomeError::Integrity {
                path: path.to_path_buf(),
                expected: expected.to_string(),
                got: got.to_string(),
            });
        }
    }
    Ok(())
}

/// Where to stream a download before decompression: `<fa>.gz` for gzipped sources, else a
/// neutral `<fa>.dl` (decompress_and_index renames a non-gzip download into place).
fn download_target(fa: &Path, url: &str) -> PathBuf {
    let suffix = if url.ends_with(".gz") { "gz" } else { "dl" };
    let mut s: OsString = fa.as_os_str().to_os_string();
    s.push(".");
    s.push(suffix);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-gw-{}-{}", std::process::id(), tag));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    fn gw(base: &Path) -> ReferenceGateway {
        ReferenceGateway::new(base.to_path_buf(), reqwest::Client::new())
    }

    #[test]
    fn status_reports_cache_state_without_network() {
        let base = scratch("status");
        let g = gw(&base);
        // Unknown build.
        assert!(matches!(g.reference_status("nope"), RefStatus::Unknown));
        // Missing → needs download.
        assert!(matches!(g.reference_status("chm13v2.0"), RefStatus::NeedsDownload { .. }));
        assert!(g.cached_reference("chm13v2.0").is_none());

        // Seed a cached reference (.fa + .fai).
        let refs = base.join("references");
        std::fs::create_dir_all(&refs).unwrap();
        std::fs::write(refs.join("chm13v2.0.fa"), b">x\nACGT\n").unwrap();
        std::fs::write(refs.join("chm13v2.0.fa.fai"), b"x\t4\t3\t4\t5\n").unwrap();
        match g.reference_status("CHM13") {
            RefStatus::Cached(p) => assert!(p.ends_with("chm13v2.0.fa")),
            other => panic!("expected Cached, got {other:?}"),
        }
        assert!(g.cached_reference("chm13v2.0").is_some());
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn masked_rcrs_reuses_the_chm13_chain_cache() {
        let base = scratch("maskedchain");
        let dir = base.join("liftover");
        std::fs::create_dir_all(&dir).unwrap();
        // Only the plain-CHM13 chain file exists on disk.
        std::fs::write(
            dir.join("GRCh38-to-chm13v2.0.chain"),
            "chain 1000 chrZ 1000 + 0 100 chrZp 1000 + 0 100 1\n100\n",
        )
        .unwrap();
        let g = gw(&base);
        // A chain is "available" for the masked variant, and loading it reuses the CHM13 file
        // (normalized to nuclear coords) rather than a missing masked-named one.
        assert!(g.chain_available("GRCh38", "chm13v2.0_maskedY_rCRS"));
        let lo = g.load_liftover("GRCh38", "chm13v2.0_maskedY_rCRS").unwrap();
        assert_eq!(lo.lift("chrZ", 50), Some(("chrZp".to_string(), 50)));
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn load_liftover_parses_a_cached_chain() {
        let base = scratch("chain");
        let dir = base.join("liftover");
        std::fs::create_dir_all(&dir).unwrap();
        // A minimal UCSC chain (du-bio format): chrZ -> chrZp, one 100bp block.
        std::fs::write(
            dir.join("GRCh38-to-chm13v2.0.chain"),
            "chain 1000 chrZ 1000 + 0 100 chrZp 1000 + 0 100 1\n100\n",
        )
        .unwrap();
        let g = gw(&base);
        let lo = g.load_liftover("GRCh38", "chm13v2.0").unwrap();
        assert_eq!(lo.lift("chrZ", 50), Some(("chrZp".to_string(), 50)));
        // Not-yet-resolved pair errors clearly.
        assert!(g.load_liftover("chm13v2.0", "GRCh38").is_err());
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn lift_positions_is_one_based_in_and_out() {
        let base = scratch("liftpos");
        let dir = base.join("liftover");
        std::fs::create_dir_all(&dir).unwrap();
        // chrY t[0,100) -> chrY q[0,100) (identity over the first 100 bp).
        std::fs::write(dir.join("GRCh38-to-chm13v2.0.chain"), "chain 1 chrY 1000 + 0 100 chrY 1000 + 0 100 1\n100\n").unwrap();
        let g = gw(&base);

        // 1-based 101 -> 0-based 100 -> outside the block -> dropped.
        // 1-based 50 -> 0-based 49 -> q 49 -> 1-based 50; 1-based 100 -> 0-based 99 -> 1-based 100.
        let lifted = g.lift_positions("GRCh38", "chm13v2.0", "chrY", &[50, 100, 101]).unwrap();
        assert_eq!(
            lifted,
            vec![
                LiftedPos { tree_pos: 50, contig: "chrY".into(), pos: 50, reverse: false },
                LiftedPos { tree_pos: 100, contig: "chrY".into(), pos: 100, reverse: false },
            ]
        );
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn lift_positions_flags_minus_strand_targets() {
        let base = scratch("liftrev");
        let dir = base.join("liftover");
        std::fs::create_dir_all(&dir).unwrap();
        // chrY t[0,10) -> chrY q on the MINUS strand (q_size 100): pos 0 -> 100-1-0 = 99.
        std::fs::write(dir.join("GRCh38-to-chm13v2.0.chain"), "chain 1 chrY 1000 + 0 10 chrY 100 - 0 10 1\n10\n").unwrap();
        let g = gw(&base);
        // 1-based tree 1 -> 0-based 0 -> q 99 -> 1-based 100, flagged reverse.
        let lifted = g.lift_positions("GRCh38", "chm13v2.0", "chrY", &[1]).unwrap();
        assert_eq!(lifted, vec![LiftedPos { tree_pos: 1, contig: "chrY".into(), pos: 100, reverse: true }]);
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn genome_regions_load_from_disk_then_memory_and_reject_stale_version() {
        let base = scratch("regions");
        let dir = base.join("regions");
        std::fs::create_dir_all(&dir).unwrap();
        // Seed a parsed regions JSON (no network). Build key normalizes to chm13v2.0.
        let regions = crate::regions::GenomeRegions::from_cytoband(
            "chm13v2.0",
            "chrY\t0\t300000\tp11.32\tgneg\nchrY\t300000\t62460029\tq11\tgpos50\n",
        );
        std::fs::write(
            cache::regions_path(&base, Build::Chm13v2),
            serde_json::to_string(&regions).unwrap(),
        )
        .unwrap();

        let g = gw(&base);
        // Disk hit (any alias / the masked variant share CHM13's regions).
        let r = g.cached_genome_regions("hs1").expect("disk-cached regions");
        assert!(r.chromosome("chrY").unwrap().par.len() == 2); // PAR overlaid by the parser
        // Second call is an in-memory hit (same Arc).
        let r2 = g.cached_genome_regions("chm13v2.0_maskedY_rCRS").unwrap();
        assert!(Arc::ptr_eq(&r, &r2));

        // A wrong-version JSON is rejected (forces a refetch path), so cold-cache load is None.
        let g2 = gw(&base);
        std::fs::write(cache::regions_path(&base, Build::Chm13v2), r#"{"build":"chm13v2.0","version":"OLD","chromosomes":{}}"#).unwrap();
        assert!(g2.cached_genome_regions("chm13v2.0").is_none());
        let _ = std::fs::remove_dir_all(&base);
    }
}

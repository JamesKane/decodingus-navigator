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

#[derive(Clone)]
pub struct ReferenceGateway {
    base: PathBuf,
    http: reqwest::Client,
    registry: Registry,
    locks: Arc<StdMutex<HashMap<String, Arc<Mutex<()>>>>>,
}

impl ReferenceGateway {
    /// Build a gateway rooted at `base` (the cache dir), loading any user source overrides.
    pub fn new(base: PathBuf, http: reqwest::Client) -> Self {
        let config = UserConfig::load(&base.join("config").join("reference_sources.json"));
        ReferenceGateway {
            base,
            http,
            registry: Registry::new(config),
            locks: Arc::new(StdMutex::new(HashMap::new())),
        }
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
        if let Some(local) = self.registry.local_override(build) {
            let p = PathBuf::from(local);
            if cache::is_present(&p) {
                return RefStatus::LocalOverride(p);
            }
        }
        let fa = cache::reference_path(&self.base, build);
        if cache::is_present(&fa) && cache::is_present(&cache::reference_fai(&self.base, build)) {
            return RefStatus::Cached(fa);
        }
        let src = self.registry.reference_source(build);
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

        let src = self.registry.reference_source(build);
        let fa = cache::reference_path(&self.base, build);
        let dl = download_target(&fa, &src.url);
        download::download(&self.http, &src.url, &dl, progress).await?;

        let fa_out = fa.clone();
        tokio::task::spawn_blocking(move || index::decompress_and_index(&dl, &fa_out))
            .await
            .map_err(|e| RefgenomeError::Message(format!("indexing join error: {e}")))??;
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
            .registry
            .chain_source(from, to)
            .ok_or_else(|| RefgenomeError::NoChain { from: from.as_str().into(), to: to.as_str().into() })?;
        download::download(&self.http, &src.url, &path, progress).await?;
        Ok(path)
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

    fn chain_builds(&self, from: &str, to: &str) -> Result<(Build, Build), RefgenomeError> {
        let f = canonical_build(from).ok_or_else(|| RefgenomeError::UnknownBuild(from.to_string()))?;
        let t = canonical_build(to).ok_or_else(|| RefgenomeError::UnknownBuild(to.to_string()))?;
        Ok((f, t))
    }
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
}

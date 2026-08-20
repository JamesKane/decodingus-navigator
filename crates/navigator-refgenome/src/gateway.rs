//! The reference gateway. It resolves a build name to a cached FASTA, in plain text, with an
//! index, and it fetches that file on a cache miss. It also caches liftover chains, for `du-bio`
//! to parse.
//!
//! It is low-cost to clone, and the app holds one. A lock on each key stops two downloads of the
//! same file at the same time.

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

/// The result of [`ReferenceGateway::verify_reference`], which hashes a cached reference again and
/// compares it against the integrity sidecar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VerifyOutcome {
    /// The cached file's SHA-256 matches its recorded sidecar.
    Verified,
    /// The cached file's hash is not the one in the sidecar. The file on disk is probably
    /// corrupt.
    Mismatch { expected: String, got: String },
    /// Cached, but no sidecar to check against (e.g. a user-pinned local FASTA, or pre-dates this).
    NoSidecar,
    /// Nothing cached for this build.
    NotCached,
}

/// A tree position, lifted to another build. It holds the original `tree_pos`, and the
/// `(contig, pos)` in the target build. All of these are 1-based.
///
/// `reverse` is true when the target chain is on the minus strand. The caller must then
/// reverse-complement the base that it reads there. Large tracts of the CHM13 Y run in the
/// opposite direction from GRCh38.
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

    /// A registry over the **current** overrides on disk. Every call reads them again, so an edit
    /// in the Settings UI applies with no rebuild of the gateway. One resolve happens for each
    /// analysis, so the small JSON read costs nothing that matters.
    fn registry(&self) -> Registry {
        Registry::new(UserConfig::load(&self.config_path()))
    }

    fn lock_for(&self, key: &str) -> Arc<Mutex<()>> {
        let mut m = self.locks.lock().unwrap();
        m.entry(key.to_string())
            .or_insert_with(|| Arc::new(Mutex::new(())))
            .clone()
    }

    /// The cache and override status of a build. It does no I/O beyond a stat, and never
    /// downloads.
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
        RefStatus::NeedsDownload {
            url: src.url,
            est_bytes: src.est_bytes,
        }
    }

    /// The reference path from the cache or the override. It gives `None` when the app would have
    /// to download the file.
    pub fn cached_reference(&self, build_name: &str) -> Option<PathBuf> {
        match self.reference_status(build_name) {
            RefStatus::Cached(p) | RefStatus::LocalOverride(p) => Some(p),
            _ => None,
        }
    }

    /// Resolve a build to an indexed `.fa` that the app can use. On a cache miss it downloads the
    /// file, decompresses it, and indexes it. It calls `progress(received, total)` as the bytes
    /// arrive during any download.
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

    /// Resolve a liftover chain to a cached `.chain` file. It downloads on a cache miss.
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
            .ok_or_else(|| RefgenomeError::NoChain {
                from: from.as_str().into(),
                to: to.as_str().into(),
            })?;
        let sha = download::download(&self.http, &src.url, &path, progress).await?;
        verify_pinned(&path, src.sha256.as_deref(), &sha)?; // verify the artifact exactly as served

        // The cache stores a chain as plain text, and `load_liftover` reads it with
        // `read_to_string`. Every chain takes the same path here. A downloaded artifact in gzip
        // form decompresses in place. UCSC serves `.over.chain.gz`, and the curated bucket serves a
        // plain `.chain`. The first two bytes of the file say which one it is, so the source URL
        // and its extension do not matter.
        maybe_gunzip_in_place(&path)?;
        write_sidecar(&path, &sha);
        Ok(path)
    }

    /// Resolve a named annotation mask (see [`registry::Y_STRUCTURAL_MASKS`]) to a cached BED. It
    /// downloads on a cache miss. The cache holds it at `<base>/masks/<name>.bed`, so the app
    /// fetches it once and reads it many times.
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

    /// Resolve a published **ancestry/IBD asset** to a cached file under `<base>/ancestry/<name>`.
    /// Such an asset is a prebuilt panel, a PCA, a manifest, or a genetic map. It downloads from
    /// `url` on a cache miss, under a lock for that name.
    ///
    /// There is **no pinned sha** here, unlike a reference or a mask. The app checks the asset
    /// against the asset manifest that it downloaded. This is how a user gets the prebuilt panels,
    /// and nobody has to run the offline `panelbuild` tool.
    pub async fn resolve_ancestry_asset(
        &self,
        name: &str,
        url: &str,
        progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
    ) -> Result<PathBuf, RefgenomeError> {
        let path = self.base.join("ancestry").join(name);
        if cache::is_present(&path) {
            return Ok(path);
        }
        let lock = self.lock_for(&format!("asset:{name}"));
        let _guard = lock.lock().await;
        if cache::is_present(&path) {
            return Ok(path); // another caller finished while we waited
        }
        download::download(&self.http, url, &path, progress).await?; // creates parent, retries, streams
        Ok(path)
    }

    /// Whether a named annotation mask is already cached (no I/O beyond a stat).
    pub fn cached_mask(&self, name: &str) -> Option<PathBuf> {
        let path = cache::mask_path(&self.base, name);
        cache::is_present(&path).then_some(path)
    }

    /// Resolve a build's genome-region metadata through a cache of two layers. That metadata is
    /// the centromere, the telomere, the cytoband, and the PAR.
    ///
    /// The first layer is a parsed copy in memory, over a JSON on disk at
    /// `<base>/regions/<build>.json`. On a miss, or when the copy expires, it refreshes from the
    /// UCSC `cytoBand` table.
    ///
    /// If that refresh fails, and a disk copy exists, the code uses that copy, even a stale one.
    /// Region data does not change, so a stale copy is better than none.
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

    /// Cached genome regions, with no network. It takes the copy in memory, else a copy from disk,
    /// of any age. It gives `None` when it has neither.
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
        let gz = self
            .base
            .join("regions")
            .join(format!("{}.cytoband.txt.gz", build.as_str()));
        download::download(&self.http, &url, &gz, progress).await?;
        let text = read_gz_to_string(&gz)?;
        Ok(crate::regions::GenomeRegions::from_cytoband(build.as_str(), &text))
    }

    /// Parse the cached chain for a pair of builds into a `du-bio` `Liftover`. Call
    /// [`resolve_chain`](Self::resolve_chain) first, to make sure the chain is there.
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

    /// Whether the registry holds a liftover chain for this pair of builds. Both names must
    /// canonicalize, and a chain source must exist. It does no I/O.
    pub fn chain_available(&self, from: &str, to: &str) -> bool {
        match (canonical_build(from), canonical_build(to)) {
            (Some(f), Some(t)) => self.registry().chain_source(f, t).is_some(),
            _ => false,
        }
    }

    /// Lift the 1-based `positions` on `contig` from build `from` to build `to`, with the cached
    /// chain. Call [`resolve_chain`](Self::resolve_chain) first. The code drops a position that
    /// falls in a gap, or in a region with no synteny.
    ///
    /// A UCSC chain is 0-based half-open, and a genomic position is 1-based. So the code lifts
    /// `p - 1` and returns `q + 1`.
    pub fn lift_positions(
        &self,
        from: &str,
        to: &str,
        contig: &str,
        positions: &[i64],
    ) -> Result<Vec<LiftedPos>, RefgenomeError> {
        let lo = self.load_liftover(from, to)?;
        // Walk the chains directly, and do not call Liftover::lift. The walk gives the target
        // strand, and the base-reader needs that to reverse-complement an inverted lift.
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

    /// Lift the `[start, end)` intervals on `contig` from build `from` to build `to`, through the
    /// cached chain. Call [`resolve_chain`](Self::resolve_chain) first.
    ///
    /// An interval gets through only when **both endpoints** land on the same target contig, and
    /// the lifted span is 0.5 to 2 times the source span.
    ///
    /// [`lift_hipstr_bed`](Self::lift_hipstr_bed) applies the same check, for the same reason. A
    /// chain can map the two ends of a repeat to places far apart. A structural mask that stretched
    /// across a chromosome, with no warning, is worse than no mask at all, because it would
    /// suppress a real call.
    ///
    /// The code drops an interval that fails, and counts it. The caller gets `(intervals,
    /// dropped)`, so it can say how much of the mask got through. Without that count it would show
    /// a thin mask as the whole one.
    pub fn lift_intervals(
        &self,
        from: &str,
        to: &str,
        contig: &str,
        intervals: &[(i64, i64)],
    ) -> Result<(Vec<(i64, i64)>, usize), RefgenomeError> {
        let lo = self.load_liftover(from, to)?;
        let lift1 = |p: i64| -> Option<(String, i64)> {
            lo.chains
                .iter()
                .filter(|c| c.t_name == contig)
                .find_map(|c| c.lift(p).map(|q| (c.q_name.clone(), q)))
        };
        let (mut out, mut dropped) = (Vec::with_capacity(intervals.len()), 0usize);
        for &(s, e) in intervals {
            if e <= s {
                dropped += 1;
                continue;
            }
            // The endpoints first. That is the exact case, and the only one that can prove the
            // full span.
            if let (Some(a), Some(b)) = (lift1(s), lift1(e)) {
                if a.0 == b.0 {
                    // An inverted lift comes back reversed; the interval is the same stretch.
                    let (lo_p, hi_p) = if a.1 <= b.1 { (a.1, b.1) } else { (b.1, a.1) };
                    let (src, dst) = ((e - s) as f64, (hi_p - lo_p) as f64);
                    if dst >= src * 0.5 && dst <= src * 2.0 {
                        out.push((lo_p, hi_p));
                        continue;
                    }
                }
            }
            // Fall back to the interior. An interval whose *ends* sit in a gap can still have a
            // body that maps. That is the usual case for an amplicon, because an amplicon is
            // exactly where the two assemblies disagree.
            //
            // To recover the part that maps matters, because these masks *suppress* a call. A mask
            // that is too small lets in hundreds of false novel variants. A mask that is too large
            // costs a few true calls, in sequence that is a known paralog.
            const SAMPLES: i64 = 64;
            let step = ((e - s) / SAMPLES).max(1);
            let mut by_contig: std::collections::HashMap<String, (i64, i64, usize)> = Default::default();
            let mut tried = 0usize;
            let mut p = s;
            while p <= e {
                tried += 1;
                if let Some((c, q)) = lift1(p) {
                    let ent = by_contig.entry(c).or_insert((q, q, 0));
                    ent.0 = ent.0.min(q);
                    ent.1 = ent.1.max(q);
                    ent.2 += 1;
                }
                p += step;
            }
            // The main target contig, and only when enough of the interval mapped. One stray point
            // is not evidence of a region.
            let best = by_contig.into_values().max_by_key(|v| v.2);
            match best {
                Some((lo_p, hi_p, hits)) if hits * 4 >= tried && hi_p > lo_p => {
                    // No lower span bound here: a partial recovery is *expected* to be shorter than
                    // the source. The upper bound still guards against a lift smeared across the
                    // chromosome.
                    if (hi_p - lo_p) as f64 <= (e - s) as f64 * 2.0 {
                        out.push((lo_p, hi_p));
                    } else {
                        dropped += 1;
                    }
                }
                _ => dropped += 1,
            }
        }
        Ok((out, dropped))
    }

    /// Lift a HipSTR-format reference BED from build `from` to build `to`, through the cached
    /// chain, and write a new gzipped BED in the target coordinates. Call
    /// [`resolve_chain`](Self::resolve_chain) first.
    ///
    /// It lifts the endpoints of each tract, which are `[start+1, end+1]`, 1-based and
    /// end-inclusive, as the HipSTR convention says. It keeps a locus only when both endpoints map
    /// to the **same** target contig, with a plausible span of 0.5 to 2 times the source span. That
    /// span check guards against a bad lift through the repeat.
    ///
    /// It computes `ref_copies` again from the lifted span, which gives the target assembly's own
    /// repeat count. The period, the name, and the motif carry over.
    ///
    /// `only_contig` limits the lift, and the match ignores a `chr` prefix. An example is
    /// `Some("chrY")`, for a reference that holds Y alone. It returns the count of the loci it
    /// lifted and the count it dropped.
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
        let strip = navigator_domain::contig::bare_upper;
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
            let (Ok(start), Ok(end)) = (f[1].parse::<i64>(), f[2].parse::<i64>()) else {
                continue;
            };
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
            let ref_copies = if period_n > 0.0 {
                dst_span as f64 / period_n
            } else {
                0.0
            };
            // Back to BED (0-based-inclusive [lo-1, hi-1]); bare contig, matching the HipSTR format.
            writeln!(
                enc,
                "{}\t{}\t{}\t{period}\t{ref_copies}\t{name}\t{motif}",
                strip(&a.0),
                lo_pos - 1,
                hi_pos - 1
            )
            .map_err(|e| RefgenomeError::io(out_bed_gz, e))?;
            stats.lifted += 1;
        }
        enc.finish().map_err(|e| RefgenomeError::io(out_bed_gz, e))?;
        Ok(stats)
    }

    /// Hash a cached reference again, and compare the result to its integrity sidecar. The
    /// download writes that sidecar, on a trust-on-first-use basis. This finds corruption of the
    /// cached `.fa` on disk.
    ///
    /// It reads the whole FASTA again, so call it from a context that can block. It is an explicit
    /// check that the user starts, and not part of the hot path. A local FASTA that the user pinned
    /// has no sidecar, and gives [`VerifyOutcome::NoSidecar`].
    pub fn verify_reference(&self, build_name: &str) -> Result<VerifyOutcome, RefgenomeError> {
        let fa = match self.reference_status(build_name) {
            RefStatus::Cached(p) | RefStatus::LocalOverride(p) => p,
            _ => return Ok(VerifyOutcome::NotCached),
        };
        let Some(expected) = read_sidecar(&fa) else {
            return Ok(VerifyOutcome::NoSidecar);
        };
        let got = index::hash_file(&fa)?;
        Ok(if expected.eq_ignore_ascii_case(&got) {
            VerifyOutcome::Verified
        } else {
            VerifyOutcome::Mismatch { expected, got }
        })
    }

    /// Resolve a `(from, to)` pair of build names for a **chain**, normalized to nuclear
    /// coordinates. So the masked and rCRS variant resolves to CHM13's chains, and reuses that
    /// cache. There is no second copy under its own name.
    fn chain_builds(&self, from: &str, to: &str) -> Result<(Build, Build), RefgenomeError> {
        let f = canonical_build(from).ok_or_else(|| RefgenomeError::UnknownBuild(from.to_string()))?;
        let t = canonical_build(to).ok_or_else(|| RefgenomeError::UnknownBuild(to.to_string()))?;
        Ok((f.nuclear(), t.nuclear()))
    }
}

/// How long a genome-region cache entry stays fresh. The region metadata, which is the cytoband
/// and the centromere, does not change inside one assembly. So a long TTL keeps the app off the
/// network. A stale copy still serves when a fetch fails.
const REGIONS_TTL_DAYS: f64 = 90.0;

/// Load a cached genome-regions JSON and deserialize it. It drops a file from before the current
/// schema version. So a change to the parser, or to an overlay, invalidates a stale cache.
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

/// If `path` holds gzip-compressed data, decompress it to plain text in place. The first two
/// bytes, `1f 8b`, say whether it does. In all other cases, leave the file as it is. So the cache
/// can hold every chain as plain text, whether its source served a `.chain` or a `.chain.gz`.
fn maybe_gunzip_in_place(path: &Path) -> Result<(), RefgenomeError> {
    let bytes = std::fs::read(path).map_err(|e| RefgenomeError::io(path, e))?;
    if bytes.len() < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b {
        return Ok(()); // not gzipped — already plain text
    }
    let text = read_gz_to_string(path)?;
    std::fs::write(path, text).map_err(|e| RefgenomeError::io(path, e))
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

/// Compare the digest of an artifact that the app has just downloaded to a pinned hash from the
/// publisher, when there is one. On a mismatch it removes the partial download at `path` and gives
/// an [`RefgenomeError::Integrity`]. A pinned hash of `None` means there is nothing to check, and
/// the app trusts the file on first use.
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
        assert!(matches!(
            g.reference_status("chm13v2.0"),
            RefStatus::NeedsDownload { .. }
        ));
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
        // A chain is "available" for the masked variant. A load of it reuses the CHM13 file,
        // normalized to nuclear coordinates. It does not look for a file under the masked name,
        // which does not exist.
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
        // A small UCSC chain, in the du-bio format: chrZ to chrZp, with one 100bp block.
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
        std::fs::write(
            dir.join("GRCh38-to-chm13v2.0.chain"),
            "chain 1 chrY 1000 + 0 100 chrY 1000 + 0 100 1\n100\n",
        )
        .unwrap();
        let g = gw(&base);

        // 1-based 101 -> 0-based 100 -> outside the block -> dropped.
        // 1-based 50 -> 0-based 49 -> q 49 -> 1-based 50; 1-based 100 -> 0-based 99 -> 1-based 100.
        let lifted = g
            .lift_positions("GRCh38", "chm13v2.0", "chrY", &[50, 100, 101])
            .unwrap();
        assert_eq!(
            lifted,
            vec![
                LiftedPos {
                    tree_pos: 50,
                    contig: "chrY".into(),
                    pos: 50,
                    reverse: false
                },
                LiftedPos {
                    tree_pos: 100,
                    contig: "chrY".into(),
                    pos: 100,
                    reverse: false
                },
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
        std::fs::write(
            dir.join("GRCh38-to-chm13v2.0.chain"),
            "chain 1 chrY 1000 + 0 10 chrY 100 - 0 10 1\n10\n",
        )
        .unwrap();
        let g = gw(&base);
        // 1-based tree 1 -> 0-based 0 -> q 99 -> 1-based 100, flagged reverse.
        let lifted = g.lift_positions("GRCh38", "chm13v2.0", "chrY", &[1]).unwrap();
        assert_eq!(
            lifted,
            vec![LiftedPos {
                tree_pos: 1,
                contig: "chrY".into(),
                pos: 100,
                reverse: true
            }]
        );
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

        // The code refuses a JSON of the wrong version, which forces a fetch. So a load from a
        // cold cache gives None.
        let g2 = gw(&base);
        std::fs::write(
            cache::regions_path(&base, Build::Chm13v2),
            r#"{"build":"chm13v2.0","version":"OLD","chromosomes":{}}"#,
        )
        .unwrap();
        assert!(g2.cached_genome_regions("chm13v2.0").is_none());
        let _ = std::fs::remove_dir_all(&base);
    }
}

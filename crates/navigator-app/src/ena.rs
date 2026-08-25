//! This module gets the files of a Grid work unit from ENA.
//!
//! A node receives a full manifest with its lease. The manifest holds the URL, the md5 value and
//! the byte size of each file. The AppView makes that list. So a node does not ask ENA to find
//! anything, and this module only gets the files that the manifest names. See
//! `documents/design/distributed-compute-grid.md` §7.1.
//!
//! # Why this module is not `refgenome::download`
//!
//! §7.1 first told us to copy that function. Three properties make a copy impossible. Each one is
//! more important here than for a reference genome.
//!
//! 1. **That function can not continue a transfer.** It sends no `Range` header. An interrupted
//!    transfer must start again at zero.
//!
//!    A reference genome is about 900 MB. An ENA run file is 10 to 30 GB on the connection of a
//!    volunteer. The ability to continue decides if a unit ever completes.
//! 2. **That function calculates SHA-256.** ENA publishes md5. A checksum that you can not compare
//!    gives no integrity.
//! 3. **Its one retry is blind.** It does the whole transfer again after each error. Some errors
//!    stay after a second try.
//!
//! This module does copy the `.part` file and the rename at the end. A file that is not complete
//! must never look like a complete file. The atomic rename makes that sure.
//!
//! # How a transfer continues, and what the hash must survive
//!
//! When a transfer continues, an earlier process calculated the hash of the bytes on the disk. That
//! process is gone. So the module reads its own `.part` prefix again and builds the md5 state
//! again. Then it asks for the remainder.
//!
//! That costs one sequential read of the bytes that the disk already holds. It is much cheaper than
//! a second download of those bytes. It also occurs while the module waits for the server.
//!
//! The other method is to calculate the hash of the full file at the end. That method costs the
//! same read. But it can not start before the transfer stops. The method here also finds a bad
//! prefix. And the usual case, with no interruption, costs nothing.

use crate::error::AppError;
use md5::{Digest, Md5};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use tokio::io::{AsyncReadExt, AsyncWriteExt};

/// How much effort to give one file. A test can supply its own policy. Then the test can examine
/// the failure path and does not wait for the full delay schedule of the application.
/// [`RetryPolicy::default`] is the policy that the application uses.
#[derive(Debug, Clone, Copy)]
pub struct RetryPolicy {
    /// How many tries the module makes before it stops. Each try continues at the point where the
    /// last try stopped. So this value limits *stalls* and not the total transfer time. The module
    /// does not try again while a file continues to make progress.
    pub attempts: u32,
    /// If the module waits between two tries. `false` removes the delay. A test of the failure
    /// path does not need the delay.
    pub backoff: bool,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        RetryPolicy {
            attempts: 5,
            backoff: true,
        }
    }
}

/// The delay before try *n*, in seconds: 2, 4, 8, 16. The delay has a limit. A node of a volunteer
/// must give a unit back. It must not hold a lease in a retry loop that has no end.
fn backoff_secs(attempt: u32) -> u64 {
    1u64 << attempt.min(4)
}

/// How much disk space a unit with an alignment needs, as a factor on the size of the manifest.
/// The file arrives ready to read, so the node adds only its own analysis output.
const SPACE_MULTIPLE_ALIGNED: u64 = 3;

/// How much disk space a unit with reads needs, as a factor on the size of the manifest.
///
/// The factor is much larger here. The manifest names **compressed** reads, and the node then
/// writes three files that hold the same data in a different form. `mapped.bam` comes from the
/// reads. `sorted.bam` exists while `mapped.bam` is still on the disk, and the sort also spills to
/// the disk. For 30 GB of compressed reads, the peak is far above 90 GB.
///
/// A value that is too small gives the exact failure that this check prevents. The disk fills in
/// the middle of a unit, after hours of work.
const SPACE_MULTIPLE_READS: u64 = 10;

/// One file in a work unit's manifest, exactly as the AppView curated it.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ManifestFile {
    #[serde(default)]
    pub run_accession: String,
    pub url: String,
    #[serde(default)]
    pub index_url: Option<String>,
    #[serde(default)]
    pub md5: Option<String>,
    #[serde(default)]
    pub bytes: Option<i64>,
    #[serde(default)]
    pub format: String,
}

impl ManifestFile {
    /// The file name this entry lands under, taken from the URL's last segment.
    pub fn file_name(&self) -> &str {
        self.url.rsplit('/').next().unwrap_or(&self.url)
    }
}

/// ENA gives a location with no scheme, such as `ftp.sra.ebi.ac.uk/vol1/...`. Get it with HTTPS.
///
/// Do not use FTP. It is more difficult to continue an FTP transfer. A firewall on the network of a
/// volunteer frequently stops FTP. And ENA gives the same paths with HTTPS. This function does not
/// change a URL that already has a scheme. So a manifest can point to a different host.
pub fn to_https(url: &str) -> String {
    let u = url.trim();
    if u.starts_with("http://") || u.starts_with("https://") {
        u.to_string()
    } else {
        format!("https://{}", u.trim_start_matches("ftp://"))
    }
}

fn part_path(dest: &Path) -> PathBuf {
    let mut s = dest.as_os_str().to_os_string();
    s.push(".part");
    PathBuf::from(s)
}

fn hex(digest: &[u8]) -> String {
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

/// Whether a downloaded file's checksum is acceptable.
///
/// An entry with **no** md5 value is acceptable. ENA does not always publish one. If the module
/// refused such work, it would refuse most of the catalogue. The comparison ignores the letter
/// case. Hex letters have the same value in each case, and archives do not use one case only.
pub fn checksum_ok(expected: Option<&str>, actual: &str) -> bool {
    match expected.map(str::trim).filter(|s| !s.is_empty()) {
        Some(want) => want.eq_ignore_ascii_case(actual),
        None => true,
    }
}

/// Refuse a unit that is too large for the disk. The check occurs before the first byte arrives.
///
/// A full disk in the middle of a run is the worst result. The unit fails. The node holds the lease
/// until the lease ends. And the owner of the machine must remove the files.
///
/// `free_space` gives zero when it can not measure the disk. A zero lets the try continue. A
/// refusal after a failed measurement is worse than a write that fails.
pub fn preflight_space(dir: &Path, manifest: &[ManifestFile]) -> Result<(), AppError> {
    let total: u64 = manifest.iter().filter_map(|f| f.bytes).map(|b| b.max(0) as u64).sum();
    // A unit of reads needs much more room than a unit with an alignment. See the two constants.
    let reads = manifest.iter().any(|f| f.format == "FASTQ");
    let multiple = if reads {
        SPACE_MULTIPLE_READS
    } else {
        SPACE_MULTIPLE_ALIGNED
    };
    let needed = total.saturating_mul(multiple);
    let free = crate::realign_job::free_space(dir);
    if !crate::realign_job::has_room(needed, free) {
        return Err(AppError::Import(format!(
            "not enough room for this work unit: about {} GB is needed and {} GB is free on {}",
            needed / 1_000_000_000,
            free / 1_000_000_000,
            dir.display()
        )));
    }
    Ok(())
}

/// What one try must ask the server for. The bytes on the disk decide this.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Resume {
    /// The disk holds nothing that the module can use. Get the full file.
    FromStart,
    /// Continue from this offset.
    From(u64),
    /// The `.part` file already has the expected size. Check it, and get no more bytes.
    AlreadyComplete,
}

/// Decide how to continue. The `.part` size and the expected total give the answer.
///
/// A `.part` file that is **larger** than the expected total starts again at zero. The module does
/// not cut it to the correct length. Such a file is evidence that the disk holds a different file
/// from the one in the manifest. The cause can be a new revision of the file, or two files with the
/// same name.
///
/// A cut to the correct length gives a file with the correct size but the wrong md5 value. The
/// module finds that only after a second full download.
pub fn resume_from(part_len: u64, expected: Option<u64>) -> Resume {
    match expected {
        Some(total) if part_len == total && total > 0 => Resume::AlreadyComplete,
        Some(total) if part_len > total => Resume::FromStart,
        _ if part_len == 0 => Resume::FromStart,
        _ => Resume::From(part_len),
    }
}

/// Rebuild md5 state over an existing `.part` prefix.
async fn hash_prefix(path: &Path) -> Result<(Md5, u64), AppError> {
    let mut file = tokio::fs::File::open(path).await.map_err(|e| io_err(path, e))?;
    let mut hasher = Md5::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut total = 0u64;
    loop {
        let n = file.read(&mut buf).await.map_err(|e| io_err(path, e))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
        total += n as u64;
    }
    Ok((hasher, total))
}

fn io_err(path: &Path, e: std::io::Error) -> AppError {
    AppError::Import(format!("{}: {e}", path.display()))
}

/// Get one file of the manifest into `dir`. This continues an interrupted transfer, and it checks
/// the md5 value.
///
/// `progress` receives `(received, total)` as the bytes arrive. The `received` count includes the
/// prefix from the earlier try. A progress bar that starts again at zero would tell the user the
/// opposite of the truth.
pub async fn fetch_file(
    client: &reqwest::Client,
    dir: &Path,
    entry: &ManifestFile,
    cancel: &navigator_analysis::CancelToken,
    progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
) -> Result<PathBuf, AppError> {
    fetch_file_with(client, dir, entry, RetryPolicy::default(), cancel, progress).await
}

/// [`fetch_file`] with an explicit retry policy.
pub async fn fetch_file_with(
    client: &reqwest::Client,
    dir: &Path,
    entry: &ManifestFile,
    retry: RetryPolicy,
    cancel: &navigator_analysis::CancelToken,
    progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
) -> Result<PathBuf, AppError> {
    let dest = dir.join(entry.file_name());
    if dest.exists() {
        return Ok(dest); // a completed file is never re-fetched; the rename is what makes it final
    }
    tokio::fs::create_dir_all(dir).await.map_err(|e| io_err(dir, e))?;
    let part = part_path(&dest);
    let url = to_https(&entry.url);
    let expected = entry.bytes.filter(|b| *b > 0).map(|b| b as u64);

    let mut last: Option<AppError> = None;
    for attempt in 0..retry.attempts {
        if cancel.is_cancelled() {
            return Err(AppError::Import("cancelled".into()));
        }
        if attempt > 0 && retry.backoff {
            tokio::time::sleep(std::time::Duration::from_secs(backoff_secs(attempt))).await;
        }
        match fetch_once(client, &url, &part, expected, cancel, progress).await {
            Ok(actual) => {
                if !checksum_ok(entry.md5.as_deref(), &actual) {
                    // The module does not try again from the bytes on the disk. Those bytes are
                    // wrong. A transfer that continues from them gives the same wrong result. So
                    // the next try must start at zero.
                    let _ = tokio::fs::remove_file(&part).await;
                    last = Some(AppError::Import(format!(
                        "checksum mismatch for {}: expected {}, got {actual}",
                        entry.file_name(),
                        entry.md5.as_deref().unwrap_or("?")
                    )));
                    continue;
                }
                tokio::fs::rename(&part, &dest).await.map_err(|e| io_err(&dest, e))?;
                return Ok(dest);
            }
            Err(e) => last = Some(e),
        }
    }
    Err(last.unwrap_or_else(|| AppError::Import(format!("could not fetch {}", entry.file_name()))))
}

/// One try. Returns the md5 value of the complete file, in lowercase hex.
async fn fetch_once(
    client: &reqwest::Client,
    url: &str,
    part: &Path,
    expected: Option<u64>,
    cancel: &navigator_analysis::CancelToken,
    progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
) -> Result<String, AppError> {
    let part_len = tokio::fs::metadata(part).await.map(|m| m.len()).unwrap_or(0);
    let plan = resume_from(part_len, expected);

    let (mut hasher, mut received) = match plan {
        Resume::FromStart => (Md5::new(), 0),
        Resume::From(_) | Resume::AlreadyComplete => hash_prefix(part).await?,
    };
    if plan == Resume::AlreadyComplete {
        return Ok(hex(&hasher.finalize()));
    }

    let mut req = client.get(url);
    if let Resume::From(offset) = plan {
        req = req.header(reqwest::header::RANGE, format!("bytes={offset}-"));
    }
    let resp = req
        .send()
        .await
        .map_err(|e| AppError::Import(format!("{url}: {e}")))?
        .error_for_status()
        .map_err(|e| AppError::Import(format!("{url}: {e}")))?;

    // A server that ignores `Range` answers 200 and sends the full file. Accept that answer. Do
    // not add those bytes to the bytes on the disk. Such a file has the correct size only by
    // accident, and it fails the checksum.
    let restart = matches!(plan, Resume::From(_)) && resp.status() != reqwest::StatusCode::PARTIAL_CONTENT;
    if restart {
        hasher = Md5::new();
        received = 0;
    }
    let total = expected.or_else(|| resp.content_length().map(|c| c + received));

    let mut file = if received > 0 && !restart {
        tokio::fs::OpenOptions::new()
            .append(true)
            .open(part)
            .await
            .map_err(|e| io_err(part, e))?
    } else {
        tokio::fs::File::create(part).await.map_err(|e| io_err(part, e))?
    };

    let mut resp = resp;
    while let Some(chunk) = resp
        .chunk()
        .await
        .map_err(|e| AppError::Import(format!("{url}: {e}")))?
    {
        if cancel.is_cancelled() {
            // Keep the `.part` file. A later try needs those bytes. If the module removed the
            // file, a cancel would discard all the work that the user already paid for.
            file.flush().await.map_err(|e| io_err(part, e))?;
            return Err(AppError::Import("cancelled".into()));
        }
        file.write_all(&chunk).await.map_err(|e| io_err(part, e))?;
        hasher.update(&chunk);
        received += chunk.len() as u64;
        progress(received, total);
    }
    file.flush().await.map_err(|e| io_err(part, e))?;
    Ok(hex(&hasher.finalize()))
}

/// Get each file of a unit manifest into `dir`, one file after the other.
///
/// The order is sequential by design. The connection of the volunteer is the limit. So parallel
/// transfers do not finish earlier. They also increase the peak disk use, and they put more load on
/// a public archive that helps us at no cost.
pub async fn fetch_unit(
    client: &reqwest::Client,
    dir: &Path,
    manifest: &[ManifestFile],
    cancel: &navigator_analysis::CancelToken,
    progress: &mut (dyn FnMut(&str, u64, Option<u64>) + Send),
) -> Result<Vec<PathBuf>, AppError> {
    preflight_space(dir, manifest)?;
    let mut out = Vec::with_capacity(manifest.len());
    for entry in manifest {
        let name = entry.file_name().to_string();
        let mut per_file = |recv: u64, total: Option<u64>| progress(&name, recv, total);
        out.push(fetch_file(client, dir, entry, cancel, &mut per_file).await?);

        // Get the index file beside the alignment, when ENA has one.
        //
        // The index is some MB, and the alignment is 10 to 30 GB. Without the index, the node
        // reads the whole alignment one more time to make its own. So this small transfer removes
        // a full pass over the largest file of the unit.
        //
        // A failure here is not a failure of the unit. The node makes the index itself, which
        // costs time and gives the same result.
        if let Some(index_url) = entry.index_url.clone().filter(|u| !u.trim().is_empty()) {
            let sidecar = ManifestFile {
                run_accession: entry.run_accession.clone(),
                url: index_url,
                index_url: None,
                // ENA publishes no checksum for the index file, so there is nothing to compare.
                md5: None,
                bytes: None,
                format: "INDEX".to_string(),
            };
            let name = sidecar.file_name().to_string();
            let mut per_file = |recv: u64, total: Option<u64>| progress(&name, recv, total);
            let _ = fetch_file(client, dir, &sidecar, cancel, &mut per_file).await;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ena_paths_become_https_and_explicit_schemes_are_left_alone() {
        assert_eq!(
            to_https("ftp.sra.ebi.ac.uk/vol1/run/ERR/x.cram"),
            "https://ftp.sra.ebi.ac.uk/vol1/run/ERR/x.cram"
        );
        assert_eq!(
            to_https("ftp://ftp.sra.ebi.ac.uk/vol1/x.cram"),
            "https://ftp.sra.ebi.ac.uk/vol1/x.cram"
        );
        assert_eq!(to_https("https://example.org/x.cram"), "https://example.org/x.cram");
        assert_eq!(to_https("http://example.org/x.cram"), "http://example.org/x.cram");
    }

    #[test]
    fn the_file_name_comes_from_the_last_url_segment() {
        let f = ManifestFile {
            run_accession: "ERR1".into(),
            url: "ftp.sra.ebi.ac.uk/vol1/run/ERR1/sample.cram".into(),
            index_url: None,
            md5: None,
            bytes: None,
            format: "CRAM".into(),
        };
        assert_eq!(f.file_name(), "sample.cram");
    }

    #[test]
    fn resume_continues_from_what_is_already_there() {
        assert_eq!(resume_from(0, Some(100)), Resume::FromStart);
        assert_eq!(resume_from(40, Some(100)), Resume::From(40));
        assert_eq!(resume_from(100, Some(100)), Resume::AlreadyComplete);
    }

    /// A `.part` file larger than the manifest total is evidence that the disk holds a different
    /// file. A cut to the correct length gives the correct size and the wrong md5 value. The module
    /// finds that only after a second full download.
    #[test]
    fn an_oversized_part_starts_over_rather_than_being_trimmed() {
        assert_eq!(resume_from(140, Some(100)), Resume::FromStart);
    }

    /// With no expected size, the module has no value to compare. So it continues the transfer of
    /// a partial file. The checksum gives the final answer.
    #[test]
    fn an_unknown_total_still_resumes() {
        assert_eq!(resume_from(40, None), Resume::From(40));
        assert_eq!(resume_from(0, None), Resume::FromStart);
    }

    #[test]
    fn a_missing_checksum_is_not_a_failure() {
        assert!(checksum_ok(None, "d41d8cd98f00b204e9800998ecf8427e"));
        assert!(checksum_ok(Some(""), "d41d8cd98f00b204e9800998ecf8427e"));
    }

    #[test]
    fn checksums_compare_without_regard_to_hex_case() {
        assert!(checksum_ok(
            Some("D41D8CD98F00B204E9800998ECF8427E"),
            "d41d8cd98f00b204e9800998ecf8427e"
        ));
        assert!(!checksum_ok(
            Some("d41d8cd98f00b204e9800998ecf8427e"),
            "0bad0bad0bad0bad0bad0bad0bad0bad"
        ));
    }

    #[test]
    fn backoff_grows_and_then_stops_growing() {
        assert_eq!((1..=5).map(backoff_secs).collect::<Vec<_>>(), vec![2, 4, 8, 16, 16]);
    }

    /// The manifest comes from the curation query of the AppView. A FASTQ entry has no
    /// `index_url`, because `jsonb_strip_nulls` removes that key.
    #[test]
    fn a_curated_manifest_entry_decodes() {
        let json = r#"{"run_accession":"ERR2000001","url":"ftp.sra.ebi.ac.uk/vol1/s1.cram",
                       "index_url":"ftp.sra.ebi.ac.uk/vol1/s1.cram.crai",
                       "md5":"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","bytes":12000000000,"format":"CRAM"}"#;
        let f: ManifestFile = serde_json::from_str(json).unwrap();
        assert_eq!(f.file_name(), "s1.cram");
        assert_eq!(f.bytes, Some(12_000_000_000));

        let stripped = r#"{"run_accession":"ERR2","url":"ftp/r_1.fastq.gz","md5":"b","bytes":9,"format":"FASTQ"}"#;
        let f: ManifestFile = serde_json::from_str(stripped).unwrap();
        assert!(f.index_url.is_none(), "an absent sidecar must not fail to decode");
    }
}

//! Streaming HTTP download to a `.part` file with an atomic rename on completion, a coarse
//! progress callback, and one retry on transient failure. Async (reqwest); the caller runs
//! it on the app's tokio runtime.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};
use tokio::io::AsyncWriteExt;

use crate::error::RefgenomeError;

fn part_path(dest: &Path) -> PathBuf {
    let mut s: OsString = dest.as_os_str().to_os_string();
    s.push(".part");
    PathBuf::from(s)
}

/// Download `url` to `dest`, reporting `(received, total)` as bytes arrive (`total` is the
/// `Content-Length`, if the server sent one). Streams to `dest.part` and renames on success.
/// Retries once on a transient error. Returns the SHA-256 (lowercase hex) of the **downloaded
/// bytes** — computed on the fly so callers can verify against a pinned hash without a re-read.
pub async fn download(
    client: &reqwest::Client,
    url: &str,
    dest: &Path,
    progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
) -> Result<String, RefgenomeError> {
    match download_once(client, url, dest, progress).await {
        Ok(digest) => Ok(digest),
        Err(_) => download_once(client, url, dest, progress).await, // single retry
    }
}

async fn download_once(
    client: &reqwest::Client,
    url: &str,
    dest: &Path,
    progress: &mut (dyn FnMut(u64, Option<u64>) + Send),
) -> Result<String, RefgenomeError> {
    let http_err = |source: reqwest::Error| RefgenomeError::Http { url: url.to_string(), source };

    if let Some(parent) = dest.parent() {
        tokio::fs::create_dir_all(parent).await.map_err(|e| RefgenomeError::io(parent, e))?;
    }
    let part = part_path(dest);

    let mut resp = client.get(url).send().await.map_err(http_err)?.error_for_status().map_err(http_err)?;
    let total = resp.content_length();

    let mut file = tokio::fs::File::create(&part).await.map_err(|e| RefgenomeError::io(&part, e))?;
    let mut received = 0u64;
    let mut hasher = Sha256::new();
    while let Some(chunk) = resp.chunk().await.map_err(http_err)? {
        file.write_all(&chunk).await.map_err(|e| RefgenomeError::io(&part, e))?;
        hasher.update(&chunk);
        received += chunk.len() as u64;
        progress(received, total);
    }
    file.flush().await.map_err(|e| RefgenomeError::io(&part, e))?;
    drop(file);

    tokio::fs::rename(&part, dest).await.map_err(|e| RefgenomeError::io(dest, e))?;
    Ok(hasher.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

//! Decompress a downloaded reference to a plain FASTA and build its `.fai` — all in-Rust
//! (`flate2` for gzip/bgzip, `noodles::fasta` for indexing), so no samtools/GATK is needed.
//! Blocking; call from `spawn_blocking`.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use noodles::fasta;

use crate::error::RefgenomeError;

/// Append `.<ext>` to a path (e.g. `foo.fa` → `foo.fa.part`).
fn with_suffix(path: &Path, ext: &str) -> PathBuf {
    let mut s: OsString = path.as_os_str().to_os_string();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

/// Whether `path` starts with the gzip magic bytes (`1f 8b`). bgzip is gzip-compatible, so
/// `MultiGzDecoder` handles both.
fn is_gzip(path: &Path) -> Result<bool, RefgenomeError> {
    let mut f = File::open(path).map_err(|e| RefgenomeError::io(path, e))?;
    let mut magic = [0u8; 2];
    let n = read_up_to(&mut f, &mut magic).map_err(|e| RefgenomeError::io(path, e))?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
}

/// Read until the buffer is full or EOF; returns bytes read (a short final read isn't EOF).
fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

/// Turn `src` (gzip/bgzip or plain FASTA) into a plain FASTA at `fa_out` and write
/// `fa_out.fai`. On a gzip input, `src` is decompressed then deleted; on a plain input it is
/// moved into place. The `.fai` is built by indexing the plain FASTA.
pub fn decompress_and_index(src: &Path, fa_out: &Path) -> Result<(), RefgenomeError> {
    if let Some(parent) = fa_out.parent() {
        std::fs::create_dir_all(parent).map_err(|e| RefgenomeError::io(parent, e))?;
    }

    if is_gzip(src)? {
        let part = with_suffix(fa_out, "part");
        let input = File::open(src).map_err(|e| RefgenomeError::io(src, e))?;
        let mut dec = MultiGzDecoder::new(BufReader::new(input));
        let mut out = BufWriter::new(File::create(&part).map_err(|e| RefgenomeError::io(&part, e))?);
        io::copy(&mut dec, &mut out).map_err(|e| RefgenomeError::io(&part, e))?;
        out.into_inner().map_err(|e| RefgenomeError::io(&part, e.into_error()))?.sync_all().ok();
        std::fs::rename(&part, fa_out).map_err(|e| RefgenomeError::io(fa_out, e))?;
        let _ = std::fs::remove_file(src);
    } else if src != fa_out {
        std::fs::rename(src, fa_out).map_err(|e| RefgenomeError::io(fa_out, e))?;
    }

    let index = fasta::fs::index(fa_out).map_err(|e| RefgenomeError::io(fa_out, e))?;
    let fai = with_suffix(fa_out, "fai");
    let mut writer = fasta::fai::Writer::new(BufWriter::new(File::create(&fai).map_err(|e| RefgenomeError::io(&fai, e))?));
    writer.write_index(&index).map_err(|e| RefgenomeError::io(&fai, e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn scratch(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-refidx-{}-{}", std::process::id(), tag));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn decompresses_gz_and_builds_fai() {
        let dir = scratch("gz");
        let fasta = ">chr1\nACGTACGTAC\nACGT\n>chr2\nTTTTNNNNAA\n";
        // bgzip-compatible plain gzip via flate2.
        let gz = dir.join("ref.fa.gz");
        {
            let mut enc = flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::default());
            enc.write_all(fasta.as_bytes()).unwrap();
            enc.finish().unwrap();
        }
        let fa = dir.join("chm13v2.0.fa");
        decompress_and_index(&gz, &fa).unwrap();

        // Plain FASTA materialized, .gz removed, .fai built and matches a fresh index.
        assert_eq!(std::fs::read_to_string(&fa).unwrap(), fasta);
        assert!(!gz.exists());
        let fai = with_suffix(&fa, "fai");
        assert!(fai.exists());
        let expected = fasta::fs::index(&fa).unwrap();
        assert_eq!(expected.as_ref().len(), 2); // two contigs
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn plain_fasta_is_indexed_in_place() {
        let dir = scratch("plain");
        let fa = dir.join("GRCh38.fa");
        std::fs::write(&fa, ">x\nACGT\n").unwrap();
        decompress_and_index(&fa, &fa).unwrap();
        assert!(with_suffix(&fa, "fai").exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

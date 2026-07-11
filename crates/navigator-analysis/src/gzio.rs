//! Transparent decompression for line-oriented genomics text (VCF, BED, CompleteGenomics
//! masterVar).
//!
//! BGZF — the bgzip format used for `.vcf.gz` / `.bed.gz` — is a *concatenation* of
//! independent gzip members (each ≤64 KiB of payload). `flate2::read::GzDecoder` decodes
//! only the first member and then reports EOF, silently truncating any multi-block
//! bgzipped file; `MultiGzDecoder` decodes every member, so it reads plain gzip and BGZF
//! whole. bzip2 (the `.tsv.bz2` CompleteGenomics ships) is handled the same way via
//! `MultiBzDecoder`, which spans concatenated streams (pbzip2 output). Detection is by the
//! leading magic bytes, not the extension, so a compressed file is handled even when it is
//! misnamed (e.g. a `.vcf` that is really bgzf).

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Detected on-disk compression, by leading magic bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Compression {
    None,
    Gzip,
    Bzip2,
}

/// Open `path` for buffered line reading, transparently decoding gzip/BGZF when the file
/// begins with the gzip magic bytes. Plain (uncompressed) text is read directly.
///
/// See [`open_maybe_compressed`] to additionally decode bzip2.
pub fn open_maybe_gz(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    match detect_compression(&mut file)? {
        Compression::Gzip => Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))),
        // Callers of this legacy entry point only expect gzip/plain; treat bzip2 as plain
        // (they never pass a `.bz2`). Use `open_maybe_compressed` for bzip2 support.
        Compression::Bzip2 | Compression::None => Ok(Box::new(BufReader::new(file))),
    }
}

/// Open `path` for buffered line reading, transparently decoding gzip/BGZF **or** bzip2 by
/// content. Plain text is read directly. Used by importers that accept the compressed dumps
/// vendors ship (e.g. a CompleteGenomics `var-*-ASM.tsv.bz2`).
pub fn open_maybe_compressed(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    match detect_compression(&mut file)? {
        Compression::Gzip => Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))),
        Compression::Bzip2 => Ok(Box::new(BufReader::new(bzip2::read::MultiBzDecoder::new(file)))),
        Compression::None => Ok(Box::new(BufReader::new(file))),
    }
}

/// Peek the leading magic bytes for gzip (`1f 8b`) or bzip2 (`BZh`), then rewind to the start
/// so the returned reader sees the whole file.
fn detect_compression(file: &mut File) -> io::Result<Compression> {
    let mut magic = [0u8; 3];
    let n = read_up_to(file, &mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if n >= 2 && magic[..2] == [0x1f, 0x8b] {
        Ok(Compression::Gzip)
    } else if n >= 3 && &magic == b"BZh" {
        Ok(Compression::Bzip2)
    } else {
        Ok(Compression::None)
    }
}

/// Fill `buf` from `file`, tolerating short reads; returns the number of bytes read (may be
/// fewer than `buf.len()` only at EOF).
fn read_up_to(file: &mut File, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match file.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(filled)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp_dir() -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-gzio-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn gzip_member(bytes: &[u8]) -> Vec<u8> {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    #[test]
    fn reads_plain_text() {
        let path = tmp_dir().join("plain.txt");
        std::fs::write(&path, b"line-1\nline-2\n").unwrap();
        let reader = open_maybe_gz(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(Result::unwrap).collect();
        assert_eq!(lines, vec!["line-1", "line-2"]);
    }

    #[test]
    fn reads_single_member_gzip() {
        let path = tmp_dir().join("single.txt.gz");
        std::fs::write(&path, gzip_member(b"a\nb\nc\n")).unwrap();
        let reader = open_maybe_gz(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(Result::unwrap).collect();
        assert_eq!(lines, vec!["a", "b", "c"]);
    }

    #[test]
    fn reads_bzip2_by_content() {
        // A `.txt`-named bzip2 stream must be decoded by content (magic `BZh`), not extension.
        let path = tmp_dir().join("cg.txt");
        let mut enc = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
        enc.write_all(b"x\ny\nz\n").unwrap();
        let blob = enc.finish().unwrap();
        std::fs::write(&path, blob).unwrap();
        let reader = open_maybe_compressed(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(Result::unwrap).collect();
        assert_eq!(lines, vec!["x", "y", "z"]);
    }

    #[test]
    fn reads_concatenated_members_without_truncation() {
        // BGZF is a concatenation of independent gzip members. A single-member decoder
        // (GzDecoder) would stop after the first and drop the rest; MultiGzDecoder must
        // read all of them. Named `.txt` (not `.gz`) to prove detection is by content.
        let path = tmp_dir().join("multi.txt");
        let mut blob = gzip_member(b"first\n");
        blob.extend(gzip_member(b"second\n"));
        blob.extend(gzip_member(b"third\n"));
        std::fs::write(&path, blob).unwrap();
        let reader = open_maybe_gz(&path).unwrap();
        let lines: Vec<String> = reader.lines().map(Result::unwrap).collect();
        assert_eq!(lines, vec!["first", "second", "third"]);
    }
}

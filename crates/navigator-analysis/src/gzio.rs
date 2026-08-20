//! Transparent decompression for line-oriented genomics text (VCF, BED, CompleteGenomics
//! masterVar).
//!
//! BGZF is the bgzip format of a `.vcf.gz` or a `.bed.gz`. It is a *chain* of independent gzip
//! members, and each one holds 64 KiB of payload or less.
//!
//! `flate2::read::GzDecoder` decodes the first member alone, and it then reports an end of file.
//! So it cuts short any bgzipped file of more than one block, and nobody sees it happen.
//! `MultiGzDecoder` decodes every member, so it reads a plain gzip file and a BGZF file whole.
//!
//! bzip2 gets the same treatment, through `MultiBzDecoder`, which reads a chain of streams that
//! pbzip2 wrote. CompleteGenomics ships a `.tsv.bz2`.
//!
//! The code finds the compression from the first bytes of the file, and not from the extension. So
//! it handles a compressed file even when its name is wrong, such as a `.vcf` that holds bgzf.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// The compression that the file on disk uses, which the code finds from the first bytes of that
/// file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Compression {
    None,
    Gzip,
    Bzip2,
}

/// Open `path` to read its lines through a buffer. It decodes gzip and BGZF when the file starts
/// with the first bytes of a gzip stream, and the caller sees no difference. It reads plain text
/// directly.
///
/// See [`open_maybe_compressed`] to decode bzip2 as well.
pub fn open_maybe_gz(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    match detect_compression(&mut file)? {
        Compression::Gzip => Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))),
        // Callers of this legacy entry point only expect gzip/plain; treat bzip2 as plain
        // (they never pass a `.bz2`). Use `open_maybe_compressed` for bzip2 support.
        Compression::Bzip2 | Compression::None => Ok(Box::new(BufReader::new(file))),
    }
}

/// Open `path` to read its lines through a buffer. It decodes gzip and BGZF, **or** bzip2, from
/// the content of the file, and the caller sees no difference. It reads plain text directly. An
/// importer that accepts a compressed dump from a vendor uses this, such as a CompleteGenomics
/// `var-*-ASM.tsv.bz2`.
pub fn open_maybe_compressed(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    match detect_compression(&mut file)? {
        Compression::Gzip => Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))),
        Compression::Bzip2 => Ok(Box::new(BufReader::new(bzip2::read::MultiBzDecoder::new(file)))),
        Compression::None => Ok(Box::new(BufReader::new(file))),
    }
}

/// Look at the first bytes of the file, for gzip (`1f 8b`) or for bzip2 (`BZh`). Then go back to
/// the start, so that the reader that comes out sees the whole file.
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

/// Fill `buf` from `file`. It accepts a short read. It returns the count of bytes that it read,
/// and that count is below `buf.len()` only at the end of the file.
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
        // A bzip2 stream whose name ends in `.txt` must decode from its content, where the first
        // bytes are `BZh`. The extension must not decide.
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

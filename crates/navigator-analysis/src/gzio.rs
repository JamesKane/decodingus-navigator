//! Transparent decompression for line-oriented genomics text (VCF, BED).
//!
//! BGZF — the bgzip format used for `.vcf.gz` / `.bed.gz` — is a *concatenation* of
//! independent gzip members (each ≤64 KiB of payload). `flate2::read::GzDecoder` decodes
//! only the first member and then reports EOF, silently truncating any multi-block
//! bgzipped file; `MultiGzDecoder` decodes every member, so it reads plain gzip and BGZF
//! whole. Detection is by the two-byte gzip magic (`1f 8b`), not the extension, so a
//! bgzipped file is handled even when it is misnamed (e.g. a `.vcf` that is really bgzf).

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Open `path` for buffered line reading, transparently decoding gzip/BGZF when the file
/// begins with the gzip magic bytes. Plain (uncompressed) text is read directly.
pub fn open_maybe_gz(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    if is_gzip(&mut file)? {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Peek the first two bytes for the gzip magic (`1f 8b`), then rewind to the start so the
/// returned reader sees the whole file.
fn is_gzip(file: &mut File) -> io::Result<bool> {
    let mut magic = [0u8; 2];
    let n = read_up_to(file, &mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
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

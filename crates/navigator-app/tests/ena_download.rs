//! Tests of the ENA download module against a true HTTP server.
//!
//! The module exists because `refgenome::download` can not continue an interrupted transfer. A
//! mock server can not test that ability with enough care. The important behaviour is the
//! answer of the code to a `Range` request, to a `206` status and to a `200` status. It is also the
//! md5 value across a prefix that this process did not get.
//!
//! A stub that answers as the caller expects agrees with the code at each point where both are
//! wrong. So these tests use a temporary HTTP/1.1 server on the loopback address. It is about forty
//! lines, it adds no dependency, and it can give a wrong answer on purpose.

use navigator_analysis::CancelToken;
use navigator_app::ena::{self, ManifestFile, RetryPolicy};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::TcpListener;

/// A scratch directory with a unique name. It removes itself.
///
/// The workspace has no `tempfile` dependency. The usual method here is a fixed name under
/// `std::env::temp_dir()`. That is the cause of the flake in `a_present_file_resolves` when two
/// test runs occur together. A unique name for each test costs nothing and prevents that fault.
struct Scratch(std::path::PathBuf);

impl Scratch {
    fn new(tag: &str) -> Self {
        static N: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let unique = format!(
            "navigator-ena-{tag}-{}-{}",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        );
        let path = std::env::temp_dir().join(unique);
        let _ = std::fs::remove_dir_all(&path);
        std::fs::create_dir_all(&path).expect("scratch dir");
        Scratch(path)
    }
    fn path(&self) -> &std::path::Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// How the test server should behave.
#[derive(Clone, Copy, PartialEq)]
enum Mode {
    /// Honour `Range` properly: `206` plus the requested tail.
    Honest,
    /// Ignore `Range` and always send the full body with `200`. Many true servers do this.
    IgnoresRange,
    /// Serve a body that does not match the advertised md5.
    Corrupt,
}

/// Serve `body` until told to stop. Returns the bound address.
async fn serve(body: Vec<u8>, mode: Mode, stop: Arc<AtomicBool>) -> String {
    let listener = TcpListener::bind("127.0.0.1:0").await.expect("bind");
    let addr = listener.local_addr().expect("addr").to_string();
    tokio::spawn(async move {
        while !stop.load(Ordering::Relaxed) {
            let Ok((mut sock, _)) = listener.accept().await else {
                break;
            };
            let body = body.clone();
            tokio::spawn(async move {
                let mut buf = vec![0u8; 4096];
                let n = sock.read(&mut buf).await.unwrap_or(0);
                // Change to lowercase before the match. `reqwest` sends header *names* in
                // lowercase. So a server that looks for "Range:" never finds one, and gives no
                // message about it. Because of that fault, this stub tested the ignore-range
                // path when its purpose was to test a transfer that continues.
                let req = String::from_utf8_lossy(&buf[..n]).to_lowercase();

                // `range: bytes=N-`
                let start = req
                    .lines()
                    .find_map(|l| l.strip_prefix("range: bytes="))
                    .and_then(|r| r.trim().trim_end_matches('-').parse::<usize>().ok())
                    .unwrap_or(0);

                let send_partial = mode == Mode::Honest && start > 0 && start < body.len();
                let payload = if send_partial { &body[start..] } else { &body[..] };
                let status = if send_partial {
                    "HTTP/1.1 206 Partial Content"
                } else {
                    "HTTP/1.1 200 OK"
                };
                let head = format!(
                    "{status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    payload.len()
                );
                let _ = sock.write_all(head.as_bytes()).await;
                let _ = sock.write_all(payload).await;
                let _ = sock.flush().await;
            });
        }
    });
    addr
}

fn md5_hex(bytes: &[u8]) -> String {
    use md5::{Digest, Md5};
    let mut h = Md5::new();
    h.update(bytes);
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}

fn entry(addr: &str, name: &str, md5: Option<String>, bytes: Option<i64>) -> ManifestFile {
    ManifestFile {
        run_accession: "ERR000001".into(),
        url: format!("http://{addr}/{name}"),
        index_url: None,
        md5,
        bytes,
        format: "CRAM".into(),
        instrument: None,
    }
}

fn body() -> Vec<u8> {
    // The body is large enough that a transfer continues across some chunks, not inside one.
    (0..200_000u32).flat_map(|i| i.to_le_bytes()).collect()
}

#[tokio::test]
async fn a_whole_file_downloads_and_verifies() {
    let data = body();
    let stop = Arc::new(AtomicBool::new(false));
    let addr = serve(data.clone(), Mode::Honest, stop.clone()).await;
    let dir = Scratch::new("whole");
    let client = reqwest::Client::new();

    let e = entry(&addr, "x.cram", Some(md5_hex(&data)), Some(data.len() as i64));
    let mut seen: u64 = 0;
    let got = ena::fetch_file(&client, dir.path(), &e, &CancelToken::none(), &mut |r, _| seen = r)
        .await
        .expect("download");

    assert_eq!(std::fs::read(&got).unwrap(), data);
    assert_eq!(seen, data.len() as u64, "progress ends at the full size");
    assert!(
        !got.with_extension("cram.part").exists(),
        "the .part is renamed away, never left behind"
    );
    stop.store(true, Ordering::Relaxed);
}

/// This is the purpose of the module. An interrupted transfer continues from the bytes on the
/// disk. The md5 value is still correct, but this process did not get the prefix.
#[tokio::test]
async fn an_interrupted_transfer_resumes_and_still_verifies() {
    let data = body();
    let stop = Arc::new(AtomicBool::new(false));
    let addr = serve(data.clone(), Mode::Honest, stop.clone()).await;
    let dir = Scratch::new("resume");
    let client = reqwest::Client::new();

    // Simulate a killed download: the first third is already on disk as a `.part`.
    let split = data.len() / 3;
    std::fs::write(dir.path().join("x.cram.part"), &data[..split]).unwrap();

    let e = entry(&addr, "x.cram", Some(md5_hex(&data)), Some(data.len() as i64));
    let mut first_report: Option<u64> = None;
    let got = ena::fetch_file(&client, dir.path(), &e, &CancelToken::none(), &mut |r, _| {
        first_report.get_or_insert(r);
    })
    .await
    .expect("resumed download");

    assert_eq!(
        std::fs::read(&got).unwrap(),
        data,
        "the resumed file is byte-identical to the source"
    );
    assert!(
        first_report.unwrap() > split as u64,
        "progress counts the resumed prefix; a bar that restarts at zero tells the user the opposite \
         of what happened"
    );
    stop.store(true, Ordering::Relaxed);
}

/// Many servers ignore `Range` and send the full body with `200`. If the module added those bytes
/// to the prefix, the file would have the correct size only by accident. The checksum would fail.
#[tokio::test]
async fn a_server_that_ignores_range_is_handled_rather_than_trusted() {
    let data = body();
    let stop = Arc::new(AtomicBool::new(false));
    let addr = serve(data.clone(), Mode::IgnoresRange, stop.clone()).await;
    let dir = Scratch::new("ignores-range");
    let client = reqwest::Client::new();

    std::fs::write(dir.path().join("x.cram.part"), &data[..data.len() / 3]).unwrap();

    let e = entry(&addr, "x.cram", Some(md5_hex(&data)), Some(data.len() as i64));
    let got = ena::fetch_file(&client, dir.path(), &e, &CancelToken::none(), &mut |_, _| {})
        .await
        .expect("download restarts cleanly");
    assert_eq!(
        std::fs::read(&got).unwrap(),
        data,
        "restarted from zero, not appended to the prefix"
    );
    stop.store(true, Ordering::Relaxed);
}

/// A bad checksum must fail loudly and leave nothing behind that a later run could mistake for a
/// good file.
#[tokio::test]
async fn a_checksum_mismatch_fails_and_leaves_no_part_behind() {
    let data = body();
    let stop = Arc::new(AtomicBool::new(false));
    let addr = serve(data.clone(), Mode::Corrupt, stop.clone()).await;
    let dir = Scratch::new("corrupt");
    let client = reqwest::Client::new();

    let e = entry(
        &addr,
        "x.cram",
        Some(md5_hex(b"something else entirely")),
        Some(data.len() as i64),
    );
    // No delay. This test examines the failure path, and not the true delay schedule.
    let policy = RetryPolicy {
        attempts: 2,
        backoff: false,
    };
    let err = ena::fetch_file_with(&client, dir.path(), &e, policy, &CancelToken::none(), &mut |_, _| {})
        .await
        .expect_err("must not accept a file that fails its checksum");
    assert!(format!("{err}").contains("checksum mismatch"), "{err}");
    assert!(!dir.path().join("x.cram").exists(), "no finished file");
    assert!(
        !dir.path().join("x.cram.part").exists(),
        "and no partial one to resume from"
    );
    stop.store(true, Ordering::Relaxed);
}

/// The module never gets a file two times. The rename marks a file as complete. So a second run of
/// a unit after a crash costs nothing for the files that are already complete.
#[tokio::test]
async fn a_completed_file_is_not_downloaded_twice() {
    let dir = Scratch::new("existing");
    std::fs::write(dir.path().join("x.cram"), b"already here").unwrap();
    let client = reqwest::Client::new();

    // The URL points to no server. A request on the network would make this test fail.
    let e = entry("127.0.0.1:1", "x.cram", Some("ignored".into()), Some(999));
    let got = ena::fetch_file(&client, dir.path(), &e, &CancelToken::none(), &mut |_, _| {})
        .await
        .expect("an existing file short-circuits");
    assert_eq!(std::fs::read(got).unwrap(), b"already here");
}

//! Time each phase of a region read that uses an index. A slow path over a BAM or a CRAM then goes
//! to the part that caused it. Those parts are the open, the header, the first query, a warm
//! query, and the walk over the whole region. Without this, the blame lands on whichever call the
//! wall clock happened to be inside.
//!
//! Written to diagnose a CRAM that took ~500x longer than an equivalent BAM for one region query.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example cram_query_probe -- \
//!   <bam/cram> <ref.fa> <contig> <pos> [span_bp=1000000]
//! ```

use std::path::Path;
use std::time::Instant;

use navigator_analysis::reader::open_indexed;
use noodles::core::Region;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: cram_query_probe <bam/cram> <ref.fa> <contig> <pos> [span_bp=1000000]");
        std::process::exit(2);
    }
    let (path, refp, contig) = (Path::new(&a[1]), Path::new(&a[2]), a[3].clone());
    let pos: usize = a[4].parse().expect("pos");
    let span: usize = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(1_000_000);

    let t0 = Instant::now();
    let (header, mut reader) = open_indexed(path, Some(refp)).expect("open");
    println!("open + header            : {:>8.2?}", t0.elapsed());

    // A region of one base. The cost here is the overhead of one query, and not the work at each
    // record.
    let one = |p: usize| -> Region { format!("{contig}:{p}-{p}").parse().expect("region") };

    let t = Instant::now();
    let n: usize = reader.query(&header, &one(pos)).expect("q1").count();
    println!(
        "first 1bp query ({n:>4} rec) : {:>8.2?}   <- includes any lazy setup",
        t.elapsed()
    );

    let t = Instant::now();
    let n: usize = reader.query(&header, &one(pos)).expect("q2").count();
    println!(
        "same query again ({n:>4} rec): {:>8.2?}   <- warm: is the cost per-query or one-off?",
        t.elapsed()
    );

    let t = Instant::now();
    let n: usize = reader.query(&header, &one(pos + 5_000_000)).expect("q3").count();
    println!(
        "distant 1bp query ({n:>4} rec): {:>8.2?}   <- new container: does it re-decode?",
        t.elapsed()
    );

    let region: Region = format!("{contig}:{pos}-{}", pos + span).parse().expect("region");
    let t = Instant::now();
    let n: usize = reader.query(&header, &region).expect("bulk").count();
    let el = t.elapsed();
    println!(
        "{:.1} Mb region ({n} rec)  : {:>8.2?}  = {:.2?}/Mb",
        span as f64 / 1e6,
        el,
        el / (span as u32 / 1_000_000).max(1)
    );

    // With `VERIFY=1` the probe checks the query that skips containers against the own `Query` of
    // noodles, which walks the whole contig, on REAL data.
    //
    // The fixture in the repo holds one container. Only a large CRAM with many containers can
    // catch a container that the code skips wrongly. Such a fault would look like a faster caller,
    // and not a broken one. This check is slow by construction, because the oracle is the code
    // that we replaced.
    if std::env::var("VERIFY").is_ok_and(|v| v == "1") {
        use noodles::cram;

        let key = |r: &noodles::sam::alignment::RecordBuf| {
            (
                r.name().map(|n| n.to_vec()),
                r.alignment_start().map(usize::from),
                r.flags().bits(),
                r.sequence().as_ref().to_vec(),
            )
        };
        let mine: Vec<_> = reader
            .query(&header, &region)
            .expect("mine")
            .map(|r| key(&r.expect("rec")))
            .collect();

        let repo = navigator_analysis::reader::build_repository(refp).expect("repo");
        let mut oracle = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_path(path)
            .expect("noodles open");
        let oh = oracle.read_header().expect("noodles header");
        let t = Instant::now();
        let theirs: Vec<_> = oracle
            .query(&oh, &region)
            .expect("noodles query")
            .map(|r| key(&r.expect("rec")))
            .collect();
        println!(
            "\nVERIFY: noodles' own Query took {:?} for the same region",
            t.elapsed()
        );
        println!("        ours {} records, noodles {} records", mine.len(), theirs.len());
        assert_eq!(mine, theirs, "container skipping changed the records returned");
        println!("        IDENTICAL — container skipping is lossless");
    }
}

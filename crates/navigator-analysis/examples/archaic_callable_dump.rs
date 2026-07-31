//! Dump the Tier B callability mask as BED, so what the segment caller can and cannot see is
//! checkable against an external callset rather than assumed.
//!
//! Windows below `min_frac` of `window_bp` callable are excluded by the caller itself, so the same
//! threshold is applied here — the output is the territory a segment could actually be called in.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_callable_dump -- \
//!   ~/.decodingus/ancestry/archaic_callable_chm13v2.0.bin 0.5 chr21 chr22 > callable.bed
//! ```

use navigator_analysis::archaic::ArchaicCallable;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let path = a.next().expect("usage: archaic_callable_dump <callable.bin> [min_frac] [contig ...]");
    let min_frac: f64 = a.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let want: Vec<String> = a.collect();

    let cal = ArchaicCallable::from_bytes(&std::fs::read(&path)?).map_err(|e| e.to_string())?;
    eprintln!(
        "callable track: {} contigs, window {} bp, {:.1} Mb callable total",
        cal.contigs.len(),
        cal.window_bp,
        cal.callable_mb()
    );

    let w = cal.window_bp;
    for c in &cal.contigs {
        if !want.is_empty() && !want.contains(&c.contig) {
            continue;
        }
        let (mut run_start, mut in_run) = (0i64, false);
        for (i, &bp) in c.callable_bp.iter().enumerate() {
            let pos = c.start + (i as i64) * w;
            let ok = (bp as f64) / (w as f64) >= min_frac;
            if ok && !in_run {
                run_start = pos;
                in_run = true;
            } else if !ok && in_run {
                println!("{}\t{}\t{}", c.contig, run_start, pos);
                in_run = false;
            }
        }
        if in_run {
            let end = c.start + (c.callable_bp.len() as i64) * w;
            println!("{}\t{}\t{}", c.contig, run_start, end);
        }
    }
    Ok(())
}

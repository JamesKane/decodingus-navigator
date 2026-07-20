//! Standalone de-novo caller profiling harness — runs [`caller::call_denovo`] on one contig
//! of a BAM/CRAM with no async/test wrapper, so a sampling profiler sees only the hot path.
//!
//! ```sh
//! BAM=/Users/jkane/Genomics/WGS229/WGS229.bwa-mem.chm13v2.cram \
//! REF=$HOME/.decodingus/references/chm13v2.0.fa CONTIG=chrY \
//! cargo run --release -p navigator-analysis --example denovo_profile
//! # profile (after `cargo install samply`):
//! #   samply record target/release/examples/denovo_profile
//! ```

use std::path::Path;
use std::time::Instant;

use navigator_analysis::caller::{call_denovo, HaploidCallerParams};

fn main() {
    let bam = std::env::var("BAM").expect("set BAM=path/to.{bam,cram}");
    let reference = std::env::var("REF").expect("set REF=path/to.fa");
    let contig = std::env::var("CONTIG").unwrap_or_else(|_| "chrY".to_string());

    let params = HaploidCallerParams::default();
    let t = Instant::now();
    let calls = call_denovo(
        Path::new(&bam),
        Path::new(&reference),
        &contig,
        &params,
        &navigator_analysis::CancelToken::none(),
    ).expect("call_denovo");
    eprintln!(
        "call_denovo({contig}): {} variants in {:.1}s (realign={})",
        calls.len(),
        t.elapsed().as_secs_f64(),
        params.local_realign,
    );
}

//! A profile harness for the de-novo caller, on its own. It runs [`caller::call_denovo`] over one
//! contig of a BAM or a CRAM. There is no async wrapper and no test wrapper, so a profiler that
//! samples sees the hot path alone.
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
    )
    .expect("call_denovo");
    eprintln!(
        "call_denovo({contig}): {} variants in {:.1}s (realign={})",
        calls.len(),
        t.elapsed().as_secs_f64(),
        params.local_realign,
    );
}

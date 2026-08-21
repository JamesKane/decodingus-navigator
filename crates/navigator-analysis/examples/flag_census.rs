//! Count the records of a BAM by their class of SAM flag: primary, secondary, supplementary, and
//! unmapped.
//!
//! A realigned alignment must hold about one primary record for each input read. When it holds some
//! times that count, the extra records are the alternative placements that the mapper chose. This
//! census says so in one line. Without it, a record count that looks wrong leads to an argument.
//!
//! ```sh
//! BAM=~/.decodingus/realigned/alignment-8.chm13v2.0.bam \
//!   cargo run --release -p navigator-analysis --example flag_census
//! ```

use std::path::Path;
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::var("BAM")?;
    let (header, mut reader) = navigator_analysis::reader::open_seq(Path::new(&path), None)?;

    let (mut total, mut primary, mut secondary, mut supplementary, mut unmapped, mut no_seq) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    let started = Instant::now();

    for result in reader.records(&header) {
        let record = result?;
        let flags = record.flags();
        total += 1;
        if flags.is_secondary() {
            secondary += 1;
        } else if flags.is_supplementary() {
            supplementary += 1;
        } else {
            primary += 1;
            if flags.is_unmapped() {
                unmapped += 1;
            }
        }
        if record.sequence().as_ref().is_empty() {
            no_seq += 1;
        }
        if total % 50_000_000 == 0 {
            eprintln!("  {total} records…");
        }
    }

    let pct = |n: u64| 100.0 * n as f64 / total.max(1) as f64;
    println!("{path}");
    println!("  total          {total:>13}");
    println!("  primary        {primary:>13}  ({:.1}%)", pct(primary));
    println!("    of which unmapped {unmapped:>8}");
    println!("  secondary      {secondary:>13}  ({:.1}%)", pct(secondary));
    println!("  supplementary  {supplementary:>13}  ({:.1}%)", pct(supplementary));
    println!("  no SEQ         {no_seq:>13}  ({:.1}%)", pct(no_seq));
    println!("  read in {:.1}s", started.elapsed().as_secs_f64());
    Ok(())
}

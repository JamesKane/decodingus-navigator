//! Validate the phase-2 reassembly wiring end-to-end on a real CRAM: for each position, run the
//! bounded de-novo caller with reassembly OFF (pileup only) vs ON, and report whether the position
//! is now called. Confirms the misaligned-ref truth privates the paralog gate drops are recovered.
//!
//!   cargo run --release --example reassembly_validate -p navigator-analysis -- \
//!       <bam/cram> <ref.fa> chrY <pos[,pos,...]>

use std::path::Path;

use navigator_analysis::caller::{call_denovo_region, HaploidCallerParams};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: reassembly_validate <bam/cram> <ref.fa> <contig> <pos[,pos,...]>");
        std::process::exit(2);
    }
    let bam = Path::new(&args[1]);
    let refp = Path::new(&args[2]);
    let contig = &args[3];
    let positions: Vec<i64> = args[4].split(',').filter_map(|s| s.trim().parse().ok()).collect();

    println!("{:>10}  {:>10}  {:>10}  recovered?", "pos", "pileup", "reassembly");
    for &pos in &positions {
        let lo = (pos - 5).max(1) as usize;
        let hi = (pos + 5) as usize;

        let called = |reassembly: bool| -> Option<(char, char, u32, u32, Option<f64>)> {
            let params = HaploidCallerParams { reassembly, ..HaploidCallerParams::default() };
            let calls = call_denovo_region(bam, refp, contig, lo, hi, &params).expect("call_denovo_region");
            calls
                .into_iter()
                .find(|c| c.position == pos)
                .map(|c| (c.reference_allele, c.alternate_allele, c.alt_depth, c.depth, c.quality))
        };

        let off = called(false);
        let on = called(true);
        let fmt = |c: &Option<(char, char, u32, u32, Option<f64>)>| match c {
            Some((r, a, ad, dp, q)) => {
                let gq = q.map(|v| format!(" GQ{v:.0}")).unwrap_or_default();
                format!("{r}>{a} {ad}/{dp}{gq}")
            }
            None => "-".to_string(),
        };
        let recovered = if off.is_none() && on.is_some() {
            "YES (reassembly recovered)"
        } else if on.is_some() {
            "called by both"
        } else {
            "still missing"
        };
        println!("{pos:>10}  {:>10}  {:>10}  {recovered}", fmt(&off), fmt(&on));
    }
}

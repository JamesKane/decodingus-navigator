//! Dump the **private** variant positions the Tier B HMM actually sees — the subject's derived
//! variants after the African-outgroup strip — so the input to the model can be checked against an
//! external truth set independently of the model.
//!
//! The segment caller is a density model over exactly these positions. If they are not enriched
//! inside known archaic tracts, no amount of HMM tuning can help, and the fault is upstream in the
//! variant calls or the outgroup strip rather than in the model. That question is unanswerable from
//! the caller's own output, which is why this exists.
//!
//! ```sh
//! cargo run --release -p navigator-analysis --example archaic_private_dump -- \
//!   calls.json ~/.decodingus/ancestry/archaic_outgroup_af_chm13v2.0.bin > private.tsv
//! ```

use navigator_analysis::archaic::ArchaicOutgroup;
use navigator_analysis::caller::SiteGenotype;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut a = std::env::args().skip(1);
    let calls_path = a
        .next()
        .expect("usage: archaic_private_dump <calls.json> <outgroup.bin>");
    let og_path = a
        .next()
        .expect("usage: archaic_private_dump <calls.json> <outgroup.bin>");

    let calls: Vec<SiteGenotype> = serde_json::from_str(&std::fs::read_to_string(&calls_path)?)?;
    let og = ArchaicOutgroup::from_bytes(&std::fs::read(&og_path)?).map_err(|e| e.to_string())?;

    // Group by contig, mirroring what the caller does before it strips.
    let mut by_contig: std::collections::BTreeMap<String, Vec<&SiteGenotype>> = Default::default();
    for c in &calls {
        by_contig.entry(c.contig.clone()).or_default().push(c);
    }

    // Quality columns come out too: whether the background's excess variance is real biology or
    // this caller's own error rate varying by region is not answerable without them.
    println!("contig\tposition\tdosage\tgq\tdepth");
    for (contig, mut sites) in by_contig {
        sites.sort_by_key(|s| s.position);
        let carried: Vec<&SiteGenotype> = sites.iter().copied().filter(|s| s.dosage > 0).collect();
        let positions: Vec<i64> = carried.iter().map(|s| s.position).collect();
        let keep: std::collections::HashSet<i64> = og.retain_private(&contig, &positions).into_iter().collect();
        let mut kept = 0usize;
        for s in &carried {
            if keep.contains(&s.position) {
                println!("{contig}\t{}\t{}\t{}\t{}", s.position, s.dosage, s.gq, s.depth);
                kept += 1;
            }
        }
        eprintln!("{contig}: {kept} private of {} carried variants", carried.len());
    }
    Ok(())
}

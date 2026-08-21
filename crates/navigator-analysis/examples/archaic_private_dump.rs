//! Write out the **private** variant positions that the Tier B HMM sees. Those are the derived
//! variants of the subject, after the code removes the ones that the African outgroup also carries.
//! Somebody can then check the input of the model against an external truth set, and that check
//! does not depend on the model.
//!
//! The segment caller is a density model over exactly these positions. If they show no enrichment
//! inside a known archaic tract, then no change to the HMM can help. The fault would lie earlier,
//! in the variant calls or in the removal of the outgroup sites, and not in the model. The own
//! output of the caller can not answer that question, and that is why this tool exists.
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

    // Put them into groups by contig, as the caller does before it removes the outgroup sites.
    let mut by_contig: std::collections::BTreeMap<String, Vec<&SiteGenotype>> = Default::default();
    for c in &calls {
        by_contig.entry(c.contig.clone()).or_default().push(c);
    }

    // The quality columns come out too. Without them, nobody can answer one question. Is the
    // excess variance of the background real biology, or is it the own error rate of this caller,
    // which changes from region to region?
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

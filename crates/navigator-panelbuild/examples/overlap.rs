use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashSet;
fn main() -> anyhow::Result<()> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let p1 = AncestryPanel::from_bytes(&std::fs::read(&a[0])?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let p2 = AncestryPanel::from_bytes(&std::fs::read(&a[1])?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let s1: HashSet<(String, i64)> = p1.sites.iter().map(|s| (s.contig.clone(), s.position)).collect();
    let s2: HashSet<(String, i64)> = p2.sites.iter().map(|s| (s.contig.clone(), s.position)).collect();
    println!("{} sites={} pops={:?}", a[0], p1.sites.len(), p1.populations);
    println!("{} sites={} pops={:?}", a[1], p2.sites.len(), p2.populations);
    println!("overlap = {}", s1.intersection(&s2).count());
    println!("p1 sample contig: {:?}", p1.sites.first().map(|s| (&s.contig, s.position)));
    println!("p2 sample contig: {:?}", p2.sites.first().map(|s| (&s.contig, s.position)));
    Ok(())
}

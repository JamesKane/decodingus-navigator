//! Apply an ascertainment floor to an already-built ancient AF panel: keep only sites whose CHM13
//! (contig,pos) appears in a `contig<TAB>pos` set. Equivalent to rebuilding the panel with
//! `ancient-panel --ascertain-sites` on the same inputs (the per-site frequencies are unchanged;
//! only the site set is restricted), for when the source matrices aren't at hand but the full panel
//! is. Usage: filter_sites <ancient.bin> <sites.tsv> <out.bin>
use navigator_analysis::ancestry::AncestryPanel;
use std::collections::HashSet;

fn main() -> anyhow::Result<()> {
    let ancient = std::env::args().nth(1).expect("usage: filter_sites <ancient.bin> <sites.tsv> <out.bin>");
    let sites_tsv = std::env::args().nth(2).expect("sites.tsv");
    let out = std::env::args().nth(3).expect("out.bin");
    let keep: HashSet<(String, i64)> = std::fs::read_to_string(&sites_tsv)?
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .filter_map(|l| {
            let mut it = l.split('\t');
            let c = it.next()?.trim();
            let p: i64 = it.next()?.trim().parse().ok()?;
            (!c.eq_ignore_ascii_case("contig")).then(|| (c.to_string(), p))
        })
        .collect();
    let mut panel = AncestryPanel::from_bytes(&std::fs::read(&ancient)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let before = panel.sites.len();
    panel.sites.retain(|s| keep.contains(&(s.contig.clone(), s.position)));
    let after = panel.sites.len();
    std::fs::write(&out, panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    println!("kept {after}/{before} sites (ascertainment set {}) -> {out}", keep.len());
    Ok(())
}

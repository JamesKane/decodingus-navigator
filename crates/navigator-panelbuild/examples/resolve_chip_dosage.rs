//! Resolve a consumer-chip raw-data file (23andMe / AncestryDNA) to CHM13-oriented panel dosages,
//! for the panel-depth sweep's WGS-vs-chip STABILITY gate (sibling of `score_modern_from_tsv`).
//! Reuses the *production* chip path — `chipprofile::{detect_build, autosomal_calls}` +
//! `IbdPanel::resolve_chip` — so the dosages are byte-for-byte what the app would feed the ancestry
//! estimators (palindromes dropped, alleles oriented to the CHM13 ref/alt). Emits
//! `contig<TAB>pos<TAB>dosage` (0/1/2) that `score_modern_from_tsv` reads directly.
//!   resolve_chip_dosage <ibd_panel.bin> <chip.txt> <out.tsv>
use navigator_analysis::ibd_panel::IbdPanel;
use navigator_domain::chipprofile;
use std::io::{BufWriter, Write};

fn main() -> anyhow::Result<()> {
    let ibd_path = std::env::args().nth(1).expect("usage: resolve_chip_dosage <ibd_panel.bin> <chip.txt> <out.tsv>");
    let chip_path = std::env::args().nth(2).expect("chip.txt");
    let out = std::env::args().nth(3).expect("out.tsv");

    let text = std::fs::read_to_string(&chip_path)?;
    let build = chipprofile::detect_build(&text);
    let calls = chipprofile::autosomal_calls(&text);
    eprintln!("chip build {build}: {} autosomal calls", calls.len());

    let ibd = IbdPanel::from_bytes(&std::fs::read(&ibd_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;
    let tuples: Vec<(String, i64, char, char)> = calls.into_iter().map(|c| (c.contig, c.position, c.a1, c.a2)).collect();
    let gts = ibd.resolve_chip(&build, &tuples);

    let mut w = BufWriter::new(std::fs::File::create(&out)?);
    let mut called = 0usize;
    for g in &gts {
        if g.dosage >= 0 {
            called += 1;
        }
        writeln!(w, "{}\t{}\t{}", g.contig, g.position, g.dosage)?;
    }
    w.flush()?;
    eprintln!("wrote {out}: {} panel dosages ({} called)", gts.len(), called);
    Ok(())
}

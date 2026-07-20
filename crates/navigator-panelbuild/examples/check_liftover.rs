//! Throwaway: does the IBD panel's per-build locus actually point at the SNP on that build's genome?
//! For a spread of sites, read the reference base at the panel's stored (contig,pos) on both CHM13
//! and GRCh38, and check it equals the panel's stored REF for that build. A correct locus has
//! genome[pos] == panel.REF (or panel.ALT if that build's strand is flipped). Wholesale mismatch on
//! GRCh38 but not CHM13 == the GRCh38 liftover coordinates are wrong.
use navigator_analysis::ibd_panel::IbdPanel;
use noodles::core::Region;
use noodles::fasta;

fn base_at<R: std::io::BufRead + std::io::Seek>(
    reader: &mut fasta::io::IndexedReader<R>,
    contig: &str,
    pos: i64,
) -> Option<char> {
    let region: Region = format!("{contig}:{pos}-{pos}").parse().ok()?;
    let rec = reader.query(&region).ok()?;
    rec.sequence().as_ref().first().map(|&b| (b as char).to_ascii_uppercase())
}

fn main() -> anyhow::Result<()> {
    let panel_path = std::env::args().nth(1).expect("usage: check_liftover <ibd_panel.bin> <chm13.fa> <grch38.fa>");
    let chm13_fa = std::env::args().nth(2).expect("chm13.fa");
    let grch38_fa = std::env::args().nth(3).expect("grch38.fa");
    let panel = IbdPanel::from_bytes(&std::fs::read(&panel_path)?).map_err(|e| anyhow::anyhow!("{e}"))?;

    let open = |p: &str| fasta::io::indexed_reader::Builder::default().build_from_path(p);
    let mut chm = open(&chm13_fa)?;
    let mut g38 = open(&grch38_fa)?;

    // Try the stored contig name, then chr/bare variants, so a naming mismatch doesn't masquerade
    // as a coordinate error.
    let g38_base = |g38: &mut _, l: &navigator_analysis::ibd_panel::Locus| -> Option<char> {
        let bare = l.contig.strip_prefix("chr").unwrap_or(&l.contig).to_string();
        base_at(g38, &l.contig, l.position)
            .or_else(|| base_at(g38, &bare, l.position))
            .or_else(|| base_at(g38, &format!("chr{bare}"), l.position))
    };

    let (mut chm_ok, mut chm_n) = (0u64, 0u64);
    let (mut g38_ref, mut g38_alt, mut g38_other, mut g38_n) = (0u64, 0u64, 0u64, 0u64);
    let mut examples = Vec::new();

    let step = (panel.sites.len() / 6000).max(1);
    for s in panel.sites.iter().step_by(step) {
        if let Some(b) = base_at(&mut chm, &s.chm13.contig, s.chm13.position) {
            chm_n += 1;
            if b == s.chm13.reference.to_ascii_uppercase() || b == s.chm13.alternate.to_ascii_uppercase() {
                chm_ok += 1;
            }
        }
        if let Some(l) = s.grch38.as_ref() {
            if let Some(b) = g38_base(&mut g38, l) {
                g38_n += 1;
                if b == l.reference.to_ascii_uppercase() {
                    g38_ref += 1;
                } else if b == l.alternate.to_ascii_uppercase() {
                    g38_alt += 1;
                } else {
                    g38_other += 1;
                    if examples.len() < 12 {
                        examples.push(format!(
                            "{} chm13 {}:{} {}/{}  ->  grch38 {}:{} REF={} ALT={}  genome={}",
                            s.rsid, s.chm13.contig, s.chm13.position, s.chm13.reference, s.chm13.alternate,
                            l.contig, l.position, l.reference, l.alternate, b
                        ));
                    }
                }
            }
        }
    }

    println!("sampled every {step}th of {} sites\n", panel.sites.len());
    println!(
        "CHM13 control : genome base == panel REF|ALT at {}/{} ({:.1}%)",
        chm_ok, chm_n, 100.0 * chm_ok as f64 / chm_n.max(1) as f64
    );
    println!("GRCh38 locus  : n={g38_n}");
    println!("  genome == build REF : {g38_ref} ({:.1}%)", 100.0 * g38_ref as f64 / g38_n.max(1) as f64);
    println!("  genome == build ALT : {g38_alt} ({:.1}%)", 100.0 * g38_alt as f64 / g38_n.max(1) as f64);
    println!(
        "  genome == NEITHER   : {g38_other} ({:.1}%)  <- wrong coordinate",
        100.0 * g38_other as f64 / g38_n.max(1) as f64
    );
    if !examples.is_empty() {
        println!("\nexamples where GRCh38 genome base matches neither stored allele:");
        for e in &examples {
            println!("  {e}");
        }
    }
    Ok(())
}

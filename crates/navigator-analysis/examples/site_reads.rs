//! Dump per-read detail at a site so we can see *why* a marginal pileup is ambiguous — the input to
//! designing v2 active-region read selection. For each spanning read: query name, mate/pair flags,
//! MAPQ, the base + base-quality it carries at the site, and its mismatch count vs the reference
//! window (excluding the site itself).
//!
//!   cargo run --release --example site_reads -p navigator-analysis -- <bam/cram> <ref.fa> chrY <pos> [win=40]

use std::collections::HashMap;
use std::path::Path;

use navigator_analysis::reader::{open_indexed, read_contig_sequence};
use noodles::core::Region;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 5 {
        eprintln!("usage: site_reads <bam/cram> <ref.fa> <contig> <pos> [win=40]");
        std::process::exit(2);
    }
    let (bam, refp, contig) = (Path::new(&a[1]), Path::new(&a[2]), &a[3]);
    let pos: i64 = a[4].parse().expect("pos");
    let win: i64 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(40);
    let (lo, hi) = ((pos - win).max(1), pos + win);

    let mut refseq = read_contig_sequence(refp, contig).expect("ref");
    refseq.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
    let ref_base = refseq[(pos - 1) as usize] as char;

    let (header, mut reader) = open_indexed(bam, Some(refp)).expect("open");
    let region: Region = format!("{contig}:{lo}-{hi}").parse().expect("region");

    println!("site chrY:{pos} ref={ref_base}  window ±{win}");
    println!("{:<32} {:>4} {:>4} {:>5} {:>4} {:>3} {:>8}", "qname", "pair", "mapq", "site", "bq", "nm", "flags");
    let mut base_tally: HashMap<char, u32> = HashMap::new();
    for result in reader.query(&header, &region).expect("query") {
        let rec = result.expect("rec");
        let f = rec.flags();
        if f.is_secondary() || f.is_supplementary() || f.is_duplicate() || f.is_unmapped() {
            continue;
        }
        let Some(start) = rec.alignment_start().map(|p| p.get() as i64) else { continue };
        let seq = rec.sequence();
        let quals = rec.quality_scores();
        let qb = quals.as_ref();
        let name = rec.name().map(|n| String::from_utf8_lossy(n).into_owned()).unwrap_or_default();
        let mapq = rec.mapping_quality().map_or(255, |m| m.get());

        // Walk the CIGAR: capture the base at `pos` and count mismatches vs ref (excluding `pos`).
        let mut ref_pos = start;
        let mut qoff = 0usize;
        let mut site_base = '.';
        let mut site_bq = 0u8;
        let mut nm = 0u32;
        let mut mm_pos: Vec<i64> = Vec::new();
        for op in rec.cigar().as_ref() {
            let (cr, cq) = (op.kind().consumes_reference(), op.kind().consumes_read());
            let len = op.len();
            if cr && cq {
                for i in 0..len {
                    let rp = ref_pos + i as i64;
                    if rp >= lo && rp <= hi {
                        if let Some(b) = seq.get(qoff + i) {
                            let b = b.to_ascii_uppercase();
                            let rb = refseq[(rp - 1) as usize];
                            if rp == pos {
                                site_base = b as char;
                                site_bq = qb.get(qoff + i).copied().unwrap_or(0);
                            } else if b != rb {
                                nm += 1;
                                mm_pos.push(rp);
                            }
                        }
                    }
                }
                ref_pos += len as i64;
                qoff += len;
            } else if cr {
                ref_pos += len as i64;
            } else if cq {
                qoff += len;
            }
        }
        if site_base == '.' {
            continue; // doesn't span the site
        }
        *base_tally.entry(site_base).or_default() += 1;
        let pair = if f.is_first_segment() { "R1" } else if f.is_last_segment() { "R2" } else { "?" };
        let mm: Vec<String> = mm_pos.iter().map(|p| (p - pos).to_string()).collect();
        println!("{name:<38} {pair:>4} {mapq:>4} {site_base:>5} {site_bq:>4} {nm:>3}  mm@[{}]", mm.join(","));
    }
    let mut keys: Vec<_> = base_tally.keys().copied().collect();
    keys.sort();
    let summary: Vec<String> = keys.iter().map(|k| format!("{k}={}", base_tally[k])).collect();
    println!("site pileup: {}", summary.join(" "));
}

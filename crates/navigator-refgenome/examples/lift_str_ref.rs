//! Lift a HipSTR reference BED to another build via the cached chain.
//! Usage: cargo run -p navigator-refgenome --example lift_str_ref -- <in.bed.gz> <out.bed.gz> [from] [to] [contig]
use navigator_refgenome::{cache, gateway::ReferenceGateway};

#[tokio::main]
async fn main() {
    let mut a = std::env::args().skip(1);
    let in_bed = std::path::PathBuf::from(a.next().expect("in.bed.gz"));
    let out_bed = std::path::PathBuf::from(a.next().expect("out.bed.gz"));
    let from = a.next().unwrap_or_else(|| "GRCh38".into());
    let to = a.next().unwrap_or_else(|| "chm13v2.0".into());
    let contig = a.next();

    let gw = ReferenceGateway::new(cache::base_dir(), reqwest::Client::new());
    eprintln!("resolving {from} -> {to} chain (downloads on miss)…");
    gw.resolve_chain(&from, &to, &mut |d, t| {
        if d % (32 * 1024 * 1024) < 1_048_576 {
            eprintln!(
                "  chain {} MB{}",
                d / 1_048_576,
                t.map(|t| format!(" / {} MB", t / 1_048_576)).unwrap_or_default()
            );
        }
    })
    .await
    .expect("resolve chain");

    let t = std::time::Instant::now();
    let stats = gw
        .lift_hipstr_bed(&from, &to, &in_bed, &out_bed, contig.as_deref())
        .expect("lift");
    eprintln!("done in {:.1}s", t.elapsed().as_secs_f32());
    println!(
        "{:?}\n  lifted {}/{} ({:.0}%)  dropped: {} unmapped, {} split, {} span",
        stats,
        stats.lifted,
        stats.total,
        100.0 * stats.lifted as f64 / stats.total.max(1) as f64,
        stats.dropped_unmapped,
        stats.dropped_split,
        stats.dropped_span
    );
}

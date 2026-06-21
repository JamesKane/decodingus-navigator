//! One-off placement debug: genotype a BAM at every Y-tree position and dump the ranked
//! placement + what `deepen_terminal` returns + the calls along a named lineage.
//!
//!   cargo run --release --example place_y -p navigator-analysis -- <bam> <ref.fa> <tree.json> [FOCUS_NODE]

use std::collections::HashSet;
use std::path::Path;

use std::collections::HashMap;

use navigator_analysis::caller::{call_bases_at, HaploidCallerParams};
use navigator_analysis::haplo::{self, Locus};

/// (derived, ancestral, nocall) for a node's loci against the calls.
fn counts(loci: &[Locus], calls: &HashMap<i64, char>) -> (usize, usize, usize) {
    let (mut d, mut a, mut n) = (0, 0, 0);
    for l in loci {
        let der = l.derived.chars().next().map(|c| c.to_ascii_uppercase());
        let anc = l.ancestral.chars().next().map(|c| c.to_ascii_uppercase());
        match (der, calls.get(&l.position).map(|c| c.to_ascii_uppercase())) {
            (None, _) => n += 1,
            (Some(_), None) => n += 1,
            (Some(dd), Some(b)) if b == dd => d += 1,
            (Some(_), Some(b)) if Some(b) == anc => a += 1,
            (Some(_), Some(_)) => a += 1,
        }
    }
    (d, a, n)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: place_y <bam> <ref.fa> <tree.json> [FOCUS_NODE]");
        std::process::exit(2);
    }
    let (bam, reference, tree_path) = (Path::new(&a[1]), Path::new(&a[2]), &a[3]);
    let focus = a.get(4).cloned().unwrap_or_else(|| "R-CTS4466".to_string());

    let tree_json = std::fs::read_to_string(tree_path).unwrap();
    let tree = haplo::parse_ftdna_json(&tree_json).unwrap();
    let targets: HashSet<i64> = tree
        .nodes
        .values()
        .flat_map(|n| n.loci.iter().map(|l| l.position))
        .collect();
    eprintln!("tree nodes={} target positions={}", tree.nodes.len(), targets.len());

    let params = HaploidCallerParams::default();
    let calls = call_bases_at(bam, "chrY", &targets, &params, Some(reference)).unwrap();
    eprintln!("called {} / {} positions", calls.len(), targets.len());

    let ranked = haplo::score(&tree, &calls);
    println!("\n=== top 15 by Kulczynski (score, depth) ===");
    for r in ranked.iter().take(15) {
        println!(
            "  {:<16} score={:.4} depth={:>2} matched={}/{} found={} admissible={}",
            r.name,
            r.score,
            r.depth,
            r.matched,
            r.expected,
            r.found,
            haplo::path_admissible(&tree, &calls, r.id),
        );
    }

    let start = ranked
        .iter()
        .find(|r| haplo::path_admissible(&tree, &calls, r.id))
        .map(|r| (r.id, r.name.clone()));
    if let Some((sid, sname)) = start {
        let term = haplo::deepen_terminal(&tree, &calls, sid);
        let tname = tree.nodes.get(&term).map(|n| n.name.as_str()).unwrap_or("?");
        println!("\nfirst admissible start = {sname}; deepen_terminal -> {tname}");
    }

    // Walk root->focus and print the call state of each node's defining SNPs.
    let byname: std::collections::HashMap<&str, i64> = tree.nodes.values().map(|n| (n.name.as_str(), n.id)).collect();
    if let Some(&fid) = byname.get(focus.as_str()) {
        let mut parent = std::collections::HashMap::new();
        for n in tree.nodes.values() {
            for &c in &n.children {
                parent.insert(c, n.id);
            }
        }
        let mut path = Vec::new();
        let mut cur = Some(fid);
        while let Some(id) = cur {
            path.push(id);
            cur = parent.get(&id).copied();
        }
        path.reverse();
        println!("\n=== call states along root -> {focus} ===");
        for id in path {
            let n = &tree.nodes[&id];
            let (d, anc, nc) = counts(&n.loci, &calls);
            println!(
                "  {:<16} loci={:>2} derived={} ancestral={} nocall={}",
                n.name,
                n.loci.len(),
                d,
                anc,
                nc
            );
        }
    }
}

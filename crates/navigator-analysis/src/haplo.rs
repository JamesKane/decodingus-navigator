//! mtDNA/Y haplogroup assignment over an FTDNA haplotree using the **Kulczynski measure**
//! (HaploGrep, Weissensteiner et al.): rank each haplogroup by the set similarity between
//! its *expected* mutations (the union of branch-defining loci from root to the node) and
//! the sample's *found* polymorphisms. Higher fidelity than a flat derived/ancestral count.
//!
//! `score = ½·(|F∩E| / |E| + |F∩E| / |F|)` per node, equal site weights (a published
//! per-site weight table can be layered on later). Pure: callers supply the parsed tree
//! and the found set; fetching the FTDNA JSON lives in the app layer.
//!
//! Caveat: `found` is the sample's variants vs rCRS, so rCRS's own haplogroup-defining
//! sites aren't polymorphisms and the H backbone is under-supported (the classic
//! rCRS-vs-RSRS anchoring issue). Adequate for ranking within the rCRS-similar space; an
//! RSRS re-anchoring would be the next fidelity step.

use std::collections::HashMap;

use serde::Deserialize;

/// A branch-defining locus: a position and its ancestral/derived alleles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub name: String,
}

/// One haplotree node.
#[derive(Debug, Clone)]
pub struct HaploNode {
    pub id: i64,
    pub name: String,
    pub is_root: bool,
    pub loci: Vec<Locus>,
    pub children: Vec<i64>,
}

/// A parsed haplotree (nodes keyed by id).
#[derive(Debug, Clone)]
pub struct HaploTree {
    pub nodes: HashMap<i64, HaploNode>,
}

/// A scored candidate haplogroup, best-first after [`score`].
#[derive(Debug, Clone, PartialEq)]
pub struct ScoredHaplogroup {
    pub name: String,
    pub score: f64,
    pub depth: usize,
    /// Root→node lineage of haplogroup names.
    pub lineage: Vec<String>,
    /// Expected mutations on the path that the sample carries.
    pub matched: usize,
    /// Total expected mutations on the path.
    pub expected: usize,
    /// Total found polymorphisms in the sample.
    pub found: usize,
}

// ---- FTDNA tree JSON (subset we use) -----------------------------------------

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct FtdnaVariant {
    variant: Option<String>,
    position: Option<i64>,
    ancestral: Option<String>,
    derived: Option<String>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct FtdnaNode {
    haplogroup_id: i64,
    name: String,
    is_root: bool,
    #[serde(default)]
    variants: Vec<FtdnaVariant>,
    #[serde(default)]
    children: Option<Vec<i64>>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct FtdnaTreeJson {
    all_nodes: HashMap<String, FtdnaNode>,
}

/// Parse an FTDNA haplotree JSON document into a [`HaploTree`]. Positions are abs-valued
/// (the FTDNA data carries some negatives); variants without a position are dropped.
pub fn parse_ftdna_json(data: &str) -> Result<HaploTree, String> {
    let raw: FtdnaTreeJson = serde_json::from_str(data).map_err(|e| e.to_string())?;
    let nodes = raw
        .all_nodes
        .into_values()
        .map(|n| {
            let loci = n
                .variants
                .into_iter()
                .filter_map(|v| {
                    let pos = v.position?;
                    Some(Locus {
                        position: pos.abs(),
                        ancestral: v.ancestral.unwrap_or_default(),
                        derived: v.derived.unwrap_or_default(),
                        name: v.variant.unwrap_or_default(),
                    })
                })
                .collect();
            (n.haplogroup_id, HaploNode {
                id: n.haplogroup_id,
                name: n.name,
                is_root: n.is_root,
                loci,
                children: n.children.unwrap_or_default(),
            })
        })
        .collect();
    Ok(HaploTree { nodes })
}

/// Rank every haplogroup in `tree` against the sample's `found` polymorphisms (position →
/// uppercase base) by the Kulczynski measure. Best-first (highest score; shallower wins
/// ties — a child that adds no matched mutation shouldn't outrank its parent).
pub fn score(tree: &HaploTree, found: &HashMap<i64, char>) -> Vec<ScoredHaplogroup> {
    let mut out = Vec::new();
    let mut expected: HashMap<i64, char> = HashMap::new();
    let mut matched: usize = 0usize;
    let mut lineage: Vec<String> = Vec::new();

    let mut roots: Vec<i64> = tree.nodes.values().filter(|n| n.is_root).map(|n| n.id).collect();
    roots.sort_unstable();
    for r in roots {
        dfs(tree, r, found, &mut expected, &mut matched, 0, &mut lineage, &mut out);
    }

    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.depth.cmp(&b.depth))
    });
    out
}

#[allow(clippy::too_many_arguments)]
fn dfs(
    tree: &HaploTree,
    id: i64,
    found: &HashMap<i64, char>,
    expected: &mut HashMap<i64, char>,
    matched: &mut usize,
    depth: usize,
    lineage: &mut Vec<String>,
    out: &mut Vec<ScoredHaplogroup>,
) {
    let Some(node) = tree.nodes.get(&id) else { return };
    lineage.push(node.name.clone());

    // Add this node's loci to the path's expected set (skip recurrent positions).
    let mut added: Vec<i64> = Vec::new();
    for locus in &node.loci {
        let Some(d) = locus.derived.chars().next().map(|c| c.to_ascii_uppercase()) else { continue };
        if expected.contains_key(&locus.position) {
            continue;
        }
        expected.insert(locus.position, d);
        added.push(locus.position);
        if found.get(&locus.position).is_some_and(|b| b.eq_ignore_ascii_case(&d)) {
            *matched += 1;
        }
    }

    let (e, f) = (expected.len(), found.len());
    let kulczynski = if e > 0 && f > 0 {
        let m = *matched as f64;
        0.5 * (m / e as f64 + m / f as f64)
    } else {
        0.0
    };
    out.push(ScoredHaplogroup {
        name: node.name.clone(),
        score: kulczynski,
        depth,
        lineage: lineage.clone(),
        matched: *matched,
        expected: e,
        found: f,
    });

    let mut children = node.children.clone();
    children.sort_unstable();
    for c in children {
        dfs(tree, c, found, expected, matched, depth + 1, lineage, out);
    }

    // Backtrack.
    for p in added {
        if let Some(d) = expected.remove(&p) {
            if found.get(&p).is_some_and(|b| b.eq_ignore_ascii_case(&d)) {
                *matched -= 1;
            }
        }
    }
    lineage.pop();
}

#[cfg(test)]
mod tests {
    use super::*;

    // A tiny tree:  root --A(146)--> H --B(263)--> H2 --C(750)--> H2a
    const TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "H2", "isRoot": false,
              "variants": [{"variant":"A263G","position":263,"ancestral":"A","derived":"G"}], "children": [4]},
        "4": {"haplogroupId": 4, "name": "H2a", "isRoot": false,
              "variants": [{"variant":"C750T","position":750,"ancestral":"C","derived":"T"}], "children": []}
      }
    }"#;

    fn found(pairs: &[(i64, char)]) -> HashMap<i64, char> {
        pairs.iter().copied().collect()
    }

    #[test]
    fn parses_and_drops_positionless_variants() {
        let t = parse_ftdna_json(TREE).unwrap();
        assert_eq!(t.nodes.len(), 4);
        assert_eq!(t.nodes[&2].loci[0].position, 146);
        assert_eq!(t.nodes[&2].loci[0].derived, "G");
    }

    #[test]
    fn perfect_match_picks_the_deepest_node() {
        // sample carries all three derived alleles -> H2a is the best (matched 3 of 3).
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &found(&[(146, 'G'), (263, 'G'), (750, 'T')]));
        assert_eq!(ranked[0].name, "H2a");
        assert_eq!(ranked[0].matched, 3);
        assert_eq!(ranked[0].expected, 3);
        assert!((ranked[0].score - 1.0).abs() < 1e-9); // |F∩E|=3, |E|=3, |F|=3
    }

    #[test]
    fn partial_match_stops_at_the_supported_node() {
        // only the first two derived alleles present -> H2 wins, H2a scores lower.
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &found(&[(146, 'G'), (263, 'G')]));
        assert_eq!(ranked[0].name, "H2");
        assert!((ranked[0].score - 1.0).abs() < 1e-9); // matched 2, |E|=2, |F|=2
        let h2a = ranked.iter().find(|r| r.name == "H2a").unwrap();
        assert!(h2a.score < ranked[0].score); // H2a: matched 2, |E|=3 -> 0.5*(2/3+2/2) < 1
    }

    #[test]
    fn no_variants_yields_root() {
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &found(&[]));
        assert_eq!(ranked[0].score, 0.0);
    }
}

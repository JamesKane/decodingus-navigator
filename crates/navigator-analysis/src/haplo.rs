//! mtDNA/Y haplogroup assignment over an FTDNA haplotree using the **Kulczynski measure**
//! (HaploGrep, Weissensteiner et al.): rank each haplogroup by the set similarity between
//! its *expected* mutations (the union of branch-defining loci from root to the node) and
//! the sample's *found* polymorphisms. Higher fidelity than a flat derived/ancestral count.
//!
//! `score = ½·(|F∩E| / |E| + |F∩E| / |F|)` per node, equal site weights (a published
//! per-site weight table can be layered on later). Pure: callers supply the parsed tree
//! and the sample's base calls; fetching the FTDNA JSON lives in the app layer.
//!
//! **RSRS-anchored, reference-free.** Rather than diffing the sample against rCRS (which
//! would hide rCRS's own backbone mutations — the classic rCRS-vs-RSRS problem), we read
//! the sample's *actual base* at each tree position and compare it to the node's derived
//! allele. The FTDNA tree is RSRS-rooted, so a base equal to a node's derived allele is a
//! genuine carried mutation, backbone included — no reference subtraction needed. `found`
//! is then the set of tree sites where the sample carries the derived allele. (Assumes the
//! sample is on rCRS coordinates, i.e. ~16,569 bp; indels would shift later positions.)

use std::collections::{HashMap, HashSet};

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
    /// Tree node id (for follow-up queries like [`child_evidence`]).
    pub id: i64,
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

/// Rank every haplogroup in `tree` by the Kulczynski measure, given the sample's base at
/// each position (`calls`: 1-based position → uppercase base, from the full sequence).
/// `found` is the set of tree sites where the sample carries the derived allele; expected
/// is the root→node derived loci. Best-first (highest score; shallower wins ties — a child
/// that adds no matched mutation shouldn't outrank its parent).
pub fn score(tree: &HaploTree, calls: &HashMap<i64, char>) -> Vec<ScoredHaplogroup> {
    // |F| — distinct tree sites whose derived allele the sample carries.
    let mut carried: HashSet<i64> = HashSet::new();
    for node in tree.nodes.values() {
        for locus in &node.loci {
            if locus_carried(locus, calls) {
                carried.insert(locus.position);
            }
        }
    }
    let total_found = carried.len();

    let mut out = Vec::new();
    let mut on_path: HashSet<i64> = HashSet::new();
    let mut matched: usize = 0;
    let mut lineage: Vec<String> = Vec::new();

    let mut roots: Vec<i64> = tree.nodes.values().filter(|n| n.is_root).map(|n| n.id).collect();
    roots.sort_unstable();
    for r in roots {
        dfs(tree, r, calls, total_found, &mut on_path, &mut matched, 0, &mut lineage, &mut out);
    }

    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.depth.cmp(&b.depth))
    });
    out
}

/// Does the sample carry this locus's derived allele?
fn locus_carried(locus: &Locus, calls: &HashMap<i64, char>) -> bool {
    match locus.derived.chars().next() {
        Some(d) => calls.get(&locus.position).is_some_and(|b| b.eq_ignore_ascii_case(&d)),
        None => false,
    }
}

#[allow(clippy::too_many_arguments)]
fn dfs(
    tree: &HaploTree,
    id: i64,
    calls: &HashMap<i64, char>,
    total_found: usize,
    on_path: &mut HashSet<i64>,
    matched: &mut usize,
    depth: usize,
    lineage: &mut Vec<String>,
    out: &mut Vec<ScoredHaplogroup>,
) {
    let Some(node) = tree.nodes.get(&id) else { return };
    lineage.push(node.name.clone());

    // Add this node's loci to the path (skip positions already seen on the path).
    let mut added: Vec<(i64, bool)> = Vec::new();
    for locus in &node.loci {
        if locus.derived.is_empty() || !on_path.insert(locus.position) {
            continue;
        }
        let carried = locus_carried(locus, calls);
        if carried {
            *matched += 1;
        }
        added.push((locus.position, carried));
    }

    let expected = on_path.len();
    let kulczynski = if expected > 0 && total_found > 0 {
        let m = *matched as f64;
        0.5 * (m / expected as f64 + m / total_found as f64)
    } else {
        0.0
    };
    out.push(ScoredHaplogroup {
        id: node.id,
        name: node.name.clone(),
        score: kulczynski,
        depth,
        lineage: lineage.clone(),
        matched: *matched,
        expected,
        found: total_found,
    });

    let mut children = node.children.clone();
    children.sort_unstable();
    for c in children {
        dfs(tree, c, calls, total_found, on_path, matched, depth + 1, lineage, out);
    }

    // Backtrack.
    for (pos, carried) in added {
        on_path.remove(&pos);
        if carried {
            *matched -= 1;
        }
    }
    lineage.pop();
}

/// Every position that defines some branch → the name of a haplogroup that uses it (for
/// annotating off-path private variants). Recurrent positions keep one name.
pub fn tree_positions(tree: &HaploTree) -> HashMap<i64, String> {
    let mut m = HashMap::new();
    for n in tree.nodes.values() {
        for l in &n.loci {
            m.entry(l.position).or_insert_with(|| n.name.clone());
        }
    }
    m
}

/// The defining-SNP positions on the root→`node_id` path (the placement's backbone).
pub fn path_positions(tree: &HaploTree, node_id: i64) -> HashSet<i64> {
    let mut parent: HashMap<i64, i64> = HashMap::new();
    for n in tree.nodes.values() {
        for &c in &n.children {
            parent.insert(c, n.id);
        }
    }
    let mut positions = HashSet::new();
    let mut cur = Some(node_id);
    while let Some(id) = cur {
        match tree.nodes.get(&id) {
            Some(node) => {
                positions.extend(node.loci.iter().map(|l| l.position));
                cur = parent.get(&id).copied();
            }
            None => break,
        }
    }
    positions
}

/// The sample's state at a defining SNP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallState {
    /// Carries the branch's derived allele.
    Derived,
    /// Carries the ancestral allele (this branch's split is not supported).
    Ancestral,
    /// No confident base call at this position.
    NoCall,
}

/// One defining SNP of a branch, with the sample's state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpEvidence {
    pub name: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub state: CallState,
}

/// A child branch below the reported terminal, with the per-SNP evidence that explains why
/// descent did or didn't continue into it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchEvidence {
    pub name: String,
    pub snps: Vec<SnpEvidence>,
    /// How many of this branch's defining SNPs the sample carries (derived).
    pub derived: usize,
}

/// For the node `node_id` (typically the reported terminal), evaluate each child branch's
/// defining SNPs against the sample `calls` — `Derived` / `Ancestral` / `NoCall` per SNP.
/// This explains a stop: a child with all-`Ancestral` SNPs is an unsupported split; one
/// with `NoCall` SNPs is unresolved for lack of coverage. Children without defining SNPs
/// are omitted.
pub fn child_evidence(tree: &HaploTree, calls: &HashMap<i64, char>, node_id: i64) -> Vec<BranchEvidence> {
    let Some(node) = tree.nodes.get(&node_id) else { return Vec::new() };
    let mut out = Vec::new();
    let mut children = node.children.clone();
    children.sort_unstable();
    for cid in children {
        let Some(child) = tree.nodes.get(&cid) else { continue };
        if child.loci.is_empty() {
            continue;
        }
        let mut snps = Vec::with_capacity(child.loci.len());
        let mut derived = 0;
        for l in &child.loci {
            let d = l.derived.chars().next().map(|c| c.to_ascii_uppercase());
            let a = l.ancestral.chars().next().map(|c| c.to_ascii_uppercase());
            let state = match calls.get(&l.position).map(|c| c.to_ascii_uppercase()) {
                Some(b) if Some(b) == d => CallState::Derived,
                Some(b) if Some(b) == a => CallState::Ancestral,
                Some(_) => CallState::Ancestral, // a third allele — not this branch's derived
                None => CallState::NoCall,
            };
            if state == CallState::Derived {
                derived += 1;
            }
            snps.push(SnpEvidence {
                name: l.name.clone(),
                position: l.position,
                ancestral: l.ancestral.clone(),
                derived: l.derived.clone(),
                state,
            });
        }
        out.push(BranchEvidence { name: child.name.clone(), snps, derived });
    }
    out
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

    /// Sample base calls by position (the bases the sample carries at these positions).
    fn calls(pairs: &[(i64, char)]) -> HashMap<i64, char> {
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
        let ranked = score(&t, &calls(&[(146, 'G'), (263, 'G'), (750, 'T')]));
        assert_eq!(ranked[0].name, "H2a");
        assert_eq!(ranked[0].matched, 3);
        assert_eq!(ranked[0].expected, 3);
        assert!((ranked[0].score - 1.0).abs() < 1e-9); // |F∩E|=3, |E|=3, |F|=3
    }

    #[test]
    fn partial_match_stops_at_the_supported_node() {
        // only the first two derived alleles present -> H2 wins, H2a scores lower.
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &calls(&[(146, 'G'), (263, 'G')]));
        assert_eq!(ranked[0].name, "H2");
        assert!((ranked[0].score - 1.0).abs() < 1e-9); // matched 2, |E|=2, |F|=2
        let h2a = ranked.iter().find(|r| r.name == "H2a").unwrap();
        assert!(h2a.score < ranked[0].score); // H2a: matched 2, |E|=3 -> 0.5*(2/3+2/2) < 1
    }

    #[test]
    fn child_evidence_explains_an_unsupported_split() {
        // H2 has a child H2a (derived T@750). Sample is ancestral (C) at 750 -> the split
        // into H2a is not supported, shown per-SNP.
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &calls(&[(146, 'G'), (263, 'G'), (750, 'C')]));
        assert_eq!(ranked[0].name, "H2"); // stops at H2 (750 ancestral)
        let ev = child_evidence(&t, &calls(&[(146, 'G'), (263, 'G'), (750, 'C')]), ranked[0].id);
        assert_eq!(ev.len(), 1);
        assert_eq!(ev[0].name, "H2a");
        assert_eq!(ev[0].derived, 0);
        assert_eq!(ev[0].snps[0].position, 750);
        assert_eq!(ev[0].snps[0].state, CallState::Ancestral);
    }

    #[test]
    fn no_variants_yields_root() {
        let t = parse_ftdna_json(TREE).unwrap();
        let ranked = score(&t, &calls(&[]));
        assert_eq!(ranked[0].score, 0.0);
    }
}

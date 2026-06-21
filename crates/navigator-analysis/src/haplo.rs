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
            (
                n.haplogroup_id,
                HaploNode {
                    id: n.haplogroup_id,
                    name: n.name,
                    is_root: n.is_root,
                    loci,
                    children: n.children.unwrap_or_default(),
                },
            )
        })
        .collect();
    Ok(HaploTree { nodes })
}

// ---- DecodingUs tree JSON (the AppView `/api/v1/y-tree/full` shape) -----------

#[derive(Deserialize)]
struct DuCoord {
    position: i64,
    #[serde(default)]
    ancestral: Option<String>,
    #[serde(default)]
    derived: Option<String>,
}

#[derive(Deserialize)]
struct DuVariant {
    #[serde(default)]
    canonical_name: String,
    /// Coordinates keyed by build label (`"hs1"`, `"GRCh38"`, `"GRCh37"`).
    #[serde(default)]
    coordinates: HashMap<String, DuCoord>,
}

#[derive(Deserialize)]
struct DuNode {
    id: i64,
    name: String,
    #[serde(default)]
    variants: Vec<DuVariant>,
    #[serde(default)]
    children: Vec<DuNode>,
}

#[derive(Deserialize)]
struct DuTreeJson {
    roots: Vec<DuNode>,
}

/// Parse the DecodingUs AppView Y-tree (`/api/v1/y-tree/full`) into a [`HaploTree`], taking
/// each variant's coordinate for `build_key` (`"hs1"` for CHM13, `"GRCh38"`, `"GRCh37"`).
/// Because positions are read in the *alignment's own build*, no liftover is needed —
/// variants without a coordinate on `build_key` are dropped (they can't be placed there).
/// Node ids come from the AppView (unique); the nested `children` flatten into child-id lists.
pub fn parse_decodingus_json(data: &str, build_key: &str) -> Result<HaploTree, String> {
    let raw: DuTreeJson = serde_json::from_str(data).map_err(|e| e.to_string())?;
    let mut nodes = HashMap::new();
    for root in &raw.roots {
        flatten_du_node(root, true, build_key, &mut nodes);
    }
    Ok(HaploTree { nodes })
}

fn flatten_du_node(n: &DuNode, is_root: bool, build_key: &str, out: &mut HashMap<i64, HaploNode>) {
    let loci = n
        .variants
        .iter()
        .filter_map(|v| {
            let c = v.coordinates.get(build_key)?;
            Some(Locus {
                position: c.position.abs(),
                ancestral: c.ancestral.clone().unwrap_or_default(),
                derived: c.derived.clone().unwrap_or_default(),
                name: v.canonical_name.clone(),
            })
        })
        .collect();
    let children = n.children.iter().map(|c| c.id).collect();
    out.insert(
        n.id,
        HaploNode {
            id: n.id,
            name: n.name.clone(),
            is_root,
            loci,
            children,
        },
    );
    for c in &n.children {
        flatten_du_node(c, false, build_key, out);
    }
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
        dfs(
            tree,
            r,
            calls,
            total_found,
            &mut on_path,
            &mut matched,
            0,
            &mut lineage,
            &mut out,
        );
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

/// The sample's state at one defining SNP: carries the derived allele, carries the
/// ancestral allele (or any non-derived base — a contradiction of the branch), or has no
/// confident call. A locus with no derived allele (indel-only / marker-less) is `NoCall`.
fn locus_state(locus: &Locus, calls: &HashMap<i64, char>) -> CallState {
    let d = locus.derived.chars().next().map(|c| c.to_ascii_uppercase());
    let a = locus.ancestral.chars().next().map(|c| c.to_ascii_uppercase());
    if d.is_none() {
        return CallState::NoCall;
    }
    match calls.get(&locus.position).map(|c| c.to_ascii_uppercase()) {
        Some(b) if Some(b) == d => CallState::Derived,
        Some(b) if Some(b) == a => CallState::Ancestral,
        Some(_) => CallState::Ancestral, // a third allele — not this branch's derived
        None => CallState::NoCall,
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

/// child → parent map, used to walk any node back to the root.
fn build_parent_map(tree: &HaploTree) -> HashMap<i64, i64> {
    let mut parent: HashMap<i64, i64> = HashMap::new();
    for n in tree.nodes.values() {
        for &c in &n.children {
            parent.insert(c, n.id);
        }
    }
    parent
}

/// The defining-SNP positions on the root→`node_id` path (the placement's backbone).
pub fn path_positions(tree: &HaploTree, node_id: i64) -> HashSet<i64> {
    let parent = build_parent_map(tree);
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
    let Some(node) = tree.nodes.get(&node_id) else {
        return Vec::new();
    };
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
            let state = locus_state(l, calls);
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
        out.push(BranchEvidence {
            name: child.name.clone(),
            snps,
            derived,
        });
    }
    out
}

/// Per-SNP evidence along the lineage root→`terminal_id`: every defining SNP of every node on
/// the path, with the sample's `Derived`/`Ancestral`/`NoCall` state. Used to compare exactly
/// which defining mutations a sample carries (e.g. GRCh38 vs a lifted CHM13 call).
pub fn lineage_evidence(tree: &HaploTree, calls: &HashMap<i64, char>, terminal_id: i64) -> Vec<SnpEvidence> {
    // child → parent, to walk the terminal back to the root.
    let parent = build_parent_map(tree);
    let mut path = Vec::new();
    let mut cur = Some(terminal_id);
    while let Some(id) = cur {
        path.push(id);
        cur = parent.get(&id).copied();
    }
    path.reverse();

    let mut out = Vec::new();
    for id in path {
        let Some(node) = tree.nodes.get(&id) else { continue };
        for l in &node.loci {
            let state = locus_state(l, calls);
            out.push(SnpEvidence {
                name: l.name.clone(),
                position: l.position,
                ancestral: l.ancestral.clone(),
                derived: l.derived.clone(),
                state,
            });
        }
    }
    out
}

// ---- path-supported parsimony guard ------------------------------------------
//
// The Kulczynski `score` ranks every node by *proportional* set-similarity, and on real
// data that places the terminal well (validated: GFX0457637 → R-FGC29071). Its one weakness
// is the distal-Y paralog artifact: a deep node reached only by *tunnelling through a branch
// the sample contradicts* can still score highly off a few coincidental matches.
//
// Parsimony guards exactly that failure — it rejects any candidate whose root→node lineage
// crosses a contradicted branch — without disturbing the proportional ranking that gets the
// clean case right. (A descent-style "follow the most-derived subtree" router was tried and
// derailed onto a bushier wrong fork on the 4× GFX sample: absolute derived count favours
// long/bushy paths, where Kulczynski's proportion does not. The proportional rank + this
// guard is the validated combination. The remaining paralog *false-positive* defence — when
// the wrong branch carries spurious derived calls rather than honest ancestral ones — is the
// haploid allele-balance filter, a separate Phase-1 item.) See PangenomeExpansion.md.

/// Per-node tally over evaluable defining SNPs (loci with a derived allele): how many the
/// sample calls derived, ancestral (a contradiction), or has no confident base for.
fn node_counts(node: &HaploNode, calls: &HashMap<i64, char>) -> (usize, usize, usize) {
    let (mut d, mut a, mut n) = (0usize, 0usize, 0usize);
    for l in &node.loci {
        if l.derived.is_empty() {
            continue; // marker-less locus — not evaluable by a SNP caller
        }
        match locus_state(l, calls) {
            CallState::Derived => d += 1,
            CallState::Ancestral => a += 1,
            CallState::NoCall => n += 1,
        }
    }
    (d, a, n)
}

/// A node is *contradicted* when the sample carries the ancestral allele at more of the
/// node's defining SNPs than the derived allele — it confidently does **not** belong to this
/// branch. A no-evidence node (all no-call, `d == a == 0`) is *not* contradicted: it is a
/// pass-through, so low coverage never blocks a lineage. A stray ancestral at an otherwise
/// well-supported node (`d >= a`) is tolerated for the same reason.
fn is_contradicted(node: &HaploNode, calls: &HashMap<i64, char>) -> bool {
    let (d, a, _) = node_counts(node, calls);
    a > d
}

/// Derived defining-SNPs that must appear *below* a contradicted ancestor (on the path toward the
/// candidate) before that contradiction stops vetoing the lineage. A lone ancestral at a sparse
/// intermediate node is usually a genotyping artifact — common on targeted-Y data (FTDNA Big Y),
/// whose coverage gaps turn most intermediate SNPs into no-calls — and would otherwise veto an
/// entire deep lineage that the terminal overwhelmingly supports (one stray ancestral at R-Z16250
/// blocking R-CTS4466 with 10 derived / 0 ancestral). A real off-branch *tunnel* artifact, by
/// contrast, carries only a coincidental hit or two below the contradicted branch-point, so it
/// stays vetoed. The threshold sits above the coincidental noise and well below a genuine clade's
/// derived count.
const REDEEM_DERIVED: usize = 4;

/// Is the root→`node_id` lineage free of any *unredeemed* contradicted branch? An off-path paralog
/// artifact sits below a branch the sample is ancestral for, so it fails this guard; the genuine
/// lineage (derived or merely no-call along its length) passes. Used to veto otherwise high-scoring
/// tunnel artifacts from the [`score`] ranking.
///
/// A contradicted ancestor only vetoes when it isn't *redeemed* by derived support further down the
/// path: a single stray ancestral at a sparse intermediate node (a Big Y miscall) is overridden
/// when ≥[`REDEEM_DERIVED`] derived SNPs below it confirm the branch, while a coincidental tunnel —
/// a contradicted branch-point with only a hit or two beneath it — stays vetoed.
pub fn path_admissible(tree: &HaploTree, calls: &HashMap<i64, char>, node_id: i64) -> bool {
    let parent = build_parent_map(tree);
    // Root→node path (root first), with each node's derived-call count.
    let mut path: Vec<i64> = Vec::new();
    let mut cur = Some(node_id);
    while let Some(id) = cur {
        path.push(id);
        cur = parent.get(&id).copied();
    }
    path.reverse();
    let derived: Vec<usize> = path
        .iter()
        .map(|id| tree.nodes.get(id).map_or(0, |n| node_counts(n, calls).0))
        .collect();
    for (i, id) in path.iter().enumerate() {
        let Some(node) = tree.nodes.get(id) else { continue };
        if is_contradicted(node, calls) {
            // Derived SNPs strictly below this contradicted node, toward the candidate.
            let derived_below: usize = derived[i + 1..].iter().sum();
            if derived_below < REDEEM_DERIVED {
                return false;
            }
        }
    }
    true
}

/// Minimum derived defining-SNPs of a child the sample must carry for [`deepen_terminal`] to
/// descend into it. Two independent shared-derived mutations confirm membership while staying
/// robust to a lone recurrent/artefactual match.
const MIN_DERIVED_TO_DEEPEN: usize = 2;

/// From the guard-selected `start`, descend further into any child the sample has *clearly
/// entered* — carries at least [`MIN_DERIVED_TO_DEEPEN`] of its derived SNPs and is not
/// contradicted (`ancestral ≤ derived`). Routes by derived count (ties → lower id).
///
/// This corrects under-calling at **unsplit tree nodes** — a common case in published trees
/// (FTDNA especially): when a node's SNP block has not yet been divided into sub-branches, a
/// sample on one sub-lineage is derived for the SNPs defining its own line and ancestral for
/// the SNPs of the *other* (not-yet-split) sub-lineages. The node then looks "half ancestral",
/// so its proportional [`score`] falls just below its parent's and the guard stops one node
/// too shallow — even though the sample genuinely carries several of the node's mutations. The
/// ancestral SNPs are an unresolved downstream split, not a contradiction.
pub fn deepen_terminal(tree: &HaploTree, calls: &HashMap<i64, char>, start: i64) -> i64 {
    let mut current = start;
    while let Some(node) = tree.nodes.get(&current) {
        let mut children = node.children.clone();
        children.sort_unstable();
        let mut best: Option<(i64, usize)> = None; // (child id, derived count)
        for cid in children {
            let Some(child) = tree.nodes.get(&cid) else { continue };
            if is_contradicted(child, calls) {
                continue; // the sample is net-ancestral here — not below this branch
            }
            let (d, _, _) = node_counts(child, calls);
            if d < MIN_DERIVED_TO_DEEPEN {
                continue;
            }
            if best.map_or(true, |(_, bd)| d > bd) {
                best = Some((cid, d));
            }
        }
        match best {
            Some((cid, _)) => current = cid,
            None => break,
        }
    }
    current
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

    // The DecodingUs AppView `/api/v1/y-tree/full` shape (snake_case, nested children,
    // multi-build coordinates). R-M207 → R1 with an `hs1` (CHM13) and a GRCh38 coordinate.
    const DU_TREE: &str = r#"{
      "roots": [
        {"id": 10, "name": "R-M207", "haplogroup_type": "Y_DNA", "variants": [
            {"canonical_name": "M207", "coordinates": {
                "hs1": {"contig":"chrY","position":2800000,"ancestral":"A","derived":"G"},
                "GRCh38": {"contig":"chrY","position":2900000,"ancestral":"A","derived":"G"}}}],
         "children": [
            {"id": 11, "name": "R-M173", "haplogroup_type": "Y_DNA", "variants": [
                {"canonical_name": "M173", "coordinates": {
                    "hs1": {"contig":"chrY","position":2810000,"ancestral":"C","derived":"T"}}},
                {"canonical_name": "GRCh38only", "coordinates": {
                    "GRCh38": {"contig":"chrY","position":2999999,"ancestral":"G","derived":"A"}}}],
             "children": []}]}
      ]
    }"#;

    #[test]
    fn parse_decodingus_picks_target_build_and_flattens() {
        // hs1: both M207 and M173 resolve; the GRCh38-only variant is dropped.
        let t = parse_decodingus_json(DU_TREE, "hs1").unwrap();
        assert_eq!(t.nodes.len(), 2);
        assert!(t.nodes[&10].is_root && !t.nodes[&11].is_root);
        assert_eq!(t.nodes[&10].children, vec![11]);
        assert_eq!(t.nodes[&10].loci[0].position, 2800000);
        assert_eq!(t.nodes[&10].loci[0].name, "M207");
        // R-M173 keeps only the hs1-coordinated M173 (GRCh38only dropped).
        assert_eq!(t.nodes[&11].loci.len(), 1);
        assert_eq!(t.nodes[&11].loci[0].position, 2810000);

        // GRCh38: M207 uses its GRCh38 position; M173's hs1-only locus drops, GRCh38only stays.
        let g = parse_decodingus_json(DU_TREE, "GRCh38").unwrap();
        assert_eq!(g.nodes[&10].loci[0].position, 2900000);
        let names: Vec<&str> = g.nodes[&11].loci.iter().map(|l| l.name.as_str()).collect();
        assert_eq!(names, vec!["GRCh38only"]);
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

    // ---- parsimony admissibility guard ----

    /// The terminal `assemble_assignment` reports: the best-ranked candidate whose lineage
    /// the parsimony guard admits (mirrors the app-layer integration).
    fn guarded_terminal(t: &HaploTree, c: &HashMap<i64, char>) -> String {
        let ranked = score(t, c);
        ranked
            .iter()
            .find(|r| path_admissible(t, c, r.id))
            .map(|r| r.name.clone())
            .unwrap()
    }

    fn id_of(t: &HaploTree, name: &str) -> i64 {
        t.nodes.values().find(|n| n.name == name).unwrap().id
    }

    #[test]
    fn reference_polarity_comes_from_the_tree_not_the_reference() {
        // The CHM13 trap: at a Y-SNP the tree calls ancestral=A, derived=G, and the analysis
        // reference (CHM13 chrY = HG002, haplogroup J) carries the DERIVED base G. Polarity
        // must come from comparing the SAMPLE's base to the tree — never from the reference.
        let locus = Locus {
            position: 146,
            ancestral: "A".into(),
            derived: "G".into(),
            name: "M-test".into(),
        };

        // Sample carries the ANCESTRAL allele (A). A "reference base = ancestral" assumption
        // (ref here is G) would see A ≠ G and wrongly flip this to Derived. Tree-driven: Ancestral.
        assert_eq!(locus_state(&locus, &calls(&[(146, 'A')])), CallState::Ancestral);
        // Sample carries the DERIVED allele (G, == the reference here): Derived, from the tree.
        assert_eq!(locus_state(&locus, &calls(&[(146, 'G')])), CallState::Derived);

        // End-to-end: a sample ANCESTRAL at the J-derived backbone site 146 does not carry H's
        // defining mutation, so it must not be placed into H — even though the CHM13 reference
        // base there is the derived G. A REF-as-ancestral assumption would flip 146 and wrongly
        // descend; tree-driven, H is contradicted and the call stays at root.
        let t = parse_ftdna_json(TREE).unwrap(); // root→H(146 A→G)→H2(263)→H2a(750)
        let c = calls(&[(146, 'A'), (263, 'G'), (750, 'T')]);
        assert!(
            !path_admissible(&t, &c, id_of(&t, "H")),
            "H is contradicted (sample ancestral at 146)"
        );
        assert_eq!(guarded_terminal(&t, &c), "root");
    }

    #[test]
    fn guard_admits_a_clean_lineage() {
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[(146, 'G'), (263, 'G'), (750, 'T')]);
        assert!(path_admissible(&t, &c, id_of(&t, "H2a")));
        assert_eq!(guarded_terminal(&t, &c), "H2a");
    }

    #[test]
    fn guard_rejects_a_contradicted_terminal_but_admits_its_parent() {
        // Ancestral (C) at 750 -> H2a is contradicted; the report falls back to H2.
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[(146, 'G'), (263, 'G'), (750, 'C')]);
        assert!(!path_admissible(&t, &c, id_of(&t, "H2a")));
        assert!(path_admissible(&t, &c, id_of(&t, "H2")));
        assert_eq!(guarded_terminal(&t, &c), "H2");
    }

    #[test]
    fn no_calls_admit_the_whole_tree() {
        // Empty calls: nothing is contradicted, so every lineage is admissible (the guard is
        // a veto, not a selector — Kulczynski still picks root for lack of matches).
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[]);
        assert!(path_admissible(&t, &c, id_of(&t, "H2a")));
        assert_eq!(guarded_terminal(&t, &c), "root");
    }

    // root -> H(146) -> B(500, contradicted) -> Bdeep(900, coincidental derived).
    // Kulczynski is lured to Bdeep (matches 146 + 900); the guard must veto it (tunnels
    // through the contradicted B) and fall back to H.
    const TUNNEL_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "B", "isRoot": false,
              "variants": [{"variant":"C500T","position":500,"ancestral":"C","derived":"T"}], "children": [4]},
        "4": {"haplogroupId": 4, "name": "Bdeep", "isRoot": false,
              "variants": [{"variant":"G900A","position":900,"ancestral":"G","derived":"A"}], "children": []}
      }
    }"#;

    #[test]
    fn guard_vetoes_the_tunnel_artifact() {
        let t = parse_ftdna_json(TUNNEL_TREE).unwrap();
        // Carries 146 (H) and a coincidental 900 (Bdeep) but is ANCESTRAL (C) at 500.
        let c = calls(&[(146, 'G'), (500, 'C'), (900, 'A')]);
        // Kulczynski alone is lured deeper by the coincidental match...
        assert_eq!(score(&t, &c)[0].name, "Bdeep");
        // ...but Bdeep tunnels through the contradicted B, so the guard reports H.
        assert!(!path_admissible(&t, &c, id_of(&t, "Bdeep")));
        assert_eq!(guarded_terminal(&t, &c), "H");
    }

    // root -> H(146) -> M(marker-less / no SNPs) -> D(263). The guard must pass through M.
    const MARKERLESS_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "M", "isRoot": false, "variants": [], "children": [4]},
        "4": {"haplogroupId": 4, "name": "D", "isRoot": false,
              "variants": [{"variant":"A263G","position":263,"ancestral":"A","derived":"G"}], "children": []}
      }
    }"#;

    #[test]
    fn guard_passes_through_marker_less_and_no_call_nodes() {
        let t = parse_ftdna_json(MARKERLESS_TREE).unwrap();
        // 146 + 263 derived: D is admissible through the marker-less M, and is the call.
        let full = calls(&[(146, 'G'), (263, 'G')]);
        assert!(path_admissible(&t, &full, id_of(&t, "D")));
        assert_eq!(guarded_terminal(&t, &full), "D");
        // 263 no-call (low coverage): D is *still* admissible (a no-call is not a
        // contradiction) — the guard never blocks for lack of coverage. Kulczynski stops at H.
        let sparse = calls(&[(146, 'G')]);
        assert!(path_admissible(&t, &sparse, id_of(&t, "D")));
        assert_eq!(guarded_terminal(&t, &sparse), "H");
    }

    // root -> H(146) -> D with three defining SNPs, used to exercise the net contradiction rule.
    const NET_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "D", "isRoot": false, "variants": [
                {"variant":"A263G","position":263,"ancestral":"A","derived":"G"},
                {"variant":"A600G","position":600,"ancestral":"A","derived":"G"},
                {"variant":"C500T","position":500,"ancestral":"C","derived":"T"}
              ], "children": []}
      }
    }"#;

    // root → P (2 derived SNPs) → C, an UNSPLIT node: 3 SNPs that define it + 3 SNPs of a
    // not-yet-split sub-branch. A sample on C's trunk is derived for the first 3 and ancestral
    // for the other 3 — so C looks half-ancestral and Kulczynski can rank it below P.
    const UNSPLIT_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "P", "isRoot": false, "variants": [
                {"variant":"P1","position":100,"ancestral":"A","derived":"G"},
                {"variant":"P2","position":200,"ancestral":"A","derived":"G"}
              ], "children": [3]},
        "3": {"haplogroupId": 3, "name": "C", "isRoot": false, "variants": [
                {"variant":"C1","position":300,"ancestral":"A","derived":"G"},
                {"variant":"C2","position":400,"ancestral":"A","derived":"G"},
                {"variant":"C3","position":500,"ancestral":"A","derived":"G"},
                {"variant":"M1","position":600,"ancestral":"A","derived":"G"},
                {"variant":"M2","position":700,"ancestral":"A","derived":"G"},
                {"variant":"M3","position":800,"ancestral":"A","derived":"G"}
              ], "children": []}
      }
    }"#;

    #[test]
    fn deepen_enters_an_unsplit_node_the_sample_clearly_carries() {
        let t = parse_ftdna_json(UNSPLIT_TREE).unwrap();
        // Derived for P (100,200) and 3 of C's SNPs (300,400,500); ancestral for the other 3.
        let c = calls(&[
            (100, 'G'),
            (200, 'G'),
            (300, 'G'),
            (400, 'G'),
            (500, 'G'),
            (600, 'A'),
            (700, 'A'),
            (800, 'A'),
        ]);
        // Deepen enters C from P: it carries 3 derived (≥2) and isn't contradicted (3 anc ≤ 3 der).
        // (The "Kulczynski stops at the parent" condition needs a long backbone — validated on
        // the real WGS229 short-read sample, where the guard stops at R-FGC29067 and deepen
        // recovers R-FGC29071.)
        assert_eq!(deepen_terminal(&t, &c, id_of(&t, "P")), id_of(&t, "C"));
    }

    #[test]
    fn deepen_does_not_enter_on_a_lone_match_or_a_net_ancestral_child() {
        let t = parse_ftdna_json(UNSPLIT_TREE).unwrap();
        // Only one of C's SNPs derived → below the ≥2 threshold (and net-ancestral): stay at P.
        let lone = calls(&[
            (100, 'G'),
            (200, 'G'),
            (300, 'G'),
            (400, 'A'),
            (500, 'A'),
            (600, 'A'),
            (700, 'A'),
            (800, 'A'),
        ]);
        assert_eq!(deepen_terminal(&t, &lone, id_of(&t, "P")), id_of(&t, "P"));
        // 2 derived but 4 ancestral → contradicted (a > d), don't enter even at ≥2 derived.
        let net_anc = calls(&[
            (100, 'G'),
            (200, 'G'),
            (300, 'G'),
            (400, 'G'),
            (500, 'A'),
            (600, 'A'),
            (700, 'A'),
            (800, 'A'),
        ]);
        assert_eq!(deepen_terminal(&t, &net_anc, id_of(&t, "P")), id_of(&t, "P"));
    }

    #[test]
    fn deepen_is_a_no_op_at_a_true_terminal() {
        // On the clean linear TREE, a perfect sample already reaches H2a; deepen adds nothing.
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[(146, 'G'), (263, 'G'), (750, 'T')]);
        assert_eq!(deepen_terminal(&t, &c, id_of(&t, "H2a")), id_of(&t, "H2a"));
    }

    #[test]
    fn guard_tolerates_a_stray_contradiction_but_blocks_a_net_one() {
        let t = parse_ftdna_json(NET_TREE).unwrap();
        // d=2 (263,600), a=1 (500 ancestral): derived outweighs -> D admitted (stray error).
        let tolerated = calls(&[(146, 'G'), (263, 'G'), (600, 'G'), (500, 'C')]);
        assert!(path_admissible(&t, &tolerated, id_of(&t, "D")));
        // d=1 (263), a=2 (600,500 ancestral): contradictions dominate -> D blocked.
        let blocked = calls(&[(146, 'G'), (263, 'G'), (600, 'A'), (500, 'C')]);
        assert!(!path_admissible(&t, &blocked, id_of(&t, "D")));
    }
}

//! mtDNA and Y haplogroup assignment over an FTDNA haplotree, with the **Kulczynski measure**
//! (HaploGrep, Weissensteiner and others). It ranks each haplogroup by the set similarity
//! between two sets. The first set is the *expected* mutations of that haplogroup. That is the
//! union of the loci that define a branch, from the root down to the node. The second set is the
//! *found* polymorphisms of the sample. This measure is more accurate than a flat count of
//! derived and ancestral sites.
//!
//! The score at each node is `score = ½·(|F∩E| / |E| + |F∩E| / |F|)`, and the site weights are
//! equal. A published table of weights for each site can go on top of this later.
//!
//! This module is pure. The caller gives the parsed tree and the base calls of the sample. The
//! app layer gets the FTDNA JSON.
//!
//! **The RSRS is the anchor, and the code uses no reference.** It does not take the difference
//! of the sample against rCRS, because that would hide the backbone mutations of rCRS itself.
//! That is the classic rCRS-against-RSRS problem. It instead reads the *actual base* of the
//! sample at each tree position, and compares it to the derived allele of the node.
//!
//! The FTDNA tree has RSRS at its root. A base that equals the derived allele of a node is a
//! true mutation of the sample, and the backbone is part of that. The code needs no subtraction
//! of a reference.
//!
//! `found` is then the set of tree sites where the sample carries the derived allele. This
//! assumes that the sample is on rCRS coordinates, which is about 16,569 bp. An indel would move
//! the positions after it.

use std::collections::{HashMap, HashSet};

use serde::Deserialize;

/// A locus that defines a branch: a position and its ancestral and derived alleles.
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

/// Parse an FTDNA haplotree JSON document into a [`HaploTree`]. The code takes the absolute
/// value of each position, because the FTDNA data carries some negative ones. It drops a variant
/// that has no position.
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
    /// The **authoritative** phylogenetic polarity at the variant level. It does not depend on
    /// the build, and every variant carries it. About 1.4% of the variants carry a *swapped*
    /// `coordinates.ancestral/derived` in a build, which is a clean exchange of the two roles.
    /// Trust this field, and not the coordinate alleles.
    #[serde(default)]
    link_ancestral: Option<String>,
    #[serde(default)]
    link_derived: Option<String>,
}

impl DuVariant {
    /// The trustworthy `(ancestral, derived)` alleles from `link_*`, or `None` when absent/empty.
    fn link_alleles(&self) -> Option<(String, String)> {
        match (self.link_ancestral.as_deref(), self.link_derived.as_deref()) {
            (Some(a), Some(d)) if !a.is_empty() && !d.is_empty() => Some((a.to_string(), d.to_string())),
            _ => None,
        }
    }

    /// Polarity for a locus: `link_*` (authoritative) with the build coordinate's alleles as a
    /// last-resort fallback if `link_*` is somehow missing.
    fn polarity(&self, coord: &DuCoord) -> (String, String) {
        if let Some(p) = self.link_alleles() {
            return p;
        }
        (
            coord.ancestral.clone().unwrap_or_default(),
            coord.derived.clone().unwrap_or_default(),
        )
    }
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

/// Parse the DecodingUs AppView Y-tree (`/api/v1/y-tree/full`) into a [`HaploTree`]. It takes the
/// coordinate of each variant for `build_key`, which is `"hs1"` for CHM13, `"GRCh38"` or
/// `"GRCh37"`.
///
/// The code reads the positions in the *build of the alignment itself*, so it needs no liftover.
/// It drops a variant that has no coordinate on `build_key`, because it can not place that
/// variant there. The node ids come from the AppView and are unique. The nested `children`
/// flatten into lists of child ids.
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
            // The position and the contig come from the build coordinate. The polarity comes
            // from the authoritative `link_*`. The ancestral and derived alleles of the
            // coordinate itself are in the wrong order on about 1.4% of the variants.
            let (ancestral, derived) = v.polarity(c);
            Some(Locus {
                position: c.position.abs(),
                ancestral,
                derived,
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

/// The **ancestral and derived polarity** of each SNP **name**, from the DecodingUs tree JSON,
/// as `name → (ancestral, derived)`. It does not depend on the build.
///
/// The DecodingUs tree carries the true phylogenetic polarity. FTDNA instead records the GRCh38
/// *reference* base as the "ancestral" one. At a site where the reference carries the derived
/// allele, the FTDNA polarity is then the wrong way round. [`normalize_polarity`] uses this map
/// to repair an FTDNA tree.
///
/// The code takes the alleles from the coordinate of any one build, because the alleles of a
/// mutation are the same in every build. The names are universal.
pub fn decodingus_polarity_map(data: &str) -> Result<HashMap<String, (String, String)>, String> {
    let raw: DuTreeJson = serde_json::from_str(data).map_err(|e| e.to_string())?;
    let mut out = HashMap::new();
    // The `link_*` at the variant level is the authoritative polarity, and it does not depend on
    // the build. Use it directly. Fall back to a build coordinate only if `link_*` is absent.
    // Choose that coordinate in a **deterministic** way: hs1 first, then the keys in sorted
    // order. A walk over `coordinates.values()` follows the HashMap order, which is not
    // deterministic. Where a build records a swapped polarity, such a walk would take a different
    // orientation on each run.
    fn pick_polarity(v: &DuVariant) -> Option<(String, String)> {
        if let Some(p) = v.link_alleles() {
            return Some(p);
        }
        let alleles = |c: &DuCoord| match (c.ancestral.as_deref(), c.derived.as_deref()) {
            (Some(a), Some(d)) if !a.is_empty() && !d.is_empty() => Some((a.to_string(), d.to_string())),
            _ => None,
        };
        let coords = &v.coordinates;
        for b in ["hs1", "GRCh38", "GRCh37"] {
            if let Some(p) = coords.get(b).and_then(alleles) {
                return Some(p);
            }
        }
        let mut keys: Vec<&String> = coords.keys().collect();
        keys.sort_unstable();
        keys.into_iter().find_map(|k| coords.get(k).and_then(alleles))
    }
    fn walk(n: &DuNode, out: &mut HashMap<String, (String, String)>) {
        for v in &n.variants {
            if v.canonical_name.is_empty() {
                continue;
            }
            if let Some(p) = pick_polarity(v) {
                out.entry(v.canonical_name.clone()).or_insert(p);
            }
        }
        for c in &n.children {
            walk(c, out);
        }
    }
    for r in &raw.roots {
        walk(r, &mut out);
    }
    Ok(out)
}

/// Repair the polarity of an FTDNA tree in place, against a `reference` polarity map from
/// [`decodingus_polarity_map`].
///
/// Take any locus whose alleles are the **same two nucleotides** as those of the reference, but
/// in the opposite roles. That is the FTDNA inversion that makes the reference the ancestral
/// allele. Put the ancestral and derived alleles of that locus back into the orientation of the
/// reference.
///
/// A locus on the other strand, which holds different nucleotides, does not change. The lift path
/// already takes the complement of those. Returns the count of the loci that changed.
pub fn normalize_polarity(tree: &mut HaploTree, reference: &HashMap<String, (String, String)>) -> usize {
    let mut flipped = 0;
    for node in tree.nodes.values_mut() {
        for l in &mut node.loci {
            if l.name.is_empty() {
                continue;
            }
            let Some((ra, rd)) = reference.get(&l.name) else {
                continue;
            };
            let la = l.ancestral.to_ascii_uppercase();
            let ld = l.derived.to_ascii_uppercase();
            let ra = ra.to_ascii_uppercase();
            let rd = rd.to_ascii_uppercase();
            // Same two nucleotides, opposite polarity → flip to the reference orientation.
            if la != ld && la == rd && ld == ra {
                std::mem::swap(&mut l.ancestral, &mut l.derived);
                flipped += 1;
            }
        }
    }
    flipped
}

/// Rank every haplogroup in `tree` by the Kulczynski measure. `calls` gives the base of the
/// sample at each position, as a 1-based position to an uppercase base, from the full sequence.
///
/// `found` is the set of tree sites where the sample carries the derived allele. The expected set
/// is the derived loci from the root down to the node. The result comes back best first, which
/// is the highest score. A node that is nearer to the root wins a tie, because a child that adds
/// no matched mutation must not rank above its parent.
pub fn score(tree: &HaploTree, calls: &HashMap<i64, char>) -> Vec<ScoredHaplogroup> {
    // |F|, the count of distinct tree sites whose derived allele the sample carries.
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

/// Watson-Crick complement of a single base (non-ACGT passes through unchanged).
fn complement_base(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// True when the two alleles of a SNP are ambiguous about the strand, which holds for an `A↔T`
/// or a `C↔G` transversion. The complement of one allele equals the other, so the observed base
/// does not tell you the strand. The code must not try a match on the complement for these.
fn strand_ambiguous(a: char, d: char) -> bool {
    let mut pair = [a.to_ascii_uppercase(), d.to_ascii_uppercase()];
    pair.sort_unstable();
    pair == ['A', 'T'] || pair == ['C', 'G']
}

/// True when a locus is a single-base SNP that a genotype at the base level can read. That is,
/// it is NOT an indel and NOT an MNP, which carry an allele of more than one character.
///
/// An insertion (`G`→`GAGC…`) or a deletion (`GAGC`→`G`) holds the same anchor base in both of
/// its alleles. A comparison of a base against an allele can then not separate the two, and it
/// reads every sample as derived. The DecodingUs Y tree carries about 12.7k such indel loci, and
/// all of them share an anchor. To count them by base turned a node with many indels into a
/// collector of homoplasy, and such a node then took the placement.
///
/// [`crate::caller::call_indels_at`] genotypes the indel loci instead. It writes a resolved
/// **sentinel** into the genotype map, which is [`INDEL_DERIVED`] or [`INDEL_ANCESTRAL`].
/// `locus_state` and `locus_carried` read that sentinel at a locus that is not a SNP.
fn is_snp_locus(locus: &Locus) -> bool {
    locus.derived.chars().count() == 1 && locus.ancestral.chars().count() <= 1
}

/// The sentinel that goes at the anchor position of an indel locus when the sample **carries**
/// the derived insertion or deletion. It is not a nucleotide, so it never collides with the base
/// call of a SNP. If an indel and a SNP do share a position, which is rare, both fall to a
/// no-call. Neither becomes a wrong call.
pub const INDEL_DERIVED: char = '+';
/// The sentinel that goes at an indel locus that the sample does **not** carry. Such a locus
/// covers the reference, and it is ancestral.
pub const INDEL_ANCESTRAL: char = '-';

/// Does the sample carry the derived allele of this locus? The code accepts the strand
/// complement of the derived base, except at a SNP that is ambiguous about the strand. Some tree
/// variants record their alleles on the other strand from the reference that the caller used on
/// the alignment. See [`locus_state`].
///
fn locus_carried(locus: &Locus, calls: &HashMap<i64, char>) -> bool {
    if !is_snp_locus(locus) {
        // Indel locus: carried iff the indel genotyper resolved it to the derived sentinel.
        return calls.get(&locus.position) == Some(&INDEL_DERIVED);
    }
    let Some(d) = locus.derived.chars().next().map(|c| c.to_ascii_uppercase()) else {
        return false;
    };
    let Some(b) = calls.get(&locus.position).map(|c| c.to_ascii_uppercase()) else {
        return false;
    };
    if b == d {
        return true;
    }
    let ambiguous = locus.ancestral.chars().next().is_some_and(|a| strand_ambiguous(a, d));
    !ambiguous && complement_base(b) == d
}

/// The state of the sample at one SNP that defines a branch. The sample carries the derived
/// allele, or it carries the ancestral allele, or it has no confident call. A locus with no
/// derived allele, which holds only an indel or no marker, gives `NoCall`.
///
/// Some haplotree variants record their ancestral and derived alleles on the **other strand**
/// from the reference that the caller genotyped the alignment against. FTDNA and YBrowse report
/// on the discovery strand of the SNP. A clean read then shows the complement of both tree
/// alleles, and it matches neither the literal ancestral allele nor the literal derived one. The
/// code also accepts a match on the strand complement, which is what the chip reconciliation
/// does.
///
/// A SNP that is ambiguous about the strand (`A↔T` or `C↔G`) is the exception. There the
/// complement of one allele *is* the other, and the data does not tell you the strand. Those
/// keep a strict literal match. A base that matches neither strand of either allele is a true
/// third allele, and it gives `NoCall`. It does not contradict the branch.
fn locus_state(locus: &Locus, calls: &HashMap<i64, char>) -> CallState {
    // The code can not read an indel or an MNP locus as a single base. It genotypes those
    // separately, and their resolved state comes in as a sentinel at the anchor, which is
    // [`INDEL_DERIVED`] or [`INDEL_ANCESTRAL`].
    if !is_snp_locus(locus) {
        return match calls.get(&locus.position) {
            Some(&INDEL_DERIVED) => CallState::Derived,
            Some(&INDEL_ANCESTRAL) => CallState::Ancestral,
            _ => CallState::NoCall,
        };
    }
    let Some(d) = locus.derived.chars().next().map(|c| c.to_ascii_uppercase()) else {
        return CallState::NoCall;
    };
    let a = locus.ancestral.chars().next().map(|c| c.to_ascii_uppercase());
    let Some(b) = calls.get(&locus.position).map(|c| c.to_ascii_uppercase()) else {
        return CallState::NoCall;
    };
    // Reference-strand match first.
    if b == d {
        return CallState::Derived;
    }
    if Some(b) == a {
        return CallState::Ancestral;
    }
    // Opposite-strand match (non-ambiguous SNPs only).
    let ambiguous = a.is_some_and(|a| strand_ambiguous(a, d));
    if !ambiguous {
        let bc = complement_base(b);
        if bc == d {
            return CallState::Derived;
        }
        if Some(bc) == a {
            return CallState::Ancestral;
        }
    }
    CallState::NoCall // a genuine third allele — not a confident call for this branch
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
        // An indel or MNP locus counts only when the indel genotyper resolved it, which is when
        // a sentinel is present. An indel with no call must not make the expected set or the
        // matched set of a node larger. A SNP locus counts as before, and a SNP with no call
        // still goes into `expected`.
        if !is_snp_locus(locus) {
            match calls.get(&locus.position) {
                Some(&INDEL_DERIVED) | Some(&INDEL_ANCESTRAL) => {}
                _ => continue,
            }
        }
        if !on_path.insert(locus.position) {
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

/// A map from every position that defines a branch to the name of a haplogroup that uses it. It
/// puts a note on a private variant that is off the path. A recurrent position keeps one name.
pub fn tree_positions(tree: &HaploTree) -> HashMap<i64, String> {
    let mut m = HashMap::new();
    for n in tree.nodes.values() {
        for l in &n.loci {
            m.entry(l.position).or_insert_with(|| n.name.clone());
        }
    }
    m
}

/// The polarity map for the consensus interpreter: **SNP name → (ancestral, derived)**, over
/// every locus in the tree that defines a branch. This is the polarity of each SNP in the tree of
/// record. `navigator_domain::consensus::interpret` applies it when it reads the data, so a
/// corrected tree changes the states and no sample needs a new genotype.
///
/// Use this for any parsed [`HaploTree`], such as the mtDNA rCRS tree or an FTDNA tree, where
/// there is no JSON polarity map. For the DecodingUs Y JSON, use [`decodingus_polarity_map`]
/// instead, which gives the true phylogenetic polarity.
///
/// The code skips a locus that has no name or no derived allele. A name that occurs more than
/// once keeps the polarity of its first occurrence.
pub fn polarity_from_tree(tree: &HaploTree) -> std::collections::BTreeMap<String, (String, String)> {
    let mut m: std::collections::BTreeMap<String, (String, String)> = std::collections::BTreeMap::new();
    for n in tree.nodes.values() {
        for l in &n.loci {
            if l.name.trim().is_empty() || l.derived.is_empty() {
                continue;
            }
            m.entry(l.name.trim().to_uppercase())
                .or_insert_with(|| (l.ancestral.clone(), l.derived.clone()));
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

/// The positions of the SNPs that define a branch, on the path from the root to `node_id`. They
/// are the backbone of the placement.
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

/// The state of the sample at a SNP that defines a branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallState {
    /// Carries the branch's derived allele.
    Derived,
    /// Carries the ancestral allele (this branch's split is not supported).
    Ancestral,
    /// No confident base call at this position.
    NoCall,
}

/// One SNP that defines a branch, with the state of the sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpEvidence {
    pub name: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub state: CallState,
    /// The sample's **observed base** at this position (`None` = no call). Carried alongside the
    /// imputed `state` so downstream consumers can persist the raw base and re-impute later.
    pub base: Option<char>,
}

/// A child branch below the terminal that the report names, with the evidence at each SNP. That
/// evidence explains why the descent went into that branch, or why it did not.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchEvidence {
    pub name: String,
    pub snps: Vec<SnpEvidence>,
    /// The count of the SNPs that define this branch, and that the sample carries as derived.
    pub derived: usize,
}

/// For the node `node_id`, which is usually the terminal that the report names, examine the SNPs
/// that define each child branch against the sample `calls`. Each SNP gets `Derived`,
/// `Ancestral` or `NoCall`.
///
/// This explains why the descent stopped. A child whose SNPs are all `Ancestral` is a split that
/// the data does not support. A child with `NoCall` SNPs has too little coverage to resolve. The
/// result leaves out a child that has no SNP to define it.
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
                base: calls.get(&l.position).copied(),
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

/// One node on the path from the root to the terminal. It holds the SNPs that define that node,
/// and the state of the sample at each of those SNPs. It is the grouped form of
/// [`lineage_evidence`], and it draws a descent report in the YFull style, one row for each
/// node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NodeEvidence {
    pub name: String,
    /// True for the reported terminal (the deepest node on the path).
    pub is_terminal: bool,
    pub snps: Vec<SnpEvidence>,
}

/// Group the path from the root to `terminal_id` into the evidence at each node, in root to
/// terminal order. The code walks the tree from the terminal up to the root. At each node it
/// attaches the loci of that node, with the state of the sample from `state_by_name`. An
/// equivalent that the sample did not call gets `NoCall`.
///
/// The key is the **SNP name**, which does not depend on the build. A name such as `M269` is the
/// same in every coordinate system. A cached variant profile that the code placed under any build
/// can colour a path in an FTDNA tree.
pub fn descent_by_node(
    tree: &HaploTree,
    terminal_id: i64,
    state_by_name: &HashMap<String, CallState>,
) -> Vec<NodeEvidence> {
    let parent = build_parent_map(tree);
    let mut ids = Vec::new();
    let mut cur = Some(terminal_id);
    while let Some(id) = cur {
        ids.push(id);
        cur = parent.get(&id).copied();
    }
    ids.reverse();

    ids.iter()
        .filter_map(|id| {
            let node = tree.nodes.get(id)?;
            let snps = node
                .loci
                .iter()
                .map(|l| SnpEvidence {
                    name: l.name.clone(),
                    position: l.position,
                    ancestral: l.ancestral.clone(),
                    derived: l.derived.clone(),
                    state: state_by_name.get(&l.name).copied().unwrap_or(CallState::NoCall),
                    base: None, // name-keyed (from a cached profile) — no raw call available here
                })
                .collect();
            Some(NodeEvidence {
                name: node.name.clone(),
                is_terminal: *id == terminal_id,
                snps,
            })
        })
        .collect()
}

/// The evidence at each SNP along the lineage from the root to `terminal_id`. It holds every SNP
/// that defines a node on the path, with the `Derived`, `Ancestral` or `NoCall` state of the
/// sample. Use it to compare which of those mutations a sample carries, for example a GRCh38
/// call against a CHM13 call that came through a liftover.
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
                base: calls.get(&l.position).copied(),
            });
        }
    }
    out
}

// ---- path-supported parsimony guard ------------------------------------------
//
// The Kulczynski `score` ranks every node by *proportional* set similarity. On real data that
// places the terminal well, and a check on GFX0457637 gave R-FGC29071. It has one weakness: the
// paralog artifact on the distal Y. A deep node can still score high on a few matches that are
// only coincidence. That happens when the path reaches the node only by a tunnel *through a
// branch that the sample contradicts*.
//
// Parsimony guards that exact failure. It refuses any candidate whose lineage from the root to
// the node crosses a branch that the sample contradicts. It does not change the proportional
// rank, which gets the clean case right.
//
// An earlier try used a router in the descent style, which follows the most derived subtree. It
// went onto a wrong fork with more branches on the 4x GFX sample. An absolute count of derived
// sites prefers a long path with many branches, where the Kulczynski proportion does not. The
// proportional rank plus this guard is the combination that the checks support.
//
// One paralog defence is still open: the *false-positive* case, where the wrong branch carries
// derived calls that are not real, and not honest ancestral ones. The haploid allele-balance
// filter covers that, and it is a separate Phase-1 item. See PangenomeExpansion.md.

/// The tally at each node. It covers the SNPs that define a branch and that the code can
/// examine, which are the loci with a derived allele. It counts three things: how many the sample
/// calls derived, how many it calls ancestral, and how many it has no confident base for. An
/// ancestral call is a contradiction.
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

/// The public view of the `(derived, ancestral, no-call)` tally of a node against `calls`, over
/// the SNPs that define the node. Use it for diagnostics, and to follow a placement path.
pub fn node_call_counts(tree: &HaploTree, calls: &HashMap<i64, char>, node_id: i64) -> (usize, usize, usize) {
    tree.nodes
        .get(&node_id)
        .map(|n| node_counts(n, calls))
        .unwrap_or((0, 0, 0))
}

/// Find a node by name, for the branch-report tool. It matches a **haplogroup name** such as
/// `R-FGC29071`, or the name of any **marker that defines the node**, such as `FGC29071`. The
/// case does not matter. It returns the node id.
///
/// A match on a marker name wins over no match at all. When more than one node matches, the
/// function returns the first one in a stable order of the ids. In practice a marker that defines
/// a node belongs to that node alone, so the result is deterministic.
pub fn find_node(tree: &HaploTree, query: &str) -> Option<i64> {
    let q = query.trim().to_ascii_uppercase();
    let mut ids: Vec<i64> = tree.nodes.keys().copied().collect();
    ids.sort_unstable();
    // Try a match on the haplogroup name first, then on the marker name. Both go over the
    // stable order of the ids.
    ids.iter()
        .find(|&&id| tree.nodes.get(&id).is_some_and(|n| n.name.to_ascii_uppercase() == q))
        .or_else(|| {
            ids.iter().find(|&&id| {
                tree.nodes
                    .get(&id)
                    .is_some_and(|n| n.loci.iter().any(|l| l.name.to_ascii_uppercase() == q))
            })
        })
        .copied()
}

/// One row of a branch report. It holds a marker that defines a node in the subtree that the
/// report covers, with the state of the sample. `node` and `parent` are haplogroup names. `snp`
/// carries the marker, the observed base, and the derived, ancestral or no-call state, which
/// [`locus_state`] gives.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchReportRow {
    pub node: String,
    pub parent: String,
    pub snp: SnpEvidence,
}

/// Every marker that defines a node in the subtree below `root_id`, in **pre-order**. The own
/// markers of a node come before the markers of its descendants. The code scores each marker
/// against `calls`.
///
/// The walk goes over the children in a stable order, by name and then by id, so the report is
/// deterministic. `max_depth` limits how far the walk goes down from the root. `None` sets no
/// limit, and `Some(0)` gives the markers of the root node alone.
///
/// This function covers the subtree of the descendants. [`lineage_evidence`] is its counterpart,
/// and it walks the ancestors from the root down to the node.
pub fn subtree_report(
    tree: &HaploTree,
    calls: &HashMap<i64, char>,
    root_id: i64,
    max_depth: Option<usize>,
) -> Vec<BranchReportRow> {
    let parent = build_parent_map(tree);
    let mut out = Vec::new();
    let mut stack = vec![(root_id, 0usize)];
    while let Some((id, depth)) = stack.pop() {
        let Some(node) = tree.nodes.get(&id) else { continue };
        let parent_name = parent
            .get(&id)
            .and_then(|p| tree.nodes.get(p))
            .map(|p| p.name.clone())
            .unwrap_or_default();
        for locus in &node.loci {
            out.push(BranchReportRow {
                node: node.name.clone(),
                parent: parent_name.clone(),
                snp: SnpEvidence {
                    name: locus.name.clone(),
                    position: locus.position,
                    ancestral: locus.ancestral.clone(),
                    derived: locus.derived.clone(),
                    state: locus_state(locus, calls),
                    base: calls.get(&locus.position).copied(),
                },
            });
        }
        let descend = match max_depth {
            Some(m) => depth < m,
            None => true,
        };
        if descend {
            // Push the children in the reverse of the stable order, so that they come off the
            // stack in pre-order, from the lowest up.
            let mut kids = node.children.clone();
            kids.sort_by(|a, b| {
                let na = tree.nodes.get(a).map(|n| n.name.as_str()).unwrap_or("");
                let nb = tree.nodes.get(b).map(|n| n.name.as_str()).unwrap_or("");
                na.cmp(nb).then(a.cmp(b))
            });
            for c in kids.into_iter().rev() {
                stack.push((c, depth + 1));
            }
        }
    }
    out
}

/// A node is *contradicted* when the sample carries the ancestral allele at more of its SNPs
/// than it carries the derived allele. The sample then clearly does **not** belong to this
/// branch.
///
/// A node with no evidence, where every SNP is a no-call and `d == a == 0`, is *not*
/// contradicted. The walk passes through it, so low coverage never blocks a lineage. One stray
/// ancestral at a node that the data otherwise supports, where `d >= a`, is acceptable for the
/// same reason.
fn is_contradicted(node: &HaploNode, calls: &HashMap<i64, char>) -> bool {
    let (d, a, _) = node_counts(node, calls);
    a > d
}

/// The count of derived branch SNPs that must occur *below* a contradicted ancestor, on the
/// path toward the candidate. Below that count, the contradiction keeps its veto on the
/// lineage.
///
/// One ancestral call at a thin intermediate node is usually an artifact of the genotype. It is
/// common on targeted-Y data such as FTDNA Big Y, whose coverage gaps turn most intermediate SNPs
/// into no-calls. Without this threshold, such a call would veto a whole deep lineage that the
/// terminal strongly supports. One stray ancestral at R-Z16250 blocked R-CTS4466, which had 10
/// derived and 0 ancestral.
///
/// A true off-branch *tunnel* artifact is different. It carries only one or two coincidental hits
/// below the contradicted branch point, so the veto stays. The threshold sits above the noise of
/// coincidence, and well below the derived count of a true clade.
const REDEEM_DERIVED: usize = 4;

/// Is the lineage from the root to `node_id` free of every contradicted branch that nothing
/// below it redeems?
///
/// An off-path paralog artifact lies below a branch for which the sample is ancestral, so it
/// fails this guard. The true lineage clears the guard, because along its length the sample is
/// derived or has only a no-call. [`score`] uses this to veto a tunnel artifact that would
/// otherwise rank high.
///
/// A contradicted ancestor vetoes only when derived support further down the path does not
/// *redeem* it. Below it, [`REDEEM_DERIVED`] derived SNPs or more confirm the branch. They
/// override one stray ancestral at a thin intermediate node, which is a Big Y miscall. A tunnel
/// of coincidence keeps its veto, because it is a contradicted branch point with only one or two
/// hits below it.
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
            // Nothing can redeem a **confident divergence**, where the sample carries the
            // ancestral allele here and *none* of the derived SNPs of the node. The sample left
            // the lineage above this node, so everything below is a different branch. It is a
            // clade beside this one, or a block of homoplasy or indels that happens to score.
            //
            // To redeem such a node lets a large block of false derived calls tunnel through
            // to a terminal that the sample can not reach. That block lies on a parallel
            // branch.
            //
            // Only a *mixed block* can take redemption from strong derived support below. A
            // mixed block has `d > 0`. The sample carries some of the derived SNPs of the node,
            // and the ancestral ones are a downstream split that the tree has not resolved.
            if derived[i] == 0 {
                return false;
            }
            // Derived SNPs strictly below this contradicted node, toward the candidate.
            let derived_below: usize = derived[i + 1..].iter().sum();
            if derived_below < REDEEM_DERIVED {
                return false;
            }
        }
    }
    true
}

/// The count of derived SNPs of a child that the sample must carry before [`deepen_terminal`]
/// goes down into that child. Two shared derived mutations that are independent confirm that the
/// sample belongs there, and they stay robust against one recurrent or false match.
const MIN_DERIVED_TO_DEEPEN: usize = 2;

/// From the `start` that the guard chose, go further down into any child that the sample has
/// *clearly entered*. Such a child meets two conditions: the sample carries at least
/// [`MIN_DERIVED_TO_DEEPEN`] of its derived SNPs, and nothing contradicts it, which means
/// `ancestral ≤ derived`. The route follows the derived count, and a tie goes to the lower id.
///
/// This corrects a low call at a **tree node that nobody has split yet**. Published trees hold
/// many of those, and the FTDNA tree most of all. Take a node whose SNP block has no
/// sub-branches yet. A sample on one sub-lineage is then derived for the SNPs of its own line.
/// It is ancestral for the SNPs of the *other* sub-lineages, which the tree has not split.
///
/// The node then looks half ancestral. Its proportional [`score`] falls just below the score of
/// its parent, and the guard stops one node too high. The sample truly carries some of the
/// mutations of the node. Those ancestral SNPs are a downstream split that nobody has resolved,
/// and they are not a contradiction.
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

// ---- induced subtree (the block-tree substrate) ------------------------------

/// One node of an [`induced_subtree`]. It holds the equivalent SNPs that define the branch. It
/// is a *block* in the sense of the FTDNA "block tree". Every sample below this node carries all
/// of `loci`, and no observation separates them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InducedNode {
    pub id: i64,
    pub name: String,
    /// Parent **within the induced subtree** (`None` at an induced root).
    pub parent: Option<i64>,
    /// The depth **inside the induced subtree**, where the root is 0. This is the coordinate for
    /// the layout. The induced root is the common ancestor of the members, and not the root of
    /// the tree. A depth from the full tree would waste the whole left margin on nodes that hold
    /// no sample.
    pub depth: usize,
    /// The equivalent SNPs that define this branch, which are the own loci of the node.
    pub loci: Vec<Locus>,
}

/// A map from a haplogroup name to a node id, over the whole tree. A caller that places *many*
/// samples must build this once. The path for one subject, which is the caller of
/// [`crate::haplo::descent_by_node`], does a linear scan for its one terminal. That is correct
/// for one subject, and quadratic for a cohort.
///
/// The code takes the names to be unique inside a tree. If two are the same, the **lowest id**
/// wins, so the result does not depend on the iteration order of a `HashMap`.
pub fn name_index(tree: &HaploTree) -> HashMap<&str, i64> {
    let mut idx: HashMap<&str, i64> = HashMap::with_capacity(tree.nodes.len());
    for n in tree.nodes.values() {
        idx.entry(n.name.as_str())
            .and_modify(|id| *id = (*id).min(n.id))
            .or_insert(n.id);
    }
    idx
}

/// The **induced subtree** that covers `terminals`. It holds every node that lies on a path from
/// the root to a terminal, for one terminal or more. It comes out in pre-order, so a parent
/// always comes before its children.
///
/// This is the skeleton of a cohort block tree. It is the union of the descent paths of the
/// members, which is exactly the set of branches that any of them share. The code ignores an id
/// that `tree` does not hold. A caller may give terminals that came from a different provider or
/// build, and it does not have to remove them first.
///
/// The order of the nodes beside each other goes by `(name, id)`, and the order of the roots does
/// the same. The output order is deterministic, and it does not depend on the iteration order of
/// a `HashMap`. The layout tests and the snapshot tests need that.
pub fn induced_subtree(tree: &HaploTree, terminals: &[i64]) -> Vec<InducedNode> {
    let parent = build_parent_map(tree);

    // Every node on a path from the root to a terminal. `seen` also breaks a cycle in a tree that
    // is not correct. A node that the set already holds means that the code also kept the rest of
    // its path, so the walk up can stop there.
    let mut kept: HashSet<i64> = HashSet::new();
    for &t in terminals {
        if !tree.nodes.contains_key(&t) {
            continue; // terminal not in this tree (provider/build skew) — the caller reports it
        }
        let mut cur = Some(t);
        while let Some(id) = cur {
            if !kept.insert(id) {
                break;
            }
            cur = parent.get(&id).copied();
        }
    }
    if kept.is_empty() {
        return Vec::new();
    }

    // The induced roots are the nodes that the code kept, and whose parent it did not keep.
    // There is usually exactly one, which is the common ancestor of the members. But a tree with
    // more than one root can give more than one, and the DecodingUs document does have a `roots`
    // array. A cohort that reaches across those roots can do the same.
    let sort_key = |id: &i64| tree.nodes.get(id).map(|n| (n.name.clone(), n.id));
    let mut roots: Vec<i64> = kept
        .iter()
        .copied()
        .filter(|id| !parent.get(id).is_some_and(|p| kept.contains(p)))
        .collect();
    roots.sort_by_key(sort_key);

    let mut out = Vec::with_capacity(kept.len());
    let mut stack: Vec<(i64, Option<i64>, usize)> = Vec::new();
    // In the reverse order, so that the nodes come off the stack in the sorted order.
    stack.extend(roots.iter().rev().map(|&id| (id, None, 0)));
    while let Some((id, par, depth)) = stack.pop() {
        let Some(node) = tree.nodes.get(&id) else { continue };
        out.push(InducedNode {
            id,
            name: node.name.clone(),
            parent: par,
            depth,
            loci: node.loci.clone(),
        });
        let mut children: Vec<i64> = node.children.iter().copied().filter(|c| kept.contains(c)).collect();
        children.sort_by_key(sort_key);
        stack.extend(children.into_iter().rev().map(|c| (c, Some(id), depth + 1)));
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

    /// A tree with more than one branch, for the induced-subtree and block cases:
    ///
    /// ```text
    /// root ──> R ──> R1 ──> R1a
    ///            └─> R2      └─> R1b
    /// ```
    ///
    /// `R1` carries two equivalent SNPs, which is the block case.
    const BRANCHY: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "R", "isRoot": false,
              "variants": [{"variant":"M207","position":100,"ancestral":"A","derived":"G"}], "children": [3, 6]},
        "3": {"haplogroupId": 3, "name": "R1", "isRoot": false,
              "variants": [{"variant":"M173","position":200,"ancestral":"C","derived":"T"},
                           {"variant":"M306","position":201,"ancestral":"G","derived":"A"}], "children": [4, 5]},
        "4": {"haplogroupId": 4, "name": "R1a", "isRoot": false,
              "variants": [{"variant":"M420","position":300,"ancestral":"A","derived":"T"}], "children": []},
        "5": {"haplogroupId": 5, "name": "R1b", "isRoot": false,
              "variants": [{"variant":"M343","position":400,"ancestral":"C","derived":"A"}], "children": []},
        "6": {"haplogroupId": 6, "name": "R2", "isRoot": false,
              "variants": [{"variant":"M479","position":500,"ancestral":"T","derived":"C"}], "children": []}
      }
    }"#;

    #[test]
    fn name_index_maps_every_haplogroup_name() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        let idx = name_index(&t);
        assert_eq!(idx.len(), 6);
        assert_eq!(idx.get("R1b"), Some(&5));
        assert_eq!(idx.get("root"), Some(&1));
        assert_eq!(idx.get("nope"), None);
    }

    #[test]
    fn induced_subtree_spans_only_the_members_paths() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        // Two members: one at R1a, one at R1b. R2 is off every path and must not appear.
        let nodes = induced_subtree(&t, &[4, 5]);
        let names: Vec<&str> = nodes.iter().map(|n| n.name.as_str()).collect();
        assert_eq!(names, vec!["root", "R", "R1", "R1a", "R1b"]);
        assert!(!names.contains(&"R2"));
    }

    #[test]
    fn induced_subtree_is_preorder_with_induced_parent_and_depth() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        let nodes = induced_subtree(&t, &[4, 5, 6]);
        // Pre-order: a parent is always emitted before its children.
        let pos = |name: &str| nodes.iter().position(|n| n.name == name).unwrap();
        for n in &nodes {
            if let Some(p) = n.parent {
                let parent_name = &nodes.iter().find(|x| x.id == p).unwrap().name;
                assert!(pos(parent_name) < pos(&n.name), "{} preceded its parent", n.name);
            }
        }
        // The depth counts from the induced root, and not from the full tree.
        assert_eq!(nodes[0].depth, 0);
        assert_eq!(nodes[0].parent, None);
        assert_eq!(nodes.iter().find(|n| n.name == "R1a").unwrap().depth, 3);
    }

    #[test]
    fn induced_subtree_carries_the_equivalent_snp_block() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        let nodes = induced_subtree(&t, &[4]);
        let r1 = nodes.iter().find(|n| n.name == "R1").unwrap();
        // The two SNPs of R1 are equivalent in the phylogeny. That pair *is* the block.
        let mut markers: Vec<&str> = r1.loci.iter().map(|l| l.name.as_str()).collect();
        markers.sort_unstable();
        assert_eq!(markers, vec!["M173", "M306"]);
    }

    #[test]
    fn induced_subtree_ignores_terminals_absent_from_the_tree() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        // 999 does not exist (provider/build skew); the real terminal still resolves.
        let nodes = induced_subtree(&t, &[4, 999]);
        let names: Vec<&str> = nodes.iter().map(|n| n.name.as_str()).collect();
        assert_eq!(names, vec!["root", "R", "R1", "R1a"]);
        // ...and a cohort of unknowns alone gives nothing. It does not panic.
        assert!(induced_subtree(&t, &[999]).is_empty());
        assert!(induced_subtree(&t, &[]).is_empty());
    }

    #[test]
    fn induced_subtree_is_deterministic_across_runs() {
        let t = parse_ftdna_json(BRANCHY).unwrap();
        // The iteration order of a HashMap changes from one process to the next. The output
        // order must not.
        let first = induced_subtree(&t, &[4, 5, 6]);
        for _ in 0..8 {
            assert_eq!(induced_subtree(&t, &[6, 5, 4]), first);
        }
    }

    #[test]
    fn parses_and_drops_positionless_variants() {
        let t = parse_ftdna_json(TREE).unwrap();
        assert_eq!(t.nodes.len(), 4);
        assert_eq!(t.nodes[&2].loci[0].position, 146);
        assert_eq!(t.nodes[&2].loci[0].derived, "G");
    }

    #[test]
    fn descent_by_node_buckets_path_with_state() {
        let t = parse_ftdna_json(TREE).unwrap();
        // The sample carries the derived allele at H (A146G) and at H2 (A263G). Nobody called
        // the SNP of H2a, which is C750T.
        let state: HashMap<String, CallState> = [
            ("A146G".to_string(), CallState::Derived),
            ("A263G".to_string(), CallState::Derived),
        ]
        .into_iter()
        .collect();
        let grouped = descent_by_node(&t, 4, &state); // terminal H2a

        // root → H → H2 → H2a. The root carries no locus that defines a branch.
        let names: Vec<&str> = grouped.iter().map(|n| n.name.as_str()).collect();
        assert_eq!(names, vec!["root", "H", "H2", "H2a"]);
        assert!(grouped[0].snps.is_empty()); // root has no loci
        assert!(!grouped[0].is_terminal && grouped[3].is_terminal);

        assert_eq!(grouped[1].snps[0].name, "A146G");
        assert_eq!(grouped[1].snps[0].state, CallState::Derived);
        assert_eq!(grouped[2].snps[0].state, CallState::Derived);
        // Nobody called the SNP that defines H2a, so it gives NoCall, which the UI draws grey.
        assert_eq!(grouped[3].snps[0].state, CallState::NoCall);
    }

    #[test]
    fn find_node_matches_haplogroup_or_marker_name_case_insensitively() {
        let t = parse_ftdna_json(TREE).unwrap();
        assert_eq!(find_node(&t, "H2"), Some(3)); // haplogroup name
        assert_eq!(find_node(&t, "h2"), Some(3)); // case-insensitive
        assert_eq!(find_node(&t, "A263G"), Some(3)); // defining-marker name → owning node
        assert_eq!(find_node(&t, "c750t"), Some(4)); // marker, case-insensitive
        assert_eq!(find_node(&t, "nope"), None);
    }

    #[test]
    fn subtree_report_walks_descendants_preorder_with_state() {
        let t = parse_ftdna_json(TREE).unwrap();
        // Derived at H(146) and H2(263); H2a(750) never called.
        let c = calls(&[(146, 'G'), (263, 'G')]);
        let rows = subtree_report(&t, &c, 2, None); // subtree rooted at H

        let seen: Vec<(&str, &str, &str, CallState)> = rows
            .iter()
            .map(|r| (r.node.as_str(), r.parent.as_str(), r.snp.name.as_str(), r.snp.state))
            .collect();
        assert_eq!(
            seen,
            vec![
                ("H", "root", "A146G", CallState::Derived),
                ("H2", "H", "A263G", CallState::Derived),
                ("H2a", "H2", "C750T", CallState::NoCall),
            ]
        );
        assert_eq!(rows[0].snp.base, Some('G'));
        assert_eq!(rows[2].snp.base, None);
    }

    #[test]
    fn subtree_report_depth_zero_is_root_node_only() {
        let t = parse_ftdna_json(TREE).unwrap();
        let rows = subtree_report(&t, &calls(&[]), 2, Some(0));
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].node, "H");
        assert_eq!(rows[0].snp.name, "A146G");
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
    fn decodingus_polarity_map_prefers_hs1_deterministically() {
        // A SNP whose hs1 and GRCh38 coords record *swapped* polarity: the map must always pick the
        // hs1 (native/authoritative) orientation, never a HashMap-order-dependent one. A second SNP
        // has only a GRCh38 coord → falls back to it.
        let json = r#"{
          "roots": [
            {"id": 1, "name": "R", "variants": [
                {"canonical_name": "A2627", "coordinates": {
                    "hs1": {"contig":"chrY","position":100,"ancestral":"C","derived":"T"},
                    "GRCh38": {"contig":"chrY","position":200,"ancestral":"T","derived":"C"}}},
                {"canonical_name": "GRCh38only", "coordinates": {
                    "GRCh38": {"contig":"chrY","position":300,"ancestral":"G","derived":"A"}}}],
             "children": []}
          ]
        }"#;
        let pol = decodingus_polarity_map(json).unwrap();
        assert_eq!(pol["A2627"], ("C".to_string(), "T".to_string())); // hs1 wins, every run
        assert_eq!(pol["GRCh38only"], ("G".to_string(), "A".to_string())); // fallback build
    }

    #[test]
    fn decodingus_link_polarity_overrides_swapped_coordinate() {
        // This is a real property of the DecodingUs data. About 1.4% of the variants carry a
        // `coordinates.ancestral/derived` in a build that is the *exchange* of the authoritative
        // `link_ancestral/link_derived` at the variant level. Both the tree parse and the
        // polarity map must trust `link_*`. If they do not, a backbone SNP that the sample
        // carries as derived reads as ancestral. That was the huF98AFD bug, where the report was
        // full of ancestral calls.
        let json = r#"{
          "roots": [
            {"id": 1, "name": "A0-T", "variants": [
                {"canonical_name": "A2614",
                 "link_ancestral": "G", "link_derived": "A",
                 "coordinates": {
                    "hs1": {"contig":"chrY","position":6964116,"ancestral":"A","derived":"G"},
                    "GRCh38": {"contig":"chrY","position":7311793,"ancestral":"A","derived":"G"}}}],
             "children": []}
          ]
        }"#;
        // parse: position from the build coord, polarity from link_* (G>A, not the coord's A>G).
        let t = parse_decodingus_json(json, "hs1").unwrap();
        let locus = &t.nodes[&1].loci[0];
        assert_eq!(locus.position, 6964116);
        assert_eq!((locus.ancestral.as_str(), locus.derived.as_str()), ("G", "A"));
        // A sample that carries the derived allele A now genotypes as Derived. Before the fix it
        // gave Ancestral.
        let calls: HashMap<i64, char> = [(6964116, 'A')].into_iter().collect();
        assert_eq!(locus_state(locus, &calls), CallState::Derived);
        // The polarity map agrees.
        let pol = decodingus_polarity_map(json).unwrap();
        assert_eq!(pol["A2614"], ("G".to_string(), "A".to_string()));
    }

    #[test]
    fn parse_decodingus_picks_target_build_and_flattens() {
        // In hs1, both M207 and M173 resolve. The code drops the variant that exists only in
        // GRCh38.
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
    fn a_node_whose_every_variant_lacks_the_build_survives_with_no_snps() {
        // This is the shape behind the short descent of `1087`. Most of the SNPs that DecodingUs
        // found, which carry a `DU` name, exist in CHM13 coordinates alone. Only a few hundred of
        // them went back to the older references. A terminal that one of them defines has
        // *nothing* under GRCh38.
        //
        // To drop the loci does not drop the node. The node stays, with its name and on the path,
        // and its `loci` is empty. `descent_by_node` then reports it correctly, with no SNPs. But
        // a renderer that hides an empty block, which is correct because the root is truly empty,
        // shows the lineage one branch short. The *name* of the terminal stays right.
        //
        // `DECODINGUS_NATIVE_BUILD` answers this: parse in hs1 wherever the join goes by SNP
        // name.
        let json = r#"{
          "roots": [
            {"id": 1, "name": "R-BY57568", "variants": [
                {"canonical_name": "BY57568", "link_ancestral": "C", "link_derived": "A",
                 "coordinates": {
                    "hs1":    {"contig":"chrY","position":3150089,"ancestral":"C","derived":"A"},
                    "GRCh38": {"contig":"chrY","position":3472892,"ancestral":"C","derived":"A"}}}],
             "children": [
               {"id": 2, "name": "R-DU17762", "variants": [
                   {"canonical_name": "DU17762", "link_ancestral": "G", "link_derived": "A",
                    "coordinates": {
                       "hs1": {"contig":"chrY","position":27785335,"ancestral":"G","derived":"A"}}}],
                "children": []}
             ]}
          ]
        }"#;
        let states: HashMap<String, CallState> = [
            ("BY57568".to_string(), CallState::Derived),
            ("DU17762".to_string(), CallState::Derived),
        ]
        .into_iter()
        .collect();

        // In GRCh38 the terminal is on the path, it has its name, and it is empty. That is the
        // bug.
        let g38 = parse_decodingus_json(json, "GRCh38").unwrap();
        assert!(g38.nodes.contains_key(&2), "the node itself is never dropped");
        assert!(g38.nodes[&2].loci.is_empty(), "its only locus has no GRCh38 coordinate");
        let d = descent_by_node(&g38, 2, &states);
        let terminal = d.last().expect("terminal is reported");
        assert_eq!(terminal.name, "R-DU17762");
        assert!(
            terminal.snps.is_empty(),
            "so a hide-empty-blocks renderer stops at R-BY57568"
        );

        // hs1: the same descent carries the marker, so the terminal block renders.
        let hs1 = parse_decodingus_json(json, DECODINGUS_NATIVE_BUILD_FOR_TEST).unwrap();
        let d = descent_by_node(&hs1, 2, &states);
        let terminal = d.last().unwrap();
        assert_eq!(terminal.name, "R-DU17762");
        let named: Vec<&str> = terminal.snps.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(named, vec!["DU17762"]);
        assert_eq!(terminal.snps[0].position, 27785335);
        assert_eq!(terminal.snps[0].state, CallState::Derived);
    }

    /// Mirrors `navigator_app::DECODINGUS_NATIVE_BUILD`, which this crate sits below and so can not
    /// import. The test above is the reason that constant exists.
    const DECODINGUS_NATIVE_BUILD_FOR_TEST: &str = "hs1";

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
        // H2 has a child H2a, whose derived allele is T at 750. The sample is ancestral, C, at
        // 750. The data does not support the split into H2a, and the report shows that at each
        // SNP.
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
        // This is the CHM13 trap. At a Y-SNP the tree gives ancestral=A and derived=G. The
        // analysis reference carries the DERIVED base G, because CHM13 chrY is HG002, which is a
        // haplogroup-J Y. The polarity must come from a comparison of the base of the SAMPLE
        // against the tree. It must never come from the reference.
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

        // This is the end-to-end case. A sample that is ANCESTRAL at the J-derived backbone site
        // 146 does not carry the mutation that defines H. The code must not place it into H.
        // That holds even though the CHM13 reference base there is the derived G. An assumption
        // that the REF base is the ancestral one would turn 146 around and descend wrongly. The
        // tree controls the polarity, so the data contradicts H and the call stays at the
        // root.
        let t = parse_ftdna_json(TREE).unwrap(); // root→H(146 A→G)→H2(263)→H2a(750)
        let c = calls(&[(146, 'A'), (263, 'G'), (750, 'T')]);
        assert!(
            !path_admissible(&t, &c, id_of(&t, "H")),
            "H is contradicted (sample ancestral at 146)"
        );
        assert_eq!(guarded_terminal(&t, &c), "root");
    }

    #[test]
    fn opposite_strand_reads_match_via_the_complement() {
        // A SNP that is not ambiguous. Its tree alleles sit on the other strand from the
        // reference that the caller genotyped the alignment against: ancestral=A, derived=C. A
        // derived sample, read on the reference strand, shows G, which is the complement of C. An
        // ancestral one shows T, the complement of A. Neither one matches a tree allele
        // literally, but the strand complement does.
        let locus = Locus {
            position: 146,
            ancestral: "A".into(),
            derived: "C".into(),
            name: "CTS-test".into(),
        };
        assert_eq!(locus_state(&locus, &calls(&[(146, 'G')])), CallState::Derived);
        assert_eq!(locus_state(&locus, &calls(&[(146, 'T')])), CallState::Ancestral);
        assert!(locus_carried(&locus, &calls(&[(146, 'G')])));
        assert!(!locus_carried(&locus, &calls(&[(146, 'T')])));

        // A SNP that is ambiguous about the strand (C↔G). The complement of the derived G is the
        // ancestral C, so the data does not tell you the strand. Keep a strict literal match, and
        // do not turn the base to its complement.
        let palindrome = Locus {
            position: 200,
            ancestral: "C".into(),
            derived: "G".into(),
            name: "PF-pal".into(),
        };
        assert_eq!(locus_state(&palindrome, &calls(&[(200, 'C')])), CallState::Ancestral);
        assert_eq!(locus_state(&palindrome, &calls(&[(200, 'G')])), CallState::Derived);
        // A genuine third allele (A/T at a C/G site) is unresolvable → NoCall, not Ancestral.
        assert_eq!(locus_state(&palindrome, &calls(&[(200, 'A')])), CallState::NoCall);
    }

    #[test]
    fn normalize_polarity_flips_ftdna_reference_as_ancestral_inversion() {
        // FTDNA records the GRCh38 reference base as the "ancestral" one. At PF1016 the
        // reference carries the derived allele. FTDNA lists T>C, where the true polarity, from
        // DecodingUs, is C>T.
        let mut tree = HaploTree { nodes: HashMap::new() };
        tree.nodes.insert(
            1,
            HaploNode {
                id: 1,
                name: "CT".into(),
                is_root: false,
                loci: vec![
                    Locus {
                        position: 100,
                        ancestral: "T".into(),
                        derived: "C".into(),
                        name: "PF1016".into(),
                    },
                    // This SNP already agrees with the reference map. It must not change.
                    Locus {
                        position: 200,
                        ancestral: "A".into(),
                        derived: "G".into(),
                        name: "M168".into(),
                    },
                    // The alleles sit on different strands: G>A against C>T. That is not a
                    // clean exchange of the two roles, so it must not change.
                    Locus {
                        position: 300,
                        ancestral: "G".into(),
                        derived: "A".into(),
                        name: "S3".into(),
                    },
                ],
                children: vec![],
            },
        );
        let reference: HashMap<String, (String, String)> = [
            ("PF1016".to_string(), ("C".to_string(), "T".to_string())),
            ("M168".to_string(), ("A".to_string(), "G".to_string())),
            ("S3".to_string(), ("C".to_string(), "T".to_string())),
        ]
        .into_iter()
        .collect();

        let flipped = normalize_polarity(&mut tree, &reference);
        assert_eq!(flipped, 1, "only the swapped PF1016 should flip");
        let loci = &tree.nodes[&1].loci;
        let pf = loci.iter().find(|l| l.name == "PF1016").unwrap();
        assert_eq!((pf.ancestral.as_str(), pf.derived.as_str()), ("C", "T"));
        let m168 = loci.iter().find(|l| l.name == "M168").unwrap();
        assert_eq!((m168.ancestral.as_str(), m168.derived.as_str()), ("A", "G"));
        let s3 = loci.iter().find(|l| l.name == "S3").unwrap();
        assert_eq!((s3.ancestral.as_str(), s3.derived.as_str()), ("G", "A"));
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
        // The sample is ancestral, C, at 750. That contradicts H2a, and the report falls back to
        // H2.
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[(146, 'G'), (263, 'G'), (750, 'C')]);
        assert!(!path_admissible(&t, &c, id_of(&t, "H2a")));
        assert!(path_admissible(&t, &c, id_of(&t, "H2")));
        assert_eq!(guarded_terminal(&t, &c), "H2");
    }

    #[test]
    fn no_calls_admit_the_whole_tree() {
        // With no calls at all, nothing contradicts a node, so every lineage passes the guard.
        // The guard is a veto and not a selector. Kulczynski still takes the root, because it
        // finds no match.
        let t = parse_ftdna_json(TREE).unwrap();
        let c = calls(&[]);
        assert!(path_admissible(&t, &c, id_of(&t, "H2a")));
        assert_eq!(guarded_terminal(&t, &c), "root");
    }

    // root -> H(146) -> B(500, contradicted) -> Bdeep(900, derived by coincidence).
    // Kulczynski goes to Bdeep, because that node matches 146 and 900. The guard must veto it,
    // because the path tunnels through the contradicted B, and fall back to H.
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
        // Kulczynski alone goes deeper, because of the match by coincidence...
        assert_eq!(score(&t, &c)[0].name, "Bdeep");
        // ...but Bdeep tunnels through the contradicted B, so the guard reports H.
        assert!(!path_admissible(&t, &c, id_of(&t, "Bdeep")));
        assert_eq!(guarded_terminal(&t, &c), "H");
    }

    // root -> H(146) -> Ins(an insertion G->GAGC at 200). A single-base genotype can't evaluate the
    // insertion, so the sample must never place onto the indel-defined node.
    const INDEL_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "Ins", "isRoot": false,
              "variants": [{"variant":"ins200","position":200,"ancestral":"G","derived":"GAGC"}], "children": []}
      }
    }"#;

    #[test]
    fn indel_locus_is_not_evaluable_from_a_raw_base() {
        let t = parse_ftdna_json(INDEL_TREE).unwrap();
        // The sample carries the derived allele at H (146). At the insertion position it carries
        // the anchor base G. A simple comparison of the first base reads that G as the derived
        // allele of the insertion. A *raw base* at an indel position, which is not the indel
        // sentinel, must NOT place the sample onto the node that the indel defines.
        let c = calls(&[(146, 'G'), (200, 'G')]);
        assert!(!locus_carried(&t.nodes[&3].loci[0], &c));
        assert_eq!(locus_state(&t.nodes[&3].loci[0], &c), CallState::NoCall);
        assert_eq!(node_call_counts(&t, &c, 3), (0, 0, 1)); // the sole locus is an un-called indel (no-call)
        assert_eq!(score(&t, &c)[0].name, "H");
        assert_eq!(guarded_terminal(&t, &c), "H");
    }

    #[test]
    fn indel_sentinel_drives_placement_and_state() {
        let t = parse_ftdna_json(INDEL_TREE).unwrap();
        let indel = &t.nodes[&3].loci[0];
        // Genotyper resolved the insertion as PRESENT: the derived sentinel places onto Ins.
        let derived = calls(&[(146, 'G'), (200, INDEL_DERIVED)]);
        assert!(locus_carried(indel, &derived));
        assert_eq!(locus_state(indel, &derived), CallState::Derived);
        assert_eq!(node_call_counts(&t, &derived, 3), (1, 0, 0));
        assert_eq!(guarded_terminal(&t, &derived), "Ins");

        // The genotyper resolved the indel as ABSENT. The ancestral sentinel keeps the terminal
        // at H, and the data contradicts Ins.
        let ancestral = calls(&[(146, 'G'), (200, INDEL_ANCESTRAL)]);
        assert!(!locus_carried(indel, &ancestral));
        assert_eq!(locus_state(indel, &ancestral), CallState::Ancestral);
        assert_eq!(node_call_counts(&t, &ancestral, 3), (0, 1, 0));
        assert_eq!(guarded_terminal(&t, &ancestral), "H");

        // No sentinel (un-called indel): does not inflate Ins's expected set → terminal stays H.
        let uncalled = calls(&[(146, 'G')]);
        assert_eq!(guarded_terminal(&t, &uncalled), "H");
    }

    // root -> H(146) -> B(500) -> Bdeep(five derived SNPs). The sample is ANCESTRAL at B, so it
    // carries NONE of the derived alleles of B. But by coincidence it matches all five SNPs of
    // Bdeep. That is a large block of homoplasy or indels on a clade beside this one, which
    // diverged earlier. The redeem clause, which needs REDEEM_DERIVED derived alleles below or
    // more, would tunnel to Bdeep. Nothing must ever redeem a *confident* divergence, where
    // d == 0.
    const CONFIDENT_DIVERGENCE_TREE: &str = r#"{
      "allNodes": {
        "1": {"haplogroupId": 1, "name": "root", "isRoot": true, "variants": [], "children": [2]},
        "2": {"haplogroupId": 2, "name": "H", "isRoot": false,
              "variants": [{"variant":"A146G","position":146,"ancestral":"A","derived":"G"}], "children": [3]},
        "3": {"haplogroupId": 3, "name": "B", "isRoot": false,
              "variants": [{"variant":"C500T","position":500,"ancestral":"C","derived":"T"}], "children": [4]},
        "4": {"haplogroupId": 4, "name": "Bdeep", "isRoot": false, "variants": [
                {"variant":"G900A","position":900,"ancestral":"G","derived":"A"},
                {"variant":"G901A","position":901,"ancestral":"G","derived":"A"},
                {"variant":"G902A","position":902,"ancestral":"G","derived":"A"},
                {"variant":"G903A","position":903,"ancestral":"G","derived":"A"},
                {"variant":"G904A","position":904,"ancestral":"G","derived":"A"}
              ], "children": []}
      }
    }"#;

    #[test]
    fn guard_never_redeems_a_confident_divergence() {
        let t = parse_ftdna_json(CONFIDENT_DIVERGENCE_TREE).unwrap();
        // Derived at H(146), and ancestral at B(500), so it carries none of the derived alleles
        // of B and d == 0. But it matches all five SNPs of Bdeep. That homoplasy block gives 5
        // derived alleles below B, which is past REDEEM_DERIVED.
        let c = calls(&[
            (146, 'G'),
            (500, 'C'),
            (900, 'A'),
            (901, 'A'),
            (902, 'A'),
            (903, 'A'),
            (904, 'A'),
        ]);
        // Kulczynski goes to Bdeep, because of the five matches by coincidence...
        assert_eq!(score(&t, &c)[0].name, "Bdeep");
        // ...but B is a confident divergence, with zero derived alleles, and nothing redeems it.
        // The derived block below it does not help. The terminal is H and not Bdeep. That is the
        // fix that stops an indel tunnel on a parallel branch.
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
        // 263 gives a no-call, because the coverage is low. D is *still* admissible, because a
        // no-call is not a contradiction. The guard never blocks a path for lack of coverage.
        // Kulczynski stops at H.
        let sparse = calls(&[(146, 'G')]);
        assert!(path_admissible(&t, &sparse, id_of(&t, "D")));
        assert_eq!(guarded_terminal(&t, &sparse), "H");
    }

    // root -> H(146) -> D, with three SNPs that define D. This covers the rule about the net
    // contradiction.
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

    // root → P (2 derived SNPs) → C, where C is a node that nobody has SPLIT yet. C holds 3 SNPs
    // that define it, and 3 SNPs of a sub-branch that the tree has not split. A sample on the
    // trunk of C carries the derived allele at the first 3, and the ancestral allele at the
    // other 3. C then looks half ancestral, and Kulczynski can rank it below P.
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
        // Deepen goes into C from P. C carries 3 derived alleles, which is 2 or more, and
        // nothing contradicts it, because 3 ancestral is not more than 3 derived.
        //
        // The condition where Kulczynski stops at the parent needs a long backbone. A check on
        // the real WGS229 short-read sample showed it: the guard stops at R-FGC29067, and deepen
        // recovers R-FGC29071.
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
        // 2 derived but 4 ancestral → contradicted (a > d), do not enter even at ≥2 derived.
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

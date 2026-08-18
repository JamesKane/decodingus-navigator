//! The project **block tree**. It is the cohort form of the descent report of one subject.
//!
//! For a project, this module builds the part of the haplotree that covers the terminal haplogroups
//! of its members. That subtree holds each branch that a member lies on.
//!
//! Each branch carries a *block*, which is the run of SNPs that define it and that the tree treats
//! as equivalent. Each member appears below its own terminal branch.
//!
//! This surface is the FTDNA "Block Tree", and it uses the placements that Navigator already made.
//!
//! Two rules control this module:
//!
//! - **It reads a placement, and it never makes one.** Each terminal comes from
//!   `haplogroup_terminals`. The subjects table and the project report use the same
//!   reconciliation. No code here can move a subject.
//! - **The report names a member with no placement** ([`UnplacedMember`]). It does not remove that
//!   member. In a cohort from many laboratories, a difference between providers and builds is
//!   normal. Without those members, the tree looks complete when it is not.
//!
//! The design is in `documents/design/project-block-tree.md`.

use std::collections::BTreeSet;

use super::*;

/// The minimum length of a run of branches that the code joins into one. Each branch in the run has
/// one child and no member.
///
/// One branch between two splits has a name that the reader needs. A run of two or more such
/// branches only fills space between the splits that the cohort resolves.
pub const COLLAPSE_MIN_RUN: usize = 2;

impl App {
    /// Build the [`ProjectBlockTree`] of `project_id`.
    ///
    /// The method returns `Ok(None)` only when the project has no member.
    ///
    /// A project with no placed member still gives a tree. That tree has an empty `blocks` list,
    /// and it holds each member in `unplaced`. This answer is a useful one, because it tells the
    /// user that the app placed nothing yet.
    ///
    /// That answer also costs nothing. When no member has a terminal, the method does not download
    /// the tree and does not parse it. That document is many MB.
    pub async fn project_block_tree(
        &self,
        project_id: i64,
        dna: DnaType,
    ) -> Result<Option<ProjectBlockTree>, AppError> {
        let members = biosample::list_members_for_project(self.store.pool(), project_id).await?;
        if members.is_empty() {
            return Ok(None);
        }

        // One reconciliation covers the full workspace. The code does not send one query for each
        // member. `project_report` and `project_str_overview` call the same function.
        let terminals = self.haplogroup_terminals().await?;

        // The guid, the display name, and the terminal of each member, for the lineage that the
        // caller requested.
        let wanted: Vec<(SampleGuid, String, Option<String>)> = members
            .iter()
            .map(|b| {
                let t = terminals.get(&b.guid).and_then(|(y, mt)| match dna {
                    DnaType::Y => y.clone(),
                    DnaType::Mt => mt.clone(),
                });
                (b.guid, b.donor_identifier.clone(), t)
            })
            .collect();

        // Do the fast test first, as `descent_report` does. With no placement there is no tree to
        // draw. A download and a parse of a document of many MB, only to learn that, is work with
        // no result. This state is common in a project that a user imported a moment ago.
        if wanted.iter().all(|(_, _, t)| t.is_none()) {
            let mut unplaced: Vec<UnplacedMember> = wanted
                .into_iter()
                .map(|(guid, name, terminal)| UnplacedMember { guid, name, terminal })
                .collect();
            unplaced.sort_by(|a, b| (&a.name, a.guid.0).cmp(&(&b.name, b.guid.0)));
            return Ok(Some(ProjectBlockTree {
                dna,
                blocks: Vec::new(),
                unplaced,
                // This value is the provider in the settings. The code downloaded no tree, so it
                // also used no second provider. `build_key` stays empty for the same reason,
                // because the code parsed nothing.
                provider: match y_tree_provider() {
                    YTreeProvider::DecodingUs => "decodingus".to_string(),
                    YTreeProvider::Ftdna => "ftdna".to_string(),
                },
                build_key: String::new(),
                candidate_conflicts: 0,
                candidate_recurrent: 0,
            }));
        }

        // The provider and the build key come from the download, and not from the settings. The
        // mtDNA path uses the FTDNA tree when the code can not remap the DecodingUs tree. The loci
        // of that path are rCRS loci in each case, and they are not in the Y coordinate space.
        let (tree, provider, build_key) = match dna {
            DnaType::Y => match y_tree_provider() {
                YTreeProvider::DecodingUs => {
                    let build_key = self.project_build_key(&members).await;
                    let json = self.fetch_decodingus_y_tree().await?;
                    let tree =
                        navigator_analysis::haplo::parse_decodingus_json(&json, build_key).map_err(AppError::Import)?;
                    (tree, "decodingus", build_key)
                }
                YTreeProvider::Ftdna => {
                    let json = self.fetch_ftdna_y_tree().await?;
                    let tree = navigator_analysis::haplo::parse_ftdna_json(&json).map_err(AppError::Import)?;
                    // FTDNA publishes its Y tree on GRCh38. The build of each member does not
                    // change that.
                    (tree, "ftdna", "GRCh38")
                }
            },
            DnaType::Mt => {
                let (tree, provider) = self.mt_tree_rcrs().await?;
                (tree, provider, "rCRS")
            }
        };

        // One name index covers the full cohort. The path for one subject reads the node map from
        // start to end to find its terminal. Here that method would cost O(n²).
        let index = navigator_analysis::haplo::name_index(&tree);
        // Find the placements first. The private-Y read below then covers only the members that
        // appear on the tree. In a real cohort that group is small: 243 of 1,881 members on
        // R1b-CTS4466Plus. The view can draw no private variant of a member with no placement.
        let mut placed: Vec<(SampleGuid, String, i64)> = Vec::new();
        let mut unplaced = Vec::new();
        for (guid, name, terminal) in wanted {
            match terminal.as_deref().and_then(|t| index.get(t).copied()) {
                Some(id) => placed.push((guid, name, id)),
                None => unplaced.push(UnplacedMember { guid, name, terminal }),
            }
        }

        // One bulk load for the placed members' private-Y. Absent = never computed (`None`), which
        // the view must not confuse with "computed, none found" (`Some(0)`).
        let placed_guids: Vec<SampleGuid> = placed.iter().map(|(g, _, _)| *g).collect();
        let private = self.private_y_for_biosamples(&placed_guids).await?;
        let mut at_node: HashMap<i64, Vec<BlockMember>> = HashMap::new();
        for (guid, name, id) in placed {
            at_node.entry(id).or_default().push(BlockMember {
                guid,
                name,
                private_novel: private.get(&guid).map(|b| b.novel_in_unique_sequence()),
                private_publishable: private
                    .get(&guid)
                    .map(|b| b.publishable_count(crate::PublishGate::default())),
                private_total: private.get(&guid).map(|b| b.variants.len()),
            });
        }

        let terminal_ids: Vec<i64> = at_node.keys().copied().collect();
        let induced = navigator_analysis::haplo::induced_subtree(&tree, &terminal_ids);
        let mut blocks: Vec<Block> = induced
            .into_iter()
            .map(|n| Block {
                members: at_node.remove(&n.id).unwrap_or_default(),
                node_id: n.id,
                name: n.name,
                parent: n.parent,
                depth: n.depth,
                loci: n.loci,
                subtree_members: 0, // filled by `roll_up_subtree_members`
                collapsed: Vec::new(),
                candidate: false,
                evidence: Vec::new(),
            })
            .collect();
        // Stable leaf order, so the layout does not reshuffle between opens.
        for b in &mut blocks {
            b.members.sort_by(|x, y| (&x.name, x.guid.0).cmp(&(&y.name, y.guid.0)));
        }
        roll_up_subtree_members(&mut blocks);
        let blocks = collapse_blocks(blocks, COLLAPSE_MIN_RUN);
        // The code adds each candidate *after* the collapse step. A candidate is a leaf with
        // members, so the collapse can never absorb it. An earlier insert only makes the collapse
        // examine nodes that the code made.
        let (blocks, candidate_conflicts, candidate_recurrent) = insert_candidate_branches(blocks, &private);
        unplaced.sort_by(|a, b| (&a.name, a.guid.0).cmp(&(&b.name, b.guid.0)));

        Ok(Some(ProjectBlockTree {
            dna,
            blocks,
            unplaced,
            provider: provider.to_string(),
            build_key: build_key.to_string(),
            candidate_conflicts,
            candidate_recurrent,
        }))
    }

    /// The one coordinate space for the tree of the cohort. It is the most frequent DecodingUs
    /// build key across the alignments of the members. The default is `hs1`.
    ///
    /// A cohort holds more than one build, so there is no answer for one subject, as there is in
    /// `descent_report`.
    ///
    /// One choice is safe, because the node names and the shape of the tree do not depend on the
    /// build. Only the *positions* of the loci depend on it. The aggregate carries the key, so the
    /// view can name the space that it shows.
    ///
    /// Two keys with the same count give the key that sorts first. So the choice does not depend on
    /// the order of a map.
    async fn project_build_key(&self, members: &[Biosample]) -> &'static str {
        let guids: Vec<SampleGuid> = members.iter().map(|b| b.guid).collect();
        let Ok(alns) = alignment::list_for_biosamples(self.store.pool(), &guids).await else {
            return "hs1";
        };
        let mut votes: HashMap<&'static str, usize> = HashMap::new();
        for (_, a) in &alns {
            if let Some(k) = decodingus_build_key(&a.reference_build) {
                *votes.entry(k).or_default() += 1;
            }
        }
        votes.into_iter().max_by_key(|&(k, n)| (n, k)).map_or("hs1", |(k, _)| k)
    }
}

/// Fill `subtree_members`, which holds the members at each block and below it.
///
/// The `blocks` list is in pre-order. So a read from the end to the start reaches each child before
/// its parent, and one pass is enough.
fn roll_up_subtree_members(blocks: &mut [Block]) {
    let mut from_children: HashMap<i64, usize> = HashMap::with_capacity(blocks.len());
    for i in (0..blocks.len()).rev() {
        let total = blocks[i].members.len() + from_children.get(&blocks[i].node_id).copied().unwrap_or(0);
        blocks[i].subtree_members = total;
        if let Some(p) = blocks[i].parent {
            *from_children.entry(p).or_default() += total;
        }
    }
}

/// Make a [`Locus`] value for a private variant, which has no name. The view then draws the shared
/// variants of a candidate branch through the same code as the SNPs of a named branch.
fn private_locus(v: &PrivateVariant) -> Locus {
    Locus {
        position: v.position,
        ancestral: v.reference.to_string(),
        derived: v.alternate.to_string(),
        name: format!("chrY:{}", v.position),
    }
}

/// The positions that a subject carries as **new-branch candidates with high confidence**. Such a
/// position is new, so the tree does not hold it, and it is in unique sequence.
///
/// The function removes a *known* variant that is off the path, by design. That variant supports a
/// finer branch that already exists. It asks a question about the placement, and not about a new
/// branch.
///
/// The function also removes a call in a structural region. A palindrome and an amplicon on chrY
/// hold paralogs. Two men with the "same" call in such a region share a mapping artefact more often
/// than they share an ancestor.
///
/// Shared noise makes a branch that does not exist. That result is the one fault that this feature
/// must not produce.
fn candidate_positions(bucket: &PrivateBucket) -> BTreeSet<i64> {
    let novel: BTreeSet<i64> = bucket
        .variants
        .iter()
        .filter(|v| v.class == PrivateClass::Novel && v.region.is_none())
        .map(|v| v.position)
        .collect();
    drop_clustered(&novel)
}

/// The minimum distance between two new calls that come from separate mutations.
///
/// A real Y mutation is far from the next one, at a distance of megabases. A group of "new" calls
/// inside tens of bases comes from one read that the mapper placed wrongly, and that read gives
/// some false SNVs. The GVCF path applies a depth limit for the same reason.
///
/// In the CTS4466 cohort, almost every first candidate branch came from such a group. One group
/// held six positions inside 32 bp, with gaps of 5 bp to 8 bp.
const CANDIDATE_MIN_SEPARATION_BP: i64 = 100;

/// Remove each position that has another candidate inside [`CANDIDATE_MIN_SEPARATION_BP`].
///
/// The function removes the full group, and not only the extra positions. When some calls come from
/// one mapping event, no rule can select one of them as the real mutation.
fn drop_clustered(positions: &BTreeSet<i64>) -> BTreeSet<i64> {
    let ordered: Vec<i64> = positions.iter().copied().collect();
    ordered
        .iter()
        .enumerate()
        .filter(|(i, &p)| {
            let prev_far = *i == 0 || p - ordered[i - 1] > CANDIDATE_MIN_SEPARATION_BP;
            let next_far = *i + 1 == ordered.len() || ordered[i + 1] - p > CANDIDATE_MIN_SEPARATION_BP;
            prev_far && next_far
        })
        .map(|(_, &p)| p)
        .collect()
}

/// The share of the members with private-Y data above which the code treats a position as
/// **shared by the population**, and not as private.
///
/// A variant that most of a cohort carries did not start on one branch of that cohort. On
/// R1b-CTS4466Plus, *all* 111 donors with private-Y data carried five positions. Those positions
/// are differences between the reference and the population. They are real, but they are not
/// private.
///
/// The blocklist in the application bundle does not name them. That list comes from a CHM13 cohort
/// of 3,352 samples, and that cohort is older than this collection. A rule that reads the cohort in
/// the workspace finds what a fixed list can not.
const COHORT_SHARED_FRACTION: f64 = 0.25;

/// The count of donors that the frequency rule needs before it applies.
///
/// The rule examines what a *large* population shares. A candidate branch needs two carriers, by
/// definition. So in a small cohort, two carriers are already a large share. With four donors, each
/// true branch goes past a limit of 25% and the code removes it.
///
/// Below this count of donors there is no population for the rule to examine. So the rule does
/// nothing, and it makes no estimate.
const COHORT_SHARED_MIN_DONORS: usize = 20;

/// Positions carried by more than [`COHORT_SHARED_FRACTION`] of the members that have private-Y.
fn population_shared_positions(blocks: &[Block], private: &HashMap<SampleGuid, PrivateBucket>) -> BTreeSet<i64> {
    let mut carriers: HashMap<i64, usize> = HashMap::new();
    let mut donors = 0usize;
    for block in blocks {
        for m in &block.members {
            let Some(bucket) = private.get(&m.guid) else { continue };
            donors += 1;
            for pos in candidate_positions(bucket) {
                *carriers.entry(pos).or_default() += 1;
            }
        }
    }
    if donors < COHORT_SHARED_MIN_DONORS {
        return BTreeSet::new();
    }
    let ceiling = (donors as f64 * COHORT_SHARED_FRACTION).ceil() as usize;
    carriers
        .into_iter()
        .filter(|&(_, n)| n > ceiling)
        .map(|(pos, _)| pos)
        .collect()
}

/// The window in which the code treats two positions that define a candidate as one mapping event.
///
/// This window is wider than [`CANDIDATE_MIN_SEPARATION_BP`], which applies to one donor, because
/// the question is different. That rule asks whether the calls of one donor come from one read.
/// This rule asks whether separate branches, with different member sets, are very close together.
///
/// On R1b-CTS4466Plus, three of nine candidates were inside a **window of 567 bp** at 56.83 Mb.
/// Three separate lineage events inside one kilobase do not occur. One repeat unit that the mapper
/// places wrongly across some donors does occur.
///
/// A kilobase is the length of a sequencing fragment. So it is the scale at which one mapping fault
/// can give branches that look separate.
const CANDIDATE_CLUSTER_WINDOW_BP: i64 = 1_000;

/// Every position that would define a candidate branch, anywhere in the cohort.
fn candidate_defining_positions(blocks: &[Block], private: &HashMap<SampleGuid, PrivateBucket>) -> BTreeSet<i64> {
    let mut out = BTreeSet::new();
    for block in blocks {
        if block.members.len() < 2 {
            continue;
        }
        let mut carriers: HashMap<i64, usize> = HashMap::new();
        for m in &block.members {
            let Some(bucket) = private.get(&m.guid) else { continue };
            for pos in candidate_positions(bucket) {
                *carriers.entry(pos).or_default() += 1;
            }
        }
        out.extend(carriers.into_iter().filter(|&(_, n)| n >= 2).map(|(pos, _)| pos));
    }
    out
}

/// The positions that define a candidate and that are inside [`CANDIDATE_CLUSTER_WINDOW_BP`] of
/// another such position.
///
/// The code refuses the full group, as the rule for one donor does. When some branches come from one
/// mapping fault, no rule can select one of them as real.
fn clustered_candidate_positions(blocks: &[Block], private: &HashMap<SampleGuid, PrivateBucket>) -> BTreeSet<i64> {
    let defining: Vec<i64> = candidate_defining_positions(blocks, private).into_iter().collect();
    defining
        .iter()
        .enumerate()
        .filter(|(i, &p)| {
            let near_prev = *i > 0 && p - defining[i - 1] <= CANDIDATE_CLUSTER_WINDOW_BP;
            let near_next = *i + 1 < defining.len() && defining[i + 1] - p <= CANDIDATE_CLUSTER_WINDOW_BP;
            near_prev || near_next
        })
        .map(|(_, &p)| p)
        .collect()
}

/// The positions that would define a candidate branch below **more than one** named block.
///
/// A variant that defines a branch below two different parents did not occur one time. It occurs
/// again in the tree, or the caller makes the same error at that site. In each case it is the one
/// thing that a *new-branch* candidate must not be.
///
/// The laminar test can not find this state, because it examines one block only. So this function
/// finds a repeat across two blocks, before the code accepts any group.
fn recurrent_positions(blocks: &[Block], private: &HashMap<SampleGuid, PrivateBucket>) -> BTreeSet<i64> {
    let mut blocks_per_position: HashMap<i64, BTreeSet<i64>> = HashMap::new();
    for block in blocks {
        if block.members.len() < 2 {
            continue;
        }
        let mut carriers: HashMap<i64, usize> = HashMap::new();
        for m in &block.members {
            let Some(bucket) = private.get(&m.guid) else { continue };
            for pos in candidate_positions(bucket) {
                *carriers.entry(pos).or_default() += 1;
            }
        }
        for (pos, n) in carriers {
            if n >= 2 {
                blocks_per_position.entry(pos).or_default().insert(block.node_id);
            }
        }
    }
    blocks_per_position
        .into_iter()
        .filter(|(_, blocks)| blocks.len() > 1)
        .map(|(pos, _)| pos)
        .collect()
}

/// Add the **candidate branches**. For each block, the function groups the members by the private
/// variants that they share. It then adds a new child block below each group of two members or
/// more.
///
/// Variants that exactly the same set of members carries are equivalent inside this cohort. The
/// SNPs of a named node form a block for the same reason. So each distinct member set becomes one
/// candidate block, and that block carries each of those variants.
///
/// The function accepts a group at once, and it takes the largest group first. It accepts a group
/// only while the set stays **laminar**. Two accepted sets must share no member, or one set must
/// hold the other.
///
/// A set that shares only some members with an accepted set is a conflict. The cause is a recurrent
/// call, or a real disagreement in the phylogeny. The function counts such a set. It does not force
/// the set into a shape that does not fit. The function returns the blocks and that count of
/// conflicts.
///
/// Each member then goes to the smallest accepted set that holds it. That set is unique, because a
/// laminar family of sets is a tree. A member in no set stays on its named block.
pub(crate) fn insert_candidate_branches(
    blocks: Vec<Block>,
    private: &HashMap<SampleGuid, PrivateBucket>,
) -> (Vec<Block>, usize, usize) {
    // The code calculates this set across each block, before it accepts any group. A position
    // that defines a branch below two parents fails everywhere. It does not fail only at the
    // second place where the code reads it.
    let recurrent: BTreeSet<i64> = recurrent_positions(&blocks, private)
        .into_iter()
        .chain(population_shared_positions(&blocks, private))
        .chain(clustered_candidate_positions(&blocks, private))
        .collect();
    let mut out: Vec<Block> = Vec::with_capacity(blocks.len());
    let mut conflicts = 0;
    let mut next_id: i64 = -1;

    for mut block in blocks {
        if block.members.len() < 2 {
            out.push(block);
            continue;
        }
        // A map from a position to the members of *this block* that hold it. The key of a member
        // is its index in `block.members`.
        let mut carriers: HashMap<i64, BTreeSet<usize>> = HashMap::new();
        for (i, m) in block.members.iter().enumerate() {
            let Some(bucket) = private.get(&m.guid) else { continue };
            for pos in candidate_positions(bucket) {
                if recurrent.contains(&pos) {
                    continue;
                }
                carriers.entry(pos).or_default().insert(i);
            }
        }
        // Equivalent variants = same carrier set. Sets of one carrier are private to that member and
        // define no branch.
        let mut groups: BTreeMap<BTreeSet<usize>, Vec<i64>> = BTreeMap::new();
        for (pos, set) in carriers {
            if set.len() >= 2 {
                groups.entry(set).or_default().push(pos);
            }
        }
        if groups.is_empty() {
            out.push(block);
            continue;
        }

        // Take the largest group first, so the code accepts a wide branch before the finer
        // branches inside it. Two groups of the same size compare on the count of shared variants,
        // and then on the lowest position. So the order is always the same.
        let mut ordered: Vec<(BTreeSet<usize>, Vec<i64>)> = groups.into_iter().collect();
        for (_, positions) in &mut ordered {
            positions.sort_unstable();
        }
        ordered.sort_by(|a, b| {
            b.0.len()
                .cmp(&a.0.len())
                .then(b.1.len().cmp(&a.1.len()))
                .then(a.1.first().cmp(&b.1.first()))
        });

        let mut accepted: Vec<(BTreeSet<usize>, Vec<i64>, i64)> = Vec::new(); // (members, positions, id)
        for (set, positions) in ordered {
            let laminar = accepted.iter().all(|(other, _, _)| {
                let shared = set.intersection(other).count();
                shared == 0 || shared == set.len() || shared == other.len()
            });
            if !laminar {
                conflicts += 1;
                continue;
            }
            accepted.push((set, positions, next_id));
            next_id -= 1;
        }

        // Parent = the smallest accepted strict superset; otherwise the named block itself.
        let parent_of = |i: usize, accepted: &[(BTreeSet<usize>, Vec<i64>, i64)]| -> Option<i64> {
            accepted
                .iter()
                .enumerate()
                .filter(|(j, (other, _, _))| {
                    *j != i && other.len() > accepted[i].0.len() && accepted[i].0.is_subset(other)
                })
                .min_by_key(|(_, (other, _, _))| other.len())
                .map(|(_, (_, _, id))| *id)
        };

        // Each member goes to the smallest accepted set that holds it.
        let owner: HashMap<usize, i64> = (0..block.members.len())
            .filter_map(|m| {
                accepted
                    .iter()
                    .filter(|(set, _, _)| set.contains(&m))
                    .min_by_key(|(set, _, _)| set.len())
                    .map(|(_, _, id)| (m, *id))
            })
            .collect();

        // Look up each shared position's variant for its alleles; any carrier's bucket will do,
        // since they all carry the same call at that position.
        let variant_at = |pos: i64, members: &BTreeSet<usize>| -> Option<PrivateVariant> {
            members.iter().find_map(|&i| {
                private
                    .get(&block.members[i].guid)?
                    .variants
                    .iter()
                    .find(|v| v.position == pos)
                    .cloned()
            })
        };

        let mut kids: Vec<Block> = Vec::with_capacity(accepted.len());
        for (i, (set, positions, id)) in accepted.iter().enumerate() {
            let parent = parent_of(i, &accepted).unwrap_or(block.node_id);
            let members: Vec<BlockMember> = set
                .iter()
                .filter(|m| owner.get(m) == Some(id))
                .map(|&m| block.members[m].clone())
                .collect();
            let evidence: Vec<CandidateEvidence> = set
                .iter()
                .flat_map(|&m| {
                    let member = &block.members[m];
                    let bucket = private.get(&member.guid);
                    positions.iter().filter_map(move |&pos| {
                        let v = bucket?.variants.iter().find(|v| v.position == pos)?;
                        Some(CandidateEvidence {
                            guid: member.guid,
                            member: member.name.clone(),
                            position: pos,
                            reference: v.reference,
                            alternate: v.alternate,
                            depth: v.depth,
                            alt_depth: v.alt_depth,
                            allele_fraction: v.allele_fraction,
                            publishable: PublishGate::default().admits(v),
                        })
                    })
                })
                .collect();
            kids.push(Block {
                node_id: *id,
                name: String::new(), // the view localizes a candidate's label
                parent: Some(parent),
                depth: 0, // set below, once the parent chain is known
                loci: positions
                    .iter()
                    .filter_map(|&p| variant_at(p, set))
                    .map(|v| private_locus(&v))
                    .collect(),
                subtree_members: set.len(),
                members,
                collapsed: Vec::new(),
                candidate: true,
                evidence,
            });
        }
        // Depth: walk up the synthetic parents to the named block.
        let depth_of: HashMap<i64, usize> = {
            let by_id: HashMap<i64, &Block> = kids.iter().map(|k| (k.node_id, k)).collect();
            kids.iter()
                .map(|k| {
                    let mut d = block.depth + 1;
                    let mut cur = k.parent;
                    while let Some(p) = cur.filter(|p| *p < 0) {
                        d += 1;
                        cur = by_id.get(&p).and_then(|b| b.parent);
                    }
                    (k.node_id, d)
                })
                .collect()
        };
        for k in &mut kids {
            k.depth = depth_of[&k.node_id];
        }
        // Pre-order: shallowest first, so a parent always precedes its children.
        kids.sort_by_key(|k| (k.depth, -k.node_id));

        // Members that joined a candidate branch leave the named block; the rest stay.
        block.members = block
            .members
            .iter()
            .enumerate()
            .filter(|(i, _)| !owner.contains_key(i))
            .map(|(_, m)| m.clone())
            .collect();
        out.push(block);
        out.extend(kids);
    }
    (out, conflicts, recurrent.len())
}

/// Fold runs of **member-less single-child** branches into the branch below them, when the run is at
/// least `min_run` long.
///
/// A subtree over a deep haplotree holds mostly such chains. Those are the branches that no member
/// sits on, and that divide nothing inside this cohort.
///
/// The join is not only a change to the display. Inside this cohort those branches are one block
/// that nothing divides. So the loci of an absorbed branch go to the branch that stays, with the
/// loci nearest the root first. [`Block::collapsed`] keeps the names of the absorbed branches.
///
/// The function does not change `subtree_members`. An absorbed node has no member of its own and
/// has one child. So its count is already the count of the branch that stays.
///
/// The function is pure. It reads no tree and does no I/O, so a unit test can call it with a
/// `Vec<Block>` that the test builds.
pub(crate) fn collapse_blocks(blocks: Vec<Block>, min_run: usize) -> Vec<Block> {
    if blocks.is_empty() || min_run == 0 {
        return blocks;
    }
    let mut children: HashMap<i64, Vec<i64>> = HashMap::new();
    for b in &blocks {
        if let Some(p) = b.parent {
            children.entry(p).or_default().push(b.node_id);
        }
    }
    let only_child = |id: i64| match children.get(&id) {
        Some(c) if c.len() == 1 => Some(c[0]),
        _ => None,
    };
    let absorbable = |b: &Block| b.members.is_empty() && only_child(b.node_id).is_some();

    let mut by_id: HashMap<i64, Block> = blocks.iter().map(|b| (b.node_id, b.clone())).collect();
    let mut absorbed: HashSet<i64> = HashSet::new();

    // The `blocks` list is in pre-order, so the code always reaches the node of a run that is
    // nearest the root first. A node whose parent the code can also absorb is in the middle of a
    // run. So the head of that run already covers it.
    for b in &blocks {
        if !absorbable(b) {
            continue;
        }
        if b.parent.and_then(|p| by_id.get(&p)).is_some_and(absorbable) {
            continue;
        }
        // Walk down the chain of absorbable nodes; `cur` ends on the run's last one.
        let mut run: Vec<i64> = Vec::new();
        let mut cur = b.node_id;
        loop {
            run.push(cur);
            match only_child(cur).and_then(|c| by_id.get(&c)) {
                Some(n) if absorbable(n) => cur = n.node_id,
                _ => break,
            }
        }
        if run.len() < min_run {
            continue;
        }
        let Some(target_id) = only_child(cur) else { continue };

        // Merge root-most first, so the survivor's loci read in descent order.
        let mut folded_loci = Vec::new();
        let mut folded_names = Vec::new();
        for id in &run {
            let n = &by_id[id];
            folded_loci.extend(n.loci.iter().cloned());
            folded_names.extend(n.collapsed.iter().cloned());
            folded_names.push(n.name.clone());
            absorbed.insert(*id);
        }
        let head_parent = by_id[&run[0]].parent;
        let target = by_id
            .get_mut(&target_id)
            .expect("a child of an induced node is induced");
        folded_loci.append(&mut target.loci);
        target.loci = folded_loci;
        folded_names.append(&mut target.collapsed);
        target.collapsed = folded_names;
        target.parent = head_parent;
    }

    // Write the list again in the first pre-order, with no absorbed node. The code calculates each
    // depth against the parents that stay. A change to a link only moves a node *up* to an
    // ancestor, so the list stays in pre-order.
    let mut depth: HashMap<i64, usize> = HashMap::new();
    let mut out = Vec::with_capacity(blocks.len() - absorbed.len());
    for b in &blocks {
        if absorbed.contains(&b.node_id) {
            continue;
        }
        let mut b = by_id
            .remove(&b.node_id)
            .expect("each surviving node is emitted exactly once");
        b.depth = b.parent.and_then(|p| depth.get(&p)).map_or(0, |d| d + 1);
        depth.insert(b.node_id, b.depth);
        out.push(b);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus(name: &str, position: i64) -> Locus {
        Locus {
            position,
            ancestral: "A".into(),
            derived: "G".into(),
            name: name.into(),
        }
    }

    fn member(name: &str) -> BlockMember {
        BlockMember {
            guid: SampleGuid(uuid::Uuid::new_v4()),
            name: name.into(),
            private_novel: None,
            private_publishable: None,
            private_total: None,
        }
    }

    /// One block. `parent` of `0` means "root" (no parent), since node ids here start at 1.
    fn block(id: i64, name: &str, parent: i64, members: &[&str]) -> Block {
        Block {
            node_id: id,
            name: name.into(),
            parent: (parent != 0).then_some(parent),
            depth: 0,
            loci: vec![locus(&format!("M{id}"), id * 100)],
            members: members.iter().map(|m| member(m)).collect(),
            subtree_members: 0,
            collapsed: Vec::new(),
            candidate: false,
            evidence: Vec::new(),
        }
    }

    fn pvar(position: i64, class: PrivateClass, region: Option<YRegionClass>) -> PrivateVariant {
        PrivateVariant {
            position,
            reference: 'A',
            alternate: 'G',
            depth: 30,
            alt_depth: 30,
            allele_fraction: 1.0,
            class,
            region,
        }
    }

    /// A bucket of high-confidence (novel, unique-sequence) private calls at `positions`.
    fn bucket(positions: &[i64]) -> PrivateBucket {
        PrivateBucket {
            terminal: "R-X".into(),
            variants: positions.iter().map(|&p| pvar(p, PrivateClass::Novel, None)).collect(),
        }
    }

    /// Map each of `block`'s members, in order, to a private bucket.
    fn privates(block: &Block, buckets: &[PrivateBucket]) -> HashMap<SampleGuid, PrivateBucket> {
        block
            .members
            .iter()
            .zip(buckets)
            .map(|(m, b)| (m.guid, b.clone()))
            .collect()
    }

    /// ```text
    /// root ─┬─> A ──> B ──> C(kane)
    ///       └─> D(smith)
    /// ```
    ///
    /// A run of two branches, `A` and `B`, above a branch with a member. Each branch in the run
    /// has one child and no member.
    ///
    /// The member `D` keeps `root` a point where the tree divides. So the head of the run is `A`.
    /// Without `D`, the code could also absorb `root`, and the full line would become one block.
    /// That result is correct, and `collapse_stops_a_run_at_a_placed_branch` covers it.
    fn chain() -> Vec<Block> {
        let mut b = vec![
            block(1, "root", 0, &[]),
            block(2, "A", 1, &[]),
            block(3, "B", 2, &[]),
            block(4, "C", 3, &["kane"]),
            block(5, "D", 1, &["smith"]),
        ];
        roll_up_subtree_members(&mut b);
        b
    }

    #[test]
    fn roll_up_counts_members_at_and_below() {
        // root → A ─┬→ B(x, y)
        //           └→ C(z)
        let mut blocks = vec![
            block(1, "root", 0, &[]),
            block(2, "A", 1, &[]),
            block(3, "B", 2, &["x", "y"]),
            block(4, "C", 2, &["z"]),
        ];
        roll_up_subtree_members(&mut blocks);
        let n = |name: &str| blocks.iter().find(|b| b.name == name).unwrap().subtree_members;
        assert_eq!(n("root"), 3);
        assert_eq!(n("A"), 3);
        assert_eq!(n("B"), 2);
        assert_eq!(n("C"), 1);
    }

    #[test]
    fn collapse_folds_a_run_into_the_branch_below_it() {
        let out = collapse_blocks(chain(), 2);
        let names: Vec<&str> = out.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, vec!["root", "C", "D"], "A and B should have been folded into C");

        let c = out.iter().find(|b| b.name == "C").unwrap();
        assert_eq!(c.collapsed, vec!["A".to_string(), "B".to_string()], "root-most first");
        assert_eq!(c.parent, Some(1), "C should re-attach to the run head's parent");
        // Loci read root-most first: A's, then B's, then C's own.
        let markers: Vec<&str> = c.loci.iter().map(|l| l.name.as_str()).collect();
        assert_eq!(markers, vec!["M2", "M3", "M4"]);
    }

    #[test]
    fn collapse_respects_min_run() {
        // With a limit of 3, the code does not change the run of 2.
        let out = collapse_blocks(chain(), 3);
        let names: Vec<&str> = out.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, vec!["root", "A", "B", "C", "D"]);
        assert!(out.iter().all(|b| b.collapsed.is_empty()));
    }

    #[test]
    fn collapse_recomputes_depth_against_surviving_parents() {
        let out = collapse_blocks(chain(), 2);
        let d = |name: &str| out.iter().find(|b| b.name == name).unwrap().depth;
        assert_eq!(d("root"), 0);
        assert_eq!(d("C"), 1, "C sat at depth 3 before its two ancestors folded away");
    }

    #[test]
    fn collapse_preserves_subtree_member_counts() {
        let before: HashMap<String, usize> = chain().iter().map(|b| (b.name.clone(), b.subtree_members)).collect();
        for b in collapse_blocks(chain(), 2) {
            assert_eq!(b.subtree_members, before[&b.name], "{} changed count", b.name);
        }
    }

    #[test]
    fn collapse_never_absorbs_a_branch_point_or_a_placed_branch() {
        // root → A ─┬→ B(x)      A branches, so it is not a single-child run;
        //           └→ C(y)      B and C carry members and are leaves.
        let mut blocks = vec![
            block(1, "root", 0, &[]),
            block(2, "A", 1, &[]),
            block(3, "B", 2, &["x"]),
            block(4, "C", 2, &["y"]),
        ];
        roll_up_subtree_members(&mut blocks);
        let out = collapse_blocks(blocks, 2);
        let names: Vec<&str> = out.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, vec!["root", "A", "B", "C"], "a branch point must survive");
    }

    #[test]
    fn collapse_stops_a_run_at_a_placed_branch() {
        // root → A → B(x) → C(y): B carries a member, so the run above it is just [root, A]
        // (length 2) and B itself is never absorbed.
        let mut blocks = vec![
            block(1, "root", 0, &[]),
            block(2, "A", 1, &[]),
            block(3, "B", 2, &["x"]),
            block(4, "C", 3, &["y"]),
        ];
        roll_up_subtree_members(&mut blocks);
        let out = collapse_blocks(blocks, 2);
        let names: Vec<&str> = out.iter().map(|b| b.name.as_str()).collect();
        assert_eq!(names, vec!["B", "C"], "root and A fold into B; B keeps its member");
        let b = out.iter().find(|b| b.name == "B").unwrap();
        assert_eq!(b.collapsed, vec!["root".to_string(), "A".to_string()]);
        assert_eq!(b.parent, None, "B becomes the new root");
        assert_eq!(b.depth, 0);
        assert_eq!(out.iter().find(|b| b.name == "C").unwrap().depth, 1);
    }

    // ---- candidate branches (phase 3) ----------------------------------------

    #[test]
    fn two_members_sharing_a_novel_variant_become_a_candidate_branch() {
        let b = block(1, "R-X", 0, &["kane", "smith", "jones"]);
        // kane + smith share 100 and 200; jones shares nothing.
        let p = privates(
            &b,
            &[
                bucket(&[100_000, 200_000]),
                bucket(&[100_000, 200_000]),
                bucket(&[900_000]),
            ],
        );
        let (out, conflicts, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(conflicts, 0);
        assert_eq!(out.len(), 2, "the named block plus one candidate");

        let cand = out.iter().find(|x| x.candidate).unwrap();
        assert!(cand.node_id < 0, "candidate ids are synthetic and negative");
        assert!(cand.name.is_empty(), "the view localizes a candidate's label");
        assert_eq!(cand.parent, Some(1));
        assert_eq!(cand.depth, 1);
        // The two shared positions are equivalent on this branch. They give one block with two
        // loci.
        let mut pos: Vec<i64> = cand.loci.iter().map(|l| l.position).collect();
        pos.sort_unstable();
        assert_eq!(pos, vec![100_000, 200_000]);

        let mut names: Vec<&str> = cand.members.iter().map(|m| m.name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(names, vec!["kane", "smith"]);
        // The two who moved down leave the named block; the third stays.
        let named = out.iter().find(|x| !x.candidate).unwrap();
        assert_eq!(
            named.members.iter().map(|m| m.name.as_str()).collect::<Vec<_>>(),
            vec!["jones"]
        );
    }

    #[test]
    fn a_variant_only_one_member_carries_defines_no_branch() {
        let b = block(1, "R-X", 0, &["kane", "smith"]);
        let p = privates(&b, &[bucket(&[100]), bucket(&[200])]);
        let (out, _, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(out.len(), 1, "no shared variant → no candidate");
        assert!(!out[0].candidate);
        assert_eq!(out[0].members.len(), 2, "members stay on the named block");
    }

    #[test]
    fn structural_region_and_off_path_calls_never_form_a_branch() {
        let b = block(1, "R-X", 0, &["kane", "smith"]);
        // The two men "share" a call in a palindrome and a known SNP that is off the path.
        // Neither call shows a shared ancestor. The first is a paralog artefact. The second marks
        // a branch that already exists.
        let shared = PrivateBucket {
            terminal: "R-X".into(),
            variants: vec![
                pvar(100, PrivateClass::Novel, Some(YRegionClass::Palindrome)),
                pvar(200, PrivateClass::OffPathKnown("M269".into()), None),
            ],
        };
        let p = privates(&b, &[shared.clone(), shared]);
        let (out, _, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(out.len(), 1, "noise must not manufacture a branch");
    }

    #[test]
    fn nested_sharing_nests_the_candidate_branches() {
        let b = block(1, "R-X", 0, &["a", "b", "c"]);
        // Each of the three members shares position 100. Members a and b also share position 200,
        // which gives a finer branch inside the wider one.
        let p = privates(
            &b,
            &[
                bucket(&[100_000, 200_000]),
                bucket(&[100_000, 200_000]),
                bucket(&[100_000]),
            ],
        );
        let (out, conflicts, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(conflicts, 0);

        let cands: Vec<&Block> = out.iter().filter(|x| x.candidate).collect();
        assert_eq!(cands.len(), 2);
        let broad = cands.iter().find(|x| x.depth == 1).unwrap();
        let fine = cands.iter().find(|x| x.depth == 2).unwrap();
        assert_eq!(broad.parent, Some(1), "the broad branch hangs off the named block");
        assert_eq!(fine.parent, Some(broad.node_id), "the finer branch nests inside it");
        assert_eq!(broad.loci[0].position, 100_000);
        assert_eq!(fine.loci[0].position, 200_000);
        // c sits on the broad branch; a and b moved down to the finer one.
        assert_eq!(
            broad.members.iter().map(|m| m.name.as_str()).collect::<Vec<_>>(),
            vec!["c"]
        );
        let mut fine_names: Vec<&str> = fine.members.iter().map(|m| m.name.as_str()).collect();
        fine_names.sort_unstable();
        assert_eq!(fine_names, vec!["a", "b"]);
        // The list is in pre-order, so a parent comes before its child.
        let pos = |id: i64| out.iter().position(|x| x.node_id == id).unwrap();
        assert!(pos(broad.node_id) < pos(fine.node_id));
    }

    #[test]
    fn overlapping_non_nested_sharing_is_counted_as_a_conflict_not_forced() {
        let b = block(1, "R-X", 0, &["a", "b", "c"]);
        // The set {a,b} shares position 100, and the set {b,c} shares position 200. Neither set
        // holds the other. So one tree can not hold both as a branch. The code removes the set with
        // the lower rank and counts it.
        let p = privates(
            &b,
            &[bucket(&[100_000]), bucket(&[100_000, 200_000]), bucket(&[200_000])],
        );
        let (out, conflicts, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(conflicts, 1, "the conflicting group is reported, not silently dropped");
        assert_eq!(
            out.iter().filter(|x| x.candidate).count(),
            1,
            "only the laminar one survives"
        );
    }

    #[test]
    fn a_member_with_no_computed_private_y_is_simply_not_grouped() {
        let b = block(1, "R-X", 0, &["kane", "smith"]);
        // Only `kane` has a bucket. No analysis ran for `smith`, and that state is not the same
        // as a subject with no private variant.
        let p: HashMap<SampleGuid, PrivateBucket> = [(b.members[0].guid, bucket(&[100]))].into_iter().collect();
        let (out, _, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].members.len(), 2);
    }

    // ---- artefact filters ------------------------------------------------------

    #[test]
    fn clustered_calls_are_dropped_whole() {
        // Six "new" calls inside 32 bp come from one read that the mapper placed wrongly. They are
        // not six mutations. This shape produced most of the first candidate branches on the
        // CTS4466 cohort.
        let cluster = bucket(&[16342231, 16342238, 16342245, 16342253, 16342258, 16342263]);
        assert!(
            candidate_positions(&cluster).is_empty(),
            "the whole cluster goes; there is no basis for electing one call the real one"
        );

        // A lone call keeps its place, and a distant neighbour does not drag it down.
        let spread = bucket(&[1_000_000, 2_000_000, 2_000_050]);
        let kept: Vec<i64> = candidate_positions(&spread).into_iter().collect();
        assert_eq!(kept, vec![1_000_000], "only the pair within 100 bp is dropped");
    }

    #[test]
    fn a_position_defining_branches_under_two_parents_is_rejected() {
        // Two members below one block shared position 11311865, and two members below another
        // block also shared it. A variant that occurred two times can not mark a new branch. The
        // laminar test can not find this state, because it examines one block only.
        let mut left = block(1, "R-A", 0, &["a", "b"]);
        left.subtree_members = 2;
        let mut right = block(2, "R-B", 0, &["c", "d"]);
        right.subtree_members = 2;

        let mut private = HashMap::new();
        for m in left.members.iter().chain(right.members.iter()) {
            private.insert(m.guid, bucket(&[11311865]));
        }
        let (out, _, recurrent) = insert_candidate_branches(vec![left, right], &private);
        assert_eq!(recurrent, 1, "the shared position is counted as recurrent");
        assert_eq!(
            out.iter().filter(|b| b.candidate).count(),
            0,
            "and defines no candidate under either parent"
        );
    }

    #[test]
    fn a_position_confined_to_one_parent_still_defines_a_branch() {
        // The guard must not reject an ordinary shared variant just because two blocks exist.
        let mut left = block(1, "R-A", 0, &["a", "b"]);
        left.subtree_members = 2;
        let mut right = block(2, "R-B", 0, &["c", "d"]);
        right.subtree_members = 2;

        let mut private = HashMap::new();
        for m in &left.members {
            private.insert(m.guid, bucket(&[500_000]));
        }
        for m in &right.members {
            private.insert(m.guid, bucket(&[900_000]));
        }
        let (out, _, recurrent) = insert_candidate_branches(vec![left, right], &private);
        assert_eq!(recurrent, 0);
        assert_eq!(out.iter().filter(|b| b.candidate).count(), 2, "one candidate per block");
    }

    #[test]
    fn candidates_clustering_across_the_cohort_are_all_rejected() {
        // The case at 56.83 Mb. Three branches with *different* member sets are inside 567 bp.
        // Three separate lineage events inside one kilobase do not occur. One repeat unit that the
        // mapper places wrongly across some donors does occur.
        let mut a = block(1, "R-A", 0, &["a", "b"]);
        a.subtree_members = 2;
        let mut b = block(2, "R-B", 0, &["c", "d"]);
        b.subtree_members = 2;
        let mut private = HashMap::new();
        for m in &a.members {
            private.insert(m.guid, bucket(&[56_832_495]));
        }
        for m in &b.members {
            private.insert(m.guid, bucket(&[56_833_062]));
        }
        let (out, _, dropped) = insert_candidate_branches(vec![a, b], &private);
        assert_eq!(dropped, 2, "both positions go — neither can be elected the real one");
        assert_eq!(out.iter().filter(|x| x.candidate).count(), 0);
    }

    #[test]
    fn candidates_far_apart_are_left_alone() {
        // Same shape, megabases apart: two ordinary branches, and the rule must not touch them.
        let mut a = block(1, "R-A", 0, &["a", "b"]);
        a.subtree_members = 2;
        let mut b = block(2, "R-B", 0, &["c", "d"]);
        b.subtree_members = 2;
        let mut private = HashMap::new();
        for m in &a.members {
            private.insert(m.guid, bucket(&[10_756_695]));
        }
        for m in &b.members {
            private.insert(m.guid, bucket(&[18_055_974]));
        }
        let (out, _, dropped) = insert_candidate_branches(vec![a, b], &private);
        assert_eq!(dropped, 0);
        assert_eq!(out.iter().filter(|x| x.candidate).count(), 2);
    }

    // ---- export ---------------------------------------------------------------

    /// A tree with two blocks, one candidate branch, and one member with no placement. This shape
    /// covers each column that the export must write correctly.
    fn exportable() -> ProjectBlockTree {
        let mut named = block(1, "R-X", 0, &["kane"]);
        named.subtree_members = 3;
        let mut cand = block(-1, "", 1, &["a", "b"]);
        cand.name = String::new();
        cand.candidate = true;
        cand.depth = 1;
        cand.subtree_members = 2;
        cand.loci = vec![locus("chrY:100", 100)];
        ProjectBlockTree {
            dna: DnaType::Y,
            blocks: vec![named, cand],
            unplaced: vec![
                UnplacedMember {
                    guid: SampleGuid(uuid::Uuid::new_v4()),
                    name: "nolabel".into(),
                    terminal: None,
                },
                UnplacedMember {
                    guid: SampleGuid(uuid::Uuid::new_v4()),
                    name: "skewed".into(),
                    terminal: Some("F-M89".into()),
                },
            ],
            provider: "decodingus".into(),
            build_key: "hs1".into(),
            candidate_conflicts: 2,
            candidate_recurrent: 0,
        }
    }

    #[test]
    fn tsv_export_marks_candidates_and_keeps_the_unplaced() {
        let tsv = crate::export::block_tree_tsv(&exportable());
        assert!(tsv.contains("decodingus"), "header names the tree");
        assert!(tsv.contains("hs1"), "header names the coordinate space");
        assert!(
            tsv.contains("2 shared-variant grouping(s) dropped"),
            "conflicts reported"
        );
        // A candidate must never read as a published haplogroup name.
        assert!(tsv.contains("\tcandidate\tcandidate\t"));
        assert!(tsv.contains("\tbranch\tR-X\t"));
        // Both unplaced members appear, each with why.
        assert!(tsv.contains("nolabel\t\tno placement"));
        assert!(tsv.contains("skewed\tF-M89\tterminal absent from this tree"));
    }

    #[test]
    fn html_export_is_self_contained_and_escapes_names() {
        let mut tree = exportable();
        tree.blocks[0].name = "R-<script>".into();
        let html = crate::export::block_tree_html(&tree, "Smith & Sons");
        assert!(html.starts_with("<!doctype html>"));
        assert!(
            !html.contains("http://") && !html.contains("https://"),
            "no external assets"
        );
        assert!(html.contains("Smith &amp; Sons"), "project name is escaped");
        assert!(!html.contains("<script>"), "block names are escaped");
        assert!(html.contains("candidate branch"));
        assert!(html.contains("Not on this tree (2)"));
    }

    #[test]
    fn collapse_is_a_no_op_on_an_empty_tree_or_a_zero_threshold() {
        assert!(collapse_blocks(Vec::new(), 2).is_empty());
        let names: Vec<String> = collapse_blocks(chain(), 0).iter().map(|b| b.name.clone()).collect();
        assert_eq!(names, vec!["root", "A", "B", "C", "D"]);
    }
}

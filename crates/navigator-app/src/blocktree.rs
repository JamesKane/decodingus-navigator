//! Project **block tree** — the cohort counterpart to the per-subject descent report.
//!
//! Given a project, build the induced subtree of the haplotree spanning its members' terminal
//! haplogroups: every branch that any member lies on, each carrying the run of defining SNPs that
//! are phylogenetically equivalent on it (a *block*), with members hanging off their own terminal.
//! This is the FTDNA "Block Tree" surface, over placements Navigator already computed.
//!
//! Two rules shape everything here:
//!
//! - **It reads placements, never re-places.** Terminals come from `haplogroup_terminals` — the same
//!   reconciliation the subjects table and project report use. Nothing here can move a subject.
//! - **A member that can't be placed is reported, not dropped** ([`UnplacedMember`]). On a multi-lab
//!   cohort provider/build skew is expected; silently omitting those members would make the tree
//!   look like it accounts for the whole project when it doesn't.
//!
//! Design: `documents/design/project-block-tree.md`.

use std::collections::BTreeSet;

use super::*;

/// Collapse a run of member-less single-child branches only when it is at least this long. A lone
/// intermediate branch is worth naming; a run of two or more is noise between the splits the cohort
/// actually resolves.
pub const COLLAPSE_MIN_RUN: usize = 2;

impl App {
    /// Build the [`ProjectBlockTree`] for `project_id`.
    ///
    /// `Ok(None)` only when the project has no members at all. A project whose members are all
    /// unplaced still yields a tree — empty `blocks`, everyone in `unplaced` — because that is a
    /// meaningful answer ("nothing here is placed yet"), and it costs no tree fetch: the multi-MB
    /// download + parse is skipped entirely when no member has a terminal.
    pub async fn project_block_tree(
        &self,
        project_id: i64,
        dna: DnaType,
    ) -> Result<Option<ProjectBlockTree>, AppError> {
        let members = biosample::list_members_for_project(self.store.pool(), project_id).await?;
        if members.is_empty() {
            return Ok(None);
        }

        // One bulk reconciliation for the whole workspace rather than a query per member — the same
        // call `project_report` and `project_str_overview` make.
        let terminals = self.haplogroup_terminals().await?;

        // (guid, display name, terminal) per member, for the requested lineage.
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

        // Cheap first, as `descent_report` does: with nothing placed there is no tree to draw, and
        // fetching + parsing a multi-MB document to discover that is pure waste. Common on a
        // freshly imported project.
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
                // The *configured* provider: no tree was fetched, so no runtime fallback happened
                // either. `build_key` stays empty for the same reason — nothing was parsed.
                provider: match y_tree_provider() {
                    YTreeProvider::DecodingUs => "decodingus".to_string(),
                    YTreeProvider::Ftdna => "ftdna".to_string(),
                },
                build_key: String::new(),
                candidate_conflicts: 0,
                candidate_recurrent: 0,
            }));
        }

        // Provider and build key are whatever the fetch actually resolved to, not what was
        // configured: the mtDNA path falls back to FTDNA at runtime when the DecodingUs tree can't be
        // remapped, and its loci are rCRS either way — not the Y coordinate space.
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
                    // The FTDNA Y tree is published on GRCh38, whatever the members are aligned to.
                    (tree, "ftdna", "GRCh38")
                }
            },
            DnaType::Mt => {
                let (tree, provider) = self.mt_tree_rcrs().await?;
                (tree, provider, "rCRS")
            }
        };

        // One name index for the whole cohort. The per-subject path scans the node map linearly for
        // its single terminal, which would be quadratic here.
        let index = navigator_analysis::haplo::name_index(&tree);
        // Resolve placement first, so the private-Y load below covers only members who actually
        // appear on the tree — on a real cohort that is a fraction of the roster (243 of 1881 on
        // R1b-CTS4466Plus), and an unplaced member's private variants can't be drawn anywhere.
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
        // Stable leaf order, so the layout doesn't reshuffle between opens.
        for b in &mut blocks {
            b.members.sort_by(|x, y| (&x.name, x.guid.0).cmp(&(&y.name, y.guid.0)));
        }
        roll_up_subtree_members(&mut blocks);
        let blocks = collapse_blocks(blocks, COLLAPSE_MIN_RUN);
        // Candidates go in *after* the collapse: they are leaves with members, so they could never
        // be absorbed, and inserting them earlier would only make the collapse reason about
        // synthetic nodes.
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

    /// The one coordinate space the cohort's tree is parsed under: the **modal** DecodingUs build key
    /// across the members' alignments, falling back to `hs1`.
    ///
    /// A cohort spans builds, so there is no per-subject answer as there is in `descent_report`.
    /// Picking one is safe because node names and topology are build-independent — only the loci
    /// *positions* are, and the aggregate carries the key so the view can say which it means. Ties
    /// break on the key name, so the choice doesn't depend on map iteration order.
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

/// Fill `subtree_members` — members at or below each block.
///
/// `blocks` is in pre-order, so walking it **backwards** visits every child before its parent and one
/// pass suffices.
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

/// Synthesize a [`Locus`] standing for a private (unnamed) variant, so a candidate branch's shared
/// variants render through exactly the same path as a named branch's defining SNPs.
fn private_locus(v: &PrivateVariant) -> Locus {
    Locus {
        position: v.position,
        ancestral: v.reference.to_string(),
        derived: v.alternate.to_string(),
        name: format!("chrY:{}", v.position),
    }
}

/// The positions a subject carries as **high-confidence new-branch candidates**: novel (not in the
/// tree at all) *and* in unique sequence.
///
/// Off-path-*known* variants are excluded on purpose — those support an existing finer branch, which
/// is a placement question, not a new one. Structural-region calls are excluded because chrY
/// palindromes and amplicons are paralog-prone: two men "sharing" a call there are far more likely to
/// share a mapping artefact than an ancestor. Sharing noise would manufacture branches, which is the
/// one failure mode this feature must not have.
fn candidate_positions(bucket: &PrivateBucket) -> BTreeSet<i64> {
    let novel: BTreeSet<i64> = bucket
        .variants
        .iter()
        .filter(|v| v.class == PrivateClass::Novel && v.region.is_none())
        .map(|v| v.position)
        .collect();
    drop_clustered(&novel)
}

/// How far apart two novel calls must be to count as independent mutations.
///
/// Real Y mutations are scattered across megabases; a handful of "novel" calls within tens of bases
/// is one misaligned read smearing several false SNVs, which is why the GVCF path already imposes a
/// depth floor for the same reason. On the CTS4466 cohort the first candidate branches were built
/// almost entirely from such clusters — six positions inside 32 bp, gaps of 5–8 bp.
const CANDIDATE_MIN_SEPARATION_BP: i64 = 100;

/// Drop every position that has another candidate within [`CANDIDATE_MIN_SEPARATION_BP`].
///
/// The whole cluster goes, not the extras: when several calls share one mapping event there is no
/// basis for electing one of them the real mutation.
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

/// Share of the cohort's private-Y-bearing members above which a position is treated as
/// **population-shared** rather than private.
///
/// A variant carried by most of a cohort did not arise on one branch of it. On R1b-CTS4466Plus five
/// positions were carried by *all* 111 donors with private-Y — those are reference-vs-population
/// differences, real but not private, and the bundled cohort-shared blocklist (derived from a
/// 3,352-sample CHM13 cohort that predates this collection) does not list them. Deriving the
/// exclusion from the cohort in hand catches what a bundled list cannot anticipate.
const COHORT_SHARED_FRACTION: f64 = 0.25;

/// Donors required before the frequency rule engages at all.
///
/// The rule reasons about what a *large* population shares. A candidate branch needs two carriers by
/// definition, so in a small cohort two carriers are already a large share — at four donors, every
/// genuine branch would exceed a 25% ceiling and be thrown away. Below this many donors there is no
/// population to argue from, so the rule abstains rather than guessing.
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

/// Window within which two *candidate-defining* positions are treated as one mapping event.
///
/// Wider than the per-donor [`CANDIDATE_MIN_SEPARATION_BP`] because this is a different question.
/// That rule asks whether one donor's calls smear across a read; this asks whether separate branches,
/// with different member sets, land suspiciously close together. On R1b-CTS4466Plus three of nine
/// candidates fell inside a **567 bp window** at 56.83 Mb — three independent lineage events in less
/// than a kilobase is not a thing that happens, whereas one repeat unit mis-mapping across several
/// donors is. A kilobase spans a sequencing fragment, so it is the scale at which one mis-mapping
/// event can produce apparently separate branches.
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

/// Candidate-defining positions lying within [`CANDIDATE_CLUSTER_WINDOW_BP`] of another.
///
/// The whole cluster is rejected, as in the per-donor rule: when several apparent branches share one
/// mis-mapping there is no basis for electing one of them real.
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

/// Positions that would define a candidate branch under **more than one** named block.
///
/// A variant defining a branch below two different parents did not arise once: it is recurrent, or a
/// systematic call error. Either way it is the one thing a *new-branch* candidate must not be. The
/// laminar check cannot see this — it reasons within a single block — so cross-block recurrence is
/// caught here, before any group is accepted.
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

/// Insert **candidate branches**: for every block, group its members by the private variants they
/// share, and hang a synthetic child block off each group of two or more.
///
/// Variants shared by exactly the same set of members are phylogenetically equivalent within this
/// cohort — the same reasoning that makes a named node's SNPs a block — so each distinct member set
/// becomes one candidate block carrying all of them.
///
/// Groups are accepted greedily, largest first, and only while they stay **laminar**: any two
/// accepted sets must be disjoint or nested. A set that partly overlaps an accepted one is a
/// conflict (a recurrent call, or real phylogenetic disagreement) and is counted, not forced into a
/// shape it doesn't fit. Returns the blocks plus that conflict count.
///
/// Each member then lands in the smallest accepted set containing it — unambiguous, because a
/// laminar family is a tree — and members in no set stay on their named block.
pub(crate) fn insert_candidate_branches(
    blocks: Vec<Block>,
    private: &HashMap<SampleGuid, PrivateBucket>,
) -> (Vec<Block>, usize, usize) {
    // Computed across all blocks before any group is accepted — a position defining branches under
    // two parents is disqualified everywhere, not just wherever it happens to be seen second.
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
        // position → the members of *this block* carrying it, keyed by index into `block.members`.
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

        // Largest first, so a broader branch is accepted before the finer ones nested inside it.
        // Ties broken deterministically: more shared variants, then lowest position.
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

        // Each member goes to the smallest accepted set containing it.
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
/// An induced subtree over a deep haplotree is mostly such chains: intermediate branches no member
/// sits on that split nothing within this cohort. Merging them is not a display trick — within this
/// cohort those branches genuinely are one undivided block, so the absorbed loci join the survivor's
/// own (root-most first) and the absorbed names are kept in [`Block::collapsed`].
///
/// `subtree_members` survives untouched: an absorbed node has no members of its own and exactly one
/// child, so its count already equals the survivor's.
///
/// Pure — no tree, no I/O — so it is unit-testable on a hand-built `Vec<Block>`.
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

    // `blocks` is pre-order, so a run's root-most node is always reached first. A node whose parent
    // is itself absorbable is therefore mid-run and already handled by the run's head.
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

    // Re-emit in the original pre-order minus the absorbed nodes, depth recomputed against the
    // surviving parents. Rewiring only ever moves a node *up* to an ancestor, so pre-order holds.
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
    /// A run of exactly two member-less single-child branches (`A`, `B`) above a placed one. `D`
    /// keeps `root` a branch point, so the run's head is `A` — otherwise `root` would be absorbable
    /// too and the whole spine would fold (which is correct, and is what
    /// `collapse_stops_a_run_at_a_placed_branch` covers).
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
        // With a threshold of 3, the run of 2 is left alone.
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
        // Both shared positions are equivalent on this branch — one block, two loci.
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
        // Both men "share" a palindrome call and a known off-path SNP. Neither is evidence of a
        // shared ancestor — the first is a paralog artefact, the second an existing branch.
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
        // All three share 100; a and b additionally share 200 — a finer branch inside the broader one.
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
        // Pre-order: a parent is emitted before its child.
        let pos = |id: i64| out.iter().position(|x| x.node_id == id).unwrap();
        assert!(pos(broad.node_id) < pos(fine.node_id));
    }

    #[test]
    fn overlapping_non_nested_sharing_is_counted_as_a_conflict_not_forced() {
        let b = block(1, "R-X", 0, &["a", "b", "c"]);
        // {a,b} share 100; {b,c} share 200. Neither set contains the other, so they cannot both be
        // branches of one tree — the smaller-ranked one is dropped and counted.
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
        // Only `kane` has a bucket at all — `smith` was never analyzed, which is not the same as
        // having no private variants.
        let p: HashMap<SampleGuid, PrivateBucket> = [(b.members[0].guid, bucket(&[100]))].into_iter().collect();
        let (out, _, _) = insert_candidate_branches(vec![b], &p);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].members.len(), 2);
    }

    // ---- artefact filters ------------------------------------------------------

    #[test]
    fn clustered_calls_are_dropped_whole() {
        // Six "novel" calls inside 32 bp is one misaligned read, not six mutations — the shape that
        // produced most of the first candidate branches on the CTS4466 cohort.
        let cluster = bucket(&[16342231, 16342238, 16342245, 16342253, 16342258, 16342263]);
        assert!(
            candidate_positions(&cluster).is_empty(),
            "the whole cluster goes; there is no basis for electing one call the real one"
        );

        // A lone call keeps its place, and a distant neighbour doesn't drag it down.
        let spread = bucket(&[1_000_000, 2_000_000, 2_000_050]);
        let kept: Vec<i64> = candidate_positions(&spread).into_iter().collect();
        assert_eq!(kept, vec![1_000_000], "only the pair within 100 bp is dropped");
    }

    #[test]
    fn a_position_defining_branches_under_two_parents_is_rejected() {
        // 11311865 was shared by two members under one block *and* two under another. A variant that
        // arose twice cannot mark a new branch, and the laminar check can't see it — it reasons
        // inside a single block.
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
        // The 56.83 Mb case: three branches with *different* member sets inside 567 bp. Three
        // independent lineage events in under a kilobase is not a thing; one repeat unit
        // mis-mapping across several donors is.
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

    /// A two-block tree with one candidate branch and one unplaced member — enough to exercise
    /// every column the export has to get right.
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

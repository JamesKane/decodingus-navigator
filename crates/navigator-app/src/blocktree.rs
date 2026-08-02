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
        let mut at_node: HashMap<i64, Vec<BlockMember>> = HashMap::new();
        let mut unplaced = Vec::new();
        for (guid, name, terminal) in wanted {
            match terminal.as_deref().and_then(|t| index.get(t).copied()) {
                Some(id) => at_node.entry(id).or_default().push(BlockMember {
                    guid,
                    name,
                    // Phase 3 populates these from `donor_private_y`; `None` = not computed.
                    private_novel: None,
                    private_total: None,
                }),
                None => unplaced.push(UnplacedMember { guid, name, terminal }),
            }
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
            })
            .collect();
        // Stable leaf order, so the layout doesn't reshuffle between opens.
        for b in &mut blocks {
            b.members.sort_by(|x, y| (&x.name, x.guid.0).cmp(&(&y.name, y.guid.0)));
        }
        roll_up_subtree_members(&mut blocks);
        let blocks = collapse_blocks(blocks, COLLAPSE_MIN_RUN);
        unplaced.sort_by(|a, b| (&a.name, a.guid.0).cmp(&(&b.name, b.guid.0)));

        Ok(Some(ProjectBlockTree {
            dna,
            blocks,
            unplaced,
            provider: provider.to_string(),
            build_key: build_key.to_string(),
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
        }
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

    #[test]
    fn collapse_is_a_no_op_on_an_empty_tree_or_a_zero_threshold() {
        assert!(collapse_blocks(Vec::new(), 2).is_empty());
        let names: Vec<String> = collapse_blocks(chain(), 0).iter().map(|b| b.name.clone()).collect();
        assert_eq!(names, vec!["root", "A", "B", "C", "D"]);
    }
}

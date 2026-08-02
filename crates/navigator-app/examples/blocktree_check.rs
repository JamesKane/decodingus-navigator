//! Throwaway validation for the project block tree: run `App::project_block_tree` against the live
//! workspace and print the result as an indented tree. Proves the real data path — tree fetch, name
//! index, induced subtree, collapse — on an actual multi-thousand-member cohort.
//!
//! ```bash
//! cargo run -p navigator-app --example blocktree_check -- <project_id>
//! ```

use navigator_app::{App, DnaType};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let project_id: i64 = std::env::args().nth(1).unwrap_or_else(|| "4".into()).parse()?;
    let db = std::path::PathBuf::from(std::env::var("HOME")?).join(".decodingus/navigator-rs.db");
    let app = App::open(&db).await?;

    let started = std::time::Instant::now();
    let Some(tree) = app.project_block_tree(project_id, DnaType::Y).await? else {
        println!("project {project_id}: no members");
        return Ok(());
    };
    let elapsed = started.elapsed();

    let placed: usize = tree.blocks.iter().map(|b| b.members.len()).sum();
    println!(
        "project {project_id} · {} · {} · {} blocks · {placed} placed · {} unplaced · built in {:.1}s",
        tree.provider,
        tree.build_key,
        tree.blocks.len(),
        tree.unplaced.len(),
        elapsed.as_secs_f64(),
    );
    println!(
        "max depth {} · collapsed branches folded {}",
        tree.blocks.iter().map(|b| b.depth).max().unwrap_or(0),
        tree.blocks.iter().map(|b| b.collapsed.len()).sum::<usize>(),
    );

    // Split the unplaced: "no placement at all" is expected (STR-only kits), but "has a terminal
    // this tree doesn't carry" is provider/build skew worth naming.
    let (skew, unplaced_none): (Vec<_>, Vec<_>) = tree.unplaced.iter().partition(|u| u.terminal.is_some());
    println!(
        "unplaced: {} with no Y placement · {} with a terminal absent from this tree",
        unplaced_none.len(),
        skew.len(),
    );
    let mut names: Vec<&str> = skew.iter().filter_map(|u| u.terminal.as_deref()).collect();
    names.sort_unstable();
    names.dedup();
    if !names.is_empty() {
        println!("  unresolved terminals ({}): {}", names.len(), {
            let head: Vec<&str> = names.iter().copied().take(15).collect();
            head.join(", ")
        });
    }

    // Pre-order with depth as indent is exactly how the aggregate is ordered, so this prints itself.
    for b in tree.blocks.iter().take(60) {
        let indent = "  ".repeat(b.depth);
        let folded = if b.collapsed.is_empty() {
            String::new()
        } else {
            format!(" (+{} folded)", b.collapsed.len())
        };
        let members = if b.members.is_empty() {
            String::new()
        } else {
            let names: Vec<&str> = b.members.iter().map(|m| m.name.as_str()).take(4).collect();
            format!(
                "  ← {}{}",
                names.join(", "),
                if b.members.len() > 4 { ", …" } else { "" }
            )
        };
        println!(
            "{indent}{} [{} SNP{}]{folded} ({} below){members}",
            b.name,
            b.loci.len(),
            if b.loci.len() == 1 { "" } else { "s" },
            b.subtree_members,
        );
    }
    if tree.blocks.len() > 60 {
        println!("… {} more blocks", tree.blocks.len() - 60);
    }
    Ok(())
}

//! DUNavigator desktop front end (egui/eframe). The UI is a thin renderer + command
//! dispatcher; a worker thread owns the async `App` and streams events back (see
//! `worker`). Business logic, DB access, and domain decisions live below the UI.

use std::path::PathBuf;

use clap::Parser;

mod charts;
mod cli;
mod i18n;
mod ui;
mod widgets;
mod worker;

use ui::NavigatorApp;

/// `~/.decodingus/navigator-rs.db` (separate from the legacy H2 file). Shared by the GUI and
/// the CLI subcommands so scripted ingestion lands in the same workbench the GUI shows.
pub(crate) fn default_db_path() -> PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    let dir = PathBuf::from(home).join(".decodingus");
    let _ = std::fs::create_dir_all(&dir);
    dir.join("navigator-rs.db")
}

fn main() -> eframe::Result<()> {
    // First-run setup (both GUI and headless): seed bundled ancestry/IBD assets shipped inside the
    // installer image into ~/.decodingus/ancestry/ if missing. No-op on a dev build (no bundle) and
    // on later runs. Done before the CLI dispatch so a scripted analysis also gets the assets.
    let seeded = navigator_app::seed_bundled_assets();
    if seeded.copied > 0 {
        eprintln!("seeded {} bundled asset(s) into the cache", seeded.copied);
    }
    // Likewise seed the chrY private-Y filtering masks (callable mask + cohort-shared exclude) into
    // ~/.decodingus/masks/. These ship gzipped in the repo `assets/masks/`, so a dev build seeds too.
    let seeded_masks = navigator_app::seed_bundled_masks();
    if seeded_masks.copied > 0 {
        eprintln!("seeded {} bundled mask(s) into the cache", seeded_masks.copied);
    }

    // With a subcommand, run headless (ingest/probe) and exit; with none, launch the GUI.
    let parsed = cli::Cli::parse();
    if let Some(command) = parsed.command {
        std::process::exit(cli::run(command));
    }

    let db_path = default_db_path();
    // The reworked layout (subjects table + detail panel + action bar) needs room; the eframe
    // default window is far too small for it. Open at a comfortable size with a sane floor.
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1360.0, 900.0])
            .with_min_inner_size([1024.0, 680.0]),
        ..Default::default()
    };
    eframe::run_native(
        "DUNavigator",
        options,
        Box::new(move |cc| Ok(Box::new(NavigatorApp::new(cc, db_path)))),
    )
}

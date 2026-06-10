//! DUNavigator desktop front end (egui/eframe). The UI is a thin renderer + command
//! dispatcher; a worker thread owns the async `App` and streams events back (see
//! `worker`). Business logic, DB access, and domain decisions live below the UI.

use std::path::PathBuf;

use clap::Parser;

mod cli;
mod i18n;
mod ui;
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
    // With a subcommand, run headless (ingest/probe) and exit; with none, launch the GUI.
    let parsed = cli::Cli::parse();
    if let Some(command) = parsed.command {
        std::process::exit(cli::run(command));
    }

    let db_path = default_db_path();
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "DUNavigator",
        options,
        Box::new(move |cc| Ok(Box::new(NavigatorApp::new(cc, db_path)))),
    )
}

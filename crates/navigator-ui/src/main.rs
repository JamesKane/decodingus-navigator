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
    let dir = navigator_domain::paths::decodingus_dir();
    let _ = std::fs::create_dir_all(&dir);
    dir.join("navigator-rs.db")
}

fn main() -> eframe::Result<()> {
    // Opt this process into the OS keychain — sessions and device keys must survive a restart.
    // This is the *only* place it may be called: everything else (tests, CI, examples) keeps the
    // in-memory default and so can never read or write the user's real credentials. Must run
    // before any `App` is built, since `App::new` reloads the active account.
    navigator_app::use_os_keychain();

    // With a subcommand, run headless (ingest/probe) and exit; with none, launch the GUI.
    let parsed = cli::Cli::parse();

    // First-run setup: seed the bundled ancestry/IBD assets, chrY masks, and HipSTR reference BEDs
    // shipped inside the installer image into ~/.decodingus/ if missing. No-op on later runs.
    //
    // Headless seeds synchronously — the analysis starts at once, and there is no window whose
    // appearance the copy could delay. The GUI seeds on a background thread instead, because the
    // GRCh38 HipSTR BED alone is ~20 MB and that copy would otherwise sit in front of the first
    // frame on exactly the run a new user is watching. `App::open` (on the worker thread) waits for
    // it, so nothing can read a half-seeded cache.
    if let Some(command) = parsed.command {
        let seeded = navigator_app::seed_bundled_all();
        if seeded.copied > 0 {
            eprintln!("seeded {} bundled asset(s) into the cache", seeded.copied);
        }
        std::process::exit(cli::run(command));
    }
    navigator_app::spawn_bundled_seed();

    let db_path = default_db_path();
    // Open at the remembered size (falling back to a comfortable default), with a sane floor. This is
    // only the initial hint — the builder sizes before the UI scale is applied, so `NavigatorApp`
    // re-asserts the remembered size and fits it to the current screen on its first frames.
    let initial_size = navigator_app::AppSettings::load()
        .window_size
        .unwrap_or(ui::DEFAULT_WINDOW);
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size(initial_size)
            .with_min_inner_size(ui::MIN_WINDOW),
        ..Default::default()
    };
    eframe::run_native(
        "DUNavigator",
        options,
        Box::new(move |cc| Ok(Box::new(NavigatorApp::new(cc, db_path)))),
    )
}

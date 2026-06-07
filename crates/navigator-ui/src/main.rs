//! DUNavigator desktop front end (egui/eframe). The UI is a thin renderer + command
//! dispatcher; a worker thread owns the async `App` and streams events back (see
//! `worker`). Business logic, DB access, and domain decisions live below the UI.

use std::path::PathBuf;

mod i18n;
mod ui;
mod worker;

use ui::NavigatorApp;

/// `~/.decodingus/navigator-rs.db` (separate from the legacy H2 file).
fn default_db_path() -> PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    let dir = PathBuf::from(home).join(".decodingus");
    let _ = std::fs::create_dir_all(&dir);
    dir.join("navigator-rs.db")
}

fn main() -> eframe::Result<()> {
    let db_path = default_db_path();
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "DUNavigator",
        options,
        Box::new(move |cc| Ok(Box::new(NavigatorApp::new(cc, db_path)))),
    )
}

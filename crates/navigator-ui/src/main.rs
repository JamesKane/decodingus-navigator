//! DUNavigator desktop front end (egui/eframe). Thin by design: renders immutable
//! view-state and dispatches commands to `navigator-app`; long-running analysis runs
//! off-thread and streams progress back over a channel the repaint loop reads. No
//! business logic, DB calls, or domain decisions here. egui shell lands in phase 5.

fn main() {
    println!("DUNavigator (Rust rewrite) — workspace skeleton. egui UI arrives in roadmap phase 5.");
}

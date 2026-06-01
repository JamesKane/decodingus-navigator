//! Navigator application/command layer — the antidote to the `WorkbenchViewModel`
//! god object. Orchestrates `navigator-store`, `navigator-analysis`, and
//! `navigator-sync` behind a command/query API the UI dispatches to. Holds the
//! policy that today's dialogs embed (e.g. fingerprint-match resolution); the UI only
//! renders view-state and prompts when genuinely needed.
//!
//! No UI types, no widget code below here. Implemented in roadmap phase 4+.

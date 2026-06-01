//! Navigator domain types — the desktop-only aggregates that `du-domain` does not
//! cover: `SequenceRun`, `Alignment`, `AnalysisArtifact`, `YProfile`, IBD, and the
//! `Workspace`/`Project` aggregate. Pure types, zero IO; this is the bottom of the
//! dependency graph (`ui → app → {analysis, store, sync} → domain`).
//!
//! Re-exports `du-domain` so every higher crate sees a single domain surface and the
//! `HaplogroupResult`/`ScoredHaplogroup`/`HaplogroupAssignments` triplication collapses
//! to one shared type. Types land in roadmap phase 4; this is the phase-1 skeleton.

pub use du_domain;

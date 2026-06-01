//! Navigator local persistence — SQLite via `sqlx`, replacing the H2/Slick layer.
//! Query-module-per-aggregate (the `du-db` pattern); complex children modelled as
//! proper rows, not 22-tuple JSONB blobs. `Json<T>` reserved for AT Proto record
//! snapshots. Versioned, reversible `sqlx migrate` migrations.
//!
//! Persisted state is authoritative; the UI reads a projection — no imperative
//! `_workspace.value = …` mutation racing async writes. Implemented in roadmap phase 4.

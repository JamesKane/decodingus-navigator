//! Navigator sync — completes the stubbed `AsyncSyncService`: PDS push/pull with
//! explicit retry/backoff, conflict policy, and an offline indicator. Also publishes
//! per-sample coverage summaries as public PDS records and submits variant/branch
//! proposals to the AppView curation API (post-AppView-scope-reduction work).
//!
//! Authenticates via `du-atproto` OAuth as a public/native client: PKCE only, native
//! `client-metadata.json`, loopback redirect (`http://127.0.0.1:<port>/callback`),
//! tokens in the OS keychain. Implemented in roadmap phase 6.

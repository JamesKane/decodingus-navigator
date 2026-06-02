//! Navigator sync — AT Proto OAuth (public/native client) + PDS push/pull.
//!
//! Authenticates as a **public/native client** (plan §7): PKCE only, a loopback redirect
//! (`http://127.0.0.1:<port>/callback`), DPoP-bound tokens in the OS keychain. Reuses the
//! shared `du-atproto` OAuth primitives. PDS publishing (coverage summaries, variant
//! proposals) and the completed AsyncSync build on the authenticated [`tokens::Session`].

pub mod error;
pub mod oauth;
pub mod publish;
pub mod records;
pub mod tokens;

pub use error::SyncError;
pub use oauth::{login, OAuthConfig};
pub use publish::{dev_http_client, PdsClient, RecordRef};
pub use records::{
    CoverageSummaryRecord, PrivateVariantsRecord, VariantCallEntry, COVERAGE_SUMMARY_COLLECTION,
    PRIVATE_VARIANTS_COLLECTION,
};
pub use tokens::{Session, TokenStore};

//! Navigator sync — AT Proto OAuth (public/native client) + PDS push/pull.
//!
//! Authenticates as a **public/native client** (plan §7): PKCE only, a loopback redirect
//! (`http://127.0.0.1:<port>/callback`), DPoP-bound tokens in the OS keychain. Reuses the
//! shared `du-atproto` OAuth primitives. PDS publishing (coverage summaries, variant
//! proposals) and the completed AsyncSync build on the authenticated [`tokens::Session`].

pub mod device_key;
pub mod error;
pub mod oauth;
pub mod publish;
pub mod records;
pub mod sync;
pub mod tokens;

pub use device_key::{DeviceKey, DEVICE_KEY_COLLECTION};
pub use error::SyncError;
pub use oauth::{login, login_default, refresh, OAuthConfig};
pub use publish::{dev_http_client, PdsClient, RecordRef};
pub use sync::{AsyncSync, RetryPolicy};
pub use records::{
    AuditEntryRecord, HaplogroupReconciliationRecord, HeteroplasmyObservationRecord,
    IdentityVerificationRecord, ManualOverrideRecord, PrivateVariantsRecord, RecordMeta,
    ReconciliationStatusRecord, RunHaplogroupCallRecord, VariantCallEntry,
    HAPLOGROUP_RECONCILIATION_COLLECTION, PRIVATE_VARIANTS_COLLECTION,
};
pub use tokens::{Session, TokenStore};

// Federated atproto wire records — the single source of truth lives in the shared
// `du-domain::fed` module so the AppView's Jetstream consumer cannot drift from us.
// (Its `RecordMeta` is intentionally not re-exported to avoid colliding with the
// reconciliation record's `RecordMeta`; `::new` builds it internally.)
pub use du_domain::fed::{
    AlignmentRecord, BiosampleRecord, CoverageMetrics, PopulationBreakdownRecord,
    PopulationComponent as FedPopulationComponent, SequenceRunRecord,
    SuperPopulationSummary as FedSuperPopulationSummary, WireF64, NS_ALIGNMENT, NS_BIOSAMPLE,
    NS_POPULATION_BREAKDOWN, NS_SEQUENCERUN,
};

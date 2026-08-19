//! Navigator sync: AT Proto OAuth for a public/native client, plus PDS push and pull.
//!
//! It authenticates as a **public/native client** (plan §7), with PKCE alone, a loopback redirect
//! (`http://127.0.0.1:<port>/callback`), and DPoP-bound tokens in the OS keychain. It reuses the
//! shared `du-atproto` OAuth primitives.
//!
//! The PDS push, which sends coverage summaries and variant proposals, and the completed
//! AsyncSync, both build on the authenticated [`tokens::Session`].

pub mod device_key;
pub mod error;
pub mod exchange;
pub mod oauth;
pub mod publish;
pub mod records;
pub mod recruitment;
pub mod secret_store;
pub mod social;
pub mod sync;
pub mod tokens;

pub use device_key::{DeviceKey, DEVICE_KEY_COLLECTION};
pub use error::SyncError;
pub use exchange::{Envelope, EphemeralKey, ExchangeKey};
pub use oauth::{login, login_default, refresh, OAuthConfig};
pub use publish::{dev_http_client, PdsClient, RecordRef, RemoteRecord};
pub use records::{
    AncestralOriginRecord, AuditEntryRecord, HaplogroupReconciliationRecord, HeteroplasmyObservationRecord,
    IdentityVerificationRecord, ManualOverrideRecord, OriginExternalId, PrivateVariantsRecord,
    ReconciliationStatusRecord, RecordMeta, RunHaplogroupCallRecord, VariantCallEntry, ANCESTRAL_ORIGIN_COLLECTION,
    HAPLOGROUP_RECONCILIATION_COLLECTION, PRIVATE_VARIANTS_COLLECTION,
};
pub use secret_store::{os_keychain_enabled, use_os_keychain};
pub use sync::{AsyncSync, RetryPolicy};
pub use tokens::{Session, TokenStore};

// Federated atproto wire records. The single source of truth is the shared `du-domain::fed`
// module, so the AppView's Jetstream consumer can not become different from us.
//
// This module does not re-export that module's `RecordMeta`, on purpose. It would then have the
// same name as the reconciliation record's `RecordMeta`. `::new` builds it inside the module.
pub use du_domain::fed::{
    AlignmentRecord, BiosampleRecord, ContigMetrics, CoverageMetrics, FeedPostRecord, PopulationBreakdownRecord,
    PopulationComponent as FedPopulationComponent, PostRef, ReplyRef, SequenceRunRecord,
    SuperPopulationSummary as FedSuperPopulationSummary, WireF64, NS_ALIGNMENT, NS_BIOSAMPLE, NS_FEED_POST,
    NS_POPULATION_BREAKDOWN, NS_SEQUENCERUN,
};

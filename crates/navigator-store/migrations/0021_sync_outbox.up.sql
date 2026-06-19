-- Sync durability (gap §5): a persistent outbox so PDS publishes survive restart and
-- offline→online. A publish enqueues a fully-built record here; a background drain pushes it with
-- exponential backoff. A transient/offline failure is rescheduled (not lost); a successful push is
-- recorded in sync_history and the outbox row removed.

CREATE TABLE sync_outbox (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    -- Destination PDS account (DID). Only drained while signed in as this DID.
    account_did   TEXT    NOT NULL,
    -- Human label for the UI / history (coverage, ancestry, biosample, seqrun, variants,
    -- reconciliation).
    kind          TEXT    NOT NULL,
    -- Stable id of the published thing (dedup key), e.g. "alignment:12",
    -- "ancestry:12:ADMIXTURE", "variants:12:chrY".
    entity_ref    TEXT    NOT NULL,
    -- AT-Proto collection NSID.
    collection    TEXT    NOT NULL,
    -- Explicit record key for idempotent singletons; NULL = server-assigned TID.
    rkey          TEXT,
    -- The JSON record to POST.
    payload       TEXT    NOT NULL,
    -- PENDING (awaiting a drain) | FAILED (non-transient error, not auto-retried). A successful
    -- push removes the row (its outcome lands in sync_history).
    status        TEXT    NOT NULL DEFAULT 'PENDING',
    attempt_count INTEGER NOT NULL DEFAULT 0,
    -- ISO-8601; NULL = ready now. Set to a future time on a transient failure (backoff).
    next_retry_at TEXT,
    last_error    TEXT,
    created_at    TEXT    NOT NULL,
    updated_at    TEXT    NOT NULL,
    -- Re-publishing the same thing coalesces onto one row (newest payload wins).
    UNIQUE(account_did, collection, entity_ref)
);

CREATE INDEX idx_sync_outbox_ready ON sync_outbox(account_did, status, next_retry_at);

-- Append-only audit trail of completed push attempts (success or terminal failure).
CREATE TABLE sync_history (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    account_did   TEXT    NOT NULL,
    kind          TEXT    NOT NULL,
    entity_ref    TEXT    NOT NULL,
    collection    TEXT    NOT NULL,
    direction     TEXT    NOT NULL DEFAULT 'PUSH', -- PULL deferred (no remote→local sync yet)
    status        TEXT    NOT NULL,                -- SUCCESS | FAILED
    at_uri        TEXT,
    at_cid        TEXT,
    attempt_count INTEGER NOT NULL DEFAULT 0,
    error         TEXT,
    created_at    TEXT    NOT NULL
);

CREATE INDEX idx_sync_history_time ON sync_history(account_did, created_at);

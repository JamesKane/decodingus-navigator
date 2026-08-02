-- The federated-IBD request ledger: one durable row per matching conversation, from the moment we
-- ask for (or receive) an introduction until the encrypted exchange produces a result. Before this
-- table the only persisted state was the *completed* exchange (`ibd_exchange_result`), so a restart
-- lost every in-flight request — "I asked X and I'm waiting on their consent" existed only in UI
-- memory.
--
-- Keyed by the broker's `request_uri` (`urn:ibd:<sha256>` for a suggestion-mediated introduction,
-- `exchange:<uuid>` for a direct one), which is stable and idempotent per (caller, candidate).
--
-- `my_sample_ref` / `partner_sample_ref` are the **AppView** `core.biosample` guids (from the
-- suggestion's `target_sample_guid` / `suggested_sample_guid`) — not local subject guids. Both are
-- needed to attest a completed comparison, which the AppView gates on sample ownership.
-- `biosample_guid` is the *local* subject whose dosages we exchange.
--
-- `consent_given` records only OUR OWN decision (NULL = undecided). The broker is symmetric-blind:
-- a partner's refusal is never reported, it simply never becomes READY.
--
-- PII-free: DIDs, opaque sample handles, and a lifecycle status.
CREATE TABLE ibd_request (
    request_uri        TEXT PRIMARY KEY,
    direction          TEXT NOT NULL,
    purpose            TEXT NOT NULL DEFAULT '',
    status             TEXT NOT NULL,
    partner_did        TEXT,
    session_id         TEXT,
    biosample_guid     TEXT REFERENCES biosample(guid),
    my_sample_ref      TEXT,
    partner_sample_ref TEXT,
    consent_given      INTEGER,
    consent_at         TEXT,
    attested_at        TEXT,
    last_error         TEXT,
    created_at         TEXT NOT NULL,
    updated_at         TEXT NOT NULL
);

CREATE INDEX ix_ibd_request_status ON ibd_request (status);
CREATE INDEX ix_ibd_request_biosample ON ibd_request (biosample_guid);

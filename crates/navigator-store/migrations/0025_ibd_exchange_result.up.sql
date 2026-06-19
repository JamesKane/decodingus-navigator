-- Federated IBD exchange results (gap §4). One row per completed exchange session: the locally
-- computed match summary + both peers' signed attestations + the segment list (JSON). Keyed by the
-- session id (one result per session). `biosample_guid` is the local subject whose dosages were
-- exchanged; `partner_did` is the peer. `agreed` ⇒ the partner's signature verified AND both summary
-- hashes matched. PII-free: DIDs + opaque sample refs + segment cM only — no identifiers.
CREATE TABLE ibd_exchange_result (
    session_id          TEXT PRIMARY KEY,
    request_uri         TEXT NOT NULL,
    my_did              TEXT NOT NULL,
    partner_did         TEXT NOT NULL,
    biosample_guid      TEXT NOT NULL REFERENCES biosample(guid),
    partner_sample_ref  TEXT,
    total_shared_cm     REAL NOT NULL DEFAULT 0,
    segment_count       INTEGER NOT NULL DEFAULT 0,
    longest_segment_cm  REAL NOT NULL DEFAULT 0,
    relationship        TEXT NOT NULL DEFAULT '',
    agreed              INTEGER NOT NULL DEFAULT 0,
    segments            TEXT NOT NULL DEFAULT '[]',
    my_attestation      TEXT NOT NULL DEFAULT '{}',
    partner_attestation TEXT NOT NULL DEFAULT '{}',
    created_at          TEXT NOT NULL
);

CREATE INDEX ix_ibd_exchange_result_biosample ON ibd_exchange_result (biosample_guid);

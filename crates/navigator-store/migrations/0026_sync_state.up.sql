-- The PDS-assigned identity of each published entity, so re-publishing UPDATES that record
-- (putRecord at the kept rkey) instead of creating a duplicate. `rkey` is the TID the PDS assigned on
-- first publish (parsed from the returned at-uri). `payload_hash` is the sha256 of the published JSON
-- at push time — lets a PULL detect local edits since the last push (for conflict resolution). One row
-- per (account, entity) keyed by the same stable `entity_ref` the outbox coalesces on.
CREATE TABLE sync_state (
    account_did   TEXT NOT NULL,
    entity_ref    TEXT NOT NULL,
    kind          TEXT NOT NULL,
    collection    TEXT NOT NULL,
    rkey          TEXT NOT NULL,
    at_uri        TEXT NOT NULL,
    at_cid        TEXT NOT NULL,
    payload_hash  TEXT NOT NULL,
    pushed_at     TEXT NOT NULL,
    PRIMARY KEY (account_did, entity_ref)
);

CREATE INDEX ix_sync_state_collection ON sync_state (account_did, collection);

-- Cached chromosome painting (local-ancestry segments) for a subject, keyed to the autosomal
-- consensus it was computed from. `consensus_sig` is the consensus's `last_reconciled_at`; when the
-- consensus is rebuilt the signature changes and the cached painting is considered stale (the app
-- compares them on read and recomputes on mismatch). `segments` is the painted segment list as JSON.
CREATE TABLE consensus_painting (
    biosample_guid TEXT PRIMARY KEY REFERENCES biosample(guid),
    consensus_sig  TEXT NOT NULL,
    segments       TEXT NOT NULL,
    painted_at     TEXT NOT NULL
);

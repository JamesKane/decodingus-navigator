-- Cached archaic (Neanderthal / Denisovan) Tier-A marker count for a subject, keyed to the
-- autosomal consensus it was computed from. `consensus_sig` is the consensus's `last_reconciled_at`;
-- when the consensus is rebuilt the signature changes and the cached result is stale (the app
-- compares on read and recomputes on mismatch). `archaic` is the full ArchaicMarkerResult as JSON.
-- Mirrors consensus_roh / consensus_painting.
CREATE TABLE consensus_archaic (
    biosample_guid TEXT PRIMARY KEY REFERENCES biosample(guid),
    consensus_sig  TEXT NOT NULL,
    archaic        TEXT NOT NULL,
    computed_at    TEXT NOT NULL
);

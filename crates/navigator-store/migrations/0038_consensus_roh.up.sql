-- Cached runs-of-homozygosity (ROH) result for a subject, keyed to the autosomal consensus it was
-- computed from. `consensus_sig` is the consensus's `last_reconciled_at`; when the consensus is
-- rebuilt the signature changes and the cached ROH is considered stale (the app compares them on
-- read and recomputes on mismatch). `roh` is the full RohResult (segments + summary) as JSON.
CREATE TABLE consensus_roh (
    biosample_guid TEXT PRIMARY KEY REFERENCES biosample(guid),
    consensus_sig  TEXT NOT NULL,
    roh            TEXT NOT NULL,
    computed_at    TEXT NOT NULL
);

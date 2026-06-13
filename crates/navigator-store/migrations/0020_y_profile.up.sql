-- Persisted snapshot of a subject's multi-source Y-variant profile (one per biosample). The full
-- reconciled profile (variants + per-source provenance + summary) lives in `payload` as JSON; the
-- scalar columns mirror the header for quick listing without decoding. Rebuilt by "Build/Refresh"
-- in the Y-DNA tab — re-genotyping each source is expensive, so this snapshot lets the tab load
-- instantly. Cross-profile SNP queries (a future match list) would want a relational variant table.
CREATE TABLE y_profile (
    biosample_guid       TEXT PRIMARY KEY REFERENCES biosample(guid),
    consensus_haplogroup TEXT,
    overall_confidence   REAL NOT NULL DEFAULT 0,
    source_count         INTEGER NOT NULL DEFAULT 0,
    total                INTEGER NOT NULL DEFAULT 0,
    confirmed            INTEGER NOT NULL DEFAULT 0,
    novel                INTEGER NOT NULL DEFAULT 0,
    conflict             INTEGER NOT NULL DEFAULT 0,
    single_source        INTEGER NOT NULL DEFAULT 0,
    tree_provider        TEXT,
    payload              TEXT NOT NULL,
    last_reconciled_at   TEXT NOT NULL
);

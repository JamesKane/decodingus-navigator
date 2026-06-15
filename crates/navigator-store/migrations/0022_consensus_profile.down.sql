-- Reverse 0022: restore the Y-only y_profile table (0020) and copy the Y rows back.
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

INSERT INTO y_profile (biosample_guid, consensus_haplogroup, overall_confidence, source_count,
    total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at)
SELECT biosample_guid, consensus_label, overall_confidence, source_count,
    total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at
FROM consensus_profile WHERE dna_type = 'Y';

DROP TABLE consensus_profile;

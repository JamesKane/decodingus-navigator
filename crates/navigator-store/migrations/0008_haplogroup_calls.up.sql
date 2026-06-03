-- Per-source Y/mtDNA haplogroup calls for a subject (donor-level reconciliation input).
-- One row per (biosample, dna_type, source); lineage stored tab-joined.
CREATE TABLE haplogroup_call (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    dna_type       TEXT NOT NULL,   -- 'Y' | 'Mt'
    source_key     TEXT NOT NULL,   -- dedup key, e.g. 'aln:5'
    source_label   TEXT NOT NULL,
    haplogroup     TEXT NOT NULL,
    lineage        TEXT NOT NULL,   -- tab-joined root->terminal
    score          REAL NOT NULL,
    matched        INTEGER NOT NULL,
    expected       INTEGER NOT NULL,
    UNIQUE(biosample_guid, dna_type, source_key)
);
CREATE INDEX idx_haplogroup_call_biosample ON haplogroup_call(biosample_guid, dna_type);

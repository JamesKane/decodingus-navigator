-- Vendor mtDNA FASTA sequences for a subject. The sequence is ~16,569 bp, small enough
-- to store inline as TEXT.
CREATE TABLE mtdna_sequence (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid   TEXT NOT NULL REFERENCES biosample(guid),
    defline          TEXT,
    sequence         TEXT NOT NULL,
    n_count          INTEGER NOT NULL,
    source_file_name TEXT
);

CREATE INDEX idx_mtdna_sequence_biosample ON mtdna_sequence(biosample_guid);

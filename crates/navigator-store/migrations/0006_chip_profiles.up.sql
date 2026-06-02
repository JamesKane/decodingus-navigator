-- Genotyping-array (chip) profiles: the QC summary of a vendor raw-data export. We store
-- the summary (call/no-call/het + per-region counts), not the ~600k individual genotypes.
CREATE TABLE chip_profile (
    id                       INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid           TEXT NOT NULL REFERENCES biosample(guid),
    provider                 TEXT NOT NULL,
    chip_version             TEXT,
    total_markers_possible   INTEGER NOT NULL,
    total_markers_called     INTEGER NOT NULL,
    no_call_rate             REAL NOT NULL,
    het_rate                 REAL,
    y_markers_called         INTEGER NOT NULL,
    mt_markers_called        INTEGER NOT NULL,
    autosomal_markers_called INTEGER NOT NULL,
    source_file_name         TEXT
);

CREATE INDEX idx_chip_profile_biosample ON chip_profile(biosample_guid);

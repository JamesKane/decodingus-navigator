-- Initial workspace schema. Proper relational tables (not Slick 22-tuple JSONB blobs);
-- the only JSON is analysis_artifact.payload (a versioned, recomputable result cache).

CREATE TABLE project (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    name          TEXT NOT NULL,
    description   TEXT,
    administrator TEXT NOT NULL
);

CREATE TABLE biosample (
    guid             TEXT PRIMARY KEY NOT NULL,
    sample_accession TEXT,
    donor_identifier TEXT NOT NULL,
    description      TEXT,
    center_name      TEXT,
    sex              TEXT,
    project_id       INTEGER REFERENCES project(id)
);

CREATE TABLE sequence_run (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid   TEXT NOT NULL REFERENCES biosample(guid),
    platform_name    TEXT NOT NULL,
    instrument_model TEXT,
    test_type        TEXT NOT NULL,
    library_layout   TEXT,
    total_reads      INTEGER,
    pf_reads_aligned INTEGER,
    mean_read_length REAL,
    mean_insert_size REAL
);

CREATE TABLE alignment (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_run_id INTEGER NOT NULL REFERENCES sequence_run(id),
    reference_build TEXT NOT NULL,
    aligner         TEXT NOT NULL,
    variant_caller  TEXT
);

CREATE TABLE analysis_artifact (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    alignment_id      INTEGER NOT NULL REFERENCES alignment(id),
    kind              TEXT NOT NULL,
    algorithm_version TEXT NOT NULL,
    created_at        TEXT NOT NULL,
    payload           TEXT NOT NULL,
    UNIQUE (alignment_id, kind, algorithm_version)
);

CREATE INDEX idx_biosample_project   ON biosample(project_id);
CREATE INDEX idx_sequence_run_sample ON sequence_run(biosample_guid);
CREATE INDEX idx_alignment_run       ON alignment(sequence_run_id);
CREATE INDEX idx_artifact_alignment  ON analysis_artifact(alignment_id);

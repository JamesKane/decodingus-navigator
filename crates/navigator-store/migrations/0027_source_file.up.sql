-- Content-hash file identity (gap §5-p2): imported files are tracked by their SHA-256, not their path,
-- so moving/renaming a BAM just updates the path instead of orphaning its analyses. Ports the legacy
-- `SourceFileRepository`. `content_sha256` is the stable key; `is_accessible` tracks moved/deleted
-- files; `alignment_id` links the file to its analysis (nullable until linked).
CREATE TABLE source_file (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    content_sha256  TEXT NOT NULL UNIQUE,
    file_path       TEXT,
    file_size       INTEGER,
    file_format     TEXT,
    alignment_id    INTEGER REFERENCES alignment(id),
    is_accessible   INTEGER NOT NULL DEFAULT 1,
    last_verified_at TEXT,
    created_at      TEXT NOT NULL,
    updated_at      TEXT NOT NULL
);

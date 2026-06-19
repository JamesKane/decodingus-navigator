-- A subject's imported SNP variant calls, grouped into a named set (one per VCF/CSV
-- import). SNP-only: ref/alt are single bases.
CREATE TABLE variant_set (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    source_label   TEXT NOT NULL
);

CREATE TABLE variant_call (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    variant_set_id INTEGER NOT NULL REFERENCES variant_set(id),
    contig         TEXT NOT NULL,
    position       INTEGER NOT NULL,
    reference      TEXT NOT NULL,
    alternate      TEXT NOT NULL,
    rs_id          TEXT,
    genotype       TEXT
);

CREATE INDEX idx_variant_set_biosample ON variant_set(biosample_guid);
CREATE INDEX idx_variant_call_set ON variant_call(variant_set_id);

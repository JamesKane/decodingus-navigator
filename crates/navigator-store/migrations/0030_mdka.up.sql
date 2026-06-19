-- Most Distant Known Ancestor (FTDNA project-import design §4.3). Genealogy, not genetics: the
-- earliest documented ancestor on a lineage, with name, dates, place, and optional geocoded
-- coordinates. Vendor-agnostic (FTDNA import, manual entry, or another vendor), one per lineage.
--
-- PII / project-shared-private: the most sensitive data in the importer — it names real people and
-- places. NEVER federated/published, NEVER stored in AppView; it may only move admin-to-admin over
-- the encrypted Edge-to-Edge channel.
CREATE TABLE mdka (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid  TEXT NOT NULL REFERENCES biosample(guid),
    lineage         TEXT NOT NULL,        -- 'Y' | 'Mt' | 'Auto'
    ancestor_name   TEXT,
    birth_year      INTEGER,
    death_year      INTEGER,
    origin_place    TEXT,                 -- free-text place as entered (e.g. "Creegh South, Co. Clare, Ireland")
    origin_country  TEXT,                 -- normalized country (for grouping/maps)
    latitude        REAL,                 -- geocoded; nullable
    longitude       REAL,
    source          TEXT,                 -- 'FTDNA' | 'MANUAL' | ...
    notes           TEXT,
    updated_at      TEXT NOT NULL,        -- ISO-8601
    UNIQUE (biosample_guid, lineage)      -- one MDKA per line per Subject
);
CREATE INDEX idx_mdka_biosample ON mdka(biosample_guid);

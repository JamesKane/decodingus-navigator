-- Y-STR profiles: a subject's short-tandem-repeat marker calls, grouped by the panel/
-- test that produced them. Markers are kept as text (multi-copy markers carry several
-- alleles, e.g. "16-17").
CREATE TABLE str_profile (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    panel_name     TEXT NOT NULL,
    provider       TEXT,
    source         TEXT
);

CREATE TABLE str_marker (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    str_profile_id INTEGER NOT NULL REFERENCES str_profile(id),
    marker         TEXT NOT NULL,
    value          TEXT NOT NULL
);

CREATE INDEX idx_str_profile_biosample ON str_profile(biosample_guid);
CREATE INDEX idx_str_marker_profile ON str_marker(str_profile_id);

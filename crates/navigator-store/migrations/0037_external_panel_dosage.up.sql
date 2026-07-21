-- Per-source autosomal 1240K panel dosages imported from a trusted external caller (a GATK4 / 1240K
-- EIGENSTRAT call set), resolved once to canonical CHM13 sites and stored so the autosomal consensus
-- can pool them WITHOUT decoding a CRAM — the autosomal counterpart to the Y/mt sidecar fast path.
--
-- `dosages` is the resolved `Vec<SiteGenotype>` (CHM13-oriented, dosage 0/1/2) as JSON, the same
-- shape `ibd_panel_dosages` produces for a chip/alignment. One row per (biosample, source_label);
-- `provenance` mirrors the haplogroup_call tier ('external' here) for future precedence work.
CREATE TABLE external_panel_dosage (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL,
    source_label   TEXT NOT NULL,
    provenance     TEXT NOT NULL DEFAULT 'external',
    panel_sig      TEXT,
    site_count     INTEGER NOT NULL DEFAULT 0,
    dosages        TEXT NOT NULL,
    created_at     TEXT NOT NULL,
    UNIQUE (biosample_guid, source_label)
);

CREATE INDEX idx_external_panel_dosage_biosample ON external_panel_dosage (biosample_guid);

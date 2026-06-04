-- Ancestry estimate for an alignment (one row per biosample + alignment + panel). The
-- ranked components and super-population summary are stored as JSON; pca_json is reserved
-- for the phase-2 PCA coordinates.
CREATE TABLE ancestry_result (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid    TEXT NOT NULL REFERENCES biosample(guid),
    alignment_id      INTEGER NOT NULL,
    panel_type        TEXT NOT NULL,   -- 'aims' | 'genome-wide'
    reference_version TEXT NOT NULL,
    confidence_level  REAL NOT NULL,
    snps_analyzed     INTEGER NOT NULL,
    snps_with_genotype INTEGER NOT NULL,
    components_json   TEXT NOT NULL,   -- Vec<PopulationComponent>
    super_pop_json    TEXT NOT NULL,   -- Vec<SuperPopulationSummary>
    pca_json          TEXT,            -- Option<Vec<f64>> (phase 2)
    UNIQUE(biosample_guid, alignment_id, panel_type)
);
CREATE INDEX idx_ancestry_result_biosample ON ancestry_result(biosample_guid);
CREATE INDEX idx_ancestry_result_alignment ON ancestry_result(alignment_id);

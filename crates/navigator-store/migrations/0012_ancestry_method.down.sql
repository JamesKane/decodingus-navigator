-- Revert to the panel_type-keyed ancestry_result (drops the method column).
DROP TABLE IF EXISTS ancestry_result;
CREATE TABLE ancestry_result (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid    TEXT NOT NULL REFERENCES biosample(guid),
    alignment_id      INTEGER NOT NULL,
    panel_type        TEXT NOT NULL,
    reference_version TEXT NOT NULL,
    confidence_level  REAL NOT NULL,
    snps_analyzed     INTEGER NOT NULL,
    snps_with_genotype INTEGER NOT NULL,
    components_json   TEXT NOT NULL,
    super_pop_json    TEXT NOT NULL,
    pca_json          TEXT,
    UNIQUE(biosample_guid, alignment_id, panel_type)
);
CREATE INDEX idx_ancestry_result_biosample ON ancestry_result(biosample_guid);
CREATE INDEX idx_ancestry_result_alignment ON ancestry_result(alignment_id);

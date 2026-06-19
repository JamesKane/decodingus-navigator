-- Add the estimator `method` to ancestry_result and make it the per-alignment
-- discriminator, so multiple estimates (e.g. ADMIXTURE + PCA_PROJECTION_GMM) can
-- coexist for one alignment. `panel_type` reverts to its true meaning (panel kind:
-- 'aims' | 'genome-wide'). No data continuity (pre-beta), so the table is recreated.
DROP TABLE IF EXISTS ancestry_result;
CREATE TABLE ancestry_result (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid    TEXT NOT NULL REFERENCES biosample(guid),
    alignment_id      INTEGER NOT NULL,
    method            TEXT NOT NULL,   -- 'AF_LIKELIHOOD' | 'ADMIXTURE' | 'PCA_PROJECTION_GMM'
    panel_type        TEXT NOT NULL,   -- 'aims' | 'genome-wide'
    reference_version TEXT NOT NULL,
    confidence_level  REAL NOT NULL,
    snps_analyzed     INTEGER NOT NULL,
    snps_with_genotype INTEGER NOT NULL,
    components_json   TEXT NOT NULL,   -- Vec<PopulationComponent>
    super_pop_json    TEXT NOT NULL,   -- Vec<SuperPopulationSummary>
    pca_json          TEXT,            -- Option<Vec<f64>>
    UNIQUE(biosample_guid, alignment_id, method)
);
CREATE INDEX idx_ancestry_result_biosample ON ancestry_result(biosample_guid);
CREATE INDEX idx_ancestry_result_alignment ON ancestry_result(alignment_id);

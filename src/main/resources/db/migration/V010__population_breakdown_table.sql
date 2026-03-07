-- ============================================
-- Population Breakdown Table
-- Stores ancestry composition analysis results as first-class Atmosphere Lexicon records.
-- NSID: com.decodingus.atmosphere.populationBreakdown
-- ============================================

CREATE TABLE population_breakdown (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,

    -- Analysis parameters
    analysis_method VARCHAR(50) NOT NULL DEFAULT 'PCA_PROJECTION_GMM',
    panel_type VARCHAR(20) NOT NULL,
    reference_populations VARCHAR(50) DEFAULT '1000G_HGDP_v1',

    -- SNP quality metrics
    snps_analyzed INT NOT NULL,
    snps_with_genotype INT NOT NULL,
    snps_missing INT NOT NULL,
    confidence_level DOUBLE NOT NULL,

    -- Population components (JSON array of PopulationComponent)
    components JSON NOT NULL,

    -- Super-population summary (JSON array of SuperPopulationSummary)
    super_population_summary JSON NOT NULL,

    -- Optional PCA coordinates for visualization (JSON array of doubles)
    pca_coordinates JSON,

    -- Pipeline info
    analysis_date TIMESTAMP,
    pipeline_version VARCHAR(50),
    reference_version VARCHAR(50),

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- One breakdown per biosample per panel type
    CONSTRAINT uq_breakdown_biosample_panel UNIQUE (biosample_id, panel_type),
    CONSTRAINT chk_panel_type CHECK (panel_type IN ('aims', 'genome-wide'))
);

CREATE INDEX idx_population_breakdown_biosample ON population_breakdown(biosample_id);
CREATE INDEX idx_population_breakdown_panel ON population_breakdown(panel_type);
CREATE INDEX idx_population_breakdown_sync ON population_breakdown(sync_status);
